use std::collections::HashMap;
use std::io::{self, BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// A parsed variant record: sort key + raw line bytes.
struct Record {
    contig_rank: usize,
    pos: i64,
    ref_allele: Vec<u8>,
    alt_allele: Vec<u8>,
    line: Vec<u8>,
}

/// Byte index of the `n`-th tab (1-based) in `line`, or None.
fn nth_tab(line: &[u8], n: usize) -> Option<usize> {
    line.iter()
        .enumerate()
        .filter(|&(_, &b)| b == b'\t')
        .map(|(i, _)| i)
        .nth(n - 1)
}

/// Extract the contig ID from a `##contig=<ID=...,…>` header line.
/// Returns `None` if the line is not a contig meta-info line.
fn parse_contig_id(line: &[u8]) -> Option<&[u8]> {
    let prefix = b"##contig=<";
    if !line.starts_with(prefix) {
        return None;
    }
    let inner = &line[prefix.len()..];
    let id_marker = b"ID=";
    let id_pos = inner
        .windows(id_marker.len())
        .position(|w| w == id_marker)?;
    let after_id = &inner[id_pos + id_marker.len()..];
    let end = after_id
        .iter()
        .position(|&b| b == b',' || b == b'>')
        .unwrap_or(after_id.len());
    Some(&after_id[..end])
}

/// Sort a VCF by (contig header order, POS, REF, ALT) — matching `bcftools sort` ordering.
///
/// Chromosome order is determined by the `##contig=<ID=...>` header lines. Records on
/// contigs absent from the header are sorted after all declared contigs, in appearance
/// order. Within a contig, records are sorted by POS ascending then REF then ALT
/// lexicographically (matching bcftools sort tie-breaking).
pub fn sort_vcf(input: &Path, output: &mut dyn io::Write) -> Result<u64> {
    let raw = std::fs::read(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    let data = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut d = Vec::new();
        flate2::read::MultiGzDecoder::new(&raw[..])
            .read_to_end(&mut d)
            .map_err(RsomicsError::Io)?;
        d
    } else {
        raw
    };

    let mut header_lines: Vec<Vec<u8>> = Vec::new();
    // contig ID → rank (0-based, insertion order)
    let mut contig_rank: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut records: Vec<Record> = Vec::new();
    // rank assigned to contigs not in header (assigned on first appearance)
    let mut unknown_contigs: HashMap<Vec<u8>, usize> = HashMap::new();

    for raw_line in data.split(|&b| b == b'\n') {
        let line = match raw_line.last() {
            Some(b'\r') => &raw_line[..raw_line.len() - 1],
            _ => raw_line,
        };
        if line.is_empty() {
            continue;
        }
        if line[0] == b'#' {
            if let Some(id) = parse_contig_id(line) {
                let rank = contig_rank.len();
                contig_rank.entry(id.to_vec()).or_insert(rank);
            }
            header_lines.push(line.to_vec());
            continue;
        }

        // Minimal parse: only CHROM, POS, REF, ALT needed for sort key; full raw line is re-emitted.
        let t1 = nth_tab(line, 1)
            .ok_or_else(|| RsomicsError::InvalidInput("VCF record missing POS column".into()))?;
        let t2 = nth_tab(line, 2)
            .ok_or_else(|| RsomicsError::InvalidInput("VCF record missing ID column".into()))?;
        let t3 = nth_tab(line, 3)
            .ok_or_else(|| RsomicsError::InvalidInput("VCF record missing REF column".into()))?;
        let t4 = nth_tab(line, 4)
            .ok_or_else(|| RsomicsError::InvalidInput("VCF record missing ALT column".into()))?;
        let t5 = nth_tab(line, 5).unwrap_or(line.len());

        let chrom = &line[..t1];
        let pos_bytes = &line[t1 + 1..t2];
        let ref_allele = line[t3 + 1..t4].to_vec();
        let alt_allele = line[t4 + 1..t5].to_vec();

        let pos = std::str::from_utf8(pos_bytes)
            .ok()
            .and_then(|s| s.parse::<i64>().ok())
            .ok_or_else(|| {
                RsomicsError::InvalidInput(format!(
                    "invalid POS: {}",
                    String::from_utf8_lossy(pos_bytes)
                ))
            })?;

        let rank = if let Some(&r) = contig_rank.get(chrom) {
            r
        } else {
            // Contig not in header: assign a stable rank beyond all declared contigs,
            // preserving first-appearance order among unknown contigs.
            let next = unknown_contigs.len();
            *unknown_contigs.entry(chrom.to_vec()).or_insert(next)
        };
        records.push(Record {
            contig_rank: rank,
            pos,
            ref_allele,
            alt_allele,
            line: line.to_vec(),
        });
    }

    // Unknown-contig ranks were assigned as 0-based appearance indices during parsing;
    // shift them past all declared contigs so they sort last.
    let declared_count = contig_rank.len();
    if !unknown_contigs.is_empty() {
        for rec in &mut records {
            let t1 = nth_tab(&rec.line, 1).unwrap();
            let chrom = &rec.line[..t1];
            if !contig_rank.contains_key(chrom)
                && let Some(&unk_idx) = unknown_contigs.get(chrom)
            {
                rec.contig_rank = declared_count + unk_idx;
            }
        }
    }

    // Stable sort: (contig_rank, pos, ref, alt)
    records.sort_by(|a, b| {
        a.contig_rank
            .cmp(&b.contig_rank)
            .then_with(|| a.pos.cmp(&b.pos))
            .then_with(|| a.ref_allele.cmp(&b.ref_allele))
            .then_with(|| a.alt_allele.cmp(&b.alt_allele))
    });

    let mut out = BufWriter::new(output);
    for h in &header_lines {
        out.write_all(h).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    let n = records.len() as u64;
    for rec in &records {
        out.write_all(&rec.line).map_err(RsomicsError::Io)?;
        out.write_all(b"\n").map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;

    Ok(n)
}
