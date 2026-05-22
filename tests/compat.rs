use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-vcf-sort"))
}

fn fixture() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden/unsorted.vcf")
}

fn bcftools_path() -> Option<String> {
    // bcftools may live in a conda env rather than $PATH on this machine
    let candidates = [
        "bcftools",
        "/opt/homebrew/Caskroom/miniforge/base/envs/imotif-pipeline/bin/bcftools",
        "/usr/bin/bcftools",
        "/usr/local/bin/bcftools",
    ];
    for candidate in &candidates {
        let ok = Command::new(candidate)
            .arg("--version")
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        if ok {
            return Some(candidate.to_string());
        }
    }
    None
}

/// Data (non-header) records only — bcftools adds its own header lines.
fn records(vcf: &[u8]) -> Vec<String> {
    String::from_utf8_lossy(vcf)
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// Sort order must match `bcftools sort` on the golden fixture.
///
/// bcftools sort orders by:
///   1. Contig rank (##contig header order)
///   2. POS ascending
///   3. REF lexicographic (tie-break)
///   4. ALT lexicographic (tie-break)
///
/// bcftools also injects a `##FILTER=<ID=PASS,...>` header line that we do not
/// emit (intentional — we keep the header as-is). Comparison is on data records only.
#[test]
fn sort_matches_bcftools() {
    let Some(bcftools) = bcftools_path() else {
        eprintln!("skipping: bcftools not found");
        return;
    };

    let version = Command::new(&bcftools)
        .arg("--version")
        .output()
        .unwrap()
        .stdout;
    eprintln!(
        "bcftools: {}",
        String::from_utf8_lossy(&version)
            .lines()
            .next()
            .unwrap_or("")
    );

    let vcf = fixture();

    let ours_out = ours().arg(&vcf).output().unwrap();
    assert!(
        ours_out.status.success(),
        "rsomics-vcf-sort failed: {}",
        String::from_utf8_lossy(&ours_out.stderr)
    );

    let theirs = Command::new(&bcftools)
        .args(["sort"])
        .arg(&vcf)
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .output()
        .unwrap();
    assert!(theirs.status.success(), "bcftools sort failed");

    let ours_records = records(&ours_out.stdout);
    let their_records = records(&theirs.stdout);

    assert_eq!(
        ours_records, their_records,
        "record order differs from bcftools sort"
    );
}
