use flate2::{write::GzEncoder, Compression};
use serde_json::Value;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Output};

fn sequences(start: usize, count: usize, dark: bool) -> Vec<String> {
    let white: Vec<_> = include_str!("../assets/barcode_whitelist.txt")
        .lines()
        .collect();
    (start..start + count)
        .map(|i| {
            if dark {
                format!("{}{}", white[i], white[200 + i])
            } else {
                format!("TCTGCG{}CCTTCC{}", white[i], white[200 + i])
            }
        })
        .collect()
}

fn fastq(path: &Path, sequences: &[String], gzip: bool, quality: char) {
    let mut body = String::new();
    for (i, seq) in sequences.iter().enumerate() {
        for duplicate in 0..3 {
            body.push_str(&format!(
                "@read{i}_{duplicate}\n{seq}\n+\n{}\n",
                quality.to_string().repeat(seq.len())
            ));
        }
    }
    let bytes = if gzip {
        let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
        encoder.write_all(body.as_bytes()).unwrap();
        encoder.finish().unwrap()
    } else {
        body.into_bytes()
    };
    fs::write(path, bytes).unwrap();
}

fn write_library(dir: &Path, files: Vec<(&str, Vec<String>, bool, char)>) {
    fs::create_dir_all(dir).unwrap();
    for (name, seqs, gzip, quality) in &files {
        fastq(&dir.join(name), seqs, *gzip, *quality);
    }
}

fn run(args: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_library-pair-check"))
        .args(args)
        .output()
        .unwrap()
}

fn base_args(root: &Path, output: &Path) -> Vec<String> {
    vec![
        "--input".into(),
        root.display().to_string(),
        "--output".into(),
        output.display().to_string(),
    ]
}

fn pairing<'a>(report: &'a Value, sample: &str) -> &'a Value {
    report["pairings"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["sample"] == sample)
        .unwrap()
}

#[test]
fn sample_directories_expose_swapped_and_correct_pairings() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let mut b_oligo = sequences(0, 12, true);
    b_oligo.extend(sequences(12, 12, true));
    write_library(
        &root.join("SampleA/cDNA"),
        vec![
            ("S1_L001_R1_001.fq", sequences(0, 12, false), false, 'I'),
            ("S1_L002_R1_001.fq.gz", sequences(12, 12, true), true, 'I'),
            ("S1_L001_R2_001.fq", sequences(0, 1, false), false, 'I'),
        ],
    );
    write_library(
        &root.join("SampleA/oligo"),
        vec![("S1_L001_R1_001.fq.gz", sequences(100, 24, true), true, 'I')],
    );
    write_library(
        &root.join("SampleB/cDNA"),
        vec![("S2_L001_R1_001.fq", sequences(100, 24, false), false, 'I')],
    );
    write_library(
        &root.join("SampleB/oligo"),
        vec![("S2_L001_R1_001.fq.gz", b_oligo, true, 'I')],
    );
    write_library(
        &root.join("SampleC/cDNA"),
        vec![("S3_L001_R1_001.fq", sequences(200, 24, false), false, 'I')],
    );
    write_library(
        &root.join("SampleC/oligo"),
        vec![("S3_L001_R1_001.fq", sequences(200, 24, true), false, 'I')],
    );
    let output = root.join("report.json");
    let args = vec![
        "--input".into(),
        root.display().to_string(),
        "--min-distinct-barcodes".into(),
        "10".into(),
        "--output".into(),
        output.display().to_string(),
    ];
    let result = run(&args);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(report["schema_version"], 3);
    assert_eq!(report["samples"].as_array().unwrap().len(), 3);
    let cdna_a = report["cdna"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "SampleA/cDNA")
        .unwrap();
    assert_eq!(cdna_a["supported_barcodes"], 24);
    assert_eq!(cdna_a["files"][0]["layout"], "standard");
    assert_eq!(cdna_a["files"][1]["layout"], "dark");
    assert_eq!(
        report["pairs"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|p| p["expected_pair"].as_bool().unwrap())
            .count(),
        3
    );
    assert_eq!(
        pairing(&report, "SampleA")["assessment"],
        "suspected_cross_sample_pairing"
    );
    assert_eq!(
        pairing(&report, "SampleA")["preferred_oligo"],
        "SampleB/oligo"
    );
    assert_eq!(
        pairing(&report, "SampleB")["assessment"],
        "suspected_cross_sample_pairing"
    );
    assert_eq!(
        pairing(&report, "SampleC")["assessment"],
        "expected_pair_supported"
    );
    assert_eq!(
        pairing(&report, "SampleC")["preferred_oligo"],
        "SampleC/oligo"
    );
    let swapped = report["pairs"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["cdna"] == "SampleA/cDNA" && p["oligo"] == "SampleB/oligo")
        .unwrap();
    assert_eq!(swapped["cdna_read_coverage"], 1.0);
    let original = fs::read(&output).unwrap();
    assert!(!run(&args).status.success());
    assert_eq!(fs::read(&output).unwrap(), original);
}

#[test]
fn low_quality_and_zero_denominators_do_not_produce_a_match() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    write_library(
        &root.join("SampleA/cDNA"),
        vec![("S1_L001_R1_001.fq", sequences(0, 24, true), false, '!')],
    );
    write_library(
        &root.join("SampleA/oligo"),
        vec![("S1_L001_R1_001.fq", sequences(0, 24, true), false, 'I')],
    );
    let output = dir.path().join("report.json");
    let mut args = base_args(root, &output);
    let failed = run(&args);
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("layout detection inconclusive"));
    assert!(!output.exists());
    args.extend([
        "--cdna-layout".into(),
        "dark".into(),
        "--oligo-layout".into(),
        "dark".into(),
    ]);
    assert!(run(&args).status.success());
    let report: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(report["cdna"][0]["files"][0]["low_quality_barcode"], 72);
    assert!(report["pairs"][0]["cdna_read_coverage"].is_null());
    assert_eq!(report["rankings"][0]["status"], "insufficient_data");
    assert_eq!(report["pairings"][0]["assessment"], "unresolved");
}

#[test]
fn duplicate_inputs_and_truncated_gzip_fail_without_overwriting_results() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("samples");
    write_library(
        &root.join("SampleA/cDNA"),
        vec![("S1_L001_R1_001.fq", sequences(0, 24, true), false, 'I')],
    );
    let oligo_dir = root.join("SampleA/oligo");
    write_library(
        &oligo_dir,
        vec![("S1_L001_R1_001.fq.gz", sequences(0, 24, true), true, 'I')],
    );
    let output = dir.path().join("report.json");
    fs::write(&output, "preserve me").unwrap();
    let make_args = |inputs: &[&Path]| -> Vec<String> {
        let mut args = Vec::new();
        for input in inputs {
            args.extend(["--input".into(), input.display().to_string()]);
        }
        args.extend([
            "--output".into(),
            output.display().to_string(),
            "--force".into(),
            "--max-records-per-file".into(),
            "0".into(),
        ]);
        args
    };
    let duplicate = run(&make_args(&[&root, &root]));
    assert!(!duplicate.status.success());
    assert!(String::from_utf8_lossy(&duplicate.stderr).contains("duplicate sample name"));
    let oligo = oligo_dir.join("S1_L001_R1_001.fq.gz");
    let mut gzip = fs::read(&oligo).unwrap();
    gzip.truncate(gzip.len() - 8);
    fs::write(&oligo, gzip).unwrap();
    assert!(!run(&make_args(&[&root])).status.success());
    assert_eq!(fs::read_to_string(&output).unwrap(), "preserve me");
}

#[test]
fn reservoir_reduces_prefix_order_bias_and_is_reproducible() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("samples");
    write_library(
        &root.join("SampleA/cDNA"),
        vec![("S1_L001_R1_001.fq", sequences(0, 200, false), false, 'I')],
    );
    let mut reordered = sequences(100, 100, true);
    reordered.extend(sequences(0, 100, true));
    write_library(
        &root.join("SampleA/oligo"),
        vec![("S1_L001_R1_001.fq", reordered.clone(), false, 'I')],
    );
    let output = dir.path().join("result.json");
    let mut base = base_args(&root, &output);
    base.extend([
        "--force".into(),
        "--max-records-per-file".into(),
        "300".into(),
        "--min-support".into(),
        "1".into(),
    ]);
    let mut prefix = base.clone();
    prefix.extend(["--sampling".into(), "prefix".into()]);
    assert!(run(&prefix).status.success());
    let p: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(p["pairs"][0]["cdna_read_coverage"], 0.0);
    assert_eq!(p["cdna"][0]["files"][0]["input_fully_scanned"], false);
    assert!(run(&base).status.success());
    let first = fs::read(&output).unwrap();
    let r: Value = serde_json::from_slice(&first).unwrap();
    assert!(r["pairs"][0]["cdna_read_coverage"].as_f64().unwrap() > 0.6);
    assert_eq!(r["cdna"][0]["files"][0]["scanned_records"], 600);
    assert_eq!(r["cdna"][0]["files"][0]["sampled_records"], 300);
    assert_eq!(r["cdna"][0]["files"][0]["input_fully_scanned"], true);
    assert_ne!(
        r["cdna"][0]["files"][0]["sampling_seed"],
        r["oligo"][0]["files"][0]["sampling_seed"]
    );
    assert!(run(&base).status.success());
    assert_eq!(fs::read(&output).unwrap(), first);
    // Adding another sample must not change an existing file's reservoir draw.
    write_library(
        &root.join("SampleB/cDNA"),
        vec![("S2_L001_R1_001.fq", sequences(0, 200, false), false, 'I')],
    );
    write_library(
        &root.join("SampleB/oligo"),
        vec![("S2_L001_R1_001.fq", reordered, false, 'I')],
    );
    assert!(run(&base).status.success());
    let expanded: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    let same = expanded["cdna"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "SampleA/cDNA")
        .unwrap();
    assert_eq!(same, &r["cdna"][0]);
}

#[test]
fn reservoir_checks_tail_after_sample_limit_but_explicit_prefix_does_not() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("samples");
    let cdna_dir = root.join("SampleA/cDNA");
    write_library(
        &cdna_dir,
        vec![("S1_L001_R1_001.fq", sequences(0, 20, true), false, 'I')],
    );
    write_library(
        &root.join("SampleA/oligo"),
        vec![("S1_L001_R1_001.fq", sequences(0, 20, true), false, 'I')],
    );
    let broken = cdna_dir.join("S1_L001_R1_001.fq");
    let mut body = fs::read(&broken).unwrap();
    body.extend(b"@broken\nACGT\n+\n");
    fs::write(&broken, body).unwrap();
    let output = dir.path().join("report.json");
    let mut base = base_args(&root, &output);
    base.extend([
        "--cdna-layout".into(),
        "dark".into(),
        "--oligo-layout".into(),
        "dark".into(),
        "--max-records-per-file".into(),
        "10".into(),
    ]);
    assert!(!run(&base).status.success());
    assert!(!output.exists());
    let mut prefix = base;
    prefix.extend(["--sampling".into(), "prefix".into()]);
    assert!(run(&prefix).status.success());
}

#[test]
fn sample_directory_discovery_follows_dnbc4tools_naming() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("samples");
    // Aliases from dnbc4tools_py, multiple lanes, and R2 files (not analysed).
    write_library(
        &root.join("SampleA/cdna"),
        vec![
            ("S1_L001_R1_001.fq", sequences(0, 24, false), false, 'I'),
            ("S1_L002_R1_001.fq.gz", sequences(24, 24, false), true, 'I'),
            ("S1_L001_R2_001.fq", sequences(0, 1, false), false, 'I'),
        ],
    );
    write_library(
        &root.join("SampleA/Oligo"),
        vec![
            ("S1_L001_R1_001.fq", sequences(0, 24, true), false, 'I'),
            ("S1_L001_R2_001.fq", sequences(0, 1, true), false, 'I'),
        ],
    );
    let output = dir.path().join("report.json");
    let mut args = base_args(&root, &output);
    args.extend(["--min-distinct-barcodes".into(), "10".into()]);
    let result = run(&args);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(report["samples"][0]["name"], "SampleA");
    assert_eq!(
        report["samples"][0]["cdna"]["r1"].as_array().unwrap().len(),
        2
    );
    assert_eq!(
        report["samples"][0]["cdna"]["r2"].as_array().unwrap().len(),
        1
    );
    assert_eq!(
        report["samples"][0]["oligo"]["r1"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert!(report["pairs"][0]["expected_pair"].as_bool().unwrap());
    assert_eq!(report["pairs"][0]["shared_barcodes"], 24);

    // A single sample directory can be passed directly; an unlabelled lone
    // FASTQ is treated as R1.
    let single = dir.path().join("SingleSample");
    write_library(
        &single.join("cDNA"),
        vec![("all.fq", sequences(0, 24, true), false, 'I')],
    );
    write_library(
        &single.join("oligo"),
        vec![("all.fq", sequences(0, 24, true), false, 'I')],
    );
    let single_output = dir.path().join("single.json");
    assert!(run(&base_args(&single, &single_output)).status.success());

    let only_r2 = dir.path().join("only_r2");
    write_library(
        &only_r2.join("SampleA/cDNA"),
        vec![("S1_L001_R1_001.fq", sequences(0, 24, true), false, 'I')],
    );
    write_library(
        &only_r2.join("SampleA/oligo"),
        vec![("S1_L001_R2_001.fq", sequences(0, 24, true), false, 'I')],
    );
    let missing_r1 = run(&base_args(&only_r2, &dir.path().join("missing.json")));
    assert!(!missing_r1.status.success());
    assert!(String::from_utf8_lossy(&missing_r1.stderr).contains("no R1 FASTQ found"));

    let empty_root = dir.path().join("empty_root");
    fs::create_dir(&empty_root).unwrap();
    let no_samples = run(&base_args(&empty_root, &dir.path().join("no_samples.json")));
    assert!(!no_samples.status.success());
    assert!(String::from_utf8_lossy(&no_samples.stderr).contains("no sample directories under"));

    let duplicate_name = dir.path().join("duplicate_name");
    write_library(
        &duplicate_name.join("root1/SampleA/cDNA"),
        vec![("S1_L001_R1_001.fq", sequences(0, 24, true), false, 'I')],
    );
    write_library(
        &duplicate_name.join("root1/SampleA/oligo"),
        vec![("S1_L001_R1_001.fq", sequences(0, 24, true), false, 'I')],
    );
    write_library(
        &duplicate_name.join("root2/SampleA/cDNA"),
        vec![("S2_L001_R1_001.fq", sequences(40, 24, true), false, 'I')],
    );
    write_library(
        &duplicate_name.join("root2/SampleA/oligo"),
        vec![("S2_L001_R1_001.fq", sequences(40, 24, true), false, 'I')],
    );
    let duplicate_args: Vec<String> = vec![
        "--input".into(),
        duplicate_name.join("root1").display().to_string(),
        "--input".into(),
        duplicate_name.join("root2").display().to_string(),
        "--output".into(),
        dir.path().join("duplicate.json").display().to_string(),
    ];
    let duplicate = run(&duplicate_args);
    assert!(!duplicate.status.success());
    assert!(String::from_utf8_lossy(&duplicate.stderr).contains("duplicate sample name"));
}

#[test]
fn input_model_rejects_legacy_flags() {
    let dir = tempfile::tempdir().unwrap();
    let legacy = vec![
        "--cdna".into(),
        "RNA=a.fq".into(),
        "--oligo".into(),
        "O=b.fq".into(),
        "--output".into(),
        dir.path().join("legacy.json").display().to_string(),
    ];
    assert!(!run(&legacy).status.success());
}

#[cfg(unix)]
#[test]
fn symlinked_fastq_are_discovered_by_link_name() {
    use std::os::unix::fs::symlink;
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("samples");
    let raw = dir.path().join("raw");
    let cdna = root.join("SampleA/cDNA");
    let oligo = root.join("SampleA/oligo");
    for path in [&raw, &cdna, &oligo] {
        fs::create_dir_all(path).unwrap();
    }
    let seqs = sequences(0, 24, true);
    // Real files have opaque names and are plain text, while links carry the
    // R1/R2 tokens and .fq.gz suffix used on the cluster.
    fastq(&raw.join("cdna_data.fq"), &seqs, false, 'I');
    fastq(&raw.join("cdna_part2.fq"), &seqs, false, 'I');
    fastq(&raw.join("oligo_data.fq"), &seqs, false, 'I');
    symlink(
        raw.join("cdna_data.fq"),
        cdna.join("E100062882_L01_9_1.fq.gz"),
    )
    .unwrap();
    symlink(
        raw.join("cdna_part2.fq"),
        cdna.join("E100062882_L01_9_2.fq.gz"),
    )
    .unwrap();
    symlink(
        raw.join("oligo_data.fq"),
        oligo.join("E100062882_L01_10_1.fq.gz"),
    )
    .unwrap();
    let output = dir.path().join("report.json");
    let result = run(&base_args(&root, &output));
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(
        report["samples"][0]["cdna"]["r1"].as_array().unwrap().len(),
        1
    );
    assert_eq!(
        report["samples"][0]["cdna"]["r2"].as_array().unwrap().len(),
        1
    );
    assert_eq!(
        report["samples"][0]["oligo"]["r1"]
            .as_array()
            .unwrap()
            .len(),
        1
    );
    assert!(report["samples"][0]["cdna"]["r1"][0]
        .as_str()
        .unwrap()
        .ends_with("raw/cdna_data.fq"));
    assert_eq!(report["pairs"][0]["shared_barcodes"], 24);
}

#[test]
fn parallel_scans_reproduce_the_serial_report() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path().join("samples");
    // Three samples with uneven, multi-lane inputs to exercise the work queue.
    for (index, start) in [0usize, 60, 120].into_iter().enumerate() {
        let name = format!("Sample{}", index + 1);
        write_library(
            &root.join(&name).join("cDNA"),
            vec![
                ("L001_R1_001.fq", sequences(start, 30, false), false, 'I'),
                (
                    "L002_R1_001.fq.gz",
                    sequences(start + 30, 30, true),
                    true,
                    'I',
                ),
            ],
        );
        let mut oligo = sequences(start, 30, true);
        oligo.extend(sequences(start + 30, 30, true));
        write_library(
            &root.join(&name).join("oligo"),
            vec![("L001_R1_001.fq", oligo, false, 'I')],
        );
    }
    for (case, extra) in [
        ("full", vec!["--max-records-per-file", "0"]),
        (
            "sampled",
            vec!["--max-records-per-file", "200", "--min-support", "1"],
        ),
    ] {
        let serial = dir.path().join(format!("{case}-serial.json"));
        let parallel = dir.path().join(format!("{case}-parallel.json"));
        let mut serial_args = base_args(&root, &serial);
        serial_args.extend([
            "--min-distinct-barcodes".into(),
            "10".into(),
            "--threads".into(),
            "1".into(),
        ]);
        serial_args.extend(extra.iter().map(|value| (*value).into()));
        let mut parallel_args = base_args(&root, &parallel);
        parallel_args.extend([
            "--min-distinct-barcodes".into(),
            "10".into(),
            "--threads".into(),
            "4".into(),
        ]);
        parallel_args.extend(extra.iter().map(|value| (*value).into()));
        let serial_result = run(&serial_args);
        assert!(
            serial_result.status.success(),
            "{}",
            String::from_utf8_lossy(&serial_result.stderr)
        );
        let parallel_result = run(&parallel_args);
        assert!(
            parallel_result.status.success(),
            "{}",
            String::from_utf8_lossy(&parallel_result.stderr)
        );
        assert!(
            String::from_utf8_lossy(&parallel_result.stderr).contains("with 4 worker thread(s)")
        );
        assert_eq!(
            fs::read(&serial).unwrap(),
            fs::read(&parallel).unwrap(),
            "{case}: parallel report differs from serial"
        );
        let default_report = dir.path().join(format!("{case}-default.json"));
        let mut default_args = base_args(&root, &default_report);
        default_args.extend(["--min-distinct-barcodes".into(), "10".into()]);
        default_args.extend(extra.iter().map(|value| (*value).into()));
        let default_result = run(&default_args);
        assert!(
            default_result.status.success(),
            "{}",
            String::from_utf8_lossy(&default_result.stderr)
        );
        assert_eq!(
            fs::read(&serial).unwrap(),
            fs::read(&default_report).unwrap(),
            "{case}: default-thread report differs from serial"
        );
    }
    let mut zero_args = base_args(&root, &dir.path().join("zero.json"));
    zero_args.extend(["--threads".into(), "0".into()]);
    let zero = run(&zero_args);
    assert!(!zero.status.success());
    assert!(String::from_utf8_lossy(&zero.stderr).contains("--threads must be at least 1"));
}
