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

fn report(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}

/// Compare reports while ignoring the diagnostic timings field.
fn without_timings(value: &Value) -> Value {
    let mut copy = value.clone();
    copy.as_object_mut().unwrap().remove("timings");
    copy
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
fn all_reads_by_default_confirm_expected_and_swapped_pairings() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    // Three samples: A/B oligo files are swapped, C is paired correctly.
    write_library(
        &root.join("SampleA/cDNA"),
        vec![
            ("S1_L001_R1_001.fq", sequences(0, 60, false), false, 'I'),
            ("S1_L002_R1_001.fq.gz", sequences(60, 60, true), true, 'I'),
            ("S1_L001_R2_001.fq", sequences(0, 1, false), false, 'I'),
        ],
    );
    write_library(
        &root.join("SampleA/oligo"),
        vec![("S1_L001_R1_001.fq.gz", sequences(200, 120, true), true, 'I')],
    );
    write_library(
        &root.join("SampleB/cDNA"),
        vec![("S2_L001_R1_001.fq", sequences(200, 120, false), false, 'I')],
    );
    write_library(
        &root.join("SampleB/oligo"),
        vec![("S2_L001_R1_001.fq.gz", sequences(0, 120, true), true, 'I')],
    );
    write_library(
        &root.join("SampleC/cDNA"),
        vec![("S3_L001_R1_001.fq", sequences(400, 120, false), false, 'I')],
    );
    write_library(
        &root.join("SampleC/oligo"),
        vec![("S3_L001_R1_001.fq", sequences(400, 120, true), false, 'I')],
    );
    let output = root.join("report.json");
    let args = base_args(root, &output);
    let result = run(&args);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report = report(&output);
    assert_eq!(report["schema_version"], 7);
    assert_eq!(report["settings"]["reads_per_file"], 0);
    assert_eq!(report["settings"]["min_support"], 3);
    let cdna_a = report["cdna"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "SampleA/cDNA")
        .unwrap();
    assert_eq!(cdna_a["supported_barcodes"], 120);
    assert_eq!(cdna_a["files"][0]["layout"], "standard");
    assert_eq!(cdna_a["files"][1]["layout"], "dark");
    assert_eq!(cdna_a["files"][0]["input_fully_scanned"], true);
    assert_eq!(cdna_a["files"][0]["record_limit_reached"], false);
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
        pairing(&report, "SampleA")["status"],
        "abundance_agree_unvalidated"
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
        pairing(&report, "SampleC")["assessment"],
        "expected_pair_supported"
    );
    // Existing output is overwritten by default, with the same content.
    let second = run(&args);
    assert!(
        second.status.success(),
        "{}",
        String::from_utf8_lossy(&second.stderr)
    );
    let overwritten: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(without_timings(&report), without_timings(&overwritten));
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
    let report = report(&output);
    assert_eq!(report["cdna"][0]["files"][0]["low_quality_barcode"], 72);
    assert!(report["pairs"][0]["cdna_read_coverage"].is_null());
    assert_eq!(report["rankings"][0]["status"], "insufficient_data");
    assert_eq!(report["pairings"][0]["assessment"], "unresolved");
}

#[test]
fn duplicate_inputs_and_truncated_gzip_fail_without_replacing_the_report() {
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
        args.extend(["--output".into(), output.display().to_string()]);
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
fn reads_limit_stops_early_and_all_reads_is_reproducible() {
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
    let mut limited = base_args(&root, &output);
    limited.extend(["--reads".into(), "300".into()]);
    assert!(run(&limited).status.success());
    let limited_report = report(&output);
    assert_eq!(
        limited_report["cdna"][0]["files"][0]["scanned_records"],
        300
    );
    assert_eq!(
        limited_report["cdna"][0]["files"][0]["input_fully_scanned"],
        false
    );
    assert_eq!(
        limited_report["cdna"][0]["files"][0]["record_limit_reached"],
        true
    );
    assert_eq!(limited_report["pairs"][0]["cdna_read_coverage"], 0.0);

    let all = base_args(&root, &output);
    assert!(run(&all).status.success());
    let first = fs::read(&output).unwrap();
    let all_report: Value = serde_json::from_slice(&first).unwrap();
    assert!(
        all_report["pairs"][0]["cdna_read_coverage"]
            .as_f64()
            .unwrap()
            > 0.6
    );
    assert_eq!(all_report["cdna"][0]["files"][0]["scanned_records"], 600);
    assert_eq!(
        all_report["cdna"][0]["files"][0]["input_fully_scanned"],
        true
    );
    assert!(run(&all).status.success());
    let second: Value = serde_json::from_slice(&fs::read(&output).unwrap()).unwrap();
    assert_eq!(without_timings(&all_report), without_timings(&second));
    let mut explicit_zero = base_args(&root, &output);
    explicit_zero.extend(["--reads".into(), "0".into()]);
    assert!(run(&explicit_zero).status.success());
    assert_eq!(
        without_timings(&all_report),
        without_timings(&report(&output)),
        "--reads 0 must match the default all-reads run"
    );

    // Adding another sample must not change an existing library's counts.
    write_library(
        &root.join("SampleB/cDNA"),
        vec![("S2_L001_R1_001.fq", sequences(0, 200, false), false, 'I')],
    );
    write_library(
        &root.join("SampleB/oligo"),
        vec![("S2_L001_R1_001.fq", reordered, false, 'I')],
    );
    assert!(run(&all).status.success());
    let expanded = report(&output);
    let same = expanded["cdna"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "SampleA/cDNA")
        .unwrap();
    assert_eq!(same, &all_report["cdna"][0]);
}

#[test]
fn tail_records_are_checked_by_default_but_not_with_a_reads_limit() {
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
    ]);
    assert!(!run(&base).status.success());
    assert!(!output.exists());
    let mut limited = base;
    limited.extend(["--reads".into(), "10".into()]);
    assert!(run(&limited).status.success());
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
    let result = run(&base_args(&root, &output));
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report = report(&output);
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

    // A single sample directory can be passed directly; a lone unlabelled
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
fn removed_parameters_are_rejected() {
    let dir = tempfile::tempdir().unwrap();
    for flag in [
        "--force",
        "--sampling",
        "--min-base-quality",
        "--top-barcodes",
        "--min-distinct-barcodes",
        "--max-records-per-file",
        "--seed",
    ] {
        let args = vec![
            flag.into(),
            "1".into(),
            "--output".into(),
            dir.path().join("unused.json").display().to_string(),
        ];
        assert!(!run(&args).status.success(), "{flag} should be rejected");
    }
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
    let report = report(&output);
    assert_eq!(
        report["samples"][0]["cdna"]["r1"].as_array().unwrap().len(),
        1
    );
    assert_eq!(
        report["samples"][0]["cdna"]["r2"].as_array().unwrap().len(),
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
    for (case, extra) in [("all", vec![]), ("limited", vec!["--reads", "50"])] {
        let serial = dir.path().join(format!("{case}-serial.json"));
        let parallel = dir.path().join(format!("{case}-parallel.json"));
        let mut serial_args = base_args(&root, &serial);
        serial_args.extend(["--threads".into(), "1".into()]);
        serial_args.extend(extra.iter().map(|value| (*value).into()));
        let mut parallel_args = base_args(&root, &parallel);
        parallel_args.extend(["--threads".into(), "4".into()]);
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
        let serial_value = report(&serial);
        let parallel_value = report(&parallel);
        assert_eq!(
            serial_value["timings"]["files"].as_array().unwrap().len(),
            9
        );
        assert!(
            serial_value["timings"]["scan_wall_seconds"]
                .as_f64()
                .unwrap()
                > 0.0
        );
        assert_eq!(
            without_timings(&serial_value),
            without_timings(&parallel_value),
            "{case}: parallel report differs from serial"
        );
        if case == "limited" {
            assert_eq!(
                serial_value["cdna"][0]["files"][0]["input_fully_scanned"],
                false
            );
        }
        let default_report = dir.path().join(format!("{case}-default.json"));
        let mut default_args = base_args(&root, &default_report);
        default_args.extend(extra.iter().map(|value| (*value).into()));
        let default_result = run(&default_args);
        assert!(
            default_result.status.success(),
            "{}",
            String::from_utf8_lossy(&default_result.stderr)
        );
        assert_eq!(
            without_timings(&serial_value),
            without_timings(&report(&default_report)),
            "{case}: default-thread report differs from serial"
        );
    }
    let mut zero_args = base_args(&root, &dir.path().join("zero.json"));
    zero_args.extend(["--threads".into(), "0".into()]);
    let zero = run(&zero_args);
    assert!(!zero.status.success());
    assert!(String::from_utf8_lossy(&zero.stderr).contains("--threads must be at least 1"));
}

#[test]
fn weak_overlap_stays_unresolved_and_single_pair_is_never_confirmed() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    let mut weak = sequences(0, 1, true);
    weak.extend(sequences(200, 119, true));
    for (sample, cdna, oligo) in [
        ("A", sequences(0, 120, false), weak),
        ("B", sequences(400, 120, false), sequences(600, 120, true)),
    ] {
        write_library(
            &root.join(sample).join("cDNA"),
            vec![("lane_R1.fq", cdna, false, 'I')],
        );
        write_library(
            &root.join(sample).join("oligo"),
            vec![("lane_R1.fq", oligo, false, 'I')],
        );
    }
    let output = root.join("report.json");
    let result = run(&base_args(root, &output));
    assert!(result.status.success());
    let value = report(&output);
    assert_eq!(
        value["rankings"][0]["status"],
        "insufficient_shared_barcodes"
    );
    assert_eq!(value["pairs"][0]["shared_barcodes"], 1);
    assert_eq!(value["pairs"][0]["sufficient_shared_barcode_count"], false);
    assert_eq!(pairing(&value, "A")["assessment"], "unresolved");
    assert!(value["rankings"][0]["preferred_candidate"].is_null());
    // Even a perfect full-depth overlap cannot establish identity without alternatives.
    fastq(
        &root.join("A/oligo/lane_R1.fq"),
        &sequences(0, 120, true),
        false,
        'I',
    );
    let result = run(&base_args(&root.join("A"), &output));
    assert!(result.status.success());
    let value = report(&output);
    assert_eq!(value["pairs"][0]["barcode_jaccard"], 1.0);
    assert_eq!(
        value["rankings"][0]["status"],
        "single_candidate_unresolved"
    );
    assert!(value["rankings"][0]["preferred_candidate"].is_null());
}

#[test]
fn saturated_sets_expose_abundance_evidence_and_tie_aware_tops() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    for sample in ["A", "B"] {
        for (role, dark) in [("cDNA", false), ("oligo", true)] {
            let mut seqs = sequences(0, 240, dark);
            let abundant = if sample == "A" { 0 } else { 120 };
            // fastq() writes each item 3 times: high=30, low=3 reads.
            for _ in 0..9 {
                seqs.extend(sequences(abundant, 120, dark));
            }
            write_library(
                &root.join(sample).join(role),
                vec![("lane_R1.fq.gz", seqs, true, 'I')],
            );
        }
    }
    let serial_path = root.join("serial.json");
    let parallel_path = root.join("parallel.json");
    let mut args = base_args(root, &serial_path);
    args.extend(["--threads".into(), "1".into()]);
    assert!(run(&args).status.success());
    let mut args = base_args(root, &parallel_path);
    args.extend(["--threads".into(), "4".into()]);
    assert!(run(&args).status.success());
    let value = report(&serial_path);
    assert_eq!(
        without_timings(&value),
        without_timings(&report(&parallel_path))
    );
    let own = &value["pairs"][0];
    let wrong = &value["pairs"][1];
    assert_eq!(own["barcode_jaccard"], 1.0);
    assert_eq!(wrong["barcode_jaccard"], 1.0);
    assert_eq!(own["log_count_cosine"], 1.0);
    assert!(own["count_cosine"].as_f64().unwrap() > wrong["count_cosine"].as_f64().unwrap());
    assert_eq!(own["abundance_overlaps"][0]["cdna_barcodes"], 120);
    assert_eq!(own["abundance_overlaps"][0]["overlap"], 1.0);
    assert_eq!(wrong["abundance_overlaps"][0]["overlap"], 0.0);
    assert_eq!(wrong["abundance_overlaps"][1]["overlap"], 1.0);
    assert_eq!(value["rankings"][0]["ranking_basis"], "abundance_only");
    assert_eq!(
        pairing(&value, "A")["assessment"],
        "expected_pair_supported"
    );
    assert!(
        value["cdna"][0]["count_effective_barcodes"]
            .as_f64()
            .unwrap()
            < 240.0
    );
    let single_path = root.join("single.json");
    let result = run(&base_args(&root.join("A"), &single_path));
    assert!(result.status.success());
    let single = report(&single_path);
    assert_eq!(single["rankings"][0]["ranking_basis"], "descriptive_only");
    assert_eq!(
        single["rankings"][0]["status"],
        "single_candidate_unresolved"
    );
    assert_eq!(single["pairs"][0]["log_count_cosine"], 1.0);
}

#[test]
fn coverage_disagreement_is_reported_without_vetoing_abundance() {
    let dir = tempfile::tempdir().unwrap();
    let root = dir.path();
    for (sample, cdna, oligo) in [
        ("A", sequences(0, 120, false), sequences(0, 100, true)),
        ("B", sequences(1000, 120, false), sequences(0, 1000, true)),
    ] {
        write_library(
            &root.join(sample).join("cDNA"),
            vec![("lane_R1.fq", cdna, false, 'I')],
        );
        write_library(
            &root.join(sample).join("oligo"),
            vec![("lane_R1.fq.gz", oligo, true, 'I')],
        );
    }
    let output = root.join("report.json");
    let result = run(&base_args(root, &output));
    assert!(result.status.success());
    let value = report(&output);
    let ranked = &value["rankings"][0];
    assert_eq!(ranked["status"], "abundance_agree_unvalidated");
    assert_eq!(ranked["preferred_candidate"], "A/oligo");
    assert_eq!(ranked["abundance_candidate"], "A/oligo");
    assert_eq!(ranked["coverage_agrees"], false);
    assert_eq!(ranked["overlap_agrees"], true);
    assert_eq!(
        ranked["diagnostic_warnings"][0],
        "read_coverage_disagrees_with_abundance"
    );
    assert_eq!(pairing(&value, "A")["coverage_agrees"], false);
    assert!(String::from_utf8_lossy(&result.stderr).contains("coverage_agrees=false"));
}
