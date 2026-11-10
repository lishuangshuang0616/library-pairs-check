mod barcode;
mod compare;
mod fastq;
mod input;
mod profile;

use anyhow::{bail, Result};
use barcode::Layout;
use clap::Parser;
use compare::Pair;
use input::{discover, library_name, Sample, CDNA_ROLE, OLIGO_ROLE};
use profile::{
    auto_threads, build_all, FileTiming, LibrarySpec, Profile, Settings, MIN_DISTINCT_BARCODES,
};
use serde::Serialize;
use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Parser)]
#[command(
    version,
    about = "Check DNBelab C4 cDNA/oligo pairing across sample directories (R1 only). Rankings do not prove sample identity."
)]
struct Args {
    /// Root of sample directories, or one sample directory; repeat for multiple roots.
    /// Each sample directory follows dnbc4tools_py: cDNA/ and oligo/ subfolders of FASTQ files.
    #[arg(long = "input", value_name = "DIR", required = true, action = clap::ArgAction::Append)]
    inputs: Vec<PathBuf>,
    #[arg(long, value_enum, default_value = "auto")]
    cdna_layout: Layout,
    #[arg(long, value_enum, default_value = "auto")]
    oligo_layout: Layout,
    /// Analyse all R1 reads (default), or only the first N reads of each FASTQ file; 0 = all.
    #[arg(long, value_name = "N", default_value_t = 0)]
    reads: usize,
    /// FASTQ files scanned concurrently; default = min(CPU cores, files, 8).
    /// Higher values multiply memory used by per-file count tables.
    #[arg(long, value_name = "N")]
    threads: Option<usize>,
    /// Minimum reads supporting a full barcode within each library.
    #[arg(long, default_value_t=3, value_parser=clap::value_parser!(u64).range(1..))]
    min_support: u64,
    /// JSON output; an existing file is overwritten atomically.
    #[arg(long)]
    output: PathBuf,
}

/// Wall-clock breakdown. Values are diagnostics only and are excluded from
/// cross-run report comparisons.
#[derive(Serialize)]
struct Timings {
    total_seconds: f64,
    discover_seconds: f64,
    scan_wall_seconds: f64,
    merge_seconds: f64,
    compare_seconds: f64,
    rank_seconds: f64,
    threads: usize,
    scanned_records: u64,
    records_per_second: f64,
    files: Vec<FileTiming>,
}

#[derive(Serialize)]
struct Report {
    schema_version: u32,
    tool_version: &'static str,
    method: &'static str,
    whitelist_segments: usize,
    settings: Settings,
    limitations: Vec<&'static str>,
    samples: Vec<Sample>,
    cdna: Vec<Profile>,
    oligo: Vec<Profile>,
    pairs: Vec<Pair>,
    rankings: Vec<compare::Ranking>,
    pairings: Vec<Pairing>,
    timings: Timings,
}

/// One row per sample: the directory-implied pair versus the barcode-evidence ranking.
#[derive(Serialize)]
struct Pairing {
    sample: String,
    cdna: String,
    expected_oligo: String,
    preferred_oligo: Option<String>,
    status: &'static str,
    assessment: &'static str,
    overlap_agrees: Option<bool>,
    coverage_agrees: Option<bool>,
    diagnostic_warnings: Vec<&'static str>,
    expected_pair_shared_barcodes: Option<usize>,
    best_cross_sample_oligo: Option<String>,
    best_cross_sample_shared_barcodes: Option<usize>,
}

fn run(args: Args) -> Result<()> {
    let run_started = Instant::now();
    let discover_started = Instant::now();
    let samples = discover(&args.inputs)?;
    let discover_seconds = discover_started.elapsed().as_secs_f64();
    let input_files: HashSet<PathBuf> = samples
        .iter()
        .flat_map(|s| {
            s.cdna
                .r1
                .iter()
                .chain(&s.cdna.r2)
                .chain(&s.oligo.r1)
                .chain(&s.oligo.r2)
        })
        .cloned()
        .collect();
    for sample in &samples {
        eprintln!(
            "Sample {}: cDNA R1={} (R2 ignored={}), oligo R1={} (R2 ignored={}) [{}]",
            sample.name,
            sample.cdna.r1.len(),
            sample.cdna.r2.len(),
            sample.oligo.r1.len(),
            sample.oligo.r2.len(),
            sample.directory.display()
        );
    }
    if args.output.exists() && input_files.contains(&args.output.canonicalize()?) {
        bail!("output must not replace an input FASTQ");
    }
    let parent = args
        .output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    if !parent.is_dir() {
        bail!(
            "output parent directory does not exist: {}",
            parent.display()
        );
    }
    let settings = Settings {
        reads_per_file: args.reads,
        min_support: args.min_support,
    };
    if args.reads > 0 {
        eprintln!(
            "Reading only the first {} reads per FASTQ file (order-sensitive; unread tails are not validated).",
            args.reads
        );
    }
    let white = barcode::whitelist()?;
    let mut libraries = Vec::with_capacity(samples.len() * 2);
    for sample in &samples {
        libraries.push(LibrarySpec {
            name: library_name(&sample.name, CDNA_ROLE),
            sample: sample.name.clone(),
            role: CDNA_ROLE,
            paths: sample.cdna.r1.clone(),
            layout: args.cdna_layout,
        });
        libraries.push(LibrarySpec {
            name: library_name(&sample.name, OLIGO_ROLE),
            sample: sample.name.clone(),
            role: OLIGO_ROLE,
            paths: sample.oligo.r1.clone(),
            layout: args.oligo_layout,
        });
    }
    let file_count: usize = libraries.iter().map(|library| library.paths.len()).sum();
    let requested = args.threads.unwrap_or_else(|| auto_threads(file_count));
    if requested == 0 {
        bail!("--threads must be at least 1");
    }
    let threads = requested.min(file_count.max(1));
    eprintln!("Scanning {file_count} FASTQ file(s) with {threads} worker thread(s).");
    let (profiles, build) = build_all(libraries, settings, &white, threads)?;
    let (mut cdna, mut oligo) = (Vec::new(), Vec::new());
    for profile in profiles {
        if profile.role == CDNA_ROLE {
            cdna.push(profile);
        } else {
            oligo.push(profile);
        }
    }
    let compare_started = Instant::now();
    let pairs: Vec<Pair> = cdna
        .iter()
        .flat_map(|c| {
            oligo
                .iter()
                .map(move |o| compare::compare(c, o, MIN_DISTINCT_BARCODES))
        })
        .collect();
    let compare_seconds = compare_started.elapsed().as_secs_f64();
    let rank_started = Instant::now();
    // Pairs are row-major: rank only the current cDNA row, not the full matrix.
    let row = |index: usize| &pairs[index * oligo.len()..(index + 1) * oligo.len()];
    let rankings: Vec<compare::Ranking> = cdna
        .iter()
        .enumerate()
        .map(|(index, c)| compare::rank(&c.name, &c.sample, row(index)))
        .collect();
    let pairings: Vec<Pairing> = cdna
        .iter()
        .zip(&rankings)
        .enumerate()
        .map(|(index, (c, ranking))| {
            let candidates = row(index);
            let expected = library_name(&c.sample, OLIGO_ROLE);
            let expected_pair = candidates
                .iter()
                .find(|p| p.cdna == c.name && p.oligo == expected);
            let cross = candidates
                .iter()
                .filter(|p| {
                    p.cdna == c.name
                        && !p.expected_pair
                        && p.sufficient_barcode_count
                        && p.shared_barcodes > 0
                })
                .max_by(|a, b| {
                    a.shared_barcodes
                        .cmp(&b.shared_barcodes)
                        .then(b.oligo.cmp(&a.oligo))
                });
            Pairing {
                sample: c.sample.clone(),
                cdna: c.name.clone(),
                expected_oligo: expected,
                preferred_oligo: ranking.preferred_candidate.clone(),
                status: ranking.status,
                assessment: ranking.expected_pair_assessment,
                overlap_agrees: ranking.overlap_agrees,
                coverage_agrees: ranking.coverage_agrees,
                diagnostic_warnings: ranking.diagnostic_warnings.clone(),
                expected_pair_shared_barcodes: expected_pair.map(|p| p.shared_barcodes),
                best_cross_sample_oligo: cross.map(|p| p.oligo.clone()),
                best_cross_sample_shared_barcodes: cross.map(|p| p.shared_barcodes),
            }
        })
        .collect();
    let rank_seconds = rank_started.elapsed().as_secs_f64();
    let scanned_records: u64 = cdna
        .iter()
        .chain(&oligo)
        .flat_map(|library| &library.files)
        .map(|file| file.scanned_records)
        .sum();
    let timings = Timings {
        total_seconds: run_started.elapsed().as_secs_f64(),
        discover_seconds,
        scan_wall_seconds: build.scan_wall_seconds,
        merge_seconds: build.merge_seconds,
        compare_seconds,
        rank_seconds,
        threads,
        scanned_records,
        records_per_second: if build.scan_wall_seconds > 0.0 {
            scanned_records as f64 / build.scan_wall_seconds
        } else {
            0.0
        },
        files: build.files,
    };
    let report = Report {
        schema_version: 7,
        tool_version: env!("CARGO_PKG_VERSION"),
        method: "exact_whitelisted_CB1_CB2_v7_abundance_decision_overlap_diagnostics",
        whitelist_segments: white.len(),
        settings,
        limitations: vec![
            "At least 100 shared supported barcodes are required before selecting a candidate; this quantity safeguard is not a calibrated match threshold.",
            "scan_wall_seconds includes merging and profile assembly; merge_seconds is an overlapping subset, not an additional stage.",
            "Raw and log1p-count cosine are correlated evidence, not independent validation; log cosine is not invariant to sequencing depth.",
            "Abundance tiers include all boundary ties; a flat count distribution has no high-abundance discrimination.",
            "Fully tied metrics do not distinguish candidates; ranking_basis and informative_metrics describe which metrics contribute.",
            "Raw and log-count cosine jointly select the relative candidate; Jaccard and read coverage are diagnostics only. No calibrated same-sample decision or p-value.",
            "Only technical library compatibility can be assessed; donor identity is not established.",
            "Sample directory names only label the expected cDNA/oligo pair; they are not biological evidence.",
            "Whitelist barcodes are reused across experiments; high overlap alone is not proof.",
            "Coverage uses supported reads, not deduplicated molecules; PCR and depth can bias it.",
            "By default every R1 read is counted; --reads N stops after the first N reads per file, is order-sensitive, and leaves unread tails unvalidated.",
            "C4 standard/dark R1 only; v1HT non-dark, custom structures and C6 are unsupported.",
            "Fixed Q20 minimum at each barcode base, exact whitelist matching, top-10000 barcodes and a 100-barcode ranking floor are built in, not configurable.",
            "Exact whitelist matches only; no mismatch correction and no oligo-derived correction of cDNA.",
            "Only R1 FASTQ files are analysed; R2 files are listed but not validated as pairs.",
        ],
        samples,
        cdna,
        oligo,
        pairs,
        rankings,
        pairings,
        timings,
    };
    let mut output = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut output, &report)?;
    writeln!(output)?;
    output.as_file().sync_all()?;
    output.persist(&args.output)?;
    println!("cDNA\toligo\texpected_pair\tshared_barcodes\tcDNA_read_coverage\tbarcode_jaccard\tcount_cosine\tlog_count_cosine\tenough_barcodes");
    for pair in &report.pairs {
        let display = |v: Option<f64>| v.map(|n| format!("{n:.6}")).unwrap_or_else(|| "NA".into());
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pair.cdna,
            pair.oligo,
            pair.expected_pair,
            pair.shared_barcodes,
            display(pair.cdna_read_coverage),
            display(pair.barcode_jaccard),
            display(pair.count_cosine),
            display(pair.log_count_cosine),
            pair.sufficient_barcode_count,
        );
    }
    for pairing in &report.pairings {
        eprintln!(
            "Pairing {}: expected={}, preferred={}, status={}, assessment={}, overlap_agrees={}, coverage_agrees={}, warnings={}",
            pairing.sample,
            pairing.expected_oligo,
            pairing.preferred_oligo.as_deref().unwrap_or("none"),
            pairing.status,
            pairing.assessment,
            pairing.overlap_agrees.map_or("NA", |v| if v { "true" } else { "false" }),
            pairing.coverage_agrees.map_or("NA", |v| if v { "true" } else { "false" }),
            if pairing.diagnostic_warnings.is_empty() { "none".to_owned() } else { pairing.diagnostic_warnings.join(",") }
        );
    }
    eprintln!(
        "Timings: discover={:.3}s scan={:.3}s merge={:.3}s compare={:.3}s rank={:.3}s total={:.3}s ({:.2} M records/s, {} threads)",
        report.timings.discover_seconds,
        report.timings.scan_wall_seconds,
        report.timings.merge_seconds,
        report.timings.compare_seconds,
        report.timings.rank_seconds,
        report.timings.total_seconds,
        report.timings.records_per_second / 1e6,
        report.timings.threads,
    );
    let mut slowest: Vec<&FileTiming> = report.timings.files.iter().collect();
    slowest.sort_by(|a, b| b.scan_seconds.total_cmp(&a.scan_seconds));
    for timing in slowest.iter().take(5) {
        eprintln!(
            "  slow file {:.3}s (detect {:.3}s, read {:.3}s) {} records: {}",
            timing.scan_seconds,
            timing.detect_seconds,
            timing.read_seconds,
            timing.scanned_records,
            timing.path.display()
        );
    }
    eprintln!(
        "Saved {}. Results are uncalibrated rankings, not sample-identity calls.",
        args.output.display()
    );
    Ok(())
}

fn main() {
    if let Err(error) = run(Args::parse()) {
        eprintln!("Error: {error:#}");
        std::process::exit(1);
    }
}
