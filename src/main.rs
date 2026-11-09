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
use profile::{auto_threads, build_all, LibrarySpec, Profile, Sampling, Settings};
use serde::Serialize;
use std::collections::HashSet;
use std::io::Write;
use std::path::{Path, PathBuf};

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
    /// Records retained per input file; 0 keeps all. Reservoir mode still scans the entire file.
    #[arg(long, default_value_t = 1_000_000)]
    max_records_per_file: usize,
    /// Reservoir scans the full stream; prefix is faster but order-sensitive.
    #[arg(long, value_enum, default_value = "reservoir")]
    sampling: Sampling,
    /// Reproducible base seed; independent streams are used for each library/file.
    #[arg(long, default_value_t = 42)]
    seed: u64,
    /// FASTQ files scanned concurrently; default = min(CPU cores, files, 8).
    /// Higher values multiply memory used by reservoirs and count tables.
    #[arg(long, value_name = "N")]
    threads: Option<usize>,
    /// Minimum Phred+33 quality at every barcode base.
    #[arg(long, default_value_t=20, value_parser=clap::value_parser!(u8).range(0..=93))]
    min_base_quality: u8,
    /// Minimum reads supporting a full barcode within each library.
    #[arg(long, default_value_t=3, value_parser=clap::value_parser!(u64).range(1..))]
    min_support: u64,
    #[arg(long, default_value_t = 10_000)]
    top_barcodes: usize,
    /// Eligibility floor for ranking only; not a sample-match threshold.
    #[arg(long, default_value_t = 100)]
    min_distinct_barcodes: usize,
    /// JSON output; an existing file is preserved unless --force is given.
    #[arg(long)]
    output: PathBuf,
    #[arg(long)]
    force: bool,
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
    expected_pair_shared_barcodes: Option<usize>,
    best_cross_sample_oligo: Option<String>,
    best_cross_sample_shared_barcodes: Option<usize>,
}

fn run(args: Args) -> Result<()> {
    if args.top_barcodes == 0 || args.min_distinct_barcodes == 0 {
        bail!("--top-barcodes and --min-distinct-barcodes must be positive");
    }
    let samples = discover(&args.inputs)?;
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
    if args.output.exists() && !args.force {
        bail!(
            "output already exists: {}; use --force to replace it",
            args.output.display()
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
        max_records_per_file: args.max_records_per_file,
        sampling: args.sampling,
        seed: args.seed,
        min_base_quality: args.min_base_quality,
        min_support: args.min_support,
        top_barcodes: args.top_barcodes,
        min_distinct_barcodes: args.min_distinct_barcodes,
    };
    if settings.sampling == Sampling::Reservoir && settings.max_records_per_file > 0 {
        eprintln!("Reservoir sampling scans complete input files; at most {} records per file are retained.", settings.max_records_per_file);
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
            stream_seed: settings.seed ^ 0x43444e41,
        });
        libraries.push(LibrarySpec {
            name: library_name(&sample.name, OLIGO_ROLE),
            sample: sample.name.clone(),
            role: OLIGO_ROLE,
            paths: sample.oligo.r1.clone(),
            layout: args.oligo_layout,
            stream_seed: settings.seed ^ 0x4f4c49474f,
        });
    }
    let file_count: usize = libraries.iter().map(|library| library.paths.len()).sum();
    let requested = args.threads.unwrap_or_else(|| auto_threads(file_count));
    if requested == 0 {
        bail!("--threads must be at least 1");
    }
    let threads = requested.min(file_count.max(1));
    eprintln!("Scanning {file_count} FASTQ file(s) with {threads} worker thread(s).");
    let profiles = build_all(libraries, settings, &white, threads)?;
    let (mut cdna, mut oligo) = (Vec::new(), Vec::new());
    for profile in profiles {
        if profile.role == CDNA_ROLE {
            cdna.push(profile);
        } else {
            oligo.push(profile);
        }
    }
    let pairs: Vec<Pair> = cdna
        .iter()
        .flat_map(|c| {
            oligo
                .iter()
                .map(move |o| compare::compare(c, o, settings.min_distinct_barcodes))
        })
        .collect();
    let rankings: Vec<compare::Ranking> = cdna.iter().map(|c| compare::rank(c, &pairs)).collect();
    let pairings: Vec<Pairing> = cdna
        .iter()
        .zip(&rankings)
        .map(|(c, ranking)| {
            let expected = library_name(&c.sample, OLIGO_ROLE);
            let expected_pair = pairs
                .iter()
                .find(|p| p.cdna == c.name && p.oligo == expected);
            let cross = pairs
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
                expected_pair_shared_barcodes: expected_pair.map(|p| p.shared_barcodes),
                best_cross_sample_oligo: cross.map(|p| p.oligo.clone()),
                best_cross_sample_shared_barcodes: cross.map(|p| p.shared_barcodes),
            }
        })
        .collect();
    let report = Report {
        schema_version: 3,
        tool_version: env!("CARGO_PKG_VERSION"),
        method: "exact_whitelisted_CB1_CB2_multi_metric_v3_sample_dirs",
        whitelist_segments: white.len(),
        settings,
        limitations: vec![
            "Multi-metric agreement only: no calibrated same-sample decision or p-value; conflicting metrics have no preferred candidate.",
            "Only technical library compatibility can be assessed; donor identity is not established.",
            "Sample directory names only label the expected cDNA/oligo pair; they are not biological evidence.",
            "Whitelist barcodes are reused across experiments; high overlap alone is not proof.",
            "Coverage uses supported reads, not deduplicated molecules; PCR and depth can bias it.",
            "Independent record sampling can reduce rare-barcode overlap; increase the cap or use all records when support is sparse.",
            "C4 standard/dark R1 only; v1HT non-dark, custom structures and C6 are unsupported.",
            "Reservoir sampling scans entire files; explicit prefix mode is order-sensitive and leaves unscanned tails unvalidated.",
            "Exact whitelist matches only; no mismatch correction and no oligo-derived correction of cDNA.",
            "Only R1 FASTQ files are analysed; R2 files are listed but not validated as pairs.",
        ],
        samples,
        cdna,
        oligo,
        pairs,
        rankings,
        pairings,
    };
    let mut output = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut output, &report)?;
    writeln!(output)?;
    output.as_file().sync_all()?;
    if args.force {
        output.persist(&args.output)?;
    } else {
        output.persist_noclobber(&args.output)?;
    }
    println!("cDNA\toligo\texpected_pair\tshared_barcodes\tcDNA_read_coverage\tbarcode_jaccard\tcount_cosine\tenough_barcodes");
    for pair in &report.pairs {
        let display = |v: Option<f64>| v.map(|n| format!("{n:.6}")).unwrap_or_else(|| "NA".into());
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pair.cdna,
            pair.oligo,
            pair.expected_pair,
            pair.shared_barcodes,
            display(pair.cdna_read_coverage),
            display(pair.barcode_jaccard),
            display(pair.count_cosine),
            pair.sufficient_barcode_count,
        );
    }
    for pairing in &report.pairings {
        eprintln!(
            "Pairing {}: expected={}, preferred={}, status={}, assessment={}",
            pairing.sample,
            pairing.expected_oligo,
            pairing.preferred_oligo.as_deref().unwrap_or("none"),
            pairing.status,
            pairing.assessment
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
