use crate::barcode::{extract, Layout, Rejection, MIN_BASE_QUALITY};
use crate::fastq::Fastq;
use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Fixed analysis thresholds. They are deliberately not CLI options.
pub const TOP_BARCODES: usize = 10_000;
pub const MIN_DISTINCT_BARCODES: usize = 100;
const MAX_AUTO_PROBE_RECORDS: usize = 10_000;

#[derive(Clone, Copy, Serialize)]
pub struct Settings {
    /// 0 = analyse every read; otherwise stop after this many reads per file.
    pub reads_per_file: usize,
    pub min_support: u64,
}

#[derive(Serialize)]
pub struct FileStats {
    pub path: PathBuf,
    pub layout: Layout,
    pub auto_probe_records: usize,
    pub auto_standard_hits: usize,
    pub auto_dark_hits: usize,
    pub scanned_records: u64,
    pub input_fully_scanned: bool,
    pub valid_records: u64,
    pub too_short: u64,
    pub ambiguous_barcode: u64,
    pub low_quality_barcode: u64,
    pub outside_whitelist: u64,
    pub record_limit_reached: bool,
}

#[derive(Serialize)]
pub struct Profile {
    pub name: String,
    pub sample: String,
    pub role: &'static str,
    pub files: Vec<FileStats>,
    pub observed_barcodes: usize,
    pub supported_barcodes: usize,
    pub supported_reads: u64,
    #[serde(skip)]
    pub counts: HashMap<u64, u64>,
    #[serde(skip)]
    pub top: HashSet<u64>,
    #[serde(skip)]
    pub norm: f64,
    #[serde(skip)]
    pub log_weights: HashMap<u64, u64>,
    #[serde(skip)]
    pub log_norm: f64,
    pub abundance_tiers: Vec<AbundanceTier>,
    /// (sum counts)^2 / sum(count^2); near barcode count means flat abundance.
    pub count_effective_barcodes: Option<f64>,
    pub largest_barcode_read_fraction: Option<f64>,
}

#[derive(Serialize)]
pub struct AbundanceTier {
    pub requested_top: usize,
    pub minimum_count: u64,
    pub included_barcodes: usize,
}

/// Cache by read count, not barcode. Fixed-point log weights make integer dot
/// products independent of HashMap/worker order (absolute weight error <= 5e-10).
pub fn log_summary(counts: &HashMap<u64, u64>) -> (HashMap<u64, u64>, f64) {
    let mut weights = HashMap::new();
    let mut norm_squared = 0u128;
    for &count in counts.values() {
        let weight = *weights
            .entry(count)
            .or_insert_with(|| ((count as f64).ln_1p() * 1_000_000_000.0).round() as u64);
        norm_squared += weight as u128 * weight as u128;
    }
    (weights, (norm_squared as f64).sqrt())
}

fn abundance_tiers(counts: &HashMap<u64, u64>, keys: &HashSet<u64>) -> Vec<AbundanceTier> {
    // One bounded heap supplies all cutoffs. Include every boundary tie rather
    // than selecting barcode IDs arbitrarily at the requested rank.
    let mut values: Vec<_> = keys.iter().map(|key| counts[key]).collect();
    values.sort_unstable_by(|a, b| b.cmp(a));
    [100, 1000, TOP_BARCODES]
        .into_iter()
        .map(|requested_top| {
            let minimum_count = values
                .get(requested_top.min(values.len()).saturating_sub(1))
                .copied()
                .unwrap_or(0);
            let included_barcodes = counts.values().filter(|&&n| n >= minimum_count).count();
            AbundanceTier {
                requested_top,
                minimum_count,
                included_barcodes,
            }
        })
        .collect()
}

fn detect(
    path: &Path,
    layout: Layout,
    settings: Settings,
    white: &HashSet<u32>,
) -> Result<(Layout, usize, usize, usize)> {
    if layout != Layout::Auto {
        return Ok((layout, 0, 0, 0));
    }
    let limit = if settings.reads_per_file == 0 {
        MAX_AUTO_PROBE_RECORDS
    } else {
        settings.reads_per_file.min(MAX_AUTO_PROBE_RECORDS)
    };
    let mut reader = Fastq::open(path)?;
    let (mut n, mut standard, mut dark) = (0, 0, 0);
    while n < limit {
        let Some((seq, qual)) = reader.next()? else {
            break;
        };
        n += 1;
        standard +=
            usize::from(extract(seq, qual, Layout::Standard, MIN_BASE_QUALITY, white).is_ok());
        dark += usize::from(extract(seq, qual, Layout::Dark, MIN_BASE_QUALITY, white).is_ok());
    }
    let best = standard.max(dark);
    let second = standard.min(dark);
    if n == 0
        || best < n.min(20)
        || (best as f64) / (n as f64) < 0.1
        || best < second.saturating_mul(3).max(1)
    {
        bail!("layout detection inconclusive: {n} records, standard={standard}, dark={dark}; check quality/chemistry or specify the layout explicitly");
    }
    Ok((
        if standard > dark {
            Layout::Standard
        } else {
            Layout::Dark
        },
        n,
        standard,
        dark,
    ))
}

// Rejection codes are kept in the same stream as barcode keys so that every
// record is accounted for exactly once.
const SHORT: u64 = u64::MAX;
const AMBIGUOUS: u64 = u64::MAX - 1;
const LOW_QUALITY: u64 = u64::MAX - 2;
const OUTSIDE: u64 = u64::MAX - 3;

fn accumulate(code: u64, stats: &mut FileStats, counts: &mut HashMap<u64, u64>) {
    match code {
        SHORT => stats.too_short += 1,
        AMBIGUOUS => stats.ambiguous_barcode += 1,
        LOW_QUALITY => stats.low_quality_barcode += 1,
        OUTSIDE => stats.outside_whitelist += 1,
        key => {
            *counts.entry(key).or_default() += 1;
            stats.valid_records += 1;
        }
    }
}

fn top_keys(counts: &HashMap<u64, u64>, n: usize) -> HashSet<u64> {
    // Larger score is better; the heap root is the weakest retained entry.
    let mut heap = BinaryHeap::new();
    for (&key, &count) in counts {
        let candidate = (count, Reverse(key));
        if heap.len() < n {
            heap.push(Reverse(candidate));
        } else if heap.peek().is_some_and(|smallest| candidate > smallest.0) {
            heap.pop();
            heap.push(Reverse(candidate));
        }
    }
    heap.into_iter()
        .map(|Reverse((_, Reverse(key)))| key)
        .collect()
}

pub struct LibrarySpec {
    pub name: String,
    pub sample: String,
    pub role: &'static str,
    pub paths: Vec<PathBuf>,
    pub layout: Layout,
}

#[derive(Clone)]
struct FileJob {
    library: usize,
    name: String,
    path: PathBuf,
    layout: Layout,
}

#[derive(Clone, Serialize)]
pub struct FileTiming {
    pub library: String,
    pub path: PathBuf,
    pub detect_seconds: f64,
    pub read_seconds: f64,
    pub scan_seconds: f64,
    pub scanned_records: u64,
}

#[derive(Serialize)]
pub struct BuildTimings {
    pub scan_wall_seconds: f64,
    pub merge_seconds: f64,
    pub files: Vec<FileTiming>,
}

type ScanResult = Result<(FileStats, HashMap<u64, u64>, FileTiming)>;

/// Scan one FASTQ file. Each file is independent, so parallel execution gives
/// the same counts as a serial run (merging is integer addition).
fn scan_file(job: &FileJob, settings: Settings, white: &HashSet<u32>) -> ScanResult {
    let name = &job.name;
    let path = &job.path;
    let layout = job.layout;
    let started = Instant::now();
    let detect_started = Instant::now();
    let (layout, n, standard, dark) = detect(path, layout, settings, white)
        .with_context(|| format!("library {name}, {}", path.display()))?;
    let detect_seconds = detect_started.elapsed().as_secs_f64();
    let read_started = Instant::now();
    let mut reader = Fastq::open(path)?;
    let mut counts = HashMap::new();
    let mut stats = FileStats {
        path: path.to_path_buf(),
        layout,
        auto_probe_records: n,
        auto_standard_hits: standard,
        auto_dark_hits: dark,
        scanned_records: 0,
        input_fully_scanned: false,
        valid_records: 0,
        too_short: 0,
        ambiguous_barcode: 0,
        low_quality_barcode: 0,
        outside_whitelist: 0,
        record_limit_reached: false,
    };
    loop {
        if settings.reads_per_file > 0 && stats.scanned_records >= settings.reads_per_file as u64 {
            break;
        }
        let Some((seq, qual)) = reader
            .next()
            .with_context(|| format!("library {name}, {}", path.display()))?
        else {
            stats.input_fully_scanned = true;
            break;
        };
        stats.scanned_records += 1;
        let code = match extract(seq, qual, layout, MIN_BASE_QUALITY, white) {
            Ok(key) => key,
            Err(Rejection::TooShort) => SHORT,
            Err(Rejection::Ambiguous) => AMBIGUOUS,
            Err(Rejection::LowQuality) => LOW_QUALITY,
            Err(Rejection::OutsideWhitelist) => OUTSIDE,
        };
        accumulate(code, &mut stats, &mut counts);
    }
    let read_seconds = read_started.elapsed().as_secs_f64();
    stats.record_limit_reached =
        settings.reads_per_file > 0 && stats.scanned_records >= settings.reads_per_file as u64;
    let timing = FileTiming {
        library: job.name.clone(),
        path: path.to_path_buf(),
        detect_seconds,
        read_seconds,
        scan_seconds: started.elapsed().as_secs_f64(),
        scanned_records: stats.scanned_records,
    };
    Ok((stats, counts, timing))
}

/// Default worker count is bounded so that per-file count tables do not
/// multiply peak memory without an explicit --threads value.
pub fn auto_threads(file_count: usize) -> usize {
    let cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    cpus.min(file_count).clamp(1, 8)
}

/// Merge a completed file, adopting its allocation when the destination is empty.
/// Integer addition keeps counts independent of completion order.
fn absorb(
    index: usize,
    outcome: ScanResult,
    jobs: &[FileJob],
    stats_slots: &mut [Option<FileStats>],
    merged: &mut [HashMap<u64, u64>],
    errors: &mut [Option<anyhow::Error>],
    timing_slots: &mut [Option<FileTiming>],
) {
    match outcome {
        Ok((stats, local, timing)) => {
            stats_slots[index] = Some(stats);
            timing_slots[index] = Some(timing);
            let counts = &mut merged[jobs[index].library];
            if counts.is_empty() {
                *counts = local;
            } else {
                for (key, value) in local {
                    *counts.entry(key).or_default() += value;
                }
            }
        }
        Err(error) => errors[index] = Some(error),
    }
}

/// Scan every FASTQ file with a bounded thread pool, then merge each library's
/// per-file counts as results arrive; reports retain the original file order. Serial and parallel runs produce
/// identical reports. Wall/merge timings are returned for the JSON report.
pub fn build_all(
    libraries: Vec<LibrarySpec>,
    settings: Settings,
    white: &HashSet<u32>,
    threads: usize,
) -> Result<(Vec<Profile>, BuildTimings)> {
    let build_started = Instant::now();
    let mut jobs: Vec<FileJob> = Vec::new();
    for (index, library) in libraries.iter().enumerate() {
        for path in &library.paths {
            jobs.push(FileJob {
                library: index,
                name: library.name.clone(),
                path: path.clone(),
                layout: library.layout,
            });
        }
    }
    let mut per_library: Vec<Vec<usize>> = (0..libraries.len()).map(|_| Vec::new()).collect();
    for (index, job) in jobs.iter().enumerate() {
        per_library[job.library].push(index);
    }
    let mut stats_slots: Vec<Option<FileStats>> = (0..jobs.len()).map(|_| None).collect();
    let mut merged: Vec<HashMap<u64, u64>> = (0..libraries.len()).map(|_| HashMap::new()).collect();
    let mut errors: Vec<Option<anyhow::Error>> = (0..jobs.len()).map(|_| None).collect();
    let mut timing_slots: Vec<Option<FileTiming>> = (0..jobs.len()).map(|_| None).collect();
    let mut merge_seconds = 0.0;
    if threads <= 1 || jobs.len() <= 1 {
        for (index, job) in jobs.iter().enumerate() {
            let outcome = scan_file(job, settings, white);
            let merge_started = Instant::now();
            absorb(
                index,
                outcome,
                &jobs,
                &mut stats_slots,
                &mut merged,
                &mut errors,
                &mut timing_slots,
            );
            merge_seconds += merge_started.elapsed().as_secs_f64();
        }
    } else {
        use std::collections::VecDeque;
        use std::sync::{mpsc, Mutex};
        let pending: Mutex<VecDeque<usize>> = Mutex::new((0..jobs.len()).collect());
        // Rendezvous: workers retain at most one local table each; completed
        // tables cannot accumulate in an unbounded result queue.
        let (tx, rx) = mpsc::sync_channel(0);
        std::thread::scope(|scope| {
            for _ in 0..threads.min(jobs.len()) {
                let tx = tx.clone();
                let pending = &pending;
                let jobs = &jobs;
                scope.spawn(move || loop {
                    let next = pending
                        .lock()
                        .unwrap_or_else(|error| error.into_inner())
                        .pop_front();
                    let Some(index) = next else { break };
                    let outcome = scan_file(&jobs[index], settings, white);
                    if tx.send((index, outcome)).is_err() {
                        break;
                    }
                });
            }
            drop(tx);
            for (index, outcome) in rx {
                let merge_started = Instant::now();
                absorb(
                    index,
                    outcome,
                    &jobs,
                    &mut stats_slots,
                    &mut merged,
                    &mut errors,
                    &mut timing_slots,
                );
                merge_seconds += merge_started.elapsed().as_secs_f64();
            }
        });
    }
    let mut profiles = Vec::with_capacity(libraries.len());
    for (index, library) in libraries.iter().enumerate() {
        let mut files = Vec::new();
        // Report failures in the original file order, as a serial run would.
        for &job_index in &per_library[index] {
            if let Some(error) = errors[job_index].take() {
                return Err(error);
            }
            files.push(
                stats_slots[job_index]
                    .take()
                    .ok_or_else(|| anyhow::anyhow!("missing scan result for job {job_index}"))?,
            );
        }
        let assemble_started = Instant::now();
        profiles.push(assemble(
            library,
            files,
            std::mem::take(&mut merged[index]),
            settings,
        ));
        merge_seconds += assemble_started.elapsed().as_secs_f64();
    }
    Ok((
        profiles,
        BuildTimings {
            scan_wall_seconds: build_started.elapsed().as_secs_f64(),
            merge_seconds,
            files: timing_slots.into_iter().flatten().collect(),
        },
    ))
}

fn assemble(
    library: &LibrarySpec,
    files: Vec<FileStats>,
    mut counts: HashMap<u64, u64>,
    settings: Settings,
) -> Profile {
    let observed_barcodes = counts.len();
    counts.retain(|_, count| *count >= settings.min_support);
    if counts.is_empty() {
        counts = HashMap::new();
    } else if counts.len() < counts.capacity() / 3 {
        counts.shrink_to_fit();
    }
    let supported_reads = counts.values().sum();
    let norm = (counts
        .values()
        .map(|&n| (n as u128) * (n as u128))
        .sum::<u128>() as f64)
        .sqrt();
    let top = top_keys(&counts, TOP_BARCODES);
    let (log_weights, log_norm) = log_summary(&counts);
    let abundance_tiers = abundance_tiers(&counts, &top);
    let count_effective_barcodes = (norm > 0.0).then(|| {
        (supported_reads as f64 / norm)
            .powi(2)
            .min(counts.len() as f64)
    });
    let largest_barcode_read_fraction = (supported_reads > 0)
        .then(|| counts.values().copied().max().unwrap_or(0) as f64 / supported_reads as f64);
    Profile {
        name: library.name.clone(),
        sample: library.sample.clone(),
        role: library.role,
        files,
        observed_barcodes,
        supported_barcodes: counts.len(),
        supported_reads,
        counts,
        top,
        norm,
        log_weights,
        log_norm,
        abundance_tiers,
        count_effective_barcodes,
        largest_barcode_read_fraction,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn abundance_cutoffs_include_all_boundary_ties_and_empty_is_safe() {
        let counts: HashMap<_, _> = (0..1200)
            .map(|i| (i, if i < 20 { 30 } else { 3 }))
            .collect();
        let tiers = abundance_tiers(&counts, &top_keys(&counts, TOP_BARCODES));
        assert_eq!(tiers[0].minimum_count, 3);
        assert_eq!(tiers[0].included_barcodes, 1200);
        let empty = HashMap::new();
        assert!(abundance_tiers(&empty, &HashSet::new())
            .iter()
            .all(|t| t.included_barcodes == 0));
        assert_eq!(log_summary(&empty).1, 0.0);
    }

    #[test]
    fn top_heap_matches_full_sort_including_ties() {
        let counts: HashMap<_, _> = (0..1000).map(|i| (i, i % 7)).collect();
        let mut expected: Vec<_> = counts.iter().map(|(&k, &v)| (k, v)).collect();
        expected.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        for n in [1, 17, 1000, 2000] {
            assert_eq!(
                top_keys(&counts, n),
                expected.iter().take(n).map(|p| p.0).collect()
            );
        }
    }
}
