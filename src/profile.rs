use crate::barcode::{extract, Layout, Rejection};
use crate::fastq::Fastq;
use anyhow::{bail, Context, Result};
use clap::ValueEnum;
use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap, HashSet};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, Serialize, PartialEq, Eq, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum Sampling {
    Reservoir,
    Prefix,
}

#[derive(Clone, Copy, Serialize)]
pub struct Settings {
    pub max_records_per_file: usize,
    pub sampling: Sampling,
    pub seed: u64,
    pub min_base_quality: u8,
    pub min_support: u64,
    pub top_barcodes: usize,
    pub min_distinct_barcodes: usize,
}

#[derive(Serialize)]
pub struct FileStats {
    pub path: PathBuf,
    pub layout: Layout,
    pub auto_probe_records: usize,
    pub auto_standard_hits: usize,
    pub auto_dark_hits: usize,
    pub scanned_records: u64,
    pub sampled_records: u64,
    pub scanned_valid_records: u64,
    pub sampling_seed: u64,
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
    let limit = if settings.max_records_per_file == 0 {
        10_000
    } else {
        settings.max_records_per_file.min(10_000)
    };
    let mut reader = Fastq::open(path)?;
    let (mut n, mut standard, mut dark) = (0, 0, 0);
    while n < limit {
        let Some((seq, qual)) = reader.next()? else {
            break;
        };
        n += 1;
        standard += usize::from(
            extract(
                seq,
                qual,
                Layout::Standard,
                settings.min_base_quality,
                white,
            )
            .is_ok(),
        );
        dark +=
            usize::from(extract(seq, qual, Layout::Dark, settings.min_base_quality, white).is_ok());
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

// SplitMix64 is used only for reproducible sampling, never as a biological hash.
struct Reservoir {
    items: Vec<u64>,
    state: u64,
}
impl Reservoir {
    fn next(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn below(&mut self, bound: u64) -> u64 {
        let threshold = bound.wrapping_neg() % bound;
        loop {
            let n = self.next();
            if n >= threshold {
                return n % bound;
            }
        }
    }
    fn offer(&mut self, value: u64, seen: u64, capacity: usize) {
        if self.items.len() < capacity {
            self.items.push(value);
        } else {
            let slot = self.below(seen);
            if slot < capacity as u64 {
                self.items[slot as usize] = value;
            }
        }
    }
}
const SHORT: u64 = u64::MAX;
const AMBIGUOUS: u64 = u64::MAX - 1;
const LOW_QUALITY: u64 = u64::MAX - 2;
const OUTSIDE: u64 = u64::MAX - 3;
fn accumulate(code: u64, stats: &mut FileStats, counts: &mut HashMap<u64, u64>) {
    stats.sampled_records += 1;
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
    pub stream_seed: u64,
}

#[derive(Clone)]
struct FileJob {
    library: usize,
    name: String,
    path: PathBuf,
    layout: Layout,
    stream_seed: u64,
}

/// A file's scan result is independent of execution order: the sampling seed
/// derives from the canonical path and merging is integer addition, so a
/// thread pool reproduces the serial report bit for bit.
type ScanResult = Result<(FileStats, HashMap<u64, u64>)>;

fn scan_file(job: &FileJob, settings: Settings, white: &HashSet<u32>) -> ScanResult {
    let name = &job.name;
    let path = &job.path;
    let layout = job.layout;
    let stream_seed = job.stream_seed;
    let mut counts = HashMap::new();
    let (layout, n, standard, dark) = detect(path, layout, settings, white)
        .with_context(|| format!("library {name}, {}", path.display()))?;
    let mut reader = Fastq::open(path)?;
    // Stable per canonical file path: adding/reordering other candidates
    // must not change this file's sample. FNV is only a stream identifier.
    let file_id = path
        .as_os_str()
        .as_encoded_bytes()
        .iter()
        .fold(0xcbf29ce484222325u64, |h, &byte| {
            (h ^ byte as u64).wrapping_mul(0x100000001b3)
        });
    let seed = stream_seed ^ file_id;
    let mut reservoir = Reservoir {
        items: Vec::new(),
        state: seed,
    };
    let mut stats = FileStats {
        path: path.to_path_buf(),
        layout,
        auto_probe_records: n,
        auto_standard_hits: standard,
        auto_dark_hits: dark,
        scanned_records: 0,
        sampled_records: 0,
        scanned_valid_records: 0,
        sampling_seed: seed,
        input_fully_scanned: false,
        valid_records: 0,
        too_short: 0,
        ambiguous_barcode: 0,
        low_quality_barcode: 0,
        outside_whitelist: 0,
        record_limit_reached: false,
    };
    let sample_limit = settings.max_records_per_file;
    loop {
        if settings.sampling == Sampling::Prefix
            && sample_limit > 0
            && stats.scanned_records >= sample_limit as u64
        {
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
        let code = match extract(seq, qual, layout, settings.min_base_quality, white) {
            Ok(key) => {
                stats.scanned_valid_records += 1;
                key
            }
            Err(Rejection::TooShort) => SHORT,
            Err(Rejection::Ambiguous) => AMBIGUOUS,
            Err(Rejection::LowQuality) => LOW_QUALITY,
            Err(Rejection::OutsideWhitelist) => OUTSIDE,
        };
        if settings.sampling == Sampling::Reservoir && sample_limit > 0 {
            reservoir.offer(code, stats.scanned_records, sample_limit);
        } else {
            accumulate(code, &mut stats, &mut counts);
        }
    }
    for code in reservoir.items {
        accumulate(code, &mut stats, &mut counts);
    }
    stats.record_limit_reached = sample_limit > 0 && stats.sampled_records >= sample_limit as u64;
    Ok((stats, counts))
}

/// Default worker count is bounded so that per-file reservoirs and count
/// tables do not multiply peak memory without an explicit --threads value.
pub fn auto_threads(file_count: usize) -> usize {
    let cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    cpus.min(file_count).clamp(1, 8)
}

/// Merge one file's result immediately so only in-flight files hold a local
/// count table. The merge is integer addition and therefore order-free.
fn absorb(
    index: usize,
    outcome: ScanResult,
    jobs: &[FileJob],
    stats_slots: &mut [Option<FileStats>],
    merged: &mut [HashMap<u64, u64>],
    errors: &mut [Option<anyhow::Error>],
) {
    match outcome {
        Ok((stats, local)) => {
            stats_slots[index] = Some(stats);
            let counts = &mut merged[jobs[index].library];
            for (key, value) in local {
                *counts.entry(key).or_default() += value;
            }
        }
        Err(error) => errors[index] = Some(error),
    }
}

/// Scan every FASTQ file with a bounded thread pool, then merge each library's
/// per-file counts in the original file order. Serial and parallel runs produce
/// identical reports.
pub fn build_all(
    libraries: Vec<LibrarySpec>,
    settings: Settings,
    white: &HashSet<u32>,
    threads: usize,
) -> Result<Vec<Profile>> {
    let mut jobs: Vec<FileJob> = Vec::new();
    for (index, library) in libraries.iter().enumerate() {
        for path in &library.paths {
            jobs.push(FileJob {
                library: index,
                name: library.name.clone(),
                path: path.clone(),
                layout: library.layout,
                stream_seed: library.stream_seed,
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
    if threads <= 1 || jobs.len() <= 1 {
        for (index, job) in jobs.iter().enumerate() {
            absorb(
                index,
                scan_file(job, settings, white),
                &jobs,
                &mut stats_slots,
                &mut merged,
                &mut errors,
            );
        }
    } else {
        use std::collections::VecDeque;
        use std::sync::{mpsc, Mutex};
        let pending: Mutex<VecDeque<usize>> = Mutex::new((0..jobs.len()).collect());
        let (tx, rx) = mpsc::channel();
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
                absorb(
                    index,
                    outcome,
                    &jobs,
                    &mut stats_slots,
                    &mut merged,
                    &mut errors,
                );
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
        profiles.push(assemble(
            library,
            files,
            std::mem::take(&mut merged[index]),
            settings,
        ));
    }
    Ok(profiles)
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
    let top = top_keys(&counts, settings.top_barcodes);
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
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
    #[test]
    fn reservoir_is_reproducible_and_draws_from_entire_stream() {
        let sample = |seed| {
            let mut r = Reservoir {
                items: Vec::new(),
                state: seed,
            };
            for i in 0..10000 {
                r.offer(i, i + 1, 100);
            }
            r.items
        };
        assert_eq!(sample(42), sample(42));
        assert_ne!(sample(42), sample(43));
        let a = sample(42);
        assert_eq!(a.len(), 100);
        assert!(a.iter().any(|&i| i > 9000));
        assert!(a.iter().any(|&i| i < 1000));
        assert_eq!(a.iter().copied().collect::<HashSet<_>>().len(), 100);
    }
}
