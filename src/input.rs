//! Sample-directory discovery following the dnbc4tools_py `--fastqs` layout.
//!
//! Each sample directory owns `cDNA/` and `oligo/` subdirectories with R1/R2
//! FASTQ files. Only R1 files are analysed; R2 files are recorded for review.

use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

pub const CDNA_ROLE: &str = "cDNA";
pub const OLIGO_ROLE: &str = "oligo";

const CDNA_DIRS: [&str; 3] = ["cDNA", "cdna", "CDNA"];
const OLIGO_DIRS: [&str; 3] = ["oligo", "Oligo", "OLIGO"];
const FASTQ_EXTENSIONS: [&str; 4] = [".fastq.gz", ".fq.gz", ".fastq", ".fq"];
const R1_TOKENS: [&str; 5] = ["_R1_", "_R1", "_1", "R1", "_read1"];
const R2_TOKENS: [&str; 5] = ["_R2_", "_R2", "_2", "R2", "_read2"];

pub fn library_name(sample: &str, role: &str) -> String {
    format!("{sample}/{role}")
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ReadType {
    R1,
    R2,
}

#[derive(Serialize)]
pub struct LibraryFiles {
    pub directory: PathBuf,
    pub r1: Vec<PathBuf>,
    pub r2: Vec<PathBuf>,
}

#[derive(Serialize)]
pub struct Sample {
    pub name: String,
    pub directory: PathBuf,
    pub cdna: LibraryFiles,
    pub oligo: LibraryFiles,
}

/// Last boundary-aware match of any read token, mirroring fastq_naming.py.
fn last_token_pos(filename: &str, tokens: &[&str]) -> Option<usize> {
    let mut last: Option<usize> = None;
    for token in tokens {
        let token = token.to_ascii_lowercase();
        for (start, _) in filename.match_indices(&token) {
            let end = start + token.len();
            let suffix_ok = token.ends_with('_')
                || filename
                    .as_bytes()
                    .get(end)
                    .is_none_or(|b| *b == b'.' || *b == b'_');
            let prefix_ok = token.starts_with('_')
                || start == 0
                || !matches!(
                    filename.as_bytes()[start - 1],
                    b'a'..=b'z' | b'0'..=b'9' | b'_'
                );
            if suffix_ok && prefix_ok {
                last = Some(last.map_or(start, |previous| previous.max(start)));
            }
        }
    }
    last
}

fn classify(filename: &str) -> Option<ReadType> {
    let lower = filename.to_ascii_lowercase();
    match (
        last_token_pos(&lower, &R1_TOKENS),
        last_token_pos(&lower, &R2_TOKENS),
    ) {
        (Some(r1), Some(r2)) => Some(if r2 > r1 { ReadType::R2 } else { ReadType::R1 }),
        (Some(_), None) => Some(ReadType::R1),
        (None, Some(_)) => Some(ReadType::R2),
        (None, None) => None,
    }
}

fn is_sample_dir(dir: &Path) -> bool {
    CDNA_DIRS
        .iter()
        .chain(&OLIGO_DIRS)
        .any(|name| dir.join(name).is_dir())
}

fn find_library_dir(sample_dir: &Path, aliases: &[&str], label: &str) -> Result<PathBuf> {
    for name in aliases {
        let candidate = sample_dir.join(name);
        if candidate.is_dir() {
            return Ok(candidate);
        }
    }
    let available: Vec<String> = fs::read_dir(sample_dir)
        .with_context(|| format!("listing {}", sample_dir.display()))?
        .filter_map(Result::ok)
        .filter(|entry| entry.path().is_dir())
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .collect();
    bail!(
        "{label} subdirectory not found in {} (searched {aliases:?}; available directories: {available:?})",
        sample_dir.display()
    );
}

fn collect_library(dir: &Path, label: &str) -> Result<LibraryFiles> {
    let mut files: Vec<(String, PathBuf)> = Vec::new();
    for entry in fs::read_dir(dir).with_context(|| format!("listing {}", dir.display()))? {
        let entry = entry?;
        let path = entry.path();
        // Cluster datasets usually use symlinks; Path::is_file follows them,
        // while DirEntry::file_type does not. Broken links are skipped.
        if !path.is_file() {
            continue;
        }
        // Classify by the directory entry name: a link target may be named
        // differently (e.g. an opaque basecall file without R1/R2 tokens).
        let name = entry.file_name().to_string_lossy().into_owned();
        let lower = name.to_ascii_lowercase();
        if FASTQ_EXTENSIONS.iter().any(|ext| lower.ends_with(ext)) {
            files.push((
                name,
                path.canonicalize()
                    .with_context(|| format!("resolving {}", path.display()))?,
            ));
        }
    }
    files.sort_by(|a, b| a.0.cmp(&b.0));
    if files.is_empty() {
        bail!(
            "no FASTQ files found in {label} directory {}",
            dir.display()
        );
    }

    let (mut r1, mut r2, mut unclassified): (Vec<PathBuf>, Vec<PathBuf>, Vec<PathBuf>) =
        (Vec::new(), Vec::new(), Vec::new());
    for (name, file) in &files {
        match classify(name) {
            Some(ReadType::R1) => r1.push(file.clone()),
            Some(ReadType::R2) => r2.push(file.clone()),
            None => unclassified.push(file.clone()),
        }
    }
    // A lone unlabelled FASTQ is necessarily the R1 of this library.
    if files.len() == 1 && unclassified.len() == 1 {
        r1 = unclassified.clone();
        unclassified.clear();
    }
    if !unclassified.is_empty() {
        let names: Vec<_> = unclassified
            .iter()
            .map(|p| {
                p.file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .into_owned()
            })
            .collect();
        bail!(
            "cannot determine R1/R2 for {names:?} in {}; use _R1_/_R2_ filename tokens",
            dir.display()
        );
    }
    if r1.is_empty() {
        bail!(
            "no R1 FASTQ found in {label} directory {} ({} R2 file(s), which are not analysed)",
            dir.display(),
            r2.len()
        );
    }
    Ok(LibraryFiles {
        directory: dir.to_path_buf(),
        r1,
        r2,
    })
}

/// Resolve input roots into samples. A root is either one sample directory or a
/// directory whose immediate children are sample directories.
pub fn discover(roots: &[PathBuf]) -> Result<Vec<Sample>> {
    let mut sample_dirs: Vec<(String, PathBuf)> = Vec::new();
    for root in roots {
        let root = root
            .canonicalize()
            .with_context(|| format!("resolving input directory {}", root.display()))?;
        if !root.is_dir() {
            bail!("input is not a directory: {}", root.display());
        }
        if is_sample_dir(&root) {
            let name = root
                .file_name()
                .context("input sample directory has no name")?
                .to_string_lossy()
                .into_owned();
            sample_dirs.push((name, root));
            continue;
        }
        let mut children: Vec<PathBuf> = fs::read_dir(&root)
            .with_context(|| format!("listing {}", root.display()))?
            .filter_map(Result::ok)
            .map(|entry| entry.path())
            .filter(|path| path.is_dir() && is_sample_dir(path))
            .collect();
        children.sort();
        if children.is_empty() {
            bail!(
                "no sample directories under {}; expected <root>/<sample>/{{cDNA,oligo}}/<..._R1_...>.fastq[.gz]",
                root.display()
            );
        }
        for child in children {
            let name = child
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned();
            sample_dirs.push((name, child));
        }
    }
    sample_dirs.sort_by(|a, b| a.0.cmp(&b.0));

    let mut seen_names: HashSet<String> = HashSet::new();
    let mut seen_files: HashSet<PathBuf> = HashSet::new();
    let mut samples = Vec::new();
    for (name, dir) in sample_dirs {
        if !seen_names.insert(name.clone()) {
            bail!("duplicate sample name {name:?}; sample directory names must be unique");
        }
        let cdna_dir = find_library_dir(&dir, &CDNA_DIRS, CDNA_ROLE)?;
        let oligo_dir = find_library_dir(&dir, &OLIGO_DIRS, OLIGO_ROLE)?;
        let cdna = collect_library(&cdna_dir, CDNA_ROLE)?;
        let oligo = collect_library(&oligo_dir, OLIGO_ROLE)?;
        for file in cdna
            .r1
            .iter()
            .chain(&cdna.r2)
            .chain(&oligo.r1)
            .chain(&oligo.r2)
        {
            if !seen_files.insert(file.clone()) {
                bail!(
                    "FASTQ is assigned more than once (possibly across samples/libraries): {}",
                    file.display()
                );
            }
        }
        samples.push(Sample {
            name,
            directory: dir,
            cdna,
            oligo,
        });
    }
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_tokens_follow_dnbc4tools_naming() {
        assert_eq!(classify("S1_L001_R1_001.fastq.gz"), Some(ReadType::R1));
        assert_eq!(classify("S1_L001_R2_001.fastq.gz"), Some(ReadType::R2));
        assert_eq!(classify("sample_1.fq.gz"), Some(ReadType::R1));
        assert_eq!(classify("sample_2.fq.gz"), Some(ReadType::R2));
        assert_eq!(classify("lane_R1.fastq"), Some(ReadType::R1));
        assert_eq!(classify("Sample_1_2.fq.gz"), Some(ReadType::R2));
        assert_eq!(classify("unlabelled.fq.gz"), None);
    }
}
