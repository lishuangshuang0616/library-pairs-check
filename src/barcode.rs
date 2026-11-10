use anyhow::{bail, Result};
use clap::ValueEnum;
use serde::Serialize;
use std::collections::HashSet;

/// Fixed minimum Phred+33 quality at every barcode base (Q20).
pub const MIN_BASE_QUALITY: u8 = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, ValueEnum)]
#[serde(rename_all = "snake_case")]
pub enum Layout {
    Auto,
    Standard,
    Dark,
}

impl Layout {
    pub fn offsets(self) -> (usize, usize) {
        match self {
            Self::Standard => (6, 22),
            Self::Dark => (0, 10),
            Self::Auto => unreachable!("auto must be resolved before counting"),
        }
    }
}

pub fn encode(seq: &[u8]) -> Option<u32> {
    seq.iter().try_fold(0, |n, b| {
        Some(
            (n << 2)
                | match b {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => return None,
                },
        )
    })
}

pub fn whitelist() -> Result<HashSet<u32>> {
    let mut result = HashSet::new();
    for line in include_str!("../assets/barcode_whitelist.txt").lines() {
        if line.len() != 10 {
            bail!("invalid bundled whitelist length");
        }
        result.insert(
            encode(line.as_bytes())
                .ok_or_else(|| anyhow::anyhow!("invalid bundled whitelist sequence"))?,
        );
    }
    Ok(result)
}

#[derive(Debug, PartialEq, Eq)]
pub enum Rejection {
    TooShort,
    Ambiguous,
    LowQuality,
    OutsideWhitelist,
}

pub fn extract(
    seq: &[u8],
    qual: &[u8],
    layout: Layout,
    min_quality: u8,
    whitelist: &HashSet<u32>,
) -> Result<u64, Rejection> {
    let (a, b) = layout.offsets();
    if seq.len() < b + 10 {
        return Err(Rejection::TooShort);
    }
    let left = encode(&seq[a..a + 10]).ok_or(Rejection::Ambiguous)?;
    let right = encode(&seq[b..b + 10]).ok_or(Rejection::Ambiguous)?;
    if qual[a..a + 10]
        .iter()
        .chain(&qual[b..b + 10])
        .any(|q| q - 33 < min_quality)
    {
        return Err(Rejection::LowQuality);
    }
    if !whitelist.contains(&left) || !whitelist.contains(&right) {
        return Err(Rejection::OutsideWhitelist);
    }
    Ok(((left as u64) << 20) | right as u64)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn layouts_agree_and_quality_and_n_are_rejected() -> anyhow::Result<()> {
        let white = whitelist()?;
        let left = "TAACAGCCAA";
        let right = "CTAAGAGTCC";
        let dark = format!("{left}{right}");
        let standard = format!("TCTGCG{left}CCTTCC{right}");
        assert_eq!(
            extract(dark.as_bytes(), &[b'I'; 20], Layout::Dark, 20, &white),
            extract(
                standard.as_bytes(),
                &[b'I'; 32],
                Layout::Standard,
                20,
                &white
            )
        );
        assert_eq!(
            extract(dark.as_bytes(), &[b'!'; 20], Layout::Dark, 20, &white),
            Err(Rejection::LowQuality)
        );
        assert_eq!(
            extract(
                b"NAACAGCCAACTAAGAGTCC",
                &[b'I'; 20],
                Layout::Dark,
                20,
                &white
            ),
            Err(Rejection::Ambiguous)
        );
        Ok(())
    }
}
