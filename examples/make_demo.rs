//! Deterministic synthetic data: three sample directories following the
//! dnbc4tools_py `--fastqs` layout. SampleA/SampleB oligo files are swapped.
use flate2::{write::GzEncoder, Compression};
use std::{
    env, fs,
    io::Write,
    path::{Path, PathBuf},
};

fn sequences(start: usize, count: usize, dark: bool, white: &[&str]) -> Vec<String> {
    (start..start + count)
        .map(|i| {
            if dark {
                format!("{}{}", white[i], white[500 + i])
            } else {
                format!(
                    "TCTGCG{}CCTTCC{}CGATGACGTACGTACTTTTTTTTTTT",
                    white[i],
                    white[500 + i]
                )
            }
        })
        .collect()
}

fn fastq(path: &Path, sequences: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let mut body = String::new();
    for (i, seq) in sequences.iter().enumerate() {
        for duplicate in 0..3 {
            body.push_str(&format!(
                "@read{i}_{duplicate}\n{seq}\n+\n{}\n",
                "I".repeat(seq.len())
            ));
        }
    }
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(body.as_bytes())?;
    fs::write(path, encoder.finish()?)?;
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let destination = PathBuf::from(env::args().nth(1).ok_or("usage: make_demo NEW_DIRECTORY")?);
    fs::create_dir(&destination)?; // Do not overwrite an existing example directory.
    let white: Vec<_> = include_str!("../assets/barcode_whitelist.txt")
        .lines()
        .collect();
    for (sample, cdna_start, oligo_start) in [
        ("SampleA", 0, 200),   // oligo actually belongs to SampleB
        ("SampleB", 200, 0),   // oligo actually belongs to SampleA
        ("SampleC", 400, 400), // correctly paired
    ] {
        let cdna_dir = destination.join(sample).join("cDNA");
        let oligo_dir = destination.join(sample).join("oligo");
        fs::create_dir_all(&cdna_dir)?;
        fs::create_dir_all(&oligo_dir)?;
        let cdna = sequences(cdna_start, 120, false, &white);
        let oligo = sequences(oligo_start, 120, true, &white);
        for (dir, sequences, read) in [
            (&cdna_dir, &cdna, "R1"),
            (&oligo_dir, &oligo, "R1"),
            (&cdna_dir, &cdna, "R2"),
            (&oligo_dir, &oligo, "R2"),
        ] {
            fastq(
                &dir.join(format!("{sample}_L001_{read}_001.fastq.gz")),
                sequences,
            )?;
        }
    }
    println!(
        "Created three synthetic sample directories in {}",
        destination.display()
    );
    Ok(())
}
