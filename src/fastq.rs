use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub struct Fastq {
    reader: Box<dyn BufRead>,
    record: u64,
    lines: [Vec<u8>; 4],
}

impl Fastq {
    pub fn open(path: &Path) -> Result<Self> {
        let mut source = BufReader::new(
            File::open(path).with_context(|| format!("opening {}", path.display()))?,
        );
        let gzip = source.fill_buf()?.starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn Read> = if gzip {
            Box::new(MultiGzDecoder::new(source))
        } else {
            Box::new(source)
        };
        Ok(Self {
            reader: Box::new(BufReader::new(reader)),
            record: 0,
            lines: std::array::from_fn(|_| Vec::new()),
        })
    }

    pub fn next(&mut self) -> Result<Option<(&[u8], &[u8])>> {
        for line in &mut self.lines {
            line.clear();
        }
        if self.reader.read_until(b'\n', &mut self.lines[0])? == 0 {
            return Ok(None);
        }
        self.record += 1;
        for line in &mut self.lines[1..] {
            if self.reader.read_until(b'\n', line)? == 0 {
                bail!("truncated FASTQ at record {}", self.record);
            }
        }
        for line in &mut self.lines {
            while matches!(line.last(), Some(b'\n' | b'\r')) {
                line.pop();
            }
        }
        if !self.lines[0].starts_with(b"@") || !self.lines[2].starts_with(b"+") {
            bail!(
                "invalid four-line FASTQ at record {} (wrapped FASTQ is unsupported)",
                self.record
            );
        }
        self.lines[1].make_ascii_uppercase();
        let seq = &self.lines[1];
        let qual = &self.lines[3];
        if seq.is_empty()
            || seq.len() != qual.len()
            || !seq.is_ascii()
            || !qual.iter().all(|&q| (33..=126).contains(&q))
        {
            bail!(
                "invalid sequence/Phred+33 quality lengths or encoding at record {}",
                self.record
            );
        }
        Ok(Some((seq, qual)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{write::GzEncoder, Compression};
    use std::io::Write;

    #[test]
    fn concatenated_gzip_is_read_by_magic_not_extension() -> Result<()> {
        let mut file = tempfile::NamedTempFile::new()?;
        for name in ["a", "b"] {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            write!(encoder, "@{name}\nACGT\n+\nIIII\n")?;
            file.write_all(&encoder.finish()?)?;
        }
        let mut reader = Fastq::open(file.path())?;
        assert!(reader.next()?.is_some());
        assert!(reader.next()?.is_some());
        assert!(reader.next()?.is_none());
        Ok(())
    }

    #[test]
    fn malformed_and_truncated_records_fail() -> Result<()> {
        for body in ["@a\nACGT\n+\n", "@a\nACGT\n+\nIII\n", "a\nACGT\n+\nIIII\n"] {
            let mut file = tempfile::NamedTempFile::new()?;
            file.write_all(body.as_bytes())?;
            assert!(Fastq::open(file.path())?.next().is_err());
        }
        Ok(())
    }
}
