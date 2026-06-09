//! Sequence statistics: base composition and GC content for FASTA/FASTQ.
//!
//! This is the CPU reference implementation. It is also the natural unit of
//! GPU-parallel work in this tool (a reduction over millions of bases), so the
//! CLI can offload it to a GPU device pinned by HydraMPP — see `cleaver-cli`.

use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::formats::{self, Format};

/// Base composition over the sequence characters of a file.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct SeqStats {
    pub seqs: u64,
    pub bases: u64,
    pub gc: u64,
    pub a: u64,
    pub c: u64,
    pub g: u64,
    pub t: u64,
    pub n: u64,
    pub other: u64,
}

impl SeqStats {
    pub fn gc_percent(&self) -> f64 {
        if self.bases == 0 {
            0.0
        } else {
            100.0 * self.gc as f64 / self.bases as f64
        }
    }

    /// Tally one sequence line's bytes (newlines ignored).
    pub fn add_seq_line(&mut self, line: &[u8]) {
        for &b in line {
            match b {
                b'\n' | b'\r' => {}
                b'A' | b'a' => {
                    self.a += 1;
                    self.bases += 1;
                }
                b'C' | b'c' => {
                    self.c += 1;
                    self.gc += 1;
                    self.bases += 1;
                }
                b'G' | b'g' => {
                    self.g += 1;
                    self.gc += 1;
                    self.bases += 1;
                }
                b'T' | b't' | b'U' | b'u' => {
                    self.t += 1;
                    self.bases += 1;
                }
                b'N' | b'n' => {
                    self.n += 1;
                    self.bases += 1;
                }
                _ => {
                    self.other += 1;
                    self.bases += 1;
                }
            }
        }
    }
}

/// Compute [`SeqStats`] for a FASTA or FASTQ file (gzip transparent), CPU path.
pub fn stats_file(input: &Path) -> Result<SeqStats> {
    let fmt = formats::detect(input)?;
    match fmt {
        Format::Fasta => stats_fasta(input),
        Format::Fastq => stats_fastq(input),
        other => bail!("stats supports FASTA/FASTQ, not {other:?}"),
    }
}

fn stats_fasta(input: &Path) -> Result<SeqStats> {
    let (mut r, _gz) = crate::open_reader(input)?;
    let mut s = SeqStats::default();
    let mut line = Vec::with_capacity(256);
    loop {
        line.clear();
        let n = r.read_until(b'\n', &mut line).context("reading FASTA")?;
        if n == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            s.seqs += 1;
        } else {
            s.add_seq_line(&line);
        }
    }
    Ok(s)
}

fn stats_fastq(input: &Path) -> Result<SeqStats> {
    let (mut r, _gz) = crate::open_reader(input)?;
    let mut s = SeqStats::default();
    let mut lines: [Vec<u8>; 4] = Default::default();
    loop {
        let mut got = 0;
        for slot in lines.iter_mut() {
            slot.clear();
            if r.read_until(b'\n', slot).context("reading FASTQ")? == 0 {
                break;
            }
            got += 1;
        }
        if got == 0 {
            break;
        }
        if got != 4 {
            bail!("truncated FASTQ record (expected 4 lines, got {got})");
        }
        s.seqs += 1;
        s.add_seq_line(&lines[1]); // sequence line only
    }
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn gc_and_composition() {
        let d = std::env::temp_dir().join(format!("cleaver-stats-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        let fa = d.join("s.fasta");
        File::create(&fa).unwrap().write_all(b">a\nACGT\n>b\nGGCC\n").unwrap();
        let s = stats_file(&fa).unwrap();
        assert_eq!(s.seqs, 2);
        assert_eq!(s.bases, 8);
        assert_eq!((s.a, s.c, s.g, s.t), (1, 3, 3, 1));
        assert_eq!(s.gc, 6);
        assert!((s.gc_percent() - 75.0).abs() < 1e-9);
    }
}
