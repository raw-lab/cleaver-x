//! Comprehensive genome / assembly statistics for FASTA & FASTQ.
//!
//! One streaming pass collects per-sequence lengths and base composition; the
//! N/L assembly metrics (N25/N50/N75/N90 and matching L-values) are then derived
//! with the standard definition — the smallest set of longest sequences whose
//! cumulative length first reaches X % of the total. Memory stays flat in the
//! input size except for one `usize` per sequence (the length vector), which N50
//! fundamentally requires.

use std::io::BufRead;
use std::path::Path;

use anyhow::Result;

use crate::compute::{count_bytes, Counts};
use crate::formats::Format;

/// Genome statistics for a single file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GenomeStats {
    pub n_seqs: u64,
    pub total_bp: u64,
    pub min_len: u64,
    pub max_len: u64,
    pub mean_len: f64,
    pub median_len: u64,
    /// N25/N50/N75/N90.
    pub nx: [u64; 4],
    /// L25/L50/L75/L90.
    pub lx: [u64; 4],
    pub gc_percent: f64,
    pub counts: Counts,
}

impl GenomeStats {
    pub fn n50(&self) -> u64 {
        self.nx[1]
    }
    pub fn l50(&self) -> u64 {
        self.lx[1]
    }
    pub fn n90(&self) -> u64 {
        self.nx[3]
    }
    pub fn l90(&self) -> u64 {
        self.lx[3]
    }
}

/// N/L statistics over sequence lengths sorted **descending**, via the vendored
/// `rustyomestats` crate (`rustyomestats::stats::compute_nl`).
///
/// `nx[k]` is the length of the sequence at which the running total first
/// reaches `frac[k]` of the sum (N25/N50/N75/N90); `lx[k]` is how many
/// sequences that took.
fn compute_nl(lengths_desc: &[usize]) -> ([u64; 4], [u64; 4]) {
    let (nx, lx) = rustyomestats::stats::compute_nl(lengths_desc);
    let to_u64 = |a: [usize; 4]| [a[0] as u64, a[1] as u64, a[2] as u64, a[3] as u64];
    (to_u64(nx), to_u64(lx))
}

#[inline]
fn seq_len(line: &[u8]) -> usize {
    line.iter().filter(|&&b| b != b'\n' && b != b'\r').count()
}

/// Stream `path` and compute its genome statistics.
///
/// If `precomputed` is supplied (e.g. base composition already tallied by the
/// GPU kernel), the pass only measures per-sequence lengths and reuses those
/// counts; otherwise composition is tallied here in the same pass.
pub fn genome_stats_file(
    path: &Path,
    fmt: Format,
    precomputed: Option<Counts>,
) -> Result<GenomeStats> {
    let count_here = precomputed.is_none();
    let mut counts = precomputed.unwrap_or_default();

    let (mut reader, _gz) = crate::open_reader(path)?;
    let mut lengths: Vec<usize> = Vec::new();

    let mut cur: usize = 0;
    let mut in_record = false;
    let mut lineno: u64 = 0;
    let mut line: Vec<u8> = Vec::with_capacity(256);

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        match fmt {
            Format::Fasta => {
                if line.first() == Some(&b'>') {
                    if in_record {
                        lengths.push(cur);
                    }
                    cur = 0;
                    in_record = true;
                } else {
                    cur += seq_len(&line);
                    if count_here {
                        counts.add(&count_bytes(&line));
                    }
                }
            }
            Format::Fastq => {
                if lineno % 4 == 1 {
                    lengths.push(seq_len(&line));
                    if count_here {
                        counts.add(&count_bytes(&line));
                    }
                }
            }
            _ => {
                cur += seq_len(&line);
                in_record = true;
                if count_here {
                    counts.add(&count_bytes(&line));
                }
            }
        }
        lineno += 1;
    }
    if in_record && !matches!(fmt, Format::Fastq) {
        lengths.push(cur);
    }

    Ok(summarize(lengths, counts))
}

fn summarize(mut lengths: Vec<usize>, counts: Counts) -> GenomeStats {
    let n_seqs = lengths.len() as u64;
    if n_seqs == 0 {
        return GenomeStats { counts, ..Default::default() };
    }
    let total_bp: u64 = lengths.iter().map(|&l| l as u64).sum();
    let min_len = *lengths.iter().min().unwrap() as u64;
    let max_len = *lengths.iter().max().unwrap() as u64;
    let mean_len = total_bp as f64 / n_seqs as f64;

    lengths.sort_unstable(); // ascending for the median
    let m = lengths.len();
    let median_len = if m % 2 == 1 {
        lengths[m / 2] as u64
    } else {
        (lengths[m / 2 - 1] as u64 + lengths[m / 2] as u64) / 2
    };

    let mut desc = lengths;
    desc.reverse();
    let (nx, lx) = compute_nl(&desc);

    GenomeStats {
        n_seqs,
        total_bp,
        min_len,
        max_len,
        mean_len,
        median_len,
        nx,
        lx,
        gc_percent: counts.gc() * 100.0,
        counts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    fn tmp(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("cleaver-genome-{tag}-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        d
    }

    #[test]
    fn fasta_genome_stats_and_n50() {
        let d = tmp("fa");
        let p = d.join("g.fasta");
        let mut data = Vec::new();
        data.extend_from_slice(b">a\n");
        for _ in 0..2 {
            data.extend_from_slice(&vec![b'G'; 50]); // 100 bp, all GC
            data.push(b'\n');
        }
        data.extend_from_slice(b">b\n");
        data.extend_from_slice(&vec![b'A'; 40]); // 40 bp
        data.push(b'\n');
        data.extend_from_slice(b">c\n");
        data.extend_from_slice(&vec![b'C'; 10]); // 10 bp, GC
        data.push(b'\n');
        File::create(&p).unwrap().write_all(&data).unwrap();

        let s = genome_stats_file(&p, Format::Fasta, None).unwrap();
        assert_eq!((s.n_seqs, s.total_bp, s.min_len, s.max_len, s.median_len), (3, 150, 10, 100, 40));
        assert_eq!((s.n50(), s.l50()), (100, 1)); // 50% = 75 -> first seq
        assert_eq!((s.n90(), s.l90()), (40, 2)); // 90% = 135 -> 100+40
        assert!((s.gc_percent - (110.0 / 150.0 * 100.0)).abs() < 1e-6);
    }

    #[test]
    fn fastq_lengths_per_record() {
        let d = tmp("fq");
        let p = d.join("r.fastq");
        File::create(&p).unwrap().write_all(b"@r1\nACGT\n+\nIIII\n@r2\nGGCCAA\n+\nJJJJJJ\n").unwrap();
        let s = genome_stats_file(&p, Format::Fastq, None).unwrap();
        assert_eq!((s.n_seqs, s.total_bp, s.min_len, s.max_len), (2, 10, 4, 6));
    }

    #[test]
    fn precomputed_counts_are_reused() {
        let d = tmp("pre");
        let p = d.join("g.fasta");
        File::create(&p).unwrap().write_all(b">a\nACGTACGT\n").unwrap();
        let fake = Counts { a: 0, c: 8, g: 0, t: 0, n: 0, other: 0 };
        let s = genome_stats_file(&p, Format::Fasta, Some(fake)).unwrap();
        assert_eq!(s.total_bp, 8);
        assert!((s.gc_percent - 100.0).abs() < 1e-9);
    }
}
