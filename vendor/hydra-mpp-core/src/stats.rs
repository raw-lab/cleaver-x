//! Length / composition statistics. Writes polars DataFrames to the output dir.

#[cfg(feature = "full")]
use anyhow::Result;
#[cfg(feature = "full")]
use bio::io::fasta;
#[cfg(feature = "full")]
use polars::prelude::*;
#[cfg(feature = "full")]
use rayon::prelude::*;
#[cfg(feature = "full")]
use std::collections::HashMap;
#[cfg(feature = "full")]
use std::path::Path;

#[derive(Debug, Clone)]
pub struct BasicStats {
    pub num_seq: usize,
    pub total_bp: usize,
    pub gc_count: usize,
    pub at_count: usize,
    pub n_count: usize,
    /// Lengths in the order records appeared.
    pub lengths: Vec<usize>,
}

/// Single parallel pass over the records: lengths and base composition.
#[cfg(feature = "full")]
pub fn compute_basic(records: &[fasta::Record]) -> BasicStats {
    let per: Vec<(usize, usize, usize, usize)> = records
        .par_iter()
        .map(|r| {
            let s = r.seq();
            let (mut gc, mut at, mut n) = (0usize, 0usize, 0usize);
            for &b in s {
                match b.to_ascii_uppercase() {
                    b'G' | b'C' => gc += 1,
                    b'A' | b'T' | b'U' => at += 1,
                    b'N' => n += 1,
                    _ => {}
                }
            }
            (s.len(), gc, at, n)
        })
        .collect();

    let mut stats = BasicStats {
        num_seq: per.len(),
        total_bp: 0,
        gc_count: 0,
        at_count: 0,
        n_count: 0,
        lengths: Vec::with_capacity(per.len()),
    };
    for (len, gc, at, n) in per {
        stats.lengths.push(len);
        stats.total_bp += len;
        stats.gc_count += gc;
        stats.at_count += at;
        stats.n_count += n;
    }
    stats
}

/// N and L statistics at the 25/50/75/90 thresholds.
///
/// NOTE: The original implementation stored the running cumulative length in
/// `n50`, which is not the N50 metric. Correct definition: `N_x` is the length
/// of the shortest contig such that the sum of contigs ≥ that length covers
/// `x%` of the total assembly. Here we store `len` when the threshold is first
/// crossed, matching the standard definition.
pub fn compute_nl(lengths_desc: &[usize]) -> ([usize; 4], [usize; 4]) {
    let total: usize = lengths_desc.iter().sum();
    let t = [
        total * 25 / 100,
        total * 50 / 100,
        total * 75 / 100,
        total * 90 / 100,
    ];
    let mut n = [0usize; 4];
    let mut l = [0usize; 4];
    let mut cum = 0usize;
    for (i, &len) in lengths_desc.iter().enumerate() {
        cum += len;
        let count = i + 1;
        for k in 0..4 {
            if n[k] == 0 && cum >= t[k] {
                n[k] = len;
                l[k] = count;
            }
        }
    }
    (n, l)
}

/// Per-sequence GC percentage (ambiguous bases excluded from the denominator).
#[cfg(feature = "full")]
fn gc_percent(seq: &[u8]) -> f64 {
    let (mut gc, mut denom) = (0usize, 0usize);
    for &b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' => {
                gc += 1;
                denom += 1;
            }
            b'A' | b'T' | b'U' => denom += 1,
            _ => {}
        }
    }
    if denom == 0 {
        0.0
    } else {
        gc as f64 / denom as f64 * 100.0
    }
}

/// Compute and persist summary stats, per-sequence lengths, and a length histogram.
///
/// Outputs (under `outdir`):
/// * `summary_stats.csv`
/// * `per_sequence.csv` (id, length, gc)
/// * `length_intervals.csv` (lower_bp, upper_bp, count)
///
/// Returns a human-readable summary for stdout.
#[cfg(feature = "full")]
pub fn compute_length_stats(
    records: &[fasta::Record],
    interval: usize,
    outdir: &Path,
) -> Result<String> {
    assert!(interval > 0, "interval must be > 0");

    let basic = compute_basic(records);

    // descending sort for N/L calculation
    let mut lengths_desc = basic.lengths.clone();
    lengths_desc.sort_unstable_by(|a, b| b.cmp(a));
    let (nx, lx) = compute_nl(&lengths_desc);

    let total_bp = basic.total_bp as f64;
    let gc_pct = if total_bp > 0.0 { basic.gc_count as f64 / total_bp * 100.0 } else { 0.0 };
    let at_pct = if total_bp > 0.0 { basic.at_count as f64 / total_bp * 100.0 } else { 0.0 };
    let n_pct  = if total_bp > 0.0 { basic.n_count  as f64 / total_bp * 100.0 } else { 0.0 };
    let max_seq = lengths_desc.first().copied().unwrap_or(0);
    let min_seq = lengths_desc.last().copied().unwrap_or(0);

    // --- summary_stats.csv (long form for easy plotting) ---
    let metrics = vec![
        "num_sequences", "total_bp", "max_length", "min_length",
        "gc_percent", "at_percent", "n_percent",
        "N25", "N50", "N75", "N90",
        "L25", "L50", "L75", "L90",
    ];
    let values: Vec<f64> = vec![
        basic.num_seq as f64, basic.total_bp as f64, max_seq as f64, min_seq as f64,
        gc_pct, at_pct, n_pct,
        nx[0] as f64, nx[1] as f64, nx[2] as f64, nx[3] as f64,
        lx[0] as f64, lx[1] as f64, lx[2] as f64, lx[3] as f64,
    ];
    let mut summary = df! {
        "metric" => metrics,
        "value"  => values,
    }?;
    let path = outdir.join("summary_stats.csv");
    let mut f = std::fs::File::create(&path)?;
    CsvWriter::new(&mut f).finish(&mut summary)?;

    // --- per-sequence csv (for length & GC plots) ---
    let ids:  Vec<String> = records.iter().map(|r| r.id().to_string()).collect();
    let lens: Vec<u64>    = records.iter().map(|r| r.seq().len() as u64).collect();
    let gcs:  Vec<f64>    = records.iter().map(|r| gc_percent(r.seq())).collect();
    let mut per_seq = df! {
        "id"     => ids,
        "length" => lens,
        "gc"     => gcs,
    }?;
    let path = outdir.join("per_sequence.csv");
    let mut f = std::fs::File::create(&path)?;
    CsvWriter::new(&mut f).finish(&mut per_seq)?;

    // --- length histogram (binned) ---
    let mut bins: HashMap<usize, u64> = HashMap::new();
    for &l in &lengths_desc {
        *bins.entry(l / interval).or_insert(0) += 1;
    }
    let mut sorted: Vec<_> = bins.into_iter().collect();
    sorted.sort_by_key(|(k, _)| *k);
    let lo: Vec<u64> = sorted.iter().map(|(k, _)| (k * interval) as u64).collect();
    let hi: Vec<u64> = sorted.iter().map(|(k, _)| ((k + 1) * interval - 1) as u64).collect();
    let ct: Vec<u64> = sorted.iter().map(|(_, c)| *c).collect();
    let mut hist = df! {
        "lower_bp" => lo,
        "upper_bp" => hi,
        "count"    => ct,
    }?;
    let path = outdir.join("length_intervals.csv");
    let mut f = std::fs::File::create(&path)?;
    CsvWriter::new(&mut f).finish(&mut hist)?;

    Ok(format!(
        "num_seq = {num_seq}\n\
         total_bp = {total_bp}\n\
         GC% = {gc:.2}   AT% = {at:.2}   N% = {n:.2}\n\
         max = {max}   min = {min}\n\
         N25 = {n25}   N50 = {n50}   N75 = {n75}   N90 = {n90}\n\
         L25 = {l25}   L50 = {l50}   L75 = {l75}   L90 = {l90}",
        num_seq = basic.num_seq,
        total_bp = basic.total_bp,
        gc = gc_pct, at = at_pct, n = n_pct,
        max = max_seq, min = min_seq,
        n25 = nx[0], n50 = nx[1], n75 = nx[2], n90 = nx[3],
        l25 = lx[0], l50 = lx[1], l75 = lx[2], l90 = lx[3],
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nl_definition() {
        // lengths sorted desc: 100, 40, 30, 20, 10 => total 200
        // N50: cumulative >= 100 at 100 -> N50 = 100, L50 = 1
        let lens = vec![100, 40, 30, 20, 10];
        let (n, l) = compute_nl(&lens);
        assert_eq!(n[1], 100);
        assert_eq!(l[1], 1);
        // N90 at cumulative >= 180: 100+40+30+20 = 190 -> N90 = 20, L90 = 4
        assert_eq!(n[3], 20);
        assert_eq!(l[3], 4);
    }
}
