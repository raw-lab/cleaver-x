//! Codon density:
//! * **Absolute** — all 6 reading frames of every input sequence.
//! * **Predicted** — counted only inside ORFs predicted by FragGeneScanRs
//!   (we read its `.ffn` file; see [`crate::fgs`]).
//!
//! All outputs are written as polars DataFrames (csv), ready for seaborn
//! plotting via the Python driver.

use anyhow::Result;
use bio::io::fasta;
use polars::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::Path;

/// Canonical 64 codons in AACG/ACGT lexicographic order.
pub const CODONS: [&str; 64] = [
    "AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT",
    "AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
    "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT",
    "CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
    "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT",
    "GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
    "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT",
    "TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT",
];

#[inline]
fn nt_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

#[inline]
fn codon_idx(c: &[u8]) -> Option<usize> {
    if c.len() != 3 { return None; }
    Some(nt_idx(c[0])? * 16 + nt_idx(c[1])? * 4 + nt_idx(c[2])?)
}

#[inline]
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A'        => b'T',
        b'T' | b'U' => b'A',
        b'G'        => b'C',
        b'C'        => b'G',
        x           => x,
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Count codons in a single reading frame starting at `offset`.
/// Codons containing a non-ACGT base contribute to nothing (silently skipped).
fn count_in_frame(seq: &[u8], offset: usize) -> [u64; 64] {
    let mut counts = [0u64; 64];
    let mut i = offset;
    while i + 3 <= seq.len() {
        if let Some(ix) = codon_idx(&seq[i..i + 3]) {
            counts[ix] += 1;
        }
        i += 3;
    }
    counts
}

/// All 6 frames: +1/+2/+3 on the forward strand, -1/-2/-3 on the reverse complement.
pub fn six_frame_codon_counts(seq: &[u8]) -> Vec<(String, [u64; 64])> {
    let rc = reverse_complement(seq);
    let mut out = Vec::with_capacity(6);
    for (label, off) in [("+1", 0usize), ("+2", 1), ("+3", 2)] {
        out.push((label.to_string(), count_in_frame(seq, off)));
    }
    for (label, off) in [("-1", 0usize), ("-2", 1), ("-3", 2)] {
        out.push((label.to_string(), count_in_frame(&rc, off)));
    }
    out
}

/// Public helper: count codons in a pre-oriented ORF (starts in its own frame +1).
pub fn count_codons_single_frame(seq: &[u8]) -> [u64; 64] {
    count_in_frame(seq, 0)
}

/// Write `codon_absolute.csv` (per-record × per-frame × 64 codons) and the
/// aggregated `codon_absolute_aggregate.csv`.
pub fn compute_absolute_codon_density(
    records: &[fasta::Record],
    outdir: &Path,
) -> Result<()> {
    // (id, frame, counts[64]) rows, in parallel
    let rows: Vec<(String, String, [u64; 64])> = records
        .par_iter()
        .flat_map_iter(|r| {
            let id = r.id().to_string();
            six_frame_codon_counts(r.seq())
                .into_iter()
                .map(move |(f, c)| (id.clone(), f, c))
        })
        .collect();

    // long-form dataframe
    let n = rows.len() * 64;
    let mut ids       = Vec::with_capacity(n);
    let mut frames    = Vec::with_capacity(n);
    let mut codons    = Vec::with_capacity(n);
    let mut counts    = Vec::with_capacity(n);
    let mut densities = Vec::with_capacity(n);
    let mut agg: HashMap<&'static str, u64> = CODONS.iter().map(|c| (*c, 0u64)).collect();

    for (id, frame, cs) in &rows {
        let total: u64 = cs.iter().sum();
        for (k, &c) in cs.iter().enumerate() {
            ids.push(id.clone());
            frames.push(frame.clone());
            codons.push(CODONS[k].to_string());
            counts.push(c);
            densities.push(if total == 0 { 0.0 } else { c as f64 / total as f64 });
            *agg.get_mut(CODONS[k]).unwrap() += c;
        }
    }

    let mut df = df! {
        "id"      => ids,
        "frame"   => frames,
        "codon"   => codons,
        "count"   => counts,
        "density" => densities,
    }?;
    let mut f = std::fs::File::create(outdir.join("codon_absolute.csv"))?;
    CsvWriter::new(&mut f).finish(&mut df)?;

    // aggregate csv (total across all records and all frames)
    write_aggregate(&agg, outdir.join("codon_absolute_aggregate.csv"))?;
    Ok(())
}

/// Read FGS `.ffn` predicted ORFs and write `codon_predicted.csv` +
/// `codon_predicted_aggregate.csv`. ORFs are assumed to start in their own
/// frame +1 (this is how FragGeneScanRs writes them).
pub fn compute_predicted_codon_density(ffn_path: &Path, outdir: &Path) -> Result<()> {
    let reader = fasta::Reader::from_file(ffn_path)?;
    let records: Vec<fasta::Record> = reader.records().filter_map(|r| r.ok()).collect();

    let rows: Vec<(String, [u64; 64])> = records
        .par_iter()
        .map(|r| (r.id().to_string(), count_codons_single_frame(r.seq())))
        .collect();

    let n = rows.len() * 64;
    let mut ids       = Vec::with_capacity(n);
    let mut codons    = Vec::with_capacity(n);
    let mut counts    = Vec::with_capacity(n);
    let mut densities = Vec::with_capacity(n);
    let mut agg: HashMap<&'static str, u64> = CODONS.iter().map(|c| (*c, 0u64)).collect();

    for (id, cs) in &rows {
        let total: u64 = cs.iter().sum();
        for (k, &c) in cs.iter().enumerate() {
            ids.push(id.clone());
            codons.push(CODONS[k].to_string());
            counts.push(c);
            densities.push(if total == 0 { 0.0 } else { c as f64 / total as f64 });
            *agg.get_mut(CODONS[k]).unwrap() += c;
        }
    }

    let mut df = df! {
        "orf_id"  => ids,
        "codon"   => codons,
        "count"   => counts,
        "density" => densities,
    }?;
    let mut f = std::fs::File::create(outdir.join("codon_predicted.csv"))?;
    CsvWriter::new(&mut f).finish(&mut df)?;

    write_aggregate(&agg, outdir.join("codon_predicted_aggregate.csv"))?;
    Ok(())
}

fn write_aggregate<P: AsRef<Path>>(
    agg: &HashMap<&'static str, u64>,
    path: P,
) -> Result<()> {
    let mut keys: Vec<&&str> = agg.keys().collect();
    keys.sort();
    let codons_v: Vec<String> = keys.iter().map(|k| k.to_string()).collect();
    let counts_v: Vec<u64>    = keys.iter().map(|k| *agg.get(**k).unwrap_or(&0)).collect();
    let total: u64 = counts_v.iter().sum();
    let dens_v:   Vec<f64> = counts_v
        .iter()
        .map(|c| if total == 0 { 0.0 } else { *c as f64 / total as f64 })
        .collect();
    let mut df = df! {
        "codon"   => codons_v,
        "count"   => counts_v,
        "density" => dens_v,
    }?;
    let mut f = std::fs::File::create(path)?;
    CsvWriter::new(&mut f).finish(&mut df)?;
    Ok(())
}

/// Join absolute and predicted aggregates and compute a per-codon enrichment
/// (`density_predicted / density_absolute`). Values > 1 indicate codons
/// overrepresented in coding regions relative to 6-frame background.
pub fn write_codon_comparison(outdir: &Path) -> Result<()> {
    let abs_path  = outdir.join("codon_absolute_aggregate.csv");
    let pred_path = outdir.join("codon_predicted_aggregate.csv");
    if !abs_path.exists() || !pred_path.exists() {
        return Ok(());
    }

    let abs_lf = CsvReadOptions::default()
        .try_into_reader_with_file_path(Some(abs_path))?
        .finish()?
        .lazy()
        .select([
            col("codon"),
            col("density").alias("density_absolute"),
            col("count").alias("count_absolute"),
        ]);

    let pred_lf = CsvReadOptions::default()
        .try_into_reader_with_file_path(Some(pred_path))?
        .finish()?
        .lazy()
        .select([
            col("codon"),
            col("density").alias("density_predicted"),
            col("count").alias("count_predicted"),
        ]);

    let joined = abs_lf
        .join(
            pred_lf,
            [col("codon")],
            [col("codon")],
            JoinArgs::new(JoinType::Inner),
        )
        .with_column(
            (col("density_predicted") / (col("density_absolute") + lit(1e-12)))
                .alias("enrichment_pred_over_abs"),
        );

    let mut out = joined.collect()?;
    let mut f = std::fs::File::create(outdir.join("codon_comparison.csv"))?;
    CsvWriter::new(&mut f).finish(&mut out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rc_simple() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
    }

    #[test]
    fn frame_plus1_counts() {
        // ATG CAT TAA -> 3 codons, each appearing once
        let c = count_in_frame(b"ATGCATTAA", 0);
        let atg = codon_idx(b"ATG").unwrap();
        let cat = codon_idx(b"CAT").unwrap();
        let taa = codon_idx(b"TAA").unwrap();
        assert_eq!(c[atg], 1);
        assert_eq!(c[cat], 1);
        assert_eq!(c[taa], 1);
    }

    #[test]
    fn six_frames_total_codons() {
        // length-18 sequence => +1/+2/+3 give 6/5/5 codons; same on RC.
        let c = six_frame_codon_counts(&vec![b'A'; 18]);
        assert_eq!(c.len(), 6);
        let sum: u64 = c[0].1.iter().sum();
        assert_eq!(sum, 6);
    }
}
