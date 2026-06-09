//! Assembly statistics after Castro et al. (2016):
//! **U50**, **UL50**, **UG50**, **ULG50**, **UG50%**, plus the standard
//! **N50**, **L50**, **NG50**, **LG50** for comparison.
//!
//! Input:
//! * reference FASTA (single-record; only the first record's length is used —
//!   matches the Castro script's design, typical for viral/bacterial refs)
//! * sorted BED of contigs mapped to that reference (0-based half-open,
//!   i.e. `chromEnd` is exclusive; standard UCSC BED).
//!
//! Output files under `outdir`:
//! * `u50_summary.csv`         — every metric in long form
//! * `u50_contigs.csv`         — per-contig start, stop, orig_length, unique_bp
//! * `u50_gap_intervals.csv`   — uncovered stretches [start, end)
//! * `u50_overlap_intervals.csv` — stretches covered by ≥ 2 contigs (with depth)
//!
//! ## Definitions (user-facing; see README for the formal spec)
//!
//! Let the mapped contigs, sorted longest → shortest on their **original**
//! lengths, be `c_1, c_2, …, c_n`. "Unique" contigs `c'_k` are the same
//! contigs after greedy masking: every position claimed by an earlier
//! (longer) contig is removed from a later one; non-overlapping portions are
//! preserved. Contigs reduced to zero unique bp drop out.
//!
//! * **N50** — shortest `c_k` whose cumulative length first reaches 50 % of
//!   the total original contig length.
//! * **L50** — the index `k` at which N50 is crossed.
//! * **NG50 / LG50** — same as N50 / L50 but with the cutoff set to 50 % of
//!   the reference genome length (instead of total contig length).
//! * **U50 / UL50** — N50-style on the unique-bp lengths `c'_k`, cutoff at
//!   50 % of total unique bp.
//! * **UG50 / ULG50** — same as U50 / UL50 but with the cutoff at 50 % of
//!   the reference genome length.
//! * **UG50%** — `100 * UG50 / reference_length`, i.e. the fraction of the
//!   reference covered by the UG50 contig.

use anyhow::{bail, Context, Result};
use bio::io::fasta;
use polars::prelude::*;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// --------------------------------------------------------------------------- //
// data model                                                                   //
// --------------------------------------------------------------------------- //

#[derive(Debug, Clone)]
pub struct Contig {
    pub start:       u64, // 0-based inclusive
    pub stop:        u64, // 0-based exclusive (BED half-open)
    pub orig_length: u64, // stop - start
    pub unique_bp:   u64, // populated after greedy masking
}

#[derive(Debug, Clone)]
pub struct U50Result {
    pub ref_name:             String,
    pub ref_length:           u64,
    pub num_contigs_in_bed:   usize,
    pub num_contigs_kept:     usize, // after dedup
    pub num_contigs_unique:   usize, // contributing ≥ 1 unique bp
    pub total_orig_length:    u64,
    pub total_unique_length:  u64,
    pub gap_bp:               u64,
    pub overlap_bp:           u64,

    // N50 / L50: original lengths vs total contig length
    pub n50: u64, pub l50: usize,
    // NG50 / LG50: original lengths vs reference length
    pub ng50: u64, pub lg50: usize,
    // U50 / UL50: unique-bp lengths vs total unique length
    pub u50: u64, pub ul50: usize,
    // UG50 / ULG50: unique-bp lengths vs reference length
    pub ug50: u64, pub ulg50: usize,
    // UG50% : 100 * UG50 / ref_length
    pub ug50_pct: f64,
}

// --------------------------------------------------------------------------- //
// I/O                                                                          //
// --------------------------------------------------------------------------- //

/// Return `(first_record_id, length)` from the reference FASTA.
pub fn load_reference_length(path: &Path) -> Result<(String, u64)> {
    let reader = fasta::Reader::from_file(path)
        .with_context(|| format!("opening reference FASTA {:?}", path))?;
    for r in reader.records() {
        let r = r?;
        return Ok((r.id().to_string(), r.seq().len() as u64));
    }
    bail!("reference FASTA {:?} contains no records", path);
}

/// Load a BED-3+ file. Lines starting with `#`, `track`, or `browser` are
/// ignored. Only columns 1-3 (chrom/start/stop) are consumed.
pub fn load_bed(path: &Path) -> Result<Vec<Contig>> {
    let f = File::open(path).with_context(|| format!("opening BED {:?}", path))?;
    let mut v = Vec::new();
    for (ln, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let t = line.trim();
        if t.is_empty()
            || t.starts_with('#')
            || t.starts_with("track")
            || t.starts_with("browser")
        {
            continue;
        }
        let parts: Vec<&str> = t.split('\t').collect();
        if parts.len() < 3 {
            bail!("BED line {} has fewer than 3 columns: {:?}", ln + 1, t);
        }
        let start: u64 = parts[1]
            .parse()
            .with_context(|| format!("BED line {}: bad start {:?}", ln + 1, parts[1]))?;
        let stop: u64 = parts[2]
            .parse()
            .with_context(|| format!("BED line {}: bad stop {:?}", ln + 1, parts[2]))?;
        if stop <= start {
            continue; // empty interval — skip silently
        }
        v.push(Contig { start, stop, orig_length: stop - start, unique_bp: 0 });
    }
    Ok(v)
}

/// Drop exact `(start, stop)` duplicates. This corrects the original
/// Castro script's `"start not in starts AND stop not in stops"` test,
/// which incidentally dropped non-duplicate contigs sharing a boundary.
fn dedup_contigs(contigs: Vec<Contig>) -> Vec<Contig> {
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(contigs.len());
    for c in contigs {
        if seen.insert((c.start, c.stop)) {
            out.push(c);
        }
    }
    out
}

// --------------------------------------------------------------------------- //
// masking                                                                      //
// --------------------------------------------------------------------------- //

/// Greedy mask: sort contigs by original length desc and let each claim every
/// position it covers that isn't already taken. Sets `unique_bp` in place and
/// returns the mask (bit per reference position: 1 = covered, 0 = gap).
///
/// Overlap depth is tracked in a separate per-position `u16` vector so we can
/// emit overlap intervals later; saturating-adds at u16::MAX if a position is
/// absurdly covered.
///
/// Memory: ~3 bytes per reference bp. For viral/bacterial genomes where U50
/// was designed this is trivial (a few MB). For mammalian references
/// consider an interval-tree approach instead.
pub fn mask_and_count(contigs: &mut [Contig], ref_length: u64) -> (Vec<u8>, Vec<u16>) {
    let mut order: Vec<usize> = (0..contigs.len()).collect();
    // longest first, ties broken deterministically by (start, stop)
    order.sort_by(|&i, &j| {
        contigs[j].orig_length
            .cmp(&contigs[i].orig_length)
            .then(contigs[i].start.cmp(&contigs[j].start))
            .then(contigs[i].stop.cmp(&contigs[j].stop))
    });

    let n = ref_length as usize;
    let mut mask  = vec![0u8;  n];
    let mut depth = vec![0u16; n];

    for &idx in &order {
        let c = &contigs[idx];
        let a = (c.start as usize).min(n);
        let b = (c.stop  as usize).min(n);
        let mut unique = 0u64;
        for pos in a..b {
            if mask[pos] == 0 {
                mask[pos] = 1;
                unique += 1;
            }
            depth[pos] = depth[pos].saturating_add(1);
        }
        contigs[idx].unique_bp = unique;
    }
    (mask, depth)
}

// --------------------------------------------------------------------------- //
// N-stat engine                                                                //
// --------------------------------------------------------------------------- //

/// Given contig lengths sorted descending and a reference total, return
/// `(Nx, Lx)` where Nx is the shortest contig whose cumulative length first
/// reaches `percent %` of `reference_total`. If cumulative never reaches the
/// cutoff (e.g. assembly smaller than half the reference), returns `(0, 0)`,
/// matching the Castro script's "running sum never reaches median" case.
///
/// Integer arithmetic only — no rounding ambiguity.
pub fn n_stat(lengths_desc: &[u64], reference_total: u64, percent: u32) -> (u64, usize) {
    if reference_total == 0 || percent == 0 {
        return (0, 0);
    }
    let mut cum: u128 = 0;
    let thresh = reference_total as u128 * percent as u128;
    for (i, &len) in lengths_desc.iter().enumerate() {
        cum += len as u128;
        if cum * 100 >= thresh {
            return (len, i + 1);
        }
    }
    (0, 0)
}

// --------------------------------------------------------------------------- //
// interval extraction                                                          //
// --------------------------------------------------------------------------- //

fn mask_intervals(mask: &[u8], want: u8) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < mask.len() {
        if mask[i] == want {
            let start = i;
            while i < mask.len() && mask[i] == want {
                i += 1;
            }
            out.push((start as u64, i as u64));
        } else {
            i += 1;
        }
    }
    out
}

/// Overlap runs with max depth per run. An "overlap" position is one covered
/// by ≥ 2 contigs (depth ≥ 2).
fn overlap_intervals(depth: &[u16]) -> Vec<(u64, u64, u16)> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < depth.len() {
        if depth[i] >= 2 {
            let start = i;
            let mut max_d = depth[i];
            while i < depth.len() && depth[i] >= 2 {
                max_d = max_d.max(depth[i]);
                i += 1;
            }
            out.push((start as u64, i as u64, max_d));
        } else {
            i += 1;
        }
    }
    out
}

// --------------------------------------------------------------------------- //
// main entry                                                                   //
// --------------------------------------------------------------------------- //

pub fn compute_u50(ref_fasta: &Path, bed: &Path, outdir: &Path) -> Result<U50Result> {
    std::fs::create_dir_all(outdir)?;

    let (ref_name, ref_length) = load_reference_length(ref_fasta)?;
    let raw = load_bed(bed)?;
    let num_total = raw.len();
    let mut contigs = dedup_contigs(raw);
    let num_kept = contigs.len();

    let (mask, depth) = mask_and_count(&mut contigs, ref_length);

    // --- lengths sorted desc in each flavor ---
    let mut orig_desc: Vec<u64> = contigs.iter().map(|c| c.orig_length).collect();
    orig_desc.sort_unstable_by(|a, b| b.cmp(a));
    let mut uniq_desc: Vec<u64> = contigs
        .iter()
        .map(|c| c.unique_bp)
        .filter(|&x| x > 0)
        .collect();
    uniq_desc.sort_unstable_by(|a, b| b.cmp(a));

    let total_orig: u64   = orig_desc.iter().sum();
    let total_unique: u64 = uniq_desc.iter().sum();
    let gap_bp = mask.iter().filter(|&&b| b == 0).count() as u64;
    let overlap_bp = depth.iter().filter(|&&d| d >= 2).count() as u64;

    let (n50, l50)    = n_stat(&orig_desc, total_orig,   50);
    let (ng50, lg50)  = n_stat(&orig_desc, ref_length,   50);
    let (u50, ul50)   = n_stat(&uniq_desc, total_unique, 50);
    let (ug50, ulg50) = n_stat(&uniq_desc, ref_length,   50);

    let ug50_pct = if ref_length > 0 {
        100.0 * ug50 as f64 / ref_length as f64
    } else {
        0.0
    };

    // --- write u50_contigs.csv (sorted unique_bp desc) ---
    {
        let mut cs = contigs.clone();
        cs.sort_by(|a, b| b.unique_bp.cmp(&a.unique_bp)
            .then(b.orig_length.cmp(&a.orig_length)));
        let start:  Vec<u64> = cs.iter().map(|c| c.start).collect();
        let stop:   Vec<u64> = cs.iter().map(|c| c.stop).collect();
        let olen:   Vec<u64> = cs.iter().map(|c| c.orig_length).collect();
        let uniq:   Vec<u64> = cs.iter().map(|c| c.unique_bp).collect();
        let mut df = df! {
            "start"       => start,
            "stop"        => stop,
            "orig_length" => olen,
            "unique_bp"   => uniq,
        }?;
        let mut f = File::create(outdir.join("u50_contigs.csv"))?;
        CsvWriter::new(&mut f).finish(&mut df)?;
    }

    // --- gap intervals ---
    {
        let gi = mask_intervals(&mask, 0);
        let s: Vec<u64> = gi.iter().map(|(a, _)| *a).collect();
        let e: Vec<u64> = gi.iter().map(|(_, b)| *b).collect();
        let l: Vec<u64> = gi.iter().map(|(a, b)| b - a).collect();
        let mut df = df! { "start" => s, "end" => e, "length" => l }?;
        let mut f = File::create(outdir.join("u50_gap_intervals.csv"))?;
        CsvWriter::new(&mut f).finish(&mut df)?;
    }

    // --- overlap intervals ---
    {
        let oi = overlap_intervals(&depth);
        let s: Vec<u64> = oi.iter().map(|(a, _, _)| *a).collect();
        let e: Vec<u64> = oi.iter().map(|(_, b, _)| *b).collect();
        let l: Vec<u64> = oi.iter().map(|(a, b, _)| b - a).collect();
        let d: Vec<u32> = oi.iter().map(|(_, _, d)| *d as u32).collect();
        let mut df = df! {
            "start"     => s,
            "end"       => e,
            "length"    => l,
            "max_depth" => d,
        }?;
        let mut f = File::create(outdir.join("u50_overlap_intervals.csv"))?;
        CsvWriter::new(&mut f).finish(&mut df)?;
    }

    // --- summary ---
    let num_contigs_unique = contigs.iter().filter(|c| c.unique_bp > 0).count();
    let metrics: Vec<String> = vec![
        "ref_name",
        "ref_length",
        "num_contigs_in_bed",
        "num_contigs_kept_after_dedup",
        "num_contigs_contributing_unique_bp",
        "total_original_contig_length",
        "total_unique_contig_length",
        "gap_bp",
        "overlap_bp",
        "N50", "L50",
        "NG50", "LG50",
        "U50", "UL50",
        "UG50", "ULG50",
        "UG50_percent",
    ].into_iter().map(String::from).collect();
    let values: Vec<String> = vec![
        ref_name.clone(),
        ref_length.to_string(),
        num_total.to_string(),
        num_kept.to_string(),
        num_contigs_unique.to_string(),
        total_orig.to_string(),
        total_unique.to_string(),
        gap_bp.to_string(),
        overlap_bp.to_string(),
        n50.to_string(),  l50.to_string(),
        ng50.to_string(), lg50.to_string(),
        u50.to_string(),  ul50.to_string(),
        ug50.to_string(), ulg50.to_string(),
        format!("{:.4}", ug50_pct),
    ];
    let mut df = df! { "metric" => metrics, "value" => values }?;
    let mut f = File::create(outdir.join("u50_summary.csv"))?;
    CsvWriter::new(&mut f).finish(&mut df)?;

    Ok(U50Result {
        ref_name, ref_length,
        num_contigs_in_bed: num_total,
        num_contigs_kept:   num_kept,
        num_contigs_unique,
        total_orig_length:    total_orig,
        total_unique_length:  total_unique,
        gap_bp,
        overlap_bp,
        n50, l50, ng50, lg50, u50, ul50, ug50, ulg50, ug50_pct,
    })
}

// --------------------------------------------------------------------------- //
// tests                                                                        //
// --------------------------------------------------------------------------- //

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn n_stat_basic() {
        // lengths: 100, 40, 30, 20, 10 → total 200
        let lens = vec![100u64, 40, 30, 20, 10];
        // N50 over total=200: cutoff 100; cum at step 1 = 100 → N50 = 100, L50 = 1
        assert_eq!(n_stat(&lens, 200, 50), (100, 1));
        // N90: cutoff 180; cum: 100,140,170,190 → N90 = 20, L90 = 4
        assert_eq!(n_stat(&lens, 200, 90), (20, 4));
    }

    #[test]
    fn n_stat_below_cutoff_returns_zero() {
        // total contig length 50, reference 1000 → NG50 cutoff = 500,
        // cumulative never reaches → (0, 0)
        let lens = vec![30u64, 20];
        assert_eq!(n_stat(&lens, 1000, 50), (0, 0));
    }

    #[test]
    fn masking_greedy() {
        // Ref length 100.
        // Contig A: 0..60 (longest, 60 bp)
        // Contig B: 40..80 (40 bp, overlaps A at 40..60)
        // Contig C: 80..100 (20 bp, no overlap)
        let mut contigs = vec![
            Contig { start: 0,  stop: 60,  orig_length: 60, unique_bp: 0 },
            Contig { start: 40, stop: 80,  orig_length: 40, unique_bp: 0 },
            Contig { start: 80, stop: 100, orig_length: 20, unique_bp: 0 },
        ];
        let (mask, depth) = mask_and_count(&mut contigs, 100);
        assert_eq!(contigs[0].unique_bp, 60); // A claims 0..60
        assert_eq!(contigs[1].unique_bp, 20); // B only gets 60..80
        assert_eq!(contigs[2].unique_bp, 20); // C gets 80..100
        assert_eq!(mask.iter().filter(|&&b| b == 0).count(), 0); // full coverage
        // overlap region is 40..60 (20 bp at depth 2)
        assert_eq!(depth.iter().filter(|&&d| d >= 2).count(), 20);
    }

    #[test]
    fn u50_whole_pipeline_numbers() {
        // 3 contigs as above, ref=100
        let mut contigs = vec![
            Contig { start: 0,  stop: 60,  orig_length: 60, unique_bp: 0 },
            Contig { start: 40, stop: 80,  orig_length: 40, unique_bp: 0 },
            Contig { start: 80, stop: 100, orig_length: 20, unique_bp: 0 },
        ];
        let _ = mask_and_count(&mut contigs, 100);

        let mut orig: Vec<u64> = contigs.iter().map(|c| c.orig_length).collect();
        orig.sort_unstable_by(|a, b| b.cmp(a));
        // N50: total 120, cutoff 60 → first contig (60) reaches → N50=60, L50=1
        assert_eq!(n_stat(&orig, 120, 50), (60, 1));

        let mut uniq: Vec<u64> = contigs.iter().map(|c| c.unique_bp).collect();
        uniq.sort_unstable_by(|a, b| b.cmp(a));
        // uniq desc: [60, 20, 20], total 100
        // U50: cutoff 50 → first contig (60) reaches → U50=60, UL50=1
        assert_eq!(n_stat(&uniq, 100, 50), (60, 1));
        // UG50: cutoff 50 (ref 100) → same → UG50=60, ULG50=1, UG50%=60%
        assert_eq!(n_stat(&uniq, 100, 50), (60, 1));
    }

    #[test]
    fn dedup_drops_exact_duplicates_only() {
        let cs = vec![
            Contig { start: 0,  stop: 10, orig_length: 10, unique_bp: 0 },
            Contig { start: 0,  stop: 10, orig_length: 10, unique_bp: 0 }, // dup
            Contig { start: 0,  stop: 20, orig_length: 20, unique_bp: 0 }, // same start, kept
            Contig { start: 5,  stop: 10, orig_length:  5, unique_bp: 0 }, // same stop, kept
        ];
        let d = dedup_contigs(cs);
        assert_eq!(d.len(), 3);
    }
}
