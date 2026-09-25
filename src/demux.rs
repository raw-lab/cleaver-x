//! A pure-Rust read demultiplexer in the spirit of
//! [barbell](https://github.com/rickbeeloo/barbell): each read's two ends are
//! scanned for the best-matching barcode using a fitting (semi-global) edit-
//! distance alignment, the read is assigned to a barcode when the match is
//! confident (identity threshold + margin over the runner-up), the barcode is
//! trimmed, and reads are split into per-barcode FASTQ files.
//!
//! This is a faithful re-implementation of barbell's annotate→trim→split idea in
//! safe, dependency-light Rust (barbell itself uses SIMD alignment via `sassy`
//! and a bundled kit database). Barcodes are supplied as a FASTA; the canonical
//! ONT layout (forward barcode at the 5' end, reverse-complement at the 3' end)
//! is searched at both ends.

use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::Path;

use anyhow::{Context, Result};

#[derive(Debug, Clone)]
pub struct DemuxParams {
    /// Bases at each read end to scan for a barcode.
    pub window: usize,
    /// Minimum identity (0..1) of a barcode match to accept it (`--min-score`).
    pub min_score: f64,
    /// Minimum identity margin of the best over the second-best barcode
    /// (`--min-score-diff`).
    pub min_score_diff: f64,
    /// Hard cap on edit distance for a match (`--max-errors`); overrides
    /// `min_score` when set.
    pub max_errors: Option<usize>,
    /// Trim the matched barcode region from the read (else annotate only).
    pub trim: bool,
    /// Also search the 3' end (reverse-complement barcode).
    pub both_ends: bool,
}

impl Default for DemuxParams {
    fn default() -> Self {
        DemuxParams {
            window: 150,
            min_score: 0.75,
            min_score_diff: 0.05,
            max_errors: None,
            trim: true,
            both_ends: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Barcode {
    pub name: String,
    pub seq: Vec<u8>,
}

#[derive(Debug, Clone, Default)]
pub struct DemuxStats {
    pub total: u64,
    pub classified: u64,
    pub unclassified: u64,
    /// Reads per barcode name (sorted).
    pub per_barcode: BTreeMap<String, u64>,
}

fn comp(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

/// Read barcodes from a FASTA file.
pub fn read_barcodes(path: &Path) -> Result<Vec<Barcode>> {
    let data = std::fs::read_to_string(path)
        .with_context(|| format!("reading barcodes '{}'", path.display()))?;
    let mut out = Vec::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    for line in data.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !seq.is_empty() {
                out.push(Barcode { name: std::mem::take(&mut name), seq: std::mem::take(&mut seq) });
            }
            name = h.split_whitespace().next().unwrap_or("barcode").to_string();
        } else {
            seq.extend(line.trim().bytes());
        }
    }
    if !seq.is_empty() {
        out.push(Barcode { name, seq });
    }
    if out.is_empty() {
        anyhow::bail!("no barcodes found in '{}'", path.display());
    }
    Ok(out)
}

/// Fitting (semi-global) alignment of `pat` into `text`: the whole pattern is
/// aligned but it may sit anywhere in `text` with free flanking gaps. Returns
/// `(min_edits, start, end)` where `text[start..end]` is the matched span.
fn fit_align(text: &[u8], pat: &[u8]) -> (usize, usize, usize) {
    let m = pat.len();
    let n = text.len();
    if m == 0 {
        return (0, 0, 0);
    }
    if n == 0 {
        return (m, 0, 0);
    }
    // dp[j] over text positions for the current pattern row; track the text
    // start column that each cell descends from (for the matched span).
    let mut prev = vec![0usize; n + 1]; // row 0: free leading gap in text
    let mut prev_start = (0..=n).collect::<Vec<usize>>();
    let mut cur = vec![0usize; n + 1];
    let mut cur_start = vec![0usize; n + 1];

    for i in 1..=m {
        cur[0] = i; // pattern prefix vs empty text
        cur_start[0] = 0;
        for j in 1..=n {
            let cost = if pat[i - 1].to_ascii_uppercase() == text[j - 1].to_ascii_uppercase() { 0 } else { 1 };
            // diagonal
            let mut best = prev[j - 1] + cost;
            let mut bstart = prev_start[j - 1];
            // deletion in text (gap in text) = pattern char unmatched
            if prev[j] + 1 < best {
                best = prev[j] + 1;
                bstart = prev_start[j];
            }
            // insertion in text (skip a text char)
            if cur[j - 1] + 1 < best {
                best = cur[j - 1] + 1;
                bstart = cur_start[j - 1];
            }
            cur[j] = best;
            cur_start[j] = bstart;
        }
        std::mem::swap(&mut prev, &mut cur);
        std::mem::swap(&mut prev_start, &mut cur_start);
    }
    // best end column = argmin over the last row (now in prev)
    let mut best = prev[0];
    let mut end = 0usize;
    for j in 0..=n {
        if prev[j] < best {
            best = prev[j];
            end = j;
        }
    }
    (best, prev_start[end], end)
}

struct EndMatch {
    edits: usize,
    identity: f64,
    start: usize,
    end: usize,
}

fn best_in_window(window: &[u8], pat: &[u8]) -> EndMatch {
    let (edits, start, end) = fit_align(window, pat);
    let identity = if pat.is_empty() { 0.0 } else { 1.0 - edits as f64 / pat.len() as f64 };
    EndMatch { edits, identity, start, end }
}

// ---- FASTQ I/O (shared shape with the fastp module) ------------------------

struct Fq {
    id: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    qual: Vec<u8>,
}

struct FastqReader<R: BufRead> {
    inner: R,
    buf: Vec<u8>,
}
impl<R: BufRead> FastqReader<R> {
    fn new(inner: R) -> Self {
        FastqReader { inner, buf: Vec::with_capacity(256) }
    }
    fn line(&mut self) -> Result<Option<Vec<u8>>> {
        self.buf.clear();
        if self.inner.read_until(b'\n', &mut self.buf)? == 0 {
            return Ok(None);
        }
        let mut v = std::mem::take(&mut self.buf);
        while matches!(v.last(), Some(b'\n') | Some(b'\r')) {
            v.pop();
        }
        Ok(Some(v))
    }
    fn next_record(&mut self) -> Result<Option<Fq>> {
        let id = match self.line()? {
            Some(l) if !l.is_empty() => l,
            _ => return Ok(None),
        };
        let seq = self.line()?.context("truncated FASTQ: sequence")?;
        let plus = self.line()?.context("truncated FASTQ: '+'")?;
        let qual = self.line()?.context("truncated FASTQ: quality")?;
        Ok(Some(Fq { id, seq, plus, qual }))
    }
}

fn write_record<W: Write>(w: &mut W, r: &Fq) -> Result<()> {
    w.write_all(&r.id)?;
    w.write_all(b"\n")?;
    w.write_all(&r.seq)?;
    w.write_all(b"\n")?;
    w.write_all(if r.plus.is_empty() { b"+" } else { &r.plus })?;
    w.write_all(b"\n")?;
    w.write_all(&r.qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

/// Assign a read to a barcode (or `None`) and compute the trim coordinates.
/// Returns `(Option<barcode_name>, left_trim, right_trim)`.
fn classify(seq: &[u8], barcodes: &[Barcode], rc: &[Vec<u8>], p: &DemuxParams) -> (Option<String>, usize, usize) {
    let len = seq.len();
    let w = p.window.min(len);
    let left = &seq[..w];
    let right = &seq[len - w..];

    // For each barcode, the best of (forward @ left end, revcomp @ right end).
    // Track index + both end matches together so extraction can't panic.
    let mut best_id = 0.0f64;
    let mut second_id = 0.0f64;
    let mut best: Option<(usize, EndMatch, Option<EndMatch>)> = None;

    for (bi, bc) in barcodes.iter().enumerate() {
        let lm = best_in_window(left, &bc.seq);
        let rm = if p.both_ends { Some(best_in_window(right, &rc[bi])) } else { None };
        let id = lm.identity.max(rm.as_ref().map(|m| m.identity).unwrap_or(0.0));
        if id > best_id {
            second_id = best_id;
            best_id = id;
            best = Some((bi, lm, rm));
        } else if id > second_id {
            second_id = id;
        }
    }

    let (idx, lm, rm) = match best {
        Some(t) => t,
        None => return (None, 0, 0),
    };

    // acceptance: identity threshold (or max_errors) + margin over runner-up
    let pass_thresh = match p.max_errors {
        Some(me) => {
            let le = lm.edits;
            let re = rm.as_ref().map(|m| m.edits).unwrap_or(usize::MAX);
            le <= me || re <= me
        }
        None => best_id >= p.min_score,
    };
    if !pass_thresh || (best_id - second_id) < p.min_score_diff {
        return (None, 0, 0);
    }

    // trim coordinates: barcode at the 5' end trims [0..end]; barcode at the 3'
    // end trims [right_start..] (mapped back to read coordinates).
    let mut left_trim = 0;
    let mut right_trim = 0;
    if !p.trim {
        return (Some(barcodes[idx].name.clone()), 0, 0);
    }
    let left_ok = match p.max_errors {
        Some(me) => lm.edits <= me,
        None => lm.identity >= p.min_score,
    };
    if left_ok {
        left_trim = lm.end; // window starts at read 0
    }
    if let Some(rm) = rm {
        let right_ok = match p.max_errors {
            Some(me) => rm.edits <= me,
            None => rm.identity >= p.min_score,
        };
        if right_ok {
            // window covers read[len-w..]; barcode starts at window pos rm.start
            right_trim = w - rm.start;
        }
    }
    (Some(barcodes[idx].name.clone()), left_trim, right_trim)
}

/// Demultiplex a FASTQ file into `outdir`, one `<barcode>.fastq` per barcode plus
/// `unclassified.fastq`.
pub fn demux_file(input: &Path, outdir: &Path, barcodes: &[Barcode], p: &DemuxParams) -> Result<DemuxStats> {
    std::fs::create_dir_all(outdir)
        .with_context(|| format!("creating output dir '{}'", outdir.display()))?;
    let rc: Vec<Vec<u8>> = barcodes.iter().map(|b| revcomp(&b.seq)).collect();

    let (reader, _gz) = crate::open_reader(input)?;
    let mut fqr = FastqReader::new(reader);

    // lazily-opened per-label writers
    let mut writers: BTreeMap<String, std::io::BufWriter<std::fs::File>> = BTreeMap::new();
    let mut stats = DemuxStats::default();

    while let Some(mut r) = fqr.next_record()? {
        stats.total += 1;
        let (label, lt, rt) = classify(&r.seq, barcodes, &rc, p);
        let out_label = match &label {
            Some(name) => {
                stats.classified += 1;
                *stats.per_barcode.entry(name.clone()).or_insert(0) += 1;
                name.clone()
            }
            None => {
                stats.unclassified += 1;
                "unclassified".to_string()
            }
        };
        if p.trim && (lt > 0 || rt > 0) {
            let len = r.seq.len();
            let end = len.saturating_sub(rt);
            let start = lt.min(end);
            r.seq = r.seq[start..end].to_vec();
            r.qual = r.qual[start..end].to_vec();
        }
        // Look up (or create) the writer for this label in one pass — the Entry
        // API returns a `&mut` with no possibility of a panicking re-lookup.
        let w = match writers.entry(out_label.clone()) {
            std::collections::btree_map::Entry::Occupied(e) => e.into_mut(),
            std::collections::btree_map::Entry::Vacant(e) => {
                let path = outdir.join(format!("{out_label}.fastq"));
                let f = std::fs::File::create(&path)
                    .with_context(|| format!("creating '{}'", path.display()))?;
                e.insert(std::io::BufWriter::new(f))
            }
        };
        write_record(w, &r)?;
    }
    for (_, mut w) in writers {
        w.flush()?;
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bc(name: &str, seq: &str) -> Barcode {
        Barcode { name: name.into(), seq: seq.as_bytes().to_vec() }
    }

    #[test]
    fn fit_align_exact_and_mismatch() {
        let (e, s, en) = fit_align(b"GGGGAAAACCCCTTTTGGGG", b"AAAACCCCTTTT");
        assert_eq!(e, 0);
        assert_eq!((s, en), (4, 16));
        // one mismatch
        let (e2, _, _) = fit_align(b"GGGGAAAACGCCTTTTGGGG", b"AAAACCCCTTTT");
        assert_eq!(e2, 1);
    }

    #[test]
    fn classify_and_trim_left_barcode() {
        let bcs = vec![bc("bc01", "AAAACCCCGGGGTTTT"), bc("bc02", "TTTTGGGGCCCCAAAA")];
        let rc: Vec<Vec<u8>> = bcs.iter().map(|b| revcomp(&b.seq)).collect();
        // read = bc01 + insert
        let read = b"AAAACCCCGGGGTTTTACGTACGTACGTACGTACGT".to_vec();
        let p = DemuxParams { both_ends: false, min_score_diff: 0.05, ..Default::default() };
        let (label, lt, rt) = classify(&read, &bcs, &rc, &p);
        assert_eq!(label.as_deref(), Some("bc01"));
        assert_eq!(lt, 16);
        assert_eq!(rt, 0);
    }

    #[test]
    fn unclassified_when_no_match() {
        let bcs = vec![bc("bc01", "AAAACCCCGGGGTTTT")];
        let rc: Vec<Vec<u8>> = bcs.iter().map(|b| revcomp(&b.seq)).collect();
        let read = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let p = DemuxParams { both_ends: false, min_score: 0.85, ..Default::default() };
        let (label, _, _) = classify(&read, &bcs, &rc, &p);
        assert!(label.is_none());
    }

    #[test]
    fn end_to_end_demux() {
        let dir = std::env::temp_dir().join(format!("cleaver-demux-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let bcf = dir.join("bc.fasta");
        std::fs::write(&bcf, b">bc01\nAAAACCCCGGGGTTTT\n>bc02\nTTTTGGGGCCCCAAAA\n").unwrap();
        let reads = dir.join("reads.fastq");
        // r1 -> bc01, r2 -> bc02 (its barcode at 5'), r3 -> unclassified
        std::fs::write(
            &reads,
            b"@r1\nAAAACCCCGGGGTTTTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
              @r2\nTTTTGGGGCCCCAAAAGTGTGTGTGTGTGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n\
              @r3\nACGTACGTACGTACGTACGTACGTACGTAC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        )
        .unwrap();
        let bcs = read_barcodes(&bcf).unwrap();
        let p = DemuxParams { both_ends: false, min_score: 0.85, ..Default::default() };
        let stats = demux_file(&reads, &dir, &bcs, &p).unwrap();
        assert_eq!(stats.total, 3);
        assert_eq!(stats.classified, 2);
        assert_eq!(stats.unclassified, 1);
        assert_eq!(stats.per_barcode.get("bc01"), Some(&1));
        // bc01 output exists and the barcode was trimmed off
        let out = std::fs::read_to_string(dir.join("bc01.fastq")).unwrap();
        assert!(out.contains("ACGTACGTACGTAC"));
        assert!(!out.contains("AAAACCCCGGGGTTTT"));
    }

    #[test]
    fn short_and_empty_reads_route_to_unclassified() {
        // A read shorter than the 16bp barcode and an empty read must not panic
        // (window = min(w, len); empty slices) and must land in `unclassified`.
        let dir = std::env::temp_dir().join(format!("cleaver-demux-edge-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let bcf = dir.join("bc.fasta");
        std::fs::write(&bcf, b">bc01\nAAAACCCCGGGGTTTT\n").unwrap();
        let reads = dir.join("reads.fastq");
        std::fs::write(
            &reads,
            b"@short\nACGT\n+\nIIII\n\
              @empty\n\n+\n\n\
              @ok\nAAAACCCCGGGGTTTTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
        )
        .unwrap();
        let bcs = read_barcodes(&bcf).unwrap();
        let p = DemuxParams { both_ends: false, min_score: 0.85, ..Default::default() };
        let stats = demux_file(&reads, &dir, &bcs, &p).unwrap();
        assert_eq!(stats.total, 3);
        assert_eq!(stats.classified, 1, "only the full bc01 read matches");
        assert_eq!(stats.unclassified, 2, "short + empty reads go to unclassified");
        assert_eq!(stats.per_barcode.get("bc01"), Some(&1));
    }
}
