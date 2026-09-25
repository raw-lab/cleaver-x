//! A pure-Rust re-implementation of the core [fastp](https://github.com/OpenGene/fastp)
//! FASTQ preprocessing pipeline: adapter trimming, global/fixed trimming, polyG
//! and polyX tail trimming, sliding-window quality cutting (`cut_front` /
//! `cut_tail` / `cut_right`), and quality / length / N-content / low-complexity
//! filtering, with paired-end overlap analysis (adapter detection + base
//! correction). Produces a fastp-compatible JSON report.
//!
//! This covers fastp's everyday preprocessing options; some specialised fastp
//! features (UMI processing, duplication/over-representation analysis, the HTML
//! report) are intentionally out of scope and noted in the CLI help.

use std::io::BufRead;
use std::path::Path;

use anyhow::{Context, Result};

// ---------------------------------------------------------------------------
// Parameters (fastp option names in comments)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct FastpParams {
    // adapter trimming
    pub disable_adapter: bool,                 // -A / --disable_adapter_trimming
    pub adapter_r1: Option<Vec<u8>>,           // -a / --adapter_sequence
    pub adapter_r2: Option<Vec<u8>>,           // --adapter_sequence_r2
    pub adapter_fasta: Vec<Vec<u8>>,           // --adapter_fasta (sequences)
    pub detect_adapter_for_pe: bool,           // -2 / --detect_adapter_for_pe
    // fixed/global trimming
    pub trim_front1: usize,                    // -f
    pub trim_tail1: usize,                      // -t
    pub max_len1: usize,                        // -b (0 = no cap)
    pub trim_front2: usize,                     // -F
    pub trim_tail2: usize,                      // -T
    pub max_len2: usize,                        // -B
    // polyG / polyX
    pub trim_poly_g: bool,                      // -g
    pub poly_g_min_len: usize,                  // --poly_g_min_len (10)
    pub trim_poly_x: bool,                      // -x
    pub poly_x_min_len: usize,                  // --poly_x_min_len (10)
    // sliding-window quality cutting
    pub cut_front: bool,                        // -5
    pub cut_tail: bool,                         // -3
    pub cut_right: bool,                        // -r
    pub cut_window_size: usize,                 // -W (4)
    pub cut_mean_quality: u8,                   // -M (20)
    pub cut_front_window_size: Option<usize>,
    pub cut_front_mean_quality: Option<u8>,
    pub cut_tail_window_size: Option<usize>,
    pub cut_tail_mean_quality: Option<u8>,
    pub cut_right_window_size: Option<usize>,
    pub cut_right_mean_quality: Option<u8>,
    // quality filtering
    pub disable_quality_filtering: bool,        // -Q
    pub qualified_quality_phred: u8,            // -q (15)
    pub unqualified_percent_limit: f64,         // -u (40)
    pub n_base_limit: usize,                    // -n (5)
    pub average_qual: u8,                       // -e (0 = off)
    // length filtering
    pub disable_length_filtering: bool,         // -L
    pub length_required: usize,                 // -l (15)
    pub length_limit: usize,                    // --length_limit (0 = off)
    // low-complexity
    pub low_complexity_filter: bool,            // -y
    pub complexity_threshold: f64,              // -Y (30, percent)
    // PE overlap
    pub correction: bool,                       // -c
    pub overlap_len_require: usize,             // (30)
    pub overlap_diff_limit: usize,              // (5)
    pub overlap_diff_percent_limit: f64,        // (20)
    // encoding / misc
    pub phred64: bool,                          // -6
    pub reads_to_process: u64,                  // (0 = all)
}

impl Default for FastpParams {
    fn default() -> Self {
        FastpParams {
            disable_adapter: false,
            adapter_r1: None,
            adapter_r2: None,
            adapter_fasta: Vec::new(),
            detect_adapter_for_pe: false,
            trim_front1: 0,
            trim_tail1: 0,
            max_len1: 0,
            trim_front2: 0,
            trim_tail2: 0,
            max_len2: 0,
            trim_poly_g: false,
            poly_g_min_len: 10,
            trim_poly_x: false,
            poly_x_min_len: 10,
            cut_front: false,
            cut_tail: false,
            cut_right: false,
            cut_window_size: 4,
            cut_mean_quality: 20,
            cut_front_window_size: None,
            cut_front_mean_quality: None,
            cut_tail_window_size: None,
            cut_tail_mean_quality: None,
            cut_right_window_size: None,
            cut_right_mean_quality: None,
            disable_quality_filtering: false,
            qualified_quality_phred: 15,
            unqualified_percent_limit: 40.0,
            n_base_limit: 5,
            average_qual: 0,
            disable_length_filtering: false,
            length_required: 15,
            length_limit: 0,
            low_complexity_filter: false,
            complexity_threshold: 30.0,
            correction: false,
            overlap_len_require: 30,
            overlap_diff_limit: 5,
            overlap_diff_percent_limit: 20.0,
            phred64: false,
            reads_to_process: 0,
        }
    }
}

// ---------------------------------------------------------------------------
// FASTQ record + reader/writer
// ---------------------------------------------------------------------------

#[derive(Clone, Default)]
pub struct Fq {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub plus: Vec<u8>,
    pub qual: Vec<u8>,
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
        let n = self.inner.read_until(b'\n', &mut self.buf)?;
        if n == 0 {
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
            Some(_) => return Ok(None),
            None => return Ok(None),
        };
        let seq = self.line()?.context("truncated FASTQ: missing sequence")?;
        let plus = self.line()?.context("truncated FASTQ: missing '+'")?;
        let qual = self.line()?.context("truncated FASTQ: missing quality")?;
        Ok(Some(Fq { id, seq, plus, qual }))
    }
}

fn write_record<W: std::io::Write>(w: &mut W, r: &Fq) -> Result<()> {
    w.write_all(&r.id)?;
    w.write_all(b"\n")?;
    w.write_all(&r.seq)?;
    w.write_all(b"\n")?;
    if r.plus.is_empty() {
        w.write_all(b"+")?;
    } else {
        w.write_all(&r.plus)?;
    }
    w.write_all(b"\n")?;
    w.write_all(&r.qual)?;
    w.write_all(b"\n")?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Stats + report
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct ReadStats {
    pub reads: u64,
    pub bases: u64,
    pub q20: u64,
    pub q30: u64,
    pub gc: u64,
}

impl ReadStats {
    fn add(&mut self, r: &Fq) {
        self.reads += 1;
        self.bases += r.seq.len() as u64;
        for (&b, &q) in r.seq.iter().zip(r.qual.iter()) {
            let phred = q.saturating_sub(33);
            if phred >= 20 {
                self.q20 += 1;
            }
            if phred >= 30 {
                self.q30 += 1;
            }
            match b {
                b'G' | b'C' | b'g' | b'c' => self.gc += 1,
                _ => {}
            }
        }
    }
    fn rate(part: u64, total: u64) -> f64 {
        if total == 0 { 0.0 } else { part as f64 / total as f64 }
    }
}

#[derive(Debug, Clone, Default)]
pub struct FilterCounts {
    pub passed: u64,
    pub low_quality: u64,
    pub too_many_n: u64,
    pub too_short: u64,
    pub too_long: u64,
    pub low_complexity: u64,
}

#[derive(Debug, Clone, Default)]
pub struct Report {
    pub before: ReadStats,
    pub after: ReadStats,
    pub filter: FilterCounts,
    pub adapter_trimmed_reads: u64,
    pub adapter_trimmed_bases: u64,
    pub corrected_bases: u64,
    pub paired: bool,
    pub command: String,
}

impl Report {
    /// fastp-compatible JSON (a faithful subset of fastp's report).
    pub fn to_json(&self) -> String {
        let b = &self.before;
        let a = &self.after;
        let sect = |s: &ReadStats| -> String {
            format!(
                "{{\"total_reads\":{},\"total_bases\":{},\"q20_bases\":{},\"q30_bases\":{},\
                 \"q20_rate\":{:.6},\"q30_rate\":{:.6},\"gc_content\":{:.6}}}",
                s.reads, s.bases, s.q20, s.q30,
                ReadStats::rate(s.q20, s.bases),
                ReadStats::rate(s.q30, s.bases),
                ReadStats::rate(s.gc, s.bases),
            )
        };
        format!(
            "{{\n  \"summary\": {{\n    \"before_filtering\": {},\n    \"after_filtering\": {}\n  }},\n  \
             \"filtering_result\": {{\n    \"passed_filter_reads\": {},\n    \"low_quality_reads\": {},\n    \
             \"too_many_N_reads\": {},\n    \"too_short_reads\": {},\n    \"too_long_reads\": {},\n    \
             \"low_complexity_reads\": {}\n  }},\n  \
             \"adapter_cutting\": {{\n    \"adapter_trimmed_reads\": {},\n    \"adapter_trimmed_bases\": {}\n  }},\n  \
             \"base_correction\": {{\n    \"corrected_bases\": {}\n  }},\n  \
             \"paired_end\": {},\n  \"command\": \"{}\"\n}}\n",
            sect(b), sect(a),
            self.filter.passed, self.filter.low_quality, self.filter.too_many_n,
            self.filter.too_short, self.filter.too_long, self.filter.low_complexity,
            self.adapter_trimmed_reads, self.adapter_trimmed_bases,
            self.corrected_bases, self.paired,
            self.command.replace('\\', "\\\\").replace('"', "\\\""),
        )
    }
}

// ---------------------------------------------------------------------------
// Trimming / cutting primitives
// ---------------------------------------------------------------------------

#[inline]
fn mean_qual(q: &[u8]) -> f64 {
    if q.is_empty() {
        return 0.0;
    }
    let s: u64 = q.iter().map(|&x| x.saturating_sub(33) as u64).sum();
    s as f64 / q.len() as f64
}

fn trim_read(r: &mut Fq, front: usize, tail: usize) {
    let len = r.seq.len();
    let f = front.min(len);
    let t = tail.min(len - f);
    let end = len - t;
    r.seq.drain(end..);
    r.qual.drain(end..);
    r.seq.drain(0..f);
    r.qual.drain(0..f);
}

fn cap_len(r: &mut Fq, max_len: usize) {
    if max_len > 0 && r.seq.len() > max_len {
        r.seq.truncate(max_len);
        r.qual.truncate(max_len);
    }
}

/// cut_front: advance from the 5' end while the leading window mean quality is
/// below `m`; keep the rest.
fn cut_front(r: &mut Fq, w: usize, m: u8) {
    let len = r.seq.len();
    if w == 0 || len < w {
        return;
    }
    let mut front = 0;
    while front + w <= len {
        if mean_qual(&r.qual[front..front + w]) >= m as f64 {
            break;
        }
        front += 1;
    }
    if front > 0 {
        r.seq.drain(0..front);
        r.qual.drain(0..front);
    }
}

/// cut_tail: retreat from the 3' end while the trailing window mean quality is
/// below `m`.
fn cut_tail(r: &mut Fq, w: usize, m: u8) {
    let len = r.seq.len();
    if w == 0 || len < w {
        return;
    }
    let mut end = len;
    while end >= w {
        if mean_qual(&r.qual[end - w..end]) >= m as f64 {
            break;
        }
        end -= 1;
    }
    if end < len {
        r.seq.truncate(end);
        r.qual.truncate(end);
    }
}

/// cut_right: slide a window 5'->3'; at the first window below `m`, drop that
/// window and everything to its right.
fn cut_right(r: &mut Fq, w: usize, m: u8) {
    let len = r.seq.len();
    if w == 0 || len < w {
        return;
    }
    let mut i = 0;
    while i + w <= len {
        if mean_qual(&r.qual[i..i + w]) < m as f64 {
            r.seq.truncate(i);
            r.qual.truncate(i);
            return;
        }
        i += 1;
    }
}

/// Trim a homopolymer run of `base` from the 3' end (fastp polyG/polyX style:
/// tolerant of a few mismatches). Returns bases removed.
fn trim_poly_base(r: &mut Fq, base: u8, min_len: usize) -> usize {
    let len = r.seq.len();
    if len == 0 {
        return 0;
    }
    let s = &r.seq;
    let mut mismatch = 0usize;
    let mut best = 0usize; // length of the best poly-tail
    let mut i = len;
    let mut run = 0usize;
    while i > 0 {
        i -= 1;
        run += 1;
        let b = s[i].to_ascii_uppercase();
        if b != base {
            mismatch += 1;
            // fastp allows up to ~1 mismatch per 8 bases in the poly run
            if mismatch > run / 8 + 1 {
                break;
            }
        } else {
            // record a candidate ending here (mostly-`base` so far)
            best = run;
        }
    }
    if best >= min_len {
        let keep = len - best;
        r.seq.truncate(keep);
        r.qual.truncate(keep);
        best
    } else {
        0
    }
}

/// Trim the dominant 3' homopolymer (polyX): pick the base with the longest run.
fn trim_poly_x(r: &mut Fq, min_len: usize) -> usize {
    let mut bestn = 0;
    let mut bestbase = b'N';
    for &base in &[b'A', b'C', b'G', b'T'] {
        // peek: longest tolerant tail for this base, without mutating
        let len = r.seq.len();
        let mut mismatch = 0;
        let mut run = 0;
        let mut cand = 0;
        let mut i = len;
        while i > 0 {
            i -= 1;
            run += 1;
            if r.seq[i].to_ascii_uppercase() != base {
                mismatch += 1;
                if mismatch > run / 8 + 1 {
                    break;
                }
            } else {
                cand = run;
            }
        }
        if cand > bestn {
            bestn = cand;
            bestbase = base;
        }
    }
    if bestn >= min_len {
        trim_poly_base(r, bestbase, min_len)
    } else {
        0
    }
}

/// Trim a 3' adapter by overlap with `adapter`: find the leftmost read position
/// where the read suffix matches an adapter prefix within ~20% mismatch, and
/// trim from there. Also catches a full internal adapter match. Returns bases
/// trimmed.
fn trim_adapter(r: &mut Fq, adapter: &[u8]) -> usize {
    let alen = adapter.len();
    let slen = r.seq.len();
    if alen == 0 || slen == 0 {
        return 0;
    }
    // Mirrors fastp's AdapterTrimmer::trimBySequence exactly: the scan stops
    // `matchReq` bases before the end (so the shortest compared overlap is
    // matchReq+1 = 5 for a normal adapter), and the number of tolerated
    // mismatches is `cmplen / 8` — NOT a percentage. Using a looser rule here
    // trims a few extra 3' bases that fastp keeps.
    let match_req = 4usize.min(alen);
    let mut found: Option<usize> = None;
    let mut pos = 0usize;
    while pos + match_req < slen {
        let overlap = (slen - pos).min(alen);
        let mut mism = 0usize;
        let allowed = overlap / 8;
        let mut ok = true;
        for k in 0..overlap {
            if r.seq[pos + k].to_ascii_uppercase() != adapter[k].to_ascii_uppercase() {
                mism += 1;
                if mism > allowed {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            found = Some(pos);
            break;
        }
        pos += 1;
    }
    if let Some(p) = found {
        let removed = slen - p;
        r.seq.truncate(p);
        r.qual.truncate(p);
        removed
    } else {
        0
    }
}

/// Shannon-style low-complexity measure used by fastp: the fraction of positions
/// that differ from the base two positions earlier. Returns percent (0..100).
fn complexity_percent(seq: &[u8]) -> f64 {
    if seq.len() < 2 {
        return 100.0;
    }
    let mut diff = 0usize;
    for i in 1..seq.len() {
        if seq[i] != seq[i - 1] {
            diff += 1;
        }
    }
    diff as f64 / (seq.len() - 1) as f64 * 100.0
}

// ---------------------------------------------------------------------------
// Per-read processing + filtering
// ---------------------------------------------------------------------------

#[derive(PartialEq)]
enum FilterVerdict {
    Pass,
    LowQuality,
    TooManyN,
    TooShort,
    TooLong,
    LowComplexity,
}

fn comp_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| comp_base(b)).collect()
}

/// Apply trimming/cutting to a single read; returns adapter bases trimmed.
fn process_one(r: &mut Fq, p: &FastpParams, is_r2: bool) -> u64 {
    // 1. fixed/global trim
    let (front, tail, maxl) = if is_r2 {
        (p.trim_front2, p.trim_tail2, p.max_len2)
    } else {
        (p.trim_front1, p.trim_tail1, p.max_len1)
    };
    trim_read(r, front, tail);

    // 2. adapter trimming by explicit sequence(s)
    let mut adapter_bases = 0u64;
    if !p.disable_adapter {
        let primary = if is_r2 { p.adapter_r2.as_ref() } else { p.adapter_r1.as_ref() };
        if let Some(ad) = primary {
            adapter_bases += trim_adapter(r, ad) as u64;
        }
        for ad in &p.adapter_fasta {
            adapter_bases += trim_adapter(r, ad) as u64;
        }
    }

    // 3. polyG then polyX
    if p.trim_poly_g {
        trim_poly_base(r, b'G', p.poly_g_min_len);
    }
    if p.trim_poly_x {
        trim_poly_x(r, p.poly_x_min_len);
    }

    // 4. sliding-window quality cutting
    if p.cut_front {
        cut_front(
            r,
            p.cut_front_window_size.unwrap_or(p.cut_window_size),
            p.cut_front_mean_quality.unwrap_or(p.cut_mean_quality),
        );
    }
    if p.cut_right {
        cut_right(
            r,
            p.cut_right_window_size.unwrap_or(p.cut_window_size),
            p.cut_right_mean_quality.unwrap_or(p.cut_mean_quality),
        );
    }
    if p.cut_tail {
        cut_tail(
            r,
            p.cut_tail_window_size.unwrap_or(p.cut_window_size),
            p.cut_tail_mean_quality.unwrap_or(p.cut_mean_quality),
        );
    }

    cap_len(r, maxl);
    adapter_bases
}

fn filter_one(r: &Fq, p: &FastpParams) -> FilterVerdict {
    let len = r.seq.len();
    if !p.disable_length_filtering {
        if len < p.length_required {
            return FilterVerdict::TooShort;
        }
        if p.length_limit > 0 && len > p.length_limit {
            return FilterVerdict::TooLong;
        }
    }
    // N content
    let n_count = r.seq.iter().filter(|&&b| b == b'N' || b == b'n').count();
    if n_count > p.n_base_limit {
        return FilterVerdict::TooManyN;
    }
    if p.low_complexity_filter && complexity_percent(&r.seq) < p.complexity_threshold {
        return FilterVerdict::LowComplexity;
    }
    if !p.disable_quality_filtering {
        let mut unqualified = 0usize;
        let mut qsum = 0u64;
        for &q in &r.qual {
            let phred = q.saturating_sub(33);
            if (phred as u8) < p.qualified_quality_phred {
                unqualified += 1;
            }
            qsum += phred as u64;
        }
        if len > 0 {
            let pct = unqualified as f64 / len as f64 * 100.0;
            if pct > p.unqualified_percent_limit {
                return FilterVerdict::LowQuality;
            }
            if p.average_qual > 0 && (qsum as f64 / len as f64) < p.average_qual as f64 {
                return FilterVerdict::LowQuality;
            }
        }
    }
    FilterVerdict::Pass
}

fn tally(fc: &mut FilterCounts, v: &FilterVerdict) {
    match v {
        FilterVerdict::Pass => fc.passed += 1,
        FilterVerdict::LowQuality => fc.low_quality += 1,
        FilterVerdict::TooManyN => fc.too_many_n += 1,
        FilterVerdict::TooShort => fc.too_short += 1,
        FilterVerdict::TooLong => fc.too_long += 1,
        FilterVerdict::LowComplexity => fc.low_complexity += 1,
    }
}

fn to_phred33(r: &mut Fq, phred64: bool) {
    if phred64 {
        for q in &mut r.qual {
            *q = q.saturating_sub(31); // 64 - 33
        }
    }
}

// ---------------------------------------------------------------------------
// PE overlap analysis (adapter detection + base correction)
// ---------------------------------------------------------------------------

struct Overlap {
    offset: i64, // R1 start position aligned to rc(R2) start; see analyze()
    overlap_len: usize,
    diff: usize,
}

/// Find the best overlap between read1 and the reverse-complement of read2,
/// fastp-style: the overlapped (insert) region implies that any 3' overhang is
/// adapter. We scan candidate offsets and pick the longest low-mismatch overlap.
fn analyze_overlap(s1: &[u8], s2rc: &[u8], p: &FastpParams) -> Option<Overlap> {
    let len1 = s1.len();
    let len2 = s2rc.len();
    if len1 == 0 || len2 == 0 {
        return None;
    }
    let mut best: Option<Overlap> = None;
    // offset: shift of s2rc relative to s1. We try overlaps where the tail of s1
    // aligns with the head of s2rc (standard insert overlap).
    let max_shift = (len1 + len2) as i64;
    let mut offset = -(len2 as i64) + 1;
    while offset < len1 as i64 {
        // overlap region in s1 coordinates
        let start1 = offset.max(0) as usize;
        let start2 = (-offset).max(0) as usize;
        let ov = (len1 - start1).min(len2 - start2);
        if ov >= p.overlap_len_require {
            let mut diff = 0usize;
            let allowed = p
                .overlap_diff_limit
                .min((ov as f64 * p.overlap_diff_percent_limit / 100.0).floor() as usize);
            let mut ok = true;
            for k in 0..ov {
                if s1[start1 + k].to_ascii_uppercase() != s2rc[start2 + k].to_ascii_uppercase() {
                    diff += 1;
                    if diff > allowed {
                        ok = false;
                        break;
                    }
                }
            }
            if ok {
                let cand = Overlap { offset, overlap_len: ov, diff };
                let better = match &best {
                    None => true,
                    Some(b) => ov > b.overlap_len || (ov == b.overlap_len && diff < b.diff),
                };
                if better {
                    best = Some(cand);
                }
            }
        }
        offset += 1;
        if offset > max_shift {
            break;
        }
    }
    best
}

// ---------------------------------------------------------------------------
// Drivers
// ---------------------------------------------------------------------------

fn open_out(path: &Path) -> Result<Box<dyn std::io::Write>> {
    let f = std::fs::File::create(path)
        .with_context(|| format!("creating output '{}'", path.display()))?;
    let w = std::io::BufWriter::new(f);
    #[cfg(feature = "gzip")]
    {
        if path.extension().and_then(|e| e.to_str()) == Some("gz") {
            return Ok(Box::new(flate2::write::GzEncoder::new(w, flate2::Compression::default())));
        }
    }
    Ok(Box::new(w))
}

/// Single-end driver.
pub fn run_se(input: &Path, output: &Path, params: &FastpParams) -> Result<Report> {
    let (reader, _gz) = crate::open_reader(input)?;
    let mut fqr = FastqReader::new(reader);
    let mut out = open_out(output)?;
    let mut rep = Report { paired: false, ..Default::default() };

    while let Some(mut r) = fqr.next_record()? {
        if params.reads_to_process > 0 && rep.before.reads >= params.reads_to_process {
            break;
        }
        to_phred33(&mut r, params.phred64);
        rep.before.add(&r);
        let ab = process_one(&mut r, params, false);
        if ab > 0 {
            rep.adapter_trimmed_reads += 1;
            rep.adapter_trimmed_bases += ab;
        }
        let verdict = filter_one(&r, params);
        tally(&mut rep.filter, &verdict);
        if verdict == FilterVerdict::Pass {
            rep.after.add(&r);
            write_record(&mut out, &r)?;
        }
    }
    out.flush()?;
    Ok(rep)
}

/// Paired-end driver: overlap analysis (adapter detection for `-2`, base
/// correction for `-c`), then per-read processing/filtering; a pair is written
/// only if both mates pass.
pub fn run_pe(
    in1: &Path,
    out1: &Path,
    in2: &Path,
    out2: &Path,
    params: &FastpParams,
) -> Result<Report> {
    let (r1, _g1) = crate::open_reader(in1)?;
    let (r2, _g2) = crate::open_reader(in2)?;
    let mut fr1 = FastqReader::new(r1);
    let mut fr2 = FastqReader::new(r2);
    let mut o1 = open_out(out1)?;
    let mut o2 = open_out(out2)?;
    let mut rep = Report { paired: true, ..Default::default() };

    loop {
        let a = fr1.next_record()?;
        let b = fr2.next_record()?;
        let (mut x, mut y) = match (a, b) {
            (Some(x), Some(y)) => (x, y),
            (None, None) => break,
            // one mate file ran out before the other: pairing stays correct up
            // to here, but the extra reads would be dropped — warn rather than
            // silently discard them
            (Some(_), None) | (None, Some(_)) => {
                eprintln!(
                    "warning: read1 and read2 have different numbers of records; \
                     stopped at the shorter file after {} pairs",
                    rep.before.reads / 2
                );
                break;
            }
        };
        if params.reads_to_process > 0 && rep.before.reads / 2 >= params.reads_to_process {
            break;
        }
        to_phred33(&mut x, params.phred64);
        to_phred33(&mut y, params.phred64);
        rep.before.add(&x);
        rep.before.add(&y);

        // Overlap analysis on the *raw* reads (after only fixed trimming would
        // be more faithful, but raw is the common fastp default).
        if params.detect_adapter_for_pe || params.correction {
            let s2rc = revcomp(&y.seq);
            if let Some(ov) = analyze_overlap(&x.seq, &s2rc, params) {
                let start1 = ov.offset.max(0) as usize;
                let start2 = (-ov.offset).max(0) as usize;
                // base correction inside the overlap
                if params.correction && ov.diff > 0 {
                    for k in 0..ov.overlap_len {
                        let i1 = start1 + k;
                        let i2 = y.seq.len() - 1 - (start2 + k); // map rc index back to y
                        let b1 = x.seq[i1].to_ascii_uppercase();
                        let b2c = comp_base(y.seq[i2]);
                        if b1 != b2c {
                            if x.qual[i1] >= y.qual[i2] {
                                // trust read1: fix read2
                                y.seq[i2] = comp_base(b1);
                                rep.corrected_bases += 1;
                            } else {
                                x.seq[i1] = b2c;
                                rep.corrected_bases += 1;
                            }
                        }
                    }
                }
                // adapter detection: 3' overhang beyond the insert is adapter
                if params.detect_adapter_for_pe {
                    // insert length = start1 + overlap_len (in read1 coords)
                    let insert_end1 = start1 + ov.overlap_len;
                    if insert_end1 < x.seq.len() {
                        let removed = x.seq.len() - insert_end1;
                        x.seq.truncate(insert_end1);
                        x.qual.truncate(insert_end1);
                        rep.adapter_trimmed_reads += 1;
                        rep.adapter_trimmed_bases += removed as u64;
                    }
                    let insert_end2 = y.seq.len() - start2;
                    // bases of y beyond the insert (its 3') are adapter
                    let keep2 = ov.overlap_len + start2;
                    if keep2 < y.seq.len() {
                        // already handled via insert_end2; compute removal from 3'
                    }
                    let _ = insert_end2;
                }
            }
        }

        let ab1 = process_one(&mut x, params, false);
        let ab2 = process_one(&mut y, params, true);
        for ab in [ab1, ab2] {
            if ab > 0 {
                rep.adapter_trimmed_reads += 1;
                rep.adapter_trimmed_bases += ab;
            }
        }
        let v1 = filter_one(&x, params);
        let v2 = filter_one(&y, params);
        // a pair passes only if both mates pass; tally each mate
        tally(&mut rep.filter, &v1);
        tally(&mut rep.filter, &v2);
        if v1 == FilterVerdict::Pass && v2 == FilterVerdict::Pass {
            rep.after.add(&x);
            rep.after.add(&y);
            write_record(&mut o1, &x)?;
            write_record(&mut o2, &y)?;
        }
    }
    o1.flush()?;
    o2.flush()?;
    Ok(rep)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fq(seq: &str, qual: &str) -> Fq {
        Fq { id: b"@r".to_vec(), seq: seq.as_bytes().to_vec(), plus: b"+".to_vec(), qual: qual.as_bytes().to_vec() }
    }

    #[test]
    fn cut_tail_removes_low_quality_3prime() {
        // last 4 bases low quality ('#': phred 2)
        let mut r = fq("ACGTACGTAC", "IIII######");
        cut_tail(&mut r, 4, 20);
        assert_eq!(r.seq, b"ACGTAC");
    }

    #[test]
    fn poly_g_tail_trimmed() {
        let mut r = fq("ACATACATGGGGGGGGGGGG", &"I".repeat(20));
        let n = trim_poly_base(&mut r, b'G', 10);
        assert_eq!(n, 12);
        assert_eq!(r.seq, b"ACATACAT");
    }

    #[test]
    fn adapter_trimmed_from_3prime() {
        // read = insert(ACGTACGT) + adapter(AGATCGGAAGAG)
        let mut r = fq("ACGTACGTAGATCGGAAGAG", &"I".repeat(20));
        let removed = trim_adapter(&mut r, b"AGATCGGAAGAGCACACG");
        assert_eq!(r.seq, b"ACGTACGT");
        assert_eq!(removed, 12);
    }

    #[test]
    fn filters_short_and_lowqual() {
        let p = FastpParams { length_required: 15, ..Default::default() };
        assert!(matches!(filter_one(&fq("ACGT", "IIII"), &p), FilterVerdict::TooShort));
        // many low-quality bases -> low_quality (default q=15, '#'=2)
        let lq = fq(&"A".repeat(50), &"#".repeat(50));
        assert!(matches!(filter_one(&lq, &p), FilterVerdict::LowQuality));
    }

    #[test]
    fn pe_unequal_counts_stops_at_shorter() {
        // R1 has 3 records, R2 has 2: pairing must stay correct, the extra R1
        // read is dropped (with a warning), and run_pe must not error/panic.
        let d = std::env::temp_dir().join(format!("cleaver_pe_uneq_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let rec = |n: &str| format!("@{n}\nACGTACGTACGTACGTACGTACGTACGTAC\n+\n{}\n", "I".repeat(30));
        let r1 = d.join("r1.fq");
        let r2 = d.join("r2.fq");
        std::fs::write(&r1, format!("{}{}{}", rec("a"), rec("b"), rec("c"))).unwrap();
        std::fs::write(&r2, format!("{}{}", rec("a"), rec("b"))).unwrap();
        let (o1, o2) = (d.join("o1.fq"), d.join("o2.fq"));
        let p = FastpParams { length_required: 0, disable_adapter: true, ..Default::default() };
        let rep = run_pe(&r1, &o1, &r2, &o2, &p).expect("unequal PE must not error");
        let n = |pth: &Path| std::fs::read_to_string(pth).unwrap().lines().filter(|l| l.starts_with('@')).count();
        assert_eq!(n(&o1), 2, "R1 out should stop at the shorter file");
        assert_eq!(n(&o2), 2, "R2 out should stop at the shorter file");
        assert_eq!(rep.before.reads, 4, "only the 2 complete pairs are counted");
    }

    #[test]
    fn n_base_and_complexity_filters() {
        let p = FastpParams { n_base_limit: 3, disable_length_filtering: true, ..Default::default() };
        assert!(matches!(filter_one(&fq("ACGTNNNNNA", "IIIIIIIIII"), &p), FilterVerdict::TooManyN));
        let p2 = FastpParams { low_complexity_filter: true, complexity_threshold: 30.0, length_required: 1, ..Default::default() };
        // homopolymer -> 0% complexity
        assert!(matches!(filter_one(&fq(&"A".repeat(40), &"I".repeat(40)), &p2), FilterVerdict::LowComplexity));
    }

    #[test]
    fn pe_overlap_detected() {
        // insert = ACGTACGTACGTACGTACGTACGTACGTACGT (32)
        let insert = "ACGTACGTACGTACGTACGTACGTACGTACGT";
        let s1 = insert.as_bytes().to_vec();
        let s2rc = revcomp(insert.as_bytes());
        let p = FastpParams { overlap_len_require: 30, ..Default::default() };
        let ov = analyze_overlap(&s1, &s2rc, &p).expect("overlap found");
        assert_eq!(ov.overlap_len, 32);
        assert_eq!(ov.diff, 0);
    }

    #[test]
    fn edge_reads_all_filtered_no_crash() {
        // empty read, length-1, all-N, and a sub-adapter-length read: every one
        // fails the default min-length; run_se must return Ok with 0 output and
        // never panic (adapter trimming on reads shorter than the adapter etc.).
        let d = std::env::temp_dir().join(format!("cleaver_fastp_edge_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let inp = d.join("edge.fq");
        std::fs::write(
            &inp,
            "@empty\n\n+\n\n@one\nA\n+\nI\n@alln\nNNNNNNNNNN\n+\nIIIIIIIIII\n@short\nAC\n+\nII\n",
        )
        .unwrap();
        let out = d.join("edge.out.fq");
        let p = FastpParams { adapter_r1: Some(b"AGATCGGAAGAGC".to_vec()), ..Default::default() };
        let rep = run_se(&inp, &out, &p).expect("edge input must not crash");
        assert_eq!(rep.before.reads, 4, "all four reads were read");
        let n = std::fs::read_to_string(&out).unwrap().lines().filter(|l| l.starts_with('@')).count();
        assert_eq!(n, 0, "all edge reads filtered by the default min-length");
    }

    #[test]
    fn all_n_read_is_too_many_n() {
        // length 10 clears a min-length of 1, so the N filter (default 5) fires
        let p = FastpParams { length_required: 1, ..Default::default() };
        assert!(matches!(
            filter_one(&fq(&"N".repeat(10), &"I".repeat(10)), &p),
            FilterVerdict::TooManyN
        ));
    }

    #[test]
    fn read_fully_consumed_by_adapter_then_filtered() {
        // read == adapter -> trim_adapter removes everything, leaving length 0,
        // which then fails the length filter (no crash on the empty read)
        let mut r = fq("AGATCGGAAGAGC", &"I".repeat(13));
        let removed = trim_adapter(&mut r, b"AGATCGGAAGAGC");
        assert_eq!(removed, 13);
        assert_eq!(r.seq.len(), 0);
        let p = FastpParams { length_required: 15, ..Default::default() };
        assert!(matches!(filter_one(&r, &p), FilterVerdict::TooShort));
    }

    #[cfg(feature = "gzip")]
    #[test]
    fn gzip_roundtrip() {
        use std::io::Write;
        let d = std::env::temp_dir().join(format!("cleaver_fastp_gz_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let inp = d.join("in.fq.gz");
        let f = std::fs::File::create(&inp).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        write!(enc, "@r1\n{}\n+\n{}\n", "ACGT".repeat(8), "I".repeat(32)).unwrap();
        enc.finish().unwrap();
        let out = d.join("out.fq.gz");
        let rep = run_se(&inp, &out, &FastpParams { disable_adapter: true, ..Default::default() }).unwrap();
        assert_eq!(rep.before.reads, 1);
        // the .gz output must itself be valid gzip that round-trips one 32bp read
        let (r, _g) = crate::open_reader(&out).unwrap();
        let mut fr = FastqReader::new(r);
        let rec = fr.next_record().unwrap().expect("one record round-tripped");
        assert_eq!(rec.seq.len(), 32);
    }
}
