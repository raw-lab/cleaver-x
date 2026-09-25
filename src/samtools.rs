//! Pure-Rust re-implementations of the most-used [samtools](https://www.htslib.org)
//! subcommands, built on `noodles`: `view` (flag/MAPQ/region/subsample filtering),
//! `sort`, `index` (BAI), `fastq`/`fasta` extraction (flag-filtered, paired
//! de-interleaving), `flagstat`, `idxstats`, `merge`, `faidx`, and `depth`.
//!
//! Flag semantics match samtools exactly: a record passes when
//! `(flags & require) == require` and `(flags & exclude) == 0`, so e.g.
//! `-f 4` keeps only unmapped reads and `-F 0x900` drops secondary/supplementary.

use std::collections::BTreeMap;
use std::io::{BufRead, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use bstr::BString;
use noodles_bam as bam;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_csi::binning_index::Indexer;
use noodles_sam::header::record::value::map::tag;
use noodles_sam::header::record::value::{map, Map};
use noodles_sam::{self as sam};

use crate::align::{open_align_reader, read_align_header, AlignReader};
use crate::{formats, Format};

// ---------------------------------------------------------------------------
// Flag / MAPQ filter (shared by view, fastq, fasta)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default)]
pub struct Filter {
    pub require: u16, // -f
    pub exclude: u16, // -F
    pub min_mapq: u8, // -q
}

impl Filter {
    #[inline]
    fn keep(&self, flags: u16, mapq: u8) -> bool {
        (flags & self.require) == self.require && (flags & self.exclude) == 0 && mapq >= self.min_mapq
    }
}

fn flag_bits(rec: &sam::alignment::RecordBuf) -> u16 {
    u16::from(rec.flags())
}
fn mapq_of(rec: &sam::alignment::RecordBuf) -> u8 {
    rec.mapping_quality().map(|m| m.get()).unwrap_or(0)
}

// ---------------------------------------------------------------------------
// Writers
// ---------------------------------------------------------------------------

fn align_writer(output: Option<&Path>, fmt: Format) -> Result<Box<dyn sam::alignment::io::Write>> {
    let w: Box<dyn Write> = match output {
        Some(p) => Box::new(std::io::BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating '{}'", p.display()))?,
        )),
        None => Box::new(std::io::BufWriter::new(std::io::stdout())),
    };
    Ok(match fmt {
        Format::Sam => Box::new(sam::io::Writer::new(w)),
        Format::Bam => Box::new(bam::io::Writer::new(w)),
        _ => bail!("alignment output must be SAM or BAM"),
    })
}

fn text_writer(output: Option<&Path>) -> Result<Box<dyn Write>> {
    match output {
        Some(p) => {
            let f = std::fs::File::create(p).with_context(|| format!("creating '{}'", p.display()))?;
            let w = std::io::BufWriter::new(f);
            #[cfg(feature = "gzip")]
            {
                if p.extension().and_then(|e| e.to_str()) == Some("gz") {
                    return Ok(Box::new(flate2::write::GzEncoder::new(w, flate2::Compression::default())));
                }
            }
            Ok(Box::new(w))
        }
        None => Ok(Box::new(std::io::BufWriter::new(std::io::stdout()))),
    }
}

// ---------------------------------------------------------------------------
// Header helpers
// ---------------------------------------------------------------------------

fn ref_name(header: &sam::Header, rid: usize) -> Option<String> {
    header
        .reference_sequences()
        .get_index(rid)
        .and_then(|(n, _)| std::str::from_utf8(n.as_ref()).ok().map(|s| s.to_string()))
}

fn ref_id(header: &sam::Header, name: &str) -> Option<usize> {
    header.reference_sequences().get_index_of(name.as_bytes())
}

fn ref_len(header: &sam::Header, rid: usize) -> u64 {
    header
        .reference_sequences()
        .get_index(rid)
        .map(|(_, m)| usize::from(m.length()) as u64)
        .unwrap_or(0)
}

/// Rewrite the `@HD SO:` sort-order tag (creating an `@HD` line if absent).
fn set_sort_order(header: &mut sam::Header, order: &str) {
    if header.header().is_none() {
        *header.header_mut() = Some(Map::<map::Header>::default());
    }
    if let Some(h) = header.header_mut() {
        if let Ok(so) = tag::Other::try_from([b'S', b'O']) {
            h.other_fields_mut().insert(so, BString::from(order));
        }
    }
}

// ---------------------------------------------------------------------------
// Regions
// ---------------------------------------------------------------------------

struct Region {
    rid: usize,
    start: u64, // 1-based inclusive; 0 = open
    end: u64,   // inclusive; u64::MAX = open
}

fn parse_region(header: &sam::Header, spec: &str) -> Result<Region> {
    let (name, range) = match spec.split_once(':') {
        Some((n, r)) => (n, Some(r)),
        None => (spec, None),
    };
    let rid = ref_id(header, name).with_context(|| format!("reference '{name}' not in header"))?;
    let (mut start, mut end) = (0u64, u64::MAX);
    if let Some(r) = range {
        let r = r.replace(',', "");
        let pos = |s: &str| -> Result<u64> {
            s.parse::<u64>()
                .with_context(|| format!("invalid coordinate '{s}' in region '{spec}'"))
        };
        match r.split_once('-') {
            Some((b, e)) => {
                start = if b.is_empty() { 1 } else { pos(b)? };
                end = if e.is_empty() { u64::MAX } else { pos(e)? };
            }
            None => {
                start = pos(&r)?;
                end = start;
            }
        }
        if start > end {
            bail!("region '{spec}' has start > end");
        }
    }
    Ok(Region { rid, start, end })
}

fn record_overlaps(rec: &sam::alignment::RecordBuf, region: &Region) -> bool {
    if rec.reference_sequence_id() != Some(region.rid) {
        return false;
    }
    let s = rec.alignment_start().map(|p| p.get() as u64).unwrap_or(0);
    let e = rec.alignment_end().map(|p| p.get() as u64).unwrap_or(s);
    s <= region.end && e >= region.start
}

// deterministic subsample keep (samtools -s: hash of read name)
fn keep_subsample(name: &[u8], frac: f64, seed: u64) -> bool {
    if frac >= 1.0 {
        return true;
    }
    let mut h = 0xcbf29ce484222325u64 ^ seed;
    for &b in name {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    ((h >> 11) as f64 / (1u64 << 53) as f64) < frac
}

// ---------------------------------------------------------------------------
// view
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct ViewParams {
    pub filter: Filter,
    pub count: bool,        // -c
    pub header_only: bool,  // -H
    pub no_header: bool,    // (default writes header)
    pub region: Option<String>,
    pub subsample: Option<f64>, // fraction 0..1
    pub seed: u64,
    pub output: Option<PathBuf>,
    pub bam: bool, // -b
}

pub fn view(input: &Path, in_fmt: Format, p: &ViewParams) -> Result<u64> {
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let region = match &p.region {
        Some(s) => Some(parse_region(&header, s)?),
        None => None,
    };
    let out_fmt = if p.bam { Format::Bam } else { Format::Sam };

    if p.count {
        let mut n = 0u64;
        macro_rules! go {
            ($r:expr) => {
                for res in $r.record_bufs(&header) {
                    let rec = res.context("reading record")?;
                    if !p.filter.keep(flag_bits(&rec), mapq_of(&rec)) {
                        continue;
                    }
                    if let Some(reg) = &region {
                        if !record_overlaps(&rec, reg) {
                            continue;
                        }
                    }
                    if let Some(frac) = p.subsample {
                        let name = rec.name().map(|x| x.as_ref()).unwrap_or(b"");
                        if !keep_subsample(name, frac, p.seed) {
                            continue;
                        }
                    }
                    n += 1;
                }
            };
        }
        match &mut reader {
            AlignReader::Sam(r) => go!(r),
            AlignReader::Bam(r) => go!(r),
        }
        println!("{n}");
        return Ok(n);
    }

    let mut writer = align_writer(p.output.as_deref(), out_fmt)?;
    writer.write_alignment_header(&header).context("writing header")?;
    if p.header_only {
        writer.finish(&header).context("finalising")?;
        return Ok(0);
    }
    let mut n = 0u64;
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                if !p.filter.keep(flag_bits(&rec), mapq_of(&rec)) {
                    continue;
                }
                if let Some(reg) = &region {
                    if !record_overlaps(&rec, reg) {
                        continue;
                    }
                }
                if let Some(frac) = p.subsample {
                    let name = rec.name().map(|x| x.as_ref()).unwrap_or(b"");
                    if !keep_subsample(name, frac, p.seed) {
                        continue;
                    }
                }
                writer.write_alignment_record(&header, &rec).context("writing record")?;
                n += 1;
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    writer.finish(&header).context("finalising")?;
    Ok(n)
}

// ---------------------------------------------------------------------------
// sort
// ---------------------------------------------------------------------------

pub fn sort(input: &Path, in_fmt: Format, output: Option<&Path>, by_name: bool, bam: bool) -> Result<u64> {
    let mut reader = open_align_reader(input, in_fmt)?;
    let mut header = read_align_header(&mut reader)?;
    let mut recs: Vec<sam::alignment::RecordBuf> = Vec::new();
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                recs.push(res.context("reading record")?);
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }

    if by_name {
        recs.sort_by(|a, b| {
            let an = a.name().map(|x| x.as_ref()).unwrap_or(b"");
            let bn = b.name().map(|x| x.as_ref()).unwrap_or(b"");
            an.cmp(bn)
        });
    } else {
        // coordinate: by (reference id, position); unplaced (no ref) sort last.
        // Ties are broken exactly as samtools does in bam_sort.c — forward
        // strand before reverse, then QNAME — so the record order matches
        // `samtools sort` byte for byte rather than merely being "a" valid sort.
        recs.sort_by(|a, b| {
            let ak = a.reference_sequence_id().unwrap_or(usize::MAX);
            let bk = b.reference_sequence_id().unwrap_or(usize::MAX);
            let ap = a.alignment_start().map(|p| p.get()).unwrap_or(0);
            let bp = b.alignment_start().map(|p| p.get()).unwrap_or(0);
            let ar = (flag_bits(a) & 0x10) != 0;
            let br = (flag_bits(b) & 0x10) != 0;
            let an = a.name().map(|x| x.as_ref()).unwrap_or(b"");
            let bn = b.name().map(|x| x.as_ref()).unwrap_or(b"");
            ak.cmp(&bk)
                .then(ap.cmp(&bp))
                .then(ar.cmp(&br))
                .then(an.cmp(bn))
        });
    }

    set_sort_order(&mut header, if by_name { "queryname" } else { "coordinate" });

    let out_fmt = if bam { Format::Bam } else { Format::Sam };
    let mut writer = align_writer(output, out_fmt)?;
    writer.write_alignment_header(&header).context("writing header")?;
    for rec in &recs {
        writer.write_alignment_record(&header, rec).context("writing record")?;
    }
    writer.finish(&header).context("finalising")?;
    Ok(recs.len() as u64)
}

// ---------------------------------------------------------------------------
// index (BAI)
// ---------------------------------------------------------------------------

pub fn index(input: &Path, output: Option<&Path>) -> Result<PathBuf> {
    let in_fmt = formats::detect(input)?;
    if in_fmt != Format::Bam {
        bail!("index requires a coordinate-sorted BAM (got {:?})", in_fmt);
    }
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let n_refs = header.reference_sequences().len();

    let r = match &mut reader {
        AlignReader::Bam(r) => r,
        AlignReader::Sam(_) => bail!("index requires BAM"),
    };

    let mut indexer = Indexer::default();
    let mut rec = sam::alignment::RecordBuf::default();
    let mut last_start = r.virtual_position();
    loop {
        let start = last_start;
        let n = r.read_record_buf(&header, &mut rec).context("reading record")?;
        if n == 0 {
            break;
        }
        let end = r.virtual_position();
        last_start = end;
        let chunk = Chunk::new(start, end);
        let ctx = match (rec.reference_sequence_id(), rec.alignment_start(), rec.alignment_end()) {
            (Some(rid), Some(s), Some(e)) => Some((rid, s, e, !rec.flags().is_unmapped())),
            _ => None,
        };
        indexer.add_record(ctx, chunk).context("indexing record")?;
    }
    let idx = indexer.build(n_refs);

    let out = output
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| PathBuf::from(format!("{}.bai", input.display())));
    let file = std::fs::File::create(&out).with_context(|| format!("creating '{}'", out.display()))?;
    let mut w = bam::bai::Writer::new(std::io::BufWriter::new(file));
    w.write_index(&idx).context("writing BAI")?;
    Ok(out)
}

/// Read a BAI back and return its reference-sequence count — used to validate
/// that a freshly written index is well-formed.
pub fn read_bai_ref_count(path: &Path) -> Result<usize> {
    let idx = bam::bai::read(path).with_context(|| format!("reading BAI '{}'", path.display()))?;
    Ok(idx.reference_sequences().len())
}

// ---------------------------------------------------------------------------
// fastq / fasta
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct FastxParams {
    pub filter: Filter,
    pub fasta: bool, // fasta output (drop quality)
    pub out1: Option<PathBuf>,
    pub out2: Option<PathBuf>,
    pub singleton: Option<PathBuf>,
    pub out0: Option<PathBuf>,
    pub output: Option<PathBuf>, // interleaved / single default (stdout)
    pub no_suffix: bool,         // -n
}

#[derive(Debug, Clone, Default)]
pub struct FastxStats {
    pub total: u64,
    pub read1: u64,
    pub read2: u64,
    pub singletons: u64,
    pub single: u64,
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

fn strip_pair_suffix(name: &[u8]) -> &[u8] {
    if name.len() >= 2 && name[name.len() - 2] == b'/' {
        let c = name[name.len() - 1];
        if c == b'1' || c == b'2' {
            return &name[..name.len() - 2];
        }
    }
    name
}

fn emit_read(
    w: &mut dyn Write,
    name: &[u8],
    suffix: &[u8],
    seq: &[u8],
    qual: &[u8],
    reverse: bool,
    fasta: bool,
) -> Result<()> {
    let (s, q): (Vec<u8>, Vec<u8>) = if reverse {
        let s = seq.iter().rev().map(|&b| comp(b)).collect();
        let q = qual.iter().rev().copied().collect();
        (s, q)
    } else {
        (seq.to_vec(), qual.to_vec())
    };
    w.write_all(if fasta { b">" } else { b"@" })?;
    w.write_all(name)?;
    w.write_all(suffix)?;
    w.write_all(b"\n")?;
    w.write_all(&s)?;
    w.write_all(b"\n")?;
    if !fasta {
        w.write_all(b"+\n")?;
        if q.is_empty() {
            // no stored qualities: emit a placeholder of matching length
            let ph = vec![b'I'; s.len()];
            w.write_all(&ph)?;
        } else {
            let ascii: Vec<u8> = q.iter().map(|&x| x.saturating_add(33)).collect();
            w.write_all(&ascii)?;
        }
        w.write_all(b"\n")?;
    }
    Ok(())
}

pub fn fastx(input: &Path, in_fmt: Format, p: &FastxParams) -> Result<FastxStats> {
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;

    let paired_split = p.out1.is_some() && p.out2.is_some();
    let mut w1 = if paired_split { Some(text_writer(p.out1.as_deref())?) } else { None };
    let mut w2 = if paired_split { Some(text_writer(p.out2.as_deref())?) } else { None };
    let mut ws = match &p.singleton {
        Some(_) => Some(text_writer(p.singleton.as_deref())?),
        None => None,
    };
    let mut w0 = match &p.out0 {
        Some(_) => Some(text_writer(p.out0.as_deref())?),
        None => None,
    };
    let mut wsingle = text_writer(p.output.as_deref())?;

    let mut st = FastxStats::default();
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                let bits = flag_bits(&rec);
                if !p.filter.keep(bits, mapq_of(&rec)) {
                    continue;
                }
                st.total += 1;
                let raw = rec.name().map(|x| x.as_ref()).unwrap_or(b"");
                let name = strip_pair_suffix(raw);
                let seq = rec.sequence().as_ref();
                let qual = rec.quality_scores().as_ref();
                let reverse = bits & 0x10 != 0;
                let paired = bits & 0x1 != 0;
                let is_r1 = bits & 0x40 != 0;
                let is_r2 = bits & 0x80 != 0;
                let mate_unmapped = bits & 0x8 != 0;

                let suffix: &[u8] = if p.no_suffix || !paired {
                    b""
                } else if is_r1 {
                    b"/1"
                } else if is_r2 {
                    b"/2"
                } else {
                    b""
                };

                if paired_split && paired {
                    if mate_unmapped && ws.is_some() {
                        st.singletons += 1;
                        emit_read(ws.as_mut().unwrap().as_mut(), name, suffix, seq, qual, reverse, p.fasta)?;
                    } else if is_r1 {
                        st.read1 += 1;
                        emit_read(w1.as_mut().unwrap().as_mut(), name, suffix, seq, qual, reverse, p.fasta)?;
                    } else if is_r2 {
                        st.read2 += 1;
                        emit_read(w2.as_mut().unwrap().as_mut(), name, suffix, seq, qual, reverse, p.fasta)?;
                    } else {
                        st.single += 1;
                        let w = w0.as_mut().map(|w| w.as_mut()).unwrap_or(wsingle.as_mut());
                        emit_read(w, name, suffix, seq, qual, reverse, p.fasta)?;
                    }
                } else if paired_split && !paired {
                    st.single += 1;
                    let w = w0.as_mut().map(|w| w.as_mut()).unwrap_or(wsingle.as_mut());
                    emit_read(w, name, suffix, seq, qual, reverse, p.fasta)?;
                } else {
                    st.single += 1;
                    emit_read(wsingle.as_mut(), name, suffix, seq, qual, reverse, p.fasta)?;
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    wsingle.flush()?;
    for w in [&mut w1, &mut w2, &mut ws, &mut w0].into_iter().flatten() {
        w.flush()?;
    }
    Ok(st)
}

// ---------------------------------------------------------------------------
// flagstat
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct FlagStat {
    pub total: [u64; 2],
    pub primary: [u64; 2],
    pub secondary: [u64; 2],
    pub supplementary: [u64; 2],
    pub duplicates: [u64; 2],
    pub primary_duplicates: [u64; 2],
    pub mapped: [u64; 2],
    pub primary_mapped: [u64; 2],
    pub paired: [u64; 2],
    pub read1: [u64; 2],
    pub read2: [u64; 2],
    pub properly_paired: [u64; 2],
    pub with_mate_mapped: [u64; 2],
    pub singletons: [u64; 2],
    pub mate_diff_chr: [u64; 2],
    pub mate_diff_chr_hq: [u64; 2],
}

pub fn flagstat(input: &Path, in_fmt: Format) -> Result<FlagStat> {
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let mut fs = FlagStat::default();
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                let b = flag_bits(&rec);
                let q = (b & 0x200 != 0) as usize; // 0 = QC pass, 1 = QC fail
                fs.total[q] += 1;
                let secondary = b & 0x100 != 0;
                let supplementary = b & 0x800 != 0;
                let primary = !secondary && !supplementary;
                let mapped = b & 0x4 == 0;
                let dup = b & 0x400 != 0;
                if secondary {
                    fs.secondary[q] += 1;
                }
                if supplementary {
                    fs.supplementary[q] += 1;
                }
                if primary {
                    fs.primary[q] += 1;
                }
                if dup {
                    fs.duplicates[q] += 1;
                    if primary {
                        fs.primary_duplicates[q] += 1;
                    }
                }
                if mapped {
                    fs.mapped[q] += 1;
                    if primary {
                        fs.primary_mapped[q] += 1;
                    }
                }
                if b & 0x1 != 0 {
                    fs.paired[q] += 1;
                    if b & 0x40 != 0 {
                        fs.read1[q] += 1;
                    }
                    if b & 0x80 != 0 {
                        fs.read2[q] += 1;
                    }
                    if b & 0x2 != 0 {
                        fs.properly_paired[q] += 1;
                    }
                    if mapped && b & 0x8 == 0 {
                        fs.with_mate_mapped[q] += 1;
                        if rec.reference_sequence_id() != rec.mate_reference_sequence_id() {
                            fs.mate_diff_chr[q] += 1;
                            if mapq_of(&rec) >= 5 {
                                fs.mate_diff_chr_hq[q] += 1;
                            }
                        }
                    }
                    if mapped && b & 0x8 != 0 {
                        fs.singletons[q] += 1;
                    }
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    Ok(fs)
}

impl FlagStat {
    /// samtools-format report.
    pub fn report(&self) -> String {
        // samtools prints a bare "N/A" (no percent sign) when the denominator
        // is zero, and "12.34%" otherwise — so the suffix belongs here, not at
        // the call sites.
        let pct = |n: u64, d: u64| -> String {
            if d == 0 {
                "N/A".to_string()
            } else {
                format!("{:.2}%", n as f64 / d as f64 * 100.0)
            }
        };
        let line = |label: &str, c: &[u64; 2]| format!("{} + {} {}\n", c[0], c[1], label);
        let mut s = String::new();
        s.push_str(&line("in total (QC-passed reads + QC-failed reads)", &self.total));
        s.push_str(&line("primary", &self.primary));
        s.push_str(&line("secondary", &self.secondary));
        s.push_str(&line("supplementary", &self.supplementary));
        s.push_str(&line("duplicates", &self.duplicates));
        s.push_str(&line("primary duplicates", &self.primary_duplicates));
        s.push_str(&format!(
            "{} + {} mapped ({} : {})\n",
            self.mapped[0], self.mapped[1],
            pct(self.mapped[0], self.total[0]), pct(self.mapped[1], self.total[1])
        ));
        s.push_str(&format!(
            "{} + {} primary mapped ({} : {})\n",
            self.primary_mapped[0], self.primary_mapped[1],
            pct(self.primary_mapped[0], self.primary[0]), pct(self.primary_mapped[1], self.primary[1])
        ));
        s.push_str(&line("paired in sequencing", &self.paired));
        s.push_str(&line("read1", &self.read1));
        s.push_str(&line("read2", &self.read2));
        s.push_str(&format!(
            "{} + {} properly paired ({} : {})\n",
            self.properly_paired[0], self.properly_paired[1],
            pct(self.properly_paired[0], self.paired[0]), pct(self.properly_paired[1], self.paired[1])
        ));
        s.push_str(&line("with itself and mate mapped", &self.with_mate_mapped));
        s.push_str(&format!(
            "{} + {} singletons ({} : {})\n",
            self.singletons[0], self.singletons[1],
            pct(self.singletons[0], self.paired[0]), pct(self.singletons[1], self.paired[1])
        ));
        s.push_str(&line("with mate mapped to a different chr", &self.mate_diff_chr));
        s.push_str(&line("with mate mapped to a different chr (mapQ>=5)", &self.mate_diff_chr_hq));
        s
    }
}

// ---------------------------------------------------------------------------
// idxstats
// ---------------------------------------------------------------------------

pub fn idxstats(input: &Path, in_fmt: Format) -> Result<Vec<(String, u64, u64, u64)>> {
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let n = header.reference_sequences().len();
    let mut mapped = vec![0u64; n];
    let mut unmapped = vec![0u64; n];
    let mut star = 0u64;
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                let b = flag_bits(&rec);
                match rec.reference_sequence_id() {
                    Some(rid) => {
                        if b & 0x4 != 0 {
                            unmapped[rid] += 1;
                        } else {
                            mapped[rid] += 1;
                        }
                    }
                    None => star += 1,
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    let mut out = Vec::with_capacity(n + 1);
    for rid in 0..n {
        out.push((ref_name(&header, rid).unwrap_or_default(), ref_len(&header, rid), mapped[rid], unmapped[rid]));
    }
    out.push(("*".to_string(), 0, 0, star));
    Ok(out)
}

// ---------------------------------------------------------------------------
// merge (concat + coordinate sort; uses the first input's header)
// ---------------------------------------------------------------------------

pub fn merge(inputs: &[PathBuf], output: Option<&Path>, bam: bool) -> Result<u64> {
    if inputs.is_empty() {
        bail!("merge needs at least one input");
    }
    let first_fmt = formats::detect(&inputs[0])?;
    let mut base_reader = open_align_reader(&inputs[0], first_fmt)?;
    let mut header = read_align_header(&mut base_reader)?;
    drop(base_reader);

    let mut recs: Vec<sam::alignment::RecordBuf> = Vec::new();
    for inp in inputs {
        let fmt = formats::detect(inp)?;
        let mut reader = open_align_reader(inp, fmt)?;
        let h = read_align_header(&mut reader)?;
        macro_rules! go {
            ($r:expr) => {
                for res in $r.record_bufs(&h) {
                    recs.push(res.context("reading record")?);
                }
            };
        }
        match &mut reader {
            AlignReader::Sam(r) => go!(r),
            AlignReader::Bam(r) => go!(r),
        }
    }
    recs.sort_by(|a, b| {
        let ak = a.reference_sequence_id().unwrap_or(usize::MAX);
        let bk = b.reference_sequence_id().unwrap_or(usize::MAX);
        let ap = a.alignment_start().map(|p| p.get()).unwrap_or(0);
        let bp = b.alignment_start().map(|p| p.get()).unwrap_or(0);
        ak.cmp(&bk).then(ap.cmp(&bp))
    });

    set_sort_order(&mut header, "coordinate");

    let out_fmt = if bam { Format::Bam } else { Format::Sam };
    let mut writer = align_writer(output, out_fmt)?;
    writer.write_alignment_header(&header).context("writing header")?;
    for rec in &recs {
        writer.write_alignment_record(&header, rec).context("writing record")?;
    }
    writer.finish(&header).context("finalising")?;
    Ok(recs.len() as u64)
}

// ---------------------------------------------------------------------------
// depth (samtools-compatible: rows for every spanned position; D/N show 0)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct DepthParams {
    pub min_mapq: u8, // -Q
    pub region: Option<String>,
    pub output: Option<PathBuf>,
}

pub fn depth(input: &Path, in_fmt: Format, p: &DepthParams) -> Result<()> {
    use sam::alignment::record::cigar::op::Kind;
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let region = match &p.region {
        Some(s) => Some(parse_region(&header, s)?),
        None => None,
    };
    // depth per (reference id, position)
    let mut cov: BTreeMap<(usize, u64), u32> = BTreeMap::new();
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                let b = flag_bits(&rec);
                // samtools depth's default --excl-flags is UNMAP,SECONDARY,
                // QCFAIL,DUP (0x704). Supplementary alignments are NOT excluded.
                if b & 0x704 != 0 {
                    continue;
                }
                if mapq_of(&rec) < p.min_mapq {
                    continue;
                }
                let rid = match rec.reference_sequence_id() {
                    Some(r) => r,
                    None => continue,
                };
                if let Some(reg) = &region {
                    if rid != reg.rid {
                        continue;
                    }
                }
                let start = match rec.alignment_start() {
                    Some(s) => s.get() as u64,
                    None => continue,
                };
                let ops: Vec<_> = rec.cigar().as_ref().to_vec();
                let in_region = |pos: u64| -> bool {
                    match &region {
                        Some(reg) => pos >= reg.start && pos <= reg.end,
                        None => true,
                    }
                };
                // samtools reports a row for every position SPANNED by the
                // alignment — including deletion (D) and reference-skip (N)
                // gaps, which show depth 0 — while only M/=/X bases add depth.
                let mut refpos = start;
                for op in &ops {
                    let l = op.len() as u64;
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            for k in 0..l {
                                let pos = refpos + k;
                                if in_region(pos) {
                                    *cov.entry((rid, pos)).or_insert(0) += 1;
                                }
                            }
                            refpos += l;
                        }
                        Kind::Deletion | Kind::Skip => {
                            for k in 0..l {
                                let pos = refpos + k;
                                if in_region(pos) {
                                    cov.entry((rid, pos)).or_insert(0);
                                }
                            }
                            refpos += l;
                        }
                        _ => {}
                    }
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }
    let mut w = text_writer(p.output.as_deref())?;
    for ((rid, pos), d) in &cov {
        let name = ref_name(&header, *rid).unwrap_or_default();
        writeln!(w, "{name}\t{pos}\t{d}")?;
    }
    w.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// faidx (FASTA index + region extraction)
// ---------------------------------------------------------------------------

struct FaiEntry {
    length: u64,
    offset: u64,    // byte offset of the first base
    linebases: u64, // bases per line
    linewidth: u64, // bytes per line (incl newline)
}

fn build_fai(fasta: &Path) -> Result<BTreeMap<String, FaiEntry>> {
    let file = std::fs::File::open(fasta).with_context(|| format!("opening '{}'", fasta.display()))?;
    let mut reader = std::io::BufReader::new(file);
    let mut map: BTreeMap<String, FaiEntry> = BTreeMap::new();
    let mut order: Vec<String> = Vec::new();
    let mut cur: Option<(String, u64, u64, u64, u64)> = None; // name, offset, len, linebases, linewidth
    let mut offset = 0u64;
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            break;
        }
        let raw_len = line.len() as u64;
        let trimmed_len = {
            let mut l = line.len();
            while l > 0 && (line[l - 1] == b'\n' || line[l - 1] == b'\r') {
                l -= 1;
            }
            l as u64
        };
        if line.first() == Some(&b'>') {
            if let Some((name, off, len, lb, lw)) = cur.take() {
                map.insert(name.clone(), FaiEntry { length: len, offset: off, linebases: lb, linewidth: lw });
                order.push(name);
            }
            let name = std::str::from_utf8(&line[1..trimmed_len as usize + 1 - 1])
                .unwrap_or("")
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            cur = Some((name, offset + raw_len, 0, 0, 0));
        } else if let Some((_, _, len, lb, lw)) = cur.as_mut() {
            if *lb == 0 {
                *lb = trimmed_len;
                *lw = raw_len;
            }
            *len += trimmed_len;
        }
        offset += raw_len;
    }
    if let Some((name, off, len, lb, lw)) = cur.take() {
        map.insert(name.clone(), FaiEntry { length: len, offset: off, linebases: lb, linewidth: lw });
        order.push(name);
    }
    let _ = order;
    Ok(map)
}

pub fn faidx(fasta: &Path, regions: &[String], output: Option<&Path>) -> Result<()> {
    let fai = build_fai(fasta)?;
    if regions.is_empty() {
        // write <fasta>.fai
        let out = PathBuf::from(format!("{}.fai", fasta.display()));
        let mut w = std::io::BufWriter::new(
            std::fs::File::create(&out).with_context(|| format!("creating '{}'", out.display()))?,
        );
        for (name, e) in &fai {
            writeln!(w, "{}\t{}\t{}\t{}\t{}", name, e.length, e.offset, e.linebases, e.linewidth)?;
        }
        w.flush()?;
        eprintln!("wrote {}", out.display());
        return Ok(());
    }
    // extract each region
    let mut file = std::fs::File::open(fasta)?;
    let mut w = text_writer(output)?;
    for spec in regions {
        let (name, range) = match spec.split_once(':') {
            Some((n, r)) => (n, Some(r)),
            None => (spec.as_str(), None),
        };
        let e = fai
            .get(name)
            .with_context(|| format!("sequence '{name}' not in FASTA index"))?;
        let (beg, end) = match range {
            Some(r) => {
                let r = r.replace(',', "");
                let pos = |s: &str| -> Result<u64> {
                    s.parse::<u64>()
                        .with_context(|| format!("invalid coordinate '{s}' in region '{spec}'"))
                };
                match r.split_once('-') {
                    Some((b, e2)) => (
                        if b.is_empty() { 1 } else { pos(b)?.max(1) },
                        if e2.is_empty() { e.length } else { pos(e2)? },
                    ),
                    None => (pos(&r)?.max(1), e.length),
                }
            }
            None => (1, e.length),
        };
        let end = end.min(e.length);
        if beg > end || e.linebases == 0 {
            writeln!(w, ">{spec}")?;
            continue;
        }
        // byte offset of base `beg` (1-based): full lines + remainder
        let idx0 = beg - 1;
        let byte = e.offset + (idx0 / e.linebases) * e.linewidth + (idx0 % e.linebases);
        let want = (end - beg + 1) as usize;
        file.seek(SeekFrom::Start(byte))?;
        let mut buf = Vec::new();
        let mut collected = 0usize;
        // read enough bytes (including newlines) to cover `want` bases
        let mut chunk = vec![0u8; want + want / e.linebases as usize + 16];
        let got = file.read(&mut chunk)?;
        for &c in &chunk[..got] {
            if c == b'\n' || c == b'\r' {
                continue;
            }
            buf.push(c);
            collected += 1;
            if collected >= want {
                break;
            }
        }
        writeln!(w, ">{}:{}-{}", name, beg, end)?;
        // wrap at 60 like samtools
        for line in buf.chunks(60) {
            w.write_all(line)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// mpileup / pileup (text pileup) + coverage
// ---------------------------------------------------------------------------

/// Load a FASTA into memory keyed by sequence name (first whitespace token).
fn load_reference(path: &Path) -> Result<std::collections::HashMap<String, Vec<u8>>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening '{}'", path.display()))?;
    let mut reader = std::io::BufReader::new(file);
    let mut map = std::collections::HashMap::new();
    let mut name = String::new();
    let mut seq: Vec<u8> = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            if !name.is_empty() {
                map.insert(std::mem::take(&mut name), std::mem::take(&mut seq));
            }
            name = String::from_utf8_lossy(&line[1..])
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
        } else {
            seq.extend_from_slice(&line);
        }
    }
    if !name.is_empty() {
        map.insert(name, seq);
    }
    Ok(map)
}

#[derive(Default)]
struct Column {
    bases: String,
    quals: String,
    depth: u32,
}

#[derive(Debug, Clone, Default)]
pub struct MpileupParams {
    pub min_mapq: u8,  // -q
    pub min_baseq: u8, // -Q
    pub region: Option<String>,
    pub reference: Option<PathBuf>, // -f
    pub output: Option<PathBuf>,
}

/// samtools-style text pileup. With a reference (`-f`), matches show as `.`/`,`
/// and mismatches as the base letter; without one the reference column is `N`
/// and every base is shown as a letter. Read starts (`^`+mapQ), read ends (`$`),
/// insertions (`+N…`), and deletions (`-N…` plus `*` placeholders) are emitted.
pub fn mpileup(input: &Path, in_fmt: Format, p: &MpileupParams) -> Result<()> {
    use sam::alignment::record::cigar::op::Kind;

    let refseq = match &p.reference {
        Some(f) => Some(load_reference(f)?),
        None => None,
    };
    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let region = match &p.region {
        Some(s) => Some(parse_region(&header, s)?),
        None => None,
    };

    let mut cols: BTreeMap<(usize, u64), Column> = BTreeMap::new();
    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                let b = flag_bits(&rec);
                // samtools mpileup's default --excl-flags is UNMAP,SECONDARY,
                // QCFAIL,DUP (0x704); supplementary alignments are kept.
                if b & 0x704 != 0 {
                    continue;
                }
                let mapq = mapq_of(&rec);
                if mapq < p.min_mapq {
                    continue;
                }
                let rid = match rec.reference_sequence_id() {
                    Some(r) => r,
                    None => continue,
                };
                if let Some(reg) = &region {
                    if rid != reg.rid {
                        continue;
                    }
                }
                let start = match rec.alignment_start() {
                    Some(s) => s.get() as u64,
                    None => continue,
                };
                let reverse = b & 0x10 != 0;
                let seq = rec.sequence().as_ref().to_vec();
                let qual = rec.quality_scores().as_ref().to_vec();
                let ops: Vec<_> = rec.cigar().as_ref().to_vec();
                let rname = ref_name(&header, rid).unwrap_or_default();
                let rbytes = refseq.as_ref().and_then(|m| m.get(&rname));

                let mapq_char = (mapq.min(93) + 33) as char;
                let mut refpos = start;
                let mut qidx = 0usize;
                let mut first = true;
                let mut prev_pos = 0u64;

                for op in &ops {
                    let l = op.len() as usize;
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            for k in 0..l {
                                let pos = refpos + k as u64;
                                let qb = seq.get(qidx + k).copied().unwrap_or(b'N');
                                let bq = qual.get(qidx + k).copied().unwrap_or(0);
                                if bq < p.min_baseq {
                                    // filtered base: not shown, so it neither
                                    // carries the read-start marker nor advances
                                    // the indel anchor
                                    continue;
                                }
                                let in_region = region
                                    .as_ref()
                                    .map(|r| pos >= r.start && pos <= r.end)
                                    .unwrap_or(true);
                                if !in_region {
                                    continue;
                                }
                                let refb = rbytes.and_then(|s| s.get((pos - 1) as usize)).copied();
                                let ch = match refb {
                                    Some(r) if qb.eq_ignore_ascii_case(&r) => {
                                        if reverse { ',' } else { '.' }
                                    }
                                    _ => {
                                        if reverse {
                                            qb.to_ascii_lowercase() as char
                                        } else {
                                            qb.to_ascii_uppercase() as char
                                        }
                                    }
                                };
                                let col = cols.entry((rid, pos)).or_default();
                                if first {
                                    col.bases.push('^');
                                    col.bases.push(mapq_char);
                                }
                                col.bases.push(ch);
                                col.quals.push((bq.min(93) + 33) as char);
                                col.depth += 1;
                                prev_pos = pos;
                                first = false;
                            }
                            refpos += l as u64;
                            qidx += l;
                        }
                        Kind::Insertion => {
                            if prev_pos > 0 {
                                if let Some(col) = cols.get_mut(&(rid, prev_pos)) {
                                    let ins: String = (0..l)
                                        .map(|k| {
                                            let q = seq.get(qidx + k).copied().unwrap_or(b'N');
                                            (if reverse { q.to_ascii_lowercase() } else { q.to_ascii_uppercase() }) as char
                                        })
                                        .collect();
                                    col.bases.push_str(&format!("+{l}{ins}"));
                                }
                            }
                            qidx += l;
                        }
                        Kind::Deletion => {
                            if prev_pos > 0 {
                                if let Some(col) = cols.get_mut(&(rid, prev_pos)) {
                                    let del: String = (0..l)
                                        .map(|k| {
                                            let r = rbytes
                                                .and_then(|s| s.get((refpos + k as u64 - 1) as usize))
                                                .copied()
                                                .unwrap_or(b'N');
                                            // samtools prints the deleted
                                            // reference bases in the read's
                                            // strand case, as it does for bases.
                                            (if reverse {
                                                r.to_ascii_lowercase()
                                            } else {
                                                r.to_ascii_uppercase()
                                            }) as char
                                        })
                                        .collect();
                                    col.bases.push_str(&format!("-{l}{del}"));
                                }
                            }
                            for k in 0..l {
                                let dp = refpos + k as u64;
                                let in_region = region
                                    .as_ref()
                                    .map(|r| dp >= r.start && dp <= r.end)
                                    .unwrap_or(true);
                                if in_region {
                                    // samtools gives the '*' placeholder the
                                    // quality of the base at the junction (qidx
                                    // is not advanced by a deletion).
                                    let bq = qual.get(qidx).copied().unwrap_or(0);
                                    let col = cols.entry((rid, dp)).or_default();
                                    col.bases.push('*');
                                    col.quals.push((bq.min(93) + 33) as char);
                                    col.depth += 1;
                                }
                            }
                            refpos += l as u64;
                        }
                        Kind::Skip => {
                            // samtools renders an N (reference skip) as '>' on
                            // the forward strand and '<' on the reverse, counted
                            // in the depth, with the quality of the base at the
                            // junction (qidx is not advanced by an N).
                            for k in 0..l {
                                let pos = refpos + k as u64;
                                let in_region = region
                                    .as_ref()
                                    .map(|r| pos >= r.start && pos <= r.end)
                                    .unwrap_or(true);
                                if !in_region {
                                    continue;
                                }
                                let bq = qual.get(qidx).copied().unwrap_or(0);
                                let col = cols.entry((rid, pos)).or_default();
                                col.bases.push(if reverse { '<' } else { '>' });
                                col.quals.push((bq.min(93) + 33) as char);
                                col.depth += 1;
                            }
                            refpos += l as u64;
                        }
                        Kind::SoftClip => qidx += l,
                        Kind::HardClip | Kind::Pad => {}
                    }
                }

                // end-of-read marker on the last base this read actually emitted
                if prev_pos > 0 {
                    if let Some(col) = cols.get_mut(&(rid, prev_pos)) {
                        col.bases.push('$');
                    }
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }

    let mut w = text_writer(p.output.as_deref())?;
    for ((rid, pos), col) in &cols {
        let rname = ref_name(&header, *rid).unwrap_or_default();
        let refb = refseq
            .as_ref()
            .and_then(|m| m.get(&rname))
            .and_then(|s| s.get((*pos - 1) as usize))
            .map(|&b| (b as char).to_ascii_uppercase())
            .unwrap_or('N');
        writeln!(w, "{rname}\t{pos}\t{refb}\t{}\t{}\t{}", col.depth, col.bases, col.quals)?;
    }
    w.flush()?;
    Ok(())
}

/// `pileup` is the modern samtools alias for the text `mpileup`.
pub fn pileup(input: &Path, in_fmt: Format, p: &MpileupParams) -> Result<()> {
    mpileup(input, in_fmt, p)
}

// ---------------------------------------------------------------------------
// coverage (per-reference summary)
// ---------------------------------------------------------------------------

/// C `printf("%g")` formatting: `sig` significant digits, trailing zeros and a
/// trailing point removed, switching to scientific notation outside
/// `10^-4 .. 10^sig`. `samtools coverage` formats its float columns this way, so
/// matching it keeps our output byte-identical.
fn fmt_g(v: f64, sig: usize) -> String {
    if !v.is_finite() {
        return format!("{v}");
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let exp = v.abs().log10().floor() as i32;
    if exp < -4 || exp >= sig as i32 {
        let mant_prec = sig.saturating_sub(1);
        let s = format!("{:.*e}", mant_prec, v);
        // trim trailing zeros in the mantissa (e.g. 1.20000e5 -> 1.2e5)
        match s.split_once('e') {
            Some((m, e)) => {
                let m = if m.contains('.') {
                    m.trim_end_matches('0').trim_end_matches('.')
                } else {
                    m
                };
                format!("{m}e{e}")
            }
            None => s,
        }
    } else {
        let prec = (sig as i32 - 1 - exp).max(0) as usize;
        let s = format!("{:.*}", prec, v);
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct CoverageParams {
    pub region: Option<String>,
    pub output: Option<PathBuf>,
}

pub fn coverage(input: &Path, in_fmt: Format, p: &CoverageParams) -> Result<()> {
    use sam::alignment::record::cigar::op::Kind;

    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let region = match &p.region {
        Some(s) => Some(parse_region(&header, s)?),
        None => None,
    };
    let n = header.reference_sequences().len();
    let mut numreads = vec![0u64; n];
    let mut sum_mapq = vec![0u64; n];
    let mut sum_baseq = vec![0u64; n];
    let mut aligned_bases = vec![0u64; n];
    let mut covered: BTreeMap<(usize, u64), ()> = BTreeMap::new();

    macro_rules! go {
        ($r:expr) => {
            for res in $r.record_bufs(&header) {
                let rec = res.context("reading record")?;
                let b = flag_bits(&rec);
                // exclude unmapped/secondary/qc-fail/duplicate (samtools coverage default)
                if b & 0x704 != 0 {
                    continue;
                }
                let rid = match rec.reference_sequence_id() {
                    Some(r) => r,
                    None => continue,
                };
                if let Some(reg) = &region {
                    if rid != reg.rid {
                        continue;
                    }
                }
                let start = match rec.alignment_start() {
                    Some(s) => s.get() as u64,
                    None => continue,
                };
                numreads[rid] += 1;
                sum_mapq[rid] += mapq_of(&rec) as u64;
                let qual = rec.quality_scores().as_ref().to_vec();
                let ops: Vec<_> = rec.cigar().as_ref().to_vec();
                let mut refpos = start;
                let mut qidx = 0usize;
                for op in &ops {
                    let l = op.len() as usize;
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            for k in 0..l {
                                let pos = refpos + k as u64;
                                let in_region = region
                                    .as_ref()
                                    .map(|r| pos >= r.start && pos <= r.end)
                                    .unwrap_or(true);
                                if in_region {
                                    aligned_bases[rid] += 1;
                                    sum_baseq[rid] += qual.get(qidx + k).copied().unwrap_or(0) as u64;
                                    covered.insert((rid, pos), ());
                                }
                            }
                            refpos += l as u64;
                            qidx += l;
                        }
                        Kind::Deletion | Kind::Skip => refpos += l as u64,
                        Kind::SoftClip | Kind::Insertion => qidx += l,
                        Kind::HardClip | Kind::Pad => {}
                    }
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => go!(r),
        AlignReader::Bam(r) => go!(r),
    }

    let mut covbases = vec![0u64; n];
    for (rid, _) in covered.keys() {
        covbases[*rid] += 1;
    }

    let mut w = text_writer(p.output.as_deref())?;
    writeln!(w, "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq")?;
    for rid in 0..n {
        if let Some(reg) = &region {
            if rid != reg.rid {
                continue;
            }
        }
        let name = ref_name(&header, rid).unwrap_or_default();
        let len = ref_len(&header, rid).max(1);
        let cov_pct = covbases[rid] as f64 / len as f64 * 100.0;
        let meandepth = aligned_bases[rid] as f64 / len as f64;
        let meanbaseq = if aligned_bases[rid] > 0 {
            sum_baseq[rid] as f64 / aligned_bases[rid] as f64
        } else {
            0.0
        };
        let meanmapq = if numreads[rid] > 0 {
            sum_mapq[rid] as f64 / numreads[rid] as f64
        } else {
            0.0
        };
        // samtools prints coverage/meandepth with %g (6 significant digits)
        // and the quality columns with %.3g.
        writeln!(
            w,
            "{name}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            ref_len(&header, rid),
            numreads[rid],
            covbases[rid],
            fmt_g(cov_pct, 6),
            fmt_g(meandepth, 6),
            fmt_g(meanbaseq, 3),
            fmt_g(meanmapq, 3)
        )?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn workdir() -> PathBuf {
        let d = std::env::temp_dir().join(format!("cleaver-samtools-{}", std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        d
    }

    fn write_sam(path: &Path) {
        // 4 records: 2 mapped (chr1), 1 unmapped, 1 mapped chr2
        let sam = "\
@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:1000
r1\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII
r2\t0\tchr1\t50\t60\t10M\t*\t0\t0\tTTTTTTTTTT\tIIIIIIIIII
r3\t4\t*\t0\t0\t*\t*\t0\t0\tGGGGGGGGGG\tIIIIIIIIII
r4\t16\tchr2\t200\t30\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII
";
        let mut f = std::fs::File::create(path).unwrap();
        f.write_all(sam.as_bytes()).unwrap();
    }

    #[test]
    fn filter_flag_semantics() {
        let f = Filter { require: 0x4, exclude: 0, min_mapq: 0 };
        assert!(f.keep(0x4, 0)); // unmapped kept with -f 4
        assert!(!f.keep(0x0, 60)); // mapped rejected
        let g = Filter { require: 0, exclude: 0x904, min_mapq: 0 };
        assert!(g.keep(0x0, 0));
        assert!(!g.keep(0x100, 0)); // secondary excluded
    }

    #[test]
    fn view_count_and_filter() {
        let d = workdir();
        let sam = d.join("a.sam");
        write_sam(&sam);
        // count unmapped (-f 4) => 1
        let p = ViewParams { filter: Filter { require: 0x4, ..Default::default() }, count: true, ..Default::default() };
        assert_eq!(view(&sam, Format::Sam, &p).unwrap(), 1);
        // count mapped (-F 4) => 3
        let p2 = ViewParams { filter: Filter { exclude: 0x4, ..Default::default() }, count: true, ..Default::default() };
        assert_eq!(view(&sam, Format::Sam, &p2).unwrap(), 3);
    }

    #[test]
    fn sort_then_index_roundtrip() {
        let d = workdir();
        let sam = d.join("s.sam");
        write_sam(&sam);
        let sorted = d.join("s.sorted.bam");
        let n = sort(&sam, Format::Sam, Some(&sorted), false, true).unwrap();
        assert_eq!(n, 4);
        // index the sorted BAM
        let bai = index(&sorted, None).unwrap();
        assert!(bai.exists());
        assert!(std::fs::metadata(&bai).unwrap().len() > 0);
    }

    #[test]
    fn fastq_extracts_unmapped() {
        let d = workdir();
        let sam = d.join("f.sam");
        write_sam(&sam);
        let out = d.join("unmapped.fastq");
        let p = FastxParams {
            filter: Filter { require: 0x4, ..Default::default() },
            output: Some(out.clone()),
            ..Default::default()
        };
        let st = fastx(&sam, Format::Sam, &p).unwrap();
        assert_eq!(st.total, 1);
        let s = std::fs::read_to_string(&out).unwrap();
        assert!(s.contains("@r3"));
        assert!(s.contains("GGGGGGGGGG"));
    }

    #[test]
    fn flagstat_counts() {
        let d = workdir();
        let sam = d.join("fs.sam");
        write_sam(&sam);
        let fs = flagstat(&sam, Format::Sam).unwrap();
        assert_eq!(fs.total[0], 4);
        assert_eq!(fs.mapped[0], 3);
    }

    #[test]
    fn idxstats_per_ref() {
        let d = workdir();
        let sam = d.join("ix.sam");
        write_sam(&sam);
        let rows = idxstats(&sam, Format::Sam).unwrap();
        // chr1: 2 mapped, chr2: 1 mapped, *: 1 unmapped
        assert_eq!(rows[0], ("chr1".to_string(), 1000, 2, 0));
        assert_eq!(rows[1], ("chr2".to_string(), 1000, 1, 0));
        assert_eq!(rows[2], ("*".to_string(), 0, 0, 1));
    }

    #[test]
    fn faidx_index_and_extract() {
        let d = workdir();
        let fa = d.join("g.fa");
        std::fs::write(&fa, b">chr1\nACGTACGTAC\nGGGGCCCCTT\n>chr2\nTTTTAAAACC\n").unwrap();
        faidx(&fa, &[], None).unwrap();
        let fai = std::fs::read_to_string(d.join("g.fa.fai")).unwrap();
        assert!(fai.contains("chr1\t20\t"));
        // extract chr1:11-14 -> GGGG
        let out = d.join("sub.fa");
        faidx(&fa, &["chr1:11-14".to_string()], Some(&out)).unwrap();
        let s = std::fs::read_to_string(&out).unwrap();
        assert!(s.contains("GGGG"), "got: {s}");
    }

    #[test]
    fn sort_sets_so_tag() {
        let d = workdir();
        let sam = d.join("so.sam");
        write_sam(&sam);
        // coordinate sort -> SAM out, header should carry SO:coordinate
        let out = d.join("so.sorted.sam");
        sort(&sam, Format::Sam, Some(&out), false, false).unwrap();
        let s = std::fs::read_to_string(&out).unwrap();
        assert!(s.contains("SO:coordinate"), "no SO tag: {}", &s[..s.find('\n').unwrap_or(s.len())]);
        // name sort -> SO:queryname
        let outn = d.join("so.byname.sam");
        sort(&sam, Format::Sam, Some(&outn), true, false).unwrap();
        assert!(std::fs::read_to_string(&outn).unwrap().contains("SO:queryname"));
    }

    #[test]
    fn index_is_valid_bai() {
        let d = workdir();
        let sam = d.join("v.sam");
        write_sam(&sam);
        let bam = d.join("v.sorted.bam");
        sort(&sam, Format::Sam, Some(&bam), false, true).unwrap();
        let bai = index(&bam, None).unwrap();
        // read the index back with noodles to prove it is a well-formed BAI
        let idx = bam::bai::read(&bai).expect("BAI should parse");
        assert_eq!(idx.reference_sequences().len(), 2);
    }

    #[test]
    fn fastq_from_bam_mapped_and_unmapped() {
        let d = workdir();
        let sam = d.join("b.sam");
        write_sam(&sam);
        // make a BAM
        let bam = d.join("b.bam");
        sort(&sam, Format::Sam, Some(&bam), false, true).unwrap();
        // unmapped (-f 4) from BAM
        let un = d.join("un.fq");
        let su = fastx(
            &bam,
            Format::Bam,
            &FastxParams { filter: Filter { require: 4, ..Default::default() }, output: Some(un.clone()), ..Default::default() },
        )
        .unwrap();
        assert_eq!(su.total, 1);
        assert!(std::fs::read_to_string(&un).unwrap().contains("GGGGGGGGGG"));
        // mapped (-F 4) from BAM
        let mp = d.join("mp.fq");
        let sm = fastx(
            &bam,
            Format::Bam,
            &FastxParams { filter: Filter { exclude: 4, ..Default::default() }, output: Some(mp.clone()), ..Default::default() },
        )
        .unwrap();
        assert_eq!(sm.total, 3);
        let s = std::fs::read_to_string(&mp).unwrap();
        assert!(s.contains("@r1") && s.contains("@r2") && s.contains("@r4"));
    }

    #[test]
    fn mpileup_matches_reference() {
        let d = workdir();
        // reference chr1: positions 10..14 (1-based) = ACGTA
        let fa = d.join("ref.fa");
        std::fs::write(&fa, b">chr1\nNNNNNNNNNACGTANNNNNN\n").unwrap();
        let sam = d.join("m.sam");
        std::fs::write(
            &sam,
            b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:20\nrd\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n",
        )
        .unwrap();
        let out = d.join("m.pileup");
        mpileup(
            &sam,
            Format::Sam,
            &MpileupParams { reference: Some(fa), output: Some(out.clone()), ..Default::default() },
        )
        .unwrap();
        let s = std::fs::read_to_string(&out).unwrap();
        // 5 covered positions, all matches (.), read start ^ at 10 and end $ at 14
        assert_eq!(s.lines().count(), 5);
        assert!(s.contains("chr1\t12\tG\t1"), "pos12: {s}");
        assert!(s.contains('^') && s.contains('$'), "markers missing: {s}");
        assert!(s.lines().nth(2).unwrap().contains('.'), "match dot missing: {s}");
    }

    #[test]
    fn coverage_summary() {
        let d = workdir();
        let sam = d.join("c.sam");
        // chr1: two 10bp reads at 50 and 100 (20 covered bases of 1000)
        write_sam(&sam);
        let out = d.join("c.cov");
        coverage(&sam, Format::Sam, &CoverageParams { output: Some(out.clone()), ..Default::default() }).unwrap();
        let s = std::fs::read_to_string(&out).unwrap();
        let chr1 = s.lines().find(|l| l.starts_with("chr1")).unwrap();
        let f: Vec<&str> = chr1.split('\t').collect();
        assert_eq!(f[3], "2", "numreads"); // 2 reads on chr1
        assert_eq!(f[4], "20", "covbases"); // 20 covered bases
    }

    #[test]
    fn malformed_region_errors() {
        let d = workdir();
        let sam = d.join("reg.sam");
        write_sam(&sam);
        // a non-numeric coordinate must error, not silently scan the whole ref
        assert!(mpileup(
            &sam,
            Format::Sam,
            &MpileupParams { region: Some("chr1:foo-bar".into()), output: Some(d.join("x")), ..Default::default() },
        )
        .is_err());
        // start > end must error
        assert!(coverage(
            &sam,
            Format::Sam,
            &CoverageParams { region: Some("chr1:200-100".into()), output: Some(d.join("y")) },
        )
        .is_err());
        // a well-formed region still succeeds
        assert!(coverage(
            &sam,
            Format::Sam,
            &CoverageParams { region: Some("chr1:1-1000".into()), output: Some(d.join("z")) },
        )
        .is_ok());
    }
}
