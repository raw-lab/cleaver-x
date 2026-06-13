//! Alignment-format support: split SAM/BAM by record count (header replicated
//! into every chunk so each is independently valid), and convert between SAM
//! and BAM with optional mapped/unmapped filtering. Conversion uses the
//! format-agnostic `RecordBuf` bridge so a record read from either format can be
//! written to either format.

use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam::{self as sam, alignment::io::Write as _};

use crate::formats::{self, Format};

/// Mapped/unmapped record selection (SAM FLAG 0x4 = segment unmapped).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MapFilter {
    All,
    Mapped,
    Unmapped,
}

impl MapFilter {
    #[inline]
    fn keep(self, is_unmapped: bool) -> bool {
        match self {
            MapFilter::All => true,
            MapFilter::Mapped => !is_unmapped,
            MapFilter::Unmapped => is_unmapped,
        }
    }
}

/// Stem of `path`, ignoring a trailing `.gz`/`.bgz` and the format extension.
fn stem(path: &Path) -> String {
    let name = path.file_name().and_then(|s| s.to_str()).unwrap_or("chunk");
    let base = name.strip_suffix(".gz").or_else(|| name.strip_suffix(".bgz")).unwrap_or(name);
    base.rsplit_once('.').map(|(s, _)| s).unwrap_or(base).to_string()
}

// ---------------------------------------------------------------------------
// SAM split (text; header replicated)
// ---------------------------------------------------------------------------

/// Split a SAM file into chunks of at most `records_per_chunk` records. The
/// `@`-prefixed header block heads every chunk, so each is a valid SAM. gzip
/// input is transparent.
pub fn split_sam(input: &Path, dest: &Path, records_per_chunk: u64) -> Result<Vec<PathBuf>> {
    if records_per_chunk == 0 {
        bail!("--records must be >= 1");
    }
    fs::create_dir_all(dest)
        .with_context(|| format!("creating output directory '{}'", dest.display()))?;
    let (mut reader, _gz) = crate::open_reader(input)?;
    let name = stem(input);
    let path_for = |i: u64| dest.join(format!("{name}.{i:05}.sam"));

    let mut header: Vec<u8> = Vec::new();
    let mut header_done = false;
    let mut chunks: Vec<PathBuf> = Vec::new();
    let mut idx = 0u64;
    let mut writer: Option<BufWriter<File>> = None;
    let mut in_chunk = 0u64;
    let mut line: Vec<u8> = Vec::with_capacity(256);

    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line).context("reading SAM")?;
        if n == 0 {
            break;
        }
        if !header_done && line.first() == Some(&b'@') {
            header.extend_from_slice(&line);
            continue;
        }
        header_done = true;
        if writer.is_none() || in_chunk >= records_per_chunk {
            if let Some(mut w) = writer.take() {
                w.flush().context("flushing SAM chunk")?;
            }
            let p = path_for(idx);
            let mut w = BufWriter::new(
                File::create(&p).with_context(|| format!("creating chunk '{}'", p.display()))?,
            );
            w.write_all(&header).context("writing SAM header")?;
            chunks.push(p);
            writer = Some(w);
            in_chunk = 0;
            idx += 1;
        }
        if let Some(w) = writer.as_mut() {
            w.write_all(&line).context("writing SAM record")?;
            in_chunk += 1;
        }
    }

    if let Some(mut w) = writer.take() {
        w.flush().context("flushing final SAM chunk")?;
    } else {
        let p = path_for(0);
        let mut w = BufWriter::new(File::create(&p)?);
        w.write_all(&header)?;
        w.flush()?;
        chunks.push(p);
    }
    Ok(chunks)
}

// ---------------------------------------------------------------------------
// BAM split (noodles; BGZF finalised per chunk)
// ---------------------------------------------------------------------------

/// Split a BAM into chunks of at most `records_per_chunk` records. Every chunk
/// carries the header and a finalised BGZF EOF block, so each is a valid BAM.
pub fn split_bam(input: &Path, dest: &Path, records_per_chunk: u64) -> Result<Vec<PathBuf>> {
    if records_per_chunk == 0 {
        bail!("--records must be >= 1");
    }
    fs::create_dir_all(dest)
        .with_context(|| format!("creating output directory '{}'", dest.display()))?;
    let file = File::open(input).with_context(|| format!("opening '{}'", input.display()))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().context("reading BAM header")?;

    let name = stem(input);
    let path_for = |i: u64| dest.join(format!("{name}.{i:05}.bam"));

    let make = |path: &PathBuf| -> Result<bam::io::Writer<bgzf::Writer<File>>> {
        let mut writer = bam::io::Writer::new(
            File::create(path).with_context(|| format!("creating chunk '{}'", path.display()))?,
        );
        writer.write_header(&header).context("writing BAM header")?;
        Ok(writer)
    };

    let mut chunks = Vec::new();
    let mut idx = 0u64;
    let mut path = path_for(idx);
    let mut writer = make(&path)?;
    chunks.push(path.clone());
    let mut in_chunk = 0u64;

    for result in reader.records() {
        let record = result.context("reading BAM record")?;
        if in_chunk >= records_per_chunk {
            writer.finish(&header).context("finalising BAM chunk")?;
            idx += 1;
            path = path_for(idx);
            writer = make(&path)?;
            chunks.push(path.clone());
            in_chunk = 0;
        }
        writer.write_record(&header, &record).context("writing BAM record")?;
        in_chunk += 1;
    }
    writer.finish(&header).context("finalising final BAM chunk")?;
    Ok(chunks)
}

/// Alignment split entry point.
pub fn split_alignment(input: &Path, dest: &Path, records_per_chunk: u64, fmt: Format) -> Result<Vec<PathBuf>> {
    match fmt {
        Format::Sam => split_sam(input, dest, records_per_chunk),
        Format::Bam => split_bam(input, dest, records_per_chunk),
        _ => bail!("split_alignment called on a non-alignment format"),
    }
}

// ---------------------------------------------------------------------------
// Conversion (SAM <-> BAM via RecordBuf; FASTQ -> FASTA)
// ---------------------------------------------------------------------------

pub(crate) enum AlignReader {
    Sam(sam::io::Reader<Box<dyn BufRead>>),
    Bam(bam::io::Reader<bgzf::Reader<File>>),
}

pub(crate) fn open_align_reader(input: &Path, fmt: Format) -> Result<AlignReader> {
    Ok(match fmt {
        Format::Sam => {
            let (br, _gz) = crate::open_reader(input)?;
            AlignReader::Sam(sam::io::Reader::new(br))
        }
        Format::Bam => {
            let file = File::open(input).with_context(|| format!("opening '{}'", input.display()))?;
            AlignReader::Bam(bam::io::Reader::new(file))
        }
        _ => bail!("expected SAM/BAM input"),
    })
}

pub(crate) fn read_align_header(reader: &mut AlignReader) -> Result<sam::Header> {
    Ok(match reader {
        AlignReader::Sam(r) => r.read_header().context("reading SAM header")?,
        AlignReader::Bam(r) => r.read_header().context("reading BAM header")?,
    })
}

fn make_align_writer(output: &Path, fmt: Format) -> Result<Box<dyn sam::alignment::io::Write>> {
    let out = BufWriter::new(
        File::create(output).with_context(|| format!("creating '{}'", output.display()))?,
    );
    Ok(match fmt {
        Format::Sam => Box::new(sam::io::Writer::new(out)),
        Format::Bam => Box::new(bam::io::Writer::new(out)),
        _ => bail!("convert output must be .sam or .bam"),
    })
}

/// Convert between alignment formats (SAM <-> BAM) with mapped/unmapped
/// filtering. Returns the number of records written.
pub fn convert_alignment(input: &Path, output: &Path, filter: MapFilter) -> Result<u64> {
    let in_fmt = formats::detect(input)?;
    let out_fmt = formats::detect_by_ext(output)
        .with_context(|| format!("unrecognised output extension for '{}'", output.display()))?;
    if !in_fmt.is_alignment() || !out_fmt.is_alignment() {
        bail!("convert_alignment expects SAM/BAM endpoints");
    }

    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let mut writer = make_align_writer(output, out_fmt)?;
    writer.write_alignment_header(&header).context("writing header")?;

    let mut n = 0u64;
    macro_rules! drain {
        ($r:expr) => {
            for result in $r.record_bufs(&header) {
                let rec = result.context("reading record")?;
                if filter.keep(rec.flags().is_unmapped()) {
                    writer.write_alignment_record(&header, &rec).context("writing record")?;
                    n += 1;
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => drain!(r),
        AlignReader::Bam(r) => drain!(r),
    }
    writer.finish(&header).context("finalising output")?;
    Ok(n)
}

/// Convert an alignment file, routing mapped records to `mapped_out` and
/// unmapped records to `unmapped_out` in a single pass. Output formats are taken
/// from each path's extension. Returns `(mapped, unmapped)` counts.
pub fn convert_alignment_split(
    input: &Path,
    mapped_out: &Path,
    unmapped_out: &Path,
) -> Result<(u64, u64)> {
    let in_fmt = formats::detect(input)?;
    if !in_fmt.is_alignment() {
        bail!("split-by-mapping expects a SAM/BAM input");
    }
    let mfmt = formats::detect_by_ext(mapped_out)
        .with_context(|| format!("unrecognised extension for '{}'", mapped_out.display()))?;
    let ufmt = formats::detect_by_ext(unmapped_out)
        .with_context(|| format!("unrecognised extension for '{}'", unmapped_out.display()))?;

    let mut reader = open_align_reader(input, in_fmt)?;
    let header = read_align_header(&mut reader)?;
    let mut wm = make_align_writer(mapped_out, mfmt)?;
    let mut wu = make_align_writer(unmapped_out, ufmt)?;
    wm.write_alignment_header(&header).context("writing mapped header")?;
    wu.write_alignment_header(&header).context("writing unmapped header")?;

    let (mut nm, mut nu) = (0u64, 0u64);
    macro_rules! drain {
        ($r:expr) => {
            for result in $r.record_bufs(&header) {
                let rec = result.context("reading record")?;
                if rec.flags().is_unmapped() {
                    wu.write_alignment_record(&header, &rec).context("writing unmapped")?;
                    nu += 1;
                } else {
                    wm.write_alignment_record(&header, &rec).context("writing mapped")?;
                    nm += 1;
                }
            }
        };
    }
    match &mut reader {
        AlignReader::Sam(r) => drain!(r),
        AlignReader::Bam(r) => drain!(r),
    }
    wm.finish(&header).context("finalising mapped")?;
    wu.finish(&header).context("finalising unmapped")?;
    Ok((nm, nu))
}

/// Convert FASTQ to FASTA (drops quality). gzip input transparent. Returns the
/// number of records written.
pub fn fastq_to_fasta(input: &Path, output: &Path) -> Result<u64> {
    let (mut reader, _gz) = crate::open_reader(input)?;
    let mut out = BufWriter::new(
        File::create(output).with_context(|| format!("creating '{}'", output.display()))?,
    );
    let mut n = 0u64;
    let mut lines: [Vec<u8>; 4] = Default::default();
    loop {
        let mut got = 0;
        for slot in lines.iter_mut() {
            slot.clear();
            let r = reader.read_until(b'\n', slot).context("reading FASTQ")?;
            if r == 0 {
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
        out.write_all(b">").context("writing FASTA")?;
        let hdr = &lines[0];
        let body = if hdr.first() == Some(&b'@') { &hdr[1..] } else { &hdr[..] };
        out.write_all(body).context("writing FASTA")?;
        out.write_all(&lines[1]).context("writing FASTA")?;
        n += 1;
    }
    out.flush().context("flushing FASTA")?;
    Ok(n)
}

/// Top-level convert dispatcher (CLI). `filter` only affects SAM/BAM endpoints.
pub fn convert(input: &Path, output: &Path, filter: MapFilter) -> Result<(u64, &'static str)> {
    let in_fmt = formats::detect(input)?;
    let out_fmt = formats::detect_by_ext(output)
        .with_context(|| format!("unrecognised output extension for '{}'", output.display()))?;
    match (in_fmt, out_fmt) {
        (a, b) if a.is_alignment() && b.is_alignment() => {
            Ok((convert_alignment(input, output, filter)?, "alignment"))
        }
        (Format::Fastq, Format::Fasta) => Ok((fastq_to_fasta(input, output)?, "fastq->fasta")),
        (a, b) => bail!("unsupported conversion {a:?} -> {b:?} (supported: SAM<->BAM, FASTQ->FASTA)"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn tmp(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("cleaver-align-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    // 2 header lines, 4 records: r1,r2,r4 mapped; r3 unmapped (FLAG 4).
    const SAM: &str = "\
@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:ref\tLN:100
r1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\tIIII
r2\t0\tref\t5\t60\t4M\t*\t0\t0\tTTTT\tIIII
r3\t4\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII
r4\t0\tref\t13\t60\t4M\t*\t0\t0\tCCCC\tIIII
";

    #[test]
    fn sam_split_replicates_header_and_preserves_records() {
        let d = tmp("sam");
        let inp = d.join("aln.sam");
        File::create(&inp).unwrap().write_all(SAM.as_bytes()).unwrap();
        let chunks = split_sam(&inp, &d.join("out"), 2).unwrap();
        assert!(chunks.len() >= 2);
        let mut total = 0;
        for c in &chunks {
            let t = fs::read_to_string(c).unwrap();
            assert!(t.starts_with("@HD\t"));
            total += t.lines().filter(|l| !l.starts_with('@')).count();
        }
        assert_eq!(total, 4);
    }

    #[test]
    fn sam_to_bam_to_sam_roundtrip() {
        let d = tmp("convert");
        let s = d.join("in.sam");
        File::create(&s).unwrap().write_all(SAM.as_bytes()).unwrap();
        let bam = d.join("mid.bam");
        assert_eq!(convert_alignment(&s, &bam, MapFilter::All).unwrap(), 4);
        let out = d.join("out.sam");
        assert_eq!(convert_alignment(&bam, &out, MapFilter::All).unwrap(), 4);
        assert_eq!(fs::read_to_string(&out).unwrap().lines().filter(|l| !l.starts_with('@')).count(), 4);
    }

    #[test]
    fn mapped_unmapped_filter() {
        let d = tmp("filter");
        let s = d.join("in.sam");
        File::create(&s).unwrap().write_all(SAM.as_bytes()).unwrap();
        let bam = d.join("all.bam");
        convert_alignment(&s, &bam, MapFilter::All).unwrap();
        // 3 mapped, 1 unmapped
        assert_eq!(convert_alignment(&bam, &d.join("m.sam"), MapFilter::Mapped).unwrap(), 3);
        assert_eq!(convert_alignment(&bam, &d.join("u.sam"), MapFilter::Unmapped).unwrap(), 1);
        // single-pass split
        let (nm, nu) = convert_alignment_split(&bam, &d.join("m2.bam"), &d.join("u2.bam")).unwrap();
        assert_eq!((nm, nu), (3, 1));
        // each split output re-reads as valid + correct counts
        assert_eq!(convert_alignment(&d.join("m2.bam"), &d.join("m3.sam"), MapFilter::All).unwrap(), 3);
        assert_eq!(convert_alignment(&d.join("u2.bam"), &d.join("u3.sam"), MapFilter::All).unwrap(), 1);
    }

    #[test]
    fn bam_split_chunks_are_valid_and_lossless() {
        let d = tmp("bamsplit");
        let s = d.join("in.sam");
        File::create(&s).unwrap().write_all(SAM.as_bytes()).unwrap();
        let bam = d.join("all.bam");
        convert_alignment(&s, &bam, MapFilter::All).unwrap();
        let chunks = split_bam(&bam, &d.join("out"), 2).unwrap();
        assert!(chunks.len() >= 2);
        let mut total = 0u64;
        for c in &chunks {
            let back = d.join(format!("{}.sam", c.file_stem().unwrap().to_str().unwrap()));
            total += convert_alignment(c, &back, MapFilter::All).unwrap();
        }
        assert_eq!(total, 4);
    }

    #[test]
    fn fastq_to_fasta_drops_quality() {
        let d = tmp("fq2fa");
        let fq = d.join("reads.fastq");
        File::create(&fq).unwrap().write_all(b"@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n").unwrap();
        let fa = d.join("reads.fasta");
        assert_eq!(fastq_to_fasta(&fq, &fa).unwrap(), 2);
        assert_eq!(fs::read_to_string(&fa).unwrap(), ">r1\nACGT\n>r2\nTTTT\n");
    }
}
