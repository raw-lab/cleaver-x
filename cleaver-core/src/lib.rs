#![forbid(unsafe_code)]
//! Cleaver core — a streaming, record-aware file chunker.
//!
//! Splits FASTA/FASTQ (or any line-oriented text) into chunks of approximately
//! a target size **without ever splitting a record**. The file is streamed one
//! line at a time through large buffered I/O, so peak memory is the buffers
//! plus a single line — independent of input size.
//!
//! Two divergences from the original Python, both deliberate bug fixes:
//!   1. Delimiter matching defaults to *start of line* (`Match::StartsWith`)
//!      rather than substring containment. `> ` only legitimately starts a
//!      FASTA header, so this avoids false boundaries. The old "anywhere"
//!      behaviour is still available via [`Match::Contains`].
//!   2. We track bytes written in memory instead of `flush()`+`fstat()` on
//!      every boundary line. Identical chunk sizing, none of the syscalls.

use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

mod size;
pub use size::parse_size;

pub mod align;
pub mod compute;
pub mod formats;
pub mod genome;
pub mod stats;
pub use formats::Format;

/// I/O buffer size (1 MiB) — large buffers keep read/write syscalls rare.
const BUF: usize = 1 << 20;

/// How a delimiter line is recognised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Match {
    /// Boundary if the line *starts with* the delimiter (correct for FASTA).
    StartsWith,
    /// Boundary if the line *contains* the delimiter anywhere
    /// (faithful to the original Python `delim in line`).
    Contains,
}

/// Splitting strategy.
#[derive(Debug, Clone)]
pub enum Mode {
    /// Break at lines matching `delim` (e.g. `>` for FASTA).
    Delimiter { delim: Vec<u8>, matching: Match },
    /// Break every `n` lines (use `4` for FASTQ).
    Lines { n: u64 },
}

/// A chunking job.
#[derive(Debug, Clone)]
pub struct Config {
    /// Approximate maximum bytes per chunk.
    pub chunk_size: u64,
    pub mode: Mode,
}

/// Result of a chunking run.
#[derive(Debug)]
pub struct Stats {
    /// Paths of the chunks actually created, in order.
    pub chunks: Vec<PathBuf>,
    pub bytes_in: u64,
    pub lines_in: u64,
    pub gzip_input: bool,
}

/// `os.path.splitext` semantics: leading dots are part of the name.
fn splitext(filename: &str) -> (&str, &str) {
    let bytes = filename.as_bytes();
    let mut start = 0;
    while start < bytes.len() && bytes[start] == b'.' {
        start += 1;
    }
    if let Some(pos) = filename.rfind('.') {
        if pos >= start {
            return (&filename[..pos], &filename[pos..]);
        }
    }
    (filename, "")
}

fn strip_gz(name: &str) -> &str {
    for suf in [".gz", ".gzip", ".GZ", ".bgz"] {
        if let Some(stripped) = name.strip_suffix(suf) {
            return stripped;
        }
    }
    name
}

/// Open `path`, transparently decompressing gzip if detected by magic bytes.
pub fn open_reader(path: &Path) -> Result<(Box<dyn BufRead>, bool)> {
    let file =
        File::open(path).with_context(|| format!("opening input '{}'", path.display()))?;
    let mut reader = BufReader::with_capacity(BUF, file);

    let is_gz = {
        let head = reader
            .fill_buf()
            .with_context(|| format!("reading '{}'", path.display()))?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };

    if is_gz {
        #[cfg(feature = "gzip")]
        {
            let dec = flate2::bufread::MultiGzDecoder::new(reader);
            return Ok((Box::new(BufReader::with_capacity(BUF, dec)), true));
        }
        #[cfg(not(feature = "gzip"))]
        {
            bail!(
                "input '{}' is gzip-compressed but cleaver was built without the `gzip` feature",
                path.display()
            );
        }
    }

    Ok((Box::new(reader), false))
}

/// Naive byte-substring search (only the `Contains` fallback uses it).
#[inline]
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    if needle.is_empty() {
        return true;
    }
    if needle.len() > haystack.len() {
        return false;
    }
    haystack.windows(needle.len()).any(|w| w == needle)
}

/// Split `input` into `dest`, respecting record boundaries.
///
/// A new chunk is started at the next boundary line once the current chunk has
/// reached `chunk_size` bytes, so chunk sizes overshoot by at most one record.
pub fn chunk_file(input: &Path, dest: &Path, cfg: &Config) -> Result<Stats> {
    if let Mode::Lines { n } = &cfg.mode {
        if *n == 0 {
            bail!("--lines must be >= 1");
        }
    }

    fs::create_dir_all(dest)
        .with_context(|| format!("creating output directory '{}'", dest.display()))?;

    let (mut reader, gzip_input) = open_reader(input)?;

    let raw = input
        .file_name()
        .and_then(|s| s.to_str())
        .with_context(|| format!("input '{}' has no valid UTF-8 filename", input.display()))?;
    let base = if gzip_input { strip_gz(raw) } else { raw };
    let (name, ext) = splitext(base);

    let path_for = |i: u64| dest.join(format!("{name}.{i:05}{ext}"));

    let mut chunks: Vec<PathBuf> = Vec::new();
    let mut idx: u64 = 0;
    let first = path_for(idx);
    let mut writer = BufWriter::with_capacity(
        BUF,
        File::create(&first).with_context(|| format!("creating chunk '{}'", first.display()))?,
    );
    chunks.push(first);

    let mut written: u64 = 0; // bytes in the current chunk so far
    let mut bytes_in: u64 = 0;
    let mut lines_in: u64 = 0;
    let mut line: Vec<u8> = Vec::with_capacity(512);

    loop {
        line.clear();
        let n = reader
            .read_until(b'\n', &mut line)
            .context("reading input")? as u64;
        if n == 0 {
            break; // EOF
        }
        bytes_in += n;

        let boundary = match &cfg.mode {
            Mode::Delimiter { delim, matching } => match matching {
                Match::StartsWith => line.starts_with(delim.as_slice()),
                Match::Contains => contains(&line, delim),
            },
            Mode::Lines { n } => lines_in % *n == 0,
        };

        if boundary && written >= cfg.chunk_size {
            writer.flush().context("flushing chunk")?;
            idx += 1;
            let p = path_for(idx);
            writer = BufWriter::with_capacity(
                BUF,
                File::create(&p).with_context(|| format!("creating chunk '{}'", p.display()))?,
            );
            chunks.push(p);
            written = 0;
        }

        writer.write_all(&line).context("writing chunk")?;
        written += n;
        lines_in += 1;
    }

    writer.flush().context("flushing final chunk")?;

    Ok(Stats { chunks, bytes_in, lines_in, gzip_input })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpdir(tag: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("cleaver-test-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn write(path: &Path, data: &[u8]) {
        let mut f = File::create(path).unwrap();
        f.write_all(data).unwrap();
    }

    fn cat(paths: &[PathBuf]) -> Vec<u8> {
        let mut out = Vec::new();
        for p in paths {
            out.extend_from_slice(&fs::read(p).unwrap());
        }
        out
    }

    #[test]
    fn splitext_matches_python() {
        assert_eq!(splitext("reads.fasta"), ("reads", ".fasta"));
        assert_eq!(splitext("a.tar.gz"), ("a.tar", ".gz"));
        assert_eq!(splitext("noext"), ("noext", ""));
        assert_eq!(splitext(".bashrc"), (".bashrc", ""));
        assert_eq!(splitext("..foo.txt"), ("..foo", ".txt"));
    }

    #[test]
    fn fasta_never_splits_records_and_is_lossless() {
        let d = tmpdir("fasta");
        let inp = d.join("seqs.fasta");
        // 50 records, each header + 4 wrapped sequence lines.
        let mut data = Vec::new();
        for i in 0..50 {
            data.extend_from_slice(format!(">seq{i} description here\n").as_bytes());
            for _ in 0..4 {
                data.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\n");
            }
        }
        write(&inp, &data);

        let cfg = Config {
            chunk_size: 1024, // tiny -> forces many chunks
            mode: Mode::Delimiter { delim: b">".to_vec(), matching: Match::StartsWith },
        };
        let stats = chunk_file(&inp, &d.join("out"), &cfg).unwrap();

        // lossless: concatenated chunks == original bytes
        assert_eq!(cat(&stats.chunks), data);
        // every chunk starts with a header (no record was split)
        for p in &stats.chunks {
            let bytes = fs::read(p).unwrap();
            assert_eq!(bytes.first(), Some(&b'>'), "chunk {p:?} did not start at a record");
        }
        assert!(stats.chunks.len() > 1);
    }

    #[test]
    fn fastq_by_four_lines_is_lossless() {
        let d = tmpdir("fastq");
        let inp = d.join("reads.fastq");
        let mut data = Vec::new();
        for i in 0..40 {
            // a quality line deliberately starting with '@' (Phred Q31) — would
            // fool a naive '@' delimiter, but line-mode handles it correctly.
            data.extend_from_slice(format!("@read{i}\n").as_bytes());
            data.extend_from_slice(b"ACGTACGTACGTACGT\n");
            data.extend_from_slice(b"+\n");
            data.extend_from_slice(b"@IIIIIIIIIIIIIII\n");
        }
        write(&inp, &data);

        let cfg = Config { chunk_size: 256, mode: Mode::Lines { n: 4 } };
        let stats = chunk_file(&inp, &d.join("out"), &cfg).unwrap();

        assert_eq!(cat(&stats.chunks), data);
        for p in &stats.chunks {
            let bytes = fs::read(p).unwrap();
            assert_eq!(bytes.first(), Some(&b'@'), "chunk {p:?} split a FASTQ record");
            // each chunk holds a whole number of 4-line records
            assert_eq!(bytes.iter().filter(|&&b| b == b'\n').count() % 4, 0);
        }
    }

    #[test]
    fn last_line_without_newline_preserved() {
        let d = tmpdir("nonl");
        let inp = d.join("x.fasta");
        let data = b">a\nACGT\n>b\nTTTT"; // no trailing newline
        write(&inp, data);
        let cfg = Config {
            chunk_size: 1,
            mode: Mode::Delimiter { delim: b">".to_vec(), matching: Match::StartsWith },
        };
        let stats = chunk_file(&inp, &d.join("out"), &cfg).unwrap();
        assert_eq!(cat(&stats.chunks), data);
    }
}
