//! `cleaver doctor` — verify the build is installed and every code path works.
//!
//! This is a real self-test, not a stub: it writes tiny inputs to a temp dir and
//! exercises the actual engine (FASTA splitting, FASTQ→FASTA, SAM↔BAM round-trip,
//! mapped/unmapped partition, base composition, genome stats/N50, and
//! featureCounts/htseq counting), checks GPU detection, and — when built
//! `--features hydra` — runs a live local HydraMPP task. It prints a pass/fail
//! line per check and exits non-zero if anything fails.

use std::fs;
use std::io::Write;

use anyhow::{anyhow, Result};

use cleaver_core::{align, annotation, chunk_file, compute, count, formats, Config, Format, Match, Mode};
use cleaver_core::genome;

use crate::gpu;

fn check(name: &str, f: impl FnOnce() -> Result<String>) -> bool {
    match f() {
        Ok(detail) => {
            println!("  [ ok ] {name:<22} {detail}");
            true
        }
        Err(e) => {
            println!("  [FAIL] {name:<22} {e:#}");
            false
        }
    }
}

/// Read the listed chunk files in order and concatenate their bytes.
fn concat(paths: &[std::path::PathBuf]) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    for p in paths {
        out.extend_from_slice(&fs::read(p)?);
    }
    Ok(out)
}

pub fn run() -> Result<()> {
    print!("{}", crate::BANNER);
    println!("cleaver doctor — self-test\n");

    // ---- build / environment ------------------------------------------------
    let cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    let backend = if cfg!(feature = "hydra") { "HydraMPP (multi-core/cross-node)" } else { "rayon (in-node)" };
    let gpuk = if gpu::compute_enabled() { "compiled" } else { "not compiled" };
    println!("  build:  cleaver v{}", env!("CARGO_PKG_VERSION"));
    println!("  exe:    {}", std::env::current_exe().map(|p| p.display().to_string()).unwrap_or_else(|_| "?".into()));
    println!("  config: backend = {backend} | gpu kernel = {gpuk} | {cores} logical core(s)\n");

    let tmp = std::env::temp_dir().join(format!("cleaver-doctor-{}", std::process::id()));
    fs::create_dir_all(&tmp)?;

    let mut pass = 0usize;
    let mut total = 0usize;
    let mut run = |name: &str, f: &mut dyn FnMut() -> Result<String>| {
        total += 1;
        if check(name, f) {
            pass += 1;
        }
    };

    // ---- temp workspace -----------------------------------------------------
    run("temp workspace", &mut || {
        let probe = tmp.join("probe.txt");
        fs::write(&probe, b"ok")?;
        fs::read(&probe)?;
        Ok(format!("writable at {}", tmp.display()))
    });

    // ---- FASTA split, lossless ---------------------------------------------
    let fasta = tmp.join("d.fasta");
    fs::File::create(&fasta)?.write_all(b">a\nACGTACGT\nACGT\n>b\nGGGGCCCC\n>c\nTTTTAAAA\n")?;
    run("fasta split", &mut || {
        let outdir = tmp.join("fa_chunks");
        fs::create_dir_all(&outdir)?;
        let cfg = Config {
            chunk_size: 8, // tiny: force multiple record-aware chunks
            mode: Mode::Delimiter { delim: b">".to_vec(), matching: Match::StartsWith },
        };
        let st = chunk_file(&fasta, &outdir, &cfg)?;
        if st.chunks.len() < 2 {
            return Err(anyhow!("expected ≥2 chunks, got {}", st.chunks.len()));
        }
        if concat(&st.chunks)? != fs::read(&fasta)? {
            return Err(anyhow!("reassembled chunks differ from input"));
        }
        Ok(format!("{} chunks, byte-identical reassembly", st.chunks.len()))
    });

    // ---- FASTQ -> FASTA -----------------------------------------------------
    let fastq = tmp.join("d.fastq");
    fs::File::create(&fastq)?.write_all(b"@r1\nACGT\n+\nIIII\n@r2\nGGCC\n+\nJJJJ\n")?;
    run("fastq -> fasta", &mut || {
        let out = tmp.join("d.from_fastq.fasta");
        let n = align::fastq_to_fasta(&fastq, &out)?;
        let got = fs::read_to_string(&out)?.lines().filter(|l| l.starts_with('>')).count();
        if n != 2 || got != 2 {
            return Err(anyhow!("expected 2 records, wrote {n}, found {got} headers"));
        }
        Ok("2 records, quality dropped".into())
    });

    // ---- SAM <-> BAM round-trip + mapped/unmapped --------------------------
    let sam = tmp.join("d.sam");
    fs::File::create(&sam)?.write_all(
        b"@HD\tVN:1.6\tSO:unsorted\n\
          @SQ\tSN:chr1\tLN:1000\n\
          r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
          r2\t16\tchr1\t20\t60\t4M\t*\t0\t0\tTTGG\tIIII\n\
          u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n",
    )?;
    run("sam -> bam -> sam", &mut || {
        let bam = tmp.join("d.bam");
        let back = tmp.join("d.back.sam");
        let to_bam = align::convert_alignment(&sam, &bam, align::MapFilter::All)?;
        let to_sam = align::convert_alignment(&bam, &back, align::MapFilter::All)?;
        if to_bam != 3 || to_sam != 3 {
            return Err(anyhow!("record count drift: SAM→BAM={to_bam}, BAM→SAM={to_sam} (expected 3)"));
        }
        Ok("3 records preserved each way (noodles/BGZF)".into())
    });
    run("mapped/unmapped split", &mut || {
        let mbam = tmp.join("d.mapped.bam");
        let ubam = tmp.join("d.unmapped.bam");
        let (nm, nu) = align::convert_alignment_split(&sam, &mbam, &ubam)?;
        if (nm, nu) != (2, 1) {
            return Err(anyhow!("expected (2 mapped, 1 unmapped), got ({nm}, {nu})"));
        }
        Ok("2 mapped + 1 unmapped, routed by 0x4 flag".into())
    });

    // ---- stats / base composition ------------------------------------------
    run("stats (base comp/GC)", &mut || {
        let c = compute::count_file(&fasta, Format::Fasta)?;
        if c.total() != 28 {
            return Err(anyhow!("expected 28 bases, counted {}", c.total()));
        }
        if (c.gc() - 0.5).abs() > 1e-9 {
            return Err(anyhow!("expected GC=0.50, got {:.3}", c.gc()));
        }
        let dev = if gpu::compute_enabled() { "gpu kernel available" } else { "cpu path" };
        Ok(format!("28 bases, GC=50.0% ({dev})"))
    });

    // ---- genome stats (N50/L50) --------------------------------------------
    run("genome stats (N50)", &mut || {
        // d.fasta sequences are 12, 8, 8 bp -> total 28, N50=8, L50=2.
        let s = genome::genome_stats_file(&fasta, Format::Fasta, None)?;
        if (s.n_seqs, s.total_bp) != (3, 28) {
            return Err(anyhow!("expected 3 seqs / 28 bp, got {}/{}", s.n_seqs, s.total_bp));
        }
        if (s.n50(), s.l50()) != (8, 2) {
            return Err(anyhow!("expected N50=8 L50=2, got N50={} L50={}", s.n50(), s.l50()));
        }
        Ok(format!("3 seqs, N50={} L50={} GC={:.1}%", s.n50(), s.l50(), s.gc_percent))
    });

    // ---- count (featureCounts/htseq) ---------------------------------------
    run("count (featureCounts)", &mut || {
        let gtf = tmp.join("d.gtf");
        fs::write(
            &gtf,
            b"chr1\ts\texon\t100\t300\t.\t+\t.\tgene_id \"gX\";\n\
              chr1\ts\texon\t200\t400\t.\t+\t.\tgene_id \"gY\";\n",
        )?;
        let csam = tmp.join("d.count.sam");
        fs::write(
            &csam,
            b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\n\
              rA\t0\tchr1\t120\t60\t50M\t*\t0\t0\t*\t*\n\
              rB\t0\tchr1\t250\t60\t20M\t*\t0\t0\t*\t*\n\
              rC\t0\tchr1\t900\t60\t20M\t*\t0\t0\t*\t*\n",
        )?;
        let ann = annotation::Annotation::from_path(&gtf, "exon", "gene_id")?;
        let fc = count::count_file(&csam, Format::Sam, &ann, &count::CountParams::default())?;
        if (fc.assigned, fc.ambiguous, fc.no_feature) != (1, 1, 1) {
            return Err(anyhow!(
                "expected 1 assigned / 1 ambiguous / 1 no-feature, got {}/{}/{}",
                fc.assigned,
                fc.ambiguous,
                fc.no_feature
            ));
        }
        if ann.n_genes() != 2 || fc.counts[0] != 1 {
            return Err(anyhow!("expected gX=1 over 2 meta-features"));
        }
        Ok("3 reads -> 1 assigned, 1 ambiguous, 1 no-feature".into())
    });

    // ---- format detection ---------------------------------------------------
    run("format detection", &mut || {
        let f1 = formats::detect(&fasta)?;
        let f2 = formats::detect(&fastq)?;
        let f3 = formats::detect(&sam)?;
        if (f1, f2, f3) != (Format::Fasta, Format::Fastq, Format::Sam) {
            return Err(anyhow!("misdetected: {f1:?}/{f2:?}/{f3:?}"));
        }
        Ok("FASTA/FASTQ/SAM detected".into())
    });

    // ---- GPU detection (informational; never fails the run) -----------------
    run("gpu detection", &mut || {
        let gpus = gpu::detect();
        if gpus.is_empty() {
            Ok("no NVIDIA GPU (CPU fallback active)".into())
        } else {
            Ok(format!("{} GPU(s): {}", gpus.len(), gpus.iter().map(|g| g.name.clone()).collect::<Vec<_>>().join(", ")))
        }
    });

    // ---- HydraMPP live local task (only when compiled in) -------------------
    #[cfg(feature = "hydra")]
    run("hydra engine", &mut || hydra_selftest());

    // ---- cleanup + summary --------------------------------------------------
    let _ = fs::remove_dir_all(&tmp);
    println!("\n{pass}/{total} checks passed.");
    if pass == total {
        println!("cleaver is healthy. ✓");
        Ok(())
    } else {
        Err(anyhow!("{} check(s) failed", total - pass))
    }
}

/// Register a trivial task, run several copies through a local HydraMPP runtime,
/// and verify the results — proving the distributed engine works end-to-end.
#[cfg(feature = "hydra")]
fn hydra_selftest() -> Result<String> {
    use hydra_mpp_core::prelude::*;

    fn square(x: u64) -> u64 {
        x * x
    }

    let h = Hydra::new();
    h.register("doctor_square", square);
    h.init(Config::local().quiet(true)).map_err(|e| anyhow!("init: {e}"))?;

    let inputs: Vec<u64> = (1..=5).collect();
    let mut ids = Vec::new();
    for x in &inputs {
        ids.push(h.task("doctor_square").cpus(1).submit(x).map_err(|e| anyhow!("submit: {e}"))?);
    }
    let _ = h.wait(&ids, None, ids.len());
    for (x, id) in inputs.iter().zip(ids) {
        let r = h.get(id).map_err(|e| anyhow!("get: {e}"))?;
        if let Some(e) = r.error() {
            h.shutdown();
            return Err(anyhow!("task errored: {e}"));
        }
        let v: u64 = r.value().map_err(|e| anyhow!("decode: {e}"))?;
        if v != x * x {
            h.shutdown();
            return Err(anyhow!("wrong result: {x}^2 = {v}"));
        }
    }
    h.shutdown();
    Ok("local runtime ran 5 tasks correctly".into())
}
