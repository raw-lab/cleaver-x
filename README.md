```
   ___ _
  / __\ | ___  __ ___   _____ _ __
 / /  | |/ _ \/ _` \ \ / / _ \ '__|
/ /___| |  __/ (_| |\ V /  __/ |
\____/|_|\___|\__,_| \_/ \___|_|
```

# Cleaver ✂️

[![Rust](https://img.shields.io/badge/rust-1.75%2B-orange.svg)](https://www.rust-lang.org)
[![Build](https://img.shields.io/badge/build-passing-brightgreen.svg)](#-build)
[![Tests](https://img.shields.io/badge/tests-passing-brightgreen.svg)](#-correctness)
[![unsafe](https://img.shields.io/badge/unsafe-forbidden%20(core)-success.svg)](#-design)
[![Formats](https://img.shields.io/badge/formats-FASTA%20%C2%B7%20FASTQ%20%C2%B7%20SAM%20%C2%B7%20BAM-informational.svg)](#-formats)
[![License](https://img.shields.io/badge/license-CC--BY--NC--4.0-blue.svg)](#-license)
[![Scaling](https://img.shields.io/badge/scaling-rayon%20%2B%20HydraMPP-blueviolet.svg)](#-parallelism-in-node--cross-node)

**A streaming, record-aware splitter and converter for the everyday sequence and
alignment formats.** Cleaver splits FASTA, FASTQ, SAM, and BAM into chunks
**without ever cutting a record in half**; converts between SAM and BAM (with
**mapped/unmapped** partitioning) and FASTQ→FASTA; reports **full genome
statistics** (N50/L50/N90, length min/max/mean/median, GC%, A/C/G/T/N) via the
vendored **[rustyomestats](https://crates.io/crates/rustyomestats)** (`stats`);
and fans work out one-file-per-task — on a work-stealing pool in-node, or across a
whole cluster via the vendored **HydraMPP** engine. The sequence path streams at
constant memory; the alignment path uses
[`noodles`](https://github.com/zaeleus/noodles), the pure-Rust htslib equivalent.

---

## 📦 Formats

| Family | Extensions | Split unit | Engine |
|---|---|---|---|
| FASTA | `.fasta .fa .fna .ffn .faa .frn .mpfa .fas` | size (`--chunk-size`) | streaming, byte-level |
| FASTQ | `.fastq .fq` | size (`--chunk-size`, 4-line aware) | streaming, byte-level |
| SAM | `.sam` | records (`--records`) | text, header replicated |
| BAM | `.bam` | records (`--records`) | `noodles` (BGZF) |

`.gz` inputs are decompressed transparently. Format is auto-detected from the
extension, falling back to a content sniff. Chunks keep the input's extension
(`genome.fna` → `genome.00000.fna`); SAM/BAM chunks are written as `.sam`/`.bam`.

---

## 📥 Install

`cargo build` only writes `target/release/cleaver` — that path is **not** on your
`PATH`, so bare `cleaver` gives `command not found`. The bundled installer builds
the binary and drops a `cleaver` command **directly onto your PATH**:

```bash
./install.sh                 # build + install (adds the dir to PATH if needed)
./install.sh hydra           # build --features hydra, then install
DEST=/usr/local/bin ./install.sh   # system-wide (uses sudo if required)

cleaver doctor               # verify the install is healthy
```

`install.sh` installs to `~/.cargo/bin` (or `~/.local/bin`) and, if that
directory isn't already on your `PATH`, appends it to your shell rc and tells you
to open a new shell. Equivalent manual routes:

```bash
cargo install --path cleaver-cli            # -> ~/.cargo/bin/cleaver (idiomatic)
# or copy the built binary onto an existing PATH dir:
cargo build --release && install -m 0755 target/release/cleaver ~/.local/bin/
```

Add `--features hydra` (cluster backend) and/or `--features gpu` (CUDA stats
kernel) to `./install.sh`, `cargo install`, or `cargo build`.

## 🩺 Self-test (`cleaver doctor`)

`cleaver doctor` verifies the install end-to-end — it writes tiny inputs to a
temp dir and exercises the real engine (FASTA split + lossless reassembly,
FASTQ→FASTA, SAM↔BAM round-trip, mapped/unmapped partition, base-composition
stats, format detection, GPU detection, and — with `--features hydra` — a live
local HydraMPP task), printing a pass/fail line per check and exiting non-zero on
any failure:

```text
$ cleaver doctor
  [ ok ] fasta split            3 chunks, byte-identical reassembly
  [ ok ] sam -> bam -> sam      3 records preserved each way (noodles/BGZF)
  [ ok ] mapped/unmapped split  2 mapped + 1 unmapped, routed by 0x4 flag
  [ ok ] hydra engine           local runtime ran 5 tasks correctly
  ...
  9/9 checks passed.  cleaver is healthy. ✓
```

## 🚀 Usage

```bash
# split — auto-detects format per file; runs files in parallel (-t, 0 = all cores)
cleaver split reads.fastq genome.fna aln.bam -o out/ -c 256M -r 1000000

# every FASTA variant, every FASTQ; gzip transparent
cleaver split contigs.ffn proteins.faa reads.fq.gz -o out/ -c 100M

# convert — direction inferred from extensions
cleaver convert aln.sam aln.bam            # SAM -> BAM (BGZF)
cleaver convert aln.bam aln.sam            # BAM -> SAM
cleaver convert reads.fastq reads.fasta    # FASTQ -> FASTA (drops quality)

# convert with mapped/unmapped selection (alignment inputs)
cleaver convert aln.sam aln.bam --mapped-only      # keep only mapped reads
cleaver convert aln.sam aln.bam --unmapped-only    # keep only unmapped reads
cleaver convert aln.sam aln.bam --partition        # -> aln.mapped.bam + aln.unmapped.bam

# stats — full genome statistics: N50/L50/N90, lengths, GC%, A/C/G/T/N (parallel across files)
cleaver stats genome.fna contigs.ffn reads.fq

# environment / accelerators / backend
cleaver info

# health check (verify the install and every code path)
cleaver doctor

# help & version (these replace --help / --version)
cleaver -h                 # full options + examples   (also: cleaver help)
cleaver help stats         # detailed help for one command
cleaver version            # ASCII banner + build config   (also: -v)
```

### Scaling across a cluster (HydraMPP)

Built `--features hydra`, the same subcommands fan tasks across a
[HydraMPP](https://github.com/raw-lab/HydraMPP) cluster — the engine is vendored
in-tree (`vendor/hydra-mpp-core`), no external dependency:

```bash
cargo build --release --features hydra

# head node (workers connect in); then split/stats fan out across the cluster
cleaver split *.fna -o out/ --head
# each worker joins the head:
cleaver split --client 10.0.0.1
# simulate GPU scheduling on a CPU-only box (devices get pinned per task):
cleaver stats *.fna --sim-gpus 4
```

Splitting bounds map to each family's natural unit: **sequence** files by size
(`-c/--chunk-size`, e.g. `1G`, `256M`, `50Mi`), **alignment** files by record
count (`-r/--records`). Every SAM/BAM chunk carries the full header (and, for
BAM, a finalised BGZF EOF block), so each chunk is an independently valid file.

---

## 📖 CLI reference

```
cleaver [GLOBAL OPTIONS] <COMMAND> [ARGS]
```

`-h` works on **every** command (`cleaver -h`, `cleaver split -h`, …) and is
equivalent to `cleaver help [COMMAND]`. `--help` and `--version` are intentionally
not provided — use `-h` / `help` and `-v` / `version`.

### Global options (accepted by every command)

| Option | Description |
|---|---|
| `-h` | Print help (top-level, or for the current command). Same as `cleaver help [COMMAND]`. |
| `-v` | Print the ASCII banner, version, and build configuration. Same as `cleaver version`. |

### Cluster options — only in the `--features hydra` build (global on every command)

| Option | Value | Default | Description |
|---|---|---|---|
| `--head` | (flag) | off | Run this process as a HydraMPP **head** node; workers connect to it. |
| `--client <ADDR>` | host or host:port | — | Join a head at `ADDR` as a worker/client. |
| `--cpus <N>` | integer | all logical cores | CPUs to advertise to the scheduler. |
| `--sim-gpus <N>` | integer | 0 | Advertise/schedule `N` GPUs; simulates GPUs on a CPU-only box to exercise device pinning. |
| `--port <PORT>` | 1–65535 | HydraMPP default | TCP port for the head/client. |

---

### `cleaver split <INPUTS>... -o <OUTDIR> [OPTIONS]`

Split files into record-aware chunks — FASTA/FASTQ by **size**, SAM/BAM by
**record count** — never cutting a record. Multiple inputs run in parallel.

| Argument / option | Default | Description |
|---|---|---|
| `<INPUTS>...` | *(required, ≥1)* | Input files: FASTA (`.fasta/.fa/.fna/.ffn/.faa/.frn/.mpfa`), FASTQ (`.fastq/.fq`), SAM (`.sam`), BAM (`.bam`). `.gz` is decompressed transparently. |
| `-o, --outdir <OUTDIR>` | *(required)* | Output directory. Chunks are written as `<stem>.NNNNN.<ext>` (zero-padded index), e.g. `genome.00000.fasta`. |
| `-c, --chunk-size <SIZE>` | `1G` | Approx. max **bytes per chunk** for FASTA/FASTQ. A unit suffix is **required**: `B K M G T P E Z Y` (case-insensitive) or IEC `Ki Mi Gi …`, or long forms (`kilo`/`mega`/`giga`/`kibi`/`mebi`/…). **All are 1024-based** (`1K`=1024, `1M`=1048576). Fractions allowed (`0.5G`). Ignored for SAM/BAM. |
| `-r, --records <N>` | `1000000` | Records per chunk for **SAM/BAM**. Ignored for FASTA/FASTQ. |
| `-t, --threads <N>` | `0` (all cores) | Worker threads for the in-node backend. *(Ignored in the `--features hydra` build — use `--cpus`.)* |
| `--format <FMT>` | auto-detect | Force input format: `fasta`, `fastq`, `sam`, or `bam`. |

```bash
cleaver split genome.fasta -o chunks/ -c 500M        # ~500 MB FASTA chunks
cleaver split reads.fastq.gz -o chunks/ -c 200M      # gzip in, FASTQ by size
cleaver split aln.bam -o chunks/ -r 500000           # 500k records / chunk
cleaver split a.fna b.fna c.fna -o chunks/ -t 8      # 3 inputs, 8 worker threads
cleaver split contigs.txt -o chunks/ --format fasta  # force the format
cleaver split genome.fasta -o chunks/ -c 1073741824B # raw bytes (suffix required)
```

---

### `cleaver convert <INPUT> <OUTPUT> [OPTIONS]`

Convert between formats; the **direction is inferred from the file extensions**.
Supports SAM↔BAM (BGZF) and FASTQ→FASTA. The three filter flags apply to
alignment inputs only.

| Argument / option | Default | Description |
|---|---|---|
| `<INPUT>` | *(required)* | Input file (format from extension/content). |
| `<OUTPUT>` | *(required)* | Output file; its **extension selects the target** (`.bam`, `.sam`, `.fasta/.fa`). |
| `--partition` | off | *(alignment only)* Write mapped and unmapped records to **two** files: `<out_stem>.mapped.<ext>` and `<out_stem>.unmapped.<ext>`. Mutually exclusive with `--mapped-only` / `--unmapped-only`. |
| `--mapped-only` | off | *(alignment only)* Keep only mapped records. Mutually exclusive with `--unmapped-only`. |
| `--unmapped-only` | off | *(alignment only)* Keep only unmapped records. |

```bash
cleaver convert aln.sam aln.bam              # SAM -> BAM (BGZF-compressed)
cleaver convert aln.bam aln.sam              # BAM -> SAM
cleaver convert aln.bam kept.bam --mapped-only
cleaver convert aln.bam drop.bam --unmapped-only
cleaver convert aln.sam out.sam --partition  # -> out.mapped.sam + out.unmapped.sam
cleaver convert reads.fastq reads.fasta      # drop quality lines
```

---

### `cleaver stats <INPUTS>... [OPTIONS]`

Full genome statistics — one row per file plus a `TOTAL` row: `seqs`, `total_bp`,
`min` / `max` / `mean` / `median` length, **N50 / L50 / N90**, `GC%`, and
`A / C / G / T / N`. N/L metrics come from rustyomestats' `compute_nl`.

| Argument / option | Default | Description |
|---|---|---|
| `<INPUTS>...` | *(required, ≥1)* | FASTA/FASTQ files; `.gz` transparent. Files run in parallel. |
| `--format <FMT>` | auto-detect | Force `fasta` or `fastq`. (`sam`/`bam` are rejected — stats is for sequence files.) |
| `-t, --threads <N>` | `0` (all cores) | Worker threads for the in-node backend. *(Ignored under `--features hydra`.)* |

```bash
cleaver stats genome.fasta                   # single file
cleaver stats *.fasta *.fastq.gz             # many files + TOTAL row
cleaver stats reads.fastq.gz --format fastq  # force format
cleaver --sim-gpus 4 stats *.fna             # (hydra build) GPU scheduling demo
```

---

### `cleaver version`  ·  `cleaver doctor`  ·  `cleaver info`  ·  `cleaver help`

These take no options.

| Command | What it does |
|---|---|
| `cleaver version` *(or `-v`)* | Print the ASCII banner, version, and build configuration (scaling backend, whether the GPU kernel is compiled, logical cores). |
| `cleaver doctor` | Self-test the install — writes tiny inputs and exercises FASTA split (+ lossless reassembly), FASTQ→FASTA, SAM↔BAM round-trip, mapped/unmapped partition, base composition/GC, **genome stats (N50)**, format detection, GPU detection, and — with `--features hydra` — a live local HydraMPP task. Exits non-zero on any failure. |
| `cleaver info` | Show version, active scaling backend, GPU-kernel status, NVIDIA GPUs detected via `nvidia-smi`, supported formats, and what `convert`/`stats` do. |
| `cleaver help [COMMAND]` | Print help; `cleaver help <command>` shows that command's full options and examples. Equivalent to `-h`. |

```bash
cleaver version          # or: cleaver -v
cleaver doctor
cleaver info
cleaver help split       # detailed help for one command (same as: cleaver split -h)
```

---
---
## 🧩 How it works

```mermaid 
flowchart LR
    A[inputs] --> D{detect format}
    D -->|FASTA/FASTQ| S[streaming engine<br/>O(1) memory, record-aware]
    D -->|SAM| T[text splitter<br/>header replicated]
    D -->|BAM| B[noodles reader/writer<br/>BGZF finalised per chunk]
    A --> P[work-stealing pool<br/>one task per file]
    P --- D
    S --> O[(chunks)]
    T --> O
    B --> O
```

* **Sequence formats** stream one line at a time through 1 MiB buffers, tracking
  bytes in memory (no `fstat` per boundary). A new chunk starts at the next
  record boundary once the current chunk passes the target size, so chunks
  overshoot by at most one record and are never cut mid-record.
* **SAM** is split as text: the leading `@` header block is replicated at the
  top of every chunk.
* **BAM** is split through `noodles`: header written per chunk, records copied,
  and `finish()` called per chunk to emit the BGZF EOF block (a missing EOF is
  the classic way to produce a truncated, unreadable BAM — Cleaver does not).

---

## ⚡ Parallelism (in-node + cross-node)

Work is **one task per file**. The default build runs tasks on a work-stealing
pool ([`rayon`](https://github.com/rayon-rs/rayon)); `-t 0` uses all logical
cores. Built `--features hydra`, the *same* task functions run on a **HydraMPP**
cluster instead — multi-core on one box, or across nodes, with no other code
change. The `cleaver-core` engine runs unchanged on whichever worker receives a
file.

HydraMPP is **vendored in-tree** (`vendor/hydra-mpp-core`, a workspace member),
so there is no external dependency to fetch. Because a scheduler cannot ship a
closure to a remote process, each unit of work is a *named* function over
`serde`-serializable jobs (`cleaver-cli/src/jobs.rs`), registered with HydraMPP
and dispatched with `--head` / `--client ADDR`; the request is clamped to each
node's advertised resources so it schedules rather than fails fast.

> HydraMPP keeps `panic = "unwind"` so one panicking task can't abort a worker;
> Cleaver's release profile matches. The default (rayon) backend remains fully
> tested in-sandbox; the HydraMPP backend is exercised here in local mode.

## 🖥️ GPUs

Splitting and conversion are I/O-bound, so they have **no GPU kernel**. Base
composition / GC (`stats`), by contrast, is a reduction over millions of bytes —
genuinely data-parallel — so it has a real GPU path:

* **Scheduling.** Under HydraMPP, `stats` reserves a GPU per task; the scheduler
  **pins** a concrete device id (race-free, via `current_gpus()`), shown in the
  `device` column. `--sim-gpus N` exercises this on CPU-only hardware (you'll see
  `cpu(gpu:0)` — device pinned, kernel run on CPU).
* **Kernel.** Built `--features gpu`, a [`cudarc`](https://github.com/coreylowman/cudarc)
  kernel streams sequence bytes to the device in 64 MiB batches and builds a
  256-bin byte histogram with one `atomicAdd` per base, then folds the bins into
  A/C/G/T/N/other — *identical arithmetic* to the CPU reference, so results match
  bit-for-bit. `cudarc` needs the CUDA toolkit and rustc ≥ 1.85 at build time, so
  it is **off by default** and built on GPU hosts; the CPU path is what the test
  suite exercises, and `stats` falls back to it automatically when no device is
  present.

`cleaver info` reports the active scaling backend, whether the GPU kernel was
compiled in, and any NVIDIA GPUs detected via `nvidia-smi` (no build dependency).

---

## 🧬 Genome statistics (rustyomestats)

`cleaver stats` reports a full genome-statistics table — one row per file, plus a
`TOTAL` row across files:

```text
  file                       seqs      total_bp      min       max       mean   median        N50   L50       N90    GC%           A           C           G           T         N     device
  genome.fna                  176      4_215_606      842    142_318    23_952    18_440    41_207    34    12_905  61.40   1_010_233   1_097_588   1_101_402   1_006_383     1_998        cpu
```

The N/L assembly metrics (**N25/N50/N75/N90** and the matching **L** indices) are
computed by **rustyomestats**' own `compute_nl` kernel — the exact function the
standalone `rustyomestats` tool uses — so the numbers agree. Lengths and base
composition are gathered in a single streaming pass (one `usize` per sequence is
the only growth in memory); when built `--features gpu`, the A/C/G/T/N tally is
offloaded to the GPU and reused by the length pass.

rustyomestats is **vendored** at `vendor/rustyomestats`. Its pure stat kernels
(`compute_nl`, the 6-frame codon counters, the Castro `n_stat`) need only `std`
and build on **rustc 1.75**, so they are wired into cleaver's default build and
exercised by the test suite (`cargo test -p rustyomestats`). Its heavier
capabilities — **codon density** CSVs, the full **Castro U50/UG50** assembly
metrics from a reference + BED, and **FragGeneScan** ORF-density — depend on
`polars` + `bio` and are gated behind rustyomestats' own `full` feature
(**rustc ≥ 1.85**); build and run them through the vendored crate's CLI:

```bash
cargo run -p rustyomestats --features full -- genome --input contigs.fna --outdir stats/
```

---

## 🏁 Benchmarks

Measured in-sandbox (single core, rustc 1.75, warm cache, best of 3). No numbers
are fabricated; field tools that are not pure-Rust or not reachable here
(`seqkit`, `samtools`) are named but not benchmarked.

**FASTA split** — 180 MB / 35,839 records, `>`-delimited, 32 MB chunks:

| tool | best (s) | peak RSS | record-safe? |
|---|---:|---:|:--:|
| **cleaver** (streaming) | **0.17** | **7 MB** | ✅ 0/6 broken |
| naive Rust (in-memory) | 0.22 | 206 MB | ✅ |
| python original | 0.49 | 13 MB | ✅ |
| GNU `split -b` | 0.06 | 1.8 MB | ❌ 5/6 mid-record |

vs a naive in-memory Rust splitter (same compiler): comparable-to-faster wall
time and **~30× less memory** — the naive RSS scales with file size, Cleaver's is
flat. vs Python: ~3× faster. `split` is faster only because it skips boundary
logic and corrupts 5 of 6 records.

**Alignment** — 300,000 records (47 MB SAM ↔ 9.7 MB BAM), single core:

| operation | time (s) | peak RSS | throughput |
|---|---:|---:|---|
| SAM → BAM convert | 1.38 | 4.4 MB | ~217k rec/s |
| BAM → SAM convert | 0.34 | 3.3 MB | ~880k rec/s |
| BAM split (`-r 50000`) | 1.46 | 3.7 MB | 6 chunks, lossless |
| SAM split (`-r 50000`) | 0.05 | 4.1 MB | 6 chunks, lossless |

Memory stays flat (3–4 MB) regardless of file size on the alignment path too.

**Stats (base composition / GC)** — CPU path, 172 MB FASTA / 175.6 M bases,
single core: **1.52 s at 4.6 MB RSS** (~115 M bases/s), memory flat. The GPU
kernel (`--features gpu`) targets this same reduction on CUDA hardware.

---
---

## 🔧 Build

```bash
# rustc 1.75+ (no rustup needed)
cargo build --release            # -> target/release/cleaver  (rayon backend)
cargo build --release --features hydra   # + HydraMPP cluster backend
cargo build --release --features gpu     # + CUDA stats kernel (needs CUDA toolkit, rustc 1.85+)
cargo test --workspace                   # unit + doctests, all passing
```

MSRV is held at **1.75** by pinning `noodles-sam=0.52`, `noodles-bam=0.55`,
`noodles-bgzf=0.26`, `indexmap=2.2.6`, `rayon=1.10`, `rayon-core=1.12.1`. A fully
static binary builds with the `x86_64-unknown-linux-musl` target. The vendored
**rustyomestats** exposes its pure stat kernels (`compute_nl`, codon counters)
with no extra dependencies on 1.75; its `full` feature (codon-density / U50 /
FragGeneScan via `polars` + `bio`) requires **rustc ≥ 1.85** and is built through
the vendored crate, not cleaver.

### ✅ Correctness

`#![forbid(unsafe_code)]` on the core crate (the CLI's GPU launch is the only
`unsafe`, behind `--features gpu`). Tests cover FASTA/FASTQ losslessness and
boundary behaviour, SAM header replication, SAM↔BAM round-trips, BAM-chunk
validity, FASTQ→FASTA, mapped/unmapped filtering and partitioning, and base
composition / GC. End-to-end (verified here): all FASTA variants, FASTQ, gzip,
and multi-file parallel split reassemble byte-identical; `--partition` routes
every record by its `0x4` flag (6 mapped + 4 unmapped, lossless); and `stats`
runs across files through the HydraMPP backend with devices pinned per task.

---

## 📚 Citation

```bibtex
@software{white_cleaver_2025,
  author  = {White III, Richard Allen},
  title   = {Cleaver: a streaming, record-aware splitter and converter for
             FASTA/FASTQ/SAM/BAM},
  year    = {2025},
  note    = {Pure-Rust; noodles for SAM/BAM; rayon in-node, HydraMPP cross-node},
  license = {CC-BY-NC-4.0}
}
```

## 📄 License

Creative Commons Attribution-NonCommercial 4.0 International (**CC-BY-NC-4.0**).
Free for academic and non-commercial use; contact the author for commercial use.
