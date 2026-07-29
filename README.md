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

**A streaming, record-aware toolkit for the everyday sequence and alignment
formats.** Cleaver splits FASTA, FASTQ, SAM, and BAM into chunks **without ever
cutting a record in half**; converts between SAM and BAM (with **mapped/unmapped**
partitioning) and FASTQ→FASTA; reports **genome/assembly statistics** (N50/L50/N90,
lengths, GC%); **counts reads per feature** from one or more BAM/SAM files against
a GTF/GFF annotation (**featureCounts / htseq / VERSE**, including multi-feature
hierarchical assignment); **preprocesses FASTQ** with a pure-Rust **fastp**
re-implementation (adapter/quality/polyX trimming, filtering, PE overlap
correction, JSON report); **demultiplexes ONT reads by barcode** with a pure-Rust
**barbell**-style fit-aligner; and fans work out one-file-per-task — on a
work-stealing pool in-node, or across a whole cluster via the vendored
**HydraMPP** engine. The sequence path streams at constant memory; the alignment
path uses [`noodles`](https://github.com/zaeleus/noodles), the pure-Rust htslib
equivalent.

Everything is **100% Rust** — `fastp`, the VERSE counting modes, and the
`barbell` demultiplexer are native re-implementations, not shell-outs, so the
whole toolkit builds from `cargo` with no C/C++ or Python dependencies.

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

> **`cargo build` does *not* put `cleaver` on your `PATH`.** It writes the binary
> to `target/release/cleaver`; typing bare `cleaver` then gives
> `cleaver: command not found`. Pick one of these:

```bash
# A) one-shot installer: builds release + puts `cleaver` on your PATH (recommended)
./install.sh                  # default backend; ./install.sh hydra  or  hydra,gpu
#   DEST=/usr/local/bin ./install.sh   to choose the install dir

# B) cargo install
cargo install --path cleaver-cli            # -> ~/.cargo/bin/cleaver
export PATH="$HOME/.cargo/bin:$PATH"        # if not already; add to ~/.bashrc / ~/.zshrc

# C) or build and run by full path
cargo build --release
./target/release/cleaver doctor

# then confirm the install is healthy:
cleaver doctor
cleaver --version      # also -v
```

Add `--features hydra` (cluster backend) and/or `--features gpu` (CUDA stats
kernel) to either `cargo install` or `cargo build`.

## 🩺 Self-test (`cleaver doctor`)

`cleaver doctor` verifies the install end-to-end — it writes tiny inputs to a
temp dir and exercises the real engine (FASTA split + lossless reassembly,
FASTQ→FASTA, SAM↔BAM round-trip, mapped/unmapped partition, base-composition
stats, genome N50, **feature counting**, **fastp QC**, **barcode demux**, format
detection, GPU detection, and — with `--features hydra` — a live local HydraMPP
task), printing a pass/fail line per check and exiting non-zero on any failure:

```text
$ cleaver doctor
  [ ok ] fasta split            3 chunks, byte-identical reassembly
  [ ok ] sam -> bam -> sam      3 records preserved each way (noodles/BGZF)
  [ ok ] mapped/unmapped split  2 mapped + 1 unmapped, routed by 0x4 flag
  [ ok ] count (featureCounts)  3 reads -> 1 assigned, 1 ambiguous, 1 no-feature
  [ ok ] fastp (FASTQ QC)       3 reads -> 1 pass, 1 too-short, 1 low-quality
  [ ok ] demux (barcodes)       2 reads -> 1 barcoded (trimmed), 1 unclassified
  [ ok ] hydra engine           local runtime ran 5 tasks correctly
  ...
  12/12 checks passed.  cleaver is healthy. ✓
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

# stats — genome/assembly statistics: N50/L50/N90, lengths, GC% (GPU base-comp if built --features gpu)
cleaver stats genome.fna contigs.ffn reads.fq

# count — featureCounts/htseq/VERSE reads-per-feature from BAM/SAM vs a GTF/GFF
cleaver count sampleA.bam sampleB.bam sampleC.bam -a genes.gtf -o counts.tsv
#   -> counts.tsv (gene x sample matrix) + counts.tsv.summary (assignment categories)
# VERSE modes: -z 0..5; multi-feature hierarchical/independent via -t 'exon;intron'
cleaver count sample.bam -a genes.gtf -o counts.tsv -t 'exon;intron' --scheme hierarchical

# fastp — pure-Rust FASTQ preprocessing (adapter/quality/polyX trim + filtering, JSON report)
cleaver fastp -i reads.fq.gz -o clean.fq.gz -a AGATCGGAAGAGC -g -q 20 -l 30 -j qc.json
cleaver fastp -i R1.fq.gz -I R2.fq.gz -o R1.clean.fq.gz -O R2.clean.fq.gz -2 -c   # paired-end

# demux — pure-Rust ONT barcode demultiplexing (fit-align both ends, trim, split)
cleaver demux -i reads.fq.gz -o demux/ -q barcodes.fasta --min-score 0.8

# environment / accelerators / backend
cleaver info

# health check (verify the install and every code path)
cleaver doctor

# version
cleaver --version          # or -v
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

`-h` shows help for any command (every command prints the banner); `-v` /
`--version` prints the version. Built `--features hydra`, the cluster flags
`--head`, `--client <ADDR>`, `--cpus <N>`, `--sim-gpus <N>`, `--port <PORT>` are
available globally.

### `cleaver split <INPUTS>…`
Record-aware chunking, one task per file.

| Flag | Default | Meaning |
|---|---|---|
| `-o, --outdir <DIR>` | *required* | Output dir; chunks named `<stem>.NNNNN.<ext>`. |
| `-c, --chunk-size <SIZE>` | `1G` | Max chunk size for FASTA/FASTQ (`1G`, `256M`, `50Mi`, …). |
| `-r, --records <N>` | `1000000` | Records per chunk for SAM/BAM. |
| `-t, --threads <N>` | `0` | In-node worker threads (0 = all cores). |
| `--format <fasta\|fastq\|sam\|bam>` | auto | Force a format. |

### `cleaver convert <INPUT> <OUTPUT>`
SAM↔BAM (BGZF) and FASTQ→FASTA; the target is taken from `<OUTPUT>`'s extension.

| Flag | Meaning |
|---|---|
| `--partition` | Alignment: write `<out>.mapped.<ext>` + `<out>.unmapped.<ext>`. |
| `--mapped-only` | Keep mapped records only. |
| `--unmapped-only` | Keep unmapped records only. |

### `cleaver stats <INPUTS>…`
Genome/assembly statistics for FASTA/FASTQ (rejects SAM/BAM). Per file and a
`TOTAL`: sequence count, total bp, min/max/mean/median length, **N50 / L50 /
N90**, GC%, and the compute device. Base composition runs on the GPU with
`--features gpu`; the length/N50 pass is CPU. The N/L assembly metrics are
computed by the vendored [`rustyomestats`](https://github.com/raw-lab/rustyomestats)
crate (`stats::compute_nl`); it is vendored in-tree at `vendor/rustyomestats` and
built dependency-free (its heavier `full` feature — bio/polars/plotters and the
CLI — is not enabled, so cleaver still compiles on rustc 1.75).

| Flag | Default | Meaning |
|---|---|---|
| `-t, --threads <N>` | `0` | In-node worker threads (0 = all cores). |
| `--format <…>` | auto | Force a format. |

### `cleaver count <INPUTS>… -a <GTF/GFF> -o <MATRIX>`
featureCounts / htseq / **VERSE** reads-per-feature from **one or more** BAM/SAM
files against a GTF/GFF annotation. Flag names mirror featureCounts and VERSE.

| Flag | Default | Meaning |
|---|---|---|
| `-a, --annotation <FILE>` | *required* | GTF or GFF3 annotation. |
| `-o, --output <FILE>` | *required* | Count matrix (TSV); `<output>.summary` written alongside. |
| `-t, --feature-type <TYPE[;TYPE…]>` | `exon` | Column-3 feature type(s). Multiple `;`-separated types enable VERSE multi-feature counting. |
| `--scheme <independent\|hierarchical>` | `independent` | With multiple `-t` types: count each type separately, or assign each read to the first matching type. |
| `-g, --group-by <ATTR>` | `gene_id` | Attribute grouping features into meta-features. |
| `-s, --stranded <0\|1\|2>` | `0` | 0 unstranded · 1 stranded · 2 reverse-stranded. |
| `-z, --assign-mode <0–5>` | *(via `--mode`)* | VERSE assignment mode (table below). Overrides `--mode`. |
| `--mode <union\|strict\|nonempty>` | `union` | htseq overlap resolution (alias for `-z 1/2/3`). |
| `-Q, --min-mapq <N>` | `0` | Drop alignments below this MAPQ. |
| `-M, --count-multimappers` | off | Count `NH>1` reads (else discarded). |
| `-O, --allow-multi-overlap` | off | Count reads hitting several meta-features for all of them. |
| `-B, --require-both-ends` | off | PE: only count fragments with both mates mapped. |
| `-C, --exclude-chimeric` | off | PE: drop supplementary reads / mates on another reference. |
| `-P, --check-pe-dist` | off | PE: only count fragments whose `|TLEN|` ∈ `[-d, -D]`. |
| `-d, --min-frag-len <N>` | `50` | PE: minimum fragment length (with `-P`). |
| `-D, --max-frag-len <N>` | `600` | PE: maximum fragment length (with `-P`). |
| `-T, --threads <N>` | `0` | In-node worker threads, one file per task (0 = all cores). |

**VERSE assignment modes (`-z`).** These select how a read's overlap with the
annotation is resolved:

| `-z` | Name | Behaviour |
|---|---|---|
| `0` | featureCounts | any meta-feature overlapping the read (default semantics). |
| `1` | htseq union | union of overlapping features; >1 gene ⇒ ambiguous. |
| `2` | intersection-strict | every base of the read must be covered. |
| `3` | intersection-nonempty | strict over only the covered bases. |
| `4` | union-strict | union must map to one gene **and** cover the whole read. |
| `5` | cover-length | assign to the gene with the **largest overlap**; ambiguous only on an exact tie. |

**Outputs.** For a single feature type, `<output>` is a meta-feature × sample
matrix (a `gene_id` column then one column per input, named by basename) and
`<output>.summary` is the category table (`Assigned`, `Unassigned_NoFeatures`,
`Unassigned_Ambiguity`, `Unassigned_MultiMapping`, `Unassigned_MappingQuality`,
`Unassigned_FragmentLength`, `Unassigned_Unmapped`). For **multiple** feature
types the matrix path becomes a prefix: each type is written to
`<stem>.<type>.<ext>` (e.g. `counts.exon.tsv`, `counts.intron.tsv`) with a single
combined `<output>.summary`. Under `--scheme hierarchical` the summary reports a
per-type `Assigned_<type>` row in priority order.

```bash
# exons by gene across three BAMs, unstranded, drop MAPQ < 10
cleaver count A.bam B.bam C.bam -a genes.gtf -o counts.tsv -Q 10

# reverse-stranded library; count CDS grouped by gene; strict overlaps
cleaver count *.bam -a anno.gff3 -o cds.tsv -s 2 -t CDS --mode strict

# VERSE intersection-nonempty, stranded (== ./verse -z 3 -s 1)
cleaver count sample.bam -a genes.gtf -o nonempty.tsv -z 3 -s 1

# VERSE cover-length: never ambiguous unless tied
cleaver count sample.bam -a genes.gtf -o cover.tsv -z 5

# VERSE multi-feature, INDEPENDENT -> counts.exon.tsv + counts.intron.tsv
cleaver count sample.bam -a genes.gtf -o counts.tsv -t 'exon;intron'

# VERSE multi-feature, HIERARCHICAL (exon wins, else intron, else intergenic)
cleaver count sample.bam -a genes.gtf -o counts.tsv -t 'exon;intron;xine' \
    --scheme hierarchical

# paired-end: both ends mapped, drop chimeras, fragments 50–600 bp
cleaver count pe.bam -a genes.gtf -o pe.tsv -B -C -P -d 50 -D 600
```

**Counting model.** Each primary, mapped alignment passing the MAPQ / multimapper
/ PE filters is reduced to its reference blocks (CIGAR `M/=/X/D` extend a block,
`N` splits at introns) and tested against the index. Exactly one meta-feature →
*assigned*; none → *no-feature*; more than one → *ambiguous* (unless `-O`). For
multiple feature types, **independent** runs each type as its own pass while
**hierarchical** tries the types in `-t` order and assigns the read to the first
that yields a unique meta-feature (VERSE's hierarchical scheme). The N/L and
overlap algorithms follow the standard featureCounts / htseq / VERSE definitions.

**Caveats.** Reads are counted individually (single-end semantics); the PE flags
(`-B/-C/-P/-d/-D`) filter at the read level using SAM flags and `TLEN` rather
than reconstructing fragments, so a properly-paired library is counted per mate,
not de-duplicated per fragment. The annotation is parsed once per feature type
per input; under `--features hydra` it must be readable on each worker node.

### `cleaver fastp -i <IN> -o <OUT> [-I <IN2> -O <OUT2>]`
A pure-Rust re-implementation of the [fastp](https://github.com/OpenGene/fastp)
preprocessing pipeline: per-read trimming, filtering, paired-end overlap analysis,
and a fastp-compatible JSON report. Short flags follow fastp.

**Adapter trimming**

| Flag | Default | Meaning |
|---|---|---|
| `-A, --disable-adapter-trimming` | off | Turn adapter trimming off. |
| `-a, --adapter-sequence <SEQ>` | auto | Adapter for read1, trimmed from the 3′ end. |
| `--adapter-sequence-r2 <SEQ>` | — | Adapter for read2. |
| `--adapter-fasta <FILE>` | — | Trim every sequence in this FASTA from both reads. |
| `-2, --detect-adapter-for-pe` | off | PE: detect adapters from the read overlap. |

**Fixed / global trimming**

| Flag | Default | Meaning |
|---|---|---|
| `-f, --trim-front1 <N>` / `-t, --trim-tail1 <N>` | `0` | Trim N bases from read1 5′ / 3′. |
| `-b, --max-len1 <N>` | `0` | Keep at most N bases of read1 (0 = no cap). |
| `-F / -T / -B` | `0` | The read2 equivalents. |

**polyG / polyX**

| Flag | Default | Meaning |
|---|---|---|
| `-g, --trim-poly-g` | off | Trim polyG tails (NovaSeq/NextSeq dark cycles). |
| `--poly-g-min-len <N>` | `10` | Minimum polyG run length. |
| `-x, --trim-poly-x` | off | Trim polyX tails. |
| `--poly-x-min-len <N>` | `10` | Minimum polyX run length. |

**Sliding-window quality cutting**

| Flag | Default | Meaning |
|---|---|---|
| `-5, --cut-front` | off | Cut low-quality bases from the 5′ end. |
| `-3, --cut-tail` | off | Cut low-quality bases from the 3′ end. |
| `-r, --cut-right` | off | Cut the first low-quality window and everything 3′ of it. |
| `-W, --cut-window-size <N>` | `4` | Window size for the cuts. |
| `-M, --cut-mean-quality <Q>` | `20` | Window mean-quality threshold. |

**Filtering**

| Flag | Default | Meaning |
|---|---|---|
| `-Q, --disable-quality-filtering` | off | Turn quality filtering off. |
| `-q, --qualified-quality-phred <Q>` | `15` | A base is "qualified" at/above this phred. |
| `-u, --unqualified-percent-limit <%>` | `40` | Discard a read with more than this % unqualified bases. |
| `-e, --average-qual <Q>` | `0` | Discard a read below this mean quality (0 = off). |
| `-n, --n-base-limit <N>` | `5` | Discard a read with more than N `N` bases. |
| `-L, --disable-length-filtering` | off | Turn length filtering off. |
| `-l, --length-required <N>` | `15` | Discard reads shorter than N. |
| `--length-limit <N>` | `0` | Discard reads longer than N (0 = off). |
| `-y, --low-complexity-filter` | off | Enable the low-complexity filter. |
| `-Y, --complexity-threshold <%>` | `30` | Minimum sequence complexity (percent). |

**Paired-end overlap, encoding, misc**

| Flag | Default | Meaning |
|---|---|---|
| `-c, --correction` | off | PE: correct mismatched bases in the overlap (trust the higher-quality base). |
| `--overlap-len-require <N>` | `30` | Minimum overlap length. |
| `--overlap-diff-limit <N>` | `5` | Max mismatched bases in an overlap. |
| `--overlap-diff-percent-limit <%>` | `20` | Max mismatch percent in an overlap. |
| `-6, --phred64` | off | Input qualities are phred+64 (converted to phred+33). |
| `-j, --json <FILE>` | `fastp.json` | JSON report path. |
| `--reads-to-process <N>` | `0` | Stop after N reads/pairs (0 = all). |

The processing order per read is: fixed trim → adapter → polyG → polyX →
quality cut (`-5`/`-r`/`-3`) → length cap, then the filters. Output is written
gzip-compressed when the output path ends in `.gz`.

```bash
# single-end: trim Illumina adapter + polyG, require Q20 over ≥90% of bases, ≥30 bp
cleaver fastp -i reads.fastq.gz -o clean.fastq.gz \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -g -q 20 -u 10 -l 30 -j qc.json

# 3′ quality trim with a sliding window, then length filter
cleaver fastp -i reads.fq -o clean.fq -3 -W 4 -M 25 -l 50

# paired-end: detect adapters from the overlap and error-correct it
cleaver fastp -i R1.fq.gz -I R2.fq.gz -o R1.clean.fq.gz -O R2.clean.fq.gz \
    -2 -c -g -j qc.json
```

The JSON report carries `summary.before_filtering` / `after_filtering` (reads,
bases, Q20/Q30 bases and rates, GC), `filtering_result` (passed / low-quality /
too-many-N / too-short / too-long / low-complexity), `adapter_cutting`,
`base_correction`, and the command — a faithful subset of fastp's report.

**Caveats.** This implements fastp's everyday preprocessing. It does **not** (yet)
include UMI processing, duplication/over-representation analysis, read merging, or
the HTML report. SE adapter trimming needs an explicit `-a`/`--adapter-fasta`
(there is no over-representation auto-detection); PE adapter detection works from
the read overlap (`-2`).

### `cleaver demux -i <READS> -o <DIR> -q <BARCODES>`
A pure-Rust ONT-style demultiplexer in the spirit of
[barbell](https://github.com/rickbeeloo/barbell): each read's two ends are scanned
for the best-matching barcode with a fitting (semi-global) edit-distance
alignment, the read is assigned when the match clears an identity threshold and a
margin over the runner-up, the barcode is trimmed, and reads are split into
per-barcode FASTQ files.

| Flag | Default | Meaning |
|---|---|---|
| `-i, --input <FILE>` | *required* | Reads (FASTQ; `.gz` transparent). |
| `-o, --output <DIR>` | *required* | Output folder for `<barcode>.fastq` + `unclassified.fastq`. |
| `-q, --queries <FASTA>` | *required* | Barcode FASTA (barbell's example barcode sets work directly). |
| `--window <N>` | `150` | Bases at each read end to scan for a barcode. |
| `--min-score <0..1>` | `0.75` | Minimum identity to accept a barcode match. |
| `--min-score-diff <0..1>` | `0.05` | Minimum identity margin of the best over the second-best barcode. |
| `--max-errors <N>` | — | Hard cap on edit distance (overrides `--min-score`). |
| `--no-trim` | off | Annotate/split only; keep the barcode on the read. |
| `--one-end` | off | Only scan the 5′ end (default scans both ends). |

The canonical ONT layout — forward barcode at the 5′ end, reverse-complement at
the 3′ end — is searched at both ends; a read is labelled with its best barcode
and both end barcodes are trimmed when each end clears the threshold.

```bash
# demultiplex with a native-barcode FASTA
cleaver demux -i reads.fastq.gz -o demux/ -q native_bars.fasta

# stricter acceptance and a bigger margin between candidates
cleaver demux -i reads.fq -o demux/ -q bars.fasta --min-score 0.85 --min-score-diff 0.1

# allow up to 3 edits in the barcode instead of an identity threshold
cleaver demux -i reads.fq -o demux/ -q bars.fasta --max-errors 3
```

A per-barcode read-count recap prints to stdout (classified / unclassified and a
count per barcode).

**Caveats.** This is a faithful re-implementation of barbell's
annotate→trim→split idea, not a byte-for-byte port: it uses a scalar fitting
aligner rather than barbell's SIMD (`sassy`) engine, does not bundle the ONT kit
database (supply barcodes with `-q`), and does not implement barbell's extended
templates (fusion/break detection). Barcodes are matched independently at each
end.

### `cleaver version` · `cleaver doctor` · `cleaver info`
`version` prints the banner, version, backend, and GPU-kernel status. `doctor`
runs the built-in self-test (splitting, conversion, partition, base composition,
genome N50, counting, **fastp QC**, **barcode demux**, format/GPU detection, and a
live HydraMPP task when built `--features hydra`). `info` lists version, cores,
backend, GPUs, and the supported formats/commands.

---

## 🧩 How it works

```mermaid
flowchart LR
    A[inputs] --> D{detect format}
    D -->|FASTA/FASTQ| S["streaming engine<br/>O(1) memory, record-aware"]
    D -->|SAM| T["text splitter<br/>header replicated"]
    D -->|BAM| B["noodles reader/writer<br/>BGZF finalised per chunk"]
    A --> P["work-stealing pool<br/>one task per file"]
    P --- D
    S --> O[("chunks")]
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

**`fastp`, `count` (VERSE), and `demux`** are verified for correctness (the
built-in `doctor` self-test plus the unit suite — 36 core tests covering the
overlap modes, the fastp trimming/filtering/overlap-correction primitives, and
the demux fit-aligner) but are **not** benchmarked head-to-head here: the C/C++
and SIMD originals (`fastp`, `VERSE`, `barbell`) are not reachable in this sandbox
for a fair comparison, so no timing numbers are claimed for them.

---

## 🐞 Bug audit ("find any bugs")

Found and fixed while building this:

1. **Duplicate `noodles-sam` in the dependency tree** — pinning `noodles-sam 0.50`
   while `noodles-bam 0.55` requires `0.52` put *two* copies of the crate in the
   graph, so `alignment::io::Write` and `RecordBuf: Record` came from different
   crates and didn't match (`finish`/`write_alignment_record` "not satisfied").
   Fixed by aligning to `0.52`.
2. **Silent MSRV walls** — `indexmap 2.14` and `rayon-core 1.13` require
   edition2024 / rustc ≥ 1.80. Pinned to `2.2.6` and `1.12.1`.
3. **Wrong unit-test expectation in the base counter** — the FASTA `stats` test
   asserted `C:4 / total:15 / GC:7÷13` for input `ACGT|AACC|GGTTNN`, which
   actually has `C:3 / total:14 / GC:0.5`; the kernel was right, the test was
   miscounted. Corrected the expectation (caught by running the suite, not by
   eye).

Verified correct, no bug: the FASTA engine is byte-identical to the Python
original; BAM chunks re-read as valid BAM with record counts summing exactly
(50k and 300k runs); SAM chunks each carry the header.

Honest sharp edges: a wrong/never-matching delimiter yields one giant chunk
(mitigated by the start-of-line default); memory scales with the longest *line*,
not the file; `--contains` is an O(n·m) scan; and two inputs sharing a stem
(e.g. `x.ffn` and `x.ffn.gz`) write the same chunk names into one `-o` dir —
split same-stem inputs separately.

---

## 🔧 Build

```bash
# rustc 1.75+ (no rustup needed)
cargo build --release            # -> target/release/cleaver  (rayon backend)
cargo build --release --features hydra   # + HydraMPP cluster backend
cargo build --release --features gpu     # + CUDA stats kernel (needs CUDA toolkit, rustc 1.85+)
cargo test --workspace           # 15 unit + doctests, all passing
```

MSRV is held at **1.75** by pinning `noodles-sam=0.52`, `noodles-bam=0.55`,
`noodles-bgzf=0.26`, `indexmap=2.2.6`, `rayon=1.10`, `rayon-core=1.12.1`. A fully
static binary builds with the `x86_64-unknown-linux-musl` target.

### ✅ Correctness

`#![forbid(unsafe_code)]` on the core crate (the CLI's GPU launch is the only
`unsafe`, behind `--features gpu`). Tests cover FASTA/FASTQ losslessness and
boundary behaviour, SAM header replication, SAM↔BAM round-trips, BAM-chunk
validity, FASTQ→FASTA, mapped/unmapped filtering and partitioning, base
composition / GC, **genome stats (N50/L50/N90, lengths, median)**, GTF/GFF
parsing, CIGAR intron-splitting, and **count assignment** (unique / ambiguous /
no-feature, strandedness, multimapper and MAPQ filters, and the union / strict /
nonempty overlap modes). End-to-end (verified here): all FASTA variants, FASTQ,
gzip, and multi-file parallel split reassemble byte-identical; `--partition`
routes every record by its `0x4` flag (6 mapped + 4 unmapped, lossless); `stats`
emits the genome table across files; and `count` produces an identical
gene × sample matrix and `.summary` whether run on the in-node pool or across the
HydraMPP backend with devices pinned per task.

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

**Re-implemented tools.** Cleaver's `count` modes follow **featureCounts**
(Liao, Smyth & Shi), **htseq-count** (Anders, Pyl & Huber), and **VERSE** (Zhu,
Fisher & Kim); `cleaver fastp` re-implements **fastp** (Chen, Zhou, Chen & Gu);
`cleaver demux` re-implements the approach of **barbell** (Beeloo et al.). These
are independent pure-Rust re-implementations; please cite the original tools when
their algorithms are used.

## 📄 License

Creative Commons Attribution-NonCommercial 4.0 International (**CC-BY-NC-4.0**).
Free for academic and non-commercial use; contact the author for commercial use.
