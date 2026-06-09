# 🐉 HydraMPP (Rust)

### *Lightweight, many-headed parallel processing — from your laptop to the whole cluster.*

<div align="center">

![Rust](https://img.shields.io/badge/Rust-1.75%2B-orange?logo=rust&logoColor=white)
![crates.io](https://img.shields.io/crates/v/hydra-mpp-core?logo=rust&logoColor=white)
![License](https://img.shields.io/badge/license-CC_BY--NC_4.0-blue)
![Build](https://img.shields.io/github/actions/workflow/status/raw-lab/HydraMPP/rust.yml?branch=main)
![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS-success)
![Dependencies](https://img.shields.io/badge/deps-tiny-purple)
![Runtime](https://img.shields.io/badge/runtime-none-blueviolet)
![GPU](https://img.shields.io/badge/GPU-aware%20scheduling-76B900?logo=nvidia&logoColor=white)
![Bioinformatics](https://img.shields.io/badge/domain-bioinformatics-green)

### ⚡ Ray-like API • 🧵 native threads • 🌐 Multi-node • 🎮 GPU-aware • 🧮 Native SLURM • 🪶 Featherweight • 🦀 100% Rust

</div>

---

# 🔬 What is HydraMPP?

**HydraMPP** (Hydra **M**assive **P**arallel **P**rocessing) is a high-performance library for **distributed parallel processing**. It scales the *same binary* seamlessly from:

* 💻 A single laptop
* 🖥️ A multi-core workstation
* 🌍 A multi-node HPC cluster

This is the **100% Rust** rewrite of the original Python [HydraMPP](https://github.com/raw-lab/HydraMPP). It keeps the familiar `register` / `remote` programming model — but compiles to a single, dependency-light, `cargo install`-able static binary with **no interpreter and no async runtime**. Register a plain Rust function, fire it with `remote()`, and collect results with `wait()` / `get()`.

---

# ✨ Features

<table>
<tr>
<td width="50%">

## 🧵 Parallel Execution

* Familiar `register` + `remote` model
* Non-blocking `remote()` task submission
* Per-call `task(..).cpus(N)` reservations
* `task(..).gpus(K)` GPU reservations
* Automatic CPU detection
* CPU-slot-aware queue and dispatch
* GPU-aware dispatch + device pinning
* One native worker thread per job
* Per-task **panic isolation**

</td>
<td width="50%">

## 🌐 Distributed Computing

* Single machine → multi-node cluster
* TCP transport with length framing
* `bincode` object passing
* Native SLURM auto-configuration
* GPU autodetect (`nvidia-smi` / `--gres`)
* Live UDP status monitor
* Three roles: local • host • client
* Type-checked task boundaries

</td>
</tr>
</table>

---

# ⚡ Why HydraMPP?

| Feature                                       | HydraMPP |
| --------------------------------------------- | :------: |
| 🪶 Featherweight (tiny dependency tree)        |    ✅     |
| 🦀 100% safe Rust (`#![forbid(unsafe_code)]`) |    ✅     |
| 🧬 Ray-style API (`register` / `remote`)       |    ✅     |
| 💻 Laptop → 🌍 cluster, **same binary**        |    ✅     |
| 🧮 Native SLURM auto-wiring                    |    ✅     |
| 🎮 GPU-aware scheduling + device pinning       |    ✅     |
| 📡 Built-in live status monitor                |    ✅     |
| ⏱️ No async runtime, no GC, no interpreter     |    ✅     |
| 🔒 Per-task panic isolation                    |    ✅     |
| 📦 `cargo install` static binary               |    ✅     |
| 🔍 Small, readable, auditable codebase         |    ✅     |

---

# 🧱 Architecture

<div align="center">
  <img src="docs/hydra-mpp-workflow.svg" alt="HydraMPP task flow" width="100%">
</div>

Your program registers typed tasks and submits them with `remote()`. Each job lands in a **pending FIFO**; the **scheduler** pops jobs and places them onto any node with a free CPU slot **and** enough free GPUs — locally on a fresh worker thread, or shipped to a connected **client node** over TCP. A GPU job is pinned to concrete device ids for its lifetime, which the task reads back race-free. Finished results flow back and wake any `wait()` callers through a condition variable (no busy-polling). A **UDP status monitor** can be queried at any time.

---

# 🐉 Tech Stack

| Component           | Technology                                   |
| ------------------- | -------------------------------------------- |
| Language            | Rust ≥ 1.75 (edition 2021)                   |
| Parallelism         | `std::thread`; thread-per-job or adaptive warm pool (`PoolMode`) |
| Resource detection  | std `available_parallelism` · `nvidia-smi` / `CUDA_VISIBLE_DEVICES` (GPUs) |
| Networking          | `std::net` TCP (4-byte big-endian framing)   |
| Serialization       | `bincode` + `serde`                          |
| Cluster integration | SLURM (`scontrol` / `srun`, `--gres` aware)  |
| Status monitor      | UDP datagram                                 |
| Errors              | hand-written `Display`/`Error` (library, zero deps) · `anyhow` (CLI) |
| CLI                 | `clap`                                       |

> **Dependency philosophy.** HydraMPP pulls in only small, well-trusted crates and uses **no async runtime**. The whole thing builds on a stock toolchain with no system libraries — unlike MPI-based stacks, there is nothing to install on the cluster but the binary.

---

# 🚀 Installation

## 1️⃣ From crates.io

```bash
# the command-line tools (hydra-status, hydra-demo, hydra-doctor)
cargo install hydra-mpp-cli

# optional: LZ4-compress large wire frames (enable on all nodes or none; needs rustc >= 1.81)
# cargo install hydra-mpp-cli --features hydra-mpp-core/compression
```

```toml
# the library, in your Cargo.toml
[dependencies]
hydra-mpp-core = "1.0"
serde = { version = "1", features = ["derive"] }
```

## 2️⃣ From source

```bash
git clone https://github.com/raw-lab/HydraMPP
cd HydraMPP
cargo build --release          # binaries land in target/release/
cargo install --path crates/hydra-mpp-cli
```

> 🪶 **Smallest possible binary?** Build with the opt-in minimal profile, which
> turns on `panic = "abort"` (giving up per-task panic isolation):
> ```bash
> cargo build --profile release-tiny
> ```

---

# ⚡ Quick Start

```rust
use hydra_mpp_core::prelude::*;
use std::{thread, time::Duration};

// 1️⃣ A plain function you want to run in parallel.
fn slow_square(x: u64) -> u64 {
    thread::sleep(Duration::from_millis(200));
    x * x
}

fn main() -> anyhow::Result<()> {
    let hydra = Hydra::new();

    // 2️⃣ Register it by name (the same name is used on every node).
    hydra.register("slow_square", slow_square);

    // 3️⃣ Initialize. Config::from_args() reads any --hydra-* flags;
    //    with none, this is plain local mode using all detected CPUs.
    hydra.init(Config::from_args())?;

    // 4️⃣ Fire off jobs — submission is non-blocking and returns a job id.
    let mut ids: Vec<JobId> = (0..20u64)
        .map(|i| hydra.remote("slow_square", &i))
        .collect::<Result<_, _>>()?;

    // 5️⃣ Reserve more resources for a heavier task.
    ids.push(hydra.task("slow_square").cpus(4).submit(&99u64)?);

    // 6️⃣ Drain results as they finish.
    let mut pending = ids;
    while !pending.is_empty() {
        let (ready, still) = hydra.wait(&pending, None, 1)?;
        for id in ready {
            let r = hydra.get(id)?;
            println!("{} -> {}", r.func_name, r.value::<u64>()?);
        }
        pending = still;
    }

    hydra.shutdown();
    Ok(())
}
```

Run it — and turn the very same program into a cluster node without touching the code:

```bash
cargo run --release --example quickstart                       # local
cargo run --release --example quickstart -- --hydra-host       # coordinator
cargo run --release --example quickstart -- --hydra-client 10.0.0.5   # worker
```

> 💡 `bincode` is **not** self-describing: a task's argument type must match the
> value you submit (e.g. register `|x: u64|` and submit `&7u64`, not `&7`).

---

# 🌐 Running Modes

HydraMPP runs in **three roles**, selected by the [`Config`] you pass to `init()` (or by `--hydra-*` flags via `Config::from_args()`):

| Mode       | How to start                                   | Role                              |
| ---------- | ---------------------------------------------- | --------------------------------- |
| 💻 local   | `Config::local()`                              | Single machine, all CPUs local    |
| 🖥️ host    | `Config::host().port(24515)`                   | Coordinator that accepts clients  |
| 🛰️ client  | `Config::client("10.0.0.5").num_cpus(36)`      | Worker node that joins a host     |

```rust
// Coordinator (host) — waits up to `timeout` for clients to connect.
hydra.init(Config::host().port(24515).timeout(Duration::from_secs(10)))?;

// Worker (client) — connects, serves tasks, then exits when the host leaves.
hydra.init(Config::client("10.0.0.5").port(24515).num_cpus(36))?;
```

> A **client** call to `init()` does not return: it serves tasks until the host
> disconnects, then exits the process with status `0` — matching the Python
> client. Drive your workload from the **host** / **local** node.

---

# 📚 API Reference

| Method                                            | Description                                              |
| ------------------------------------------------- | -------------------------------------------------------- |
| `Hydra::new()`                                    | Create an un-initialized runtime                         |
| `hydra.register::<In, Out>(name, f)`              | Register a typed task under `name`                       |
| `hydra.init(config)`                              | Start in local / host / client / SLURM mode              |
| `hydra.remote(func, &arg)`                        | Queue `func(arg)`; returns a **job id** (non-blocking)   |
| `hydra.task(func).cpus(N).gpus(K).blocking(b).submit(&a)` | Customized submission (CPU/GPU reservation, block-on-dispatch) |
| `hydra.wait(&ids, timeout, max)`                  | Split a slice of ids into `(ready, pending)`             |
| `hydra.get(id)`                                   | Retrieve a finished job's result record (pops it)        |
| `hydra.get_blocking::<R>(id)`                     | Block until `id` finishes, then return its value         |
| `hydra.put(name, &value)`                         | Place a value into the queue as a finished result        |
| `hydra.nodes()`                                   | List the nodes in the cluster (incl. GPU counts/ids)     |
| `hydra.shutdown()`                                | Tear down the scheduler and worker threads               |
| `query_status(addr, port, timeout)`               | Fetch a host's status snapshot over UDP                  |
| `current_gpus()`                                  | (in-task) the GPU device ids pinned to the running task  |
| `cuda_visible_devices()`                          | (in-task) those ids as a `CUDA_VISIBLE_DEVICES` string   |

### 📋 `JobResult` record

`hydra.get(id)` returns a [`JobResult`]. Field names mirror the Python `get()`
tuple so existing mental models carry over; the return value is decoded lazily
into your own type via `.value::<R>()`.

| Index | Field          | Meaning                                            |
| :---: | -------------- | -------------------------------------------------- |
|   0   | `finished`     | `bool` — whether the job has completed             |
|   1   | `func_name`    | name of the executed task                          |
|   2   | `.value::<R>()`| the task's return value, decoded into `R`          |
|   3   | `num_cpus`     | number of CPU slots the job reserved               |
|   –   | `num_gpus`     | number of GPUs the job reserved                    |
|   –   | `gpu_ids`      | the concrete GPU device ids it was pinned to       |
|   4   | `runtime_secs` | wall-clock run time, in seconds                    |
|   5   | `hostname`     | the node the job ran on (`Option<String>`)         |

`.is_ok()` reports success; `.error()` returns the failure message for a job
whose task returned an error or panicked.

---

# 📡 Status Monitor

`hydra-status` inspects a running cluster over UDP — the Rust replacement for
the Python `hydra-status.py`.

```bash
hydra-status [ADDRESS] [PORT]
# e.g.
hydra-status 10.0.0.1 24515
```

### `hydra-doctor` — check the install

Run after installing (or in CI) to verify HydraMPP works in your environment. It prints detected CPUs/GPUs, SLURM, and the control port, then self-tests every execution path and exits non-zero on any failure:

```text
hydra-doctor            # full report
hydra-doctor --quiet    # failures + summary only
hydra-doctor && cargo test   # one-line smoke test
```

It prints total/available CPUs, the connected nodes, and the current queue, then
exits. For continuous monitoring, wrap it with `watch`:

```bash
watch -n1 hydra-status localhost
```

---

# 🧮 SLURM Integration

HydraMPP can auto-configure host and client nodes inside a SLURM allocation. Add
`--hydra-slurm $SLURM_JOB_NODELIST` to your program's invocation and the head
node expands the nodelist (`scontrol show hostnames`), resolves its own IP, and
`srun`s a copy of your binary on every other node as a client.

> 💡 Call `hydra.init()` **after** all tasks have been registered.
> Use `--hydra-cpus` to set CPUs per node; `0` or omitted auto-detects.
> Client stdout/stderr are captured under `tmp-hydra/<date>/`.

```bash
#!/bin/bash
#SBATCH --job-name=My_Slurm_Job
#SBATCH --nodes=3
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=100G
#SBATCH --time=1-0
#SBATCH -o slurm-%x-%j.out

echo "Node List : $SLURM_JOB_NODELIST"

./target/release/hydra-demo --custom-args \
    --hydra-slurm "$SLURM_JOB_NODELIST" \
    --hydra-cpus  "$SLURM_CPUS_ON_NODE"
```

## 🎮 GPU allocations under SLURM

For GPU jobs, request GPUs with `--gres=gpu:N` (or `--gpus-per-node=N`). SLURM
sets `CUDA_VISIBLE_DEVICES` on each node/step, which HydraMPP reads
automatically — so you usually need **no GPU flag at all**. To pin an explicit
count, forward `--hydra-gpus N` (the head node passes it through to every `srun`
client).

```bash
#!/bin/bash
#SBATCH --job-name=GPU_Job
#SBATCH --nodes=4
#SBATCH --tasks-per-node=1
#SBATCH --gpus-per-node=2          # or: --gres=gpu:2
#SBATCH --cpus-per-task=16
#SBATCH --time=1-0
#SBATCH -o slurm-%x-%j.out

./target/release/my_gpu_program \
    --hydra-slurm "$SLURM_JOB_NODELIST" \
    --hydra-cpus  "$SLURM_CPUS_ON_NODE"
    # GPUs auto-detected from CUDA_VISIBLE_DEVICES; add --hydra-gpus 2 to force.
```

This gives a cluster of 4 nodes × 2 GPUs; HydraMPP schedules GPU tasks across all
8 devices and pins each task to specific ids (see below).

---

# 🎮 GPU Scheduling

HydraMPP treats **GPUs as a reservable resource**, exactly like CPU slots — the
same model Ray uses with `num_gpus`. A node advertises the device ids it owns; a
task reserves `K` of them; the scheduler pins `K` *concrete* device ids to that
job for its lifetime and hands them to the task.

> **HydraMPP does not link CUDA.** It never calls a GPU API itself — it only
> *schedules* and *pins*. Your task does the real GPU work (a Rust CUDA crate,
> PyTorch in a subprocess, a GPU-accelerated aligner, …). This keeps the
> "nothing to install but the binary" promise: detection is dependency-free and
> degrades gracefully to *0 GPUs* on CPU-only nodes, so the same binary runs
> everywhere.

### Reserving and reading devices

```rust
use hydra_mpp_core::prelude::*;

fn infer(batch: Vec<f32>) -> usize {
    // The devices HydraMPP pinned to THIS task — race-free (a thread-local),
    // so concurrent GPU jobs never see each other's ids.
    let devices = hydra_mpp_core::current_gpus();          // e.g. [2]
    let cvd     = hydra_mpp_core::cuda_visible_devices();  // e.g. "2"

    // For a subprocess, hand it the assigned devices:
    //   Command::new("my_gpu_tool").env("CUDA_VISIBLE_DEVICES", cvd).status()?;
    let _ = (devices, cvd);
    batch.len()
}

let hydra = Hydra::new();
hydra.register("infer", infer);
hydra.init(Config::from_args())?;

// Reserve 1 GPU. `.cpus(0)` makes the job GPU-gated (not throttled by a CPU
// slot) — ideal for GPU-bound work. Default is 1 CPU if you omit `.cpus()`.
let id = hydra.task("infer").cpus(0).gpus(1).submit(&vec![0.0f32; 1024])?;
let used = hydra.get_blocking::<usize>(id)?;
```

### How it schedules

* **Distinct devices, automatically.** With 4 GPUs, four 1-GPU jobs run
  concurrently pinned to ids `0,1,2,3`; as each finishes its device returns to
  the pool and the next job claims it. A 2-GPU job is pinned to two ids (e.g.
  `[0,1]`). The scheduler and the in-task `current_gpus()` always agree.
* **Across nodes.** GPU counts from every node are pooled. A 4-node × 2-GPU
  SLURM job schedules over all 8 devices, shipping GPU tasks to whichever node
  has a free device.
* **Fail-fast.** A job that reserves more GPUs than *any* node provides fails
  immediately with a clear message instead of wedging the queue.
* **Autodetection** order: explicit `--hydra-gpus N` → `CUDA_VISIBLE_DEVICES`
  (set by SLURM `--gres`) → `nvidia-smi` enumeration → none.

### Try it

```bash
# Real GPUs are autodetected; simulate a GPU box anywhere with --hydra-gpus:
cargo run --release --example gpu_tasks -- --hydra-gpus 4
```

```text
node gpu01 — 16 CPUs, 4 GPUs [0, 1, 2, 3]

submitted 9 GPU jobs (≤4 run concurrently):
  job  0 → CUDA_VISIBLE_DEVICES="0" (devices [0])   [scheduler pinned [0]]
  job  1 → CUDA_VISIBLE_DEVICES="1" (devices [1])   [scheduler pinned [1]]
  job  2 → CUDA_VISIBLE_DEVICES="2" (devices [2])   [scheduler pinned [2]]
  job  3 → CUDA_VISIBLE_DEVICES="3" (devices [3])   [scheduler pinned [3]]
  job  4 → CUDA_VISIBLE_DEVICES="0" (devices [0])   [scheduler pinned [0]]
  ...
  job 999 → CUDA_VISIBLE_DEVICES="0,1" (devices [0, 1])   [scheduler pinned [0, 1]]

fail-fast check: over-subscribed GPU job rejected → job reserves 1 CPU(s) and 104 GPU(s) but no node provides that many
```

---

# 🆚 How does it compare to other Rust libraries?

Rust has excellent parallelism and networking crates — but most solve a
*different* problem than HydraMPP. HydraMPP's niche is a **featherweight,
SLURM-native task engine with Ray-like ergonomics that runs the same binary from
a laptop to a cluster**. The table below is about *fit*, not "winning":

| Capability                                  | **HydraMPP** | [rayon](https://crates.io/crates/rayon) | [tokio](https://crates.io/crates/tokio) | [rsmpi (MPI)](https://crates.io/crates/mpi) | [timely-dataflow](https://crates.io/crates/timely) | [tonic/tarpc](https://crates.io/crates/tonic) |
| ------------------------------------------- | :----------: | :----: | :----: | :---------: | :-------------: | :---------: |
| Multi-core on one machine                   |      ✅       |   ✅    |   ✅    |     ✅       |       ✅        |   build it  |
| **Multi-node distribution**                 |      ✅       |   ❌    |  manual |     ✅       |       ✅        |   build it  |
| **Native SLURM auto-wiring**                |      ✅       |   ❌    |   ❌    |   via mpirun |       ❌        |     ❌      |
| Same binary laptop → cluster                |      ✅       |   ✅¹   |   ✅¹   |     ✅       |       ✅        |     ✅¹     |
| Ray-style task model (submit → futures)     |      ✅       |   ❌²   |   ❌    |     ❌³      |       ❌³       |     ❌      |
| CPU-slot accounting / admission             |      ✅       |   ❌    |   ❌    |     ❌       |       ❌        |     ❌      |
| **GPU-aware scheduling + device pinning**   |      ✅       |   ❌    |   ❌    |     ❌⁵      |       ❌        |     ❌      |
| No external system library required         |      ✅       |   ✅    |   ✅    |     ❌⁴      |       ✅        |     ✅      |
| No async runtime required                   |      ✅       |   ✅    |   ❌    |     ✅       |       ✅        |     ❌      |
| Built-in live status monitor                |      ✅       |   ❌    |   ❌    |     ❌       |       ❌        |     ❌      |
| Lowest-overhead fine-grained local parallel |     good⁶     | **best** |  n/a  |    good     |      good       |     n/a     |

<sub>¹ only the local part; no built-in cross-node story. ² rayon is data-parallel iterators/joins, not a submit-and-collect task queue. ³ lower-level message passing / dataflow, not a task-with-result model. ⁴ requires an MPI implementation (OpenMPI/MPICH) installed and an `mpirun`/`srun` launcher. ⁵ MPI leaves GPU placement to you (rank→device by hand); HydraMPP reserves GPUs as a resource and pins device ids, like Ray's `num_gpus`. ⁶ measured: a `criterion` bench (`bench/benches/per_task_overhead.rs`) over 64 trivial tasks puts rayon at ~6 µs, crossbeam at ~13 µs, threadpool at ~19 µs, and HydraMPP at ~53 ms (warm pool) / ~139 ms (thread-per-job) — HydraMPP's serialize→schedule→dispatch contract is far heavier by design, so for fine-grained local work call rayon *inside* a task. See [`bench/BENCHMARKS.md`](bench/BENCHMARKS.md).</sub>

> **Honesty note on "best":** rayon owns the fine-grained-local cell and the
> numbers say so — HydraMPP is the heaviest per-task of any Rust task tool tested,
> the cost of being distributed/GPU/SLURM-aware. The new opt-in warm pool
> (`PoolMode`) narrows HydraMPP's own overhead (~2.6× on a 64-task burst) but does
> not approach rayon. HydraMPP's edge is distribution: vs the only other
> distributed engine here (Ray) it is ~43× faster wall, ~33× leaner RSS, and ~27×
> smaller to install (8 MB static binary vs 217 MB) on the MerCat2 kernel.

### When to reach for what

* **rayon** — the right tool for **fine-grained data parallelism on one
  machine** (`par_iter`, parallel `map`/`reduce`, divide-and-conquer joins). Its
  work-stealing scheduler has lower per-item overhead than spawning a thread per
  task. HydraMPP does **not** replace rayon — they compose: call `rayon` *inside*
  a HydraMPP task to parallelize within a node while HydraMPP fans work out
  *across* nodes.
* **tokio** — an **async runtime for I/O-bound concurrency** (servers, many
  sockets, timers). It is not a CPU task scheduler, and HydraMPP deliberately
  needs no async runtime. Use tokio when your bottleneck is I/O concurrency, not
  CPU work that needs distributing.
* **rsmpi / MPI** — the **HPC standard for tightly-coupled, low-latency message
  passing** with explicit ranks and collectives. It is more powerful and lower
  level, but requires an MPI library on every node and an `mpirun`/`srun`
  launcher, and you express your computation in terms of ranks and messages
  rather than tasks-with-results. Choose MPI for communication-heavy numerical
  kernels; choose HydraMPP for embarrassingly-parallel, task-shaped workloads you
  want running in minutes with nothing to install but the binary.
* **timely / differential-dataflow** — a **powerful distributed dataflow** model
  for streaming and iterative computations. Far more expressive for complex
  pipelines, with a correspondingly steeper learning curve. HydraMPP is a simple
  task queue, not a dataflow engine.
* **tonic / tarpc** — **RPC frameworks** (gRPC and native-Rust RPC). Great for
  request/response between services, but they give you a transport, not a
  scheduler, queue, CPU accounting, or cluster bring-up — you'd build the
  HydraMPP parts yourself on top.

> **In one line:** if you have a pile of independent Rust tasks and a SLURM
> allocation (or just a laptop today and a cluster next week), HydraMPP gets you
> running with the least ceremony. For maximum local throughput, drop to rayon;
> for communication-heavy kernels, reach for MPI.

---

# 📈 Performance

HydraMPP runs **one OS thread per admitted job** (much lighter than the Python
original, which spawned a *process* per job), and the scheduler dispatches in
**O(jobs dispatched)** per pass via a pending FIFO — so a huge backlog does not
slow placement.

A reproducible micro-benchmark is included:

```bash
cargo run --release --example bench_local
```

On a **single constrained vCPU** (CI sandbox) it measured ≈ **10,000–11,000
trivial jobs/second** end-to-end — submit → schedule → spawn thread → run →
collect — i.e. **~90–99 µs/job**, dominated by thread creation. On a real
multi-core host, throughput and speed-up scale with the number of cores; run the
example on your own hardware for representative numbers.

> No fabricated numbers here — the figure above is what the bundled benchmark
> prints in CI. Your mileage will be better on real cores.

**Per-task overhead vs other Rust task tools.** A `criterion` micro-benchmark
(`bench/benches/per_task_overhead.rs`, run with `cargo bench --bench
per_task_overhead`) pits HydraMPP against rayon, a `threadpool`, a
`crossbeam-channel` pool, and `std::thread`. For 64 trivial tasks on one vCPU:
rayon ~6 µs, crossbeam ~13 µs, threadpool ~19 µs, `std::thread` ~2.6 ms,
HydraMPP ~53 ms (warm pool) / ~139 ms (thread-per-job). **HydraMPP is the
heaviest per task** — its serialize → registry → mutex+condvar scheduler contract
does far more work than a bare pool, which is exactly what buys cross-node
distribution, GPU pinning, and SLURM wiring. **For fine-grained local work, use
rayon — ideally *inside* a HydraMPP task.**

**Adaptive warm pool (`PoolMode`).** CPU-bound jobs can run on a persistent pool
instead of a fresh thread per job (`--hydra-pool auto|always|never`). `cpus(0)`
jobs always keep their own thread, so I/O/GPU concurrency and device pinning are
unchanged. The pool cut the 64-task burst ~2.6× (139 → 53 ms) versus thread-per-job
even on a single vCPU, and helps more on multi-core where warm workers run in
parallel. `Auto` (default) turns it on only when the node has >1 CPU.

### MerCat2-style backend benchmark

The `bench/` directory holds a reproducible harness that runs MerCat2's k-mer
counting kernel across **HydraMPP (Rust), rayon, rsmpi/MPI, HydraMPP (Python),
Ray, and a serial oracle**, verifying byte-identical output and measuring wall
time, peak memory, and disk footprint. See [`bench/BENCHMARKS.md`](bench/BENCHMARKS.md);
run it with `cd bench && ./run_benchmarks.sh` (raise the size args on a real
multi-core node for true scaling numbers).

---

# 🔁 Coming from the Python HydraMPP?

| Python                                   | Rust                                                      |
| ---------------------------------------- | --------------------------------------------------------- |
| `@hydraMPP.remote` decorator             | `hydra.register("f", f)` + `hydra.remote("f", &arg)`      |
| `func.remote(args)`                      | `hydra.remote("f", &arg)`                                 |
| `func.options(num_cpus=N).remote(...)`   | `hydra.task("f").cpus(N).submit(&arg)`                    |
| `func.options(num_gpus=K).remote(...)`   | `hydra.task("f").gpus(K).submit(&arg)` (+ `current_gpus()`) |
| `hydraMPP.init(address=...)`             | `hydra.init(Config::...)`                                 |
| `--hydraMPP-slurm`, `--hydraMPP-client`  | `--hydra-slurm`, `--hydra-client`, `--hydra-gpus`         |
| `get()` returns indexable record         | `JobResult` with `.value::<R>()`                          |
| `[func.remote(a) for a in args]` + gather | `hydra.map("f", &args)?` → `Vec<Out>` (in order, one call) |
| submit a batch, collect later             | `hydra.remote_many("f", &args, cpus)?` → `Vec<JobId>` (1 lock, 1 wakeup) |
| `[f(a) for a in big_args]` (large inputs)  | `hydra.map_owned("f", args)?` → zero-copy local, no arg serialization |
| huge fan-out / large results (bound RAM)   | `hydra.for_each("f", &args, window, |i, r| …)?` → streaming, ~window results held |
| `pickle` over TCP                        | `bincode` over TCP (Rust ↔ Rust)                          |

Behavioral parity worth knowing: dispatch is **by name** (register the same
names everywhere); a departing client is dropped without requeuing its in-flight
jobs; and a job reserving more CPUs than any node provides now **fails fast**
with a clear message instead of waiting forever.

---

# 📄 License

**Creative Commons Attribution-NonCommercial (CC BY-NC 4.0)**

See the `LICENSE` file for details. For commercial licensing inquiries, contact
Richard Allen White III (rwhit101@charlotte.edu).

---

# 📖 Citation

If you use **HydraMPP** in published work, please cite:

```text
Figueroa III JL, White III RA. 2026
HydraMPP: A lightweight library for distributed massive parallel processing — threading at scale. BioRxiv
```

---

# 🤝 Contributing

We welcome:

* 🧵 Scheduling and dispatch improvements
* 🌐 Networking and fault-tolerance enhancements
* 🔐 Authentication / transport security
* 📡 Status-monitor features
* 🧪 Tests and benchmarks

Pull requests and issues are encouraged.

---

# 📞 Support

* 🐛 **Issues:** [HydraMPP Issues](https://github.com/raw-lab/HydraMPP/issues)
* 📧 **Contact:**
  * [Dr. Richard Allen White III](mailto:rwhit101@uncc.edu)
  * [Jose Luis Figueroa III](mailto:jlfiguer@uncc.edu)

---

<div align="center">

# 🐉 HydraMPP

### *Lightweight. Distributed. Parallel.*

Built with ❤️ and 🦀 in Rust.

</div>
