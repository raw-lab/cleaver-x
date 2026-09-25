//! GPU detection and the base-composition compute path.
//!
//! Detection shells out to `nvidia-smi`, so it needs no build dependency and
//! works wherever the NVIDIA driver is present (empty list otherwise).
//!
//! Splitting and format conversion are I/O-bound, so they have no GPU kernel.
//! Base composition / GC, by contrast, is a reduction over millions of bytes —
//! genuinely data-parallel — so `stats` offloads it to the GPU when the `gpu`
//! feature is built and a device is available, and falls back to the CPU
//! reference otherwise. Both paths fold the *same* per-byte histogram into
//! A/C/G/T/N/other (see `cleaver::compute`), so results are identical.
//!
//! The CUDA kernel uses `cudarc` 0.12, whose build needs the CUDA toolkit and
//! rustc >= 1.85, so it is OFF by default and compiled on GPU hosts with
//! `--features gpu`. The CPU path is the one exercised by the test suite.

use std::path::Path;

use anyhow::Result;

use cleaver::compute::{self, Counts};
use cleaver::Format;

pub struct Gpu {
    pub name: String,
    pub mem_mib: String,
}

/// Enumerate NVIDIA GPUs via `nvidia-smi`. Empty if none / driver absent.
pub fn detect() -> Vec<Gpu> {
    let out = std::process::Command::new("nvidia-smi")
        .args(["--query-gpu=name,memory.total", "--format=csv,noheader,nounits"])
        .output();
    match out {
        Ok(o) if o.status.success() => String::from_utf8_lossy(&o.stdout)
            .lines()
            .filter_map(|l| {
                let mut it = l.splitn(2, ',');
                let name = it.next()?.trim().to_string();
                let mem = it.next().unwrap_or("").trim().to_string();
                (!name.is_empty()).then_some(Gpu { name, mem_mib: mem })
            })
            .collect(),
        _ => Vec::new(),
    }
}

/// Whether the GPU kernel was compiled in.
pub const fn compute_enabled() -> bool {
    cfg!(feature = "gpu")
}

/// Count base composition for one file, on the GPU when the kernel is compiled
/// in and `device` is `Some`, otherwise the CPU reference. Returns the counts
/// and the backend used (`"gpu"` or `"cpu"`).
pub fn count_file_accel(path: &Path, fmt: Format, device: Option<usize>) -> Result<(Counts, &'static str)> {
    #[cfg(feature = "gpu")]
    {
        if let Some(dev) = device {
            match gpu_count_file(path, fmt, dev) {
                Ok(c) => return Ok((c, "gpu")),
                Err(e) => eprintln!("  warning: GPU count failed on device {dev} ({e:#}); using CPU"),
            }
        }
    }
    #[cfg(not(feature = "gpu"))]
    let _ = device;

    Ok((compute::count_file(path, fmt)?, "cpu"))
}

/// CUDA implementation: stream sequence bytes to the device in bounded batches,
/// build a 256-bin byte histogram with one `atomicAdd` per base, then fold the
/// bins into `Counts` exactly as `compute::count_bytes` does. Memory stays flat
/// regardless of file size. Requires the `gpu` feature (CUDA toolkit, rustc 1.85+).
#[cfg(feature = "gpu")]
fn gpu_count_file(path: &Path, fmt: Format, device: usize) -> Result<Counts> {
    use std::io::BufRead;

    use anyhow::Context;
    use cudarc::driver::{CudaDevice, CudaFunction, CudaSlice, LaunchAsync, LaunchConfig};
    use cudarc::nvrtc::compile_ptx;

    const SRC: &str = r#"
extern "C" __global__ void hist_bytes(const unsigned char *data,
                                      const unsigned long long n,
                                      unsigned int *hist) {
    unsigned long long i = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
    unsigned long long stride = (unsigned long long)blockDim.x * gridDim.x;
    for (; i < n; i += stride) atomicAdd(&hist[data[i]], 1u);
}
"#;
    const BATCH: usize = 64 << 20; // 64 MiB per host->device transfer

    let dev = CudaDevice::new(device).context("opening CUDA device")?;
    dev.load_ptx(compile_ptx(SRC).context("compiling kernel")?, "cleaver", &["hist_bytes"])
        .context("loading kernel")?;
    let func = dev.get_func("cleaver", "hist_bytes").context("kernel not found")?;
    let mut d_hist = dev.alloc_zeros::<u32>(256).context("allocating histogram")?;

    // One batch: copy to device and accumulate into the persistent histogram.
    fn launch(
        dev: &std::sync::Arc<CudaDevice>,
        func: &CudaFunction,
        batch: &mut Vec<u8>,
        d_hist: &mut CudaSlice<u32>,
    ) -> Result<()> {
        if batch.is_empty() {
            return Ok(());
        }
        let n = batch.len() as u64;
        let d_in = dev.htod_sync_copy(batch).context("host->device copy")?;
        let threads = 256u32;
        let blocks = (((n + threads as u64 - 1) / threads as u64).min(65_535)) as u32;
        let cfg = LaunchConfig {
            grid_dim: (blocks, 1, 1),
            block_dim: (threads, 1, 1),
            shared_mem_bytes: 0,
        };
        unsafe {
            func.clone().launch(cfg, (&d_in, n, &mut *d_hist)).context("kernel launch")?;
        }
        batch.clear();
        Ok(())
    }

    let (mut reader, _gz) = cleaver::open_reader(path)?;
    let mut line: Vec<u8> = Vec::with_capacity(256);
    let mut batch: Vec<u8> = Vec::with_capacity(BATCH);
    let mut lineno: u64 = 0;
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let is_seq = match fmt {
            Format::Fasta => line.first() != Some(&b'>'),
            Format::Fastq => lineno % 4 == 1,
            _ => true,
        };
        if is_seq {
            batch.extend_from_slice(&line);
        }
        if batch.len() >= BATCH {
            launch(&dev, &func, &mut batch, &mut d_hist)?;
        }
        lineno += 1;
    }
    launch(&dev, &func, &mut batch, &mut d_hist)?;

    let hist = dev.dtoh_sync_copy(&d_hist).context("device->host copy")?;
    let bin = |b: u8| hist[b as usize] as u64;
    let acgtn = bin(b'A') + bin(b'a') + bin(b'C') + bin(b'c') + bin(b'G') + bin(b'g')
        + bin(b'T') + bin(b't') + bin(b'N') + bin(b'n');
    let breaks = bin(b'\n') + bin(b'\r');
    let all: u64 = hist.iter().map(|&x| x as u64).sum();
    Ok(Counts {
        a: bin(b'A') + bin(b'a'),
        c: bin(b'C') + bin(b'c'),
        g: bin(b'G') + bin(b'g'),
        t: bin(b'T') + bin(b't'),
        n: bin(b'N') + bin(b'n'),
        other: all - acgtn - breaks,
    })
}
