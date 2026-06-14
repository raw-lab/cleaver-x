//! Runtime configuration and the `--hydra-*` argument convention.
//!
//! Like the Python module-load hook, [`Config::from_args`] scans the process
//! arguments for HydraMPP's own flags, removes them, and hands back everything
//! else so your program's CLI parser never sees them:
//!
//! | Flag                         | Meaning                                        |
//! |------------------------------|------------------------------------------------|
//! | `--hydra-host`               | start as the coordinator (host)                |
//! | `--hydra-client <ip>`        | join the host at `<ip>` as a worker            |
//! | `--hydra-slurm <nodelist>`   | auto-wire a SLURM allocation (host + clients)  |
//! | `--hydra-cpus <n>`           | CPUs to offer (`0`/omitted ⇒ autodetect)       |
//! | `--hydra-gpus <n>`           | GPUs to offer (omitted ⇒ autodetect via nvidia-smi / CUDA_VISIBLE_DEVICES) |
//! | `--hydra-port <n>`           | control/status port (default 24515)            |
//! | `--hydra-pool <mode>`        | local CPU pool: `auto` (default) / `always` / `never` |
//! | `--hydra-stack <size>`       | worker thread stack (e.g. `512k`, `2m`); default OS ~2 MiB |
//! | `--hydra-quiet`              | silence HydraMPP's own log lines               |

use std::time::Duration;

/// Default control + status port (matches the Python original).
pub const DEFAULT_PORT: u16 = 24515;

/// Which role this process plays in the cluster.
#[derive(Debug, Clone)]
pub enum Mode {
    /// Single machine; all CPUs local.
    Local,
    /// Coordinator that accepts client connections.
    Host,
    /// Worker that connects to a host at this address.
    Client { address: String },
    /// Head of a SLURM allocation: become a host and `srun` clients on the
    /// other nodes in `nodelist` (e.g. `$SLURM_JOB_NODELIST`).
    SlurmHost { nodelist: String },
}

/// How local **CPU-bound** jobs are executed.
///
/// `cpus(0)` jobs are *always* run on their own thread regardless of this
/// setting, so I/O- and GPU-bound work keeps unbounded concurrency (this is why
/// GPU device pinning is unaffected). This only changes how jobs that reserve
/// **one or more** CPU slots are dispatched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PoolMode {
    /// Warm worker pool when the node has more than one CPU, else thread-per-job.
    /// The honest default: the pool avoids per-job thread creation on multi-core
    /// (a throughput win), while a single core — where a condvar hand-off costs
    /// more than a spawn — keeps the simpler thread-per-job path.
    #[default]
    Auto,
    /// Always route CPU-bound jobs through the warm pool.
    Always,
    /// Always spawn a fresh thread per job (the original 1.0 behaviour).
    Never,
}

/// Full runtime configuration. Build with the constructors/builders or parse
/// from the command line with [`Config::from_args`].
#[derive(Debug, Clone)]
pub struct Config {
    pub mode: Mode,
    pub num_cpus: Option<usize>,
    /// GPUs to offer. `None` autodetects (CUDA_VISIBLE_DEVICES, then nvidia-smi).
    pub num_gpus: Option<usize>,
    pub port: u16,
    pub timeout: Duration,
    pub quiet: bool,
    /// How local CPU-bound jobs are executed (see [`PoolMode`]).
    pub pool: PoolMode,
    /// Stack size (bytes) for worker threads. `None` uses the OS default (~2 MiB).
    /// Shrinking it lowers per-worker memory when a many-core node spawns many
    /// concurrent workers; raise it for tasks with deep recursion or big stack
    /// arrays. Set with [`Config::worker_stack`] or `--hydra-stack`.
    pub worker_stack_bytes: Option<usize>,
    /// Non-HydraMPP args, preserved so SLURM client spawns can replay them.
    pub(crate) passthrough: Vec<String>,
}

/// Parse a byte size like `1048576`, `512k`, `2m`, or `1g` (binary multipliers).
fn parse_size(s: &str) -> Option<usize> {
    let s = s.trim();
    let last = s.chars().last()?;
    let (num, mult) = match last.to_ascii_lowercase() {
        'k' => (&s[..s.len() - 1], 1024usize),
        'm' => (&s[..s.len() - 1], 1024 * 1024),
        'g' => (&s[..s.len() - 1], 1024 * 1024 * 1024),
        _ => (s, 1usize),
    };
    num.trim().parse::<usize>().ok().map(|n| n * mult)
}

impl Default for Config {
    fn default() -> Self {
        Config {
            mode: Mode::Local,
            num_cpus: None,
            num_gpus: None,
            port: DEFAULT_PORT,
            timeout: Duration::from_secs(10),
            quiet: false,
            pool: PoolMode::Auto,
            worker_stack_bytes: None,
            passthrough: Vec::new(),
        }
    }
}

impl Config {
    /// Local mode (the common single-machine case).
    pub fn local() -> Self {
        Config::default()
    }

    /// Host (coordinator) mode.
    pub fn host() -> Self {
        Config {
            mode: Mode::Host,
            ..Config::default()
        }
    }

    /// Client (worker) mode joining the host at `address`.
    pub fn client(address: impl Into<String>) -> Self {
        Config {
            mode: Mode::Client {
                address: address.into(),
            },
            ..Config::default()
        }
    }

    /// Override the number of CPUs offered by this node.
    pub fn num_cpus(mut self, n: usize) -> Self {
        self.num_cpus = if n == 0 { None } else { Some(n) };
        self
    }

    /// Override the number of GPUs offered by this node. `0` means "offer no
    /// GPUs"; leave unset to autodetect (CUDA_VISIBLE_DEVICES, then nvidia-smi).
    pub fn num_gpus(mut self, n: usize) -> Self {
        self.num_gpus = Some(n);
        self
    }

    /// Override the control/status port.
    pub fn port(mut self, port: u16) -> Self {
        self.port = port;
        self
    }

    /// Set the worker-thread stack size in bytes (default: OS default, ~2 MiB).
    /// Lower it to cut per-worker memory on many-core nodes; `0` resets to default.
    pub fn worker_stack(mut self, bytes: usize) -> Self {
        self.worker_stack_bytes = if bytes == 0 { None } else { Some(bytes) };
        self
    }

    /// Override the host accept window.
    pub fn timeout(mut self, timeout: Duration) -> Self {
        self.timeout = timeout;
        self
    }

    /// Silence HydraMPP's own log chatter.
    pub fn quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    /// Choose how local CPU-bound jobs are executed (see [`PoolMode`]).
    pub fn pool_mode(mut self, mode: PoolMode) -> Self {
        self.pool = mode;
        self
    }

    /// The non-HydraMPP arguments left after parsing (program name excluded),
    /// so you can feed them to your own argument parser.
    pub fn rest_args(&self) -> &[String] {
        &self.passthrough
    }

    /// Parse [`std::env::args`], extracting the `--hydra-*` flags.
    pub fn from_args() -> Self {
        Self::parse(std::env::args().skip(1))
    }

    /// Parse an explicit argument iterator (program name already excluded).
    pub fn parse<I, S>(args: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let args: Vec<String> = args.into_iter().map(Into::into).collect();
        let mut cfg = Config::default();
        let mut slurm_nodelist: Option<String> = None;
        let mut client_addr: Option<String> = None;
        let mut be_host = false;
        let mut passthrough = Vec::new();

        let mut i = 0;
        while i < args.len() {
            let arg = &args[i];
            // Support both `--flag value` and `--flag=value`.
            let (key, inline_val) = match arg.split_once('=') {
                Some((k, v)) => (k, Some(v.to_string())),
                None => (arg.as_str(), None),
            };
            let take_val = |inline: Option<String>, i: &mut usize| -> Option<String> {
                if let Some(v) = inline {
                    Some(v)
                } else if *i + 1 < args.len() {
                    *i += 1;
                    Some(args[*i].clone())
                } else {
                    None
                }
            };
            match key {
                "--hydra-host" => be_host = true,
                "--hydra-quiet" => cfg.quiet = true,
                "--hydra-slurm" => slurm_nodelist = take_val(inline_val, &mut i),
                "--hydra-client" => client_addr = take_val(inline_val, &mut i),
                "--hydra-cpus" => {
                    if let Some(v) = take_val(inline_val, &mut i) {
                        if let Ok(n) = v.parse::<usize>() {
                            cfg.num_cpus = if n == 0 { None } else { Some(n) };
                        }
                    }
                }
                "--hydra-gpus" => {
                    if let Some(v) = take_val(inline_val, &mut i) {
                        if let Ok(n) = v.parse::<usize>() {
                            cfg.num_gpus = Some(n);
                        }
                    }
                }
                "--hydra-port" => {
                    if let Some(v) = take_val(inline_val, &mut i) {
                        if let Ok(p) = v.parse::<u16>() {
                            cfg.port = p;
                        }
                    }
                }
                "--hydra-pool" => {
                    if let Some(v) = take_val(inline_val, &mut i) {
                        cfg.pool = match v.to_ascii_lowercase().as_str() {
                            "always" => PoolMode::Always,
                            "never" => PoolMode::Never,
                            _ => PoolMode::Auto,
                        };
                    }
                }
                "--hydra-stack" => {
                    if let Some(v) = take_val(inline_val, &mut i) {
                        if let Some(b) = parse_size(&v) {
                            cfg.worker_stack_bytes = if b == 0 { None } else { Some(b) };
                        }
                    }
                }
                _ => passthrough.push(arg.clone()),
            }
            i += 1;
        }

        cfg.mode = if let Some(nodelist) = slurm_nodelist {
            Mode::SlurmHost { nodelist }
        } else if let Some(address) = client_addr {
            Mode::Client { address }
        } else if be_host {
            Mode::Host
        } else {
            Mode::Local
        };
        cfg.passthrough = passthrough;
        cfg
    }
}
