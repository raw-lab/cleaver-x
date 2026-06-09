//! Multi-node wiring: the host accept loop, the client serve loop, and the UDP
//! status server. Mirrors the Python `host`/`client`/`main_loop` behaviour, with
//! GPU device ids carried through the handshake and task messages.

use std::collections::HashSet;
use std::net::{TcpListener, TcpStream, UdpSocket};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use crate::error::Result;
use crate::node::Node;
use crate::registry::Registry;
use crate::scheduler::{self, CoreRef};
use crate::wire::{Outcome, Wire, STATUS_REQUEST};
use crate::{gpu, log, net, printlog};

/// Disambiguate a client hostname against names already in the cluster by
/// appending `(1)`, `(2)`, … exactly like the Python original.
fn dedup_hostname(existing: &HashSet<String>, hostname: String) -> String {
    if !existing.contains(&hostname) {
        return hostname;
    }
    // Strip a trailing "(n)" if present, then count upward.
    let (base, mut n) = match hostname.rfind('(') {
        Some(idx) if hostname.ends_with(')') => {
            let inner = &hostname[idx + 1..hostname.len() - 1];
            match inner.parse::<usize>() {
                Ok(v) => (hostname[..idx].to_string(), v),
                Err(_) => (hostname.clone(), 0),
            }
        }
        _ => (hostname.clone(), 0),
    };
    loop {
        n += 1;
        let candidate = format!("{base}({n})");
        if !existing.contains(&candidate) {
            return candidate;
        }
    }
}

/// Run the host accept window: for `timeout` seconds, accept client connections,
/// read each handshake, and register the client as a remote node.
pub(crate) fn host_accept(core: &CoreRef, port: u16, timeout: Duration) -> Result<()> {
    let listener = TcpListener::bind(("0.0.0.0", port))?;
    listener.set_nonblocking(true)?;
    printlog!("HOST INFO: waiting {}s for clients on port {port}", timeout.as_secs());

    let deadline = Instant::now() + timeout;
    while Instant::now() < deadline {
        match listener.accept() {
            Ok((stream, addr)) => {
                stream.set_nonblocking(false)?;
                match read_handshake(&stream) {
                    Some((cpus, gpus, raw_host)) => {
                        let host = {
                            let st = core.state.lock().unwrap();
                            let names: HashSet<String> =
                                st.nodes.iter().map(|n| n.hostname.clone()).collect();
                            dedup_hostname(&names, raw_host)
                        };
                        printlog!(
                            "HOST accepted {host} ({addr}) offering {cpus} CPUs, {} GPUs",
                            gpus.len()
                        );

                        let read_clone = stream.try_clone()?;
                        let write = Arc::new(Mutex::new(stream));
                        {
                            let mut st = core.state.lock().unwrap();
                            st.nodes.push(Node::remote(
                                host.clone(),
                                addr.ip().to_string(),
                                cpus,
                                gpus,
                                write,
                            ));
                        }
                        scheduler::spawn_peer_reader(core.clone(), read_clone, host);
                        core.work_cv.notify_all();
                    }
                    None => printlog!("HOST ERROR: bad handshake from {addr}"),
                }
            }
            Err(e) if e.kind() == std::io::ErrorKind::WouldBlock => {
                thread::sleep(Duration::from_millis(50));
            }
            Err(e) => printlog!("HOST ERROR: accept failed: {e}"),
        }
    }
    let n = core.state.lock().unwrap().nodes.len().saturating_sub(1);
    printlog!("HOST INFO: accept window closed; {n} client(s) joined");
    Ok(())
}

/// Read the first framed [`Wire::Handshake`] from a freshly accepted stream.
fn read_handshake(stream: &TcpStream) -> Option<(usize, Vec<usize>, String)> {
    let mut s = stream.try_clone().ok()?;
    match net::recv_msg(&mut s) {
        Ok(Some(Wire::Handshake {
            cpus,
            gpus,
            hostname,
        })) => Some((cpus, gpus, hostname)),
        _ => None,
    }
}

/// Connect to a host and serve tasks until the host disconnects. This blocks for
/// the lifetime of the client node (the caller exits the process afterwards).
pub(crate) fn client_serve(
    registry: Registry,
    address: &str,
    port: u16,
    num_cpus: usize,
    gpu_ids: Vec<usize>,
) -> Result<()> {
    printlog!("CLIENT connecting to {address}:{port}");
    let stream = TcpStream::connect((address, port))?;
    let write = Arc::new(Mutex::new(stream.try_clone()?));

    // Handshake: advertise our CPUs and GPU device ids.
    {
        let mut w = write.lock().unwrap();
        net::send_msg(
            &mut w,
            &Wire::Handshake {
                cpus: num_cpus,
                gpus: gpu_ids.clone(),
                hostname: log::hostname(),
            },
        )?;
    }
    printlog!(
        "CLIENT connected; offering {num_cpus} CPUs, {} GPUs {:?}",
        gpu_ids.len(),
        gpu_ids
    );

    let mut read = stream;
    loop {
        match net::recv_msg(&mut read) {
            Ok(Some(Wire::Task {
                id,
                func,
                args,
                num_cpus: _,
                gpu_ids,
            })) => {
                let task = registry.get(&func);
                let write_task = write.clone();
                let write_err = write.clone();
                let spawned = thread::Builder::new()
                    .name(format!("hydra-client-{id}"))
                    .spawn(move || {
                        gpu::set_assigned(&gpu_ids);
                        let start = Instant::now();
                        let outcome = match task {
                            Some(t) => t(&args),
                            None => Outcome::Err(format!("unknown task '{func}'")),
                        };
                        let msg = Wire::Result {
                            id,
                            func,
                            outcome,
                            runtime_secs: start.elapsed().as_secs_f64(),
                        };
                        if let Ok(mut w) = write_task.lock() {
                            let _ = net::send_msg(&mut w, &msg);
                        }
                    });
                if let Err(e) = spawned {
                    // Report the failure to the host instead of dropping the job.
                    if let Ok(mut w) = write_err.lock() {
                        let _ = net::send_msg(
                            &mut w,
                            &Wire::Result {
                                id,
                                func: String::new(),
                                outcome: Outcome::Err(format!(
                                    "client could not spawn worker thread: {e}"
                                )),
                                runtime_secs: 0.0,
                            },
                        );
                    }
                }
            }
            Ok(Some(_)) => {} // ignore non-task frames
            Ok(None) => {
                printlog!("CLIENT: host disconnected, shutting down");
                break;
            }
            Err(e) => {
                printlog!("CLIENT: connection error: {e}");
                break;
            }
        }
    }
    Ok(())
}

/// Spawn the UDP status server. It answers each `STATUS` datagram with a bincode
/// [`crate::wire::StatusSnapshot`]. Runs until the runtime stops.
pub(crate) fn spawn_status_server(core: CoreRef, port: u16) {
    let socket = match UdpSocket::bind(("0.0.0.0", port)) {
        Ok(s) => s,
        Err(e) => {
            printlog!("WARNING: status server unavailable on port {port}: {e}");
            return;
        }
    };
    let _ = socket.set_read_timeout(Some(Duration::from_millis(250)));
    printlog!("status server listening on UDP {port}");

    let _ = thread::Builder::new()
        .name("hydra-status".into())
        .spawn(move || {
            let mut buf = [0u8; 1024];
            loop {
                {
                    let st = core.state.lock().unwrap();
                    if !st.running {
                        break;
                    }
                }
                match socket.recv_from(&mut buf) {
                    Ok((n, addr)) => {
                        if &buf[..n] == STATUS_REQUEST || n == 0 {
                            let snap = scheduler::snapshot(&core);
                            if let Ok(bytes) = bincode::serialize(&snap) {
                                let _ = socket.send_to(&bytes, addr);
                            }
                        }
                    }
                    Err(e)
                        if e.kind() == std::io::ErrorKind::WouldBlock
                            || e.kind() == std::io::ErrorKind::TimedOut => {}
                    Err(e) => printlog!("status recv error: {e}"),
                }
            }
            printlog!("status server stopped");
        });
}
