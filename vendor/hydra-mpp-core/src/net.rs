//! Length-framed message transport, mirroring the Python `net.py`:
//! a 4-byte big-endian length prefix followed by the payload bytes.

use std::io::{ErrorKind, Read, Write};
use std::net::TcpStream;

use crate::error::{HydraError, Result};
use crate::wire::Wire;

/// Hard cap on a single framed message (2 GiB). A 4-byte length prefix can
/// otherwise claim up to ~4 GiB; a corrupt or hostile prefix would then trigger
/// an enormous `vec![0u8; len]` allocation (an instant OOM / DoS). Legitimate
/// task payloads are far smaller than this, so anything larger is rejected as a
/// protocol error rather than allocated.
const MAX_FRAME_BYTES: usize = 2 * 1024 * 1024 * 1024;

/// Send one framed [`Wire`] message: `[u32 BE length][body]`. Without the
/// `compression` feature the body is the raw bincode payload (the original
/// format). With it, the body is `[1-byte flag][maybe-LZ4 payload]`.
pub fn send_msg(stream: &mut TcpStream, msg: &Wire) -> Result<()> {
    let payload = bincode::serialize(msg).map_err(HydraError::serde)?;
    if payload.len() > MAX_FRAME_BYTES {
        return Err(HydraError::Protocol(format!(
            "outgoing message of {} bytes exceeds the {MAX_FRAME_BYTES}-byte frame cap",
            payload.len()
        )));
    }
    let body = frame_payload(payload);
    let len = body.len() as u32;
    stream.write_all(&len.to_be_bytes())?;
    stream.write_all(&body)?;
    stream.flush()?;
    Ok(())
}

/// Body framing without compression: identity (the original wire format).
#[cfg(not(feature = "compression"))]
#[inline]
fn frame_payload(payload: Vec<u8>) -> Vec<u8> {
    payload
}

/// Body framing with compression: prefix a 1-byte flag (`0` = raw, `1` = LZ4)
/// and compress payloads at/above the threshold — but only keep the compressed
/// form if it is actually smaller.
#[cfg(feature = "compression")]
fn frame_payload(payload: Vec<u8>) -> Vec<u8> {
    const COMPRESS_THRESHOLD: usize = 4096;
    if payload.len() >= COMPRESS_THRESHOLD {
        let compressed = lz4_flex::compress_prepend_size(&payload);
        if compressed.len() + 1 < payload.len() {
            let mut out = Vec::with_capacity(compressed.len() + 1);
            out.push(1u8);
            out.extend_from_slice(&compressed);
            return out;
        }
    }
    let mut out = Vec::with_capacity(payload.len() + 1);
    out.push(0u8);
    out.extend_from_slice(&payload);
    out
}

/// Receive one framed [`Wire`] message. Returns `Ok(None)` on a clean EOF
/// (peer disconnected), `Err` on a real I/O or decode error.
pub fn recv_msg(stream: &mut TcpStream) -> Result<Option<Wire>> {
    let mut len_buf = [0u8; 4];
    if !read_exact_eof(stream, &mut len_buf)? {
        return Ok(None); // clean EOF before any bytes
    }
    let len = u32::from_be_bytes(len_buf) as usize;
    if len > MAX_FRAME_BYTES {
        // Refuse to allocate for an implausible/corrupt length.
        return Err(HydraError::Protocol(format!(
            "incoming frame claims {len} bytes, exceeding the {MAX_FRAME_BYTES}-byte cap; \
             dropping the connection (corrupt stream or version mismatch?)"
        )));
    }
    let mut payload = vec![0u8; len];
    if !read_exact_eof(stream, &mut payload)? {
        return Ok(None); // truncated mid-message == peer gone
    }
    let payload = unframe_payload(payload)?;
    let msg = bincode::deserialize::<Wire>(&payload).map_err(HydraError::serde)?;
    Ok(Some(msg))
}

/// Inverse of `frame_payload` without compression: identity.
#[cfg(not(feature = "compression"))]
#[inline]
fn unframe_payload(body: Vec<u8>) -> Result<Vec<u8>> {
    Ok(body)
}

/// Inverse of `frame_payload` with compression: read the flag and LZ4-decode if
/// set. The decompressed size is capped to `MAX_FRAME_BYTES` so a malicious or
/// corrupt frame cannot trigger an unbounded allocation.
#[cfg(feature = "compression")]
fn unframe_payload(body: Vec<u8>) -> Result<Vec<u8>> {
    let (&flag, rest) = body
        .split_first()
        .ok_or_else(|| HydraError::Protocol("empty frame body".into()))?;
    match flag {
        0 => Ok(rest.to_vec()),
        1 => {
            let out = lz4_flex::decompress_size_prepended(rest)
                .map_err(|e| HydraError::Protocol(format!("LZ4 decompress failed: {e}")))?;
            if out.len() > MAX_FRAME_BYTES {
                return Err(HydraError::Protocol(format!(
                    "decompressed frame of {} bytes exceeds the {MAX_FRAME_BYTES}-byte cap",
                    out.len()
                )));
            }
            Ok(out)
        }
        other => Err(HydraError::Protocol(format!(
            "unknown frame flag {other} (compression build talking to a non-compression peer?)"
        ))),
    }
}

/// Like `read_exact`, but distinguishes a clean EOF (returns `Ok(false)`) from
/// a partial read (treated as EOF too) and a genuine error (`Err`).
fn read_exact_eof(stream: &mut TcpStream, buf: &mut [u8]) -> Result<bool> {
    let mut filled = 0;
    while filled < buf.len() {
        match stream.read(&mut buf[filled..]) {
            Ok(0) => return Ok(false),
            Ok(n) => filled += n,
            Err(e) if e.kind() == ErrorKind::Interrupted => continue,
            Err(e) if e.kind() == ErrorKind::WouldBlock => {
                // Caller set a read timeout and nothing arrived yet.
                return Err(HydraError::Net(e));
            }
            Err(e) => return Err(HydraError::Net(e)),
        }
    }
    Ok(true)
}

/// Best-effort "is this socket still connected?" peek (cf. `is_connected`).
#[allow(dead_code)] // parity with the Python net API; handy for callers
pub fn is_connected(stream: &TcpStream) -> bool {
    let mut probe = [0u8; 1];
    match stream.peek(&mut probe) {
        Ok(0) => false,                                   // EOF
        Ok(_) => true,                                    // data waiting
        Err(e) if e.kind() == ErrorKind::WouldBlock => true, // open, just idle
        Err(_) => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::net::TcpListener;
    use std::thread;

    /// Round-trip a large, highly-compressible message over a real loopback
    /// socket. Passes regardless of the `compression` feature; when the feature
    /// is on, the 64 KiB payload exercises the LZ4 compress/decompress path.
    #[test]
    fn frame_roundtrip_loopback() {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let addr = listener.local_addr().unwrap();
        let args = vec![7u8; 64 * 1024]; // compresses to a few bytes under LZ4
        let sent = Wire::Task {
            id: 42,
            func: "demo".into(),
            args: args.clone(),
            num_cpus: 2,
            gpu_ids: vec![0, 1],
        };
        let server = thread::spawn(move || {
            let (mut s, _) = listener.accept().unwrap();
            recv_msg(&mut s).unwrap()
        });
        let mut client = TcpStream::connect(addr).unwrap();
        send_msg(&mut client, &sent).unwrap();
        match server.join().unwrap().expect("a framed message") {
            Wire::Task {
                id,
                func,
                args: got,
                num_cpus,
                gpu_ids,
            } => {
                assert_eq!(id, 42);
                assert_eq!(func, "demo");
                assert_eq!(got, args);
                assert_eq!(num_cpus, 2);
                assert_eq!(gpu_ids, vec![0, 1]);
            }
            other => panic!("unexpected frame: {other:?}"),
        }
    }
}
