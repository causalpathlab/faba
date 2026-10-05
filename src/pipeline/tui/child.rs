//! The child process `faba run` starts: run it to its end where it can be
//! stopped, and follow its log a line at a time.
//!
//! A child's stderr is a pseudo-terminal where there is one, so the
//! progress bars it draws (indicatif hides them from a pipe) reach us:
//! each frame becomes a [`Progress`], the rest are lines of its log.

use std::process::{Command, Stdio};

/// Stops a run of child commands: the first stop interrupts the one
/// running, the second kills it, and none after it starts. A stop asked
/// for before a child is registered is replayed on registration.
#[derive(Default)]
pub struct Stopper {
    /// Stops asked for so far (saturating); changed only under `child`'s lock.
    asked: std::sync::atomic::AtomicU8,
    child: std::sync::Mutex<Option<std::process::Child>>,
}

/// Do what the `n`th stop asks of `child`: interrupt it the first time,
/// kill it after.
fn deliver(child: &mut std::process::Child, n: u8) {
    #[cfg(unix)]
    if n == 1 {
        let pid = i32::try_from(child.id())
            .ok()
            .and_then(rustix::process::Pid::from_raw);
        if let Some(pid) = pid {
            let _ = rustix::process::kill_process(pid, rustix::process::Signal::INT);
            return;
        }
    }
    let _ = child.kill();
}

impl Stopper {
    pub fn stop(&self) {
        use std::sync::atomic::Ordering;
        let Ok(mut c) = self.child.lock() else { return };
        let n = self.asked.load(Ordering::SeqCst).saturating_add(1);
        self.asked.store(n, Ordering::SeqCst);
        if let Some(child) = c.as_mut() {
            deliver(child, n);
        }
    }

    /// Make `child` the one `stop` reaches, and replay what was asked
    /// before it was here: one request interrupts it, more kill it.
    fn register(&self, child: std::process::Child) {
        let Ok(mut c) = self.child.lock() else { return };
        let child = c.insert(child);
        let n = self.asked.load(std::sync::atomic::Ordering::SeqCst);
        if n > 0 {
            deliver(child, n.min(2));
        }
    }

    pub fn is_stopped(&self) -> bool {
        self.asked.load(std::sync::atomic::Ordering::SeqCst) > 0
    }
}

/// Why [`run_one`] did not finish well.
#[derive(Debug, PartialEq, Eq)]
pub enum Failed {
    /// `stopper` was stopped.
    Stopped,
    /// The command did not start.
    Start(String),
    /// It ended badly: the reason it gave last.
    Exit(String),
}

/// What a child said on its stderr.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Said {
    /// A line of its log.
    Line(String),
    /// A frame of a progress bar.
    Progress(Progress),
}

/// Where a progress bar stands: `pos` of `len` (0 for a spinner), and
/// what it counts.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Progress {
    pub pos: u64,
    pub len: u64,
    pub what: String,
}

/// The columns a child's terminal has: wide enough that a bar and its
/// label are not cut.
#[cfg(unix)]
const COLUMNS: u16 = 160;

/// A pseudo-terminal: our end to read, and the child's to write its
/// stderr to.
#[cfg(unix)]
fn terminal() -> std::io::Result<(std::fs::File, std::os::fd::OwnedFd)> {
    use rustix::fs::{Mode, OFlags};
    use rustix::io::{fcntl_setfd, FdFlags};
    use rustix::pty::{grantpt, openpt, ptsname, unlockpt, OpenptFlags};
    let ours = openpt(OpenptFlags::RDWR | OpenptFlags::NOCTTY)?;
    fcntl_setfd(&ours, FdFlags::CLOEXEC)?;
    grantpt(&ours)?;
    unlockpt(&ours)?;
    let name = ptsname(&ours, Vec::new())?;
    let theirs = rustix::fs::open(
        name.as_c_str(),
        OFlags::RDWR | OFlags::NOCTTY | OFlags::CLOEXEC,
        Mode::empty(),
    )?;
    rustix::termios::tcsetwinsize(
        &theirs,
        rustix::termios::Winsize {
            ws_row: 50,
            ws_col: COLUMNS,
            ws_xpixel: 0,
            ws_ypixel: 0,
        },
    )?;
    Ok((std::fs::File::from(ours), theirs))
}

/// Where the child writes its stderr, and how we read it: a terminal
/// where we can make one, else a pipe.
fn stderr_for(command: &mut Command) -> Option<std::fs::File> {
    #[cfg(unix)]
    if let Ok((ours, theirs)) = terminal() {
        command.stderr(Stdio::from(theirs));
        return Some(ours);
    }
    command.stderr(Stdio::piped());
    None
}

/// Run `command` to its end where `stopper` can kill it, what it says on
/// its stderr to `each`.
pub fn run_one(
    mut command: Command,
    stopper: &Stopper,
    each: impl FnMut(Said),
) -> Result<(), Failed> {
    let program = std::path::Path::new(command.get_program())
        .file_name()
        .map_or_else(|| "faba".into(), |n| n.to_string_lossy().into_owned());
    let terminal = stderr_for(&mut command);
    let mut child = command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .spawn()
        .map_err(|e| Failed::Start(format!("cannot run {program}: {e}")))?;
    // Our copy of the child's end closes, so its end is the end of the log.
    drop(command);
    let log: Option<Box<dyn std::io::Read + Send>> = match terminal {
        Some(t) => Some(Box::new(t)),
        None => child
            .stderr
            .take()
            .map(|e| Box::new(e) as Box<dyn std::io::Read + Send>),
    };
    // Where `stop` can reach it, with any stop asked for meanwhile.
    stopper.register(child);
    let last = follow(log, each);
    let status = stopper
        .child
        .lock()
        .ok()
        .and_then(|mut c| c.take())
        .map(|mut c| c.wait());
    match status {
        // A stop that came after a good end changes nothing.
        Some(Ok(s)) if s.success() => Ok(()),
        _ if stopper.is_stopped() => Err(Failed::Stopped),
        Some(Err(e)) => Err(Failed::Exit(e.to_string())),
        _ => Err(Failed::Exit(
            last.strip_prefix("Error: ").unwrap_or(&last).to_string(),
        )),
    }
}

/// Read a child's stderr to the end, in pieces split at `\r` and `\n`
/// (a progress bar redraws after a `\r`, with no new line), each as
/// [`log_line`] trims it: a bar frame as a [`Said::Progress`], anything
/// else as a [`Said::Line`]. Returns the last line.
pub fn follow(log: Option<impl std::io::Read>, mut each: impl FnMut(Said)) -> String {
    let mut last = String::new();
    let Some(mut log) = log else {
        return last;
    };
    let mut said = |piece: &[u8]| {
        let line = log_line(&String::from_utf8_lossy(piece));
        if line.is_empty() {
            return;
        }
        match progress_of(&line) {
            Some(p) => each(Said::Progress(p)),
            None => {
                last.clone_from(&line);
                each(Said::Line(line));
            }
        }
    };
    let mut pending: Vec<u8> = Vec::new();
    let mut buf = [0u8; 8192];
    loop {
        let n = match log.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => n,
            Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
            // A terminal whose child is gone reads as an error.
            Err(_) => break,
        };
        for &b in &buf[..n] {
            if b == b'\r' || b == b'\n' {
                said(&pending);
                pending.clear();
            } else {
                pending.push(b);
            }
        }
    }
    said(&pending);
    last
}

/// Frames of the spinners the workspace draws.
const SPINNER: &str = "⠁⠂⠄⡀⢀⠠⠐⠈";

/// A progress bar's frame, its elapsed time already gone: the bar of `#`
/// and `-`, `pos/len`, the time left in brackets, then what it counts. A
/// spinner's frame starts with one of its ticks.
pub fn progress_of(line: &str) -> Option<Progress> {
    let (head, rest) = line.split_once(char::is_whitespace)?;
    if head.chars().all(|c| SPINNER.contains(c)) {
        return Some(Progress {
            what: rest.trim().to_string(),
            ..Progress::default()
        });
    }
    if !head.chars().all(|c| c == '#' || c == '-') {
        return None;
    }
    let rest = rest.trim_start();
    let (count, rest) = rest.split_once(char::is_whitespace).unwrap_or((rest, ""));
    let (pos, len) = count.split_once('/')?;
    let what = rest.trim_start();
    let what = match what.strip_prefix('(') {
        Some(w) => w.split_once(')').map_or("", |(_, after)| after),
        None => what,
    };
    Some(Progress {
        pos: pos.parse().ok()?,
        len: len.parse().ok()?,
        what: what.trim().to_string(),
    })
}

/// A log line without terminal colours or its "[time LEVEL module] "
/// prefix: the message is what matters on a status line.
pub fn log_line(line: &str) -> String {
    let mut plain = String::with_capacity(line.len());
    let mut chars = line.chars();
    while let Some(c) = chars.next() {
        if c != '\u{1b}' {
            plain.push(c);
        } else if chars.next() == Some('[') {
            // A CSI sequence runs to its final byte.
            for c in chars.by_ref() {
                if ('@'..='~').contains(&c) {
                    break;
                }
            }
        }
    }
    let t = plain.trim();
    match t.split_once("] ") {
        Some((head, msg)) if head.starts_with('[') => msg.trim().to_string(),
        _ => t.to_string(),
    }
}

#[cfg(test)]
#[path = "tests/child.rs"]
mod tests;
