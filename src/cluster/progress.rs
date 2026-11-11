//! A fit's progress: each restart's steps (its k-means start's
//! iterations, then its epochs) counted as they are done, and its ELBO at
//! each epoch, read from another thread (by `faba cluster`: a view with the
//! traces on a terminal, else a bar on stderr), and a cancel that the fits
//! see between their steps.

use std::io::IsTerminal;
use std::sync::atomic::{AtomicBool, AtomicU8, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, OnceLock};
use std::time::{Duration, Instant};

pub mod advice;
pub mod spark;
pub mod view;

/// The fit was cancelled: nothing was written.
#[derive(Debug)]
pub struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("cancelled: the fit stopped and nothing was written")
    }
}

impl std::error::Error for Cancelled {}

/// Running, cancelled, writing (too late to cancel), or stopped early by
/// an interrupt (the restarts stop, and what they reached is written).
const RUNNING: u8 = 0;
const CANCELLED: u8 = 1;
const COMMITTED: u8 = 2;
const STOPPED: u8 = 3;

/// A restart's step count on a cache line of its own, so restarts on
/// other cores do not contend for it.
#[repr(align(128))]
struct Counter(AtomicUsize);

/// Steps done per restart of a run, its schedule once it starts, and
/// whether it was cancelled; shared between the restarts and a watcher.
/// (A restart is a chain in what the browser shares with the MCMC fits'
/// records of faba before 0.19.)
pub struct Progress {
    /// Steps done by each restart, by its index (its seed's offset).
    done: Vec<Counter>,
    /// When each restart started, by the same index.
    started: Vec<OnceLock<Instant>>,
    /// Steps per restart and its k-means start's.
    schedule: OnceLock<(usize, usize)>,
    state: AtomicU8,
    /// A cancel or a stop was asked for once the fit was writing, too late
    /// for it.
    late: AtomicBool,
    /// Draw the restarts' steps on a bar on stderr while they run.
    bar: bool,
    /// Each restart's ELBO so far, by step.
    traces: Vec<Mutex<Vec<(usize, f64)>>>,
    /// The process's interrupt flag (the first Ctrl+C or SIGTERM sets it,
    /// see `legume_numeric::matrix::stop`), when this run is the command
    /// line's: set, the restarts stop early and what they reached is
    /// written.
    interrupt: Option<Arc<AtomicBool>>,
    /// The fit's verdict, once it ends.
    optimised: OnceLock<advice::Verdict>,
}

/// The line under a fit by variational inference's traces, its k-means
/// start `warmup` steps.
pub(in crate::cluster) fn vi_caption(warmup: usize) -> String {
    format!(" │ the k-means start ends at step {warmup}; then each restart's epochs")
}

/// What a fit by variational inference traces.
pub(in crate::cluster) const ELBO: &str = "ELBO";

impl Progress {
    /// For a run of `chains` restarts, none started, watched from
    /// elsewhere.
    pub fn new(chains: usize) -> Self {
        Self {
            done: (0..chains).map(|_| Counter(AtomicUsize::new(0))).collect(),
            started: (0..chains).map(|_| OnceLock::new()).collect(),
            schedule: OnceLock::new(),
            state: AtomicU8::new(RUNNING),
            late: AtomicBool::new(false),
            bar: false,
            traces: (0..chains).map(|_| Mutex::default()).collect(),
            interrupt: None,
            optimised: OnceLock::new(),
        }
    }

    /// Restart `i`'s ELBO so far, by step.
    pub fn trace(&self, i: usize) -> Vec<(usize, f64)> {
        self.trace_of(i).clone()
    }

    /// Restart `i`'s trace, held.
    fn trace_of(&self, i: usize) -> MutexGuard<'_, Vec<(usize, f64)>> {
        self.traces[i].lock().expect("a chain's trace")
    }

    /// Every restart's [`Progress::trace`].
    pub fn traces(&self) -> Vec<Vec<(usize, f64)>> {
        (0..self.chains()).map(|i| self.trace(i)).collect()
    }

    /// The fit's verdict, once it ends (see [`Progress::optimised`]).
    pub fn verdict(&self) -> Option<advice::Verdict> {
        self.optimised.get().cloned()
    }

    /// [`Progress::new`], for the command line's own run: the
    /// process's interrupt flag stops its restarts early.
    pub fn interruptible(chains: usize) -> Self {
        Self {
            interrupt: Some(legume_numeric::matrix::stop::stop_flag()),
            ..Self::new(chains)
        }
    }

    /// [`Progress::interruptible`], its restarts drawn on a bar on stderr
    /// while they run, when stderr is a terminal.
    pub fn on_stderr(chains: usize) -> Self {
        Self {
            bar: true,
            ..Self::interruptible(chains)
        }
    }

    /// Restarts of the run.
    pub fn chains(&self) -> usize {
        self.done.len()
    }

    /// Steps done by restart `i`.
    pub fn done(&self, i: usize) -> usize {
        self.done[i].0.load(Ordering::Relaxed)
    }

    /// The time since restart `i` started; none before it does.
    pub fn running(&self, i: usize) -> Option<Duration> {
        self.started[i].get().map(Instant::elapsed)
    }

    /// Steps per restart and its k-means start's; none before the fit
    /// begins.
    pub fn schedule(&self) -> Option<(usize, usize)> {
        self.schedule.get().copied()
    }

    /// Stop the restarts at their next step. False once the fit is writing
    /// its tables, which cannot be stopped.
    pub fn cancel(&self) -> bool {
        let to =
            self.state
                .compare_exchange(RUNNING, CANCELLED, Ordering::SeqCst, Ordering::SeqCst);
        self.too_late(to)
    }

    /// Whether the change of state `to` came before the fit was writing;
    /// noted when it did not (see [`Progress::asked_late`]).
    fn too_late(&self, to: Result<u8, u8>) -> bool {
        let late = matches!(to, Err(COMMITTED));
        if late {
            self.late.store(true, Ordering::SeqCst);
        }
        !late
    }

    /// Whether a cancel or a stop was asked for too late, the fit writing
    /// already: it is written, and whoever runs another fit after it (see
    /// `faba cluster --within`) is not to.
    pub fn asked_late(&self) -> bool {
        self.late.load(Ordering::SeqCst)
    }

    /// Whether [`Progress::cancel`] stopped the run.
    pub fn cancelled(&self) -> bool {
        self.state.load(Ordering::SeqCst) == CANCELLED
    }

    /// Whether the restarts are to stop: cancelled, or stopped early (an
    /// interrupt seen here first stops them). Read between their steps,
    /// which need only see it soon, not at once.
    pub fn stopping(&self) -> bool {
        match self.state.load(Ordering::Relaxed) {
            CANCELLED | STOPPED => true,
            RUNNING if self.interrupted() => {
                self.stop_early();
                true
            }
            _ => false,
        }
    }

    /// Whether this run's interrupt flag is set: the process was
    /// interrupted, now or before this run (the flag stays set).
    pub fn interrupted(&self) -> bool {
        self.interrupt
            .as_ref()
            .is_some_and(|f| f.load(Ordering::Relaxed))
    }

    /// Stop the restarts early, keeping what they reached: what the first
    /// interrupt does (see [`Progress::stopped_early`]).
    pub fn stop_early(&self) {
        let to = self
            .state
            .compare_exchange(RUNNING, STOPPED, Ordering::SeqCst, Ordering::SeqCst);
        self.too_late(to);
    }

    /// When the restarts were stopped early, the steps each finished; the
    /// fit writes what they reached, and its record says so.
    pub fn stopped_early(&self) -> Option<Vec<usize>> {
        (self.state.load(Ordering::SeqCst) == STOPPED)
            .then(|| (0..self.chains()).map(|i| self.done(i)).collect())
    }

    /// [`Cancelled`] when the run was: a fit checks between the steps
    /// before its restarts.
    pub fn check(&self) -> anyhow::Result<()> {
        // An interrupt before the restarts run has nothing to keep: it
        // cancels.
        if self.interrupted() {
            self.cancel();
        }
        if self.cancelled() {
            return Err(Cancelled.into());
        }
        Ok(())
    }

    /// Before anything is written: [`Cancelled`] when the run was, else
    /// from here on a cancel is refused.
    fn commit(&self) -> anyhow::Result<()> {
        match self
            .state
            .compare_exchange(RUNNING, COMMITTED, Ordering::SeqCst, Ordering::SeqCst)
        {
            Err(CANCELLED) => Err(Cancelled.into()),
            _ => Ok(()),
        }
    }

    /// The fit begins: `total` steps a restart, the first `warmup` its
    /// k-means start's iterations, the rest its epochs.
    pub fn begin_vi(&self, total: usize, warmup: usize) {
        let _ = self.schedule.set((total, warmup));
    }

    /// Restart `i` at step `n` at least, started now if it was not.
    pub fn reach(&self, i: usize, n: usize) {
        let _ = self.started[i].get_or_init(Instant::now);
        self.done[i].0.fetch_max(n, Ordering::Relaxed);
    }

    /// A point of restart `i`'s trace: in place of its last when that is
    /// at the same step (an estimate made exact).
    pub fn note(&self, i: usize, at: usize, value: f64) {
        let mut trace = self.trace_of(i);
        match trace.last_mut() {
            Some(last) if last.0 == at => last.1 = value,
            _ => trace.push((at, value)),
        }
    }

    /// Once the fit ends, before anything is written: its `verdict` noted, and commit to writing ([`Cancelled`]
    /// when the run was).
    pub fn optimised(&self, verdict: advice::Verdict) -> anyhow::Result<()> {
        let _ = self.optimised.set(verdict);
        self.commit()
    }

    /// Steps done by every restart together.
    fn sum(&self) -> usize {
        (0..self.chains()).map(|i| self.done(i)).sum()
    }
}

/// `work`, with `progress` drawn on a bar on stderr meanwhile when it is
/// [`Progress::on_stderr`] (off a terminal indicatif hides the bar: no
/// thread to draw it).
pub fn barred<T: Send>(progress: &Progress, work: impl FnOnce() -> T + Send) -> T {
    if progress.bar && std::io::stderr().is_terminal() {
        let (stop, stopped) = std::sync::mpsc::channel::<()>();
        std::thread::scope(|s| {
            s.spawn(|| draw_bar(progress, stopped));
            let out = work();
            drop(stop);
            out
        })
    } else {
        work()
    }
}

/// How often the bar is brought up to date.
const BAR_TICK: Duration = Duration::from_millis(200);

/// The bar of `progress`'s steps over every restart, its k-means start or
/// its epochs, until `stopped` hangs up.
fn draw_bar(progress: &Progress, stopped: std::sync::mpsc::Receiver<()>) {
    use std::sync::mpsc::RecvTimeoutError::Timeout;
    let mut bar = None;
    // Whether the message shown says epochs: set again only when it flips.
    let mut shown = None;
    while let Err(Timeout) = stopped.recv_timeout(BAR_TICK) {
        let Some((total, warmup)) = progress.schedule() else {
            continue;
        };
        let chains = progress.chains();
        let bar =
            bar.get_or_insert_with(|| crate::common::new_progress_bar((chains * total) as u64));
        let least = (0..chains).map(|i| progress.done(i)).min().unwrap_or(0);
        let epochs = least >= warmup;
        if shown != Some(epochs) {
            shown = Some(epochs);
            bar.set_message(if epochs {
                format!("epochs over {chains} restarts")
            } else {
                format!("k-means starts of {chains} restarts")
            });
        }
        bar.set_position(progress.sum() as u64);
    }
    if let Some(bar) = bar {
        bar.finish_and_clear();
    }
}

#[cfg(test)]
mod tests;
