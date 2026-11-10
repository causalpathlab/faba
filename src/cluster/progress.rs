//! Sampling progress: each chain's sweeps counted as they are done, and
//! its log joint every so many sweeps (the blockmodel's), read from
//! another thread (by `faba cluster`: a view with the traces on a
//! terminal, else a bar on stderr), and a cancel that turns the sweeps left
//! into no-ops.
//!
//! The chains run as `run_mcmc_parallel` runs them, each model wrapped so
//! that its sweep is the model's own, then a count (and the log joint,
//! read from the state without random numbers): the random numbers and so
//! the draws are those of the bare model.

use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig, McmcModel};
use rand::rngs::SmallRng;
use rand::SeedableRng;
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
/// an interrupt (the chains stop, and what they kept is written).
const RUNNING: u8 = 0;
const CANCELLED: u8 = 1;
const COMMITTED: u8 = 2;
const STOPPED: u8 = 3;

/// A chain's sweep count on a cache line of its own, so chains on other
/// cores do not contend for it.
#[repr(align(128))]
struct Counter(AtomicUsize);

/// Sweeps done per chain of a run, its schedule once it starts, and
/// whether it was cancelled; shared between the chains and a watcher.
pub struct Progress {
    /// Sweeps done by each chain, by its index (its seed's offset), as the
    /// summary numbers the chains.
    done: Vec<Counter>,
    /// When each chain started, by the same index: chains queued behind
    /// the pool's threads start late.
    started: Vec<OnceLock<Instant>>,
    /// Sweeps per chain and burn-in sweeps.
    schedule: OnceLock<(usize, usize)>,
    state: AtomicU8,
    /// A cancel or a stop was asked for once the fit was writing, too late
    /// for it.
    late: AtomicBool,
    /// The chains' split R-hat, once summarised.
    rhat: OnceLock<f64>,
    /// What pooling the chains found short of convergence, once pooled.
    unsettled: OnceLock<Vec<String>>,
    /// Draw the chains' sweeps on a bar on stderr while they run.
    bar: bool,
    /// Each chain's log joint so far, by sweep, when traced.
    traces: Vec<Mutex<Vec<(usize, f64)>>>,
    /// Note each chain's log joint every so many sweeps.
    tracing: bool,
    /// The process's interrupt flag (the first Ctrl+C or SIGTERM sets it,
    /// see `legume_numeric::matrix::stop`), when this run is the command
    /// line's: set, the chains stop early and what they kept is written.
    interrupt: Option<Arc<AtomicBool>>,
    /// What the trace holds, when not the log joint: the count model's
    /// ELBO (see [`Progress::begin_vi`]).
    quantity: OnceLock<&'static str>,
    /// A fit by variational inference's verdict, once it ends.
    optimised: OnceLock<advice::Verdict>,
}

/// The line under a fit by variational inference's traces, its k-means
/// start `warmup` steps.
pub(in crate::cluster) fn vi_caption(warmup: usize) -> String {
    format!(" │ the k-means start ends at step {warmup}; then each restart's epochs")
}

/// What a fit by variational inference traces.
pub(in crate::cluster) const ELBO: &str = "ELBO";

/// Points a chain's trace keeps about: every k-th sweep is noted, k the
/// fewest sweeps that keep it under this many.
const TRACE_POINTS: usize = 200;

/// A model whose state gives its log joint, for a chain's trace.
pub trait Traced: McmcModel {
    /// The log joint at `st`, computed without random numbers; none for a
    /// model never traced.
    fn log_joint_of(&self, _st: &Self::State) -> Option<f64> {
        None
    }
}

impl Progress {
    /// For a run of `chains` chains, none started, watched from elsewhere.
    pub fn new(chains: usize) -> Self {
        Self {
            done: (0..chains).map(|_| Counter(AtomicUsize::new(0))).collect(),
            started: (0..chains).map(|_| OnceLock::new()).collect(),
            schedule: OnceLock::new(),
            state: AtomicU8::new(RUNNING),
            late: AtomicBool::new(false),
            rhat: OnceLock::new(),
            unsettled: OnceLock::new(),
            bar: false,
            traces: (0..chains).map(|_| Mutex::default()).collect(),
            tracing: false,
            interrupt: None,
            quantity: OnceLock::new(),
            optimised: OnceLock::new(),
        }
    }

    /// [`Progress::new`], each chain's log joint traced as it sweeps.
    pub fn traced(chains: usize) -> Self {
        Self {
            tracing: true,
            ..Self::new(chains)
        }
    }

    /// Chain `i`'s log joint so far, by sweep.
    pub fn trace(&self, i: usize) -> Vec<(usize, f64)> {
        self.trace_of(i).clone()
    }

    /// Chain `i`'s trace, held.
    fn trace_of(&self, i: usize) -> MutexGuard<'_, Vec<(usize, f64)>> {
        self.traces[i].lock().expect("a chain's trace")
    }

    /// Every chain's [`Progress::trace`].
    pub fn traces(&self) -> Vec<Vec<(usize, f64)>> {
        (0..self.chains()).map(|i| self.trace(i)).collect()
    }

    /// The verdict on the chains' traces and split R-hat in the schedule
    /// they ran (see [`advice::Verdict::of`]), and once they are pooled on
    /// what that found unsettled; none until they are summarised.
    pub fn verdict(&self) -> Option<advice::Verdict> {
        if let Some(v) = self.optimised.get() {
            return Some(v.clone());
        }
        let ((total, warmup), rhat) = (self.schedule()?, self.rhat()?);
        let verdict = advice::Verdict::of(&self.traces(), warmup, total, rhat);
        Some(verdict.with_unsettled(self.unsettled.get().cloned().unwrap_or_default()))
    }

    /// Once the chains are pooled, note what that found `unsettled` (see
    /// [`advice::Verdict::with_unsettled`]).
    pub fn pooled(&self, unsettled: Vec<String>) {
        let _ = self.unsettled.set(unsettled);
    }

    /// Whether a chain's log joint is noted at its sweep `n`: every k-th,
    /// k the fewest sweeps that keep the trace within [`TRACE_POINTS`].
    fn traces_at(&self, n: usize) -> bool {
        if !self.tracing {
            return false;
        }
        let every = self
            .schedule()
            .map_or(1, |(total, _)| total.div_ceil(TRACE_POINTS).max(1));
        n.is_multiple_of(every)
    }

    /// [`Progress::traced`], for the command line's own run: the
    /// process's interrupt flag stops its chains early.
    pub fn interruptible(chains: usize) -> Self {
        Self {
            interrupt: Some(legume_numeric::matrix::stop::stop_flag()),
            ..Self::traced(chains)
        }
    }

    /// [`Progress::interruptible`], its chains drawn on a bar on stderr
    /// while they run, when stderr is a terminal.
    pub fn on_stderr(chains: usize) -> Self {
        Self {
            bar: true,
            ..Self::interruptible(chains)
        }
    }

    /// Chains of the run.
    pub fn chains(&self) -> usize {
        self.done.len()
    }

    /// Sweeps done by chain `i`.
    pub fn done(&self, i: usize) -> usize {
        self.done[i].0.load(Ordering::Relaxed)
    }

    /// The time since chain `i` started; none while it waits for a thread.
    pub fn running(&self, i: usize) -> Option<Duration> {
        self.started[i].get().map(Instant::elapsed)
    }

    /// Sweeps per chain and burn-in sweeps; none before sampling begins.
    pub fn schedule(&self) -> Option<(usize, usize)> {
        self.schedule.get().copied()
    }

    /// Stop the chains: what sweeps are left do nothing. False once the
    /// fit is writing its tables, which cannot be stopped.
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

    /// Whether the chains are to stop sweeping: cancelled, or stopped
    /// early (an interrupt seen here first stops them). Read by the
    /// chains' sweeps, which need only see it soon, not at once.
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

    /// Stop the chains early, keeping what they drew: what the first
    /// interrupt does (see [`Progress::stopped_early`]).
    pub fn stop_early(&self) {
        let to = self
            .state
            .compare_exchange(RUNNING, STOPPED, Ordering::SeqCst, Ordering::SeqCst);
        self.too_late(to);
    }

    /// When the chains were stopped early, the sweeps each finished; the
    /// fit writes what they kept, and its record says so.
    pub fn stopped_early(&self) -> Option<Vec<usize>> {
        (self.state.load(Ordering::SeqCst) == STOPPED)
            .then(|| (0..self.chains()).map(|i| self.done(i)).collect())
    }

    /// After the chains of `config`: when they were stopped early, say so
    /// in `record` (`stopped_early`: the sweeps each chain finished and
    /// those it was to run) and on stderr, or, when a chain kept fewer than
    /// the 4 draws split R-hat needs (stopped in or just after its burn-in),
    /// refuse to write: nothing to summarise.
    pub fn note_stop(
        &self,
        config: &McmcConfig,
        record: &mut crate::run_record::RunRecord,
    ) -> anyhow::Result<()> {
        let Some(done) = self.stopped_early() else {
            return Ok(());
        };
        // Draws are kept at the burn-in's end and every `thin` sweeps on.
        let kept = |d: usize| d.saturating_sub(config.warmup).div_ceil(config.thin.max(1));
        let fewest = done.iter().map(|&d| kept(d)).min().unwrap_or(0);
        anyhow::ensure!(
            fewest >= 4,
            "stopped early at sweeps {done:?} (burn-in {}): a chain kept {fewest} draws, \
             fewer than 4, so nothing was written",
            config.warmup
        );
        log::warn!(
            "stopped early at sweeps {done:?} of {}: the tables are from the draws kept so far \
             ({fewest} or more per chain), short of the run asked for",
            sweeps(config)
        );
        record.set_result(
            "stopped_early",
            &serde_json::json!({"sweeps": done, "of": sweeps(config)}),
        );
        Ok(())
    }

    /// [`Cancelled`] when the run was: a fit checks between the steps
    /// before its chains, which do not stop part way.
    pub fn check(&self) -> anyhow::Result<()> {
        // An interrupt before the chains run has no draws to keep: it
        // cancels.
        if self.interrupted() {
            self.cancel();
        }
        if self.cancelled() {
            return Err(Cancelled.into());
        }
        Ok(())
    }

    /// Once the chains are summarised to split R-hat `rhat`, before
    /// anything is written: note it, and commit to writing, [`Cancelled`]
    /// when the run was.
    pub fn summarised(&self, rhat: f64) -> anyhow::Result<()> {
        self.set_rhat(rhat);
        self.commit()
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

    /// For a fit run outside the chain engine, the count model's by
    /// variational inference: `total` steps a restart (each restart a
    /// chain), the first `warmup` its k-means start's iterations, the rest
    /// its epochs; its trace the ELBO.
    pub fn begin_vi(&self, total: usize, warmup: usize) {
        let _ = self.schedule.set((total, warmup));
        let _ = self.quantity.set(ELBO);
    }

    /// What the trace holds: `log joint`, or `ELBO` for a fit by
    /// variational inference.
    pub fn quantity(&self) -> &'static str {
        self.quantity.get().copied().unwrap_or("log joint")
    }

    /// Whether this is a fit by variational inference (see
    /// [`Progress::begin_vi`]).
    pub fn is_vi(&self) -> bool {
        self.quantity() == ELBO
    }

    /// Chain `i` at step `n` at least, started now if it was not.
    pub fn reach(&self, i: usize, n: usize) {
        let _ = self.started[i].get_or_init(Instant::now);
        self.done[i].0.fetch_max(n, Ordering::Relaxed);
    }

    /// A point of chain `i`'s trace, when traced: in place of its last
    /// when that is at the same step (an estimate made exact).
    pub fn note(&self, i: usize, at: usize, value: f64) {
        if self.tracing {
            let mut trace = self.trace_of(i);
            match trace.last_mut() {
                Some(last) if last.0 == at => last.1 = value,
                _ => trace.push((at, value)),
            }
        }
    }

    /// Once a fit by variational inference ends, before anything is
    /// written: its `verdict` noted, and commit to writing ([`Cancelled`]
    /// when the run was).
    pub fn optimised(&self, verdict: advice::Verdict) -> anyhow::Result<()> {
        let _ = self.optimised.set(verdict);
        self.commit()
    }

    /// The chains' split R-hat, once they are summarised.
    pub fn rhat(&self) -> Option<f64> {
        self.rhat.get().copied()
    }

    /// Note the chains' split R-hat.
    fn set_rhat(&self, rhat: f64) {
        let _ = self.rhat.set(rhat);
    }

    /// Sweeps done by every chain together.
    fn sum(&self) -> usize {
        (0..self.chains()).map(|i| self.done(i)).sum()
    }
}

/// Sweeps a chain of `config` runs: the burn-in, then up to its last kept
/// draw (a `thin` of 0 read as 1). The engine keeps this schedule to
/// itself, so it is mirrored here, overflow check and all.
pub fn sweeps(config: &McmcConfig) -> usize {
    match config.n_samples {
        0 => Some(config.warmup),
        n => (n - 1)
            .checked_mul(config.thin.max(1))
            .and_then(|t| t.checked_add(config.warmup + 1)),
    }
    .expect("warmup + n_samples * thin overflows usize")
}

/// A model whose sweeps are counted on a [`Progress`], its log joint noted
/// there when traced, and skipped once it is cancelled (as are the draws:
/// a cancelled chain keeps none).
struct Watched<'m, 'p, M> {
    model: &'m M,
    progress: &'p Progress,
    /// The run's seed, which the engine offsets by each chain's index.
    seed: u64,
}

impl<M> Watched<'_, '_, M> {
    /// The index of the chain whose random numbers are `rng`, as the engine
    /// seeds them, found without drawing any; none for a stranger's.
    fn chain_of(&self, rng: &SmallRng) -> Option<usize> {
        (0..self.progress.chains())
            .find(|&i| *rng == SmallRng::seed_from_u64(self.seed.wrapping_add(i as u64)))
    }
}

impl<M: Traced> McmcModel for Watched<'_, '_, M> {
    /// The chain's index (none counts nowhere), and the model's state
    /// (none for a chain that started after a cancel, which never sweeps).
    type State = (Option<usize>, Option<M::State>);
    type Sample = Option<M::Sample>;
    type Result = M::Result;

    fn init(&self, rng: &mut SmallRng) -> Self::State {
        let chain = self.chain_of(rng);
        if let Some(at) = chain.map(|i| &self.progress.started[i]) {
            let _ = at.set(Instant::now());
        }
        let state = (!self.progress.cancelled()).then(|| self.model.init(rng));
        (chain, state)
    }

    fn sweep(&self, (chain, st): &mut Self::State, rng: &mut SmallRng) {
        let Some(st) = st.as_mut().filter(|_| !self.progress.stopping()) else {
            return;
        };
        self.model.sweep(st, rng);
        let Some(i) = *chain else {
            return;
        };
        let n = self.progress.done[i].0.fetch_add(1, Ordering::Relaxed) + 1;
        if self.progress.traces_at(n) {
            if let Some(v) = self.model.log_joint_of(st) {
                self.progress.trace_of(i).push((n, v));
            }
        }
    }

    fn collect(&self, (_, st): &Self::State) -> Self::Sample {
        st.as_ref()
            .filter(|_| !self.progress.stopping())
            .map(|st| self.model.collect(st))
    }

    fn summarize(&self, samples: Vec<Self::Sample>) -> M::Result {
        self.model
            .summarize(samples.into_iter().flatten().collect())
    }
}

/// `chains` chains of `model`, as `run_mcmc_parallel` runs them (the same
/// draws), their sweeps counted on `progress` (their log joint traced if
/// it is [`Progress::traced`]), and drawn on a bar on
/// stderr if it is [`Progress::on_stderr`]. [`Cancelled`] when `progress`
/// was cancelled.
pub fn run<M>(
    model: &M,
    config: &McmcConfig,
    chains: usize,
    progress: &Progress,
) -> anyhow::Result<Vec<M::Result>>
where
    M: Traced + Sync,
    M::State: Send,
    M::Result: Send,
{
    progress.check()?;
    let _ = progress.schedule.set((sweeps(config), config.warmup));
    let watched = Watched {
        model,
        progress,
        seed: config.seed,
    };
    let out = barred(progress, || run_mcmc_parallel(&watched, config, chains));
    progress.check()?;
    Ok(out)
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

/// The bar of `progress`'s sweeps over every chain, burn-in or sampling,
/// until `stopped` hangs up.
fn draw_bar(progress: &Progress, stopped: std::sync::mpsc::Receiver<()>) {
    use std::sync::mpsc::RecvTimeoutError::Timeout;
    let mut bar = None;
    // Whether the message shown says sampling: set again only when it flips.
    let mut shown = None;
    while let Err(Timeout) = stopped.recv_timeout(BAR_TICK) {
        let Some((total, warmup)) = progress.schedule() else {
            continue;
        };
        let chains = progress.chains();
        let bar =
            bar.get_or_insert_with(|| crate::common::new_progress_bar((chains * total) as u64));
        let least = (0..chains).map(|i| progress.done(i)).min().unwrap_or(0);
        let sampling = least >= warmup;
        if shown != Some(sampling) {
            shown = Some(sampling);
            bar.set_message(match (progress.is_vi(), sampling) {
                (true, true) => format!("epochs over {chains} restarts"),
                (true, false) => format!("k-means starts of {chains} restarts"),
                (false, true) => format!("sweeps over {chains} chains, sampling"),
                (false, false) => format!("sweeps over {chains} chains, burn-in to {warmup}"),
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
