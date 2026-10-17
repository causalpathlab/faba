//! Sampling progress: each chain's sweeps counted as they are done, read
//! from another thread (a bar on stderr for `faba cluster`, a view in
//! `faba modules`), and a cancel that turns the sweeps left into no-ops.
//!
//! The chains run as `run_mcmc_parallel` runs them, each model wrapped so
//! that its sweep is the model's own, then a count: the random numbers and
//! so the draws are those of the bare model.

use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig, McmcModel};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use std::io::IsTerminal;
use std::sync::atomic::{AtomicU8, AtomicUsize, Ordering};
use std::sync::OnceLock;
use std::time::{Duration, Instant};

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

/// Running, cancelled, or writing (too late to cancel).
const RUNNING: u8 = 0;
const CANCELLED: u8 = 1;
const COMMITTED: u8 = 2;

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
    /// The chains' split R-hat, once summarised.
    rhat: OnceLock<f64>,
    /// Draw the chains' sweeps on a bar on stderr while they run.
    bar: bool,
}

impl Progress {
    /// For a run of `chains` chains, none started, watched from elsewhere.
    pub fn new(chains: usize) -> Self {
        Self {
            done: (0..chains).map(|_| Counter(AtomicUsize::new(0))).collect(),
            started: (0..chains).map(|_| OnceLock::new()).collect(),
            schedule: OnceLock::new(),
            state: AtomicU8::new(RUNNING),
            rhat: OnceLock::new(),
            bar: false,
        }
    }

    /// [`Progress::new`], its chains drawn on a bar on stderr while they
    /// run, when stderr is a terminal.
    pub fn on_stderr(chains: usize) -> Self {
        Self {
            bar: true,
            ..Self::new(chains)
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
        !matches!(to, Err(COMMITTED))
    }

    /// Whether [`Progress::cancel`] stopped the run.
    pub fn cancelled(&self) -> bool {
        self.state.load(Ordering::SeqCst) == CANCELLED
    }

    /// [`Progress::cancelled`] for the chains' sweeps, which need only see
    /// a cancel soon, not at once.
    fn stopping(&self) -> bool {
        self.state.load(Ordering::Relaxed) == CANCELLED
    }

    /// [`Cancelled`] when the run was: a fit checks between the steps
    /// before its chains, which do not stop part way.
    pub fn check(&self) -> anyhow::Result<()> {
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

/// A model whose sweeps are counted on a [`Progress`], and skipped once it
/// is cancelled (as are the draws: a cancelled chain keeps none).
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

impl<'p, M: McmcModel> McmcModel for Watched<'_, 'p, M> {
    /// The chain's sweep count (none counts nowhere), and the model's state
    /// (none for a chain that started after a cancel, which never sweeps).
    type State = (Option<&'p AtomicUsize>, Option<M::State>);
    type Sample = Option<M::Sample>;
    type Result = M::Result;

    fn init(&self, rng: &mut SmallRng) -> Self::State {
        let chain = self.chain_of(rng);
        if let Some(at) = chain.and_then(|i| self.progress.started.get(i)) {
            let _ = at.set(Instant::now());
        }
        let done = chain.and_then(|i| self.progress.done.get(i)).map(|c| &c.0);
        let state = (!self.progress.cancelled()).then(|| self.model.init(rng));
        (done, state)
    }

    fn sweep(&self, (done, st): &mut Self::State, rng: &mut SmallRng) {
        let Some(st) = st.as_mut().filter(|_| !self.progress.stopping()) else {
            return;
        };
        self.model.sweep(st, rng);
        if let Some(done) = done {
            done.fetch_add(1, Ordering::Relaxed);
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
/// draws), their sweeps counted on `progress`, and drawn on a bar on
/// stderr if it is [`Progress::on_stderr`]. [`Cancelled`] when `progress`
/// was cancelled.
pub fn run<M>(
    model: &M,
    config: &McmcConfig,
    chains: usize,
    progress: &Progress,
) -> anyhow::Result<Vec<M::Result>>
where
    M: McmcModel + Sync,
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
    // Off a terminal indicatif hides the bar: no thread to draw it.
    let out = if progress.bar && std::io::stderr().is_terminal() {
        let (stop, stopped) = std::sync::mpsc::channel::<()>();
        std::thread::scope(|s| {
            s.spawn(|| draw_bar(progress, stopped));
            let out = run_mcmc_parallel(&watched, config, chains);
            drop(stop);
            out
        })
    } else {
        run_mcmc_parallel(&watched, config, chains)
    };
    progress.check()?;
    Ok(out)
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
            bar.set_message(if sampling {
                format!("sweeps over {chains} chains, sampling")
            } else {
                format!("sweeps over {chains} chains, burn-in to {warmup}")
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
