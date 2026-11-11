//! The schedule the fits by variational inference (the count fit's, the
//! site blockmodel's) share: a local step for every cell, then stochastic
//! epochs over minibatches of cells (when a batch is fewer than every
//! cell; Hoffman et al. 2013), then coordinate ascent to the tolerance, at
//! most `epochs` epochs in all, then a last local step for every cell, its
//! pass the tables'. Each fit says what its steps are (see [`Steps`]).

use rand::rngs::SmallRng;
use rand::seq::SliceRandom;

/// A fit's options.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Options {
    /// Stochastic epochs before coordinate ascent (none when one batch
    /// holds every cell), each of batches of `batch_size` cells.
    pub svi_epochs: usize,
    pub batch_size: usize,
    /// Epochs in all, stochastic and coordinate ascent.
    pub epochs: usize,
    /// Converged when an epoch of coordinate ascent moves the ELBO by less
    /// than this per unit (count or read).
    pub tol: f64,
    /// The step size's decay, `rho_t = (tau0 + t)^(-forget)`, `tau0` the
    /// batches of an epoch.
    pub forget: f64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            svi_epochs: 10,
            batch_size: 1000,
            epochs: 200,
            tol: 1e-5,
            forget: 0.7,
        }
    }
}

/// What a fit tells as it goes: each epoch's ELBO at its start, by index
/// (exact; an epoch told again, its local step done, replaces the value
/// before it), whether to stop (and write what it has), and whether it was
/// cancelled (nothing to write).
pub struct Hooks<'h> {
    pub epoch: &'h (dyn Fn(usize, f64) + Sync),
    pub stop: &'h (dyn Fn() -> bool + Sync),
    pub cancelled: &'h (dyn Fn() -> bool + Sync),
}

/// A fit's steps, which [`run`] schedules; each none when it was
/// cancelled.
pub trait Steps {
    /// One epoch of coordinate ascent's local steps: what its ELBO and its
    /// global step need.
    type Pass;
    /// The fit's end: its last pass, its ELBO among it.
    type End;

    /// The units (counts or reads) the tolerance is per.
    fn units(&self) -> f64;
    /// Every cell fitted to the start's globals, before any global step.
    fn first(&mut self) -> Option<()>;
    /// The local steps of `cells`, then a stochastic global step of size
    /// `rho`, their statistics scaled by `scale`.
    fn batch(&mut self, cells: &[usize], scale: f64, rho: f64) -> Option<()>;
    /// After each stochastic epoch, before it is scored.
    fn end_svi_epoch(&mut self) {}
    /// The ELBO of the state as it is, every cell's local step scored.
    fn score(&mut self) -> Option<f64>;
    /// An epoch of coordinate ascent's local steps: the ELBO at its start,
    /// and its pass.
    fn epoch(&mut self) -> Option<(f64, Self::Pass)>;
    /// The end on `pass`, where the fit converged.
    fn converged(&mut self, pass: Self::Pass) -> Self::End;
    /// The epoch's global step from `pass`.
    fn step(&mut self, pass: Self::Pass);
    /// The last local steps, from the last global step, as the end.
    fn last(&mut self) -> Option<Self::End>;
    /// The ELBO at `end`.
    fn elbo(end: &Self::End) -> f64;
}

/// How a fit ran: its epochs, the stochastic among them, and whether it
/// converged.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ran {
    pub svi_epochs: usize,
    pub epochs: usize,
    pub converged: bool,
}

/// `steps` over `cells` cells as the schedule says (see the module), the
/// batches shuffled by `rng`, each epoch's ELBO told to `hooks`; its end
/// and how it ran. None when it was cancelled, or stopped before an epoch
/// ended.
pub fn run<S: Steps>(
    steps: &mut S,
    cells: usize,
    options: &Options,
    rng: &mut SmallRng,
    hooks: &Hooks,
) -> Option<(S::End, Ran)> {
    steps.first()?;
    let batches = cells.div_ceil(options.batch_size.max(1));
    let svi = if batches > 1 {
        options.svi_epochs.min(options.epochs)
    } else {
        0
    };
    // Batches as even as the cells allow: none a scrap scaled up.
    let batch = cells.div_ceil(batches.max(1));
    let tau0 = batches as f64;
    let (mut t, mut done) = (0usize, 0);
    let mut order: Vec<usize> = (0..cells).collect();
    'svi: for e in 1..=svi {
        order.shuffle(rng);
        for chunk in order.chunks(batch) {
            if (hooks.stop)() {
                break 'svi;
            }
            let rho = (tau0 + t as f64).powf(-options.forget);
            steps.batch(chunk, cells as f64 / chunk.len() as f64, rho)?;
            t += 1;
        }
        steps.end_svi_epoch();
        // The last stochastic epoch is scored by the next epoch's (or the
        // last) local steps, which start from it.
        if e < svi {
            (hooks.epoch)(e, steps.score()?);
        }
        done = e;
    }
    // Coordinate ascent. Each epoch's local steps score the state at its
    // start (the previous epoch's end) as they go, and the state they
    // leave; converged, the epoch ends there, its global step not taken.
    let (mut last, mut converged, mut end) = (None::<f64>, false, None);
    let units = steps.units().max(1.0);
    for e in (svi + 1)..=options.epochs {
        if (hooks.stop)() || done < e - 1 {
            break;
        }
        let (f, pass) = steps.epoch()?;
        (hooks.epoch)(e - 1, f);
        if last.is_some_and(|prev| (f - prev).abs() / units < options.tol) {
            converged = true;
            end = Some(steps.converged(pass));
            break;
        }
        last = Some(f);
        steps.step(pass);
        done = e;
    }
    if done == 0 && end.is_none() {
        return None;
    }
    let end = match end {
        Some(end) => end,
        None => steps.last()?,
    };
    (hooks.epoch)(done, S::elbo(&end));
    let ran = Ran {
        svi_epochs: svi.min(done),
        epochs: done,
        converged,
    };
    Some((end, ran))
}

#[cfg(test)]
mod tests;
