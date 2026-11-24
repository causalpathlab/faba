//! What a fit's restarts say of it, wherever it is looked at: whether its
//! best restart's ELBO settled, and one line of advice naming the flags to
//! run it again with.

use crate::cluster::vi::Unit;

/// The verdict on a fit by variational inference: whether its best restart
/// settled (its ELBO within the tolerance) in `epochs` of the `asked` at
/// most, or was stopped early; how many restarts ended, and how far below
/// the best the worst ended, per `unit`.
#[derive(Debug, Clone, PartialEq)]
pub struct Verdict {
    converged: bool,
    stopped: bool,
    epochs: usize,
    asked: usize,
    restarts: usize,
    spread: f64,
    unit: Unit,
}

impl Verdict {
    /// The verdict on a fit whose best restart `converged` (its ELBO
    /// settled) in `epochs` of the `asked` epochs at most, or was `stopped`
    /// early; `restarts` in all, the worst ending `spread` below the best
    /// per `unit`.
    pub fn vi(
        converged: bool,
        stopped: bool,
        (epochs, asked): (usize, usize),
        restarts: usize,
        (spread, unit): (f64, Unit),
    ) -> Self {
        Self {
            converged,
            stopped,
            epochs,
            asked,
            restarts,
            spread,
            unit,
        }
    }

    /// Whether the fit can be read as it is: its best restart settled.
    pub fn converged(&self) -> bool {
        self.converged
    }

    /// How the best restart ended, in a few words, and the restarts'
    /// spread.
    pub fn said(&self) -> String {
        let epochs = self.epochs;
        let settled = match (self.converged, self.stopped) {
            (true, _) => format!("ELBO settled in {epochs} epochs"),
            (false, true) => format!("stopped early after {epochs} epochs"),
            (false, false) => format!("ELBO still rising after {epochs} epochs"),
        };
        match self.restarts {
            0 | 1 => settled,
            n => format!(
                "{settled}; {n} restarts, the worst {:.4} per {} below the best",
                self.spread, self.unit
            ),
        }
    }

    /// One line on what to do, naming the flags.
    pub fn advice(&self) -> String {
        let asked = self.asked;
        match (self.converged, self.stopped) {
            (true, _) => {
                "nothing to rerun: the fit can be read as it is (the best restart's tables)"
                    .to_string()
            }
            (false, true) => {
                format!("stopped before it settled: run again to let it finish (--epochs {asked})")
            }
            (false, false) => format!("run again with --epochs {}", 2 * asked.max(1)),
        }
    }

    /// `ELBO ok` or `ELBO check`, for a header.
    pub fn badge(&self) -> String {
        let word = if self.converged { "ok" } else { "check" };
        format!("ELBO {word}")
    }
}

#[cfg(test)]
mod tests;
