//! A fit's progress said in lines, for whoever reads its stderr off a
//! terminal (an agent through `faba mcp`, a batch job's log), where no bar
//! is drawn: every [`EVERY`], the restart fitting, its epochs and ELBO, how
//! far the run is and about how long it has left.

use super::{Progress, NONE};
use crate::tui::clock;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{Receiver, RecvTimeoutError};
use std::time::{Duration, Instant};

/// How often a fit off a terminal says where it stands.
pub(super) const EVERY: Duration = Duration::from_secs(15);

impl Progress {
    /// Epochs restart `i` has run, as its trace says: not the steps it is
    /// credited with once it settles, or when its start failed.
    fn epochs_run(&self, i: usize, warmup: usize) -> usize {
        self.trace_of(i)
            .last()
            .map_or(0, |&(at, _)| at.saturating_sub(warmup))
    }

    /// Where the fit stands, in a line: the k-means starts' iterations, or
    /// the restart fitting, its epochs done, its ELBO per unit and (in
    /// coordinate ascent) how much its last epoch moved it against the
    /// tolerance; then the share of every restart's steps done, the time
    /// since the first started, and the time left at the pace of the epochs
    /// run so far (at most: a restart that settles ends early). None before
    /// the fit begins.
    pub fn standing(&self) -> Option<String> {
        let (total, warmup) = self.schedule()?;
        let restarts = self.restarts();
        let ran = (0..restarts).filter_map(|i| self.running(i)).max()?;
        let share = 100.0 * self.sum() as f64 / (restarts * total).max(1) as f64;
        let i = Some(self.fitting.load(Ordering::Relaxed)).filter(|&i| i != NONE);
        let fitting = (i, self.epochs_began.get(), self.aim.get());
        let (Some(i), Some(began), Some(aim)) = fitting else {
            let least = (0..restarts).map(|i| self.done(i)).min().unwrap_or(0);
            return Some(format!(
                "k-means starts of {restarts} restarts: {least}/{warmup} iterations; {} in",
                clock(ran)
            ));
        };
        let epochs = total - warmup;
        let epoch = self.epochs_run(i, warmup);
        let mut line = format!("restart {}/{restarts}: {epoch}/{epochs} epochs", i + 1);
        let per = |f: f64| f / aim.units;
        // Its last two ELBOs, newest first.
        let tail: Vec<f64> = self.trace_of(i).iter().rev().take(2).map(|p| p.1).collect();
        if let Some(f) = tail.first() {
            line.push_str(&format!(", ELBO per {} {:.5}", aim.unit, per(*f)));
        }
        match tail[..] {
            _ if epoch < aim.svi => line.push_str(" (stochastic epochs: no test yet)"),
            [b, a] => line.push_str(&format!(
                " (its last epoch moved it {:.1e}; settled under {:.0e})",
                per(b - a).abs(),
                aim.tol
            )),
            _ => {}
        }
        line.push_str(&format!("; {share:.0}% of the steps, {} in", clock(ran)));
        // The pace of the epochs run, not of those credited without running.
        let run: usize = (0..restarts).map(|i| self.epochs_run(i, warmup)).sum();
        let left: usize = (0..restarts)
            .map(|i| total.saturating_sub(self.done(i)))
            .sum();
        let since = Instant::elapsed(began).as_secs_f64();
        if run > 0 && since > 0.0 {
            let eta = Duration::from_secs_f64(since / run as f64 * left as f64);
            line.push_str(&format!(", at most about {} left", clock(eta)));
        }
        Some(line)
    }
}

/// [`Progress::standing`] said on stderr every [`EVERY`] until `stopped`
/// hangs up, each time it changed: printed, not logged, so that it shows
/// without `-v`.
pub(super) fn report(progress: &Progress, stopped: Receiver<()>) {
    let mut said = None;
    while let Err(RecvTimeoutError::Timeout) = stopped.recv_timeout(EVERY) {
        let Some(line) = progress.standing() else {
            continue;
        };
        if said.as_ref() != Some(&line) {
            eprintln!("faba cluster: {line}");
            said = Some(line);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cluster::vi::Aim;

    const AIM: Aim = Aim {
        units: 100.0,
        unit: crate::cluster::vi::Unit::Read,
        tol: 1e-5,
        svi: 2,
    };

    /// Before the fit begins nothing is said; during the starts, their
    /// iterations; once a restart fits, its epochs, ELBO per unit, the last
    /// move against the tolerance, the share done and the time left.
    #[test]
    fn the_standing_follows_the_fit() {
        let p = Progress::new(2);
        assert_eq!(p.standing(), None);
        p.begin_vi(40, 30);
        p.reach(0, 12);
        p.reach(1, 10);
        let line = p.standing().unwrap();
        assert!(
            line.starts_with("k-means starts of 2 restarts: 10/30 iterations"),
            "{line}"
        );
        p.aim(AIM);
        p.reach(0, 30);
        p.reach(1, 30);
        p.fitting(0);
        std::thread::sleep(Duration::from_millis(5));
        p.reach(0, 34);
        p.note(0, 33, -50.0);
        p.note(0, 34, -49.0);
        let line = p.standing().unwrap();
        assert!(
            line.starts_with("restart 1/2: 4/10 epochs, ELBO per read -0.49000"),
            "{line}"
        );
        assert!(
            line.contains("moved it 1.0e-2; settled under 1e-5"),
            "{line}"
        );
        assert!(
            line.contains("80% of the steps") && line.contains("left"),
            "{line}"
        );
        // Steps credited without running (a failed start, a restart that
        // settled) do not count toward the pace.
        p.reach(1, 40);
        assert_eq!((p.epochs_run(0, 30), p.epochs_run(1, 30)), (4, 0));
    }

    /// In the stochastic epochs no test reads the ELBO's moves.
    #[test]
    fn stochastic_epochs_are_not_tested() {
        let p = Progress::new(1);
        p.begin_vi(40, 30);
        p.aim(AIM);
        p.reach(0, 31);
        p.fitting(0);
        p.note(0, 30, -50.0);
        p.note(0, 31, -49.0);
        let line = p.standing().unwrap();
        assert!(line.contains("(stochastic epochs: no test yet)"), "{line}");
    }
}
