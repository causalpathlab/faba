//! Checks against the model's definition: draws of each move with the rest
//! held, the log joint, and the swap's ratio beside the log joint's change.

use super::{Model, Scratch, State};
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;

/// What a check holds the labels at: the global weights, each block's
/// shift, and each cell's level (all 0 at finite strength, where they are
/// not read).
pub struct Held<'h> {
    pub beta: &'h [f64],
    pub delta: &'h [f64],
    pub nu: &'h [f64],
}

impl Model<'_> {
    /// The labels and global weights after each of `draws` sweeps of the
    /// label moves after `warmup`, the shifts held at `delta` and the cells'
    /// levels at `nu`.
    pub fn label_draws(
        &self,
        (delta, nu): (&[f64], &[f64]),
        warmup: usize,
        draws: usize,
        seed: u64,
    ) -> Vec<(Vec<u16>, Vec<u16>, Vec<f64>)> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.init(&mut rng);
        let beta = st.beta.clone();
        self.hold(
            &mut st,
            &Held {
                beta: &beta,
                delta,
                nu,
            },
        );
        let mut scratch = Scratch::new(self);
        (0..warmup + draws)
            .filter_map(|sweep| {
                self.label_moves(&mut st, &mut scratch, &mut rng);
                (sweep >= warmup).then(|| (st.a.clone(), st.b.clone(), st.beta.clone()))
            })
            .collect()
    }

    /// `delta[0]` after each of `draws` shift updates, the labels held at
    /// `a`, `b` and the cells' levels at `nu`.
    pub fn shift_draws(
        &self,
        (a, b): (&[u16], &[u16]),
        nu: &[f64],
        draws: usize,
        seed: u64,
    ) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.state(a.to_vec(), b.to_vec());
        if !self.finite() {
            st.nu = nu.to_vec();
        }
        (0..draws)
            .map(|_| {
                self.resample_shifts(&mut st, &mut rng);
                st.delta[0]
            })
            .collect()
    }

    /// Cell 0's level after each of `draws` level updates, the labels held
    /// at `a`, `b` and the shifts at `delta`.
    pub fn level_draws(
        &self,
        (a, b): (&[u16], &[u16]),
        delta: &[f64],
        draws: usize,
        seed: u64,
    ) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.state(a.to_vec(), b.to_vec());
        st.delta = delta.to_vec();
        for k in 0..self.blocks() {
            self.refresh(&mut st, k);
        }
        (0..draws)
            .map(|_| {
                self.resample_levels(&mut st, &mut rng);
                st.nu[0]
            })
            .collect()
    }

    /// `beta[0]` after each of `draws` updates of the global weights, the
    /// labels held at `a`, `b`.
    pub fn beta_draws(&self, a: &[u16], b: &[u16], draws: usize, seed: u64) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.state(a.to_vec(), b.to_vec());
        (0..draws)
            .map(|_| {
                self.resample_beta(&mut st, &mut rng);
                st.beta[0]
            })
            .collect()
    }

    /// The log joint of labels `a`, `b` at `held`.
    pub fn log_joint_at(&self, labels: (&[u16], &[u16]), held: &Held) -> f64 {
        self.log_joint(&self.held_state(labels, held))
    }

    /// The held-out score of a draw at labels `a`, `b` and `held`.
    pub fn heldout_at(&self, labels: (&[u16], &[u16]), held: &Held) -> (f64, f64) {
        self.draw(&self.held_state(labels, held)).heldout
    }

    /// For labels `a`, `b` at `held`: the swap's computed change in log
    /// joint for roles `r`, `r2` on cell `unit` (`cell` true) or gene
    /// `unit`, and the change the log joint of the swapped labels shows
    /// (both none when nothing would move).
    pub fn swap_check(
        &self,
        (a, b): (&[u16], &[u16]),
        held: &Held,
        unit: usize,
        cell: bool,
        (r, r2): (u16, u16),
    ) -> Option<(f64, f64)> {
        let before = self.log_joint(&self.held_state((a, b), held));
        let mut st = self.held_state((a, b), held);
        let computed = self.exchange(&mut st, (unit, cell), (r, r2), &mut Scratch::new(self))?;
        Some((
            computed,
            self.log_joint(&self.held_state((&st.a, &st.b), held)) - before,
        ))
    }

    /// The labels and the live reads per site and block (converted,
    /// unconverted; `sites x K`) after `sweeps` full sweeps from a random
    /// start.
    pub fn live_reads(&self, sweeps: usize, seed: u64) -> (Vec<u16>, Vec<u16>, Vec<f64>, Vec<f64>) {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.init(&mut rng);
        for _ in 0..sweeps {
            self.sweep(&mut st, &mut rng);
        }
        (st.a, st.b, st.on, st.off)
    }

    /// The state of labels `a`, `b` at `held`.
    fn held_state(&self, (a, b): (&[u16], &[u16]), held: &Held) -> State {
        let mut st = self.state(a.to_vec(), b.to_vec());
        self.hold(&mut st, held);
        st
    }

    /// `st` at `held` (the levels only at infinite strength, where there
    /// are any).
    fn hold(&self, st: &mut State, held: &Held) {
        st.beta = held.beta.to_vec();
        st.delta = held.delta.to_vec();
        if !self.finite() {
            st.nu = held.nu.to_vec();
        }
        for k in 0..self.blocks() {
            self.refresh(st, k);
        }
    }
}
