//! Checks against the model's definition: draws of each move with the rest
//! held, the log joint, and the swap's ratio beside the log joint's change.

use super::draw::{add_columns, Tally};
use super::{ln_normal, Alignment, Blocks, Chain, Draw, Model, Scratch, State};
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use std::sync::{Arc, Mutex};

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
        let mut ch = self.init(&mut rng);
        let beta = ch.urn.beta.clone();
        self.hold(
            &mut ch,
            &Held {
                beta: &beta,
                delta,
                nu,
            },
        );
        let mut scratch = Scratch::new(self);
        (0..warmup + draws)
            .filter_map(|sweep| {
                self.label_moves(&mut ch, &mut scratch, &mut rng);
                let Chain { urn, st } = &ch;
                (sweep >= warmup).then(|| (st.a.clone(), st.b.clone(), urn.beta.clone()))
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
        let mut st = self.chain(a.to_vec(), b.to_vec()).st;
        if !self.finite() {
            st.nu = nu.to_vec();
            (0..self.links.cells).for_each(|c| self.relevel(&mut st, c));
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
        let mut st = self.chain(a.to_vec(), b.to_vec()).st;
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

    /// The log joint before and after one translation of the shifts and
    /// levels (see [`Model::resample_translation`]), the labels held at
    /// `a`, `b`, with the shifts' and levels' Normal priors' change: what
    /// the move changes, read off the log joint and off its priors alone.
    pub fn translation_check(
        &self,
        (a, b): (&[u16], &[u16]),
        (delta, nu): (&[f64], &[f64]),
        seed: u64,
    ) -> (f64, f64) {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let Chain { urn, mut st } = self.chain(a.to_vec(), b.to_vec());
        st.delta = delta.to_vec();
        (0..self.blocks()).for_each(|k| self.refresh(&mut st, k));
        st.nu = nu.to_vec();
        (0..self.links.cells).for_each(|c| self.relevel(&mut st, c));
        let (sd, cell_sd) = (self.settings.shift_sd.unwrap(), self.settings.cell_sd);
        let prior = |st: &State| ln_normal(&st.delta, sd) + ln_normal(&st.nu, cell_sd);
        let (before, prior_before) = (self.log_joint(&urn, &st), prior(&st));
        self.resample_translation(&mut st, &mut rng);
        (
            self.log_joint(&urn, &st) - before,
            prior(&st) - prior_before,
        )
    }

    /// The translation `t` of each of `draws` translations (see
    /// [`Model::resample_translation`]), each from the same shifts `delta`
    /// and levels `nu`, the labels held at `a`, `b`.
    pub fn translation_draws(
        &self,
        (a, b): (&[u16], &[u16]),
        (delta, nu): (&[f64], &[f64]),
        draws: usize,
        seed: u64,
    ) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.chain(a.to_vec(), b.to_vec()).st;
        (0..draws)
            .map(|_| {
                st.delta = delta.to_vec();
                st.nu = nu.to_vec();
                self.resample_translation(&mut st, &mut rng);
                st.delta[0] - delta[0]
            })
            .collect()
    }

    /// `beta[0]` after each of `draws` updates of the global weights, the
    /// labels held at `a`, `b`.
    pub fn beta_draws(&self, a: &[u16], b: &[u16], draws: usize, seed: u64) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut urn = self.chain(a.to_vec(), b.to_vec()).urn;
        (0..draws)
            .map(|_| {
                self.cells.resample(&mut urn, &mut rng);
                urn.beta[0]
            })
            .collect()
    }

    /// The log joint of labels `a`, `b` at `held`.
    pub fn log_joint_at(&self, labels: (&[u16], &[u16]), held: &Held) -> f64 {
        let ch = self.held_state(labels, held);
        self.log_joint(&ch.urn, &ch.st)
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
        let joint = |ch: &Chain| self.log_joint(&ch.urn, &ch.st);
        let before = joint(&self.held_state((a, b), held));
        let Chain { mut urn, mut st } = self.held_state((a, b), held);
        let scratch = &mut Scratch::new(self);
        let computed = self.exchange(&mut urn, &mut st, (unit, cell), (r, r2), scratch)?;
        Some((
            computed,
            joint(&self.held_state((&st.a, &st.b), held)) - before,
        ))
    }

    /// The labels and the live reads per site and block (converted,
    /// unconverted; `sites x K`) after `sweeps` full sweeps from a random
    /// start.
    pub fn live_reads(&self, sweeps: usize, seed: u64) -> (Vec<u16>, Vec<u16>, Vec<f64>, Vec<f64>) {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut ch = self.init(&mut rng);
        for _ in 0..sweeps {
            self.sweep(&mut ch, &mut rng);
        }
        let st = ch.st;
        (st.a, st.b, st.on, st.off)
    }

    /// After each of `sweeps` full sweeps from a random start, a
    /// fingerprint of the chain: a hash of the labels, the log joint's
    /// bits, and a hash of the shifts' bits and of the levels'.
    pub fn fingerprints(&self, sweeps: usize, seed: u64) -> Vec<[u64; 4]> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut ch = self.init(&mut rng);
        (0..sweeps)
            .map(|_| {
                self.sweep(&mut ch, &mut rng);
                let Chain { urn, st } = &ch;
                let labels = st.a.iter().chain(&st.b).map(|&x| u64::from(x));
                let bits = |x: &[f64]| fnv(x.iter().map(|v| v.to_bits()));
                [
                    fnv(labels),
                    self.log_joint(urn, st).to_bits(),
                    bits(&st.delta),
                    bits(&st.nu),
                ]
            })
            .collect()
    }

    /// After `sweeps` full sweeps from a random start and then the label
    /// moves alone (the shifts' refresh would hide a stale entry), the
    /// cached entries checked and those not, bit for bit, what their
    /// pseudo-counts and live reads give (finite strength): each entry's
    /// `ln B` of its pseudo-counts, and at the Gamma sites its `ln Gamma`
    /// terms (both there by the fast `ln Gamma`).
    pub fn stale_caches(&self, sweeps: usize, seed: u64) -> (usize, usize) {
        use crate::cluster::fast;
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut ch = self.init(&mut rng);
        for _ in 0..sweeps {
            self.sweep(&mut ch, &mut rng);
        }
        self.label_moves(&mut ch, &mut Scratch::new(self), &mut rng);
        let st = &ch.st;
        let (mut checked, mut stale) = (0, 0);
        for (at, &(x, y)) in st.prior.iter().enumerate() {
            let s = at / self.blocks();
            let prior = self.ln_beta_at(s, x, y);
            let mut fresh = st.prior_beta[at].to_bits() == prior.to_bits();
            if self.gamma_site[s] {
                let (x, y) = (x + st.on[at], y + st.off[at]);
                let lg = [fast::ln_gamma(x), fast::ln_gamma(y), fast::ln_gamma(x + y)];
                fresh &= st.lgamma[at].map(f64::to_bits) == lg.map(f64::to_bits);
            }
            checked += 1;
            stale += usize::from(!fresh);
        }
        (checked, stale)
    }

    /// For labels `a`, `b` at `held`, link `l` counted out: its reads' log
    /// predictive in each block as the relabelling scores it, and as the
    /// sum over its sites of each one's predictive by the definition (the
    /// Beta-binomial through `ln B`, by the fast `ln Gamma` at the Gamma
    /// sites; at infinite strength the Binomial).
    pub fn link_check(
        &self,
        labels: (&[u16], &[u16]),
        held: &Held,
        l: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        use crate::cluster::crp::ln_sigmoid;
        let Chain { mut urn, mut st } = self.held_state(labels, held);
        self.put(&mut urn, &mut st, l, false);
        let kk = self.blocks();
        let mut scored = vec![0.0; kk];
        self.ln_link_blocks(&st, l, Blocks::ALL, &mut scored, &mut vec![(1.0, 0); kk]);
        let (s, m, u) = self.links.of(l);
        let nu = st.nu[self.links.cell[l] as usize];
        let summed = (0..kk)
            .map(|k| {
                (0..s.len())
                    .map(|i| {
                        let j = s[i] as usize;
                        if self.finite() {
                            let (x, y) = st.prior[j * kk + k];
                            let held = (x + st.on[j * kk + k], y + st.off[j * kk + k]);
                            self.marginal(j, held, m[i], u[i])
                        } else {
                            let z = self.mu[j] + st.delta[k] + nu;
                            m[i] * ln_sigmoid(z) + u[i] * ln_sigmoid(-z)
                        }
                    })
                    .sum()
            })
            .collect();
        (scored, summed)
    }

    /// After each of `sweeps` full sweeps from a random start, its moves
    /// one at a time (the cached exponentials checked after each), the
    /// exponentials checked and those not, bit for bit, `exp(-x)` of the
    /// current shifts and levels.
    pub fn stale_exponentials(&self, sweeps: usize, seed: u64) -> (usize, usize) {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut ch = self.init(&mut rng);
        let mut scratch = Scratch::new(self);
        let (mut checked, mut stale) = (0, 0);
        let mut check = |st: &State| {
            let pairs = st.delta.iter().zip(&st.delta_exp);
            for (x, e) in pairs.chain(st.nu.iter().zip(&st.nu_exp)) {
                checked += 1;
                stale += usize::from((-x).exp().to_bits() != e.to_bits());
            }
        };
        for _ in 0..sweeps {
            self.label_moves(&mut ch, &mut scratch, &mut rng);
            self.resample_shifts(&mut ch.st, &mut rng);
            check(&ch.st);
            self.resample_levels(&mut ch.st, &mut rng);
            check(&ch.st);
            self.resample_translation(&mut ch.st, &mut rng);
            check(&ch.st);
        }
        (checked, stale)
    }

    /// The state of labels `a`, `b` at `held`.
    fn held_state(&self, (a, b): (&[u16], &[u16]), held: &Held) -> Chain {
        let mut ch = self.chain(a.to_vec(), b.to_vec());
        self.hold(&mut ch, held);
        ch
    }

    /// `ch` at `held` (the levels only at infinite strength, where there
    /// are any).
    fn hold(&self, ch: &mut Chain, held: &Held) {
        ch.urn.beta = held.beta.to_vec();
        let st = &mut ch.st;
        st.delta = held.delta.to_vec();
        if !self.finite() {
            st.nu = held.nu.to_vec();
            (0..self.links.cells).for_each(|c| self.relevel(st, c));
        }
        for k in 0..self.blocks() {
            self.refresh(st, k);
        }
    }
}

impl Draw {
    /// `chain` with its roles renumbered by `to` (its cell role `a` to
    /// `to.cell[a]`, gene role `b` to `to.gene[b]`; roles numbered alike
    /// merged, their values added): a tally of its own, each draw's shifts
    /// and cell roles' means moved alike.
    pub fn relabelled(chain: &[Draw], to: &Alignment) -> Vec<Draw> {
        let mut tally = Tally::default();
        if let Some(d) = chain.last() {
            tally.add_aligned(&d.tally.lock().expect("a chain's tally"), to);
        }
        let tally = Arc::new(Mutex::new(tally));
        let blocks = to.blocks();
        let moved = |x: &[f64], to: &[usize]| {
            let mut out = Vec::new();
            add_columns(&mut out, x, to);
            out
        };
        chain
            .iter()
            .map(|d| Draw {
                heldout: d.heldout,
                log_joint: d.log_joint,
                shift: moved(&d.shift, &blocks),
                role_mean: moved(&d.role_mean, &to.cell),
                tally: Arc::clone(&tally),
            })
            .collect()
    }
}

/// The FNV-1a hash of `words`, stable across builds and platforms (what
/// it hashes need not be).
fn fnv(words: impl Iterator<Item = u64>) -> u64 {
    words.fold(0xcbf2_9ce4_8422_2325, |h, w| {
        (h ^ w).wrapping_mul(0x0000_0100_0000_01b3)
    })
}
