//! The conditional draws beside the labels: the global weights, the block
//! shifts and the cells' levels, the last two by slice sampling.

use super::{marginal, Model, State, Weights};
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::Distribution;

impl Model<'_> {
    /// The global weights from their conditional (HDP only): each cell's
    /// table counts in each role (an Antoniak draw), then
    /// `beta ~ Dirichlet(gamma / A + tables)`.
    pub(super) fn resample_beta(&self, st: &mut State, rng: &mut SmallRng) {
        let Weights::Hdp(gamma) = self.settings.weights else {
            return;
        };
        let aa = self.settings.cell_roles;
        let mut tables = vec![0.0f64; aa];
        for c in 0..self.links.cells {
            for (a, t) in tables.iter_mut().enumerate() {
                let ab = self.settings.alpha * st.beta[a];
                for j in 0..st.n[c * aa + a] {
                    if rng.random::<f64>() < ab / (ab + j as f64) {
                        *t += 1.0;
                    }
                }
            }
        }
        let draws: Vec<f64> = tables
            .iter()
            .map(|&t| {
                rand_distr::Gamma::new(gamma / aa as f64 + t, 1.0)
                    .unwrap()
                    .sample(rng)
                    .max(1e-300)
            })
            .collect();
        let total: f64 = draws.iter().sum();
        st.beta = draws.iter().map(|d| d / total).collect();
    }

    /// Each block's shift from its conditional (its Normal prior times its
    /// sites' marginals, at infinite strength its links' reads each at its
    /// cell's level), by slice sampling (Neal 2003). Each block's sites
    /// holding reads (finite strength) or links (infinite) are listed once
    /// a sweep, and the target reads them where they are.
    pub(super) fn resample_shifts(&self, st: &mut State, rng: &mut SmallRng) {
        let Some(sd) = self.settings.shift_sd else {
            return;
        };
        let kk = self.blocks();
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); kk];
        if self.finite() {
            for (at, &h) in st.holders.iter().enumerate() {
                if h > 0 {
                    members[at % kk].push(at);
                }
            }
        } else {
            for l in 0..self.links.len() {
                members[self.block(st, l)].push(l);
            }
        }
        for (k, members) in members.iter().enumerate() {
            let target = |d: f64| {
                -d * d / (2.0 * sd * sd)
                    + if self.finite() {
                        members
                            .iter()
                            .map(|&at| marginal(self.pseudo(at / kk, d), st.on[at], st.off[at]))
                            .sum::<f64>()
                    } else {
                        members
                            .iter()
                            .map(|&l| self.link_reads(l, d + st.nu[self.links.cell[l] as usize]))
                            .sum::<f64>()
                    }
            };
            let d = slice(st.delta[k], &target, 1.0, rng);
            st.delta[k] = d;
            self.refresh(st, k);
        }
    }

    /// Each cell's level from its conditional (its Normal prior times its
    /// links' Binomial reads at their blocks' shifts), by slice sampling;
    /// at infinite strength only, where there are levels.
    pub(super) fn resample_levels(&self, st: &mut State, rng: &mut SmallRng) {
        if self.finite() {
            return;
        }
        let var = self.settings.cell_sd * self.settings.cell_sd;
        for c in 0..self.links.cells {
            let own = self.links.cell_ptr[c]..self.links.cell_ptr[c + 1];
            let target = |v: f64| {
                -v * v / (2.0 * var)
                    + own
                        .clone()
                        .map(|l| self.link_reads(l, st.delta[self.block(st, l)] + v))
                        .sum::<f64>()
            };
            let v = slice(st.nu[c], &target, 1.0, rng);
            st.nu[c] = v;
        }
    }
}

/// One slice-sampling update of `x` under log density `f` (Neal 2003):
/// a level under `f(x)`, an interval of width `w` placed at random around
/// `x` and stepped out until both ends fall below the level, then points
/// drawn from it, shrinking it toward `x`, until one is above the level.
pub(super) fn slice(x: f64, f: &dyn Fn(f64) -> f64, w: f64, rng: &mut SmallRng) -> f64 {
    let level = f(x) + (1.0 - rng.random::<f64>()).ln();
    let mut lo = x - w * rng.random::<f64>();
    let mut hi = lo + w;
    while f(lo) > level {
        lo -= w;
    }
    while f(hi) > level {
        hi += w;
    }
    loop {
        let y = lo + (hi - lo) * rng.random::<f64>();
        if f(y) > level {
            return y;
        }
        if y < x {
            lo = y;
        } else {
            hi = y;
        }
    }
}
