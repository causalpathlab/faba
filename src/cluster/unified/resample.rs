//! The conditional draws beside the labels: the block shifts and the cells'
//! levels, by slice sampling (the global weights are the cells' side's, in
//! [`super::urn`]).

use super::{Model, State};
use rand::rngs::SmallRng;
use rand::RngExt;

impl Model<'_> {
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
                            .map(|&at| {
                                let s = at / kk;
                                self.marginal(s, self.pseudo(s, d), st.on[at], st.off[at])
                            })
                            .sum::<f64>()
                    } else {
                        let e = (-d).exp();
                        members
                            .iter()
                            .map(|&l| {
                                let c = self.links.cell[l] as usize;
                                self.link_reads(l, (st.nu_exp[c], st.nu[c]), (e, d))
                            })
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
                let e = (-v).exp();
                -v * v / (2.0 * var)
                    + own
                        .clone()
                        .map(|l| {
                            let k = self.block(st, l);
                            self.link_reads(l, (e, v), (st.delta_exp[k], st.delta[k]))
                        })
                        .sum::<f64>()
            };
            let v = slice(st.nu[c], &target, 1.0, rng);
            st.nu[c] = v;
            self.relevel(st, c);
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
