//! The conditional draws beside the labels: the block shifts, slice
//! sampled at finite strength; at infinite strength the shifts and the
//! cells' levels by random-walk Metropolis, and a translation of both
//! along the line where the reads do not change (the global weights are
//! the cells' side's, in [`super::urn`]).

use super::{ln_normal, Model, State};
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::{Distribution, StandardNormal};

/// Random-walk Metropolis steps per shift and per level each sweep.
const WALK_STEPS: usize = 2;

/// The multiples of its estimated scale a walk's step is drawn at, each
/// alike likely (independent of the offset, so the walk stays symmetric).
const WALK_SCALES: [f64; 3] = [0.5, 1.0, 2.0];

/// Which logit offset a walk moves: a block's shift (each link's level
/// held) or a cell's level (each link's block's shift held).
#[derive(Clone, Copy)]
enum Offset {
    Shift,
    Level,
}

impl Model<'_> {
    /// Each block's shift from its conditional: its Normal prior times its
    /// sites' marginals (finite strength, by slice sampling, Neal 2003) or
    /// its links' reads each at its cell's level (infinite, by
    /// [`Self::walk`]). Each block's sites holding reads, or links, are
    /// listed once a sweep, and the target reads them where they are.
    pub(super) fn resample_shifts(&self, st: &mut State, rng: &mut SmallRng) {
        let Some(sd) = self.settings.shift_sd else {
            return;
        };
        if self.finite() {
            self.slice_shifts(st, sd, rng);
        } else {
            self.walk_shifts(st, sd, rng);
        }
    }

    /// [`Self::resample_shifts`] at finite strength.
    fn slice_shifts(&self, st: &mut State, sd: f64, rng: &mut SmallRng) {
        let kk = self.blocks();
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); kk];
        for (at, &h) in st.holders.iter().enumerate() {
            if h > 0 {
                members[at % kk].push(at);
            }
        }
        for (k, members) in members.iter().enumerate() {
            let target = |d: f64| {
                ln_normal(&[d], sd)
                    + members
                        .iter()
                        .map(|&at| {
                            let s = at / kk;
                            self.marginal(s, self.pseudo(s, d), st.on[at], st.off[at])
                        })
                        .sum::<f64>()
            };
            st.delta[k] = slice(st.delta[k], &target, 1.0, rng);
            self.refresh(st, k);
        }
    }

    /// [`Self::resample_shifts`] at infinite strength.
    fn walk_shifts(&self, st: &mut State, sd: f64, rng: &mut SmallRng) {
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); self.blocks()];
        for l in 0..self.links.len() {
            members[self.block(st, l)].push(l);
        }
        for (k, members) in members.iter().enumerate() {
            let level = |l: usize| {
                let c = self.links.cell[l] as usize;
                (st.nu_exp[c], st.nu[c])
            };
            let links = members.iter().copied();
            st.delta[k] = self.walk(st.delta[k], sd, (links, level), Offset::Shift, rng);
            self.refresh(st, k);
        }
    }

    /// Each cell's level from its conditional (its Normal prior times its
    /// links' Binomial reads at their blocks' shifts), by [`Self::walk`];
    /// at infinite strength only, where there are levels.
    pub(super) fn resample_levels(&self, st: &mut State, rng: &mut SmallRng) {
        if self.finite() {
            return;
        }
        for c in 0..self.links.cells {
            let shift = |l: usize| {
                let k = self.block(st, l);
                (st.delta_exp[k], st.delta[k])
            };
            let links = self.links.cell_ptr[c]..self.links.cell_ptr[c + 1];
            st.nu[c] = self.walk(
                st.nu[c],
                self.settings.cell_sd,
                (links, shift),
                Offset::Level,
                rng,
            );
            self.relevel(st, c);
        }
    }

    /// [`WALK_STEPS`] random-walk Metropolis steps for an offset at `x`
    /// under its Normal(0, `sd`) prior and its `links`' reads, each link's
    /// other offset (and its `exp(-x)`) given by `held`. Each step is a
    /// Normal step of 2.4 posterior standard deviations (near the best for
    /// one dimension, Gelman, Roberts and Gilks 1996), those taken from
    /// the links' information with the offset itself left out, so that
    /// the scale does not depend on it and the walk stays symmetric; each
    /// step's scale is that times 1/2, 1 or 2 at random, so that a scale
    /// off by the offset's own reach still moves. The
    /// first step scores the current point and the proposal in one pass;
    /// later ones carry the current point's score.
    fn walk<I: Iterator<Item = usize> + Clone>(
        &self,
        mut x: f64,
        sd: f64,
        (links, held): (I, impl Fn(usize) -> (f64, f64)),
        offset: Offset,
        rng: &mut SmallRng,
    ) -> f64 {
        let info: f64 = links
            .clone()
            .map(|l| self.link_information(l, held(l).0))
            .sum();
        let scale = 2.4 / (info + 1.0 / (sd * sd)).sqrt();
        let score = |points: &[f64; 2], both: bool| -> [f64; 2] {
            let at = points.map(|x| ((-x).exp(), x));
            let pair = |p: (f64, f64), a: (f64, f64)| match offset {
                Offset::Shift => (p, a),
                Offset::Level => (a, p),
            };
            let mut out = points.map(|x| ln_normal(&[x], sd));
            for l in links.clone() {
                let p = held(l);
                if both {
                    let r = self.link_reads_at(l, [pair(p, at[0]), pair(p, at[1])]);
                    out[0] += r[0];
                    out[1] += r[1];
                } else {
                    let (nu, delta) = pair(p, at[1]);
                    out[1] += self.link_reads(l, nu, delta);
                }
            }
            out
        };
        let mut fx = None;
        for _ in 0..WALK_STEPS {
            let z: f64 = StandardNormal.sample(rng);
            let y = x + scale * WALK_SCALES[rng.random_range(0..WALK_SCALES.len())] * z;
            let [f_x, f_y] = match fx {
                None => score(&[x, y], true),
                Some(f) => [f, score(&[x, y], false)[1]],
            };
            fx = Some(f_x);
            if (1.0 - rng.random::<f64>()).ln() < f_y - f_x {
                (x, fx) = (y, Some(f_y));
            }
        }
        x
    }

    /// Every shift moved by `t` and every cell's level by `-t`, `t` from its
    /// conditional (infinite strength with a shift only). Each read's logit,
    /// `mu_s + nu_c + delta_ab`, is unchanged, so `t`'s conditional is the
    /// two Normal priors' along that line, a Normal drawn exactly: a group
    /// move (Liu and Sabatti 2000) along the ridge the one-at-a-time updates
    /// cross slowly, which leaves the posterior invariant. It holds only
    /// while every read's logit carries one level and one shift, each
    /// with a Normal(0, sd) prior: an offset any other way (per modality,
    /// or a shift held at 0) needs its own move.
    pub(super) fn resample_translation(&self, st: &mut State, rng: &mut SmallRng) {
        let (Some(sd), false) = (self.settings.shift_sd, self.finite()) else {
            return;
        };
        let (shift_prec, level_prec) = (1.0 / (sd * sd), 1.0 / self.settings.cell_sd.powi(2));
        let precision = st.delta.len() as f64 * shift_prec + st.nu.len() as f64 * level_prec;
        let mean = (st.nu.iter().sum::<f64>() * level_prec
            - st.delta.iter().sum::<f64>() * shift_prec)
            / precision;
        let z: f64 = StandardNormal.sample(rng);
        let t = mean + z / precision.sqrt();
        for k in 0..st.delta.len() {
            st.delta[k] += t;
            self.refresh(st, k);
        }
        for c in 0..st.nu.len() {
            st.nu[c] -= t;
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
