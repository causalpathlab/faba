//! The conditional draws beside the labels: the block shifts, slice
//! sampled at finite strength; at infinite strength the shifts and the
//! cells' levels by random-walk Metropolis, a translation of both along
//! the line where the reads do not change, and per cell role a translation
//! of its shifts against its cells' levels (the global weights are the
//! sides', in [`super::urn`]).

use super::{ln_normal, Model, State, Urn};
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::{Distribution, StandardNormal};
use std::ops::Range;

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
            let y = x + walk_step(scale, rng);
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
        let (shifts, cells) = (0..st.delta.len(), 0..st.nu.len());
        let (precision, mean) = self.line_prior(st, (shifts.clone(), cells.clone()), sd);
        let z: f64 = StandardNormal.sample(rng);
        self.translate(st, (shifts, cells), mean + z / precision.sqrt());
    }

    /// The precision and mean of a translation `t` of `shifts` by `t` and
    /// `cells`' levels by `-t` under the shifts' Normal(0, `sd`) and the
    /// levels' Normal(0, `cell_sd`) priors, a Normal in `t`.
    pub(super) fn line_prior(
        &self,
        st: &State,
        (shifts, cells): (Range<usize>, impl Iterator<Item = usize> + Clone),
        sd: f64,
    ) -> (f64, f64) {
        let (shift_prec, level_prec) = (1.0 / (sd * sd), 1.0 / self.settings.cell_sd.powi(2));
        let precision =
            shifts.len() as f64 * shift_prec + cells.clone().count() as f64 * level_prec;
        let levels: f64 = cells.map(|c| st.nu[c]).sum();
        let mean =
            (levels * level_prec - st.delta[shifts].iter().sum::<f64>() * shift_prec) / precision;
        (precision, mean)
    }

    /// `shifts` moved by `t` and `cells`' levels by `-t`.
    pub(super) fn translate(
        &self,
        st: &mut State,
        (shifts, cells): (Range<usize>, impl Iterator<Item = usize>),
        t: f64,
    ) {
        for k in shifts {
            st.delta[k] += t;
            self.refresh(st, k);
        }
        for c in cells {
            st.nu[c] -= t;
            self.relevel(st, c);
        }
    }

    /// For each cell role `a` in turn, `a`'s shifts moved by `t` and the
    /// levels of `a`'s main cells (see [`RoleLines`]) by `-t`, by a
    /// random-walk Metropolis step along that line (infinite strength with
    /// a shift and two cell roles or more). Only links outside their cell's
    /// main role change their logit, so on cells of mostly one role the
    /// line is nearly flat and the step long: it crosses the ridge between
    /// a role's shifts and its cells' levels, which the global translation
    /// does not reach. The main roles depend on the labels alone, which the
    /// move leaves, so it is its own reverse at `-t`; the target along the
    /// line is the shifts' and levels' Normal priors, a Normal in `t`, times
    /// the moved links' reads. The step's scale is set as in
    /// [`Self::walk`], from the moved links' information (each at the rest
    /// of its logit, the offset the move changes left out: a moved level's
    /// link at its shift, a moved shift's at its level) and the line's
    /// prior precision, neither of which the move changes.
    pub(super) fn resample_role_translations(&self, urn: &Urn, st: &mut State, rng: &mut SmallRng) {
        let (Some(sd), false, 2..) = (self.settings.shift_sd, self.finite(), self.cells.roles)
        else {
            return;
        };
        let lines = self.role_lines(urn, st);
        for a in 0..self.cells.roles {
            self.role_step(st, &lines, a, sd, rng);
        }
    }

    /// Cell role `a`'s step of [`Self::resample_role_translations`].
    pub(super) fn role_step(
        &self,
        st: &mut State,
        lines: &RoleLines,
        a: usize,
        sd: f64,
        rng: &mut SmallRng,
    ) {
        let (precision, mean) = self.line_prior(st, lines.line(a), sd);
        let info = self.line_information(st, lines, a);
        let t = walk_step(2.4 / (info + precision).sqrt(), rng);
        let ratio = -0.5 * precision * ((t - mean).powi(2) - mean * mean)
            + self.line_reads(st, lines, a, t);
        if (1.0 - rng.random::<f64>()).ln() < ratio {
            self.translate(st, lines.line(a), t);
        }
    }

    /// The moved links' information about cell role `a`'s translation,
    /// each at the rest of its logit (a moved level's link at its shift, a
    /// moved shift's at its level): the move leaves it, so the step's
    /// scale is the same either way.
    pub(super) fn line_information(&self, st: &State, lines: &RoleLines, a: usize) -> f64 {
        lines
            .moved(a)
            .map(|(l, level)| {
                let rest = if level {
                    st.delta_exp[self.block(st, l)]
                } else {
                    st.nu_exp[self.links.cell[l] as usize]
                };
                self.link_information(l, rest)
            })
            .sum()
    }

    /// The labels' [`RoleLines`].
    pub(super) fn role_lines(&self, urn: &Urn, st: &State) -> RoleLines {
        let aa = self.cells.roles;
        let main: Vec<usize> = urn
            .n
            .chunks_exact(aa)
            .map(|row| (0..aa).fold(0, |best, a| if row[a] > row[best] { a } else { best }))
            .collect();
        let mut lines = RoleLines {
            bb: self.settings.gene_roles,
            cells: vec![Vec::new(); aa],
            leaving: vec![Vec::new(); aa],
            entering: vec![Vec::new(); aa],
        };
        for (c, &a) in main.iter().enumerate() {
            lines.cells[a].push(c);
        }
        for l in 0..self.links.len() {
            let (home, a) = (main[self.links.cell[l] as usize], st.a[l] as usize);
            if a != home {
                lines.leaving[home].push(l);
                lines.entering[a].push(l);
            }
        }
        lines
    }

    /// The change in the reads of the links that cell role `a`'s
    /// translation by `t` moves (see [`Self::resample_role_translations`]).
    pub(super) fn line_reads(&self, st: &State, lines: &RoleLines, a: usize, t: f64) -> f64 {
        let point = |nu: f64, delta: f64| (((-nu).exp(), nu), ((-delta).exp(), delta));
        lines
            .moved(a)
            .map(|(l, level)| {
                let (c, k) = (self.links.cell[l] as usize, self.block(st, l));
                let (nu, delta) = (st.nu[c], st.delta[k]);
                let now = ((st.nu_exp[c], nu), (st.delta_exp[k], delta));
                let then = if level {
                    point(nu - t, delta)
                } else {
                    point(nu, delta + t)
                };
                let [now, then] = self.link_reads_at(l, [now, then]);
                then - now
            })
            .sum()
    }
}

/// Per cell role, its main cells (those whose links it holds most of, the
/// lowest role on a tie), their links out of it (whose levels a role
/// translation moves) and its links of other cells (whose shifts it
/// moves).
pub(super) struct RoleLines {
    bb: usize,
    cells: Vec<Vec<usize>>,
    leaving: Vec<Vec<usize>>,
    entering: Vec<Vec<usize>>,
}

impl RoleLines {
    /// The links role `a`'s translation moves, each with whether its level
    /// moves (else its shift).
    fn moved(&self, a: usize) -> impl Iterator<Item = (usize, bool)> + '_ {
        let leaving = self.leaving[a].iter().map(|&l| (l, true));
        leaving.chain(self.entering[a].iter().map(|&l| (l, false)))
    }

    /// Role `a`'s line: its shifts and its main cells.
    pub(super) fn line(
        &self,
        a: usize,
    ) -> (Range<usize>, impl Iterator<Item = usize> + Clone + '_) {
        (
            a * self.bb..(a + 1) * self.bb,
            self.cells[a].iter().copied(),
        )
    }
}

/// A random-walk step at `scale` times 1/2, 1 or 2, alike likely (see
/// [`Model::walk`]).
fn walk_step(scale: f64, rng: &mut SmallRng) -> f64 {
    let z: f64 = StandardNormal.sample(rng);
    scale * WALK_SCALES[rng.random_range(0..WALK_SCALES.len())] * z
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
