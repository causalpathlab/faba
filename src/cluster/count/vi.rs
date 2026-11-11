//! The count model fitted by variational inference: a k-means start, then
//! stochastic variational inference over minibatches of cells, then
//! coordinate ascent to convergence (after hierarchical Poisson
//! factorization, Gopalan, Hofman and Blei 2015, and stochastic variational
//! inference, Hoffman, Blei, Wang and Paisley 2013).
//!
//! The model: each count Poisson about `sum_ab theta_ca lambda_ab psi_gb`,
//! the role weights of a fixed sparse shape `eta` (HPF's), each cell's and
//! gene's own rate taking its depth.
//!
//! ```text
//! theta_ca ~ Gamma(eta, rho_c),   rho_c ~ Gamma(e0, f0)       (local, per cell)
//! psi_gb ~ Gamma(eta, rho'_g),    rho'_g ~ Gamma(e0, f0)      (global)
//! lambda_ab ~ Gamma(kappa, kappa)
//! ```
//!
//! Every factor of `q` is a Gamma, by shape and rate, each rate stored; a
//! nonzero count's units are shared over the blocks by `phi`, which
//! factorizes so that a count costs `A`, not `A B`:
//! `Z_cg = sum_a T_ca H_ga`, `H_ga = sum_b L_ab P_gb`, each capital
//! `exp(E ln)` of its factor. A cell's factors are fitted against the
//! globals (its local step); the globals then take a step toward their
//! targets from a batch's statistics, scaled up to every cell: the whole
//! step in coordinate ascent, a share `rho_t = (tau0 + t)^(-forget)` of it
//! in a stochastic epoch. Only nonzero counts are visited, by cell, and by
//! gene through each batch's own transpose: nothing cells x genes is dense.
//!
//! The ELBO `f` is taken with `phi` at its optimum given the rest; in
//! coordinate ascent it never falls, and an epoch's first passes give it at
//! the epoch's start at no extra cost.

use super::{sum_parts, Counts, Settings};
use crate::cluster::vi::{self, Hooks, Options};
use rand::rngs::SmallRng;
use rayon::prelude::*;
use special::Gamma as _;
use statrs::function::gamma::{digamma, ln_gamma};

/// Local passes of a cell's step at most, and the relative change of its
/// role shapes under which it stops.
const LOCAL_ITERS: usize = 10;
const LOCAL_TOL: f64 = 1e-4;

/// The priors, as the model has them.
#[derive(Debug, Clone, Copy)]
struct Prior {
    eta: f64,
    kappa: f64,
    e0: f64,
    f0: f64,
}

impl Prior {
    /// The shape of `q`'s rate of a cell or gene of `k` roles.
    fn rho_shape(&self, k: usize) -> f64 {
        self.e0 + k as f64 * self.eta
    }
}

/// `E ln x` of a Gamma by shape and rate.
fn e_ln(shape: f64, rate: f64) -> f64 {
    digamma(shape) - rate.ln()
}

/// The entropy of a Gamma by shape and rate.
fn entropy(shape: f64, rate: f64) -> f64 {
    shape - rate.ln() + ln_gamma(shape) + (1.0 - shape) * digamma(shape)
}

/// `E ln p(x) - E ln q(x)` for `x` of `q = Gamma(shape, rate)` under the
/// prior `Gamma(a, b)`, `b` random with mean `eb` and `E ln b` `eln_b`.
fn kl_term(a: f64, eb: f64, eln_b: f64, shape: f64, rate: f64) -> f64 {
    a * eln_b - ln_gamma(a) + (a - 1.0) * e_ln(shape, rate) - eb * shape / rate
        + entropy(shape, rate)
}

/// [`kl_term`] with the prior's rate `b` fixed.
fn kl_fixed(a: f64, b: f64, shape: f64, rate: f64) -> f64 {
    kl_term(a, b, b.ln(), shape, rate)
}

/// Each column's total of the Gamma means `shape / rate`, `k` columns a
/// row.
fn mean_totals(shape: &[f64], rate: &[f64], k: usize) -> Vec<f64> {
    let mut out = vec![0.0; k];
    for (s, r) in shape.chunks_exact(k).zip(rate.chunks_exact(k)) {
        for ((o, s), r) in out.iter_mut().zip(s).zip(r) {
            *o += s / r;
        }
    }
    out
}

/// `F_b = sum_a E lambda_ab Theta_a`: psi's exposure per unit, past its
/// gene's rate.
fn exposure(g: &Globals, theta_total: &[f64]) -> Vec<f64> {
    let bb = g.bb;
    (0..bb)
        .map(|b| {
            (0..g.aa)
                .map(|a| g.l[a * bb + b] / g.lr[a * bb + b] * theta_total[a])
                .sum()
        })
        .collect()
}

/// The global factors of `q`.
#[derive(Debug, Clone)]
pub struct Globals {
    aa: usize,
    bb: usize,
    /// `psi`, genes x B, by shape and rate.
    p: Vec<f64>,
    pq: Vec<f64>,
    /// Each gene's `rho'` rate (its shape `e0 + B eta`).
    qr: Vec<f64>,
    /// `lambda`, A x B.
    l: Vec<f64>,
    lr: Vec<f64>,
}

/// The local factors of `q`, every cell's.
#[derive(Debug, Clone)]
pub struct Locals {
    /// `theta`, cells x A, by shape and rate.
    t: Vec<f64>,
    tr: Vec<f64>,
    /// Each cell's `rho` rate (its shape `e0 + A eta`).
    pr: Vec<f64>,
}

/// What the local steps and the ELBO read off the globals.
struct Cache {
    l_exp: Vec<f64>,
    p_exp: Vec<f64>,
    /// `H_ga = sum_b L_ab P_gb`, genes x A.
    h: Vec<f64>,
    /// `E_a = sum_b E lambda_ab Psi_b`: theta's exposure per unit.
    e_a: Vec<f64>,
}

impl Globals {
    fn genes(&self) -> usize {
        self.qr.len()
    }

    fn cache(&self) -> Cache {
        let (aa, bb, gg) = (self.aa, self.bb, self.genes());
        let e_lambda: Vec<f64> = self.l.iter().zip(&self.lr).map(|(s, r)| s / r).collect();
        let l_exp: Vec<f64> = (self.l.iter().zip(&self.lr))
            .map(|(&s, &r)| e_ln(s, r).exp())
            .collect();
        let p_exp: Vec<f64> = (self.p.par_iter().zip(&self.pq))
            .map(|(&s, &r)| e_ln(s, r).exp())
            .collect();
        let psi_total = mean_totals(&self.p, &self.pq, bb);
        let mut h = vec![0.0; gg * aa];
        h.par_chunks_mut(aa).enumerate().for_each(|(g, row)| {
            let pg = &p_exp[g * bb..(g + 1) * bb];
            for (a, x) in row.iter_mut().enumerate() {
                *x = (l_exp[a * bb..(a + 1) * bb].iter().zip(pg))
                    .map(|(l, p)| l * p)
                    .sum();
            }
        });
        let e_a = (0..aa)
            .map(|a| (0..bb).map(|b| e_lambda[a * bb + b] * psi_total[b]).sum())
            .collect();
        Cache {
            l_exp,
            p_exp,
            h,
            e_a,
        }
    }

    /// The globals' part of the ELBO: each factor's prior less its `q`.
    fn kl(&self, pr: &Prior) -> f64 {
        let bb = self.bb;
        let rho_shape = pr.rho_shape(bb);
        let genes: f64 = (0..self.genes())
            .into_par_iter()
            .map(|g| {
                let (e_rho, eln_rho) = (rho_shape / self.qr[g], e_ln(rho_shape, self.qr[g]));
                let psi: f64 = (g * bb..(g + 1) * bb)
                    .map(|k| kl_term(pr.eta, e_rho, eln_rho, self.p[k], self.pq[k]))
                    .sum();
                psi + kl_fixed(pr.e0, pr.f0, rho_shape, self.qr[g])
            })
            .collect::<Vec<f64>>()
            .iter()
            .sum();
        let blocks: f64 = (self.l.iter().zip(&self.lr))
            .map(|(&s, &r)| kl_fixed(pr.kappa, pr.kappa, s, r))
            .sum();
        genes + blocks
    }
}

/// A cell's locals, by value, for its step.
#[derive(Clone)]
struct Cell {
    t: Vec<f64>,
    tr: Vec<f64>,
    pr: f64,
}

impl Cell {
    /// Cell `c`'s stored locals.
    fn of(locals: &Locals, c: usize, aa: usize) -> Self {
        Self {
            t: locals.t[c * aa..(c + 1) * aa].to_vec(),
            tr: locals.tr[c * aa..(c + 1) * aa].to_vec(),
            pr: locals.pr[c],
        }
    }

    /// Stored back as cell `c`'s.
    fn store(&self, locals: &mut Locals, c: usize, aa: usize) {
        locals.t[c * aa..(c + 1) * aa].copy_from_slice(&self.t);
        locals.tr[c * aa..(c + 1) * aa].copy_from_slice(&self.tr);
        locals.pr[c] = self.pr;
    }

    /// Its `E theta` by role, added into `out`.
    fn means_into(&self, out: &mut [f64]) {
        for ((o, t), r) in out.iter_mut().zip(&self.t).zip(&self.tr) {
            *o += t / r;
        }
    }
}

/// One pass over a cell's counts: `n`, `x / Z` per count (when kept), the
/// `exp(E ln)` of its theta, and `sum x ln Z` (when asked; else 0).
struct Pass {
    n: Vec<f64>,
    w: Vec<f64>,
    t_exp: Vec<f64>,
    data: f64,
}

/// One cell's step: its new locals, its last pass, and its terms of the
/// ELBO at its stored locals (before) and at its new ones (after).
struct Step {
    cell: Cell,
    last: Pass,
    before: f64,
    after: f64,
}

/// The data the fit reads, and its priors.
struct Fit<'a> {
    counts: &'a Counts,
    prior: Prior,
    aa: usize,
    bb: usize,
}

/// The statistics cells give the globals: `m_gb` (genes x B) and `k_ab`,
/// from `u_ga = sum_c w_cg T_ca` (see [`Fit::stats`]).
struct Stats {
    m: Vec<f64>,
    k: Vec<f64>,
}

/// What a sweep over cells gives: their ELBO terms before and after their
/// steps, their `E theta` by role summed before and after, and
/// (when asked) each one's `n_ca` in their order.
struct Swept {
    before: f64,
    after: f64,
    theta_old: Vec<f64>,
    theta_new: Vec<f64>,
    n: Vec<f64>,
}

/// Cells stepped together: their steps are held at once, so this bounds
/// a sweep's memory whatever the cells.
const CHUNK: usize = 4096;

/// The fixed parts a sum over cells or genes is split into (see
/// [`sum_parts`]).
const PARTS: usize = 32;

impl<'a> Fit<'a> {
    fn new(counts: &'a Counts, settings: &Settings) -> Self {
        let (e0, f0) = settings.rate_prior;
        Self {
            counts,
            prior: Prior {
                eta: settings.eta,
                kappa: settings.kappa,
                e0,
                f0,
            },
            aa: settings.cell_roles,
            bb: settings.gene_roles,
        }
    }

    fn rho_shape(&self) -> f64 {
        self.prior.rho_shape(self.aa)
    }

    /// A cell's own prior less its `q`: its theta and rho.
    fn cell_kl(&self, cell: &Cell) -> f64 {
        let pr = &self.prior;
        let shape = self.rho_shape();
        let (e_rho, eln_rho) = (shape / cell.pr, e_ln(shape, cell.pr));
        let theta: f64 = (cell.t.iter().zip(&cell.tr))
            .map(|(&s, &r)| kl_term(pr.eta, e_rho, eln_rho, s, r))
            .sum();
        theta + kl_fixed(pr.e0, pr.f0, shape, cell.pr)
    }

    /// The ELBO's terms that are not any cell's own: the globals' prior less
    /// their `q`, every rate's expectation summed (zeros included:
    /// `sum_a Theta_a E_a`), and `sum ln x!`; `kl` the globals' (see
    /// [`Globals::kl`]).
    fn global_terms(&self, kl: f64, cache: &Cache, theta: &[f64]) -> f64 {
        let blocks: f64 = theta.iter().zip(&cache.e_a).map(|(t, e)| t * e).sum();
        kl - blocks - self.counts.ln_factorials
    }

    /// One pass over cell `c`'s counts at `cell`, `x / Z` kept when `keep_w`,
    /// `sum x ln Z` taken when `data`.
    fn pass(&self, c: usize, cell: &Cell, cache: &Cache, keep_w: bool, data: bool) -> Pass {
        let want = data;
        let aa = self.aa;
        let (genes, xs) = self.counts.of(c);
        let t_exp: Vec<f64> = (cell.t.iter().zip(&cell.tr))
            .map(|(&s, &r)| e_ln(s, r).exp())
            .collect();
        let (mut n, mut data) = (vec![0.0; aa], 0.0);
        let mut w = Vec::with_capacity(if keep_w { genes.len() } else { 0 });
        for (&g, &x) in genes.iter().zip(xs) {
            let g = g as usize;
            let h = &cache.h[g * aa..(g + 1) * aa];
            let z = t_exp.iter().zip(h).map(|(t, h)| t * h).sum::<f64>();
            let x = x as f64;
            let wx = x / z;
            n.iter_mut()
                .zip(&t_exp)
                .zip(h)
                .for_each(|((n, t), h)| *n += wx * t * h);
            if want {
                data += x * z.ln();
            }
            if keep_w {
                w.push(wx);
            }
        }
        Pass { n, w, t_exp, data }
    }

    /// Cell `c`'s step from `cell` against `cache`: up to `iters` rounds of
    /// a pass then theta and rho, then a last pass (`x / Z` kept
    /// when `keep_w`).
    fn step(&self, c: usize, mut cell: Cell, cache: &Cache, iters: usize, keep_w: bool) -> Step {
        let pr = self.prior;
        let shape = self.rho_shape();
        let mut before = None;
        for _ in 0..iters {
            let pass = self.pass(c, &cell, cache, false, before.is_none());
            before.get_or_insert_with(|| pass.data + self.cell_kl(&cell));
            let e_rho = shape / cell.pr;
            let mut change: f64 = 0.0;
            for a in 0..self.aa {
                let t = pr.eta + pass.n[a];
                change = change.max((t - cell.t[a]).abs() / t);
                cell.t[a] = t;
                cell.tr[a] = e_rho + cache.e_a[a];
            }
            cell.pr = pr.f0
                + (cell.t.iter().zip(&cell.tr))
                    .map(|(s, r)| s / r)
                    .sum::<f64>();
            if change < LOCAL_TOL {
                break;
            }
        }
        let last = self.pass(c, &cell, cache, keep_w, true);
        let after = last.data + self.cell_kl(&cell);
        Step {
            before: before.unwrap_or(after),
            after,
            cell,
            last,
        }
    }

    /// The steps of `cells` from their stored locals against `cache`, up to
    /// `iters` rounds each, [`CHUNK`] cells at a time (each chunk's in
    /// parallel), stored back: their `u_ga` added into `u` (genes x A) when
    /// given, their `n_ca` kept when `keep_n`. None when `cancelled` says so
    /// between chunks.
    #[allow(clippy::too_many_arguments)]
    fn sweep(
        &self,
        cells: &[usize],
        locals: &mut Locals,
        cache: &Cache,
        iters: usize,
        mut u: Option<&mut [f64]>,
        keep_n: bool,
        cancelled: &(dyn Fn() -> bool + Sync),
    ) -> Option<Swept> {
        let aa = self.aa;
        let mut out = Swept {
            before: 0.0,
            after: 0.0,
            theta_old: vec![0.0; aa],
            theta_new: vec![0.0; aa],
            n: Vec::with_capacity(if keep_n { cells.len() * aa } else { 0 }),
        };
        for chunk in cells.chunks(CHUNK) {
            if cancelled() {
                return None;
            }
            let (held, keep_w): (&Locals, bool) = (locals, u.is_some());
            let steps: Vec<Step> = chunk
                .par_iter()
                .map(|&c| self.step(c, Cell::of(held, c, aa), cache, iters, keep_w))
                .collect();
            for (&c, s) in chunk.iter().zip(&steps) {
                out.before += s.before;
                out.after += s.after;
                let r = c * aa..(c + 1) * aa;
                (out.theta_old
                    .iter_mut()
                    .zip(&locals.t[r.clone()])
                    .zip(&locals.tr[r]))
                .for_each(|((o, t), r)| *o += t / r);
                s.cell.means_into(&mut out.theta_new);
                if keep_n {
                    out.n.extend_from_slice(&s.last.n);
                }
                s.cell.store(locals, c, aa);
            }
            if let Some(u) = u.as_deref_mut() {
                self.add_u(chunk, &steps, u);
            }
        }
        Some(out)
    }

    /// `u_ga = sum_c w_cg T_ca` of `cells` (`cells[i]`'s step `steps[i]`)
    /// added into `u` (genes x A), in fixed parts (see [`sum_parts`]).
    fn add_u(&self, cells: &[usize], steps: &[Step], u: &mut [f64]) {
        let aa = self.aa;
        let sum = sum_parts(cells.len(), u.len(), PARTS, |acc, i| {
            let last = &steps[i].last;
            for (&g, w) in self.counts.of(cells[i]).0.iter().zip(&last.w) {
                let row = &mut acc[g as usize * aa..(g as usize + 1) * aa];
                row.iter_mut()
                    .zip(&last.t_exp)
                    .for_each(|(u, t)| *u += w * t);
            }
        });
        u.iter_mut().zip(&sum).for_each(|(u, s)| *u += s);
    }

    /// The statistics of `u` (genes x A, see [`Fit::sweep`]) against
    /// `cache`: `m_gb = P_gb sum_a L_ab u_ga`, `k_ab = L_ab sum_g P_gb u_ga`.
    fn stats(&self, u: &[f64], cache: &Cache) -> Stats {
        let (aa, bb) = (self.aa, self.bb);
        // The units of gene `g` in block (a, b).
        let units = |g: usize, a: usize, b: usize| {
            cache.l_exp[a * bb + b] * cache.p_exp[g * bb + b] * u[g * aa + a]
        };
        let mut m = vec![0.0; self.counts.genes * bb];
        m.par_chunks_mut(bb).enumerate().for_each(|(g, row)| {
            for (b, x) in row.iter_mut().enumerate() {
                *x = (0..aa).map(|a| units(g, a, b)).sum();
            }
        });
        let k = sum_parts(self.counts.genes, aa * bb, PARTS, |k, g| {
            for (i, x) in k.iter_mut().enumerate() {
                *x += units(g, i / bb, i % bb);
            }
        });
        Stats { m, k }
    }

    /// A sweep of every cell `all` against `globals`, its statistics and
    /// the `n_ca` kept: the cache it ran against, the globals' KL there, its
    /// sums and statistics. None when cancelled.
    fn full(
        &self,
        all: &[usize],
        locals: &mut Locals,
        globals: &Globals,
        cancelled: &(dyn Fn() -> bool + Sync),
    ) -> Option<(Cache, f64, Swept, Stats)> {
        let cache = globals.cache();
        let mut u = vec![0.0; self.counts.genes * self.aa];
        let swept = self.sweep(
            all,
            locals,
            &cache,
            LOCAL_ITERS,
            Some(&mut u),
            true,
            cancelled,
        )?;
        let stats = self.stats(&u, &cache);
        Some((cache, globals.kl(&self.prior), swept, stats))
    }
}

/// A step of the globals toward the targets of `stats`, scaled by `scale`
/// (every cell over the batch's), by a share `rho` (1: the whole step), the
/// cells' totals `theta_total` exact.
fn global_step(
    g: &mut Globals,
    pr: &Prior,
    stats: &Stats,
    theta_total: &[f64],
    scale: f64,
    rho: f64,
) {
    let bb = g.bb;
    let blend = |x: &mut f64, target: f64| *x = (1.0 - rho) * *x + rho * target;
    let rho_shape = pr.rho_shape(bb);
    // psi, then rho' from it.
    let f = exposure(g, theta_total);
    (g.p.par_chunks_mut(bb).zip(g.pq.par_chunks_mut(bb)))
        .zip(g.qr.par_iter_mut())
        .enumerate()
        .for_each(|(gene, ((p, pq), qr))| {
            let e_rho = rho_shape / *qr;
            for b in 0..bb {
                blend(&mut p[b], pr.eta + scale * stats.m[gene * bb + b]);
                blend(&mut pq[b], e_rho + f[b]);
            }
            *qr = pr.f0 + p.iter().zip(pq.iter()).map(|(s, r)| s / r).sum::<f64>();
        });
    // lambda, with the new psi.
    let psi_total = mean_totals(&g.p, &g.pq, bb);
    for (k, (l, lr)) in g.l.iter_mut().zip(g.lr.iter_mut()).enumerate() {
        blend(l, pr.kappa + scale * stats.k[k]);
        blend(lr, pr.kappa + theta_total[k / bb] * psi_total[k % bb]);
    }
}

/// A fit's end: its globals, its last pass's statistics, its ELBO, and
/// how it ended.
pub struct Fitted {
    globals: Globals,
    /// Each cell's `n_ca` and the cells' statistics in the last pass, at
    /// the factors kept.
    n: Vec<f64>,
    stats: Stats,
    /// The final ELBO, exact.
    pub elbo: f64,
    pub ran: vi::Ran,
}

impl Fitted {
    /// The ELBO per count: the comparison between restarts and the
    /// tolerance's scale.
    pub fn elbo_per_count(&self, counts: &Counts) -> f64 {
        self.elbo / counts.total().max(1.0)
    }
}

/// Fit `counts` from `start` on the shared schedule (see
/// [`crate::cluster::vi`]). None when it was cancelled, or stopped before
/// an epoch ended.
pub fn fit(
    counts: &Counts,
    settings: &Settings,
    options: &Options,
    start: Start,
    rng: &mut SmallRng,
    hooks: &Hooks,
) -> Option<Fitted> {
    let fit = Fit::new(counts, settings);
    let locals = initial_locals(&fit, &start);
    let mut steps = CountSteps {
        all: (0..counts.cells()).collect(),
        theta_total: Vec::new(),
        globals: start.globals,
        locals,
        cancelled: hooks.cancelled,
        fit: &fit,
    };
    let ((elbo, n, stats), ran) = vi::run(&mut steps, counts.cells(), options, rng, hooks)?;
    Some(Fitted {
        globals: steps.globals,
        n,
        stats,
        elbo,
        ran,
    })
}

/// The count fit's steps, its state as they leave it: each role's
/// `Theta_a` summed over the cells as they go (set right once a stochastic
/// epoch).
struct CountSteps<'f, 'h> {
    fit: &'f Fit<'f>,
    all: Vec<usize>,
    locals: Locals,
    globals: Globals,
    theta_total: Vec<f64>,
    cancelled: &'h (dyn Fn() -> bool + Sync),
}

impl vi::Steps for CountSteps<'_, '_> {
    /// The epoch's ELBO at its end, each gene's counts, the statistics.
    type Pass = (f64, Vec<f64>, Stats);
    type End = (f64, Vec<f64>, Stats);

    fn units(&self) -> f64 {
        self.fit.counts.total()
    }

    fn first(&mut self) -> Option<()> {
        let cache = self.globals.cache();
        let swept = self.fit.sweep(
            &self.all,
            &mut self.locals,
            &cache,
            2 * LOCAL_ITERS,
            None,
            false,
            self.cancelled,
        )?;
        self.theta_total = swept.theta_new;
        Some(())
    }

    fn batch(&mut self, cells: &[usize], scale: f64, rho: f64) -> Option<()> {
        let fit = self.fit;
        let cache = self.globals.cache();
        let mut u = vec![0.0; fit.counts.genes * fit.aa];
        let swept = fit.sweep(
            cells,
            &mut self.locals,
            &cache,
            LOCAL_ITERS,
            Some(&mut u),
            false,
            self.cancelled,
        )?;
        for ((t, new), old) in self
            .theta_total
            .iter_mut()
            .zip(&swept.theta_new)
            .zip(&swept.theta_old)
        {
            *t += new - old;
        }
        let stats = fit.stats(&u, &cache);
        let theta = &self.theta_total;
        global_step(&mut self.globals, &fit.prior, &stats, theta, scale, rho);
        Some(())
    }

    fn end_svi_epoch(&mut self) {
        // Drift from the running sums, set right once an epoch.
        self.theta_total = totals(&self.locals, self.fit.aa);
    }

    fn score(&mut self) -> Option<f64> {
        score(
            self.fit,
            &self.all,
            &mut self.locals,
            &self.globals,
            self.cancelled,
        )
    }

    fn epoch(&mut self) -> Option<(f64, Self::Pass)> {
        let fit = self.fit;
        let (cache, kl, swept, stats) =
            fit.full(&self.all, &mut self.locals, &self.globals, self.cancelled)?;
        let f = swept.before + fit.global_terms(kl, &cache, &self.theta_total);
        self.theta_total = swept.theta_new;
        let after = swept.after + fit.global_terms(kl, &cache, &self.theta_total);
        Some((f, (after, swept.n, stats)))
    }

    fn converged(&mut self, pass: Self::Pass) -> Self::End {
        pass
    }

    fn step(&mut self, (_, _, stats): Self::Pass) {
        let prior = &self.fit.prior;
        global_step(
            &mut self.globals,
            prior,
            &stats,
            &self.theta_total,
            1.0,
            1.0,
        );
    }

    fn last(&mut self) -> Option<Self::End> {
        let fit = self.fit;
        let (cache, kl, swept, stats) =
            fit.full(&self.all, &mut self.locals, &self.globals, self.cancelled)?;
        let f = swept.after + fit.global_terms(kl, &cache, &swept.theta_new);
        Some((f, swept.n, stats))
    }

    fn elbo(end: &Self::End) -> f64 {
        end.0
    }
}

/// Every cell's locals before its first step: theta's shapes from the
/// start's units, rates from the globals, rho at mean 1.
fn initial_locals(fit: &Fit, start: &Start) -> Locals {
    let (cc, aa, pr) = (fit.counts.cells(), fit.aa, fit.prior);
    let cache = start.globals.cache();
    let e_rho = 1.0;
    let t = (0..cc)
        .into_par_iter()
        .flat_map_iter(|c| {
            let total: f64 = fit.counts.of(c).1.iter().map(|&x| x as f64).sum();
            (0..aa).map(move |a| pr.eta + start.units(c, a, aa, total))
        })
        .collect();
    Locals {
        t,
        tr: (0..cc * aa).map(|k| e_rho + cache.e_a[k % aa]).collect(),
        pr: vec![fit.rho_shape(); cc],
    }
}

/// `Theta_a`, the cells' `E theta` per role.
fn totals(locals: &Locals, aa: usize) -> Vec<f64> {
    mean_totals(&locals.t, &locals.tr, aa)
}

/// The ELBO at the stored locals and the globals, exact: a pass over every
/// cell, nothing updated (no step taken). None when cancelled.
fn score(
    fit: &Fit,
    all: &[usize],
    locals: &mut Locals,
    globals: &Globals,
    cancelled: &(dyn Fn() -> bool + Sync),
) -> Option<f64> {
    let cache = globals.cache();
    let swept = fit.sweep(all, locals, &cache, 0, None, false, cancelled)?;
    Some(swept.after + fit.global_terms(globals.kl(&fit.prior), &cache, &swept.theta_new))
}

/// What the tables give of a fit.
pub struct Summary {
    /// Each cell's share of its block units by cell role (cells x A), and
    /// each gene's by gene role (genes x B).
    pub cell: Vec<f32>,
    pub gene: Vec<f32>,
    /// Each block's `E ln lambda` and its sd.
    pub log_rate: Vec<f64>,
    pub log_rate_sd: Vec<f64>,
    /// Each block's units.
    pub units: Vec<f64>,
}

/// The tables of a fit, from its last pass.
pub fn summarize(fitted: &Fitted) -> Summary {
    let globals = &fitted.globals;
    let (aa, bb) = (globals.aa, globals.bb);
    let shares = |x: &[f64], k: usize| -> Vec<f32> {
        x.chunks_exact(k)
            .flat_map(|row| {
                let total: f64 = row.iter().sum();
                row.iter().map(move |v| {
                    if total > 0.0 {
                        (v / total) as f32
                    } else {
                        1.0 / k as f32
                    }
                })
            })
            .collect()
    };
    Summary {
        cell: shares(&fitted.n, aa),
        gene: shares(&fitted.stats.m, bb),
        log_rate: (globals.l.iter().zip(&globals.lr))
            .map(|(&s, &r)| e_ln(s, r))
            .collect(),
        log_rate_sd: globals.l.iter().map(|&s| s.trigamma().sqrt()).collect(),
        units: fitted.stats.k.clone(),
    }
}

mod start;
pub use start::{kmeans_start, Profiles, Start, ITERS as KMEANS_ITERS};
pub(in crate::cluster) use start::{soft, spherical, Rows};

#[cfg(test)]
mod tests;
