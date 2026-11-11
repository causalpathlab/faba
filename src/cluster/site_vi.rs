//! The site blockmodel (`--model blocks`) fitted by variational inference:
//! a k-means start, then stochastic variational inference over minibatches
//! of cells, then coordinate ascent to convergence, as the count fit runs
//! (see [`super::count::vi`]).
//!
//! A link (a cell's reads at one gene's sites) takes a cell role `a ~
//! theta_c` and a gene role `b ~ psi_g`; its reads at site `s` are
//! `Binomial(n, sigmoid(eta))`, `eta = mu_s + nu_c + delta_ab`:
//!
//! ```text
//! theta_c ~ Dir(alpha / A),  psi_g ~ Dir(gene_alpha / B)
//! delta_ab ~ N(0, shift_sd^2),  nu_c ~ N(0, cell_sd^2)
//! ```
//!
//! with `mu_s` the logit of the site's pooled fraction. The likelihood is
//! bounded by Polya-Gamma (Jaakkola and Jordan), one `xi` per read, site and
//! block at its optimum `xi^2 = E eta^2`, where a read's bound is
//! `kappa E eta + n (ln sigmoid(xi) - xi / 2)` (`kappa = m - n / 2`); fixed
//! at those `xi`, it is quadratic in `nu` and `delta`, with `w = n tanh(xi /
//! 2) / (2 xi)` (`E omega`). `q` is mean field: each link's pair of roles a
//! table over the blocks (`phi`, never stored: each cell's step makes its
//! links' afresh), Dirichlets for `theta`, `psi`, Gaussians for `nu`,
//! `delta`. After each global step the shifts and the levels move together
//! by the constant that best fits their priors, which leaves every read's
//! bound as it is (`nu` and `delta` trade a constant freely otherwise).

use super::count::sum_parts;
use super::count::vi::{Hooks, Options};
use super::links::Links;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use statrs::function::gamma::{digamma, ln_gamma};

/// Local passes of a cell's step at most, and the change of its role
/// shapes and level under which it stops.
const LOCAL_ITERS: usize = 10;
const LOCAL_TOL: f64 = 1e-3;

/// Cells stepped together, their steps held at once.
const CHUNK: usize = 4096;

/// The fixed parts a sum over cells is split into (see [`sum_parts`]).
const PARTS: usize = 64;

/// The model's sizes and priors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Prior {
    pub cell_roles: usize,
    pub gene_roles: usize,
    /// Each cell's role weights' concentration over `A` (`alpha`), and each
    /// gene's over `B`.
    pub alpha: f64,
    pub gene_alpha: f64,
    pub shift_sd: f64,
    pub cell_sd: f64,
}

/// `E omega` of `n` reads at `xi`: `n tanh(xi / 2) / (2 xi)`, `n / 4` at 0.
fn omega(n: f64, xi: f64) -> f64 {
    jj(n, xi).0
}

/// `E omega` of `n` reads at `xi` (`n tanh(xi / 2) / (2 xi)`, `n / 4` at 0)
/// and their bound's term
/// free of `eta`, `n (ln sigmoid(xi) - xi / 2)`, from one exponential.
fn jj(n: f64, xi: f64) -> (f64, f64) {
    let q = (-xi).exp();
    let bound = -n * (q.ln_1p() + xi / 2.0);
    let w = match xi < 1e-4 {
        true => n * (0.25 - xi * xi / 48.0),
        false => n * (1.0 - q) / ((1.0 + q) * 2.0 * xi),
    };
    (w, bound)
}

/// `-KL(N(mean, var) || N(0, prior))`, `prior` a variance.
fn gaussian_term(mean: f64, var: f64, prior: f64) -> f64 {
    -0.5 * ((var + mean * mean) / prior - 1.0 - (var / prior).ln())
}

/// `E ln p(theta) - E ln q(theta)` for `q = Dir(t)` under `Dir(a0)`, `t.len()`
/// roles, and `E ln theta` of `q`.
fn dirichlet_term(t: &[f64], a0: f64) -> (f64, Vec<f64>) {
    let total: f64 = t.iter().sum();
    let dg = digamma(total);
    let eln: Vec<f64> = t.iter().map(|&x| digamma(x) - dg).collect();
    let k = t.len() as f64;
    let mut f = ln_gamma(k * a0) - k * ln_gamma(a0) - ln_gamma(total);
    for (x, e) in t.iter().zip(&eln) {
        f += ln_gamma(*x) + (a0 - x) * e;
    }
    (f, eln)
}

/// The reads the fit reads, by link, and what is fixed of them.
pub struct Data<'a> {
    pub links: &'a Links,
    /// Each link's gene.
    pub link_gene: Vec<u32>,
    pub genes: usize,
    /// Each site's baseline: the logit of its pooled fraction.
    pub mu: Vec<f64>,
    /// `sum ln C(n, m)` over the reads (weighted reads by `ln Gamma`): the
    /// likelihood's constant.
    ln_choose: f64,
    /// Every read, weighted.
    reads: f64,
}

impl<'a> Data<'a> {
    pub fn new(links: &'a Links, gene_of: &[u32], genes: usize, mu: Vec<f64>) -> Self {
        let (ln_choose, reads) = (0..links.len())
            .into_par_iter()
            .map(|l| {
                let (_, m, u) = links.of(l);
                (m.iter().zip(u)).fold((0.0, 0.0), |(c, r), (&m, &u)| {
                    let lc = ln_gamma(m + u + 1.0) - ln_gamma(m + 1.0) - ln_gamma(u + 1.0);
                    (c + lc, r + m + u)
                })
            })
            .collect::<Vec<_>>()
            .iter()
            .fold((0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        Self {
            links,
            link_gene: links.genes_of(gene_of),
            genes,
            mu,
            ln_choose,
            reads,
        }
    }

    pub fn cells(&self) -> usize {
        self.links.cells
    }

    /// Every read, weighted: the tolerance's scale.
    pub fn reads(&self) -> f64 {
        self.reads
    }

    /// Cell `c`'s links.
    fn links_of(&self, c: usize) -> std::ops::Range<usize> {
        self.links.cell_ptr[c]..self.links.cell_ptr[c + 1]
    }
}

/// The global factors of `q`: each gene's role weights' Dirichlet (genes x
/// B), each block's shift's mean and variance.
#[derive(Debug, Clone)]
pub struct Globals {
    pub(super) p: Vec<f64>,
    pub(super) d: Vec<f64>,
    pub(super) v: Vec<f64>,
}

/// The local factors of `q`, every cell's: its role weights' Dirichlet
/// (cells x A), its level's mean and variance.
#[derive(Debug, Clone)]
pub struct Locals {
    pub(super) t: Vec<f64>,
    pub(super) e: Vec<f64>,
    pub(super) s: Vec<f64>,
}

/// What a cell's step reads off the globals: each gene's `E ln psi`.
struct Cache {
    eln_psi: Vec<f64>,
}

/// One cell's locals, by value.
#[derive(Clone)]
struct Cell {
    t: Vec<f64>,
    e: f64,
    s: f64,
}

/// A pass's working space, kept across a cell's passes: one link's sums
/// over its sites per block (the bound, `sum w`, `sum (kappa - w mu)`), its
/// `phi`, and each of the cell's (entry, block) bounds fixed at the `xi` of
/// the step's first pass (`w`, `xi^2`, the term free of `eta`).
#[derive(Default)]
struct Scratch {
    j: Vec<f64>,
    w: Vec<f64>,
    r: Vec<f64>,
    phi: Vec<f64>,
    fixed: Vec<[f64; 3]>,
}

/// Where a pass takes each bound's `xi`: at its optimum, kept for the
/// passes after (`Fresh`); or as the step's first pass left it (`Fixed`),
/// the bound then quadratic, still a bound.
#[derive(Clone, Copy, PartialEq)]
enum Xi {
    Fresh,
    Fixed,
}

/// The fit: the data and the priors.
pub struct Fit<'a> {
    pub data: &'a Data<'a>,
    pub prior: Prior,
}

/// What a sweep over cells gives: their ELBO terms at their stored locals
/// (`phi` at its optimum) and after their steps.
struct Swept {
    before: f64,
    after: f64,
}

/// The statistics cells give the globals: each gene's links by gene role
/// (genes x B), and each block's `sum phi w` and `sum phi (kappa - w (mu +
/// nu))`.
struct Stats {
    gene: Vec<f64>,
    prec: Vec<f64>,
    num: Vec<f64>,
}

impl Fit<'_> {
    fn blocks(&self) -> usize {
        self.prior.cell_roles * self.prior.gene_roles
    }

    /// Link `l`'s sums into `sc` at cell level `(e, s)` against the shifts
    /// `(d, v)`, its entries from `at` in the cell's: each bound at its
    /// optimum `xi`, kept, when `xi` is `Fresh`; else at the kept `xi`.
    fn sums(&self, l: usize, (e, s): (f64, f64), g: &Globals, xi: Xi, at: usize, sc: &mut Scratch) {
        let k = self.blocks();
        for v in [&mut sc.j, &mut sc.w, &mut sc.r] {
            v.clear();
            v.resize(k, 0.0);
        }
        let (sites, m, u) = self.data.links.of(l);
        for (i, ((&site, &m), &u)) in sites.iter().zip(m).zip(u).enumerate() {
            let (n, mu) = (m + u, self.data.mu[site as usize]);
            let kappa = m - n / 2.0;
            let fixed = &mut sc.fixed[(at + i) * k..(at + i + 1) * k];
            for (b, keep) in fixed.iter_mut().enumerate() {
                let x = mu + e + g.d[b];
                let square = x * x + s + g.v[b];
                let (w, j) = match xi {
                    Xi::Fresh => {
                        let xi = square.sqrt();
                        let (w, bound) = jj(n, xi);
                        *keep = [w, square, bound];
                        (w, kappa * x + bound)
                    }
                    Xi::Fixed => {
                        let [w, xi2, bound] = *keep;
                        (w, kappa * x - w / 2.0 * (square - xi2) + bound)
                    }
                };
                sc.j[b] += j;
                sc.w[b] += w;
                sc.r[b] += kappa - w * mu;
            }
        }
    }

    /// Link `l`'s `phi` into `sc` (one weight per block) from its sums,
    /// `E ln theta` of its cell and `E ln psi` of its gene; and `ln` of its
    /// normalizer, its term of the ELBO at that `phi`.
    fn phi(&self, l: usize, eln_theta: &[f64], cache: &Cache, sc: &mut Scratch) -> f64 {
        let bb = self.prior.gene_roles;
        let g = self.data.link_gene[l] as usize;
        let eln_psi = &cache.eln_psi[g * bb..(g + 1) * bb];
        sc.phi.clear();
        sc.phi
            .extend((0..self.blocks()).map(|k| eln_theta[k / bb] + eln_psi[k % bb] + sc.j[k]));
        let top = sc.phi.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        sc.phi.iter_mut().for_each(|x| *x = (*x - top).exp());
        let z: f64 = sc.phi.iter().sum();
        sc.phi.iter_mut().for_each(|p| *p /= z);
        top + z.ln()
    }

    /// One pass over cell `c`'s links at `cell`: its terms of the ELBO
    /// (`phi` at its optimum), its links by cell role, and its level's
    /// precision and precision times mean; with `acc`, the cell's statistics
    /// for the globals added into it (genes x B, then each block's `sum phi
    /// w`, then its `sum phi (kappa - w (mu + nu))`).
    #[allow(clippy::too_many_arguments)]
    fn pass(
        &self,
        c: usize,
        cell: &Cell,
        g: &Globals,
        cache: &Cache,
        mut acc: Option<&mut [f64]>,
        xi: Xi,
        sc: &mut Scratch,
    ) -> (f64, Vec<f64>, f64, f64) {
        let pr = self.prior;
        let (aa, bb, k) = (pr.cell_roles, pr.gene_roles, self.blocks());
        let genes = self.data.genes * bb;
        let cell_var = pr.cell_sd * pr.cell_sd;
        let (dir, eln) = dirichlet_term(&cell.t, pr.alpha / aa as f64);
        let mut f = dir + gaussian_term(cell.e, cell.s, cell_var);
        let mut n = vec![0.0; aa];
        let (mut prec, mut num) = (1.0 / cell_var, 0.0);
        let links = self.data.links_of(c);
        let first = self.data.links.ptr[links.start];
        if xi == Xi::Fresh {
            let entries = self.data.links.ptr[links.end] - first;
            sc.fixed.resize(entries * k, [0.0; 3]);
        }
        for l in links {
            self.sums(
                l,
                (cell.e, cell.s),
                g,
                xi,
                self.data.links.ptr[l] - first,
                sc,
            );
            f += self.phi(l, &eln, cache, sc);
            for (b, p) in sc.phi.iter().enumerate() {
                n[b / bb] += p;
                prec += p * sc.w[b];
                num += p * (sc.r[b] - g.d[b] * sc.w[b]);
            }
            if let Some(acc) = acc.as_deref_mut() {
                let gene = self.data.link_gene[l] as usize;
                for (b, p) in sc.phi.iter().enumerate() {
                    acc[gene * bb + b % bb] += p;
                    acc[genes + b] += p * sc.w[b];
                    acc[genes + k + b] += p * (sc.r[b] - cell.e * sc.w[b]);
                }
            }
        }
        (f, n, prec, num)
    }

    /// Cell `c`'s step from `cell`: up to `iters` rounds of its links'
    /// `phi`, then its role weights and its level, every bound's `xi` kept
    /// from the first round (each round still raises the bound, and the ELBO
    /// with it); its terms of the ELBO at `cell` as given.
    fn step(
        &self,
        c: usize,
        mut cell: Cell,
        g: &Globals,
        cache: &Cache,
        iters: usize,
    ) -> (Cell, f64) {
        let a0 = self.prior.alpha / self.prior.cell_roles as f64;
        let mut before = None;
        let mut sc = Scratch::default();
        for round in 0..iters.max(1) {
            // The first pass's bounds at their optimum (the ELBO at `cell`),
            // kept for the passes after.
            let xi = if round == 0 { Xi::Fresh } else { Xi::Fixed };
            let (f, n, prec, num) = self.pass(c, &cell, g, cache, None, xi, &mut sc);
            before.get_or_insert(f);
            if iters == 0 {
                break;
            }
            let t: Vec<f64> = n.iter().map(|x| a0 + x).collect();
            let (e, s) = (num / prec, 1.0 / prec);
            let change = (t.iter().zip(&cell.t))
                .map(|(x, y)| (x - y).abs() / x)
                .fold((e - cell.e).abs(), f64::max);
            cell = Cell { t, e, s };
            if change < LOCAL_TOL {
                break;
            }
        }
        (cell, before.expect("one pass"))
    }

    /// The steps of `cells` from their stored locals, up to `iters` rounds
    /// each (none: a score), [`CHUNK`] cells at a time (each chunk's in
    /// parallel), stored back; then, when `after` or `stats` asks, a pass at
    /// the new locals for their terms of the ELBO, their statistics added
    /// into `stats` when given. None when `cancelled` says so between
    /// chunks.
    #[allow(clippy::too_many_arguments)]
    fn sweep(
        &self,
        cells: &[usize],
        locals: &mut Locals,
        g: &Globals,
        iters: usize,
        mut stats: Option<&mut Stats>,
        after: bool,
        cancelled: &(dyn Fn() -> bool + Sync),
    ) -> Option<Swept> {
        let (aa, k) = (self.prior.cell_roles, self.blocks());
        let genes = self.data.genes * self.prior.gene_roles;
        let cache = self.cache(g);
        let mut out = Swept {
            before: 0.0,
            after: 0.0,
        };
        let cell_of = |l: &Locals, c: usize| Cell {
            t: l.t[c * aa..(c + 1) * aa].to_vec(),
            e: l.e[c],
            s: l.s[c],
        };
        for chunk in cells.chunks(CHUNK) {
            if cancelled() {
                return None;
            }
            let held: &Locals = locals;
            let steps: Vec<(Cell, f64)> = (chunk.par_iter())
                .map(|&c| self.step(c, cell_of(held, c), g, &cache, iters))
                .collect();
            for (&c, (cell, before)) in chunk.iter().zip(steps) {
                out.before += before;
                locals.t[c * aa..(c + 1) * aa].copy_from_slice(&cell.t);
                (locals.e[c], locals.s[c]) = (cell.e, cell.s);
            }
            if !(after || stats.is_some()) {
                continue;
            }
            // The cells' terms at their new locals, their statistics with
            // them: in fixed parts, the last slot the terms.
            let held: &Locals = locals;
            let want = stats.is_some();
            let len = if want { genes + 2 * k + 1 } else { 1 };
            let sum = sum_parts(chunk.len(), len, PARTS, |acc, i| {
                let (last, rest) = acc.split_last_mut().expect("a slot for the terms");
                let into = want.then_some(rest);
                let mut sc = Scratch::default();
                let cell = cell_of(held, chunk[i]);
                *last += self
                    .pass(chunk[i], &cell, g, &cache, into, Xi::Fresh, &mut sc)
                    .0;
            });
            out.after += sum[len - 1];
            if let Some(stats) = stats.as_deref_mut() {
                stats
                    .gene
                    .iter_mut()
                    .zip(&sum[..genes])
                    .for_each(|(o, x)| *o += x);
                (stats.prec.iter_mut())
                    .zip(&sum[genes..genes + k])
                    .for_each(|(o, x)| *o += x);
                (stats.num.iter_mut())
                    .zip(&sum[genes + k..genes + 2 * k])
                    .for_each(|(o, x)| *o += x);
            }
        }
        if iters == 0 {
            out.after = out.before;
        }
        Some(out)
    }

    fn stats(&self) -> Stats {
        let k = self.blocks();
        Stats {
            gene: vec![0.0; self.data.genes * self.prior.gene_roles],
            prec: vec![0.0; k],
            num: vec![0.0; k],
        }
    }

    fn cache(&self, g: &Globals) -> Cache {
        let bb = self.prior.gene_roles;
        let eln_psi = (g.p.par_chunks(bb))
            .flat_map_iter(|row| {
                let dg = digamma(row.iter().sum());
                row.iter()
                    .map(move |&x| digamma(x) - dg)
                    .collect::<Vec<_>>()
            })
            .collect();
        Cache { eln_psi }
    }

    /// The ELBO's terms that are not any cell's own: the genes' role weights'
    /// and the shifts' prior less their `q`, and `sum ln C(n, m)`.
    fn global_terms(&self, g: &Globals) -> f64 {
        let pr = self.prior;
        let b0 = pr.gene_alpha / pr.gene_roles as f64;
        let genes: f64 = (g.p.chunks_exact(pr.gene_roles))
            .map(|row| dirichlet_term(row, b0).0)
            .sum();
        let var = pr.shift_sd * pr.shift_sd;
        let shifts: f64 = (g.d.iter().zip(&g.v))
            .map(|(&d, &v)| gaussian_term(d, v, var))
            .sum();
        genes + shifts + self.data.ln_choose
    }

    /// A step of the globals toward the targets of `stats`, scaled by
    /// `scale` (every cell over the batch's), by a share `rho` of it (1: the
    /// whole step); the shifts' natural parameters blended.
    fn global_step(&self, g: &mut Globals, stats: &Stats, scale: f64, rho: f64) {
        let pr = self.prior;
        let b0 = pr.gene_alpha / pr.gene_roles as f64;
        for (p, x) in g.p.iter_mut().zip(&stats.gene) {
            *p = (1.0 - rho) * *p + rho * (b0 + scale * x);
        }
        let precision = 1.0 / (pr.shift_sd * pr.shift_sd);
        for k in 0..self.blocks() {
            let (n1, n2) = (g.d[k] / g.v[k], 1.0 / g.v[k]);
            let (t1, t2) = (scale * stats.num[k], precision + scale * stats.prec[k]);
            let (n1, n2) = ((1.0 - rho) * n1 + rho * t1, (1.0 - rho) * n2 + rho * t2);
            (g.v[k], g.d[k]) = (1.0 / n2, n1 / n2);
        }
    }

    /// The shifts and every cell's level moved together, `delta + t`, `nu -
    /// t`, by the `t` their priors favour: each read's bound unchanged.
    fn translate(&self, g: &mut Globals, locals: &mut Locals) {
        let pr = self.prior;
        let (cv, sv) = (pr.cell_sd * pr.cell_sd, pr.shift_sd * pr.shift_sd);
        let e: f64 = locals.e.iter().sum();
        let d: f64 = g.d.iter().sum();
        let t = (e / cv - d / sv) / (locals.e.len() as f64 / cv + g.d.len() as f64 / sv);
        g.d.iter_mut().for_each(|x| *x += t);
        locals.e.iter_mut().for_each(|x| *x -= t);
    }
}

/// A fit's end: its factors, its ELBO and how it ended.
pub struct Fitted {
    pub globals: Globals,
    pub locals: Locals,
    pub elbo: f64,
    pub svi_epochs: usize,
    pub epochs: usize,
    pub converged: bool,
}

impl Fitted {
    /// The ELBO per read: the comparison between restarts.
    pub fn elbo_per_read(&self, data: &Data) -> f64 {
        self.elbo / data.reads().max(1.0)
    }
}

/// The ELBO at the stored locals, nothing changed. None when cancelled.
fn score(
    fit: &Fit,
    all: &[usize],
    locals: &mut Locals,
    g: &Globals,
    cancelled: &(dyn Fn() -> bool + Sync),
) -> Option<f64> {
    let swept = fit.sweep(all, locals, g, 0, None, false, cancelled)?;
    Some(swept.after + fit.global_terms(g))
}

/// Fit `data` from `start` (see [`fit`] in the count model for the
/// schedule): a local step for every cell, stochastic epochs when a batch is
/// fewer than every cell, then coordinate ascent to the tolerance, at most
/// `options.epochs` in all. None when cancelled, or stopped before an
/// epoch ended.
pub fn fit(
    fit: &Fit,
    options: &Options,
    start: Start,
    rng: &mut SmallRng,
    hooks: &Hooks,
) -> Option<Fitted> {
    let cc = fit.data.cells();
    let cancelled = hooks.cancelled;
    let (mut globals, mut locals) = (start.globals, start.locals);
    let all: Vec<usize> = (0..cc).collect();
    fit.sweep(
        &all,
        &mut locals,
        &globals,
        2 * LOCAL_ITERS,
        None,
        false,
        cancelled,
    )?;
    let batches = cc.div_ceil(options.batch_size.max(1));
    let svi = if batches > 1 {
        options.svi_epochs.min(options.epochs)
    } else {
        0
    };
    let batch = cc.div_ceil(batches.max(1));
    let tau0 = batches as f64;
    let (mut t, mut done) = (0usize, 0);
    let mut order = all.clone();
    'svi: for e in 1..=svi {
        order.shuffle(rng);
        for cells in order.chunks(batch) {
            if (hooks.stop)() {
                break 'svi;
            }
            let mut stats = fit.stats();
            fit.sweep(
                cells,
                &mut locals,
                &globals,
                LOCAL_ITERS,
                Some(&mut stats),
                false,
                cancelled,
            )?;
            let rho = (tau0 + t as f64).powf(-options.forget);
            fit.global_step(&mut globals, &stats, cc as f64 / cells.len() as f64, rho);
            fit.translate(&mut globals, &mut locals);
            t += 1;
        }
        // The last stochastic epoch is scored by the next local steps.
        if e < svi {
            (hooks.epoch)(e, score(fit, &all, &mut locals, &globals, cancelled)?);
        }
        done = e;
    }
    // Coordinate ascent: each epoch's local steps score the state at its
    // start as they go; converged, the epoch ends there.
    let (mut last, mut converged, mut end) = (None::<f64>, false, None);
    for e in (svi + 1)..=options.epochs {
        if (hooks.stop)() || done < e - 1 {
            break;
        }
        let mut stats = fit.stats();
        let swept = fit.sweep(
            &all,
            &mut locals,
            &globals,
            LOCAL_ITERS,
            Some(&mut stats),
            false,
            cancelled,
        )?;
        let terms = fit.global_terms(&globals);
        let f = swept.before + terms;
        (hooks.epoch)(e - 1, f);
        if last.is_some_and(|prev| (f - prev).abs() / fit.data.reads().max(1.0) < options.tol) {
            converged = true;
            end = Some(swept.after + terms);
            break;
        }
        last = Some(f);
        fit.global_step(&mut globals, &stats, 1.0, 1.0);
        fit.translate(&mut globals, &mut locals);
        done = e;
    }
    if done == 0 && end.is_none() {
        return None;
    }
    let elbo = match end {
        Some(f) => f,
        None => {
            let swept = fit.sweep(
                &all,
                &mut locals,
                &globals,
                LOCAL_ITERS,
                None,
                true,
                cancelled,
            )?;
            swept.after + fit.global_terms(&globals)
        }
    };
    (hooks.epoch)(done, elbo);
    Some(Fitted {
        globals,
        locals,
        elbo,
        svi_epochs: svi.min(done),
        epochs: done,
        converged,
    })
}

mod start;
mod tables;
pub use start::{kmeans_start, Profiles, Start};

#[cfg(test)]
mod tests;
