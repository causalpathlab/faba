//! One blockmodel over (cell, gene) links, with link communities and the
//! cell and gene blockmodel as two settings of it.
//!
//! A link (a cell's reads at one gene's sites) takes a cell role
//! `a ~ theta_c` out of `A` and a gene role `b ~ psi_g` out of `B`; its
//! reads at site `s` are `Binomial(n_s, p_s,ab)`, with
//!
//! ```text
//! beta ~ Dirichlet(gamma / A)                 (or fixed at 1 / A)
//! theta_c ~ Dirichlet(alpha beta),  psi_g ~ Dirichlet(gene_alpha / B)
//! delta_ab ~ N(0, shift_sd^2)                 (or 0)
//! p_s,ab ~ Beta(S sigmoid(mu_s + delta_ab), S sigmoid(-mu_s - delta_ab))
//! ```
//!
//! `mu_s` the logit of the site's pooled fraction and `S` the profile
//! strength. At `S = inf` the profile is the shifted baseline itself, and
//! each cell has its own level `nu_c ~ N(0, cell_sd^2)` (so that roles
//! carry patterns, not levels): `logit p = mu_s + nu_c + delta_ab` and the
//! reads are Binomial. At finite `S`, `p` is integrated out per block, so a
//! cell-specific level cannot enter, and there is none. A block is a pair
//! of roles: its shift is shared by every gene in role `b`. With `theta`,
//! `psi` and `p` integrated out, a link's pair of roles is drawn from
//!
//! ```text
//! P(a, b | rest) ∝ (n_ca + alpha beta_a) (m_gb + gene_alpha / B) Π_s pred(m_s, u_s | block ab at s)
//! ```
//!
//! (counts excluding the link): collapsed Gibbs, one link at a time. Whole
//! cells exchange two cell roles and whole genes two gene roles by
//! Metropolis-Hastings swaps, each its own inverse, the ratio the change in
//! log joint (the unit's urn and the blocks' sites it touches); `beta` is
//! drawn from its conditional through each cell's table counts (Teh et al.
//! 2006); each shift by slice sampling through its block's marginals, and
//! at infinite strength each cell's level through its links' reads. Every
//! move leaves the posterior invariant.
//!
//! At a site some link reaches with fractional reads (batch weights) or
//! more than 16, its blocks' Beta functions are taken through a fast
//! `ln Gamma` (the cluster's `fast` module: f32 below 16, about `1e-3`
//! error per term) in every term a move reads, the log joint's included.
//! The target there is the Beta-binomial through that `ln Gamma`, a small
//! shift from the exact one, and the moves leave it invariant exactly;
//! elsewhere the library's `ln Gamma`.
//!
//! Gene roles need the shift: a site belongs to one gene, so without it a
//! gene's role never reaches its reads, and is refused. So does infinite
//! strength: without the shift every block's reads are Binomial at
//! `sigmoid(mu_s)`, and no role reaches them.
//!
//! `faba cluster --model link` is one gene role, no shift, finite strength
//! and HDP weights (link communities); `--model blocks` gene roles, the
//! shift, infinite strength (with the cells' levels) and symmetric weights
//! (cell and gene modules).
//! Cells outside the fit are placed against a fitted chain (see [`place`]),
//! the chains summarised together (see [`summary`]), and their roles
//! aligned (see [`align`]) and their tallies pooled where they agree (see
//! [`pool`]).

use super::crp::{
    binary_split, few_whole, ln_beta, ln_sigmoid, logit, pooled_fractions, rising_ratio, sample_log,
};
use super::fast;
use super::read::Reads;
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::RngExt;
use resample::slice;
use std::f64::consts::LN_2;
use std::sync::{Arc, Mutex};
use swap::Scratch;

mod align;
#[cfg(test)]
mod check;
mod draw;
mod links;
mod place;
mod pool;
mod resample;
mod settings;
mod summary;
mod swap;

pub use align::Alignment;
#[cfg(test)]
pub use align::{align, assign};
#[cfg(test)]
pub use check::Held;
use draw::Tally;
pub use draw::{Draw, Means};
pub use links::Links;
#[cfg(test)]
pub use pool::pool;
pub use pool::Pooling;
pub use settings::{strength, Recorded, Settings, Weights};
pub use summary::{pooled, summarize};

/// The blockmodel over one set of links, with what the moves share.
pub struct Model<'a> {
    links: &'a Links,
    /// Each site's gene, numbered from 0.
    gene_of: &'a [u32],
    /// Each link's gene.
    link_gene: Vec<u32>,
    /// Gene `g`'s links are `gene_link[gene_ptr[g]..gene_ptr[g + 1]]`.
    gene_ptr: Vec<usize>,
    gene_link: Vec<usize>,
    /// Each site's baseline, the logit of its pooled fraction, and
    /// `exp(-mu_s)`.
    mu: Vec<f64>,
    mu_exp: Vec<f64>,
    /// Each link's unconverted reads, and those times their sites'
    /// baselines summed: the linear part of its Binomial log-likelihood.
    link_off: Vec<(f64, f64)>,
    /// The genes' urn base, `1 / B` each.
    gene_base: Vec<f64>,
    /// Sites some link reaches with reads not few and whole, whose
    /// entries' `ln Gamma` terms are cached (finite strength).
    gamma_site: Vec<bool>,
    settings: Settings,
    /// Each cell's held-out reads, scored at every kept draw.
    heldout: Option<&'a Reads>,
}

/// A chain's labels and what they imply: links per cell and cell role
/// (`cells x A`), per gene and gene role (`genes x B`), and reads per site
/// and block (`sites x K`, `K = A x B` blocks, block `a B + b`) with the
/// links holding them (finite strength only: at infinite strength the
/// reads factor by link and none are kept per block).
pub struct State {
    a: Vec<u16>,
    b: Vec<u16>,
    n: Vec<u32>,
    m: Vec<u32>,
    on: Vec<f64>,
    off: Vec<f64>,
    /// Links holding reads at each site in each block: at none, the
    /// block's reads there are set to exactly 0, so adding and taking
    /// away fractional reads leaves no residue.
    holders: Vec<u32>,
    /// Each block's shift, and `exp(-delta_k)`.
    delta: Vec<f64>,
    delta_exp: Vec<f64>,
    /// Each cell's level (all 0 at finite strength, where there is none),
    /// and `exp(-nu_c)`.
    nu: Vec<f64>,
    nu_exp: Vec<f64>,
    /// Global cell-role weights.
    beta: Vec<f64>,
    /// Each site's profile pseudo-counts in each block (`sites x K`; empty
    /// at infinite strength).
    prior: Vec<(f64, f64)>,
    /// Each entry's `ln B` of its pseudo-counts (`sites x K`; empty at
    /// infinite strength).
    prior_beta: Vec<f64>,
    /// At the model's Gamma sites, each entry's `ln Gamma` of its
    /// posterior pseudo-counts `x + on`, `y + off` and their sum (`sites x
    /// K`, kept current there alone; empty at infinite strength).
    lgamma: Vec<[f64; 3]>,
    tally: Arc<Mutex<Tally>>,
}

impl<'a> Model<'a> {
    /// `gene_of`: each site's gene, numbered from 0. Refuses what
    /// [`Settings::check`] refuses.
    pub fn new(
        links: &'a Links,
        gene_of: &'a [u32],
        sites: usize,
        settings: Settings,
    ) -> anyhow::Result<Self> {
        settings.check()?;
        let mu: Vec<f64> = pooled_fractions(&links.site, &links.on, &links.off, sites)
            .into_iter()
            .map(logit)
            .collect();
        let link_off = (0..links.len())
            .map(|l| {
                let (s, _, u) = links.of(l);
                let total = u.iter().sum();
                (total, (0..s.len()).map(|i| u[i] * mu[s[i] as usize]).sum())
            })
            .collect();
        let link_gene = links.genes_of(gene_of);
        let genes = gene_of.iter().max().map_or(0, |&g| g as usize + 1);
        let mut gene_ptr = vec![0usize; genes + 1];
        for &g in &link_gene {
            gene_ptr[g as usize + 1] += 1;
        }
        for g in 0..genes {
            gene_ptr[g + 1] += gene_ptr[g];
        }
        let mut next = gene_ptr.clone();
        let mut gene_link = vec![0usize; link_gene.len()];
        for (l, &g) in link_gene.iter().enumerate() {
            gene_link[next[g as usize]] = l;
            next[g as usize] += 1;
        }
        Ok(Self {
            links,
            gene_of,
            link_gene,
            gene_ptr,
            gene_link,
            mu_exp: mu.iter().map(|&x| (-x).exp()).collect(),
            mu,
            link_off,
            gene_base: vec![1.0 / settings.gene_roles as f64; settings.gene_roles],
            gamma_site: gamma_sites(links, sites),
            settings,
            heldout: None,
        })
    }

    /// Score `heldout` (each cell's held-out reads) at every kept draw.
    pub fn with_heldout(mut self, heldout: &'a Reads) -> Self {
        self.heldout = Some(heldout);
        self
    }

    /// The sites' baselines.
    pub fn baseline(&self) -> &[f64] {
        &self.mu
    }

    fn genes(&self) -> usize {
        self.gene_ptr.len() - 1
    }

    /// Blocks, `A x B`.
    fn blocks(&self) -> usize {
        self.settings.cell_roles * self.settings.gene_roles
    }

    fn finite(&self) -> bool {
        self.settings.strength.is_finite()
    }

    /// Link `l`'s block.
    fn block(&self, st: &State, l: usize) -> usize {
        st.a[l] as usize * self.settings.gene_roles + st.b[l] as usize
    }

    /// Site `s`'s profile pseudo-counts at shift `d`: `S sigmoid(z)` and
    /// `S sigmoid(-z)` at `z = mu_s + d`, from one `exp(-|z|)`, stably.
    fn pseudo(&self, s: usize, d: f64) -> (f64, f64) {
        let (z, strength) = (self.mu[s] + d, self.settings.strength);
        let e = (-z.abs()).exp();
        let (near, far) = (strength / (1.0 + e), strength * e / (1.0 + e));
        if z >= 0.0 {
            (near, far)
        } else {
            (far, near)
        }
    }

    /// Block `k`'s `exp(-delta_k)` and pseudo-counts at every site, from
    /// its shift.
    fn refresh(&self, st: &mut State, k: usize) {
        st.delta_exp[k] = (-st.delta[k]).exp();
        if !self.finite() {
            return;
        }
        let kk = self.blocks();
        for s in 0..self.mu.len() {
            let at = s * kk + k;
            let (x, y) = self.pseudo(s, st.delta[k]);
            st.prior[at] = (x, y);
            st.prior_beta[at] = self.ln_beta_at(s, x, y);
            self.recache(st, s, at);
        }
    }

    /// Entry `at`'s (at site `s`) `ln Gamma` terms from its pseudo-counts
    /// and reads, at the model's Gamma sites (by [`fast::ln_gamma`]).
    fn recache(&self, st: &mut State, s: usize, at: usize) {
        if !self.gamma_site[s] {
            return;
        }
        let (x, y) = st.prior[at];
        let (x, y) = (x + st.on[at], y + st.off[at]);
        st.lgamma[at] = [fast::ln_gamma(x), fast::ln_gamma(y), fast::ln_gamma(x + y)];
    }

    /// `ln B(x, y)` as site `s` takes it: by [`fast::ln_gamma`] at the
    /// model's Gamma sites, the library's elsewhere.
    fn ln_beta_at(&self, s: usize, x: f64, y: f64) -> f64 {
        if self.gamma_site[s] {
            fast::ln_beta(x, y)
        } else {
            ln_beta(x, y)
        }
    }

    /// A block's Beta-binomial log marginal at site `s` with profile
    /// pseudo-counts `(x, y)` holding reads `(on, off)` there (see
    /// [`Self::ln_beta_at`]).
    fn marginal(&self, s: usize, (x, y): (f64, f64), on: f64, off: f64) -> f64 {
        self.ln_beta_at(s, x + on, y + off) - self.ln_beta_at(s, x, y)
    }

    /// Cell `c`'s `exp(-nu_c)`, from its level.
    fn relevel(&self, st: &mut State, c: usize) {
        st.nu_exp[c] = (-st.nu[c]).exp();
    }

    /// Site `s`'s Binomial log-likelihood of reads `(on, off)` at offset
    /// `d` from its baseline (infinite strength: `d` a block's shift and
    /// the reads' cell's level), written with `ln sigmoid(-z) = ln
    /// sigmoid(z) - z` (binomial coefficients left out, the same for every
    /// labelling).
    fn binomial(&self, s: usize, on: f64, off: f64, d: f64) -> f64 {
        let z = self.mu[s] + d;
        (on + off) * ln_sigmoid(z) - off * z
    }

    /// Entry `at`'s (site and block) term of the log joint at finite
    /// strength: its reads' Beta-binomial marginal at the block's current
    /// shift, its pseudo-counts as cached (none where no link holds reads).
    fn site_term(&self, st: &State, at: usize) -> f64 {
        if st.holders[at] == 0 {
            return 0.0;
        }
        if self.gamma_site[at / self.blocks()] {
            let lg = st.lgamma[at];
            return lg[0] + lg[1] - lg[2] - st.prior_beta[at];
        }
        let (x, y) = st.prior[at];
        ln_beta(x + st.on[at], y + st.off[at]) - st.prior_beta[at]
    }

    /// Link `l`'s labels counted in (`add`) or out: its cell's and gene's
    /// role counts, and at finite strength its reads in its block.
    fn put(&self, st: &mut State, l: usize, add: bool) {
        let (aa, bb, kk) = (
            self.settings.cell_roles,
            self.settings.gene_roles,
            self.blocks(),
        );
        let (c, g) = (self.links.cell[l] as usize, self.link_gene[l] as usize);
        let (n, m) = (
            &mut st.n[c * aa + st.a[l] as usize],
            &mut st.m[g * bb + st.b[l] as usize],
        );
        if add {
            *n += 1;
            *m += 1;
        } else {
            *n -= 1;
            *m -= 1;
        }
        if !self.finite() {
            return;
        }
        let k = self.block(st, l);
        let (s, on, off) = self.links.of(l);
        for i in 0..s.len() {
            let at = s[i] as usize * kk + k;
            if add {
                st.holders[at] += 1;
                st.on[at] += on[i];
                st.off[at] += off[i];
            } else {
                st.holders[at] -= 1;
                if st.holders[at] == 0 {
                    (st.on[at], st.off[at]) = (0.0, 0.0);
                } else {
                    st.on[at] -= on[i];
                    st.off[at] -= off[i];
                }
            }
            self.recache(st, s[i] as usize, at);
        }
    }

    /// Link `l`'s reads' Binomial log-likelihood in block `k` at its cell's
    /// level (infinite strength).
    fn ln_link(&self, st: &State, l: usize, k: usize) -> f64 {
        let c = self.links.cell[l] as usize;
        self.link_reads(l, (st.nu_exp[c], st.nu[c]), (st.delta_exp[k], st.delta[k]))
    }

    /// Link `l`'s reads' log predictive in each block (`out`, one per
    /// block) given the other links (`l` counted out): sites outer, each
    /// site's blocks inner over contiguous entries, the path taken once per
    /// site. Away from the model's Gamma sites (every link there has few
    /// whole reads), each block takes a ratio of rising factorials,
    /// multiplied into its running product in `rising` (kept as a mantissa
    /// and a power of 2, so it never underflows) and one log taken per
    /// block at the end. At a Gamma site each block adds its site term's
    /// change, by [`fast::ln_gamma`] from the cached terms, its arguments
    /// formed as [`Self::put`] and [`Self::recache`] would form them with
    /// the link in, so that the weights are the log joint's.
    fn ln_link_blocks(&self, st: &State, l: usize, out: &mut [f64], rising: &mut [(f64, i64)]) {
        if !self.finite() {
            return self.link_reads_blocks(st, l, out, rising);
        }
        let kk = self.blocks();
        let (s, m, u) = self.links.of(l);
        out.fill(0.0);
        rising.fill((1.0, 0));
        for i in 0..s.len() {
            let row = s[i] as usize * kk..(s[i] as usize + 1) * kk;
            let (on, off) = (&st.on[row.clone()], &st.off[row.clone()]);
            let prior = &st.prior[row.clone()];
            if !self.gamma_site[s[i] as usize] {
                for (k, r) in rising.iter_mut().enumerate() {
                    let (px, py) = prior[k];
                    scale(r, rising_ratio(px + on[k], py + off[k], m[i], u[i]));
                }
            } else {
                let lg = &st.lgamma[row];
                for (k, x) in out.iter_mut().enumerate() {
                    *x += pred_gamma(prior[k], (on[k], off[k]), lg[k], (m[i], u[i]));
                }
            }
        }
        for (x, &r) in out.iter_mut().zip(rising.iter()) {
            *x += ln_split(r);
        }
    }

    /// Link `l`'s reads' Binomial log-likelihood at its cell's level `nu`
    /// and block shift `delta`, each given with its `exp(-x)` (infinite
    /// strength). With `z = mu_s + nu + delta` a site's logit, `n` its reads
    /// and `u` its unconverted ones, each site gives `n ln sigmoid(z) - u z`
    /// with `ln sigmoid(z) = -ln (1 + exp(-mu_s) exp(-nu) exp(-delta))`: few
    /// whole reads as `(1 + ...)^n` multiplied into a running product kept
    /// as a mantissa and a power of 2 (one log at the end), others each
    /// through `ln_1p`, the linear part from the link's constants. The
    /// relabelling's every block ([`Self::link_reads_blocks`]) takes the
    /// same steps in the same order, so the two agree bit for bit.
    fn link_reads(&self, l: usize, nu: (f64, f64), delta: (f64, f64)) -> f64 {
        let (s, m, u) = self.links.of(l);
        let (mut product, mut ones) = ((1.0, 0), 0.0);
        for i in 0..s.len() {
            let e = self.mu_exp[s[i] as usize] * nu.0;
            entry_reads(e * delta.0, m[i], u[i], &mut product, &mut ones);
        }
        let (off, off_mu) = self.link_off[l];
        -ln_split(product) - ones - (off_mu + off * (nu.1 + delta.1))
    }

    /// [`Self::link_reads`] in every block at once (`out`, one per block),
    /// at link `l`'s cell's level and each block's shift: sites outer, each
    /// site's `exp(-mu_s - nu_c)` once and the blocks inner, the running
    /// products in `rising`.
    fn link_reads_blocks(&self, st: &State, l: usize, out: &mut [f64], rising: &mut [(f64, i64)]) {
        let c = self.links.cell[l] as usize;
        let (s, m, u) = self.links.of(l);
        out.fill(0.0);
        rising.fill((1.0, 0));
        for i in 0..s.len() {
            let e = self.mu_exp[s[i] as usize] * st.nu_exp[c];
            for ((x, r), &d) in out.iter_mut().zip(rising.iter_mut()).zip(&st.delta_exp) {
                entry_reads(e * d, m[i], u[i], r, x);
            }
        }
        let (off, off_mu) = self.link_off[l];
        for (k, (x, &r)) in out.iter_mut().zip(rising.iter()).enumerate() {
            *x = -ln_split(r) - *x - (off_mu + off * (st.nu[c] + st.delta[k]));
        }
    }

    /// Link `l`'s pair of roles from its conditional, `scratch`'s weights
    /// and gene logs its buffers.
    fn relabel(&self, st: &mut State, l: usize, scratch: &mut Scratch, rng: &mut SmallRng) {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let (c, g) = (self.links.cell[l] as usize, self.link_gene[l] as usize);
        self.put(st, l, false);
        let gb = self.settings.gene_alpha / bb as f64;
        let Scratch {
            w, rising, gene, ..
        } = scratch;
        for (b, x) in gene.iter_mut().enumerate() {
            *x = (st.m[g * bb + b] as f64 + gb).ln();
        }
        self.ln_link_blocks(st, l, w, rising);
        for a in 0..aa {
            let cell = (st.n[c * aa + a] as f64 + self.settings.alpha * st.beta[a]).ln();
            for b in 0..bb {
                w[a * bb + b] += cell + gene[b];
            }
        }
        let k = sample_log(w, rng);
        st.a[l] = (k / bb) as u16;
        st.b[l] = (k % bb) as u16;
        self.put(st, l, true);
    }

    /// Log joint of the labels, the shifts, the cells' levels and the reads
    /// given the global weights (role weights and profiles integrated out;
    /// the weights' prior not included).
    fn log_joint(&self, st: &State) -> f64 {
        let kk = self.blocks();
        let cells: f64 = (0..self.links.cells)
            .map(|c| self.unit_term(st, c, true))
            .sum();
        let genes: f64 = (0..self.genes())
            .map(|g| self.unit_term(st, g, false))
            .sum();
        let normal =
            |x: &[f64], sd: f64| -> f64 { x.iter().map(|d| -d * d / (2.0 * sd * sd)).sum() };
        let shifts = self
            .settings
            .shift_sd
            .map_or(0.0, |sd| normal(&st.delta, sd));
        if !self.finite() {
            let reads: f64 = (0..self.links.len())
                .map(|l| self.ln_link(st, l, self.block(st, l)))
                .sum();
            return cells + genes + reads + shifts + normal(&st.nu, self.settings.cell_sd);
        }
        let sites: f64 = (0..self.mu.len() * kk)
            .map(|at| self.site_term(st, at))
            .sum();
        cells + genes + sites + shifts
    }

    /// Relabel every link once in random order, then a swap per cell and
    /// per gene, then the global weights.
    fn label_moves(&self, st: &mut State, scratch: &mut Scratch, rng: &mut SmallRng) {
        let mut order: Vec<usize> = (0..self.links.len()).collect();
        order.shuffle(rng);
        for l in order {
            self.relabel(st, l, scratch, rng);
        }
        for c in 0..self.links.cells {
            self.swap(st, (c, true), scratch, rng);
        }
        for g in 0..self.genes() {
            self.swap(st, (g, false), scratch, rng);
        }
        self.resample_beta(st, rng);
    }

    /// The counts implied by labels `a`, `b`, at shift 0 and uniform global
    /// weights.
    fn state(&self, a: Vec<u16>, b: Vec<u16>) -> State {
        let (aa, bb, kk) = (
            self.settings.cell_roles,
            self.settings.gene_roles,
            self.blocks(),
        );
        let sites = self.mu.len();
        // Reads per site and block are kept at finite strength alone.
        let per_block = if self.finite() { sites * kk } else { 0 };
        let mut st = State {
            a,
            b,
            n: vec![0; self.links.cells * aa],
            m: vec![0; self.genes() * bb],
            on: vec![0.0; per_block],
            off: vec![0.0; per_block],
            holders: vec![0; per_block],
            delta: vec![0.0; kk],
            delta_exp: vec![1.0; kk],
            nu: vec![0.0; self.links.cells],
            nu_exp: vec![1.0; self.links.cells],
            beta: vec![1.0 / aa as f64; aa],
            prior: vec![(0.0, 0.0); per_block],
            prior_beta: vec![0.0; per_block],
            lgamma: vec![[0.0; 3]; per_block],
            tally: Arc::default(),
        };
        for k in 0..kk {
            self.refresh(&mut st, k);
        }
        for l in 0..self.links.len() {
            self.put(&mut st, l, true);
        }
        st
    }
}

/// One link entry's `m` converted and `u` unconverted reads at a site of
/// `exp(-z) = e` (see [`Model::link_reads`]): few whole reads multiply
/// `(1 + e)^n` into `product` (a mantissa and a power of 2), others add
/// `n ln (1 + e)` to `ones`, as do few whole ones from [`E_OVERFLOWS`] on
/// (or at a NaN), where the power would not be finite.
fn entry_reads(e: f64, m: f64, u: f64, product: &mut (f64, i64), ones: &mut f64) {
    if few_whole(m, u) && e < E_OVERFLOWS {
        scale(product, (1.0 + e).powi((m + u) as i32));
    } else {
        *ones += (m + u) * e.ln_1p();
    }
}

/// A running product kept as a mantissa and a power of 2, multiplied by
/// `x` and split again, so that it never underflows.
fn scale(product: &mut (f64, i64), x: f64) {
    let (mantissa, e) = binary_split(product.0 * x);
    *product = (mantissa, product.1 + e);
}

/// The log of a running product kept as a mantissa and a power of 2.
fn ln_split((mantissa, e): (f64, i64)) -> f64 {
    mantissa.ln() + e as f64 * LN_2
}

/// The `exp(-z)` from which [`entry_reads`] takes few whole reads by
/// `ln_1p`: below it `(1 + e)^16` times a mantissa under 2 is finite.
const E_OVERFLOWS: f64 = (1u64 << 60) as f64;

/// Each of `sites` sites: whether some link reaches it with reads not few
/// and whole.
fn gamma_sites(links: &Links, sites: usize) -> Vec<bool> {
    let mut out = vec![false; sites];
    for (i, &s) in links.site.iter().enumerate() {
        out[s as usize] |= !few_whole(links.on[i], links.off[i]);
    }
    out
}

/// At a Gamma site with pseudo-counts `(x, y)` holding reads `(on, off)`,
/// `lg` their cached `ln Gamma` terms: the change in the site's term (by
/// [`fast::ln_gamma`]) as a link's `m` converted and `u` unconverted reads
/// join, each side's skipped when it adds none.
fn pred_gamma((x, y): (f64, f64), (on, off): (f64, f64), lg: [f64; 3], (m, u): (f64, f64)) -> f64 {
    let (x2, y2) = (x + (on + m), y + (off + u));
    let side = |z: f64, lz: f64, n: f64| if n > 0.0 { fast::ln_gamma(z) - lz } else { 0.0 };
    side(x2, lg[0], m) + side(y2, lg[1], u) - (fast::ln_gamma(x2 + y2) - lg[2])
}

/// The log fractions `(ln p, ln (1 - p))` of the posterior-mean profile at
/// a site with pseudo-counts `(x, y)` holding reads `(on, off)` there.
fn ln_mean((x, y): (f64, f64), on: f64, off: f64) -> (f64, f64) {
    let (on, off) = (x + on, y + off);
    let total = (on + off).ln();
    (on.ln() - total, off.ln() - total)
}

/// The range of gene role `b`'s shifts across the cell roles, from block
/// shifts `shift` (`cell_roles x gene_roles`, row-major).
pub fn shift_range(shift: &[f64], cell_roles: usize, gene_roles: usize, b: usize) -> f64 {
    let s = (0..cell_roles).map(|a| shift[a * gene_roles + b]);
    s.clone().fold(f64::MIN, f64::max) - s.fold(f64::MAX, f64::min)
}

impl McmcModel for Model<'_> {
    type State = State;
    type Sample = Draw;
    type Result = Vec<Draw>;

    /// Each cell's links in one cell role, each gene's in one gene role,
    /// drawn at random.
    fn init(&self, rng: &mut SmallRng) -> State {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let cell: Vec<u16> = (0..self.links.cells)
            .map(|_| rng.random_range(0..aa) as u16)
            .collect();
        let gene: Vec<u16> = (0..self.genes())
            .map(|_| rng.random_range(0..bb) as u16)
            .collect();
        let a = self.links.cell.iter().map(|&c| cell[c as usize]).collect();
        let b = self.link_gene.iter().map(|&g| gene[g as usize]).collect();
        self.state(a, b)
    }

    /// The label moves, then the shifts, then the cells' levels.
    fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        self.label_moves(st, &mut Scratch::new(self), rng);
        self.resample_shifts(st, rng);
        self.resample_levels(st, rng);
    }

    fn collect(&self, st: &State) -> Draw {
        self.draw(st)
    }

    fn summarize(&self, draws: Vec<Draw>) -> Vec<Draw> {
        draws
    }
}

impl crate::cluster::progress::Traced for Model<'_> {
    /// The kept draws' log joint, at any sweep.
    fn log_joint_of(&self, st: &State) -> Option<f64> {
        Some(self.log_joint(st))
    }
}
