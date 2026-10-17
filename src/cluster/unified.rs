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
//! and the chains summarised together (see [`summary`]).

use super::crp::{ln_beta, ln_pred, ln_sigmoid, logit, pooled_fractions, sample_log};
use super::read::Reads;
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::RngExt;
use resample::slice;
use std::sync::{Arc, Mutex};
use swap::Scratch;

#[cfg(test)]
mod check;
mod draw;
mod links;
mod place;
mod resample;
mod settings;
mod summary;
mod swap;

#[cfg(test)]
pub use check::Held;
use draw::Tally;
pub use draw::{Draw, Means};
pub use links::Links;
pub use settings::{strength, Recorded, Settings, Weights};
pub use summary::summarize;

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
    /// Each site's baseline, the logit of its pooled fraction.
    mu: Vec<f64>,
    /// The genes' urn base, `1 / B` each.
    gene_base: Vec<f64>,
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
    /// Each block's shift.
    delta: Vec<f64>,
    /// Each cell's level (all 0 at finite strength, where there is none).
    nu: Vec<f64>,
    /// Global cell-role weights.
    beta: Vec<f64>,
    /// Each site's profile pseudo-counts in each block (`sites x K`; empty
    /// at infinite strength).
    prior: Vec<(f64, f64)>,
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
        let mu = pooled_fractions(&links.site, &links.on, &links.off, sites)
            .into_iter()
            .map(logit)
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
            mu,
            gene_base: vec![1.0 / settings.gene_roles as f64; settings.gene_roles],
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

    /// Block `k`'s pseudo-counts at every site, from its shift.
    fn refresh(&self, st: &mut State, k: usize) {
        if !self.finite() {
            return;
        }
        let kk = self.blocks();
        for s in 0..self.mu.len() {
            st.prior[s * kk + k] = self.pseudo(s, st.delta[k]);
        }
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
        marginal(st.prior[at], st.on[at], st.off[at])
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
        }
    }

    /// Link `l`'s reads' log predictive in block `k` given the other links
    /// (`l` counted out).
    fn ln_link(&self, st: &State, l: usize, k: usize) -> f64 {
        let kk = self.blocks();
        let (s, m, u) = self.links.of(l);
        if self.finite() {
            (0..s.len())
                .map(|i| {
                    let at = s[i] as usize * kk + k;
                    let (x, y) = st.prior[at];
                    ln_pred(x + st.on[at], y + st.off[at], m[i], u[i])
                })
                .sum()
        } else {
            self.link_reads(l, st.delta[k] + st.nu[self.links.cell[l] as usize])
        }
    }

    /// Link `l`'s reads' Binomial log-likelihood at offset `d` from the
    /// baseline (infinite strength).
    fn link_reads(&self, l: usize, d: f64) -> f64 {
        let (s, m, u) = self.links.of(l);
        (0..s.len())
            .map(|i| self.binomial(s[i] as usize, m[i], u[i], d))
            .sum()
    }

    /// Link `l`'s pair of roles from its conditional, `scratch`'s weights
    /// and gene logs its buffers.
    fn relabel(&self, st: &mut State, l: usize, scratch: &mut Scratch, rng: &mut SmallRng) {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let (c, g) = (self.links.cell[l] as usize, self.link_gene[l] as usize);
        self.put(st, l, false);
        let gb = self.settings.gene_alpha / bb as f64;
        let Scratch { w, gene, .. } = scratch;
        for (b, x) in gene.iter_mut().enumerate() {
            *x = (st.m[g * bb + b] as f64 + gb).ln();
        }
        for a in 0..aa {
            let cell = (st.n[c * aa + a] as f64 + self.settings.alpha * st.beta[a]).ln();
            for b in 0..bb {
                w[a * bb + b] = cell + gene[b] + self.ln_link(st, l, a * bb + b);
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
            nu: vec![0.0; self.links.cells],
            beta: vec![1.0 / aa as f64; aa],
            prior: vec![(0.0, 0.0); per_block],
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

/// A block's Beta-binomial log marginal at a site with profile
/// pseudo-counts `(x, y)` holding reads `(on, off)` there.
fn marginal((x, y): (f64, f64), on: f64, off: f64) -> f64 {
    ln_beta(x + on, y + off) - ln_beta(x, y)
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
