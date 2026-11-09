//! A Poisson blockmodel of gene counts with a background: cell roles and
//! gene roles over a cell x gene count matrix, after Ball, Karrer and
//! Newman (2011) and the Bayesian Poisson Tucker decomposition (Schein et
//! al. 2016).
//!
//! ```text
//! x_cg ~ Poisson(s_c r_g + sum_ab theta_ca lambda_ab psi_gb)
//! theta_ca ~ Gamma(w_a, rho_c),   psi_gb ~ Gamma(v_b, rho'_g),   lambda_ab ~ Gamma(kappa, kappa)
//! s_c, r_g ~ Gamma(a0, b0),       rho_c, rho'_g ~ Gamma(e0, f0)
//! w_a ~ Gamma(gamma / A, c0),     v_b ~ Gamma(gamma / B, c0)
//! ```
//!
//! (Gamma by shape and rate.) The background `s_c r_g`, a cell's size times
//! a gene's level, belongs to no role: it takes the counts that no block
//! needs, so a block is a departure from it and no role has to hold the
//! bulk of the data (the resolution limit). Each count unit (UMI) is in the
//! background or in one block `(a, b)`; a cell's and a gene's units spread
//! over roles, so cells and genes have mixed memberships. The global
//! weights `w`, `v` are a Gamma-negative-binomial process's (Zhou and
//! Carin 2015): roles not needed empty.
//!
//! Each sweep, uncollapsed Gibbs, every draw from its exact conditional:
//! each nonzero count's units shared out over the background and the
//! blocks, multinomially by their rates (none kept between sweeps); then
//! the global weights through each unit's table counts (a CRT draw, its
//! role's weights integrated out); the cells' weights, their rates, the
//! genes' likewise, the blocks' rates and the background. A zero count
//! enters only through the exposures, `sum_cg` of each rate, which factor
//! into column totals (`sum_ab lambda_ab (sum_c theta_ca)(sum_g psi_gb)`
//! and `(sum_c s_c)(sum_g r_g)`): no zero is visited.

use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::{Binomial, Distribution, Gamma};

/// A cell x gene count matrix, its nonzero entries by cell: cell `c`'s are
/// `ptr[c]..ptr[c + 1]` of `gene` and `x`.
#[derive(Debug, Clone, Default)]
pub struct Counts {
    pub genes: usize,
    pub ptr: Vec<usize>,
    pub gene: Vec<u32>,
    pub x: Vec<u32>,
}

impl Counts {
    pub fn cells(&self) -> usize {
        self.ptr.len().saturating_sub(1)
    }

    /// Cell `c`'s `(genes, counts)`.
    pub fn of(&self, c: usize) -> (&[u32], &[u32]) {
        let r = self.ptr[c]..self.ptr[c + 1];
        (&self.gene[r.clone()], &self.x[r])
    }

    /// From each cell's `(gene, count)` entries, zeros dropped.
    pub fn from_cells(genes: usize, cells: &[Vec<(u32, u32)>]) -> Self {
        let mut out = Counts {
            genes,
            ptr: vec![0],
            ..Counts::default()
        };
        for cell in cells {
            for &(g, x) in cell.iter().filter(|e| e.1 > 0) {
                out.gene.push(g);
                out.x.push(x);
            }
            out.ptr.push(out.gene.len());
        }
        out
    }
}

/// The model's sizes and priors.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Settings {
    pub cell_roles: usize,
    pub gene_roles: usize,
    /// The global weights' concentration `gamma` and rate `c0`.
    pub gamma: f64,
    pub weight_rate: f64,
    /// The blocks' rates' shape and rate `kappa` (mean 1).
    pub kappa: f64,
    /// The background's sizes' and levels' shape and rate.
    pub background: (f64, f64),
    /// The cells' and genes' rates' shape and rate.
    pub rate_prior: (f64, f64),
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            cell_roles: 8,
            gene_roles: 8,
            gamma: 1.0,
            weight_rate: 1.0,
            kappa: 1.0,
            background: (0.1, 0.1),
            rate_prior: (0.1, 0.1),
        }
    }
}

/// The model over one count matrix.
pub struct Model<'a> {
    counts: &'a Counts,
    settings: Settings,
}

/// A chain's parameters, and the units' allocation of its last sweep.
#[derive(Debug, Clone)]
pub struct State {
    /// `cells x A`, `genes x B`, `A x B`, row-major.
    pub theta: Vec<f64>,
    pub psi: Vec<f64>,
    pub lambda: Vec<f64>,
    /// Each cell's size and each gene's level in the background.
    pub size: Vec<f64>,
    pub level: Vec<f64>,
    /// Each cell's and each gene's rate in its weights' prior.
    pub cell_rate: Vec<f64>,
    pub gene_rate: Vec<f64>,
    /// The global weights over the cell roles and the gene roles.
    pub w: Vec<f64>,
    pub v: Vec<f64>,
    /// Units per cell and role, gene and role, block; and in the
    /// background per cell and per gene.
    pub n: Vec<u32>,
    pub m: Vec<u32>,
    pub blocks: Vec<u32>,
    pub n0: Vec<u32>,
    pub m0: Vec<u32>,
}

/// One draw of a Gamma by shape and rate, kept off 0.
fn gamma(shape: f64, rate: f64, rng: &mut SmallRng) -> f64 {
    Gamma::new(shape, 1.0 / rate)
        .expect("a positive shape and rate")
        .sample(rng)
        .max(1e-300)
}

/// The tables of `n` customers at concentration `r` (a Chinese
/// restaurant table count): the first opens one for certain, the `j`-th
/// after it one with probability `r / (r + j)`.
fn tables(n: u32, r: f64, rng: &mut SmallRng) -> u32 {
    if n == 0 {
        return 0;
    }
    1 + (1..n)
        .filter(|&j| rng.random::<f64>() * (r + j as f64) < r)
        .count() as u32
}

/// `n` units shared out over `weights` (all finite and at least 0, not all
/// 0), multinomially, into `out` (added to).
fn multinomial(n: u32, weights: &[f64], out: &mut [u32], rng: &mut SmallRng) {
    let mut left = n;
    let mut rest: f64 = weights.iter().sum();
    for (o, &w) in out.iter_mut().zip(weights) {
        if left == 0 {
            break;
        }
        let p = (w / rest).clamp(0.0, 1.0);
        let k = if p >= 1.0 {
            left
        } else {
            Binomial::new(left as u64, p).expect("a probability").sample(rng) as u32
        };
        *o += k;
        left -= k;
        rest -= w;
    }
}

impl<'a> Model<'a> {
    pub fn new(counts: &'a Counts, settings: Settings) -> anyhow::Result<Self> {
        anyhow::ensure!(
            settings.cell_roles >= 1 && settings.gene_roles >= 1,
            "at least one cell role and one gene role"
        );
        Ok(Self { counts, settings })
    }

    fn sizes(&self) -> (usize, usize, usize, usize) {
        let s = &self.settings;
        (self.counts.cells(), self.counts.genes, s.cell_roles, s.gene_roles)
    }

    /// Every parameter drawn from its prior.
    pub fn prior_draw(&self, rng: &mut SmallRng) -> State {
        let (cc, gg, aa, bb) = self.sizes();
        let s = &self.settings;
        let w: Vec<f64> = (0..aa).map(|_| gamma(s.gamma / aa as f64, s.weight_rate, rng)).collect();
        let v: Vec<f64> = (0..bb).map(|_| gamma(s.gamma / bb as f64, s.weight_rate, rng)).collect();
        let (e0, f0) = s.rate_prior;
        let cell_rate: Vec<f64> = (0..cc).map(|_| gamma(e0, f0, rng)).collect();
        let gene_rate: Vec<f64> = (0..gg).map(|_| gamma(e0, f0, rng)).collect();
        let theta = (0..cc * aa)
            .map(|i| gamma(w[i % aa], cell_rate[i / aa], rng))
            .collect();
        let psi = (0..gg * bb)
            .map(|i| gamma(v[i % bb], gene_rate[i / bb], rng))
            .collect();
        let lambda = (0..aa * bb).map(|_| gamma(s.kappa, s.kappa, rng)).collect();
        let (a0, b0) = s.background;
        State {
            theta,
            psi,
            lambda,
            size: (0..cc).map(|_| gamma(a0, b0, rng)).collect(),
            level: (0..gg).map(|_| gamma(a0, b0, rng)).collect(),
            cell_rate,
            gene_rate,
            w,
            v,
            n: vec![0; cc * aa],
            m: vec![0; gg * bb],
            blocks: vec![0; aa * bb],
            n0: vec![0; cc],
            m0: vec![0; gg],
        }
    }

    /// Each nonzero count's units shared out over the background and the
    /// blocks by their rates, the counts per cell and role, gene and role,
    /// block and background tallied afresh.
    pub fn allocate(&self, st: &mut State, rng: &mut SmallRng) {
        let (cc, _, aa, bb) = self.sizes();
        for t in [&mut st.n, &mut st.m, &mut st.blocks, &mut st.n0, &mut st.m0] {
            t.fill(0);
        }
        // Per link: the background's rate then each cell role's (its gene
        // roles summed), and within a role each gene role's.
        let mut role = vec![0.0; 1 + aa];
        let mut within = vec![0.0; aa * bb];
        let mut units = vec![0u32; 1 + aa];
        let mut per_b = vec![0u32; bb];
        for c in 0..cc {
            let (genes, xs) = self.counts.of(c);
            let theta = &st.theta[c * aa..(c + 1) * aa];
            for (&g, &x) in genes.iter().zip(xs) {
                let g = g as usize;
                let psi = &st.psi[g * bb..(g + 1) * bb];
                role[0] = st.size[c] * st.level[g];
                for a in 0..aa {
                    let row = &mut within[a * bb..(a + 1) * bb];
                    for b in 0..bb {
                        row[b] = st.lambda[a * bb + b] * psi[b];
                    }
                    role[1 + a] = theta[a] * row.iter().sum::<f64>();
                }
                units.fill(0);
                multinomial(x, &role, &mut units, rng);
                st.n0[c] += units[0];
                st.m0[g] += units[0];
                for a in 0..aa {
                    let k = units[1 + a];
                    if k == 0 {
                        continue;
                    }
                    st.n[c * aa + a] += k;
                    per_b.fill(0);
                    multinomial(k, &within[a * bb..(a + 1) * bb], &mut per_b, rng);
                    for (b, &kb) in per_b.iter().enumerate() {
                        st.m[g * bb + b] += kb;
                        st.blocks[a * bb + b] += kb;
                    }
                }
            }
        }
    }

    /// Everything but the allocation from its conditional, the allocation
    /// held: the cell side, the gene side, the blocks, the background.
    pub fn draw_parameters(&self, st: &mut State, rng: &mut SmallRng) {
        let (cc, gg, aa, bb) = self.sizes();
        let s = &self.settings;
        let (e0, f0) = s.rate_prior;
        // The cell side: each cell role's exposure per unit of theta.
        let psi_total: Vec<f64> = (0..bb).map(|b| (0..gg).map(|g| st.psi[g * bb + b]).sum()).collect();
        let cell_exposure: Vec<f64> = (0..aa)
            .map(|a| (0..bb).map(|b| st.lambda[a * bb + b] * psi_total[b]).sum())
            .collect();
        side(
            (&st.n, &mut st.theta, &mut st.cell_rate, &mut st.w),
            &cell_exposure,
            (cc, aa),
            (s.gamma, s.weight_rate, e0, f0),
            rng,
        );
        let theta_total: Vec<f64> = (0..aa).map(|a| (0..cc).map(|c| st.theta[c * aa + a]).sum()).collect();
        let gene_exposure: Vec<f64> = (0..bb)
            .map(|b| (0..aa).map(|a| st.lambda[a * bb + b] * theta_total[a]).sum())
            .collect();
        side(
            (&st.m, &mut st.psi, &mut st.gene_rate, &mut st.v),
            &gene_exposure,
            (gg, bb),
            (s.gamma, s.weight_rate, e0, f0),
            rng,
        );
        let psi_total: Vec<f64> = (0..bb).map(|b| (0..gg).map(|g| st.psi[g * bb + b]).sum()).collect();
        for a in 0..aa {
            for b in 0..bb {
                let k = a * bb + b;
                st.lambda[k] = gamma(
                    s.kappa + st.blocks[k] as f64,
                    s.kappa + theta_total[a] * psi_total[b],
                    rng,
                );
            }
        }
        let (a0, b0) = s.background;
        let levels: f64 = st.level.iter().sum();
        for c in 0..cc {
            st.size[c] = gamma(a0 + st.n0[c] as f64, b0 + levels, rng);
        }
        let sizes: f64 = st.size.iter().sum();
        for g in 0..gg {
            st.level[g] = gamma(a0 + st.m0[g] as f64, b0 + sizes, rng);
        }
    }

    /// Each count's log-likelihood under `st`, summed: every nonzero
    /// count's `x ln mu - ln x!` less every rate's exposure.
    pub fn log_likelihood(&self, st: &State) -> f64 {
        let (cc, gg, aa, bb) = self.sizes();
        let mut out = 0.0;
        for c in 0..cc {
            let (genes, xs) = self.counts.of(c);
            for (&g, &x) in genes.iter().zip(xs) {
                let g = g as usize;
                let mut mu = st.size[c] * st.level[g];
                for a in 0..aa {
                    for b in 0..bb {
                        mu += st.theta[c * aa + a] * st.lambda[a * bb + b] * st.psi[g * bb + b];
                    }
                }
                out += x as f64 * mu.ln() - statrs::function::gamma::ln_gamma(x as f64 + 1.0);
            }
        }
        let theta_total: Vec<f64> = (0..aa).map(|a| (0..cc).map(|c| st.theta[c * aa + a]).sum()).collect();
        let psi_total: Vec<f64> = (0..bb).map(|b| (0..gg).map(|g| st.psi[g * bb + b]).sum()).collect();
        for a in 0..aa {
            for b in 0..bb {
                out -= st.lambda[a * bb + b] * theta_total[a] * psi_total[b];
            }
        }
        out - st.size.iter().sum::<f64>() * st.level.iter().sum::<f64>()
    }
}

/// One side's draws, the cells' or the genes' (`units x roles`): the
/// global weights through table counts, each unit's weights integrated out
/// (a role's units in a unit are negative binomial, `NB(w_a, E_a / (rho +
/// E_a))`), then each unit's weights, then its rate.
fn side(
    (counts, weights, rate, global): (&[u32], &mut [f64], &mut [f64], &mut [f64]),
    exposure: &[f64],
    (units, roles): (usize, usize),
    (conc, conc_rate, e0, f0): (f64, f64, f64, f64),
    rng: &mut SmallRng,
) {
    for (a, w) in global.iter_mut().enumerate() {
        let (mut t, mut ln_q) = (0u64, 0.0);
        for u in 0..units {
            t += tables(counts[u * roles + a], *w, rng) as u64;
            ln_q += (1.0 + exposure[a] / rate[u]).ln();
        }
        *w = gamma(conc / roles as f64 + t as f64, conc_rate + ln_q, rng);
    }
    for u in 0..units {
        let row = &mut weights[u * roles..(u + 1) * roles];
        for (a, x) in row.iter_mut().enumerate() {
            *x = gamma(global[a] + counts[u * roles + a] as f64, rate[u] + exposure[a], rng);
        }
        let shape: f64 = global.iter().sum();
        rate[u] = gamma(e0 + shape, f0 + row.iter().sum::<f64>(), rng);
    }
}

impl Model<'_> {
    /// A start from the data: the background at half of each cell's and
    /// gene's total (`s_c r_g` half the counts' product over their sum),
    /// each cell's and gene's other half spread over its roles at random
    /// (each role's share a unit Gamma's), the blocks' rates and the global
    /// weights at 1, the rates at their priors' means.
    pub fn start(&self, rng: &mut SmallRng) -> State {
        let (cc, gg, aa, bb) = self.sizes();
        let mut st = self.prior_draw(rng);
        let mut cell_total = vec![0.0f64; cc];
        let mut gene_total = vec![0.0f64; gg];
        for c in 0..cc {
            let (genes, xs) = self.counts.of(c);
            for (&g, &x) in genes.iter().zip(xs) {
                cell_total[c] += x as f64;
                gene_total[g as usize] += x as f64;
            }
        }
        let total: f64 = cell_total.iter().sum::<f64>().max(1.0);
        // s_c r_g = half of N_c N_g / N, split evenly between the two.
        let root = (2.0 * total).sqrt();
        for c in 0..cc {
            st.size[c] = (cell_total[c] / root).max(1e-3);
        }
        for g in 0..gg {
            st.level[g] = (gene_total[g] / root).max(1e-3);
        }
        // Blocks at rate 1: theta_ca psi_gb summed is the other half when
        // theta sums to N_c / root over a cell's roles and psi to N_g /
        // root over a gene's.
        let spread = |out: &mut [f64], totals: &[f64], k: usize, rng: &mut SmallRng| {
            for (row, &t) in out.chunks_mut(k).zip(totals) {
                let shares: Vec<f64> = (0..k).map(|_| gamma(1.0, 1.0, rng)).collect();
                let sum: f64 = shares.iter().sum();
                for (x, w) in row.iter_mut().zip(&shares) {
                    *x = (t / root * w / sum).max(1e-6);
                }
            }
        };
        spread(&mut st.theta, &cell_total, aa, rng);
        spread(&mut st.psi, &gene_total, bb, rng);
        st.lambda.fill(1.0);
        st.w.fill(1.0);
        st.v.fill(1.0);
        let (e0, f0) = self.settings.rate_prior;
        st.cell_rate.fill(e0 / f0);
        st.gene_rate.fill(e0 / f0);
        st
    }
}

impl McmcModel for Model<'_> {
    type State = State;
    type Sample = f64;
    type Result = Vec<f64>;

    fn init(&self, rng: &mut SmallRng) -> State {
        self.start(rng)
    }

    fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        self.allocate(st, rng);
        self.draw_parameters(st, rng);
    }

    fn collect(&self, st: &State) -> f64 {
        self.log_likelihood(st)
    }

    fn summarize(&self, samples: Vec<f64>) -> Vec<f64> {
        samples
    }
}
