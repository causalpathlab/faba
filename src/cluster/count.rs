//! A Poisson blockmodel of gene counts: cell roles and gene roles over a
//! cell x gene count matrix, after Ball, Karrer and Newman (2011), the
//! Bayesian Poisson Tucker decomposition (Schein et al. 2016) and
//! hierarchical Poisson factorization (Gopalan, Hofman and Blei 2015),
//! fitted by variational inference ([`vi`]).
//!
//! ```text
//! x_cg ~ Poisson(sum_ab theta_ca lambda_ab psi_gb)
//! theta_ca ~ Gamma(eta, rho_c),   psi_gb ~ Gamma(eta, rho'_g),   lambda_ab ~ Gamma(kappa, kappa)
//! rho_c, rho'_g ~ Gamma(e0, f0)
//! ```
//!
//! (Gamma by shape and rate.) A cell's and a gene's units spread over
//! roles, so cells and genes have mixed memberships; each one's rate
//! takes its depth. Only nonzero counts are visited: a zero enters through
//! the exposures, `sum_ab lambda_ab (sum_c theta_ca)(sum_g psi_gb)`.

pub mod vi;

/// A cell x gene count matrix, its nonzero entries by cell: cell `c`'s are
/// `ptr[c]..ptr[c + 1]` of `gene` and `x`.
#[derive(Debug, Clone, Default)]
pub struct Counts {
    pub genes: usize,
    pub ptr: Vec<usize>,
    pub gene: Vec<u32>,
    pub x: Vec<u32>,
    /// `sum ln x!` over the nonzero counts: the likelihood's constant.
    ln_factorials: f64,
    /// `sum x` over the nonzero counts.
    total: f64,
}

impl Counts {
    /// Every count summed.
    pub fn total(&self) -> f64 {
        self.total
    }

    pub fn cells(&self) -> usize {
        self.ptr.len().saturating_sub(1)
    }

    /// Cell `c`'s `(genes, counts)`.
    pub fn of(&self, c: usize) -> (&[u32], &[u32]) {
        let r = self.ptr[c]..self.ptr[c + 1];
        (&self.gene[r.clone()], &self.x[r])
    }

    /// From each cell's `(gene, count)` entries, zeros dropped.
    pub fn from_cells<'c>(genes: usize, cells: impl IntoIterator<Item = &'c [(u32, u32)]>) -> Self {
        let mut out = Counts {
            genes,
            ptr: vec![0],
            ..Counts::default()
        };
        for cell in cells {
            for &(g, x) in cell.iter().filter(|e| e.1 > 0) {
                out.gene.push(g);
                out.x.push(x);
                out.ln_factorials += statrs::function::gamma::ln_gamma(x as f64 + 1.0);
                out.total += x as f64;
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
    /// The role weights' fixed shape `eta` (HPF's), on both sides.
    pub eta: f64,
    /// The blocks' rates' shape and rate `kappa` (mean 1).
    pub kappa: f64,
    /// The cells' and genes' rates' shape and rate.
    pub rate_prior: (f64, f64),
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            cell_roles: 8,
            gene_roles: 8,
            eta: 0.3,
            kappa: 1.0,
            rate_prior: (0.1, 0.1),
        }
    }
}

/// Each column's total of `x`, `k` columns per row.
pub(in crate::cluster) fn column_totals(x: &[f64], k: usize) -> Vec<f64> {
    let mut out = vec![0.0; k];
    for row in x.chunks_exact(k) {
        out.iter_mut().zip(row).for_each(|(o, v)| *o += v);
    }
    out
}

/// `f(sum, i)` over `0..n` into a sum of `len` zeros, in `parts` fixed
/// parts run in parallel and added in order: the same sum, to the bit,
/// whatever the threads.
pub(in crate::cluster) fn sum_parts<T>(
    n: usize,
    len: usize,
    parts: usize,
    f: impl Fn(&mut [T], usize) + Sync,
) -> Vec<T>
where
    T: Copy + Default + std::ops::AddAssign + Send,
{
    use rayon::prelude::*;
    let size = n.div_ceil(parts.max(1)).max(1);
    let sums: Vec<Vec<T>> = (0..n.div_ceil(size))
        .into_par_iter()
        .map(|p| {
            let mut s = vec![T::default(); len];
            (p * size..((p + 1) * size).min(n)).for_each(|i| f(&mut s, i));
            s
        })
        .collect();
    let mut out = vec![T::default(); len];
    for s in sums {
        out.iter_mut().zip(s).for_each(|(o, v)| *o += v);
    }
    out
}
