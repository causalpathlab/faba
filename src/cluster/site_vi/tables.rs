//! What the tables give of a fit: the posterior means the blockmodel's
//! writers take ([`Means`]), from one last pass over every cell's links at
//! the fit's factors; and the cells kept but not fitted, each placed by a
//! local step against the fitted globals.

use super::{dirichlet_term, Cell, Data, Fit, Fitted, Scratch, Xi, LOCAL_ITERS};

/// Cells whose links' `phi` the tables pass holds at once.
const CHUNK: usize = 512;
use crate::cluster::means::Means;
use rayon::prelude::*;

impl Fit<'_> {
    /// The tables' means of `fitted`: each cell's and gene's share of its
    /// links by role, each block's shift and its sd, each cell's level, and
    /// each gene's, cell's and site's links or reads in each block.
    pub fn means(&self, fitted: &Fitted) -> Means {
        let sites = self.data.mu.len();
        let pr = self.prior;
        let (aa, bb, k) = (pr.cell_roles, pr.gene_roles, self.blocks());
        let (g, l) = (&fitted.globals, &fitted.locals);
        let cache = self.cache(g);
        let a0 = pr.alpha / aa as f64;
        // Each cell's links' phi, a chunk of cells at a time (each chunk's in
        // parallel), added in the cells' order: its own row of links per
        // block, and into the genes' and the sites'.
        let cells = self.data.cells();
        let mut cell_block = vec![0.0; cells * k];
        let mut gene_block = vec![0.0; self.data.genes * k];
        let (mut site_on, mut site_off) = (vec![0.0; sites * k], vec![0.0; sites * k]);
        let all: Vec<usize> = (0..cells).collect();
        for chunk in all.chunks(CHUNK) {
            let rows: Vec<Vec<(usize, Vec<f64>)>> = chunk
                .par_iter()
                .map(|&c| {
                    let (_, eln) = dirichlet_term(&l.t[c * aa..(c + 1) * aa], a0);
                    let links = self.data.links_of(c);
                    let first = self.data.links.ptr[links.start];
                    let mut sc = Scratch::default();
                    sc.fixed
                        .resize((self.data.links.ptr[links.end] - first) * k, [0.0; 3]);
                    links
                        .map(|link| {
                            let at = self.data.links.ptr[link] - first;
                            self.sums(link, (l.e[c], l.s[c]), g, Xi::Fresh, at, &mut sc);
                            self.phi(link, &eln, &cache, &mut sc);
                            (link, sc.phi.clone())
                        })
                        .collect()
                })
                .collect();
            for (&c, links) in chunk.iter().zip(&rows) {
                for (link, phi) in links {
                    let gene = self.data.link_gene[*link] as usize;
                    for (b, p) in phi.iter().enumerate() {
                        cell_block[c * k + b] += p;
                        gene_block[gene * k + b] += p;
                    }
                    let (s, m, u) = self.data.links.of(*link);
                    for i in 0..s.len() {
                        let at = s[i] as usize * k;
                        for (b, p) in phi.iter().enumerate() {
                            site_on[at + b] += p * m[i];
                            site_off[at + b] += p * u[i];
                        }
                    }
                }
            }
        }
        // A row's shares of its links by role: cells over the cell roles,
        // genes over the gene roles.
        let shares = |block: &[f64], rows: usize, roles: usize, of: fn(usize, usize) -> usize| {
            let mut out = vec![0f32; rows * roles];
            for (r, row) in block.chunks_exact(k).enumerate() {
                let total: f64 = row.iter().sum();
                for (b, x) in row.iter().enumerate() {
                    out[r * roles + of(b, bb)] += match total > 0.0 {
                        true => (x / total) as f32,
                        false => 1.0 / k as f32,
                    };
                }
            }
            out
        };
        Means {
            cell: shares(&cell_block, cells, aa, |b, bb| b / bb),
            gene: shares(&gene_block, self.data.genes, bb, |b, bb| b % bb),
            shift: g.d.clone(),
            shift_sd: g.v.iter().map(|v| v.sqrt()).collect(),
            level: l.e.clone(),
            gene_block,
            cell_block,
            site_on,
            site_off,
        }
    }

    /// The cells of `data` (cells kept but not fitted, over the same sites
    /// and genes) placed against `fitted`'s globals: each a local step from
    /// its links spread evenly and its level at 0; each one's share of its
    /// links by cell role (cells x A) and its level.
    pub fn place(&self, data: &Data, fitted: &Fitted) -> (Vec<f32>, Vec<f64>) {
        let fit = Fit {
            data,
            prior: self.prior,
        };
        let aa = self.prior.cell_roles;
        let g = &fitted.globals;
        let cache = fit.cache(g);
        let a0 = self.prior.alpha / aa as f64;
        let cell_var = self.prior.cell_sd * self.prior.cell_sd;
        let placed: Vec<Cell> = (0..data.cells())
            .into_par_iter()
            .map(|c| {
                let own = data.links_of(c).len() as f64;
                let start = Cell {
                    t: vec![a0 + own / aa as f64; aa],
                    e: 0.0,
                    s: cell_var,
                };
                fit.step(c, start, g, &cache, 2 * LOCAL_ITERS).0
            })
            .collect();
        let weights = placed
            .iter()
            .flat_map(|cell| {
                let n: Vec<f64> = cell.t.iter().map(|t| (t - a0).max(0.0)).collect();
                let total: f64 = n.iter().sum();
                n.into_iter().map(move |x| match total > 0.0 {
                    true => (x / total) as f32,
                    false => 1.0 / aa as f32,
                })
            })
            .collect();
        (weights, placed.iter().map(|c| c.e).collect())
    }
}
