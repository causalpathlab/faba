//! What the site blockmodel's tables are written from: its means.

use super::count::column_totals;
use super::math::{logit, pooled_fraction};

/// The fit's means: each cell's and gene's role weights (row-major), each
/// block's shift with its sd, each cell's level, and each gene's, cell's and site's links or reads in each block
/// (`K` blocks per row; the sites' only when asked for).
pub struct Means {
    pub cell: Vec<f32>,
    pub gene: Vec<f32>,
    pub shift: Vec<f64>,
    pub shift_sd: Vec<f64>,
    pub level: Vec<f64>,
    pub gene_block: Vec<f64>,
    pub cell_block: Vec<f64>,
    pub site_on: Vec<f64>,
    pub site_off: Vec<f64>,
}

/// Each row's share of its counts in each block over the block's share of
/// all counts: above 1, the row is enriched in the block (0 for an empty
/// row).
fn enrichment(counts: &[f64], blocks: usize) -> Vec<f64> {
    let column = column_totals(counts, blocks);
    let total: f64 = column.iter().sum();
    counts
        .chunks(blocks)
        .flat_map(|row| {
            let n: f64 = row.iter().sum();
            row.iter().zip(&column).map(move |(&x, &c)| {
                if n > 0.0 && c > 0.0 {
                    (x / n) / (c / total)
                } else {
                    0.0
                }
            })
        })
        .collect()
}

impl Means {
    /// Blocks, `A x B`.
    pub fn blocks(&self) -> usize {
        self.shift.len()
    }

    /// Each block's links, over all genes.
    pub fn block_links(&self) -> Vec<f64> {
        column_totals(&self.gene_block, self.blocks())
    }

    /// Each gene's enrichment in each block.
    pub fn gene_enrichment(&self) -> Vec<f64> {
        enrichment(&self.gene_block, self.blocks())
    }

    /// Each cell's enrichment in each block.
    pub fn cell_enrichment(&self) -> Vec<f64> {
        enrichment(&self.cell_block, self.blocks())
    }

    /// Each site's methylation in each block over its baseline `mu`, on the
    /// logit scale (`NaN` where the block holds none of its reads).
    pub fn site_lift(&self, mu: &[f64]) -> Vec<f64> {
        self.site_on
            .iter()
            .zip(&self.site_off)
            .enumerate()
            .map(|(at, (&on, &off))| {
                if on + off <= 0.0 {
                    return f64::NAN;
                }
                logit(pooled_fraction(on, on + off)) - mu[at / self.blocks()]
            })
            .collect()
    }
}
