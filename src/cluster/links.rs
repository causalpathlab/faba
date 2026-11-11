//! Reads grouped into (cell, gene) links: a cell's reads at one gene's
//! sites, the unit the blockmodel labels.

use crate::cluster::read::{Reads, Weighted};

/// Reads grouped into (cell, gene) links: link `l` is cell `cell[l]`'s
/// entries `ptr[l]..ptr[l + 1]` (sites and weighted reads) at one gene.
pub struct Links {
    /// Cell `c`'s links are `cell_ptr[c]..cell_ptr[c + 1]`.
    pub cell_ptr: Vec<usize>,
    pub cell: Vec<u32>,
    pub ptr: Vec<usize>,
    pub site: Vec<u32>,
    pub on: Vec<f64>,
    pub off: Vec<f64>,
    pub cells: usize,
}

impl Links {
    /// `reads` grouped by `gene_of[site]` within each cell.
    pub fn new(reads: &Weighted, gene_of: &[u32]) -> Self {
        let mut links = Links {
            cell_ptr: vec![0],
            cell: Vec::new(),
            ptr: vec![0],
            site: Vec::with_capacity(reads.on.len()),
            on: Vec::with_capacity(reads.on.len()),
            off: Vec::with_capacity(reads.on.len()),
            cells: reads.cells(),
        };
        for c in 0..reads.cells() {
            let (s, m, u) = reads.of(c);
            let mut order: Vec<usize> = (0..s.len()).collect();
            order.sort_by_key(|&i| (gene_of[s[i] as usize], s[i]));
            for (j, &i) in order.iter().enumerate() {
                let gene = gene_of[s[i] as usize];
                if j > 0 && gene_of[s[order[j - 1]] as usize] != gene {
                    links.cell.push(c as u32);
                    links.ptr.push(links.site.len());
                }
                links.site.push(s[i]);
                links.on.push(m[i]);
                links.off.push(u[i]);
            }
            if !order.is_empty() {
                links.cell.push(c as u32);
                links.ptr.push(links.site.len());
            }
            links.cell_ptr.push(links.cell.len());
        }
        links
    }

    /// The links [`Self::new`] would group `reads` into, counted without
    /// grouping them: each cell's distinct genes.
    pub fn count(reads: &Reads, gene_of: &[u32]) -> usize {
        let genes = gene_of.iter().max().map_or(0, |&g| g as usize + 1);
        // The last cell each gene was met in, plus 1 (0: none yet).
        let mut seen = vec![0usize; genes];
        let mut links = 0;
        for c in 0..reads.cells() {
            for &s in reads.of(c).0 {
                let g = gene_of[s as usize] as usize;
                if seen[g] != c + 1 {
                    seen[g] = c + 1;
                    links += 1;
                }
            }
        }
        links
    }

    /// Each link's gene, by `gene_of[site]` of its first site.
    pub fn genes_of(&self, gene_of: &[u32]) -> Vec<u32> {
        (0..self.len())
            .map(|l| gene_of[self.site[self.ptr[l]] as usize])
            .collect()
    }

    pub fn len(&self) -> usize {
        self.cell.len()
    }

    /// Link `l`'s `(sites, converted, unconverted)`.
    pub(in crate::cluster) fn of(&self, l: usize) -> (&[u32], &[f64], &[f64]) {
        let r = self.ptr[l]..self.ptr[l + 1];
        (&self.site[r.clone()], &self.on[r.clone()], &self.off[r])
    }
}
