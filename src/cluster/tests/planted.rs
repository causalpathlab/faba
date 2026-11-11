//! Planted site reads for the fits' tests (made up, small): blocks of cell
//! roles and gene roles shifting their sites' methylation; and the links
//! they group into.

use super::draw_cells;
use crate::cluster::links::Links;
use crate::cluster::read::{self, Reads, Weighted};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Cells of three planted cell roles, genes of two planted gene roles,
/// each (cell role, gene role) block shifting its sites' methylation, and
/// each cell its own level within `±level`: the reads, each site's gene,
/// the cells' and genes' roles, and the cells' levels.
pub(super) fn planted_blocks(
    cells: usize,
    genes: usize,
    level: f64,
    seed: u64,
) -> (Reads, Vec<u32>, Vec<u16>, Vec<u16>, Vec<f64>) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let sites = genes * 4;
    let gene_of: Vec<u32> = (0..sites as u32).map(|s| s / 4).collect();
    let mu: Vec<f64> = (0..sites).map(|_| rng.random_range(-1.5..0.5)).collect();
    let shift = [[1.5, -1.5], [-1.5, 1.5], [0.0, 0.0]];
    let cell_role: Vec<u16> = (0..cells).map(|c| (c % 3) as u16).collect();
    let gene_role: Vec<u16> = (0..genes).map(|g| (g % 2) as u16).collect();
    // Each cell's own level, whatever its role.
    let level: Vec<f64> = (0..cells)
        .map(|_| level * rng.random_range(-1.0..1.0))
        .collect();
    let fracs: Vec<Vec<f64>> = cell_role
        .iter()
        .zip(&level)
        .map(|(&r, &v)| {
            (0..sites)
                .map(|s| {
                    sigmoid(mu[s] + v + shift[r as usize][gene_role[gene_of[s] as usize] as usize])
                })
                .collect()
        })
        .collect();
    let refs: Vec<&[f64]> = fracs.iter().map(Vec::as_slice).collect();
    (
        draw_cells(&refs, 0.5, &mut rng),
        gene_of,
        cell_role,
        gene_role,
        level,
    )
}

/// [`planted_blocks`] (each cell's level within `±1`) as `faba cluster`
/// reads them: `cells` cells, the first `fitted` of them fitted and the
/// rest controls, over `genes` genes; the data and the fitted cells' reads.
pub(super) fn planted_site_data(
    cells: usize,
    fitted: usize,
    genes: usize,
    seed: u64,
) -> (read::SiteData, Reads) {
    let (reads, gene_of, _, _, _) = planted_blocks(cells, genes, 1.0, seed);
    let sites = gene_of
        .iter()
        .enumerate()
        .map(|(s, &g)| read::Site {
            gene: format!("GENE{}", g + 1).into(),
            chr: "chr1".into(),
            pos: s as i64 * 10,
        })
        .collect();
    let cells = (0..cells)
        .map(|c| super::cell(&format!("B{c}"), c < fitted))
        .collect();
    let fit_reads = reads.subset(&(0..fitted).collect::<Vec<_>>());
    let data = read::SiteData {
        modality: "m6A".into(),
        sites,
        cells,
        reads,
    };
    (data, fit_reads)
}

/// The links counted ahead of the fit, for its memory estimate, are those
/// built, genes interleaved within a cell or not.
#[test]
fn links_are_counted_before_they_are_built() {
    let (reads, gene_of, _, _, _) = planted_blocks(30, 5, 0.0, 161);
    let interleaved: Vec<u32> = (0..gene_of.len() as u32).map(|s| s % 3).collect();
    for genes in [&gene_of[..], &interleaved] {
        let built = Links::new(&Weighted::new(&reads, None), genes).len();
        assert_eq!(Links::count(&reads, genes), built);
    }
}

/// A cell's reads are grouped by gene into links, and each link knows its
/// gene.
#[test]
fn links_group_a_cells_reads_by_gene() {
    // Sites 0 and 2 in gene 0, site 1 in gene 1.
    let reads = Reads::from_cells(vec![vec![(0, 1, 0), (1, 0, 1), (2, 2, 0)], vec![(1, 1, 1)]]);
    let weighted = Weighted::new(&reads, None);
    let gene_of = [0, 1, 0];
    let links = Links::new(&weighted, &gene_of);
    assert_eq!(links.len(), 3);
    assert_eq!(links.cell, vec![0, 0, 1]);
    assert_eq!(links.genes_of(&gene_of), vec![0, 1, 1]);
    assert_eq!(links.cell_ptr, vec![0, 2, 3]);
    // Cell 0's gene-0 link holds both of its gene-0 sites.
    assert_eq!(&links.site[links.ptr[0]..links.ptr[1]], &[0, 2]);
}
