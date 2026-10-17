//! What a run's chains say together: their agreement, the best chain's
//! shifts and weights, and how far their strong gene roles name the same
//! genes.

use super::draw::enrichment;
use super::{shift_range, Draw, Means, Model};
use crate::cluster::summary::{chain_scores, main_labels, pairwise, pairwise_ari, per_read};

/// What the chains say, together.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Summary {
    pub cell_roles: usize,
    pub gene_roles: usize,
    pub chains: usize,
    pub draws: usize,
    /// Split R-hat of the log joint across chains.
    pub rhat_log_joint: f64,
    /// Each chain's mean log joint.
    pub chain_log_joint: Vec<f64>,
    /// Adjusted Rand index between each pair of chains' cells' main cell
    /// roles, and genes' main gene roles (posterior-mean weights).
    pub cell_agreement: Vec<f64>,
    pub gene_agreement: Vec<f64>,
    /// The best chain (highest mean log joint) and its posterior-mean block
    /// shifts (`A x B`, row-major) with their sds.
    pub best_chain: usize,
    pub shifts: Vec<f64>,
    pub shift_sd: Vec<f64>,
    /// The best chain's posterior-mean global weights over the cell roles.
    pub global_weights: Vec<f64>,
    /// Held-out reads' log-likelihood per read, over every kept draw.
    pub heldout_per_read: Option<f64>,
    /// Each chain's strong gene roles (whose shift ranges at least
    /// [`STRONG_SHIFT`] across the cell roles), and for each pair of chains
    /// the mean overlap of a strong role's top genes with the best-matching
    /// strong role's in the other chain.
    pub strong_gene_roles: Vec<usize>,
    pub top_gene_overlap: Vec<f64>,
}

/// A gene role is strong when its shift ranges at least this much (logit
/// scale) across the cell roles.
const STRONG_SHIFT: f64 = 0.5;
/// Top genes per strong gene role, compared across chains.
const TOP_GENES: usize = 50;
/// Links a gene needs in a chain to be ranked.
const MIN_LINKS: f64 = 5.0;

/// Each strong gene role's top genes in a chain: the [`TOP_GENES`] with at
/// least [`MIN_LINKS`] links most enriched in the role.
fn strong_top_genes(m: &Means, aa: usize, bb: usize) -> Vec<Vec<usize>> {
    let k = aa * bb;
    let genes = m.gene_block.len() / k.max(1);
    // Links per gene and gene role, cell roles summed.
    let mut role = vec![0.0f64; genes * bb];
    for (at, &x) in m.gene_block.iter().enumerate() {
        role[(at / k) * bb + (at % k) % bb] += x;
    }
    let enrich = enrichment(&role, bb);
    (0..bb)
        .filter(|&b| shift_range(&m.shift, aa, bb, b) >= STRONG_SHIFT)
        .map(|b| {
            let mut ranked: Vec<usize> = (0..genes)
                .filter(|&g| role[g * bb..(g + 1) * bb].iter().sum::<f64>() >= MIN_LINKS)
                .collect();
            ranked.sort_by(|&x, &y| enrich[y * bb + b].total_cmp(&enrich[x * bb + b]));
            ranked.truncate(TOP_GENES);
            ranked
        })
        .collect()
}

/// The mean, over `a`'s top-gene sets, of the best overlap (shared share)
/// with any of `b`'s (0 when either has none).
fn best_overlap(a: &[Vec<usize>], b: &[Vec<usize>]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let shared = |x: &[usize], y: &[usize]| {
        x.iter().filter(|g| y.contains(g)).count() as f64 / x.len().max(1) as f64
    };
    a.iter()
        .map(|x| b.iter().map(|y| shared(x, y)).fold(0.0, f64::max))
        .sum::<f64>()
        / a.len() as f64
}

/// The chains of `model`, summarised together.
pub fn summarize(model: &Model, chains: &[Vec<Draw>]) -> Summary {
    let (aa, bb) = (model.settings.cell_roles, model.settings.gene_roles);
    let log_joint: Vec<Vec<f64>> = chains
        .iter()
        .map(|c| c.iter().map(|d| d.log_joint).collect())
        .collect();
    let (chain_log_joint, rhat_log_joint, best_chain) = chain_scores(&log_joint);
    let means: Vec<Means> = chains
        .iter()
        .filter_map(|c| c.last().map(|d| d.means(false)))
        .collect();
    let tops: Vec<Vec<Vec<usize>>> = means.iter().map(|m| strong_top_genes(m, aa, bb)).collect();
    let best = means.get(best_chain);
    Summary {
        cell_roles: aa,
        gene_roles: bb,
        chains: chains.len(),
        draws: chains.iter().map(Vec::len).sum(),
        rhat_log_joint,
        chain_log_joint,
        cell_agreement: pairwise_ari(
            &means
                .iter()
                .map(|m| main_labels(&m.cell, aa))
                .collect::<Vec<_>>(),
        ),
        gene_agreement: pairwise_ari(
            &means
                .iter()
                .map(|m| main_labels(&m.gene, bb))
                .collect::<Vec<_>>(),
        ),
        best_chain,
        shifts: best.map_or_else(Vec::new, |m| m.shift.clone()),
        shift_sd: best.map_or_else(Vec::new, |m| m.shift_sd.clone()),
        global_weights: best.map_or_else(Vec::new, |m| m.beta.clone()),
        heldout_per_read: per_read(chains.iter().flatten().map(|d| d.heldout)),
        strong_gene_roles: tops.iter().map(Vec::len).collect(),
        top_gene_overlap: pairwise(&tops, |x, y| best_overlap(x, y)),
    }
}
