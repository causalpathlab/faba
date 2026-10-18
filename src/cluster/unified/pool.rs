//! The chains' tallies pooled for the tables, their roles aligned to the
//! best chain's, where they agree; else the best chain's alone.
//!
//! A chain the convergence verdict finds lagging on its trace (see
//! [`crate::cluster::progress::advice::lagging`]) is left out, as is one
//! whose roles in use do not match the best chain's (some role's match
//! under [`MIN_MATCH`]), each reason that holds named; the rest pool when
//! the split R-hat of their log joints together is within [`RHAT_OK`],
//! else each is left out. Over the chains pooled, each block shift and
//! each cell role's mean membership in use is checked as well, aligned:
//! split R-hat above [`RHAT_OK`] or fewer than [`MIN_ESS`] effective draws
//! across the chains is reported as not converged, and held against the
//! fit by the verdict, as is each chain left out on its roles. The split
//! R-hat is `legume_numeric`'s, without the rank normalisation of Vehtari
//! et al. (2021).

use super::align::{align, Alignment};
use super::draw::{log_joints, Tally};
use super::{Draw, Means, Model};
use crate::cluster::summary::{centred, centred_rhat, finite_or_inf};
use crate::cluster::RHAT_OK;
use legume_numeric::mcmc::engine::ess;

/// Effective draws, summed over the chains compared, below which a
/// quantity is not converged.
const MIN_ESS: f64 = 100.0;
/// Weighted Jaccard below which a role does not match the best chain's.
const MIN_MATCH: f64 = 0.5;

/// A chain left out of the pool, and why (every reason that holds, joined
/// by `; `).
#[derive(Debug, Clone, serde::Serialize)]
pub struct Excluded {
    pub chain: usize,
    pub reason: String,
}

/// Which chains the tables pool and how far they agree: the reference
/// (the best chain) and the chains pooled, in order (the reference alone
/// when nothing pools, with why); every other chain left out, with why;
/// each chain's alignment to the reference, by chain; and over the chains
/// compared (the reference and those not left out on their own), the
/// split R-hat of the log joint; over the chains pooled, the worst split
/// R-hat and fewest effective draws (summed over the chains) of any
/// quantity checked, and the quantities short of convergence; and what
/// the convergence verdict is to hold against the fit besides its traces
/// (see [`crate::cluster::progress::advice::Verdict::with_unsettled`]):
/// those quantities, then each chain left out but not lagging on its
/// trace, which the verdict names itself.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Pooling {
    pub reference: usize,
    pub pooled: Vec<usize>,
    pub fallback: Option<String>,
    pub excluded: Vec<Excluded>,
    pub alignment: Vec<Alignment>,
    #[serde(serialize_with = "finite_or_inf")]
    pub rhat_log_joint: f64,
    #[serde(serialize_with = "finite_or_inf")]
    pub worst_rhat: f64,
    pub worst_ess: f64,
    pub not_converged: Vec<String>,
    pub unsettled: Vec<String>,
}

/// The tables' means from `model`'s `chains`, `reference` the best and
/// `lagging` those the convergence verdict finds lagging on their traces:
/// the chains' tallies summed, aligned to the reference's roles, where
/// they pool; else the reference's own. With what was pooled.
pub fn pool(
    model: &Model,
    chains: &[Vec<Draw>],
    reference: usize,
    lagging: &[usize],
) -> (Means, Pooling) {
    let (aa, bb) = (model.settings.cell_roles, model.settings.gene_roles);
    let last = |c: usize| chains[c].last().expect("a chain has draws");
    let means: Vec<Means> = (0..chains.len()).map(|c| last(c).means(false)).collect();
    let alignment: Vec<Alignment> = (0..chains.len())
        .map(|c| {
            if c == reference {
                Alignment::identity(&means[c], aa, bb)
            } else {
                align(&means[reference], &means[c], aa, bb)
            }
        })
        .collect();
    let log_joint = log_joints(chains);
    let (compared, mut excluded) = agreeing(&alignment, reference, lagging);
    let rhat_log_joint = centred_rhat(&pick(&log_joint, &compared));
    let fallback = if chains.len() < 2 {
        Some("one chain".to_string())
    } else if compared.len() < 2 {
        Some("no other chain agrees with the best".to_string())
    } else if rhat_log_joint > RHAT_OK {
        let why =
            format!("split R-hat {rhat_log_joint:.2} of the log joint across the chains kept");
        excluded.extend(
            compared
                .iter()
                .filter(|&&c| c != reference)
                .map(|&c| Excluded {
                    chain: c,
                    reason: format!("the chains kept disagree: {why}"),
                }),
        );
        excluded.sort_by_key(|e| e.chain);
        Some(why)
    } else {
        None
    };
    let pooled = if fallback.is_some() {
        vec![reference]
    } else {
        compared.clone()
    };
    let mut tally = Tally::default();
    for &c in &pooled {
        let t = last(c).tally.lock().expect("a chain's tally");
        tally.add_aligned(&t, &alignment[c]);
    }
    let checks = diagnose(model, chains, (&pooled, reference), &alignment);
    let not_converged: Vec<String> = checks
        .iter()
        .filter(|c| c.1 > RHAT_OK || c.2 < MIN_ESS)
        .map(|c| c.0.clone())
        .collect();
    let left_out = excluded
        .iter()
        .filter(|e| !lagging.contains(&e.chain))
        .map(|e| format!("chain {} left out", e.chain));
    let pooling = Pooling {
        reference,
        rhat_log_joint,
        worst_rhat: checks.iter().map(|c| c.1).fold(1.0, f64::max),
        worst_ess: checks.iter().map(|c| c.2).fold(f64::INFINITY, f64::min),
        unsettled: not_converged.iter().cloned().chain(left_out).collect(),
        not_converged,
        pooled,
        fallback,
        excluded,
        alignment,
    };
    (tally.means(true), pooling)
}

/// The chains that agree with the reference on their own (it among them,
/// in order), from each chain's alignment to it and the chains `lagging`
/// on their traces; and those that do not, with every reason that holds:
/// a trace that lags, roles that do not match its. Not a chain's split
/// R-hat with the reference alone: over two chains' few draws it often
/// passes the bound by chance, and how the chains mix is the R-hat of
/// those kept together.
fn agreeing(
    alignment: &[Alignment],
    reference: usize,
    lagging: &[usize],
) -> (Vec<usize>, Vec<Excluded>) {
    let (mut kept, mut excluded) = (Vec::new(), Vec::new());
    for (c, to) in alignment.iter().enumerate() {
        let mut reasons = Vec::new();
        if c != reference {
            if lagging.contains(&c) {
                reasons
                    .push("its trace lags the best chain's (the convergence verdict)".to_string());
            }
            let worst = to.worst();
            if worst < MIN_MATCH {
                reasons.push(format!(
                    "its roles do not match the best chain's: a role's members overlap {worst:.2}"
                ));
            }
        }
        if reasons.is_empty() {
            kept.push(c);
        } else {
            excluded.push(Excluded {
                chain: c,
                reason: reasons.join("; "),
            });
        }
    }
    (kept, excluded)
}

/// `series` of chains `chains`.
fn pick(series: &[Vec<f64>], chains: &[usize]) -> Vec<Vec<f64>> {
    chains.iter().map(|&c| series[c].clone()).collect()
}

/// Each quantity checked over chains `pooled`, in the roles of the
/// `reference` among them: its name, split R-hat, and effective
/// draws summed over the chains. The log joint, each block's shift (with
/// the shift on), each cell role's mean membership; of the roles, only
/// those the reference uses (see [`Alignment::worst`]): an empty role's
/// values are its prior's.
fn diagnose(
    model: &Model,
    chains: &[Vec<Draw>],
    (pooled, reference): (&[usize], usize),
    alignment: &[Alignment],
) -> Vec<(String, f64, f64)> {
    let (aa, bb) = (model.settings.cell_roles, model.settings.gene_roles);
    // Chain `c`'s series of `f` at each kept draw.
    let series = |f: &dyn Fn(usize, &Draw) -> f64| -> Vec<Vec<f64>> {
        pooled
            .iter()
            .map(|&c| chains[c].iter().map(|d| f(c, d)).collect())
            .collect()
    };
    let check = |name: String, x: Vec<Vec<f64>>| {
        let total = centred(&x).iter().map(|c| f64::from(ess(c))).sum();
        (name, centred_rhat(&x), total)
    };
    let mut out = vec![check("log joint".into(), series(&|_, d| d.log_joint))];
    // By the reference's block and cell role, each chain's that is it.
    let (mut block, mut role) = (
        vec![vec![0; chains.len()]; aa * bb],
        vec![vec![0; chains.len()]; aa],
    );
    for &c in pooled {
        let to = &alignment[c];
        (0..aa * bb).for_each(|k| block[to.block(k)][c] = k);
        to.cell
            .iter()
            .enumerate()
            .for_each(|(a, &r)| role[r][c] = a);
    }
    let used = &alignment[reference];
    if model.settings.shift_sd.is_some() {
        for (k, from) in block.iter().enumerate() {
            let (a, b) = (k / bb, k % bb);
            if !(used.cell_used[a] && used.gene_used[b]) {
                continue;
            }
            let x = series(&|c, d| d.shift[from[c]]);
            out.push(check(format!("shift ({a}, {b})"), x));
        }
    }
    for (a, from) in role.iter().enumerate() {
        if !used.cell_used[a] {
            continue;
        }
        let x = series(&|c, d| d.role_mean[from[c]]);
        out.push(check(format!("cell role {a}"), x));
    }
    out
}
