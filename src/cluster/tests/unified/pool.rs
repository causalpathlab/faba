//! The chains' roles aligned to the best chain's, and their tallies pooled
//! where they agree: copies and relabelled copies of one chain pool to its
//! own means, a chain the convergence verdict finds lagging or one whose
//! roles in use do not match is left out and named with every reason that
//! holds, as are chains that disagree together, but not a chain on its
//! pair's R-hat with the best alone; with nothing to pool the tables are
//! the best chain's alone, and the diagnostics its own; and the record's
//! summary gives the tables' means.

use super::{blocks_preset, config, planted_blocks, settings, tiny_links, TINY_GENES};
use crate::cluster::read::Weighted;
use crate::cluster::summary::centred_rhat;
use crate::cluster::unified::{self, align, assign, Alignment, Draw, Links, Means, Model, Weights};
use crate::cluster::RHAT_OK;
use legume_numeric::mcmc::engine::run_mcmc;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};

/// Means holding only role memberships (`aa` cell roles, `bb` gene roles)
/// and block shifts, what alignment reads.
fn memberships(cell: Vec<f32>, gene: Vec<f32>, shift: Vec<f64>) -> Means {
    Means {
        cell,
        gene,
        beta: Vec::new(),
        shift_sd: vec![0.0; shift.len()],
        shift,
        level: Vec::new(),
        gene_block: Vec::new(),
        cell_block: Vec::new(),
        site_on: Vec::new(),
        site_off: Vec::new(),
    }
}

/// `rows` random membership rows over `k` roles, each summing to 1.
fn random_rows(rows: usize, k: usize, rng: &mut SmallRng) -> Vec<f32> {
    (0..rows)
        .flat_map(|_| {
            let w: Vec<f32> = (0..k)
                .map(|_| rng.random_range(0.01..1.0f32).powi(3))
                .collect();
            let total: f32 = w.iter().sum();
            w.into_iter().map(move |x| x / total)
        })
        .collect()
}

/// `w`'s columns (`k` per row) renumbered: column `j` moved to `to[j]`.
fn renumber<T: Copy + Default>(w: &[T], to: &[usize]) -> Vec<T> {
    let k = to.len();
    let mut out = vec![T::default(); w.len()];
    for (o, row) in out.chunks_mut(k).zip(w.chunks(k)) {
        for (j, &x) in row.iter().enumerate() {
            o[to[j]] = x;
        }
    }
    out
}

/// A relabelling: cell role `a` to `cell[a]`, gene role `b` to `gene[b]`.
fn relabelling(cell: &[usize], gene: &[usize]) -> Alignment {
    Alignment {
        cell: cell.to_vec(),
        gene: gene.to_vec(),
        cell_quality: Vec::new(),
        gene_quality: Vec::new(),
        cell_used: Vec::new(),
        gene_used: Vec::new(),
    }
}

/// `m` relabelled by `to`: its cell and gene memberships and its blocks.
fn relabel(m: &Means, to: &Alignment) -> Means {
    let blocks: Vec<usize> = (0..to.cell.len() * to.gene.len())
        .map(|k| to.block(k))
        .collect();
    memberships(
        renumber(&m.cell, &to.cell),
        renumber(&m.gene, &to.gene),
        renumber(&m.shift, &blocks),
    )
}

/// Every field of `m` relabelled by `to`, by the test's own renumbering:
/// the cell roles' columns (memberships, global weights), the gene roles',
/// and each block's (shifts, their sds, links and reads).
fn relabel_all(m: &Means, to: &Alignment) -> Means {
    let blocks: Vec<usize> = (0..to.cell.len() * to.gene.len())
        .map(|k| to.block(k))
        .collect();
    Means {
        cell: renumber(&m.cell, &to.cell),
        gene: renumber(&m.gene, &to.gene),
        beta: renumber(&m.beta, &to.cell),
        shift: renumber(&m.shift, &blocks),
        shift_sd: renumber(&m.shift_sd, &blocks),
        level: m.level.clone(),
        gene_block: renumber(&m.gene_block, &blocks),
        cell_block: renumber(&m.cell_block, &blocks),
        site_on: renumber(&m.site_on, &blocks),
        site_off: renumber(&m.site_off, &blocks),
    }
}

/// The inverse of permutation `p`.
fn inverse(p: &[usize]) -> Vec<usize> {
    let mut q = vec![0; p.len()];
    p.iter().enumerate().for_each(|(i, &j)| q[j] = i);
    q
}

/// Three cell roles, two gene roles, random memberships and shifts.
fn reference(seed: u64) -> Means {
    let mut rng = SmallRng::seed_from_u64(seed);
    memberships(
        random_rows(40, 3, &mut rng),
        random_rows(20, 2, &mut rng),
        (0..6).map(|_| rng.random_range(-2.0..2.0)).collect(),
    )
}

/// The optimal assignment is the best of every permutation's total.
#[test]
fn the_assignment_is_the_best_permutation() {
    let mut rng = SmallRng::seed_from_u64(301);
    for n in 1..=5 {
        for _ in 0..20 {
            let score: Vec<f64> = (0..n * n).map(|_| rng.random_range(-1.0..1.0)).collect();
            let total = |p: &[usize]| (0..n).map(|i| score[i * n + p[i]]).sum::<f64>();
            let p = assign(&score, n);
            let mut seen = p.clone();
            seen.sort();
            assert_eq!(seen, (0..n).collect::<Vec<_>>(), "one-to-one: {p:?}");
            let best = permutations(n)
                .iter()
                .map(|q| total(q))
                .fold(f64::NEG_INFINITY, f64::max);
            assert!((total(&p) - best).abs() < 1e-12, "{} vs {best}", total(&p));
        }
    }
}

/// Every permutation of `0..n`.
fn permutations(n: usize) -> Vec<Vec<usize>> {
    if n == 0 {
        return vec![Vec::new()];
    }
    permutations(n - 1)
        .into_iter()
        .flat_map(|p| {
            (0..n).map(move |at| {
                let mut q = p.clone();
                q.insert(at, n - 1);
                q
            })
        })
        .collect()
}

/// A chain aligned to itself is the identity, every role of quality 1.
#[test]
fn identical_chains_align_as_the_identity() {
    let m = reference(311);
    let a = align(&m, &m, 3, 2);
    assert_eq!((a.cell, a.gene), (vec![0, 1, 2], vec![0, 1]));
    for q in a.cell_quality.iter().chain(&a.gene_quality) {
        assert!((q - 1.0).abs() < 1e-12, "quality {q}");
    }
}

/// Relabelled copies of a chain align back to it: every relabelling of its
/// cell roles and gene roles is undone, each role of quality 1.
#[test]
fn relabelled_copies_align_back() {
    let m = reference(313);
    for cell in permutations(3) {
        for gene in permutations(2) {
            let copy = relabel(&m, &relabelling(&cell, &gene));
            let a = align(&m, &copy, 3, 2);
            assert_eq!(a.cell, inverse(&cell), "cell roles {cell:?}");
            assert_eq!(a.gene, inverse(&gene), "gene roles {gene:?}");
            assert!(a.cell_quality.iter().all(|q| (q - 1.0).abs() < 1e-12));
            assert!(a.gene_quality.iter().all(|q| (q - 1.0).abs() < 1e-12));
        }
    }
}

/// A chain with two cell roles swapped and its memberships perturbed is
/// put back, its quality below 1 but well above a mismatched role's.
#[test]
fn a_chain_with_two_roles_swapped_is_fixed() {
    let m = reference(317);
    let mut copy = relabel(&m, &relabelling(&[1, 0, 2], &[0, 1]));
    let mut rng = SmallRng::seed_from_u64(319);
    copy.cell
        .iter_mut()
        .for_each(|w| *w *= rng.random_range(0.8..1.2));
    let a = align(&m, &copy, 3, 2);
    assert_eq!(a.cell, vec![1, 0, 2]);
    assert_eq!(a.gene, vec![0, 1]);
    assert!(
        a.cell_quality.iter().all(|&q| q > 0.7 && q < 1.0),
        "{:?}",
        a.cell_quality
    );
}

/// Gene roles whose genes cannot tell them apart are matched by their
/// blocks' shifts.
#[test]
fn block_shifts_break_ties_between_gene_roles() {
    let mut rng = SmallRng::seed_from_u64(323);
    let m = memberships(
        random_rows(40, 3, &mut rng),
        vec![0.5; 40],
        vec![1.0, -1.0, 0.5, -0.5, 2.0, -2.0],
    );
    let copy = relabel(&m, &relabelling(&[0, 1, 2], &[1, 0]));
    assert_eq!(align(&m, &copy, 3, 2).gene, vec![1, 0]);
}

/// Scores that are not finite count as the worst: the assignment is still
/// one-to-one, and optimal over the finite scores.
#[test]
fn scores_not_finite_count_as_the_worst() {
    let nan = f64::NAN;
    let score = [1.0, nan, 0.0, nan, 1.0, f64::NEG_INFINITY, 0.0, 0.0, 1.0];
    assert_eq!(assign(&score, 3), vec![0, 1, 2]);
    let mut seen = assign(&[nan; 9], 3);
    seen.sort();
    assert_eq!(seen, vec![0, 1, 2]);
}

/// The gene memberships' overlap decides the gene roles; the shifts only
/// break its ties, however far apart they are.
#[test]
fn shifts_never_overrule_a_gene_overlap_that_differs() {
    let mut rng = SmallRng::seed_from_u64(327);
    // Two gene roles told apart by a hair, eight cell roles whose shifts
    // favour swapping them.
    let d = 1e-6f32;
    let gene: Vec<f32> = (0..20).flat_map(|_| [0.5 + d, 0.5 - d]).collect();
    let cell = random_rows(40, 8, &mut rng);
    let shift = |b0: f64| (0..8).flat_map(|_| [b0, -b0]).collect::<Vec<f64>>();
    let m = memberships(cell.clone(), gene.clone(), shift(2.0));
    let chain = memberships(cell, gene, shift(-2.0));
    assert_eq!(align(&m, &chain, 8, 2).gene, vec![0, 1]);
}

/// A small chain under `settings` on the tiny fixture: the model's
/// fixture held by the caller.
fn tiny_model(links: &Links) -> Model<'_> {
    Model::new(
        links,
        &TINY_GENES,
        3,
        settings(Weights::Symmetric, true, f64::INFINITY, 2),
    )
    .unwrap()
}

/// Two means the same, to rounding: every field the tables read.
fn assert_same_means(x: &Means, y: &Means) {
    let close = |a: &[f64], b: &[f64], what: &str| {
        assert_eq!(a.len(), b.len(), "{what}");
        for (u, v) in a.iter().zip(b) {
            assert!(
                (u - v).abs() <= 1e-9 * (1.0 + v.abs()),
                "{what}: {u} vs {v}"
            );
        }
    };
    let wide = |w: &[f32]| w.iter().map(|&x| f64::from(x)).collect::<Vec<_>>();
    close(&wide(&x.cell), &wide(&y.cell), "cell");
    close(&wide(&x.gene), &wide(&y.gene), "gene");
    close(&x.beta, &y.beta, "beta");
    close(&x.shift, &y.shift, "shift");
    // The sd from pooled sums and squares, to their rounding.
    for (u, v) in x.shift_sd.iter().zip(&y.shift_sd) {
        assert!((u - v).abs() < 1e-6, "shift sd: {u} vs {v}");
    }
    close(&x.level, &y.level, "level");
    close(&x.gene_block, &y.gene_block, "gene block");
    close(&x.cell_block, &y.cell_block, "cell block");
    close(&x.site_on, &y.site_on, "site on");
    close(&x.site_off, &y.site_off, "site off");
}

/// Copies of one chain pool to that chain's own means, every copy pooled.
#[test]
fn identical_chains_pool_to_one_chains_means() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let cfg = config(200, 50, 331);
    let chains = vec![run_mcmc(&model, &cfg), run_mcmc(&model, &cfg)];
    let (pooled, pooling) = unified::pool(&model, &chains, 0, &[]);
    assert_eq!(pooling.fallback, None);
    assert_eq!(pooling.pooled, vec![0, 1]);
    assert!(pooling.excluded.is_empty());
    assert_same_means(&pooled, &chains[0].last().unwrap().means(true));
}

/// A chain and a relabelled copy of it pool, after alignment, to the
/// chain's own means.
#[test]
fn relabelled_copies_pool_to_one_chains_means() {
    let (reads, gene_of, _, _, _) = planted_blocks(90, 12, 0.5, 337);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let model = Model::new(&links, &gene_of, gene_of.len(), blocks_preset()).unwrap();
    let chain = run_mcmc(&model, &config(40, 100, 339));
    let to = relabelling(&[1, 2, 0], &[1, 0]);
    let copy = Draw::relabelled(&chain, &to);
    let one = chain.last().unwrap().means(true);
    // The copy is relabelled, every field: its cells' first role is the
    // chain's last.
    let moved = copy.last().unwrap().means(true);
    assert_same_means(&moved, &relabel_all(&one, &to));
    assert_ne!(moved.cell[0], one.cell[0]);
    let (pooled, pooling) = unified::pool(&model, &[chain, copy], 0, &[]);
    assert_eq!(pooling.fallback, None, "{:?}", pooling.excluded);
    assert_eq!(pooling.pooled, vec![0, 1]);
    assert_eq!(pooling.alignment[1].cell, inverse(&to.cell));
    assert_eq!(pooling.alignment[1].gene, inverse(&to.gene));
    assert_same_means(&pooled, &one);
}

/// A chain whose log joint the verdict finds lagging the best one's is
/// left out of the pool and named with its reason; the rest pool.
#[test]
fn a_lagging_chain_is_left_out_and_named() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let cfg = config(200, 50, 341);
    let chain = run_mcmc(&model, &cfg);
    let lagging = offset(&chain, -1e3);
    let chains = vec![chain.clone(), lagging, chain];
    let (pooled, pooling) = unified::pool(&model, &chains, 0, &[1]);
    assert_eq!(pooling.fallback, None);
    assert_eq!(pooling.pooled, vec![0, 2]);
    assert_eq!(pooling.excluded.len(), 1);
    assert_eq!(pooling.excluded[0].chain, 1);
    assert!(
        pooling.excluded[0].reason.contains("lags"),
        "{}",
        pooling.excluded[0].reason
    );
    assert!(pooling.rhat_log_joint <= RHAT_OK);
    assert_same_means(&pooled, &chains[0].last().unwrap().means(true));
}

/// With no other chain to pool, the tables are the best chain's alone,
/// and the record says why: every other chain lags, or there is one.
#[test]
fn with_nothing_to_pool_the_best_chain_stands_alone() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chain = run_mcmc(&model, &config(200, 50, 347));
    let one = chain.last().unwrap().means(true);
    for chains in [
        vec![offset(&chain, -1e3), chain.clone()],
        vec![chain.clone()],
    ] {
        let best = chains.len() - 1;
        let (pooled, pooling) = unified::pool(&model, &chains, best, &[]);
        assert_eq!(pooling.pooled, vec![best]);
        assert!(pooling.fallback.is_some());
        assert_eq!(pooling.excluded.len(), chains.len() - 1);
        assert_same_means(&pooled, &one);
    }
}

/// The diagnostics name every quantity short of convergence: with fewer
/// than 100 effective draws across the chains pooled, each quantity falls
/// short, so each is named, and the worst values are the extremes. Here
/// the cells all take the first cell role: the second, empty, is not
/// checked.
#[test]
fn quantities_short_of_convergence_are_named() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let cfg = config(30, 50, 353);
    let chains = vec![run_mcmc(&model, &cfg), run_mcmc(&model, &cfg)];
    let (_, pooling) = unified::pool(&model, &chains, 0, &[]);
    assert_eq!(pooling.alignment[0].cell_used, vec![true, false]);
    assert_eq!(
        pooling.not_converged,
        ["log joint", "shift (0, 0)", "shift (0, 1)", "cell role 0"]
    );
    assert!(pooling.not_converged.iter().any(|q| q == "log joint"));
    assert_eq!(pooling.unsettled, pooling.not_converged, "nothing left out");
    assert!(pooling.worst_ess < 100.0);
    assert!(pooling.worst_rhat >= pooling.rhat_log_joint);
}

/// The effective draws are counted across the chains pooled: four chains
/// of 60 draws reach 100 between them, though none does alone.
#[test]
fn effective_draws_are_counted_across_the_chains_pooled() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chains: Vec<Vec<Draw>> = (0..4)
        .map(|c| run_mcmc(&model, &config(60, 50, 357 + c)))
        .collect();
    let (_, pooling) = unified::pool(&model, &chains, 0, &[]);
    assert_eq!(pooling.pooled, vec![0, 1, 2, 3], "{:?}", pooling.excluded);
    assert!(
        !pooling.not_converged.iter().any(|q| q == "log joint"),
        "{:?}",
        pooling.not_converged
    );
    assert!(pooling.worst_ess > 60.0, "{}", pooling.worst_ess);
}

/// A chain whose log joint matches the best one's but whose roles do not
/// (two of its cell roles merged) is left out of the pool and named.
#[test]
fn a_chain_whose_roles_do_not_match_is_left_out() {
    let (reads, gene_of, _, _, _) = planted_blocks(90, 12, 0.5, 359);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let model = Model::new(&links, &gene_of, gene_of.len(), blocks_preset()).unwrap();
    let chain = run_mcmc(&model, &config(40, 100, 361));
    let merged = Draw::relabelled(&chain, &relabelling(&[0, 0, 2], &[0, 1]));
    let one = chain.last().unwrap().means(true);
    let (pooled, pooling) = unified::pool(&model, &[chain, merged], 0, &[]);
    assert_eq!(pooling.pooled, vec![0]);
    assert_eq!(pooling.excluded.len(), 1);
    assert_eq!(pooling.excluded[0].chain, 1);
    assert!(
        pooling.excluded[0].reason.contains("roles"),
        "{}",
        pooling.excluded[0].reason
    );
    assert!(
        pooling.unsettled.iter().any(|q| q == "chain 1 left out"),
        "{:?}",
        pooling.unsettled
    );
    assert_same_means(&pooled, &one);
}

/// A chain the convergence verdict finds lagging on its trace is left out,
/// the verdict named, though its log joint mixes with the best chain's or
/// not.
#[test]
fn a_chain_the_verdict_finds_lagging_is_left_out() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chain = run_mcmc(&model, &config(200, 50, 371));
    let chains = vec![chain.clone(), chain.clone(), offset(&chain, -1e3)];
    let (pooled, pooling) = unified::pool(&model, &chains, 0, &[1, 2]);
    assert_eq!(
        pooling.fallback.as_deref(),
        Some("no other chain agrees with the best")
    );
    assert_eq!(pooling.pooled, vec![0]);
    let reasons: Vec<(usize, &str)> = pooling
        .excluded
        .iter()
        .map(|e| (e.chain, e.reason.as_str()))
        .collect();
    assert_eq!(reasons.len(), 2, "{reasons:?}");
    assert_eq!(reasons[0].0, 1);
    assert!(reasons[0].1.contains("verdict"), "{}", reasons[0].1);
    assert_eq!(reasons[1].0, 2);
    assert!(reasons[1].1.contains("verdict"), "{}", reasons[1].1);
    // The verdict names the chains it finds lagging itself.
    let left_out = pooling.unsettled.iter().any(|q| q.starts_with("chain"));
    assert!(!left_out, "{:?}", pooling.unsettled);
    assert_same_means(&pooled, &chain.last().unwrap().means(true));
    // The verdict cannot leave out the best chain itself.
    let (_, pooling) = unified::pool(&model, &chains[..2], 0, &[0]);
    assert_eq!(pooling.pooled, vec![0, 1], "{:?}", pooling.excluded);
}

/// `chain` with its log joint moved by `by`.
fn offset(chain: &[Draw], by: f64) -> Vec<Draw> {
    chain
        .iter()
        .map(|d| {
            let mut d = d.clone();
            d.log_joint += by;
            d
        })
        .collect()
}

/// Chains each close enough to the best one, but not to each other, do not
/// pool: each is named with why, and the record keeps their split R-hat.
#[test]
fn chains_that_disagree_together_are_each_named() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chain = run_mcmc(&model, &config(200, 50, 363));
    let log_joint: Vec<f64> = chain.iter().map(|d| d.log_joint).collect();
    let moved = |by: f64| log_joint.iter().map(|x| x + by).collect::<Vec<f64>>();
    // The smallest offset at which the outer two disagree together.
    let sd = {
        let n = log_joint.len() as f64;
        let mean = log_joint.iter().sum::<f64>() / n;
        (log_joint.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt()
    };
    let by = (1..400)
        .map(|i| f64::from(i) * sd / 100.0)
        .find(|&by| {
            let pair = centred_rhat(&[log_joint.clone(), moved(by)]);
            let three = centred_rhat(&[moved(-by), log_joint.clone(), moved(by)]);
            pair <= RHAT_OK && three > RHAT_OK
        })
        .expect("an offset between");
    let chains = vec![offset(&chain, -by), chain.clone(), offset(&chain, by)];
    let (pooled, pooling) = unified::pool(&model, &chains, 1, &[]);
    assert!(pooling.fallback.is_some());
    assert_eq!(pooling.pooled, vec![1]);
    let named: Vec<usize> = pooling.excluded.iter().map(|e| e.chain).collect();
    assert_eq!(named, vec![0, 2], "{:?}", pooling.excluded);
    assert!(
        pooling.rhat_log_joint > RHAT_OK,
        "{}",
        pooling.rhat_log_joint
    );
    assert_same_means(&pooled, &chain.last().unwrap().means(true));
}

/// The summary the record keeps gives the tables' shifts, their sds and
/// the global weights: the chains' pooled, not the best chain's alone.
#[test]
fn the_summary_gives_the_tables_means() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chains = vec![
        run_mcmc(&model, &config(200, 50, 367)),
        run_mcmc(&model, &config(200, 50, 369)),
    ];
    let summary = unified::summarize(&model, &chains);
    let (summary, means, pooling) = unified::pooled(&model, &chains, summary, &[]);
    assert_eq!(pooling.pooled, vec![0, 1], "{:?}", pooling.excluded);
    assert_eq!(summary.shifts, means.shift);
    assert_eq!(summary.shift_sd, means.shift_sd);
    assert_eq!(summary.global_weights, means.beta);
    let best = chains[summary.best_chain].last().unwrap().means(false);
    assert_ne!(summary.shift_sd, best.shift_sd);
}

/// A chain that mixes with the others, its log joint just apart from the
/// best chain's (their pair's split R-hat past the bound, though not all
/// the chains' together, which is what the verdict reads), is pooled: no
/// chain is left out on its pair's R-hat alone.
#[test]
fn a_chain_is_not_left_out_on_its_pairs_r_hat_alone() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chain = run_mcmc(&model, &config(200, 50, 373));
    let log_joint: Vec<f64> = chain.iter().map(|d| d.log_joint).collect();
    let moved = |by: f64| log_joint.iter().map(|x| x + by).collect::<Vec<f64>>();
    let by = (1..400)
        .map(|i| f64::from(i) / 100.0)
        .find(|&by| {
            let pair = centred_rhat(&[log_joint.clone(), moved(by)]);
            let four = centred_rhat(&[
                log_joint.clone(),
                log_joint.clone(),
                log_joint.clone(),
                moved(by),
            ]);
            pair > RHAT_OK && four <= RHAT_OK
        })
        .expect("an offset between");
    let chains = vec![
        chain.clone(),
        chain.clone(),
        chain.clone(),
        offset(&chain, by),
    ];
    let (_, pooling) = unified::pool(&model, &chains, 0, &[]);
    assert_eq!(pooling.pooled, vec![0, 1, 2, 3], "{:?}", pooling.excluded);
    assert!(pooling.excluded.is_empty(), "{:?}", pooling.excluded);
}

/// With nothing to pool, the diagnostics are the tables' chain's own, not
/// those of the chains compared and left out.
#[test]
fn the_diagnostics_are_over_the_chains_pooled() {
    let links = tiny_links();
    let model = tiny_model(&links);
    let chain = run_mcmc(&model, &config(200, 50, 379));
    let chains = vec![offset(&chain, -50.0), chain.clone(), offset(&chain, 50.0)];
    let (_, pooling) = unified::pool(&model, &chains, 1, &[]);
    assert_eq!(pooling.pooled, vec![1], "{:?}", pooling.excluded);
    let (_, alone) = unified::pool(&model, &[chain], 0, &[]);
    assert_eq!(pooling.worst_rhat, alone.worst_rhat);
    assert_eq!(pooling.worst_ess, alone.worst_ess);
    assert_eq!(pooling.not_converged, alone.not_converged);
}

/// Roles no cell or gene has as its main role are not held to a match:
/// under global weights, an empty role's memberships are its weight's
/// smoothing alone, which differs from chain to chain.
#[test]
fn empty_roles_are_not_held_to_a_match() {
    let mut rng = SmallRng::seed_from_u64(383);
    // Two cell roles in use, a third every cell holds at `empty`.
    let cells = |empty: f32| -> Vec<f32> {
        random_rows(40, 2, &mut SmallRng::seed_from_u64(389))
            .chunks(2)
            .flat_map(|r| [r[0] * (1.0 - empty), r[1] * (1.0 - empty), empty])
            .collect()
    };
    let gene = random_rows(20, 2, &mut rng);
    let shift: Vec<f64> = (0..6).map(|_| rng.random_range(-2.0..2.0)).collect();
    let m = memberships(cells(1e-4), gene.clone(), shift.clone());
    let chain = memberships(cells(4e-4), gene, shift);
    let a = align(&m, &chain, 3, 2);
    assert_eq!(a.cell, vec![0, 1, 2]);
    assert!(a.cell_quality[2] < 0.5, "{:?}", a.cell_quality);
    assert_eq!(a.cell_used, vec![true, true, false]);
    assert!(a.worst() > 0.9, "{:?}", a.cell_quality);
}
