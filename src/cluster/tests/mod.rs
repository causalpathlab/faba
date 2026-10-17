use super::crp::{self, ln_beta, ln_pred, Chain, Model, Sample, Settings};
use super::fit::{self, OUTSIDE};
use super::read::{self, Cell, Reads, Site, SiteData, Weighted};
use super::summary::{self, ari};
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use rustc_hash::FxHashMap;
use statrs::function::gamma::ln_gamma;

mod command;

fn settings(gibbs: bool, split_merge: usize, fixed_alpha: Option<f64>) -> Settings {
    Settings {
        strength: 2.0,
        split_merge,
        restricted: 2,
        start_clusters: 3,
        gibbs,
        fixed_alpha,
    }
}

fn config(n_samples: usize, warmup: usize, seed: u64) -> McmcConfig {
    McmcConfig {
        n_samples,
        warmup,
        thin: 1,
        seed,
    }
}

fn site(p: i64) -> Site {
    Site {
        gene: "GENE1".into(),
        chr: "chr1".into(),
        pos: p,
    }
}

fn cell(barcode: &str, fit: bool) -> Cell {
    Cell {
        barcode: barcode.into(),
        batch: "B1".into(),
        fit,
    }
}

/// Cells over two sites, for the exact posteriors: the first ones mostly
/// converted at site 0, the next ones mostly not.
fn small_cells(n: usize) -> Reads {
    let all = [
        vec![(0, 3, 0), (1, 1, 1)],
        vec![(0, 2, 1), (1, 0, 2)],
        vec![(0, 0, 3), (1, 2, 0)],
        vec![(0, 1, 2), (1, 1, 1)],
        vec![(0, 2, 0), (1, 0, 1)],
        vec![(0, 0, 2), (1, 2, 1)],
    ];
    Reads::from_cells(all[..n].to_vec())
}

/// Cells drawn from a fraction per site, a few reads at `cover` of the sites.
fn draw_cells(fracs: &[&[f64]], cover: f64, rng: &mut SmallRng) -> Reads {
    Reads::from_cells(
        fracs
            .iter()
            .map(|frac| {
                let mut entries = Vec::new();
                for (s, &p) in frac.iter().enumerate() {
                    if rng.random::<f64>() < cover {
                        let n = rng.random_range(1..4u32);
                        let m = (0..n).filter(|_| rng.random::<f64>() < p).count() as u32;
                        entries.push((s as u32, m, n - m));
                    }
                }
                entries
            })
            .collect(),
    )
}

#[test]
fn predictive_as_logs_matches_beta_functions() {
    for &(x, y, m, u) in &[
        (0.7, 2.3, 1, 0),
        (3.0, 1.5, 2, 3),
        (10.2, 4.4, 9, 8),
        (0.4, 0.6, 20, 1),
    ] {
        let exact = ln_beta(x + m as f64, y + u as f64) - ln_beta(x, y);
        assert!(
            (ln_pred(x, y, m as f64, u as f64) - exact).abs() < 1e-9,
            "{x} {y} {m} {u}"
        );
    }
    // Weighted reads are not whole: the Beta-function path.
    let exact = ln_beta(1.3 + 0.7, 2.0 + 1.4) - ln_beta(1.3, 2.0);
    assert!((ln_pred(1.3, 2.0, 0.7, 1.4) - exact).abs() < 1e-12);
}

/// Every partition of `0..n`, as labels by first appearance.
fn partitions(n: usize) -> Vec<Vec<u16>> {
    let mut out = vec![vec![0u16]];
    for _ in 1..n {
        out = out
            .into_iter()
            .flat_map(|p| {
                let k = *p.iter().max().unwrap() + 1;
                (0..=k).map(move |l| {
                    let mut q = p.clone();
                    q.push(l);
                    q
                })
            })
            .collect();
    }
    out
}

/// The exact posterior over every partition of `reads`' cells, `alpha`
/// fixed or, `None`, integrated out under its Gamma(1, 1) prior: the CRP
/// times each cluster's Beta-binomial marginal at the model's prior.
fn exact_posterior(reads: &Reads, strength: f64, alpha: Option<f64>) -> Vec<(Vec<u16>, f64)> {
    let (n, sites) = (
        reads.cells(),
        *reads.site.iter().max().unwrap() as usize + 1,
    );
    let (mut on, mut tot) = (vec![0.0f64; sites], vec![0.0f64; sites]);
    for i in 0..reads.site.len() {
        on[reads.site[i] as usize] += reads.on[i] as f64;
        tot[reads.site[i] as usize] += (reads.on[i] + reads.off[i]) as f64;
    }
    let p: Vec<f64> = (0..sites)
        .map(|s| crp::pooled_fraction(on[s], tot[s]))
        .collect();
    let parts = partitions(n);
    let logp: Vec<f64> = parts
        .iter()
        .map(|part| {
            let k = *part.iter().max().unwrap() as usize + 1;
            let mut v = 0.0;
            for l in 0..k {
                let members: Vec<usize> = (0..n).filter(|&c| part[c] as usize == l).collect();
                v += ln_gamma(members.len() as f64);
                let (mut m, mut u) = (vec![0u32; sites], vec![0u32; sites]);
                for &c in &members {
                    let (ss, cm, cu) = reads.of(c);
                    for i in 0..ss.len() {
                        m[ss[i] as usize] += cm[i];
                        u[ss[i] as usize] += cu[i];
                    }
                }
                for s in 0..sites {
                    let (a, b) = (strength * p[s], strength * (1.0 - p[s]));
                    v += ln_beta(a + m[s] as f64, b + u[s] as f64) - ln_beta(a, b);
                }
            }
            // The CRP's alpha^K Gamma(alpha) / Gamma(alpha + n), at alpha or
            // integrated against exp(-alpha) on a fine grid.
            let crp = |al: f64| k as f64 * al.ln() + ln_gamma(al) - ln_gamma(al + n as f64);
            v + match alpha {
                Some(al) => crp(al),
                None => {
                    let h = 1e-3;
                    let terms: Vec<f64> = (1..60_000)
                        .map(|t| t as f64 * h)
                        .map(|al| crp(al) - al)
                        .collect();
                    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                    mx + (terms.iter().map(|t| (t - mx).exp()).sum::<f64>() * h).ln()
                }
            }
        })
        .collect();
    let mx = logp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let z: f64 = logp.iter().map(|l| (l - mx).exp()).sum();
    parts
        .into_iter()
        .zip(logp.iter().map(|l| (l - mx).exp() / z))
        .collect()
}

/// Total variation between the chains' draw frequencies and `exact`.
fn total_variation(chains: &[Chain], exact: &[(Vec<u16>, f64)]) -> f64 {
    let total = chains.iter().map(|c| c.samples.len()).sum::<usize>() as f64;
    let mut freq: FxHashMap<&[u16], f64> = FxHashMap::default();
    for s in chains.iter().flat_map(|c| &c.samples) {
        *freq.entry(&s.z).or_default() += 1.0 / total;
    }
    exact
        .iter()
        .map(|(p, e)| (freq.get(p.as_slice()).copied().unwrap_or(0.0) - e).abs())
        .sum::<f64>()
        / 2.0
}

/// Exact posterior of every partition of four cells, alpha fixed, against
/// the sampler's draw frequencies: with split-merge alone and with both moves.
#[test]
fn sampler_matches_the_exact_posterior_on_four_cells() {
    let reads = small_cells(4);
    let weighted = Weighted::new(&reads, None);
    let alpha = 0.8;
    let exact = exact_posterior(&reads, 2.0, Some(alpha));
    for (gibbs, sm) in [(false, 3), (true, 1)] {
        let model = Model::new(&weighted, 2, settings(gibbs, sm, Some(alpha)));
        let chains = run_mcmc_parallel(&model, &config(20_000, 200, 7), 2);
        let tv = total_variation(&chains, &exact);
        assert!(
            tv < 0.02,
            "gibbs {gibbs}: total variation {tv:.4} from the exact posterior"
        );
    }
}

/// Every move together, alpha resampled, against the exact posterior over
/// all 203 partitions of six cells: the chain leaves the posterior invariant
/// (detailed balance of each move) and reaches it.
#[test]
fn all_moves_with_alpha_resampled_match_the_exact_posterior() {
    let reads = small_cells(6);
    let weighted = Weighted::new(&reads, None);
    let model = Model::new(&weighted, 2, settings(true, 2, None));
    let chains = run_mcmc_parallel(&model, &config(40_000, 500, 21), 2);
    let tv = total_variation(&chains, &exact_posterior(&reads, 2.0, None));
    assert!(
        tv < 0.03,
        "total variation {tv:.4} from the exact posterior"
    );
}

/// Planted clusters at the sparsity of real cells (a few reads at a few
/// percent of the sites) are recovered, with their number.
#[test]
fn recovers_planted_clusters() {
    let (cells, sites, k) = (300usize, 400usize, 3usize);
    let mut rng = SmallRng::seed_from_u64(3);
    let base: Vec<f64> = (0..sites).map(|_| rng.random_range(0.05..0.4)).collect();
    // A third of the sites differ by cluster.
    let frac: Vec<Vec<f64>> = (0..k)
        .map(|_| {
            base.iter()
                .map(|&p| {
                    if rng.random::<f64>() < 0.33 {
                        rng.random_range(0.0..0.9)
                    } else {
                        p
                    }
                })
                .collect()
        })
        .collect();
    let truth: Vec<u16> = (0..cells).map(|c| (c % k) as u16).collect();
    let fracs: Vec<&[f64]> = truth.iter().map(|&t| frac[t as usize].as_slice()).collect();
    let reads = draw_cells(&fracs, 0.15, &mut rng);
    let weighted = Weighted::new(&reads, None);
    let model = Model::new(&weighted, sites, settings(true, 5, None));
    let chains: Vec<Chain> = run_mcmc_parallel(&model, &config(60, 60, 11), 2);
    // A side group of ten-odd cells is all but tied with folding it back
    // (the posterior barely tells them apart), so count clusters of 20 up.
    let (part, _) = fit::settle(
        &reads,
        &weighted,
        &summary::partition(&chains),
        20,
        sites,
        2.0,
    );
    let post = summary::posterior(&chains, 20);
    let mode = post.k.iter().max_by(|a, b| a.1.total_cmp(&b.1)).unwrap().0;
    assert_eq!(mode, k, "K posterior {:?}", post.k);
    assert_eq!(part.k, k);
    let score = ari(&part.z, &truth);
    assert!(score > 0.95, "ARI {score:.3}");
    // Draws park a sparse cell in a passing singleton now and then, which
    // lowers its co-clustering; most cells stay put.
    let mut conf = part.confidence.clone();
    conf.sort_by(f64::total_cmp);
    let median = conf[cells / 2];
    assert!(median > 0.9, "median confidence {median:.3}");
}

/// The contingency-table summaries equal those from an explicit similarity matrix.
#[test]
fn summaries_match_an_explicit_similarity_matrix() {
    let mut rng = SmallRng::seed_from_u64(5);
    let n = 40;
    let mut draw = || {
        let k = rng.random_range(1..5u32);
        let raw: Vec<u32> = (0..n).map(|_| rng.random_range(0..k)).collect();
        let (z, k) = crp::compact(&raw);
        Sample {
            k,
            z,
            alpha: 1.0,
            log_post: 0.0,
            heldout: (0.0, 0.0),
        }
    };
    let chains: Vec<Chain> = (0..2)
        .map(|_| Chain {
            samples: (0..15).map(|_| draw()).collect(),
            splits: 0,
            merges: 0,
        })
        .collect();
    let draws: Vec<&Sample> = chains.iter().flat_map(|c| &c.samples).collect();
    let t = draws.len() as f64;
    let mut psm = vec![0.0f64; n * n];
    for d in &draws {
        for i in 0..n {
            for j in 0..n {
                if d.z[i] == d.z[j] {
                    psm[i * n + j] += 1.0 / t;
                }
            }
        }
    }
    let loss = |z: &[u16]| -> f64 {
        (0..n * n)
            .map(|ij| {
                let same = f64::from(u8::from(z[ij / n] == z[ij % n]));
                (same - psm[ij]).powi(2)
            })
            .sum()
    };
    let best = draws
        .iter()
        .map(|d| loss(&d.z))
        .fold(f64::INFINITY, f64::min);
    let part = summary::partition(&chains);
    assert!(
        (loss(&part.z) - best).abs() < 1e-9,
        "not the least-squares draw"
    );
    for i in 0..n {
        let own: Vec<usize> = (0..n).filter(|&j| part.z[j] == part.z[i]).collect();
        let explicit = own.iter().map(|&j| psm[i * n + j]).sum::<f64>() / own.len() as f64;
        assert!((part.confidence[i] - explicit).abs() < 1e-9, "cell {i}");
    }
}

#[test]
fn ari_is_one_for_relabelled_partitions_and_near_zero_for_unrelated() {
    let a = [0u16, 0, 1, 1, 2, 2];
    assert!((ari(&a, &[1, 1, 2, 2, 0, 0]) - 1.0).abs() < 1e-12);
    assert!(ari(&a, &[0, 1, 2, 0, 1, 2]) < 0.1);
}

#[test]
fn filter_drops_thin_cells_then_thin_sites_then_cells_left_thin() {
    let data = SiteData {
        modality: "m6a".into(),
        sites: (0..3).map(site).collect(),
        cells: vec![
            cell("A", true),
            cell("B", true),
            cell("C", true),
            cell("D", false),
        ],
        reads: Reads::from_cells(vec![
            vec![(0, 1, 0), (1, 1, 1), (2, 0, 1)],
            vec![(0, 2, 0), (1, 0, 1)],
            vec![(0, 1, 0)],
            vec![(0, 1, 0), (2, 3, 0)],
        ]),
    };
    // Cells need 2 sites: C goes. Sites need 2 fitted cells: site 2 (A only;
    // D is not fitted) goes, which leaves D one site, so D goes too.
    let (out, kept) = read::filter(&data, 2, 2);
    assert_eq!(kept, vec![0, 1]);
    assert_eq!(
        out.sites.iter().map(|s| s.pos).collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(out.reads.of(1).0, &[0, 1]);
}

/// A cell that passes the cell filter but loses its sites to the site filter
/// is dropped too, not kept with nothing to place it by.
#[test]
fn a_cell_left_without_sites_is_dropped() {
    let data = SiteData {
        modality: "m6a".into(),
        sites: (0..4).map(site).collect(),
        cells: vec![cell("A", true), cell("B", true), cell("C", false)],
        reads: Reads::from_cells(vec![
            vec![(0, 1, 0), (1, 1, 1)],
            vec![(0, 2, 0), (1, 0, 1)],
            // Sites 2 and 3 are covered by no fitted cell.
            vec![(2, 1, 0), (3, 3, 0)],
        ]),
    };
    let (out, kept) = read::filter(&data, 2, 2);
    assert_eq!(kept, vec![0, 1]);
    assert_eq!((out.cells.len(), out.sites.len()), (2, 2));
}

#[test]
fn a_sites_two_channels_merge_into_one_entry() {
    let reads = Reads::from_cells(vec![vec![(1, 2, 0), (0, 1, 0), (1, 0, 3)], vec![]]);
    assert_eq!(reads.of(0), (&[0, 1][..], &[1, 2][..], &[0, 3][..]));
    assert_eq!(reads.of(1).0.len(), 0);
    // A subset keeps each cell's entries, in the order asked.
    let sub = reads.subset(&[1, 0]);
    assert_eq!((sub.of(0).0.len(), sub.of(1)), (0, reads.of(0)));
}

/// Cells drawn alike in two batches, one batch shifted at every site: the
/// shift alone makes two clusters; weighted out, one.
#[test]
fn batch_shifts_are_cancelled() {
    let (cells, sites) = (240usize, 150usize);
    let mut rng = SmallRng::seed_from_u64(9);
    let base: Vec<f64> = (0..sites).map(|_| rng.random_range(0.1..0.5)).collect();
    // Batch 1 reads more converted everywhere: logit + 1.5.
    let shifted: Vec<f64> = base
        .iter()
        .map(|&p| 1.0 / (1.0 + (-((p / (1.0 - p)).ln() + 1.5)).exp()))
        .collect();
    let batch_of: Vec<u16> = (0..cells).map(|c| (c % 2) as u16).collect();
    let fracs: Vec<&[f64]> = batch_of
        .iter()
        .map(|&b| {
            if b == 1 {
                shifted.as_slice()
            } else {
                base.as_slice()
            }
        })
        .collect();
    // Enough reads per batch and site (about 140) that the shrinkage of each
    // batch's fraction toward the site's leaves little of the shift; thinly
    // read sites keep part of it by design.
    let reads = draw_cells(&fracs, 0.6, &mut rng);
    let clusters = |weighted: &Weighted| {
        let model = Model::new(weighted, sites, settings(true, 5, None));
        let chains = run_mcmc_parallel(&model, &config(20, 20, 13), 2);
        fit::settle(
            &reads,
            weighted,
            &summary::partition(&chains),
            10,
            sites,
            2.0,
        )
        .0
    };
    let raw = clusters(&Weighted::new(&reads, None));
    let score = ari(&raw.z, &batch_of);
    assert!(
        raw.k >= 2 && score > 0.8,
        "uncorrected: K {} ARI {score:.2}",
        raw.k
    );
    let w = fit::batch_weights(&reads, &batch_of, 2, sites);
    let fixed = clusters(&Weighted::new(&reads, Some((&w.on, &w.off))));
    assert_eq!(fixed.k, 1, "corrected: K {}", fixed.k);
}

/// Weighted, each batch's reads at a site keep their total and take the
/// site's pooled fraction (up to the shrinkage of thinly read batches).
#[test]
fn batch_weights_keep_totals_and_match_the_pooled_fraction() {
    // One site; batch 0: 80 of 100 converted; batch 1: 20 of 100.
    let batch_of: Vec<u16> = (0..20).map(|c| (c / 10) as u16).collect();
    let reads = Reads::from_cells(
        batch_of
            .iter()
            .map(|&b| {
                let m = if b == 0 { 8 } else { 2 };
                vec![(0u32, m, 10 - m)]
            })
            .collect(),
    );
    let w = fit::batch_weights(&reads, &batch_of, 2, 1);
    let weighted = Weighted::new(&reads, Some((&w.on, &w.off)));
    let p = crp::pooled_fraction(100.0, 200.0);
    for b in 0..2u16 {
        let (mut won, mut wtot) = (0.0, 0.0);
        for c in (0..20).filter(|&c| batch_of[c] == b) {
            let (_, m, u) = weighted.of(c);
            won += m[0];
            wtot += m[0] + u[0];
        }
        assert!((wtot - 100.0).abs() < 1e-9, "batch {b} total {wtot}");
        // The batch's fraction, shrunk by 20 pseudo-reads toward p, is what
        // the weights take out: most of the gap to p closes, not all.
        let (m, n) = (if b == 0 { 80.0 } else { 20.0 }, 100.0);
        let pb = (m + 20.0 * p) / (n + 20.0);
        let expect = m * p / pb / (m * p / pb + (n - m) * (1.0 - p) / (1.0 - pb));
        let got = won / wtot;
        assert!((got - expect).abs() < 1e-9, "batch {b} fraction {got:.4}");
        assert!((got - p).abs() < 0.08 && (m / n - p).abs() > 0.25);
    }
    assert!(
        w.effective_share.iter().all(|&e| e > 0.5 && e <= 1.0),
        "{:?}",
        w.effective_share
    );
    assert_eq!(w.capped_share, 0.0);
}

/// A cell with one cluster's pattern but a lower level everywhere goes to
/// that cluster, not to the one whose sites are all lower.
#[test]
fn placement_reads_pattern_not_level() {
    // Cluster 0: sites 0-1 high, 2-3 low. Cluster 1: everything lowish,
    // flat. Each cluster: 20 cells, and one cell left out of both.
    let profile = [[0.8, 0.8, 0.2, 0.2], [0.15, 0.15, 0.15, 0.15]];
    let mut cells: Vec<Vec<(u32, u32, u32)>> = profile
        .iter()
        .flat_map(|prof| {
            let cell: Vec<_> = prof
                .iter()
                .enumerate()
                .map(|(s, &p)| {
                    let m = (p * 20.0f64).round() as u32;
                    (s as u32, m, 20 - m)
                })
                .collect();
            std::iter::repeat_n(cell, 20)
        })
        .collect();
    cells.push(vec![(0, 9, 11)]);
    let reads = Reads::from_cells(cells);
    let mut z: Vec<u16> = (0..40).map(|c| (c / 20) as u16).collect();
    z.push(OUTSIDE);
    let tiles = fit::Tiles::new(&reads, &Weighted::new(&reads, None), &z, 2, 4, 2.0);
    // Cluster 0's pattern (high, high, low, low), shifted down to 0.4/0.4/0.05/0.05.
    let p = tiles.place(
        &[0, 1, 2, 3],
        &[8.0, 8.0, 1.0, 1.0],
        &[12.0, 12.0, 19.0, 19.0],
        &[20, 20],
    );
    assert!(
        p.cluster == 0 && p.probability > 0.9,
        "{} {}",
        p.cluster,
        p.probability
    );
    // Its offset is the level it sits below the cluster.
    assert!(p.offset < -0.5, "offset {}", p.offset);
    // The cell left out of the tiles counts toward neither cluster.
    assert_eq!((tiles.total(0, 0), tiles.total(0, 1)), (400, 400));
}

/// Folding small clusters into large ones can reorder the sizes: the
/// result is numbered by size again, largest first.
#[test]
fn settled_clusters_are_numbered_by_size() {
    // Clusters 0 (5 cells) and 1 (4 cells) are kept; cluster 2's 3 cells
    // read like cluster 1 and join it, making it the largest.
    let a = vec![(0u32, 10u32, 0u32), (1, 0, 10)];
    let b = vec![(0u32, 0u32, 10u32), (1, 10, 0)];
    let mut cells = vec![a; 5];
    cells.extend(std::iter::repeat_n(b, 7));
    let reads = Reads::from_cells(cells);
    let weighted = Weighted::new(&reads, None);
    let z: Vec<u16> = [0; 5].into_iter().chain([1; 4]).chain([2; 3]).collect();
    let part = summary::Partition {
        confidence: vec![1.0; z.len()],
        k: 3,
        z,
    };
    let (out, _) = fit::settle(&reads, &weighted, &part, 4, 2, 2.0);
    assert_eq!(out.k, 2);
    let first: Vec<u16> = out.z.iter().copied().take(5).collect();
    assert!(first.iter().all(|&x| x == 1), "{:?}", out.z);
    assert!(out.z[5..].iter().all(|&x| x == 0), "{:?}", out.z);
}

#[test]
fn a_batch_name_in_two_folders_is_two_batches() {
    let input = |path: &str, batch: &str| read::Input {
        path: path.into(),
        folder: path.split('/').next().unwrap().into(),
        batch: batch.into(),
        fit: true,
    };
    let mut inputs = vec![
        input("dirA/B1_m6a_site.zarr", "B1"),
        input("dirB/B1_m6a_site.zarr", "B1"),
        input("dirA/B2_m6a_site.zarr", "B2"),
    ];
    super::args::qualify_batches(&mut inputs);
    let names: Vec<&str> = inputs.iter().map(|i| i.batch.as_ref()).collect();
    assert_eq!(names, vec!["dirA/B1", "dirB/B1", "B2"]);
}

/// Too few draws per chain to judge convergence, or a prior with no
/// strength, are refused before anything is read.
#[test]
fn settings_that_cannot_judge_or_fit_are_refused() {
    use clap::{Args, FromArgMatches};
    let parse = |line: &[&str]| {
        let cmd = super::ClusterArgs::augment_args(clap::Command::new("cluster"));
        let m = cmd.try_get_matches_from(line).unwrap();
        super::ClusterArgs::from_arg_matches(&m).unwrap()
    };
    let ok = parse(&["cluster", "x.zarr"]);
    assert!(super::check(&ok).is_ok());
    let few = parse(&[
        "cluster",
        "x.zarr",
        "--sweeps",
        "60",
        "--burn-in",
        "50",
        "--thin",
        "4",
    ]);
    assert!(super::check(&few).is_err());
    let flat = parse(&["cluster", "x.zarr", "--prior-strength", "0"]);
    assert!(super::check(&flat).is_err());
    for (flag, value) in [
        ("--communities", "0"),
        ("--communities", "65535"),
        ("--membership-alpha", "0"),
        ("--community-gamma", "-1"),
        ("--holdout", "1"),
        ("--holdout", "-0.1"),
    ] {
        let arg = format!("{flag}={value}");
        let bad = parse(&["cluster", "x.zarr", &arg]);
        assert!(super::check(&bad).is_err(), "{flag} {value} accepted");
    }
}

// ---------------------------------------------------------------- links

use super::link;

/// Every labelling of two-community `links` with its exact posterior
/// probability: memberships and site fractions integrated out in closed
/// form, the global weight `beta_1 ~ Beta(gamma/2, gamma/2)` numerically.
/// Also the posterior mean of the weight of the first link's community,
/// `E[beta_(z_0)]`, which checks the global weights' updates.
fn exact_link_posterior(
    links: &link::Links,
    sites: usize,
    s: &link::Settings,
) -> (Vec<(Vec<u16>, f64)>, f64) {
    assert_eq!(s.k, 2);
    let model = link::Model::new(links, sites, s.clone());
    let n = links.len();
    let half = s.gamma / 2.0;
    let ln_prior = |b: f64| (half - 1.0) * (b.ln() + (1.0 - b).ln()) - ln_beta(half, half);
    let grid: Vec<f64> = (0..2000).map(|i| (i as f64 + 0.5) / 2000.0).collect();
    let mut out = Vec::new();
    let mut weight_mean = Vec::new();
    for code in 0..1usize << n {
        let z: Vec<u16> = (0..n).map(|l| ((code >> l) & 1) as u16).collect();
        let terms: Vec<f64> = grid
            .iter()
            .map(|&b| model.log_joint_of(&z, &[b, 1.0 - b]) + ln_prior(b))
            .collect();
        let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
        let mass: f64 = w.iter().sum();
        let own = |b: f64| if z[0] == 0 { b } else { 1.0 - b };
        weight_mean.push(grid.iter().zip(&w).map(|(&b, &x)| own(b) * x).sum::<f64>() / mass);
        out.push((z, mx + mass.ln()));
    }
    let mx = out.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = out.iter().map(|x| (x.1 - mx).exp()).sum();
    let probs: Vec<(Vec<u16>, f64)> = out
        .into_iter()
        .map(|(z, l)| (z, (l - mx).exp() / total))
        .collect();
    let mean = probs
        .iter()
        .zip(&weight_mean)
        .map(|((_, p), m)| p * m)
        .sum();
    (probs, mean)
}

/// Two cells, two genes (sites 0-1 and site 2): four links over two
/// communities, every labelling enumerated.
#[test]
fn link_sampler_matches_the_exact_posterior() {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 0, 2), (2, 1, 1)],
        vec![(0, 0, 3), (1, 2, 0), (2, 2, 0)],
    ]);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &[0, 0, 1]);
    assert_eq!(links.len(), 4);
    let s = link::Settings {
        k: 2,
        alpha: 0.7,
        gamma: 2.0,
        strength: 2.0,
    };
    let (exact, own_weight) = exact_link_posterior(&links, 3, &s);
    let model = link::Model::new(&links, 3, s);
    let chains = run_mcmc_parallel(&model, &config(40_000, 500, 61), 2);
    let total = chains.iter().map(Vec::len).sum::<usize>() as f64;
    let mut freq: FxHashMap<&[u16], f64> = FxHashMap::default();
    for d in chains.iter().flatten() {
        *freq.entry(&d.z).or_default() += 1.0 / total;
    }
    let tv: f64 = exact
        .iter()
        .map(|(z, p)| (freq.get(z.as_slice()).copied().unwrap_or(0.0) - p).abs())
        .sum::<f64>()
        / 2.0;
    assert!(
        tv < 0.02,
        "total variation {tv:.4} from the exact posterior"
    );
    let draws: Vec<&link::Draw> = chains.iter().flatten().collect();
    let sampled = draws
        .iter()
        .map(|d| d.beta.get(d.z[0] as usize).copied().unwrap_or(0.0))
        .sum::<f64>()
        / draws.len() as f64;
    assert!(
        (sampled - own_weight).abs() < 0.01,
        "E[beta of link 0's community]: sampled {sampled:.4}, exact {own_weight:.4}"
    );
}

/// The cell swap's Metropolis-Hastings ratio is exactly the change in log
/// joint, for every labelling, cell and pair of communities: checked
/// directly, so an error too small to move the sampled frequencies still
/// shows.
#[test]
fn link_swap_ratio_is_the_change_in_log_joint() {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 0, 2), (2, 1, 1), (3, 2, 2)],
        vec![(0, 0, 3), (1, 2, 0), (2, 2, 0), (3, 0, 1)],
        vec![(0, 1, 1), (2, 3, 1), (3, 1, 0)],
    ]);
    let weighted = Weighted::new(&reads, None);
    let gene_of = [0, 0, 1, 2];
    let links = link::Links::new(&weighted, &gene_of);
    let s = link::Settings {
        k: 3,
        alpha: 0.7,
        gamma: 2.0,
        strength: 2.0,
    };
    let model = link::Model::new(&links, 4, s);
    let beta = [0.5, 0.3, 0.2];
    let n = links.len();
    let mut checked = 0;
    for code in 0..3usize.pow(n as u32) {
        let z: Vec<u16> = (0..n)
            .map(|l| ((code / 3usize.pow(l as u32)) % 3) as u16)
            .collect();
        for c in 0..3 {
            for (k, k2) in [(0, 1), (0, 2), (1, 2)] {
                if let Some((delta, exact)) = model.swap_check(&z, &beta, c, k, k2) {
                    assert!(
                        (delta - exact).abs() < 1e-9,
                        "z {z:?} cell {c} swap {k}<->{k2}: {delta} vs {exact}"
                    );
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > 1000, "{checked} swaps checked");
}

/// Cells of three planted within-gene patterns, genes of four sites each.
fn planted_links(cells: usize, genes: usize, seed: u64) -> (Reads, Vec<u32>, Vec<u16>) {
    let sites = genes * 4;
    let mut rng = SmallRng::seed_from_u64(seed);
    let base: Vec<f64> = (0..sites).map(|_| rng.random_range(0.1..0.4)).collect();
    let frac: Vec<Vec<f64>> = (0..3)
        .map(|_| {
            base.iter()
                .map(|&p| {
                    if rng.random::<f64>() < 0.4 {
                        rng.random_range(0.0..0.9)
                    } else {
                        p
                    }
                })
                .collect()
        })
        .collect();
    let truth: Vec<u16> = (0..cells).map(|c| (c % 3) as u16).collect();
    let fracs: Vec<&[f64]> = truth.iter().map(|&t| frac[t as usize].as_slice()).collect();
    let gene_of: Vec<u32> = (0..sites as u32).map(|s| s / 4).collect();
    (draw_cells(&fracs, 0.3, &mut rng), gene_of, truth)
}

/// Each cell's most probable community: its chain's posterior-mean
/// membership.
fn main_community(d: &link::Draw, k: usize) -> Vec<u16> {
    summary::main_labels(&d.means().membership, k)
}

#[test]
fn links_recover_planted_cell_patterns() {
    let (reads, gene_of, truth) = planted_links(240, 60, 71);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = link::Settings {
        k: 3,
        alpha: 0.5,
        gamma: 1.0,
        strength: 2.0,
    };
    let model = link::Model::new(&links, gene_of.len(), s);
    let chains = run_mcmc_parallel(&model, &config(20, 80, 73), 2);
    for c in &chains {
        let z = main_community(c.last().unwrap(), 3);
        let score = ari(&z, &truth);
        assert!(score > 0.9, "ARI {score:.3}");
    }
}

/// Cells left out of the fit, placed against the fitted communities, fall
/// with the fitted cells of their own pattern.
#[test]
fn control_cells_are_placed_by_their_pattern() {
    let (reads, gene_of, truth) = planted_links(360, 60, 109);
    let fit: Vec<usize> = (0..240).collect();
    let rest: Vec<usize> = (240..360).collect();
    let (fit_reads, rest_reads) = (reads.subset(&fit), reads.subset(&rest));
    let weighted = Weighted::new(&fit_reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = link::Settings {
        k: 3,
        alpha: 0.5,
        gamma: 1.0,
        strength: 2.0,
    };
    let model = link::Model::new(&links, gene_of.len(), s);
    let chains = run_mcmc_parallel(&model, &config(20, 80, 111), 1);
    let last = chains[0].last().unwrap();
    let means = last.means();
    let profile = model.profile(&means.z);
    let rest_weighted = Weighted::new(&rest_reads, None);
    let rest_links = link::Links::new(&rest_weighted, &gene_of);
    let placed = model.place(&rest_links, &profile, &last.beta, 113);
    let fitted = summary::main_labels(&means.membership, 3);
    let placed = summary::main_labels(&placed, 3);
    // Placed cells and fitted cells of one pattern share a community.
    let all: Vec<u16> = fitted.iter().chain(&placed).copied().collect();
    let score = ari(&all, &truth);
    assert!(score > 0.9, "ARI {score:.3} over fitted and placed cells");
}

/// Every cell alike: the reads condense into one community instead of
/// being split by their own outcomes.
#[test]
fn links_condense_without_structure() {
    let (reads, gene_of, _) = planted_links(240, 60, 79);
    // Re-draw every cell from the first pattern only.
    let mut rng = SmallRng::seed_from_u64(80);
    let same: Vec<f64> = (0..gene_of.len())
        .map(|_| rng.random_range(0.1..0.6))
        .collect();
    let fracs: Vec<&[f64]> = (0..reads.cells()).map(|_| same.as_slice()).collect();
    let reads = draw_cells(&fracs, 0.3, &mut rng);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = link::Settings {
        k: 4,
        alpha: 0.5,
        gamma: 1.0,
        strength: 2.0,
    };
    let model = link::Model::new(&links, gene_of.len(), s);
    let chains = run_mcmc_parallel(&model, &config(20, 150, 81), 2);
    for c in &chains {
        let share = link::read_shares(&links, &c.last().unwrap().z, 4);
        let top = share.iter().cloned().fold(0.0, f64::max);
        assert!(top > 0.9, "read shares {share:?}");
    }
}

/// A smaller Dirichlet concentration gives cells more concentrated
/// memberships (the per-cell urn condenses harder).
#[test]
fn smaller_alpha_concentrates_memberships() {
    let (reads, gene_of, _) = planted_links(150, 40, 83);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let top = |alpha: f64| {
        let s = link::Settings {
            k: 3,
            alpha,
            gamma: 1.0,
            strength: 2.0,
        };
        let model = link::Model::new(&links, gene_of.len(), s);
        let chains = run_mcmc_parallel(&model, &config(10, 40, 85), 1);
        let m = &chains[0].last().unwrap().means().membership;
        m.chunks(3)
            .map(|row| row.iter().cloned().fold(0.0f32, f32::max) as f64)
            .sum::<f64>()
            / 150.0
    };
    let (sharp, flat) = (top(0.1), top(5.0));
    assert!(
        sharp > flat + 0.1,
        "mean top membership: alpha 0.1 {sharp:.3}, alpha 5 {flat:.3}"
    );
}

/// The global weights' update, labels held fixed, leaves `p(beta | z)`
/// invariant: draws of `beta_1` match its exact conditional mean
/// (numerical integral) for a large labelling, where a wrong table count
/// would show.
#[test]
fn global_weights_follow_their_conditional() {
    let (reads, gene_of, truth) = planted_links(300, 40, 91);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = link::Settings {
        k: 2,
        alpha: 0.5,
        gamma: 2.0,
        strength: 2.0,
    };
    let model = link::Model::new(&links, gene_of.len(), s);
    // Cells of pattern 0 in community 0, the rest in 1: unequal sizes.
    let z: Vec<u16> = links
        .cell
        .iter()
        .map(|&c| u16::from(truth[c as usize] != 0))
        .collect();
    let grid: Vec<f64> = (1..4000).map(|i| i as f64 / 4000.0).collect();
    let terms: Vec<f64> = grid
        .iter()
        .map(|&b| model.log_joint_of(&z, &[b, 1.0 - b]))
        .collect();
    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
    let mass: f64 = w.iter().sum();
    let exact = grid.iter().zip(&w).map(|(b, x)| b * x).sum::<f64>() / mass;
    let exact_sd = (grid
        .iter()
        .zip(&w)
        .map(|(b, x)| (b - exact).powi(2) * x)
        .sum::<f64>()
        / mass)
        .sqrt();
    let draws: Vec<f64> = model
        .beta_draws(&z, 4000, 93)
        .into_iter()
        .skip(200)
        .collect();
    let sampled = draws.iter().sum::<f64>() / draws.len() as f64;
    let sd = (draws.iter().map(|b| (b - sampled).powi(2)).sum::<f64>() / draws.len() as f64).sqrt();
    // The mean alone misses an update off by a common factor: the spread
    // shows it.
    assert!(
        (sampled - exact).abs() < 0.01 && (sd / exact_sd - 1.0).abs() < 0.15,
        "beta_1 | z: sampled mean {sampled:.4} sd {sd:.4}, exact mean {exact:.4} sd {exact_sd:.4}"
    );
}

/// A cell's reads are grouped by gene into links, and each link knows its
/// gene.
#[test]
fn links_group_a_cells_reads_by_gene() {
    // Sites 0 and 2 in gene 0, site 1 in gene 1.
    let reads = Reads::from_cells(vec![vec![(0, 1, 0), (1, 0, 1), (2, 2, 0)], vec![(1, 1, 1)]]);
    let weighted = Weighted::new(&reads, None);
    let gene_of = [0, 1, 0];
    let links = link::Links::new(&weighted, &gene_of);
    assert_eq!(links.len(), 3);
    assert_eq!(links.cell, vec![0, 0, 1]);
    assert_eq!(links.genes_of(&gene_of), vec![0, 1, 1]);
    assert_eq!(links.cell_ptr, vec![0, 2, 3]);
    // Cell 0's gene-0 link holds both of its gene-0 sites.
    assert_eq!(&links.site[links.ptr[0]..links.ptr[1]], &[0, 2]);
}

use super::blocks;

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// A link labelling under the block model: cell roles and gene roles.
type Labels = (Vec<u16>, Vec<u16>);

/// Two cells, two genes (sites 0-1 and site 2): four links, two cell roles
/// and two gene roles, the block shifts and cell baselines held fixed.
/// Every labelling's probability is computed here from the model's
/// definition (Dirichlet-multinomial urns for cells and genes, Binomial
/// reads at `logit p = site baseline + cell baseline + shift`), not from
/// the sampler's own log joint.
#[test]
fn block_labels_match_the_exact_posterior_given_shifts() {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 0, 2), (2, 1, 1)],
        vec![(0, 0, 3), (1, 2, 0), (2, 2, 0)],
    ]);
    let weighted = Weighted::new(&reads, None);
    let gene_of = [0, 0, 1];
    let links = link::Links::new(&weighted, &gene_of);
    let shifts = vec![1.2, -0.4, -1.0, 0.7];
    let cell_base = vec![0.5, -0.3];
    let s = blocks::Settings {
        cell_roles: 2,
        gene_roles: 2,
        alpha: 0.8,
        gene_alpha: 1.3,
        shift_sd: 1.0,
        cell_sd: 1.0,
    };
    let model = blocks::Model::new(&links, &gene_of, 3, s.clone());
    let mu = model.baseline().to_vec();
    let link_gene = links.genes_of(&gene_of);
    let n = links.len();
    let urn = |counts: &[f64], conc: f64| {
        let r = counts.len() as f64;
        ln_gamma(conc) - ln_gamma(conc + counts.iter().sum::<f64>())
            + counts
                .iter()
                .map(|&x| ln_gamma(x + conc / r) - ln_gamma(conc / r))
                .sum::<f64>()
    };
    // Each labelling, (cell roles, gene roles), with its log probability.
    let mut exact: Vec<(Labels, f64)> = Vec::new();
    for code in 0..1usize << (2 * n) {
        let a: Vec<u16> = (0..n).map(|l| ((code >> (2 * l)) & 1) as u16).collect();
        let b: Vec<u16> = (0..n).map(|l| ((code >> (2 * l + 1)) & 1) as u16).collect();
        let mut v = 0.0;
        for c in 0..2 {
            let mut counts = [0.0; 2];
            (0..n)
                .filter(|&l| links.cell[l] as usize == c)
                .for_each(|l| counts[a[l] as usize] += 1.0);
            v += urn(&counts, s.alpha);
        }
        for g in 0..2 {
            let mut counts = [0.0; 2];
            (0..n)
                .filter(|&l| link_gene[l] as usize == g)
                .for_each(|l| counts[b[l] as usize] += 1.0);
            v += urn(&counts, s.gene_alpha);
        }
        for l in 0..n {
            let shift =
                shifts[a[l] as usize * 2 + b[l] as usize] + cell_base[links.cell[l] as usize];
            let (sites, m, u) = links.of(l);
            for i in 0..sites.len() {
                let p = sigmoid(mu[sites[i] as usize] + shift);
                v += m[i] * p.ln() + u[i] * (1.0 - p).ln();
            }
        }
        exact.push(((a, b), v));
    }
    let mx = exact.iter().map(|e| e.1).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = exact.iter().map(|e| (e.1 - mx).exp()).sum();
    let labels: Vec<Labels> = [115, 116]
        .iter()
        .flat_map(|&seed| model.label_draws(&shifts, &cell_base, 500, 60_000, seed))
        .collect();
    let draws = labels.len() as f64;
    let mut freq: FxHashMap<Labels, f64> = FxHashMap::default();
    for z in labels {
        *freq.entry(z).or_default() += 1.0 / draws;
    }
    let tv: f64 = exact
        .iter()
        .map(|(z, l)| (freq.get(z).copied().unwrap_or(0.0) - (l - mx).exp() / total).abs())
        .sum::<f64>()
        / 2.0;
    assert!(
        tv < 0.02,
        "total variation {tv:.4} from the exact posterior"
    );
}

/// Cells of three planted cell roles, genes of two planted gene roles,
/// each (cell role, gene role) block shifting its sites' methylation.
fn planted_blocks(
    cells: usize,
    genes: usize,
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
    let level: Vec<f64> = (0..cells).map(|_| rng.random_range(-1.0..1.0)).collect();
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

/// The shift update, labels held fixed, leaves its conditional invariant:
/// draws of `delta[0][0]` match the exact conditional's mean and spread
/// (numerical integral of its prior times its block's Binomial reads),
/// with many reads in the block and with few, where the prior counts.
#[test]
fn block_shifts_follow_their_conditional() {
    check_shift_conditional(90, 20, 1.0);
    check_shift_conditional(6, 2, 0.5);
}

fn check_shift_conditional(cells: usize, genes: usize, shift_sd: f64) {
    let (reads, gene_of, cell_role, gene_role, _) = planted_blocks(cells, genes, 117);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = blocks::Settings {
        cell_roles: 3,
        gene_roles: 2,
        alpha: 1.0,
        gene_alpha: 1.0,
        shift_sd,
        cell_sd: 1.0,
    };
    let model = blocks::Model::new(&links, &gene_of, gene_of.len(), s.clone());
    let mu = model.baseline().to_vec();
    let link_gene = links.genes_of(&gene_of);
    let a: Vec<u16> = links.cell.iter().map(|&c| cell_role[c as usize]).collect();
    let b: Vec<u16> = link_gene.iter().map(|&g| gene_role[g as usize]).collect();
    let entries: Vec<(f64, f64, f64)> = (0..links.len())
        .filter(|&l| a[l] == 0 && b[l] == 0)
        .flat_map(|l| {
            let (sites, m, u) = links.of(l);
            (0..sites.len())
                .map(|i| (mu[sites[i] as usize], m[i], u[i]))
                .collect::<Vec<_>>()
        })
        .collect();
    let grid: Vec<f64> = (0..8000)
        .map(|i| -4.0 + 8.0 * (i as f64 + 0.5) / 8000.0)
        .collect();
    let terms: Vec<f64> = grid
        .iter()
        .map(|&d| {
            -d * d / (2.0 * s.shift_sd * s.shift_sd)
                + entries
                    .iter()
                    .map(|&(mu, m, u)| {
                        let p = sigmoid(mu + d);
                        m * p.ln() + u * (1.0 - p).ln()
                    })
                    .sum::<f64>()
        })
        .collect();
    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
    let mass: f64 = w.iter().sum();
    let exact = grid.iter().zip(&w).map(|(d, x)| d * x).sum::<f64>() / mass;
    let exact_sd = (grid
        .iter()
        .zip(&w)
        .map(|(d, x)| (d - exact).powi(2) * x)
        .sum::<f64>()
        / mass)
        .sqrt();
    let draws: Vec<f64> = model
        .shift_draws(&a, &b, 4000, 119)
        .into_iter()
        .skip(200)
        .collect();
    let sampled = draws.iter().sum::<f64>() / draws.len() as f64;
    let sd = (draws.iter().map(|d| (d - sampled).powi(2)).sum::<f64>() / draws.len() as f64).sqrt();
    assert!(
        (sampled - exact).abs() < 0.1 * exact_sd && (sd / exact_sd - 1.0).abs() < 0.1,
        "{cells} cells: delta[0][0] | labels: sampled mean {sampled:.4} sd {sd:.4}, exact mean {exact:.4} sd {exact_sd:.4}"
    );
}

/// Planted cell roles and gene roles are both recovered (cells by their
/// main cell role, genes by their main gene role) though every cell has its
/// own level, which its baseline takes.
#[test]
fn blocks_recover_planted_cell_and_gene_roles() {
    let (reads, gene_of, cell_role, gene_role, level) = planted_blocks(240, 40, 121);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = blocks::Settings {
        cell_roles: 3,
        gene_roles: 2,
        alpha: 0.5,
        gene_alpha: 0.5,
        shift_sd: 1.5,
        cell_sd: 1.0,
    };
    let model = blocks::Model::new(&links, &gene_of, gene_of.len(), s);
    let chains = run_mcmc_parallel(&model, &config(20, 100, 123), 2);
    // Both planted gene roles shift strongly, and both chains name the same
    // genes as each one's most enriched.
    let summary = blocks::summarize(&model, &chains);
    assert_eq!(summary.strong_gene_roles, vec![2, 2]);
    assert!(
        summary.top_gene_overlap.iter().all(|&o| o > 0.9),
        "top-gene overlap between chains {:?}",
        summary.top_gene_overlap
    );
    for c in &chains {
        let means = c.last().unwrap().means(false);
        let cells = summary::main_labels(&means.cell, 3);
        let genes = summary::main_labels(&means.gene, 2);
        let (ac, ag) = (ari(&cells, &cell_role), ari(&genes, &gene_role));
        assert!(ac > 0.9 && ag > 0.9, "ARI cells {ac:.3}, genes {ag:.3}");
        // Each cell's level is recovered as its baseline, not as a role.
        let r = correlation(&means.cell_base, &level);
        assert!(r > 0.8, "baseline vs planted level: r {r:.3}");
    }
}

fn correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let cov: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let vx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    let vy: f64 = y.iter().map(|b| (b - my).powi(2)).sum();
    cov / (vx * vy).sqrt()
}

/// A cell baseline's update, labels and shifts held fixed, leaves its
/// conditional invariant: draws of cell 0's baseline match the exact
/// conditional's mean and spread (its Normal prior times its links'
/// Binomial reads at their blocks' shifts), with many reads and with few.
#[test]
fn cell_baselines_follow_their_conditional() {
    check_baseline_conditional(0.5, 1.0);
    check_baseline_conditional(0.05, 0.5);
}

fn check_baseline_conditional(cover: f64, cell_sd: f64) {
    let mut rng = SmallRng::seed_from_u64(125);
    let fracs: Vec<Vec<f64>> = (0..4)
        .map(|_| (0..40).map(|_| rng.random_range(0.1..0.6)).collect())
        .collect();
    let refs: Vec<&[f64]> = fracs.iter().map(Vec::as_slice).collect();
    let reads = draw_cells(&refs, cover, &mut rng);
    let gene_of: Vec<u32> = (0..40u32).map(|s| s / 4).collect();
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let shifts = vec![0.8, -0.5, 0.3, -1.1];
    let s = blocks::Settings {
        cell_roles: 2,
        gene_roles: 2,
        alpha: 1.0,
        gene_alpha: 1.0,
        shift_sd: 1.0,
        cell_sd,
    };
    let model = blocks::Model::new(&links, &gene_of, 40, s);
    let mu = model.baseline().to_vec();
    let link_gene = links.genes_of(&gene_of);
    let a: Vec<u16> = (0..links.len()).map(|l| (l % 2) as u16).collect();
    let b: Vec<u16> = link_gene.iter().map(|&g| (g % 2) as u16).collect();
    let entries: Vec<(f64, f64, f64)> = (0..links.len())
        .filter(|&l| links.cell[l] == 0)
        .flat_map(|l| {
            let shift = shifts[a[l] as usize * 2 + b[l] as usize];
            let (sites, m, u) = links.of(l);
            (0..sites.len())
                .map(|i| (mu[sites[i] as usize] + shift, m[i], u[i]))
                .collect::<Vec<_>>()
        })
        .collect();
    let grid: Vec<f64> = (0..8000)
        .map(|i| -4.0 + 8.0 * (i as f64 + 0.5) / 8000.0)
        .collect();
    let terms: Vec<f64> = grid
        .iter()
        .map(|&v| {
            -v * v / (2.0 * cell_sd * cell_sd)
                + entries
                    .iter()
                    .map(|&(x, m, u)| {
                        let p = sigmoid(x + v);
                        m * p.ln() + u * (1.0 - p).ln()
                    })
                    .sum::<f64>()
        })
        .collect();
    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
    let mass: f64 = w.iter().sum();
    let exact = grid.iter().zip(&w).map(|(v, x)| v * x).sum::<f64>() / mass;
    let exact_sd = (grid
        .iter()
        .zip(&w)
        .map(|(v, x)| (v - exact).powi(2) * x)
        .sum::<f64>()
        / mass)
        .sqrt();
    let draws: Vec<f64> = model
        .baseline_draws(&a, &b, &shifts, 4000, 127)
        .into_iter()
        .skip(200)
        .collect();
    let sampled = draws.iter().sum::<f64>() / draws.len() as f64;
    let sd = (draws.iter().map(|v| (v - sampled).powi(2)).sum::<f64>() / draws.len() as f64).sqrt();
    assert!(
        (sampled - exact).abs() < 0.1 * exact_sd && (sd / exact_sd - 1.0).abs() < 0.1,
        "cover {cover}: cell 0 baseline | rest: sampled mean {sampled:.4} sd {sd:.4}, exact mean {exact:.4} sd {exact_sd:.4}"
    );
}

/// Each block's enriched genes and cells are the planted ones: a gene's
/// most enriched block lies in its planted gene role, a cell's in its
/// planted cell role, and a planted block's sites carry its shift.
#[test]
fn blocks_report_what_is_enriched_in_each_block() {
    let (reads, gene_of, cell_role, gene_role, _) = planted_blocks(240, 40, 129);
    let weighted = Weighted::new(&reads, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = blocks::Settings {
        cell_roles: 3,
        gene_roles: 2,
        alpha: 0.5,
        gene_alpha: 0.5,
        shift_sd: 1.5,
        cell_sd: 1.0,
    };
    let model = blocks::Model::new(&links, &gene_of, gene_of.len(), s);
    let chains = run_mcmc_parallel(&model, &config(20, 100, 131), 1);
    let means = chains[0].last().unwrap().means(true);
    let blocks = 6;
    let top = |enrichment: &[f64]| -> Vec<u16> {
        enrichment
            .chunks(blocks)
            .map(|row| {
                (0..blocks)
                    .max_by(|&a, &b| row[a].total_cmp(&row[b]))
                    .unwrap() as u16
            })
            .collect()
    };
    // A gene's top block names its gene role, a cell's its cell role.
    let genes: Vec<u16> = top(&means.gene_enrichment())
        .iter()
        .map(|k| k % 2)
        .collect();
    let cells: Vec<u16> = top(&means.cell_enrichment())
        .iter()
        .map(|k| k / 2)
        .collect();
    let (ag, ac) = (ari(&genes, &gene_role), ari(&cells, &cell_role));
    assert!(
        ag > 0.9 && ac > 0.9,
        "top-block ARI: genes {ag:.3}, cells {ac:.3}"
    );
    // Sites' lift in a block follows the block's shift, where the block
    // holds enough of the site's reads to say.
    let lift = means.site_lift(model.baseline());
    let (mut x, mut y) = (Vec::new(), Vec::new());
    for (at, &l) in lift.iter().enumerate() {
        if means.site_on[at] + means.site_off[at] >= 20.0 {
            x.push(l);
            y.push(means.shift[at % blocks]);
        }
    }
    let r = correlation(&x, &y);
    assert!(r > 0.8, "site lift vs block shift: r {r:.3}");
}

/// The held-out split: fit and held reads partition the input entry by
/// entry, the held share is near the fraction asked, and a seed gives the
/// same split every time (another seed, another split).
#[test]
fn hold_out_partitions_the_entries_by_seed() {
    let mut rng = SmallRng::seed_from_u64(133);
    let fracs: Vec<Vec<f64>> = (0..200).map(|_| vec![0.4; 100]).collect();
    let refs: Vec<&[f64]> = fracs.iter().map(Vec::as_slice).collect();
    let reads = draw_cells(&refs, 0.5, &mut rng);
    let (fit, held) = read::hold_out(&reads, 0.2, 9);
    assert_eq!((fit.cells(), held.cells()), (reads.cells(), reads.cells()));
    for c in 0..reads.cells() {
        let mut both: Vec<(u32, u32, u32)> = [fit.of(c), held.of(c)]
            .iter()
            .flat_map(|(s, m, u)| (0..s.len()).map(|i| (s[i], m[i], u[i])).collect::<Vec<_>>())
            .collect();
        both.sort_unstable();
        let (s, m, u) = reads.of(c);
        let all: Vec<(u32, u32, u32)> = (0..s.len()).map(|i| (s[i], m[i], u[i])).collect();
        assert_eq!(
            both, all,
            "cell {c}: fit and held do not partition its entries"
        );
    }
    let share = held.site.len() as f64 / reads.site.len() as f64;
    assert!((share - 0.2).abs() < 0.02, "held share {share:.3}");
    let again = read::hold_out(&reads, 0.2, 9).1;
    assert_eq!(again.site, held.site);
    assert_ne!(read::hold_out(&reads, 0.2, 10).1.site, held.site);
}

/// Held-out reads are scored as finite log-likelihoods by every model,
/// even at a site so deeply read and so fully converted that its fraction
/// rounds against 1.
#[test]
fn held_out_scores_are_finite_at_extreme_sites() {
    // Site 0: thousands of reads, all converted; sites 1-3 ordinary.
    let cells: Vec<Vec<(u32, u32, u32)>> = (0..12)
        .map(|c| vec![(0, 5000, 0), (1, 2, 1), (2, (c % 3) as u32, 2), (3, 1, 1)])
        .collect();
    let reads = Reads::from_cells(cells);
    let (fit, held) = read::hold_out(&reads, 0.5, 3);
    let gene_of = vec![0u32, 0, 1, 1];
    let weighted = Weighted::new(&fit, None);
    let links = link::Links::new(&weighted, &gene_of);
    let s = blocks::Settings {
        cell_roles: 2,
        gene_roles: 2,
        alpha: 1.0,
        gene_alpha: 1.0,
        shift_sd: 3.0,
        cell_sd: 3.0,
    };
    let model = blocks::Model::new(&links, &gene_of, 4, s).with_heldout(&held, &gene_of);
    let chains = run_mcmc_parallel(&model, &config(30, 20, 135), 1);
    let summary = blocks::summarize(&model, &chains);
    assert!(
        summary.heldout_per_read.is_some_and(f64::is_finite),
        "blocks held-out {:?}",
        summary.heldout_per_read
    );
    let s = link::Settings {
        k: 2,
        alpha: 0.5,
        gamma: 1.0,
        strength: 2.0,
    };
    let model = link::Model::new(&links, 4, s).with_heldout(&held);
    let chains = run_mcmc_parallel(&model, &config(30, 20, 137), 1);
    let summary = link::summarize(&model, &links, &chains);
    assert!(
        summary.heldout_per_read.is_some_and(f64::is_finite),
        "link held-out {:?}",
        summary.heldout_per_read
    );
}

#[test]
fn site_stats_count_every_cells_reads_and_the_cells_with_any() {
    // Two cells: one with reads at sites 0 and 2, one at site 0 only, an
    // entry with no reads not counted as a cell.
    let reads = super::read::Reads {
        ptr: vec![0, 2, 4],
        site: vec![0, 2, 0, 1],
        on: vec![3, 1, 2, 0],
        off: vec![1, 0, 4, 0],
    };
    let stats = super::read::SiteStats::of(&reads, 3);
    assert_eq!(stats.on, vec![5, 0, 1]);
    assert_eq!(stats.off, vec![5, 0, 0]);
    assert_eq!(stats.cells, vec![2, 0, 1]);
}
