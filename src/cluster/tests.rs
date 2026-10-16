use super::crp::{self, ln_beta, ln_pred, Chain, Model, Sample, Settings};
use super::fit::{self, OUTSIDE};
use super::read::{self, Cell, Reads, Site, SiteData, Weighted};
use super::summary::{self, ari};
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use rustc_hash::FxHashMap;
use statrs::function::gamma::ln_gamma;

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
    super::qualify_batches(&mut inputs);
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
}
