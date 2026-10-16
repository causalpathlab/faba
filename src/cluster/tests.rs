use super::crp::{self, ln_beta, ln_pred, Chain, Model, Sample, Settings};
use super::read::{self, Cell, Reads, Site, SiteData};
use super::summary;
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig};
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use rustc_hash::FxHashMap;

fn settings(gibbs: bool, split_merge: usize, fixed_alpha: Option<f64>) -> Settings {
    Settings {
        strength: 2.0,
        split_merge,
        restricted: 2,
        start_clusters: 3,
        neighbours: 2,
        gibbs,
        fixed_alpha,
    }
}

/// Adjusted Rand index.
fn ari(x: &[u16], y: &[u16]) -> f64 {
    let mut n: FxHashMap<(u16, u16), f64> = FxHashMap::default();
    let (mut a, mut b): (FxHashMap<u16, f64>, FxHashMap<u16, f64>) = Default::default();
    for (&p, &q) in x.iter().zip(y) {
        *n.entry((p, q)).or_default() += 1.0;
        *a.entry(p).or_default() += 1.0;
        *b.entry(q).or_default() += 1.0;
    }
    let c2 = |v: f64| v * (v - 1.0) / 2.0;
    let index: f64 = n.values().map(|&v| c2(v)).sum();
    let (sa, sb): (f64, f64) = (
        a.values().map(|&v| c2(v)).sum(),
        b.values().map(|&v| c2(v)).sum(),
    );
    let expected = sa * sb / c2(x.len() as f64);
    (index - expected) / ((sa + sb) / 2.0 - expected)
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
    // Reweighted reads are not whole: the Beta-function path.
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

/// Exact posterior of every partition of four cells over two sites, against
/// the sampler's draw frequencies: with split-merge alone and with both moves.
#[test]
fn sampler_matches_the_exact_posterior_on_four_cells() {
    // Cells 0 and 1 mostly converted at site 0, cells 2 and 3 mostly not.
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 1, 1)],
        vec![(0, 2, 1), (1, 0, 2)],
        vec![(0, 0, 3), (1, 2, 0)],
        vec![(0, 1, 2), (1, 1, 1)],
    ]);
    let alpha = 0.8;
    for (gibbs, sm) in [(false, 3), (true, 1)] {
        let model = Model::new(&reads, 2, settings(gibbs, sm, Some(alpha)), None);
        // Exact: alpha^K prod Gamma(n_k) prod_k prod_s B(a + M, b + U) / B(a, b).
        let (a, b) = model_prior(&reads, 2.0);
        let parts = partitions(4);
        let logp: Vec<f64> = parts
            .iter()
            .map(|p| {
                let k = *p.iter().max().unwrap() as usize + 1;
                let mut v = k as f64 * f64::ln(alpha);
                for l in 0..k {
                    let members: Vec<usize> = (0..4).filter(|&c| p[c] as usize == l).collect();
                    v += statrs::function::gamma::ln_gamma(members.len() as f64);
                    for s in 0..2 {
                        let (mut m, mut u) = (0u32, 0u32);
                        for &c in &members {
                            let (ss, on, off) = reads.of(c);
                            for i in 0..ss.len() {
                                if ss[i] as usize == s {
                                    m += on[i];
                                    u += off[i];
                                }
                            }
                        }
                        v += ln_beta(a[s] + m as f64, b[s] + u as f64) - ln_beta(a[s], b[s]);
                    }
                }
                v
            })
            .collect();
        let mx = logp.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let z: f64 = logp.iter().map(|l| (l - mx).exp()).sum();
        let exact: Vec<f64> = logp.iter().map(|l| (l - mx).exp() / z).collect();

        let config = McmcConfig {
            n_samples: 20_000,
            warmup: 200,
            thin: 1,
            seed: 7,
        };
        let chains = run_mcmc_parallel(&model, &config, 2);
        let mut freq: FxHashMap<Vec<u16>, f64> = FxHashMap::default();
        let total = chains.iter().map(|c| c.samples.len()).sum::<usize>() as f64;
        for s in chains.iter().flat_map(|c| &c.samples) {
            *freq.entry(s.z.clone()).or_default() += 1.0 / total;
        }
        let tv: f64 = parts
            .iter()
            .zip(&exact)
            .map(|(p, e)| (freq.get(p).copied().unwrap_or(0.0) - e).abs())
            .sum::<f64>()
            / 2.0;
        assert!(
            tv < 0.02,
            "gibbs {gibbs}: total variation {tv:.4} from the exact posterior"
        );
    }
}

/// The prior `Beta(a_s, b_s)` the model uses, recomputed.
fn model_prior(reads: &Reads, strength: f64) -> (Vec<f64>, Vec<f64>) {
    let (mut on, mut tot) = ([0.0f64; 2], [0.0f64; 2]);
    for i in 0..reads.site.len() {
        on[reads.site[i] as usize] += reads.on[i] as f64;
        tot[reads.site[i] as usize] += (reads.on[i] + reads.off[i]) as f64;
    }
    let p: Vec<f64> = (0..2).map(|s| (on[s] + 0.5) / (tot[s] + 1.0)).collect();
    (
        p.iter().map(|p| strength * p).collect(),
        p.iter().map(|p| strength * (1.0 - p)).collect(),
    )
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
    let per_cell = truth
        .iter()
        .map(|&t| {
            let mut entries = Vec::new();
            for (s, &p) in frac[t as usize].iter().enumerate() {
                if rng.random::<f64>() < 0.15 {
                    let n = rng.random_range(1..4u32);
                    let m = (0..n).filter(|_| rng.random::<f64>() < p).count() as u32;
                    entries.push((s as u32, m, n - m));
                }
            }
            entries
        })
        .collect();
    let reads = Reads::from_cells(per_cell);
    let model = Model::new(&reads, sites, settings(true, 5, None), None);
    let config = McmcConfig {
        n_samples: 60,
        warmup: 60,
        thin: 1,
        seed: 11,
    };
    let chains: Vec<Chain> = run_mcmc_parallel(&model, &config, 2);
    // A side group of ten-odd cells is all but tied with folding it back
    // (the posterior barely tells them apart), so count clusters of 20 up.
    let (part, _) = super::fit::settle(&reads, &summary::partition(&chains), 20, sites, 2.0);
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
    let draw = |rng: &mut SmallRng| {
        let k = rng.random_range(1..5usize);
        let raw: Vec<u32> = (0..n).map(|_| rng.random_range(0..k) as u32).collect();
        let mut map: FxHashMap<u32, u16> = FxHashMap::default();
        let z: Vec<u16> = raw
            .iter()
            .map(|&x| {
                let next = map.len() as u16;
                *map.entry(x).or_insert(next)
            })
            .collect();
        Sample {
            k: map.len(),
            z,
            alpha: 1.0,
            log_post: 0.0,
        }
    };
    let chains: Vec<Chain> = (0..2)
        .map(|_| Chain {
            samples: (0..15).map(|_| draw(&mut rng)).collect(),
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
fn filter_drops_thin_cells_then_thin_sites() {
    let site = |p: i64| Site {
        gene: "GENE1".into(),
        chr: "chr1".into(),
        pos: p,
    };
    let cell = |b: &str, fit: bool| Cell {
        barcode: b.into(),
        batch: "B1".into(),
        fit,
    };
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
    // D is not fitted) goes.
    let (out, kept) = read::filter(&data, 2, 2);
    assert_eq!(kept.cells, vec![0, 1, 3]);
    assert_eq!(
        out.sites.iter().map(|s| s.pos).collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(out.reads.of(2).0, &[0]);
    let _ = crp::Settings::clone(&settings(true, 0, None));
}

#[test]
fn a_sites_two_channels_merge_into_one_entry() {
    let reads = Reads::from_cells(vec![vec![(1, 2, 0), (0, 1, 0), (1, 0, 3)], vec![]]);
    assert_eq!(reads.of(0), (&[0, 1][..], &[1, 2][..], &[0, 3][..]));
    assert_eq!(reads.of(1).0.len(), 0);
}

/// Cells drawn alike in two batches, one batch shifted at every site: the
/// shift alone makes two clusters; weighted out, one.
#[test]
fn batch_shifts_are_cancelled() {
    let (cells, sites) = (240usize, 300usize);
    let mut rng = SmallRng::seed_from_u64(9);
    let base: Vec<f64> = (0..sites).map(|_| rng.random_range(0.1..0.5)).collect();
    let batch_of: Vec<u16> = (0..cells).map(|c| (c % 2) as u16).collect();
    let per_cell = batch_of
        .iter()
        .map(|&b| {
            let mut entries = Vec::new();
            for (s, &p) in base.iter().enumerate() {
                if rng.random::<f64>() < 0.2 {
                    // Batch 1 reads more converted everywhere: logit + 1.5.
                    let p = if b == 1 {
                        let logit = (p / (1.0 - p)).ln() + 1.5;
                        1.0 / (1.0 + (-logit).exp())
                    } else {
                        p
                    };
                    let n = rng.random_range(1..4u32);
                    let m = (0..n).filter(|_| rng.random::<f64>() < p).count() as u32;
                    entries.push((s as u32, m, n - m));
                }
            }
            entries
        })
        .collect();
    let reads = Reads::from_cells(per_cell);
    let config = McmcConfig {
        n_samples: 20,
        warmup: 20,
        thin: 1,
        seed: 13,
    };
    let clusters = |weights: Option<(&[f64], &[f64])>| {
        let model = Model::new(&reads, sites, settings(true, 5, None), weights);
        let chains = run_mcmc_parallel(&model, &config, 2);
        super::fit::settle(&reads, &summary::partition(&chains), 10, sites, 2.0).0
    };
    let raw = clusters(None);
    assert!(
        raw.k >= 2 && ari(&raw.z, &batch_of) > 0.8,
        "uncorrected: K {} ARI {:.2}",
        raw.k,
        ari(&raw.z, &batch_of)
    );
    let w = super::fit::batch_weights(&reads, &batch_of, 2, sites);
    let fixed = clusters(Some((&w.on, &w.off)));
    assert_eq!(fixed.k, 1, "corrected: K {}", fixed.k);
}

/// A cell with one cluster's pattern but a lower level everywhere goes to
/// that cluster, not to the one whose sites are all lower.
#[test]
fn assignment_reads_pattern_not_level() {
    // Cluster 0: sites 0-1 high, 2-3 low. Cluster 1: everything lowish,
    // flat. Each cluster: 20 cells.
    let profile = [[0.8, 0.8, 0.2, 0.2], [0.15, 0.15, 0.15, 0.15]];
    let mut cells = Vec::new();
    for prof in &profile {
        for _ in 0..20 {
            cells.push(
                prof.iter()
                    .enumerate()
                    .map(|(s, &p)| {
                        let m = (p * 20.0f64).round() as u32;
                        (s as u32, m, 20 - m)
                    })
                    .collect::<Vec<_>>(),
            );
        }
    }
    let reads = Reads::from_cells(cells);
    let z: Vec<u16> = (0..40).map(|c| (c / 20) as u16).collect();
    let tiles = super::fit::Tiles::new(&reads, &z, 2, 4, 2.0);
    // Cluster 0's pattern (high, high, low, low), shifted down to 0.4/0.4/0.05/0.05.
    let probe = tiles.assign(&[0, 1, 2, 3], &[8, 8, 1, 1], &[12, 12, 19, 19], &[20, 20]);
    assert!(probe[0] > 0.9, "{probe:?}");
}

#[test]
fn ari_is_one_for_relabelled_partitions_and_near_zero_for_unrelated() {
    let draw = |z: Vec<u16>| Sample {
        k: *z.iter().max().unwrap() as usize + 1,
        z,
        alpha: 1.0,
        log_post: 0.0,
    };
    let a = draw(vec![0, 0, 1, 1, 2, 2]);
    let b = draw(vec![1, 1, 2, 2, 0, 0]);
    assert!((summary::ari(&a, &b) - 1.0).abs() < 1e-12);
    let c = draw(vec![0, 1, 2, 0, 1, 2]);
    assert!(summary::ari(&a, &c) < 0.1);
}

/// The exact posterior over every partition of `reads`' cells (two sites),
/// `alpha` fixed or, `None`, integrated out under its Gamma(1, 1) prior.
fn exact_posterior(reads: &Reads, strength: f64, alpha: Option<f64>) -> Vec<(Vec<u16>, f64)> {
    let n = reads.cells();
    let (a, b) = model_prior(reads, strength);
    let ln_gamma = statrs::function::gamma::ln_gamma;
    let parts = partitions(n);
    let logp: Vec<f64> = parts
        .iter()
        .map(|p| {
            let k = *p.iter().max().unwrap() as usize + 1;
            let mut v = 0.0;
            for l in 0..k {
                let members: Vec<usize> = (0..n).filter(|&c| p[c] as usize == l).collect();
                v += ln_gamma(members.len() as f64);
                for s in 0..2 {
                    let (mut m, mut u) = (0u32, 0u32);
                    for &c in &members {
                        let (ss, on, off) = reads.of(c);
                        for i in 0..ss.len() {
                            if ss[i] as usize == s {
                                m += on[i];
                                u += off[i];
                            }
                        }
                    }
                    v += ln_beta(a[s] + m as f64, b[s] + u as f64) - ln_beta(a[s], b[s]);
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

/// Every move together, alpha resampled, against the exact posterior over
/// all 203 partitions of six cells: the chain leaves the posterior invariant
/// (detailed balance of each move) and reaches it.
#[test]
fn all_moves_with_alpha_resampled_match_the_exact_posterior() {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 1, 1)],
        vec![(0, 2, 1), (1, 0, 2)],
        vec![(0, 0, 3), (1, 2, 0)],
        vec![(0, 1, 2), (1, 1, 1)],
        vec![(0, 2, 0), (1, 0, 1)],
        vec![(0, 0, 2), (1, 2, 1)],
    ]);
    let exact = exact_posterior(&reads, 2.0, None);
    let model = Model::new(&reads, 2, settings(true, 2, None), None);
    let config = McmcConfig {
        n_samples: 40_000,
        warmup: 500,
        thin: 1,
        seed: 21,
    };
    let chains = run_mcmc_parallel(&model, &config, 2);
    let total = chains.iter().map(|c| c.samples.len()).sum::<usize>() as f64;
    let mut freq: FxHashMap<Vec<u16>, f64> = FxHashMap::default();
    for s in chains.iter().flat_map(|c| &c.samples) {
        *freq.entry(s.z.clone()).or_default() += 1.0 / total;
    }
    let tv: f64 = exact
        .iter()
        .map(|(p, e)| (freq.get(p).copied().unwrap_or(0.0) - e).abs())
        .sum::<f64>()
        / 2.0;
    assert!(
        tv < 0.03,
        "total variation {tv:.4} from the exact posterior"
    );
}

/// Weighted, each batch's reads at a site keep their total and take the
/// site's pooled fraction (up to the shrinkage of thinly read batches).
#[test]
fn batch_weights_keep_totals_and_match_the_pooled_fraction() {
    // One site; batch 0: 80 of 100 converted; batch 1: 20 of 100.
    let mut cells = Vec::new();
    for (b, m) in [(0u16, 80u32), (1, 20)] {
        for _ in 0..10 {
            cells.push((b, vec![(0u32, m / 10, 10 - m / 10)]));
        }
    }
    let batch_of: Vec<u16> = cells.iter().map(|c| c.0).collect();
    let reads = Reads::from_cells(cells.into_iter().map(|c| c.1).collect());
    let w = super::fit::batch_weights(&reads, &batch_of, 2, 1);
    let p = (100.0 + 0.5) / (200.0 + 1.0);
    for b in 0..2u16 {
        let (mut won, mut wtot) = (0.0, 0.0);
        for (c, &bc) in batch_of.iter().enumerate().filter(|(_, &x)| x == b) {
            let (_, m, u) = reads.of(c);
            won += w.on[c] * m[0] as f64;
            wtot += w.on[c] * m[0] as f64 + w.off[c] * u[0] as f64;
            assert_eq!(bc, b);
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
