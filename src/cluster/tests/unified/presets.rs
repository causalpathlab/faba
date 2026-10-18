//! The presets against the models they replaced, what they refuse,
//! recover, place and write, and the links they are fitted over.

use super::*;

/// The link preset (one gene role, no shift, finite strength, HDP) is
/// the link model it replaced: the same enumerated posterior on the tiny
/// fixture as that model's, and its moves (with no shift and no levels,
/// every move of a sweep) reach it.
#[test]
fn the_link_preset_is_the_link_model() {
    let links = tiny_links();
    let s = Settings {
        weights: Weights::Hdp(2.0),
        ..settings(Weights::Symmetric, false, 2.0, 1)
    };
    let model = Model::new(&links, &TINY_GENES, 3, s.clone()).expect("an allowed setting");
    let none = [0.0, 0.0];
    let (new, new_weight) = exact(&links, &TINY_GENES, model.baseline(), &s, (&none, &none));
    let (old, old_weight) = LINK_REFERENCE;
    for (code, (((a, _), p), q)) in new.iter().zip(old).enumerate() {
        let z: Vec<u16> = (0..links.len()).map(|l| (code >> l & 1) as u16).collect();
        assert_eq!(a, &z);
        assert!(
            (p - q).abs() < 1e-12,
            "z {z:?}: {p} here, {q} in the link model"
        );
    }
    assert!((new_weight - old_weight).abs() < 1e-12);
    let all: Vec<_> = [61, 62]
        .into_iter()
        .flat_map(|seed| model.label_draws((&none, &none), 500, 40_000, seed))
        .collect();
    let draws: Vec<Labels> = all.iter().map(|(a, b, _)| (a.clone(), b.clone())).collect();
    let tv = label_tv(draws, &new);
    assert!(tv < 0.02, "total variation {tv:.4} from the link posterior");
    let sampled = all
        .iter()
        .map(|(a, _, beta)| beta[a[0] as usize])
        .sum::<f64>()
        / all.len() as f64;
    assert!(
        (sampled - old_weight).abs() < 0.01,
        "E[beta of link 0's role]: sampled {sampled:.4}, link model {old_weight:.4}"
    );
}

/// The blocks preset (gene roles, the shift, infinite strength with the
/// cells' own levels, symmetric weights) is the blocks model it replaced:
/// the same enumerated posterior on the tiny fixture as that model's, at
/// unequal levels, and the label moves reach it.
#[test]
fn the_blocks_preset_is_the_blocks_model() {
    let links = tiny_links();
    let (shifts, levels) = (vec![1.2, -0.4, -1.0, 0.7], vec![0.5, -0.3]);
    let s = Settings {
        alpha: 0.8,
        ..settings(Weights::Symmetric, true, f64::INFINITY, 2)
    };
    let model = Model::new(&links, &TINY_GENES, 3, s.clone()).expect("an allowed setting");
    let old = exact_block_labels(
        &links,
        &TINY_GENES,
        model.baseline(),
        (s.alpha, s.gene_alpha),
        (&shifts, &levels),
    );
    let (new, _) = exact(
        &links,
        &TINY_GENES,
        model.baseline(),
        &s,
        (&shifts, &levels),
    );
    let by_labels: rustc_hash::FxHashMap<&Labels, f64> = old.iter().map(|(z, p)| (z, *p)).collect();
    assert_eq!(by_labels.len(), new.len());
    for (z, p) in &new {
        let q = by_labels[z];
        assert!(
            (p - q).abs() < 1e-9,
            "{z:?}: {p} here, {q} in the blocks model"
        );
    }
    let draws: Vec<Labels> = [117, 118]
        .iter()
        .flat_map(|&seed| model.label_draws((&shifts, &levels), 500, 40_000, seed))
        .map(|(a, b, _)| (a, b))
        .collect();
    let tv = label_tv(draws, &old);
    assert!(
        tv < 0.02,
        "total variation {tv:.4} from the blocks posterior"
    );
}

/// Settings under which no role reaches the reads are refused: gene roles
/// without the shift (a site belongs to one gene), and infinite strength
/// without it (every block's reads Binomial at the same `sigmoid(mu_s)`).
/// With the shift, or at finite strength, they are not.
#[test]
fn settings_whose_roles_never_reach_the_reads_are_refused() {
    let links = tiny_links();
    let new = |s| Model::new(&links, &TINY_GENES, 3, s);
    for s in [
        settings(Weights::Symmetric, false, 2.0, 2),
        settings(Weights::Hdp(1.0), false, f64::INFINITY, 1),
    ] {
        let refused = new(s.clone())
            .err()
            .unwrap_or_else(|| panic!("{s:?} accepted"));
        assert!(refused.to_string().contains("shift"), "{refused}");
    }
    for s in [
        settings(Weights::Symmetric, true, 2.0, 2),
        settings(Weights::Hdp(1.0), true, f64::INFINITY, 1),
        settings(Weights::Hdp(1.0), false, 2.0, 1),
    ] {
        assert!(new(s.clone()).is_ok(), "{s:?} refused");
    }
}

/// The link preset recovers planted cell patterns from the posterior-mean
/// memberships, scores held-out reads, and its tallies account for every
/// link and read.
#[test]
fn the_link_preset_recovers_planted_patterns_and_scores_held_out_reads() {
    let (reads, gene_of, truth) = planted_links(240, 60, 71);
    let (fit, held) = read::hold_out(&reads, 0.1, 5);
    let weighted = Weighted::new(&fit, None);
    let links = Links::new(&weighted, &gene_of);
    let s = Settings {
        cell_roles: 3,
        gene_roles: 1,
        alpha: 0.5,
        gene_alpha: 1.0,
        weights: Weights::Hdp(1.0),
        shift_sd: None,
        strength: 2.0,
        cell_sd: 1.0,
    };
    let model = Model::new(&links, &gene_of, gene_of.len(), s)
        .expect("an allowed setting")
        .with_heldout(&held);
    let chains = run_mcmc_parallel(&model, &config(20, 80, 74), 2);
    let reads_in: f64 = weighted.on.iter().chain(&weighted.off).sum();
    let held_reads: f64 = held.on.iter().chain(&held.off).map(|&x| x as f64).sum();
    for c in &chains {
        let d = c.last().unwrap();
        let means = d.means(true);
        let score = ari(&summary::main_labels(&means.cell, 3), &truth);
        assert!(score > 0.9, "ARI {score:.3}");
        assert!(means.gene.iter().all(|&w| w == 1.0));
        assert!(means.shift.iter().chain(&means.shift_sd).all(|&x| x == 0.0));
        assert!(means.level.is_empty(), "no level at finite strength");
        let sum = |v: &[f64]| v.iter().sum::<f64>();
        let n = links.len() as f64;
        assert!((sum(&means.gene_block) - n).abs() < 1e-6);
        assert!((sum(&means.cell_block) - n).abs() < 1e-6);
        let tallied = sum(&means.site_on) + sum(&means.site_off);
        assert!((tallied - reads_in).abs() < 1e-6, "{tallied} of {reads_in}");
        assert!(d.means(false).site_on.is_empty());
        let (ll, count) = d.heldout;
        assert_eq!(count, held_reads);
        assert!(ll.is_finite() && ll < 0.0, "held-out {ll}");
        assert!(d.log_joint.is_finite());
        assert!((sum(&means.beta) - 1.0).abs() < 1e-9);
    }
}

/// The blocks preset recovers planted cell and gene roles (cells by their
/// main cell role, genes by their main gene role) though every cell has its
/// own level within ±1, which its level takes, not a role; both planted
/// gene roles shift strongly, and both chains name the same genes as each
/// one's most enriched.
#[test]
fn the_blocks_preset_recovers_planted_cell_and_gene_roles() {
    let (reads, gene_of, cell_role, gene_role, level) = planted_blocks(240, 40, 1.0, 121);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let model = Model::new(&links, &gene_of, gene_of.len(), blocks_preset()).unwrap();
    let chains = run_mcmc_parallel(&model, &config(20, 100, 123), 2);
    let summary = unified::summarize(&model, &chains);
    assert_eq!(summary.strong_gene_roles, vec![2, 2]);
    assert!(
        summary.top_gene_overlap.iter().all(|&o| o > 0.9),
        "top-gene overlap between chains {:?}",
        summary.top_gene_overlap
    );
    assert!(summary
        .global_weights
        .iter()
        .all(|w| (w - 1.0 / 3.0).abs() < 1e-12));
    for c in &chains {
        let means = c.last().unwrap().means(false);
        let cells = summary::main_labels(&means.cell, 3);
        let genes = summary::main_labels(&means.gene, 2);
        let (ac, ag) = (ari(&cells, &cell_role), ari(&genes, &gene_role));
        assert!(ac >= 0.95 && ag >= 0.95, "ARI cells {ac:.3}, genes {ag:.3}");
        let r = correlation(&means.level, &level);
        assert!(r > 0.8, "level vs planted level: r {r:.3}");
    }
}

/// Each block's enriched genes and cells are the planted ones: a gene's
/// most enriched block lies in its planted gene role, a cell's in its
/// planted cell role, and a planted block's sites carry its shift.
#[test]
fn the_blocks_preset_reports_what_is_enriched_in_each_block() {
    let (reads, gene_of, cell_role, gene_role, _) = planted_blocks(240, 40, 0.0, 129);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let model = Model::new(&links, &gene_of, gene_of.len(), blocks_preset()).unwrap();
    let chains = run_mcmc_parallel(&model, &config(20, 100, 131), 1);
    let means = chains[0].last().unwrap().means(true);
    let blocks = 6;
    assert_eq!(means.blocks(), blocks);
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

/// Cells left out of the fit, placed against the fitted chain, fall with
/// the fitted cells of their own pattern: under the link preset, and under
/// the blocks preset, their gene roles summed out and each its own level
/// drawn, which follows the level planted.
#[test]
fn control_cells_are_placed_by_their_pattern() {
    let split = |reads: &Reads| {
        let fit: Vec<usize> = (0..240).collect();
        let rest: Vec<usize> = (240..360).collect();
        (reads.subset(&fit), reads.subset(&rest))
    };
    let (reads, gene_of, truth) = planted_links(360, 60, 109);
    let (r, gene_of_b, cell_role, _, level) = planted_blocks(360, 40, 1.0, 139);
    for (reads, gene_of, truth, s) in [
        (reads, gene_of, truth, link_preset(3, 0.5)),
        (r, gene_of_b, cell_role, blocks_preset()),
    ] {
        let (fit_reads, rest_reads) = split(&reads);
        let links = Links::new(&Weighted::new(&fit_reads, None), &gene_of);
        let model = Model::new(&links, &gene_of, gene_of.len(), s.clone()).unwrap();
        let chains = run_mcmc_parallel(&model, &config(20, 80, 111), 1);
        let means = chains[0].last().unwrap().means(true);
        let rest_links = Links::new(&Weighted::new(&rest_reads, None), &gene_of);
        let (placed, levels) = model.place(&rest_links, &means, 113);
        assert_eq!(placed.len(), 120 * 3);
        if s.strength.is_finite() {
            assert!(levels.is_empty(), "no level at finite strength");
        } else {
            let r = correlation(&levels, &level[240..]);
            assert!(r > 0.8, "placed level vs planted level: r {r:.3}");
        }
        assert!(placed
            .chunks(3)
            .all(|w| (w.iter().sum::<f32>() - 1.0).abs() < 1e-4));
        let fitted = summary::main_labels(&means.cell, 3);
        let placed = summary::main_labels(&placed, 3);
        // Placed cells and fitted cells of one pattern share a role.
        let all: Vec<u16> = fitted.iter().chain(&placed).copied().collect();
        let score = ari(&all, &truth);
        assert!(
            score > 0.9,
            "{s:?}: ARI {score:.3} over fitted and placed cells"
        );
    }
}

/// Every cell alike: under the link preset the reads condense into one
/// cell role instead of being split by their own outcomes.
#[test]
fn the_link_preset_condenses_without_structure() {
    let (reads, gene_of, _) = planted_links(240, 60, 79);
    // Re-draw every cell from one pattern only.
    let mut rng = SmallRng::seed_from_u64(80);
    let same: Vec<f64> = (0..gene_of.len())
        .map(|_| rng.random_range(0.1..0.6))
        .collect();
    let fracs: Vec<&[f64]> = (0..reads.cells()).map(|_| same.as_slice()).collect();
    let reads = draw_cells(&fracs, 0.3, &mut rng);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let model = Model::new(&links, &gene_of, gene_of.len(), link_preset(4, 0.5)).unwrap();
    let chains = run_mcmc_parallel(&model, &config(20, 150, 81), 2);
    for c in &chains {
        let means = c.last().unwrap().means(true);
        // Each cell role's share of the reads, its blocks' summed.
        let mut share = [0.0f64; 4];
        for (at, (on, off)) in means.site_on.iter().zip(&means.site_off).enumerate() {
            share[at % 4] += on + off;
        }
        let total: f64 = share.iter().sum();
        let top = share.iter().cloned().fold(0.0, f64::max) / total;
        assert!(top > 0.9, "read shares {share:?}");
    }
}

/// A smaller Dirichlet concentration gives cells more concentrated
/// memberships (the per-cell urn condenses harder).
#[test]
fn smaller_alpha_concentrates_memberships() {
    let (reads, gene_of, _) = planted_links(150, 40, 83);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let top = |alpha: f64| {
        let model = Model::new(&links, &gene_of, gene_of.len(), link_preset(3, alpha)).unwrap();
        let chains = run_mcmc_parallel(&model, &config(10, 40, 85), 1);
        let m = &chains[0].last().unwrap().means(false).cell;
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

/// Held-out reads are scored as finite log-likelihoods under both presets,
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
    let links = Links::new(&Weighted::new(&fit, None), &gene_of);
    let blocks = Settings {
        cell_roles: 2,
        shift_sd: Some(3.0),
        ..blocks_preset()
    };
    for s in [blocks, link_preset(2, 0.5)] {
        let model = Model::new(&links, &gene_of, 4, s.clone())
            .unwrap()
            .with_heldout(&held);
        let chains = run_mcmc_parallel(&model, &config(30, 20, 135), 1);
        let summary = unified::summarize(&model, &chains);
        assert!(
            summary.heldout_per_read.is_some_and(f64::is_finite),
            "{s:?}: held-out {:?}",
            summary.heldout_per_read
        );
    }
}

/// Both presets, run as `faba cluster` runs them on read input, write the
/// one table set `faba modules` reads: the link preset one gene module, its
/// shifts all 0, the blocks preset its gene roles ordered by shift range
/// and each cell's level; the placed cells' rows are in the cell roles, not
/// fitted.
#[test]
fn both_presets_write_the_tables_faba_modules_reads() {
    let (data, fit_reads) = planted_site_data(90, 80, 12, 141);
    let fitted: Vec<usize> = (0..80).collect();
    let weighted = Weighted::new(&fit_reads, None);
    let kept: Vec<usize> = (0..90).collect();
    let controls: Vec<usize> = (80..90).collect();
    for (model, gene_roles) in [("link", 1), ("blocks", 2)] {
        let tmp = tempfile::tempdir().unwrap();
        let args = crate::cluster::tests::parse_args(&[
            "cluster",
            "x.zarr",
            "--model",
            model,
            "--cell-roles",
            "3",
            "--gene-roles",
            &gene_roles.to_string(),
            "--chains",
            "2",
        ]);
        let fit = crate::cluster::FitData {
            data: &data,
            fitted: &fitted,
            weighted: &weighted,
            held: None,
            all_cells: &data.cells,
            kept: &kept,
            controls: &controls,
            progress: &crate::cluster::progress::Progress::new(2),
        };
        let mut record =
            crate::run_record::RunRecord::start("cluster", &tmp.path().to_string_lossy());
        let written = crate::cluster::run::blockmodel::blockmodel_fit(
            &args,
            &config(10, 20, 143),
            tmp.path(),
            &mut record,
            fit,
        )
        .unwrap();
        assert_eq!(
            written.record, "m6a_blockmodel.run.json",
            "{model}: its own record"
        );
        assert_eq!(written.modality, "m6a", "{model}");
        // The record says which chains the tables pool and how well.
        let saved = record.finish(&Ok(())).expect("the record written");
        let saved: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(saved).unwrap()).unwrap();
        let pooling = &saved["results"]["pooling"];
        let pooled = pooling["pooled"].as_array().expect("pooled chains");
        assert!(!pooled.is_empty(), "{model}: {pooling}");
        let excluded = pooling["excluded"].as_array().expect("excluded chains");
        assert_eq!(pooled.len() + excluded.len(), 2, "{model}: {pooling}");
        if !pooling["fallback"].is_null() {
            assert_eq!(pooled.len(), 1, "{model}: the best chain alone");
        }
        // The record's shifts are the tables' own.
        let blocks =
            crate::common::read_record_batch(tmp.path().join("m6a_blocks.parquet")).unwrap();
        let sd = blocks
            .column_by_name("sd")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow::array::Float32Array>()
            .unwrap()
            .values()
            .to_vec();
        let recorded = saved["results"]["blockmodel"]["shift_sd"]
            .as_array()
            .unwrap();
        assert_eq!(recorded.len(), sd.len(), "{model}");
        for (r, t) in recorded.iter().zip(&sd) {
            let r = r.as_f64().unwrap();
            assert!((r - f64::from(*t)).abs() < 1e-5, "{model}: {r} vs {t}");
        }
        let alignment = pooling["alignment"].as_array().unwrap();
        assert_eq!(alignment.len(), 2, "{model}: each chain's alignment");
        assert_eq!(alignment[0]["cell_quality"].as_array().unwrap().len(), 3);
        assert!(pooling["worst_rhat"].is_number() || pooling["worst_rhat"] == "inf");
        assert!(pooling["worst_ess"].is_number(), "{model}: {pooling}");
        assert!(pooling["not_converged"].is_array(), "{model}: {pooling}");
        let m = crate::cluster::browse::Modules::load(tmp.path(), "m6a").unwrap();
        assert_eq!((m.cell_roles, m.gene_roles), (3, gene_roles), "{model}");
        assert!(!m.approximate.sizes && !m.approximate.pileup, "{model}");
        let fitted_size: f64 = m.cell_sizes.iter().sum();
        assert!((fitted_size - 80.0).abs() < 1e-3, "{model}: {fitted_size}");
        // The cells' own levels, every kept cell's, under the blocks preset
        // (infinite strength) alone.
        let roles =
            crate::common::read_record_batch(tmp.path().join("m6a_cell_roles.parquet")).unwrap();
        let level = roles.column_by_name("level");
        if gene_roles == 1 {
            assert!(level.is_none(), "link: no level column");
        } else {
            let level = level.expect("blocks: a level column");
            assert_eq!(level.null_count(), 0, "every kept cell's level");
        }
        if gene_roles == 1 {
            assert_eq!(m.order, vec![0]);
            assert_eq!(m.shift_range(0), 0.0);
        } else {
            let ranges: Vec<f64> = m.order.iter().map(|&g| m.shift_range(g)).collect();
            assert!(ranges.windows(2).all(|w| w[0] >= w[1]), "{ranges:?}");
        }
    }
}

/// Under fractional read weights, the live reads per site and block stay
/// those the labels imply after many moves, and a block's site no link
/// reaches holds exactly none (no residue left by adding and taking away).
#[test]
fn live_reads_stay_those_the_labels_imply() {
    let (reads, gene_of, _, _, _) = planted_blocks(40, 6, 0.0, 151);
    let mut rng = SmallRng::seed_from_u64(152);
    let mut weight =
        |n: usize| -> Vec<f64> { (0..n).map(|_| rng.random_range(0.1..1.9)).collect() };
    let (on, off) = (weight(reads.on.len()), weight(reads.off.len()));
    let links = Links::new(&Weighted::new(&reads, Some((&on, &off))), &gene_of);
    let sites = gene_of.len();
    let s = Settings {
        cell_roles: 3,
        ..settings(Weights::Hdp(1.0), true, 2.0, 2)
    };
    let model = Model::new(&links, &gene_of, sites, s.clone()).unwrap();
    let (a, b, live_on, live_off) = model.live_reads(60, 153);
    let (on, off) = block_reads(&links, &s, (&a, &b), sites);
    let mut empty = 0;
    for at in 0..on.len() {
        if on[at] == 0.0 && off[at] == 0.0 {
            assert_eq!((live_on[at], live_off[at]), (0.0, 0.0), "entry {at}");
            empty += 1;
        } else {
            assert!(
                (live_on[at] - on[at]).abs() < 1e-9 && (live_off[at] - off[at]).abs() < 1e-9,
                "entry {at}: live ({}, {}), from the labels ({}, {})",
                live_on[at],
                live_off[at],
                on[at],
                off[at]
            );
        }
    }
    assert!(empty > 0, "some block's site is empty");
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
