use super::read::{self, Cell, Reads, Site, SiteData, Weighted};
use super::{math, weights};
use rand::rngs::SmallRng;
use rand::RngExt;

mod command;
mod planted;
mod progress;

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
    let w = weights::batch_weights(&reads, &batch_of, 2, 1);
    let weighted = Weighted::new(&reads, Some((&w.on, &w.off)));
    let p = math::pooled_fraction(100.0, 200.0);
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

/// `line` parsed as `faba cluster` parses it.
pub(in crate::cluster) fn parse_args(line: &[&str]) -> super::ClusterArgs {
    use clap::{Args, FromArgMatches};
    let cmd = super::ClusterArgs::augment_args(clap::Command::new("cluster"));
    let m = cmd.try_get_matches_from(line).unwrap();
    super::ClusterArgs::from_arg_matches(&m).unwrap()
}

/// Settings that would make the fit meaningless are refused before
/// anything is read: roles out of range, concentrations and sds not
/// positive, no restart, no epoch.
#[test]
fn settings_that_cannot_fit_are_refused() {
    assert!(super::check(&parse_args(&["cluster", "x.zarr"])).is_ok());
    for (flag, value) in [
        ("--communities", "0"),
        ("--cell-roles", "65535"),
        ("--gene-roles", "0"),
        ("--membership-alpha", "0"),
        ("--gene-alpha", "-1"),
        ("--shift-sd", "0"),
        ("--cell-sd", "0"),
        ("--chains", "0"),
        ("--epochs", "0"),
        ("--tol", "-1"),
    ] {
        let arg = format!("{flag}={value}");
        let bad = parse_args(&["cluster", "x.zarr", &arg]);
        assert!(super::check(&bad).is_err(), "{flag} {value} accepted");
    }
}

/// The MCMC fits' flags (faba before 0.19) are refused by name, each, with
/// what to use instead; `--model blocks`, the model left, is taken.
#[test]
fn the_mcmc_flags_are_refused_by_name() {
    for line in [
        &["--model", "link"][..],
        &["--model", "crp"],
        &["--sweeps", "100"],
        &["--burn-in", "10"],
        &["--thin", "2"],
        &["--holdout", "0.1"],
        &["--profile-strength", "inf"],
        &["--prior-strength", "5"],
        &["--shift", "on"],
        &["--global", "hdp"],
        &["--community-gamma", "2"],
        &["--split-merge", "3"],
        &["--restricted-scans", "1"],
        &["--start-clusters", "4"],
        &["--min-cluster-cells", "4"],
    ] {
        let words: Vec<&str> = ["cluster", "x.zarr"].iter().chain(line).copied().collect();
        let refused = super::check(&parse_args(&words)).unwrap_err().to_string();
        let named = if line[0] == "--prior-strength" {
            "--profile-strength"
        } else {
            line[0]
        };
        assert!(
            refused.contains(named) && refused.contains("removed in faba 0.19"),
            "{line:?}: {refused}"
        );
    }
    let both = parse_args(&["cluster", "x.zarr", "--sweeps", "9", "--model", "crp"]);
    let refused = super::check(&both).unwrap_err().to_string();
    assert!(refused.starts_with("--model, --sweeps "), "{refused}");
    let blocks = parse_args(&["cluster", "x.zarr", "--model", "blocks"]);
    assert!(super::check(&blocks).is_ok());
    let count = parse_args(&["cluster", "x/", "--modality", "count", "--sweeps", "9"]);
    assert!(super::check(&count).is_err());
}

/// The site blockmodel `faba cluster` fits: cell roles and gene roles (8
/// unless given), the shift and each cell's own level, symmetric weights;
/// recorded as the blocks preset unless `--gene-roles` changes it.
#[test]
fn the_site_blockmodel_is_the_blocks_preset() {
    use super::settings::{Settings, Weights};
    let parse = |line: &[&str]| {
        let words: Vec<&str> = ["cluster", "x.zarr"].iter().chain(line).copied().collect();
        parse_args(&words)
    };
    let blocks = Settings {
        cell_roles: 8,
        gene_roles: 8,
        alpha: 0.5,
        gene_alpha: 0.5,
        weights: Weights::Symmetric,
        shift_sd: Some(1.0),
        strength: f64::INFINITY,
        cell_sd: 1.0,
    };
    assert_eq!(super::args::blockmodel(&parse(&[])), blocks);
    let set = parse(&[
        "--communities",
        "4",
        "--shift-sd",
        "0.5",
        "--cell-sd",
        "0.3",
    ]);
    assert_eq!(
        super::args::blockmodel(&set),
        Settings {
            cell_roles: 4,
            shift_sd: Some(0.5),
            cell_sd: 0.3,
            ..blocks.clone()
        }
    );
    let recorded = |line: &[&str]| super::args::recorded(&parse(line)).preset;
    assert_eq!(recorded(&[]).as_deref(), Some("blocks"));
    assert_eq!(recorded(&["--gene-roles", "8"]).as_deref(), Some("blocks"));
    assert_eq!(recorded(&["--gene-roles", "4"]), None);
}

/// The record of a fit is read whether faba wrote it before 0.19 (a
/// Gibbs-sampled link fit: drawn weights, finite strength, no shift) or
/// after (the blocks preset by variational inference).
#[test]
fn old_and_new_records_of_the_model_read() {
    use super::settings::{Recorded, Weights};
    let old: Recorded = serde_json::from_str(
        r#"{"preset":"link","settings":{"cell_roles":8,"gene_roles":1,"alpha":0.5,
            "gene_alpha":0.5,"weights":{"hdp":1.0},"shift_sd":null,"strength":5.0,
            "cell_sd":1.0}}"#,
    )
    .unwrap();
    assert_eq!(old.settings.weights, Weights::Hdp(1.0));
    assert_eq!((old.settings.shift_sd, old.settings.strength), (None, 5.0));
    assert_eq!(old.label(), "link preset");
    let new: Recorded = serde_json::from_str(
        r#"{"preset":"blocks","settings":{"cell_roles":8,"gene_roles":8,"alpha":0.5,
            "gene_alpha":0.5,"weights":"symmetric","shift_sd":1.0,"strength":"inf",
            "cell_sd":1.0}}"#,
    )
    .unwrap();
    assert_eq!(
        new,
        super::args::recorded(&parse_args(&["cluster", "x.zarr"]))
    );
    assert_eq!(new.label(), "blocks preset");
}
