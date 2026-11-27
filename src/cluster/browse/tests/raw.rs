//! The raw view's reads: pooled into tiles, cut under too few reads, the
//! cell table's batches matched to the matrices read again, and the fit's
//! command line read back.

use super::*;

/// Two cell roles (three cells and one), two gene roles (two features
/// and one): each cell's reads at its features.
fn raw() -> Raw {
    Raw {
        unit: RawUnit::Gene,
        cells: vec![
            vec![
                vec![(0, 1, 4), (2, 5, 5)],
                vec![(0, 2, 2), (1, 0, 3)],
                vec![(1, 3, 3)],
            ],
            vec![vec![(2, 1, 1)]],
        ],
        starts: vec![0, 2, 3],
    }
}

/// The raw view's options, each block's cells and features as they are.
fn as_they_are() -> RawArgs {
    RawArgs {
        sort: RawSort::None,
        ..RawArgs::default()
    }
}

#[test]
fn a_block_sorted_by_value_puts_the_cells_and_features_with_most_dots_first() {
    let shown = RawArgs {
        min_reads: 1,
        ..RawArgs::default()
    };
    let t = raw().tile(0, 0, (10, 10), &shown);
    let at = |r: usize, q: usize| t.pixels[r * t.cols + q];
    // Dots at half methylated: the second cell's first feature (2 of 2)
    // and the third's second (3 of 3). Cells: the second and the third (a
    // dot each, a tie kept as they were), then the first (1 of 4, none);
    // features: a dot each, as they were.
    let row = |r: usize| {
        (0..t.cols)
            .map(|q| (at(r, q).on, at(r, q).total))
            .collect::<Vec<_>>()
    };
    assert_eq!(row(0), vec![(2.0, 2.0), (0.0, 3.0)]);
    assert_eq!(row(1), vec![(0.0, 0.0), (3.0, 3.0)]);
    assert_eq!(row(2), vec![(1.0, 4.0), (0.0, 0.0)]);
    let dots: Vec<bool> = t.pixels.iter().map(|p| p.dot()).collect();
    assert_eq!(dots, vec![true, false, false, true, false, false]);
}

#[test]
fn a_block_sorted_by_reads_puts_its_most_read_cells_first() {
    let shown = RawArgs {
        sort: RawSort::Reads,
        ..RawArgs::default()
    };
    let t = raw().tile(0, 0, (10, 1), &shown);
    let totals: Vec<f64> = t.pixels.iter().map(|p| p.total).collect();
    assert_eq!(totals, vec![5.0, 4.0, 3.0]);
}

#[test]
fn a_tile_has_a_pixel_for_each_cell_and_feature_when_there_is_room() {
    let t = raw().tile(0, 0, (10, 10), &as_they_are());
    assert_eq!((t.rows, t.cols), (3, 2));
    let at = |r: usize, q: usize| t.pixels[r * t.cols + q];
    assert_eq!(
        (at(0, 0).on, at(0, 0).total, at(0, 0).pairs),
        (1.0, 4.0, 1.0)
    );
    assert_eq!((at(1, 1).on, at(1, 1).total), (0.0, 3.0));
    // No reads: none to say.
    assert_eq!(at(0, 1).total, 0.0);
    // The other gene role's reads stay out.
    assert_eq!(t.pixels.iter().map(|p| p.total).sum::<f64>(), 12.0);
}

#[test]
fn a_small_tile_pools_the_cells_and_features_it_covers() {
    let t = raw().tile(0, 0, (1, 1), &as_they_are());
    assert_eq!((t.rows, t.cols), (1, 1));
    let p = t.pixels[0];
    assert_eq!((p.on, p.total, p.pairs), (6.0, 12.0, 6.0));
    assert_eq!(p.value(RawValue::Fraction, 3), Some(0.5));
    assert_eq!(p.value(RawValue::Methylated, 3), Some(1.0));
    assert_eq!(p.value(RawValue::Log1p, 3), Some(2.0f64.ln_1p()));
}

#[test]
fn a_pixel_under_the_least_reads_is_blank() {
    let p = Pooled {
        on: 1.0,
        total: 2.0,
        pairs: 1.0,
        ..Pooled::default()
    };
    assert_eq!(p.value(RawValue::Fraction, 3), None);
    assert_eq!(p.value(RawValue::Fraction, 2), Some(0.5));
    // At least a read, whatever is asked.
    assert_eq!(Pooled::default().value(RawValue::Fraction, 0), None);
}

#[test]
fn the_least_reads_is_one_at_least() {
    let parsed = cluster_raw(&["--raw-min-reads", "0"]);
    assert!(parsed.is_err());
    assert_eq!(cluster_raw(&[]).unwrap().min_reads, MIN_READS);
}

/// `--raw` options parsed from `words`.
fn cluster_raw(words: &[&str]) -> Result<RawArgs, clap::Error> {
    #[derive(clap::Parser)]
    struct Raw {
        #[command(flatten)]
        args: RawArgs,
    }
    let words = std::iter::once("raw").chain(words.iter().copied());
    Ok(<Raw as clap::Parser>::try_parse_from(words)?.args)
}

#[test]
fn a_qualified_batch_is_matched_to_the_one_read_by_its_name() {
    let read = ["rep1_wt", "rep2_wt"];
    assert_eq!(matched_batch("rep1_wt", &read), Some("rep1_wt"));
    assert_eq!(matched_batch("out/a/rep2_wt", &read), Some("rep2_wt"));
    assert_eq!(matched_batch("rep3_wt", &read), None);
    // Qualified by a folder the fit named from where it ran.
    let two = ["/work/out1/rep1", "/work/out2/rep1"];
    assert_eq!(matched_batch("out2/rep1", &two), Some("/work/out2/rep1"));
    assert_eq!(matched_batch("t2/rep1", &two), None);
    // Two read under one name: neither is guessed.
    let twice = ["a/rep1", "b/rep1"];
    assert_eq!(matched_batch("c/rep1", &twice), None);
}

#[test]
fn a_fits_inputs_are_read_from_its_record_with_their_batches() {
    let dir = tempfile::tempdir().unwrap();
    let matrix = dir.path().join("rep1_wt_m6a_site.zarr.zip");
    std::fs::write(&matrix, b"").unwrap();
    let path = matrix.to_string_lossy();
    let record = serde_json::json!({
        "job": "cluster",
        "command_line": ["faba", "cluster", "elsewhere.zarr.zip"],
        "inputs": {"matrices": [{"path": "/gone/elsewhere.zarr.zip"}]},
        "results": {"inputs": [{"path": path, "batch": "out1/rep1_wt"}]},
    });
    let name = crate::cluster::record_name("m6a");
    std::fs::write(dir.path().join(name), record.to_string()).unwrap();
    let inputs = fitted_inputs(dir.path(), "m6a", &[]).unwrap();
    assert_eq!(inputs.len(), 1);
    assert_eq!(
        (&*inputs[0].path, &*inputs[0].batch),
        (&*path, "out1/rep1_wt")
    );
}

#[test]
fn a_record_without_its_inputs_listed_says_to_name_them() {
    let dir = tempfile::tempdir().unwrap();
    let record = serde_json::json!({"job": "cluster", "inputs": {}});
    let name = crate::cluster::record_name("m6a");
    std::fs::write(dir.path().join(name), record.to_string()).unwrap();
    let Err(e) = fitted_inputs(dir.path(), "m6a", &[]) else {
        panic!("refused");
    };
    let e = e.to_string();
    assert!(e.contains("--raw-data"), "{e}");
}

#[test]
fn the_threshold_is_a_beta_and_the_colour_a_colour() {
    assert!(cluster_raw(&["--raw-threshold", "50"]).is_err());
    assert_eq!(
        cluster_raw(&["--raw-threshold", "0.3"]).unwrap().threshold,
        0.3
    );
    assert!(cluster_raw(&["--raw-colour", "x1"]).is_err());
    let red = cluster_raw(&["--raw-color", "#ff0000"]).unwrap();
    assert_eq!(red.colour, Some((255, 0, 0)));
}
