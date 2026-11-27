use super::*;

#[test]
fn col_of_endpoints_and_midpoint() {
    let e = BinEdges::new(100, 200, 10);
    // first site -> first column, last site -> last column
    assert_eq!(e.col_of(100), 0);
    assert_eq!(e.col_of(200), 9);
    assert_eq!(e.col_of(150), 5);
    // out-of-range low clamps to col 0
    assert_eq!(e.col_of(50), 0);
    // degenerate single-position extent collapses to column 0
    assert_eq!(BinEdges::new(100, 100, 10).col_of(100), 0);
}

#[test]
fn bin_conserves_mass_without_log() {
    let e = BinEdges::new(0, 99, 10);
    let pts = [(0, 1.0), (5, 2.0), (50, 3.0), (99, 4.0)];
    let bins = e.bin(&pts, false);
    let total: f64 = bins.iter().sum();
    assert!((total - 10.0).abs() < 1e-9, "got {total}");
    // positions outside the extent are dropped
    let bins2 = e.bin(&[(0, 1.0), (1000, 9.0)], false);
    assert!((bins2.iter().sum::<f64>() - 1.0).abs() < 1e-9);
}

#[test]
fn x_px_spans_plot_band() {
    let e = BinEdges::new(100, 200, 10);
    assert!((e.x_px(100, 10.0, 80.0) - 10.0).abs() < 1e-4);
    assert!((e.x_px(200, 10.0, 80.0) - 90.0).abs() < 1e-4);
    assert!((e.x_px(150, 10.0, 80.0) - 50.0).abs() < 1e-4);
}

#[test]
fn robust_max_ignores_zeros() {
    assert_eq!(robust_max([0.0, 0.0].into_iter()), 0.0);
    let m = robust_max([1.0, 2.0, 3.0, 4.0].into_iter());
    assert!(m >= 3.0);
}

/// On a sqrt x the left end is spread: the first tenth of the columns
/// holds the first hundredth of the bases, and a column's range maps back
/// to that column.
#[test]
fn a_sqrt_x_spreads_the_left_end() {
    let e = BinEdges::new(1_000, 11_000, 10).on(Scale::Sqrt);
    assert_eq!(e.col_of(1_000), 0);
    assert_eq!(e.col_of(1_099), 0);
    assert_eq!(e.col_of(1_101), 1);
    assert_eq!(e.col_of(11_000), 9);
    for col in 0..10 {
        let (start, stop) = e.col_range(col);
        assert_eq!(e.col_of(start + (stop - start) / 2), col);
    }
    assert!((e.x_px(3_500, 0.0, 100.0) - 50.0).abs() < 1e-3);
    // Linear as it was.
    assert_eq!(BinEdges::new(1_000, 11_000, 10).col_of(1_101), 0);
}

/// On a log x over a short window the left columns are narrower than a
/// base: those hold none, and every other column's first base is binned
/// in it.
#[test]
fn a_log_xs_columns_start_where_they_are_binned() {
    let e = BinEdges::new(0, 1_000, 80).on(Scale::Log);
    let mut held = 0;
    for col in 0..80 {
        let (start, stop) = e.col_range(col);
        assert!(stop >= start, "{col}: {start}..{stop}");
        if stop > start {
            held += 1;
            assert_eq!(e.col_of(start), col, "{col}: {start}..{stop}");
        }
    }
    assert!(held < 80 && held > 40, "{held}");
}
