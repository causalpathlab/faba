//! The thresholds two at a time: the axes a table has, the pairs Auto
//! picks, and a panel's squares and box.

use super::*;
use crate::qc::pairs::{self, Axis};

/// A table with a control arm and cells has every axis; without cells,
/// none of them.
#[test]
fn a_table_has_the_axes_of_its_criteria() {
    let t = table("m6a", 500);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(500)), &f);
    assert_eq!(pairs::axes(&v.sites()), Axis::ALL.to_vec());
    let v = SiteView::new(&t, None, &f);
    assert!(!pairs::axes(&v.sites()).contains(&Axis::Cells));
}

/// Auto picks four distinct pairs, the first the best scored, each split
/// counted as the sites failing exactly one of its thresholds, and spread
/// over more axes than any one pair holds.
#[test]
fn auto_picks_distinct_pairs_spread_over_the_axes() {
    let t = table("m6a", 2000);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(2000)), &f);
    let s = v.sites();
    let picked = pairs::auto_pairs(&s, pairs::PANELS);
    assert_eq!(picked.len(), 4);
    let best = (pairs::score_pairs(&s).iter().map(|p| p.score())).fold(f64::MIN, f64::max);
    assert_eq!(picked[0].score(), best);
    for p in &picked {
        let split = (0..t.len())
            .filter(|&i| (v.fails[i] & p.x.bits() != 0) != (v.fails[i] & p.y.bits() != 0))
            .count();
        assert_eq!(p.split, split, "{p:?}");
        assert!(p.rho.abs() <= 1.0 + 1e-9);
    }
    let shown: rustc_hash::FxHashSet<Axis> = picked.iter().flat_map(|p| [p.x, p.y]).collect();
    assert!(shown.len() >= 5, "{shown:?}");
    let keys: rustc_hash::FxHashSet<(Axis, Axis)> = picked.iter().map(|p| (p.x, p.y)).collect();
    assert_eq!(keys.len(), 4);
}

/// A panel bins every site with both values, kept apart from removed, and
/// its box holds the sites passing both drawn thresholds: those of them
/// removed fail some threshold not drawn there.
#[test]
fn a_panel_bins_every_site_and_counts_its_box() {
    let t = table("m6a", 1500);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(1500)), &f);
    let s = v.sites();
    let g = pairs::grid(&s, (Axis::Coverage, Axis::Pvalue), &f, pairs::BINS);
    let binned: u32 = g.squares.kept.iter().chain(&g.squares.other).sum();
    assert_eq!(binned as usize, t.len());
    let kept = v.fails.iter().filter(|&&m| m == 0).count();
    assert_eq!(g.squares.kept.iter().sum::<u32>() as usize, kept);
    let drawn = Axis::Coverage.bits() | Axis::Pvalue.bits();
    let in_box = |want_kept: bool| {
        (0..t.len())
            .filter(|&i| v.fails[i] & drawn == 0 && (v.fails[i] == 0) == want_kept)
            .count()
    };
    assert_eq!((g.in_kept, g.in_removed), (in_box(true), in_box(false)));
    assert!(
        g.in_box(false).starts_with("in the box: "),
        "{}",
        g.in_box(false)
    );
}

/// A p-value's upper bound is the box's lower end on -log10 p; an edit
/// ratio holds both ends; an off threshold leaves its end open.
#[test]
fn the_box_ends_follow_the_thresholds() {
    let t = table("m6a", 200);
    let mut f = SiteFilterArgs::default_values();
    f.site_max_pv = 0.01;
    f.site_min_edit_ratio = 0.1;
    f.site_max_edit_ratio = 0.9;
    let v = SiteView::new(&t, None, &f);
    let s = v.sites();
    let (lo, hi) = Axis::Pvalue.passes(&f, &s);
    assert!((lo.unwrap() - 2.0).abs() < 1e-6 && hi.is_none());
    assert_eq!(
        Axis::EditRatio.passes(&f, &s),
        (Some(0.1f32 as f64), Some(0.9f32 as f64))
    );
    let mut off = f.clone();
    Criterion::MinCoverage.set(&mut off, Criterion::MinCoverage.permissive());
    assert_eq!(Axis::Coverage.passes(&off, &s), (None, None));
}

/// The figure draws a panel per pair with its box and why it was picked.
#[test]
fn the_figure_draws_each_pair() {
    let t = table("m6a", 800);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(800)), &f);
    let s = v.sites();
    let picked: Vec<_> = (pairs::auto_pairs(&s, 4).into_iter())
        .map(|p| (p.x, p.y, p.why()))
        .collect();
    let svg = pairs::figure(
        &s,
        &picked,
        &f,
        (720.0, 640.0),
        "m6a sites, two thresholds at a time",
        &Default::default(),
    );
    assert!(
        svg.starts_with("<svg") && svg.contains("stroke-dasharray"),
        "{}",
        &svg[..200]
    );
    assert_eq!(svg.matches("in the box: ").count(), 4);
    assert!(svg.contains("auto: "));
}

/// Every fill of a drawing is a gray (or white).
fn all_gray(svg: &str) -> bool {
    svg.split("fill=\"#")
        .skip(1)
        .chain(svg.split("stroke=\"#").skip(1))
        .all(|rest| {
            let hex = &rest[..6];
            hex[0..2] == hex[2..4] && hex[2..4] == hex[4..6]
        })
}

/// The panels alone: no title, colour key or notes, only the squares,
/// their boxes, axes and axis titles; in gray, darker where more sites.
#[test]
fn the_pairs_panels_alone_in_gray() {
    use crate::figure::style::{Palette, Sheet, StyleArgs};
    let t = table("m6a", 800);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(800)), &f);
    let s = v.sites();
    let picked: Vec<_> = (pairs::auto_pairs(&s, 4).into_iter())
        .map(|p| (p.x, p.y, p.why()))
        .collect();
    let look = StyleArgs {
        font_size: None,
        panel_only: true,
        palette: Palette::Gray,
    };
    let sheet = Sheet {
        page: (Some(360.0), Some(320.0)),
        look,
    };
    let style = sheet.style((360.0, 320.0));
    let svg = pairs::figure(&s, &picked, &f, (360.0, 320.0), "TITLE1", &style);
    let first = format!("{} by {}", picked[0].1.short(), picked[0].0.short());
    for gone in [
        "TITLE1",
        "in the box",
        "auto:",
        "kept by every",
        "darker:",
        &first,
    ] {
        assert!(!svg.contains(gone), "{gone}");
    }
    assert!(svg.contains(picked[0].0.name()), "axis titles stay");
    assert!(all_gray(&svg), "{svg}");
    assert!(!all_gray(&pairs::figure(
        &s,
        &picked,
        &f,
        (360.0, 320.0),
        "T",
        &Default::default()
    )));
}

/// The cells figure alone: no title or key; in gray.
#[test]
fn the_cells_figure_alone_in_gray() {
    use crate::figure::style::{Palette, StyleArgs};
    let report = data_beans::qc_lib::QcReport {
        train_keep: vec![true, true, false, true],
        near_empty: vec![false, false, false, true],
        feature_keep: vec![],
        n_genes: vec![900.0, 700.0, 9000.0, 1.0],
        total_counts: vec![3000.0, 2000.0, 60000.0, 2.0],
        mito_frac: None,
        ribo_frac: None,
        n_cells_dropped: 1,
        n_features_dropped: 0,
    };
    let c = crate::qc::cells_plot::Cells::of("B1", &report);
    let look = StyleArgs {
        font_size: Some(6.0),
        panel_only: true,
        palette: Palette::Gray,
    };
    let style = look.style((300.0, 150.0), true);
    let svg = c.figure((300.0, 150.0), &style);
    for gone in ["B1 cells", "MAD outlier", "near empty", "right: "] {
        assert!(!svg.contains(gone), "{gone}");
    }
    assert!(svg.contains("barcode rank") && svg.contains("genes detected"));
    let labels = crate::figure::labels(&svg);
    assert!(labels.iter().all(|l| l.0 == "6"), "{labels:?}");
    assert!(all_gray(&svg), "{svg}");
}

/// The screen as text, drawn at `w` by `h`.
fn drawn(p: &mut SitePicker, (w, h): (u16, u16)) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let buf = term.backend().buffer().clone();
    (0..h)
        .map(|y| (0..w).map(|x| buf[(x, y)].symbol()).collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

/// `v` puts the four panels in place of the histogram, each with its box's
/// counts and why it is shown; `v` again brings the histogram back.
#[test]
fn v_shows_the_pairs_in_place_of_the_histogram() {
    let t = table("m6a", 600);
    let mut p = picker(&t, Some(cells(600)), SiteFilterArgs::default_values());
    let before = drawn(&mut p, (160, 50));
    assert!(!before.contains("in the box"), "{before}");
    press(&mut p, KeyCode::Char('v'));
    let s = drawn(&mut p, (160, 50));
    assert_eq!(s.matches("in the box").count(), 4, "{s}");
    assert!(s.contains("auto: "), "{s}");
    assert!(s.contains("(by it alone)"), "the table's counts: {s}");
    press(&mut p, KeyCode::Char('v'));
    assert!(!drawn(&mut p, (160, 50)).contains("in the box"));
}

/// Tab reaches the panels; there `x` and `y` step the focused panel's axes
/// (never onto the other one), and `A` puts back Auto's pairs; the
/// thresholds still move the panels' counts.
#[test]
fn the_panels_axes_step_and_auto_puts_them_back() {
    let t = table("m6a", 900);
    let mut p = picker(&t, Some(cells(900)), SiteFilterArgs::default_values());
    press(&mut p, KeyCode::Char('v'));
    let auto = p.view().pairs.clone();
    press(&mut p, KeyCode::Tab);
    press(&mut p, KeyCode::Tab);
    assert_eq!(p.panel, Panel::Pairs);
    press(&mut p, KeyCode::Right);
    assert_eq!(p.pair_at, 1);
    let (x0, y0, _) = p.view().pairs[1].clone();
    press(&mut p, KeyCode::Char('x'));
    press(&mut p, KeyCode::Char('Y'));
    let (x1, y1, why) = p.view().pairs[1].clone();
    assert!(x1 != x0 && y1 != y0 && x1 != y1, "{x1:?} {y1:?}");
    assert_eq!(why, "chosen by hand");
    press(&mut p, KeyCode::Char('A'));
    assert_eq!(p.view().pairs, auto);
    // Back at the table, a threshold moves the counts in the boxes.
    press(&mut p, KeyCode::Tab);
    assert_eq!(p.panel, Panel::Thresholds);
    let (x, y, _) = p.view().pairs[0];
    let before = pairs::grid(&p.view().sites(), (x, y), &p.filter, pairs::BINS).in_kept;
    focus_on(&mut p, Criterion::MinCoverage);
    for _ in 0..30 {
        press(&mut p, KeyCode::Char('+'));
    }
    let after = pairs::grid(&p.view().sites(), (x, y), &p.filter, pairs::BINS).in_kept;
    assert!(after < before, "{before} -> {after}");
}

/// A p-value of 0 is drawn at the axis' end, not stretching it.
#[test]
fn a_zero_p_value_does_not_stretch_the_axis() {
    let t = table("m6a", 600);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(600)), &f);
    let s = v.sites();
    assert!(t.pv.contains(&0.0));
    let (_, hi) = s.ruler(Axis::Pvalue, &f).span();
    assert!(hi < 20.0, "{hi}");
    let g = pairs::grid(&s, (Axis::LogOdds, Axis::Pvalue), &f, pairs::BINS);
    assert_eq!(
        g.squares.kept.iter().chain(&g.squares.other).sum::<u32>(),
        600
    );
}

/// The knob, pairs and cells figures at their own size, no font asked for,
/// are drawn as they always were; the knob and the pairs on each
/// quantity's one scale (see `crate::qc::scale`), as the site panels draw
/// it: the p-value's knob on linear bins.
#[test]
fn the_default_qc_figures_are_unchanged() {
    let t = table("m6a", 800);
    let f = SiteFilterArgs::default_values();
    let p = picker(&t, Some(cells(800)), f.clone());
    let knob = p.figure(crate::figure::style::Sheet::OWN).unwrap();
    let pairs = p.pairs_figure(0, crate::figure::style::Sheet::OWN);
    let report = data_beans::qc_lib::QcReport {
        train_keep: vec![true, true, false, true],
        near_empty: vec![false, false, false, true],
        feature_keep: vec![],
        n_genes: vec![900.0, 700.0, 9000.0, 1.0],
        total_counts: vec![3000.0, 2000.0, 60000.0, 2.0],
        mito_frac: None,
        ribo_frac: None,
        n_cells_dropped: 1,
        n_features_dropped: 0,
    };
    let cells =
        crate::qc::cells_plot::Cells::of("B1", &report).figure((720.0, 320.0), &Default::default());
    let shown = [&knob, &pairs, &cells].map(|s| crate::figure::fingerprint(s));
    assert_eq!(
        shown,
        [
            1047802977818497795,
            12215480049523947226,
            12252359045459261505
        ]
    );
}

/// Each axis of the 2D view is drawn on its quantity's ruler, on the scale
/// the cache holds (`--scale`), reaching the thresholds on: its range, its
/// values and its ticks.
#[test]
fn each_axis_is_drawn_on_its_quantitys_ruler() {
    use crate::qc::scale::{spread, Ruler, Scales};
    let t = table("m6a", 600);
    let f = SiteFilterArgs::default_values();
    let mut v = SiteView::new(&t, Some(cells(600)), &f);
    let scales = Scales::parse(&["pv=log", "coverage=linear"]).unwrap();
    v.cache = pairs::Cache::new(scales, 30);
    let s = v.sites();
    for a in Axis::ALL {
        let shown: Vec<f32> = (0..t.len()).map(|i| a.value(&s, i) as f32).collect();
        let on = |c: Criterion| s.criteria.contains(&c) && !c.is_off(&f);
        let (r, _) = Ruler::cut((a, scales.of(a)), spread(&shown), &f, on, 30);
        assert_eq!(s.ruler(a, &f), r, "{a:?}");
        let i = shown.iter().position(|v| v.is_finite()).unwrap();
        let col = s.column(a);
        assert_eq!(col.values[i], r.along(shown[i] as f64) as f32);
    }
    assert!(s.ruler(Axis::Pvalue, &f).axis.log);
    assert!(!s.ruler(Axis::Coverage, &f).axis.log);
    let g = pairs::grid(&s, (Axis::Coverage, Axis::Pvalue), &f, pairs::BINS);
    let ticks = [Axis::Coverage, Axis::Pvalue].map(|a| s.ruler(a, &f).ticks());
    assert_eq!(g.ticks, ticks);
    let ranges = [Axis::Coverage, Axis::Pvalue].map(|a| s.ruler(a, &f).span());
    assert_eq!([g.xr, g.yr], ranges);
}

/// A site whose p-value is not a number has no place on the -log10 p
/// axis: it is left out of a panel's squares, not drawn at 0.
#[test]
fn a_nan_p_value_is_left_out_of_the_panels() {
    let mut t = table("m6a", 300);
    t.pv[3] = f32::NAN;
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(300)), &f);
    let s = v.sites();
    assert!(Axis::Pvalue.value(&s, 3).is_nan());
    assert!(Axis::Pvalue.shown(f64::NAN).is_nan());
    let g = pairs::grid(&s, (Axis::Coverage, Axis::Pvalue), &f, pairs::BINS);
    let binned: u32 = g.squares.kept.iter().chain(&g.squares.other).sum();
    assert_eq!(binned as usize, t.len() - 1);
}
