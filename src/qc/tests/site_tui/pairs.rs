//! The thresholds two at a time: the axes a table has, the pairs Auto
//! picks, and a panel's squares and box.

use super::*;
use crate::qc::pairs::{self, Axis, Sites};

/// The view's sites as [`Sites`].
fn sites<'a>(v: &'a SiteView) -> Sites<'a> {
    Sites {
        table: v.table,
        n_cells: v.n_cells.as_deref(),
        fails: &v.fails,
        criteria: &v.criteria,
    }
}

/// A table with a control arm and cells has every axis; without cells,
/// none of them.
#[test]
fn a_table_has_the_axes_of_its_criteria() {
    let t = table("m6a", 500);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(500)), &f);
    assert_eq!(pairs::axes(&sites(&v)), Axis::ALL.to_vec());
    let v = SiteView::new(&t, None, &f);
    assert!(!pairs::axes(&sites(&v)).contains(&Axis::Cells));
}

/// Auto picks four distinct pairs, the first the best scored, each split
/// counted as the sites failing exactly one of its thresholds, and spread
/// over more axes than any one pair holds.
#[test]
fn auto_picks_distinct_pairs_spread_over_the_axes() {
    let t = table("m6a", 2000);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(2000)), &f);
    let s = sites(&v);
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
    let s = sites(&v);
    let g = pairs::grid(&s, (Axis::Coverage, Axis::Pvalue), &f, pairs::BINS);
    let binned: u32 = g.kept.iter().chain(&g.removed).sum();
    assert_eq!(binned as usize, t.len());
    let kept = v.fails.iter().filter(|&&m| m == 0).count();
    assert_eq!(g.kept.iter().sum::<u32>() as usize, kept);
    let drawn = Axis::Coverage.bits() | Axis::Pvalue.bits();
    let in_box = |want_kept: bool| {
        (0..t.len())
            .filter(|&i| v.fails[i] & drawn == 0 && (v.fails[i] == 0) == want_kept)
            .count()
    };
    assert_eq!((g.in_kept, g.in_removed), (in_box(true), in_box(false)));
    assert!(g.in_box().starts_with("in the box: "), "{}", g.in_box());
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
    let s = sites(&v);
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
    let s = sites(&v);
    let picked: Vec<_> = (pairs::auto_pairs(&s, 4).into_iter())
        .map(|p| (p.x, p.y, p.why()))
        .collect();
    let svg = pairs::figure(
        &s,
        &picked,
        &f,
        (720.0, 640.0),
        "m6a sites, two thresholds at a time",
    );
    assert!(
        svg.starts_with("<svg") && svg.contains("stroke-dasharray"),
        "{}",
        &svg[..200]
    );
    assert_eq!(svg.matches("in the box: ").count(), 4);
    assert!(svg.contains("auto: "));
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

/// `d` puts the four panels in place of the histogram, each with its box's
/// counts and why it is shown; `d` again brings the histogram back.
#[test]
fn d_shows_the_pairs_in_place_of_the_histogram() {
    let t = table("m6a", 600);
    let mut p = picker(&t, Some(cells(600)), SiteFilterArgs::default_values());
    let before = drawn(&mut p, (160, 50));
    assert!(!before.contains("in the box"), "{before}");
    press(&mut p, KeyCode::Char('d'));
    let s = drawn(&mut p, (160, 50));
    assert_eq!(s.matches("in the box").count(), 4, "{s}");
    assert!(s.contains("auto: "), "{s}");
    assert!(s.contains("(by it alone)"), "the table's counts: {s}");
    press(&mut p, KeyCode::Char('d'));
    assert!(!drawn(&mut p, (160, 50)).contains("in the box"));
}

/// Tab reaches the panels; there `x` and `y` step the focused panel's axes
/// (never onto the other one), and `A` puts back Auto's pairs; the
/// thresholds still move the panels' counts.
#[test]
fn the_panels_axes_step_and_auto_puts_them_back() {
    let t = table("m6a", 900);
    let mut p = picker(&t, Some(cells(900)), SiteFilterArgs::default_values());
    press(&mut p, KeyCode::Char('d'));
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
    let before = pairs::grid(&sites(p.view()), (x, y), &p.filter, pairs::BINS).in_kept;
    focus_on(&mut p, Criterion::MinCoverage);
    for _ in 0..30 {
        press(&mut p, KeyCode::Char('+'));
    }
    let after = pairs::grid(&sites(p.view()), (x, y), &p.filter, pairs::BINS).in_kept;
    assert!(after < before, "{before} -> {after}");
}

/// A p-value of 0 is drawn at the axis' end, not stretching it.
#[test]
fn a_zero_p_value_does_not_stretch_the_axis() {
    let t = table("m6a", 600);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(600)), &f);
    let s = sites(&v);
    assert!(t.pv.contains(&0.0));
    let (_, hi) = pairs::range(&s, Axis::Pvalue);
    assert!(hi < 20.0, "{hi}");
    let g = pairs::grid(&s, (Axis::LogOdds, Axis::Pvalue), &f, pairs::BINS);
    assert_eq!(g.kept.iter().chain(&g.removed).sum::<u32>(), 600);
}
