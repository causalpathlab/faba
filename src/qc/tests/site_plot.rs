//! The site panels: what is asked, one ruler per quantity across every
//! panel, infinities at the ends, kept and removed as the thresholds cut,
//! the size and font asked, and cells per site off the `_site` matrices.

use super::*;
use crate::figure::style::{Palette, Sheet, StyleArgs};
use crate::qc::fixtures::{cells, site_table, Row};

/// An m6A table of `rows`, its derived columns as `read_site_table` makes
/// them.
fn table(rows: &[Row]) -> SiteTable {
    site_table("m6a", 3, rows, |i| 100 * (i + 1))
}

/// `n` sites spread over every quantity: a clean control every fifth, a
/// p-value of 0 every 31st.
fn many(n: usize) -> SiteTable {
    let rows: Vec<Row> = (0..n)
        .map(|i| {
            let cov = 2 + (i * 37 % 150) as u64;
            let conv = (i * 13 % 40) as u64 % (cov + 1);
            let pv = match i % 31 {
                0 => 0.0,
                _ => 10f32.powf(-((i * 7919 % 600) as f32) / 100.0),
            };
            let (ctl_cov, ctl_conv) = if i % 5 == 0 {
                (12, 0)
            } else {
                (20, (i % 7) as u64)
            };
            (pv, cov, conv, ctl_cov, ctl_conv)
        })
        .collect();
    table(&rows)
}

fn args(hist: &[&str], pairs: &[&str], marginals: &[&str], scale: &[&str]) -> PanelArgs {
    let own = |v: &[&str]| v.iter().map(|s| s.to_string()).collect();
    PanelArgs {
        hist: own(hist),
        pairs: own(pairs),
        marginals: own(marginals),
        scale: own(scale),
        bins: None,
        xlab: Vec::new(),
        ylab: Vec::new(),
    }
}

fn sheet(w: f64, h: f64, font: Option<f64>, panel_only: bool) -> Sheet {
    Sheet {
        page: (Some(w), Some(h)),
        look: StyleArgs {
            font_size: font,
            panel_only,
            palette: Palette::Colour,
        },
    }
}

/// `--pairs` reads `x:y`, `--hist` and `--marginals` their names; an
/// unknown name is refused, naming the ones there are; `none` asks for
/// none, nothing given for Auto's (or every histogram).
#[test]
fn the_panels_asked_for_are_read_by_name() {
    let p = args(&[], &["pv:cells", "coverage:cells", "fold:pv"], &[], &[])
        .plan()
        .unwrap();
    use Quantity::*;
    assert_eq!(
        p.pairs,
        Choice::Given(vec![(Pvalue, Cells), (Coverage, Cells), (Fold, Pvalue)])
    );
    assert_eq!(p.hist, Choice::Auto);
    let p = args(
        &["log-odds", "edit_ratio"],
        &["none"],
        &["top", "right"],
        &[],
    )
    .plan()
    .unwrap();
    assert_eq!(p.hist, Choice::Given(vec![LogOdds, EditRatio]));
    assert_eq!(p.pairs, Choice::Given(vec![]));
    assert_eq!(
        p.marginals,
        Marginals {
            top: true,
            right: true
        }
    );
    assert_eq!(
        args(&["none"], &[], &[], &[]).plan().unwrap().hist,
        Choice::Given(vec![])
    );
    for bad in [
        args(&[], &["pv-cells"], &[], &[]),
        args(&[], &["pv:pv"], &[], &[]),
        args(&["pvx"], &[], &[], &[]),
        args(&[], &[], &["left"], &[]),
        args(&[], &[], &[], &["log_odds=log"]),
    ] {
        assert!(bad.plan().is_err(), "{:?}", bad.plan().err());
    }
    let e = args(&[], &["pv:bogus"], &[], &[])
        .plan()
        .unwrap_err()
        .to_string();
    assert!(e.contains("bogus") && e.contains("edit_ratio"), "{e}");
}

/// `--xlab` and `--ylab` title a quantity's axis, `q=TEXT`, in every
/// panel of it: on the axis named, else on the other; `--ylab sites=`
/// titles the histograms' counts. The rest keep their own titles; a title
/// that names no quantity is refused.
#[test]
fn axis_titles_are_given_by_quantity() {
    use Quantity::*;
    let mut a = args(&["pv", "cells"], &["pv:cells", "cells:coverage"], &[], &[]);
    a.xlab = vec!["pv=-log10 p (wt vs. mut)".into()];
    a.ylab = vec![
        "cells=cells with a converted read".into(),
        "sites=number of sites".into(),
    ];
    let plan = a.plan().unwrap();
    let t = many(300);
    let n = cells(t.len());
    let s = sheet(144.0, 130.0, Some(7.0), true);
    let p = Panels::new(&t, Some(&n), &Cut::none(), &plan, s, true).unwrap();
    let h = p.hist(Pvalue);
    assert_eq!(
        (&*h.x_title, &*h.y_title),
        ("-log10 p (wt vs. mut)", "number of sites")
    );
    assert_eq!(p.hist(Cells).x_title, "cells with a converted read");
    let d = p.pair(Pvalue, Cells);
    assert_eq!(
        (&*d.x_title, &*d.y_title),
        ("-log10 p (wt vs. mut)", "cells with a converted read")
    );
    let d = p.pair(Cells, Coverage);
    assert_eq!(
        (&*d.x_title, &*d.y_title),
        ("cells with a converted read", "coverage")
    );
    let svg = p.figures().into_iter().map(|f| f.svg).collect::<String>();
    assert!(svg.contains(">-log10 p (wt vs. mut)<"), "{svg}");
    assert!(svg.contains(">cells with a converted read<"), "{svg}");
    // Without any, each its own.
    let plain = args(&["pv"], &["none"], &[], &[]).plan().unwrap();
    let p = Panels::new(&t, Some(&n), &Cut::none(), &plain, s, true).unwrap();
    let h = p.hist(Pvalue);
    assert_eq!((&*h.x_title, &*h.y_title), ("-log10 p", "sites"));
    for (x, y) in [
        (&["-log10 p"][..], &[][..]),
        (&["bogus=TEXT1"], &[]),
        (&["sites=TEXT1"], &[]),
        (&[], &["pv"]),
    ] {
        let mut bad = args(&[], &[], &[], &[]);
        bad.xlab = x.iter().map(|s| s.to_string()).collect();
        bad.ylab = y.iter().map(|s| s.to_string()).collect();
        let e = bad.plan().unwrap_err().to_string();
        assert!(e.contains("--xlab") || e.contains("--ylab"), "{e}");
    }
}

/// The sites `faba qc`'s 2D view reads, with every criterion and the fail
/// bits of `f`, on the axes in `cache`.
fn qc_sites<'a>(
    t: &'a SiteTable,
    n: &'a [usize],
    criteria: &'a [Criterion],
    fails: &'a mut Vec<u8>,
    f: &SiteFilterArgs,
    cache: &'a crate::qc::pairs::Cache,
) -> crate::qc::pairs::Sites<'a> {
    *fails = (0..t.len())
        .map(|i| {
            (criteria.iter())
                .filter(|c| c.fails(f, t, i, Some(n[i])))
                .fold(0, |b, &c| b | 1 << c as u8)
        })
        .collect();
    crate::qc::pairs::Sites {
        table: t,
        n_cells: Some(n),
        fails,
        criteria,
        cache,
    }
}

/// A quantity's ruler is one: its histogram's axis, its axis in every 2D
/// panel and in their marginals (the panel's own), on the scale `--scale`
/// set; and what `faba qc`'s 2D view draws it on, a threshold beyond the
/// data's range and all: the axis reaches it in both.
#[test]
fn one_ruler_serves_a_histogram_a_2d_axis_and_its_marginal() {
    let t = many(800);
    let n = cells(800);
    let plan = args(
        &["pv", "fold", "coverage"],
        &["pv:cells", "fold:pv", "coverage:pv"],
        &["top", "right"],
        &["pv=log", "fold=linear", "coverage=linear"],
    )
    .plan()
    .unwrap();
    // Coverage under 160, its threshold far over it.
    let cut = Cut::of_qc(&SiteFilterArgs {
        site_min_coverage: 5000,
        ..SiteFilterArgs::default_values()
    });
    let sheet = sheet(252.0, 200.0, None, false);
    let p = Panels::new(&t, Some(&n), &cut, &plan, sheet, true).unwrap();
    use Quantity::*;
    let (pv, fold, cov) = (p.hist(Pvalue), p.hist(Fold), p.hist(Coverage));
    let (a, b) = (p.pair(Pvalue, Cells), p.pair(Fold, Pvalue));
    assert!(pv.axis.log && !fold.axis.log && !cov.axis.log);
    assert!(cov.axis.hi >= 5000.0, "{:?}", cov.axis);
    assert_eq!((a.x, a.x_bins), (pv.axis, pv.bins));
    assert_eq!((b.y, b.y_bins), (pv.axis, pv.bins));
    assert_eq!((b.x, b.x_bins), (fold.axis, fold.bins));
    let c = p.pair(Coverage, Pvalue);
    assert_eq!((c.x, c.x_bins), (cov.axis, cov.bins));
    assert_eq!(a.marginals, plan.marginals);
    assert_eq!(a.x.ticks(), pv.axis.ticks());
    // `faba qc`'s 2D view reads the same rulers for them.
    let criteria: Vec<Criterion> = Criterion::ALL.to_vec();
    let cache = crate::qc::pairs::Cache::new(plan.scales, plan.bins_on(&sheet));
    let mut fails = Vec::new();
    let s = qc_sites(&t, &n, &criteria, &mut fails, &cut.filter, &cache);
    for (q, h) in [(Pvalue, &pv), (Fold, &fold), (Coverage, &cov)] {
        let r = s.ruler(q, &cut.filter);
        assert_eq!((r.axis, r.bins), (h.axis, h.bins), "{q:?}");
    }
}

/// A count's threshold on a linear scale is drawn between the bar of the
/// count below it and the bar of its own (bars are centred on counts),
/// in a histogram, a 2D panel and `faba qc`'s 2D view alike; on a log
/// scale, at the count.
#[test]
fn a_counts_threshold_sits_between_its_bars() {
    let rows: Vec<Row> = (2..=30).map(|cov| (0.01, cov, 1, 10, 0)).collect();
    let t = table(&rows);
    let n: Vec<usize> = (0..t.len()).collect();
    let cut = Thresholds {
        site_min_coverage: Some(10),
        ..Thresholds::default()
    }
    .cut();
    let plan = args(&["coverage"], &["coverage:pv"], &[], &["coverage=linear"])
        .plan()
        .unwrap();
    let sheet = sheet(252.0, 180.0, None, false);
    let p = Panels::new(&t, Some(&n), &cut, &plan, sheet, true).unwrap();
    let h = p.hist(Quantity::Coverage);
    assert_eq!(h.edges, vec![9.5]);
    assert_eq!(
        p.pair(Quantity::Coverage, Quantity::Pvalue).x_edges,
        h.edges
    );
    // Bars one count wide from 0, centred on their counts: the edge is
    // where bar 9 ends, bar 10 begins.
    let width = (h.axis.hi - h.axis.lo) / h.bins as f64;
    assert_eq!((h.axis.lo, width), (-0.5, 1.0));
    // Drawn there: the dashed line on the boundary of the two bars.
    let look = data_beans::figure::panel::Look::default();
    let (w, ht) = (252.0, 180.0);
    let svg = h.svg(&look, w, ht);
    let (px, _, pw, _) = h.plot_area(&look, (0.0, 0.0, w, ht));
    let at = px + h.axis.frac(9.5) * pw;
    let dashed = (svg.lines())
        .find(|l| l.contains("stroke-dasharray"))
        .expect("a dashed edge");
    let x1: f64 = dashed
        .split("x1=\"")
        .nth(1)
        .unwrap()
        .split('"')
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!((x1 - at).abs() < 0.01, "{dashed} not at {at}");
    assert!(
        (px + 10.0 * pw / h.bins as f64 - x1).abs() < 0.01,
        "{dashed}"
    );
    // `faba qc`'s 2D view puts its box there too.
    let criteria: Vec<Criterion> = Criterion::ALL.to_vec();
    let cache = crate::qc::pairs::Cache::new(plan.scales, plan.bins_on(&sheet));
    let mut fails = Vec::new();
    let s = qc_sites(&t, &n, &criteria, &mut fails, &cut.filter, &cache);
    assert_eq!(
        Quantity::Coverage.passes(&cut.filter, &s),
        (Some(9.5), None)
    );
    // On a log scale, at the count.
    let plan = args(&["coverage"], &["none"], &[], &[]).plan().unwrap();
    let p = Panels::new(&t, Some(&n), &cut, &plan, sheet, true).unwrap();
    assert_eq!(p.hist(Quantity::Coverage).edges, vec![10.0]);
}

/// A clean control's infinite fold and a p-value of 0 are drawn at their
/// axis' end, counted in its last bin: no site is lost.
#[test]
fn an_infinite_fold_and_a_zero_p_land_at_the_end() {
    let t = table(&[
        (0.0, 20, 8, 20, 0),
        (0.01, 20, 4, 20, 2),
        (0.2, 30, 3, 30, 2),
        (0.5, 10, 1, 10, 1),
    ]);
    let plan = args(&["pv", "fold"], &["none"], &[], &[]).plan().unwrap();
    let p = Panels::new(
        &t,
        None,
        &Cut::none(),
        &plan,
        sheet(252.0, 180.0, None, false),
        true,
    )
    .unwrap();
    for q in [Quantity::Pvalue, Quantity::Fold] {
        let h = p.hist(q);
        assert_eq!(h.values[0], h.axis.hi as f32, "{q:?}");
        let st = data_beans::figure::panel::stacks(&h.axis, h.bins, h.values, h.kept);
        assert_eq!(st.iter().map(|(k, r)| k + r).sum::<usize>(), 4, "{q:?}");
        assert!(st[h.bins - 1].0 >= 1, "{q:?}");
    }
}

/// With no threshold given nothing is removed: no mask, no edge, no key.
#[test]
fn no_threshold_removes_nothing() {
    let t = many(300);
    let plan = args(&["pv"], &["pv:coverage"], &[], &[]).plan().unwrap();
    let p = Panels::new(
        &t,
        None,
        &Cut::none(),
        &plan,
        sheet(252.0, 200.0, None, false),
        true,
    )
    .unwrap();
    let h = p.hist(Quantity::Pvalue);
    assert!(h.kept.is_none() && h.edges.is_empty());
    let d = p.pair(Quantity::Pvalue, Quantity::Coverage);
    assert!(d.kept.is_none() && d.x_edges.is_empty() && d.y_edges.is_empty());
    for f in p.figures() {
        assert!(!f.svg.contains(">removed<"), "{}", f.stem);
    }
}

/// Thresholds given cut as `faba qc` would: a site is kept when none of
/// them fails it, and each is drawn as an edge on its quantity's axis.
#[test]
fn thresholds_split_the_sites_as_their_criteria_fail() {
    let t = many(500);
    let n = cells(500);
    let given = Thresholds {
        site_max_pv: Some(0.05),
        site_min_coverage: Some(10),
        site_min_cells: Some(20),
        ..Thresholds::default()
    };
    let cut = given.cut();
    let plan = args(&["pv", "coverage", "fold"], &["pv:cells"], &[], &[])
        .plan()
        .unwrap();
    let p = Panels::new(
        &t,
        Some(&n),
        &cut,
        &plan,
        sheet(252.0, 200.0, None, false),
        true,
    )
    .unwrap();
    let mask = p.hist(Quantity::Pvalue).kept.expect("a mask");
    let on = [
        Criterion::MaxPv,
        Criterion::MinCoverage,
        Criterion::MinCells,
    ];
    for i in 0..t.len() {
        let fails = on.iter().any(|c| c.fails(&cut.filter, &t, i, Some(n[i])));
        assert_eq!(mask[i], !fails, "site {i}");
    }
    assert!(mask.iter().any(|k| *k) && mask.iter().any(|k| !k));
    assert_eq!(
        p.hist(Quantity::Pvalue).edges,
        vec![-(0.05f32 as f64).log10()]
    );
    assert_eq!(p.hist(Quantity::Coverage).edges, vec![10.0]);
    assert!(p.hist(Quantity::Fold).edges.is_empty());
    assert_eq!(
        p.pair(Quantity::Pvalue, Quantity::Cells).y_edges,
        vec![20.0]
    );
}

/// Auto's pairs are `faba qc`'s: those it picks at its default thresholds,
/// with any given over them.
#[test]
fn auto_pairs_are_qcs_at_its_defaults() {
    let t = many(900);
    let n = cells(900);
    let plan = args(&[], &[], &[], &[]).plan().unwrap();
    let p = Panels::new(
        &t,
        Some(&n),
        &Cut::none(),
        &plan,
        sheet(252.0, 200.0, None, false),
        true,
    )
    .unwrap();
    let f = SiteFilterArgs::default_values();
    let fails: Vec<u8> = (0..t.len())
        .map(|i| {
            Criterion::ALL
                .iter()
                .filter(|c| c.fails(&f, &t, i, Some(n[i])))
                .fold(0, |b, &c| b | 1 << c as u8)
        })
        .collect();
    let cache = crate::qc::pairs::Cache::default();
    let criteria = Criterion::ALL.to_vec();
    let s = crate::qc::pairs::Sites {
        table: &t,
        n_cells: Some(&n),
        fails: &fails,
        criteria: &criteria,
        cache: &cache,
    };
    let auto: Vec<_> = (crate::qc::pairs::auto_pairs(&s, crate::qc::pairs::PANELS).into_iter())
        .map(|p| (p.x, p.y))
        .collect();
    assert_eq!(p.pairs(), auto.as_slice());
    assert_eq!(
        p.hists(),
        crate::qc::pairs::axes(&s).as_slice(),
        "every quantity"
    );
    // Without cells, none of cells.
    let p = Panels::new(
        &t,
        None,
        &Cut::none(),
        &plan,
        sheet(252.0, 200.0, None, false),
        true,
    )
    .unwrap();
    assert!(!p.hists().contains(&Quantity::Cells));
    let asked = args(&["cells"], &[], &[], &[]).plan().unwrap();
    let e = Panels::new(
        &t,
        None,
        &Cut::none(),
        &asked,
        sheet(252.0, 200.0, None, false),
        true,
    );
    assert!(e.err().unwrap().to_string().contains("cells"));
}

/// Every figure is the size asked for exactly, every label at the font
/// asked; with a title, it alone a little larger. Named by modality and
/// quantities.
#[test]
fn figures_are_the_size_asked_with_the_font_on_every_label() {
    let t = many(400);
    let n = cells(400);
    let plan = args(&["pv", "cells"], &["pv:cells"], &["top", "right"], &[])
        .plan()
        .unwrap();
    let (w, h) = (1.8 * 72.0, 1.4 * 72.0);
    for alone in [true, false] {
        let p = Panels::new(
            &t,
            Some(&n),
            &Cut::none(),
            &plan,
            sheet(w, h, Some(6.0), alone),
            true,
        )
        .unwrap();
        let figures = p.figures();
        let stems: Vec<&str> = figures.iter().map(|f| f.stem.as_str()).collect();
        assert_eq!(stems, ["m6a_hist_pv", "m6a_hist_cells", "m6a_cells_vs_pv"]);
        for f in &figures {
            let head = f.svg.lines().next().unwrap();
            assert!(head.contains(r#"width="129.6" height="100.8""#), "{head}");
            let labels = crate::figure::labels(&f.svg);
            assert!(!labels.is_empty());
            for (size, text, _) in &labels {
                let title = !alone && text.starts_with("m6a");
                let want = if title { "7.5" } else { "6" };
                assert_eq!(size, want, "{}: {text}", f.stem);
            }
            assert_eq!(labels.iter().any(|l| l.1.starts_with("m6a")), !alone);
        }
    }
}

/// Write `{dir}/{name}.zarr` from triplets.
fn write_matrix(
    dir: &std::path::Path,
    name: &str,
    rows: &[&str],
    cols: &[&str],
    e: &[(u64, u64, f32)],
) {
    use data_beans::sparse_io::{create_sparse_from_triplets_owned, SparseIoBackend};
    let path = dir
        .join(format!("{name}.zarr"))
        .to_string_lossy()
        .into_owned();
    let mut data = create_sparse_from_triplets_owned(
        e.to_vec(),
        (rows.len(), cols.len(), e.len()),
        Some(&path),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let names = |v: &[&str]| -> Vec<Box<str>> { v.iter().map(|s| (*s).into()).collect() };
    data.register_row_names_vec(&names(rows));
    data.register_column_names_vec(&names(cols));
}

/// Cells per site are the cells with a converted read there, over every
/// cell (none cut) and every batch's `_site` matrix.
#[test]
fn cells_per_site_count_every_cell_with_a_converted_read() {
    let dir = tempfile::tempdir().unwrap();
    let rows = [
        "GENE1/m6a/chr1:100/methylated",
        "GENE1/m6a/chr1:100/unmethylated",
        "GENE2/m6a/chr1:200/methylated",
        "GENE2/m6a/chr1:200/unmethylated",
    ];
    let cols = ["C1", "C2", "C3", "C4"];
    // Site 1 converted in C1, C2, C4; site 2 only unconverted reads but C3.
    let b1 = [
        (0, 0, 1.0),
        (0, 1, 2.0),
        (0, 3, 1.0),
        (1, 2, 5.0),
        (3, 0, 4.0),
        (2, 2, 1.0),
    ];
    write_matrix(dir.path(), "b1_m6a_site", &rows, &cols, &b1);
    write_matrix(
        dir.path(),
        "b2_m6a_site",
        &rows,
        &cols,
        &[(0, 1, 1.0), (1, 1, 1.0)],
    );
    let layout = crate::qc::layout::scan_input_dir(&dir.path().to_string_lossy()).unwrap();
    let t = table(&[(0.01, 10, 3, 10, 0), (0.5, 10, 1, 10, 1), (0.5, 5, 1, 5, 1)]);
    let n = site_cells(&layout, &t).unwrap().expect("site matrices");
    assert_eq!(n, vec![4, 1, 0]);
    let none = tempfile::tempdir().unwrap();
    std::fs::write(none.path().join("m6a_sites.parquet"), b"").unwrap();
    let layout = crate::qc::layout::scan_input_dir(&none.path().to_string_lossy()).unwrap();
    assert_eq!(site_cells(&layout, &t).unwrap(), None);
}
