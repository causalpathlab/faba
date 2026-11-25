//! `--style view`: the site panels drawn as the view draws its own, one
//! drawing for both: the focused knob's histogram and a pair's squares.

use super::*;
use crate::figure::style::{Palette, StyleArgs};
use crate::qc::site_plot::{Cut, PanelArgs, PanelStyle, Panels};

/// A 2.4 by 1.8 inch panel at 7 points, in `palette`, alone or not.
fn sheet(palette: Palette, alone: bool) -> Sheet {
    Sheet {
        page: (Some(2.4 * 72.0), Some(1.8 * 72.0)),
        look: StyleArgs {
            font_size: Some(7.0),
            panel_only: alone,
            palette,
        },
    }
}

/// What the view draws into a figure `sheet`'s size: `draw` on its page.
fn drawn(sheet: Sheet, draw: impl FnOnce(&mut Canvas, (f64, f64, f64, f64))) -> String {
    let (w, h) = sheet.fill(|| unreachable!("both sides asked"));
    let mut c = Canvas::styled(w, h, &sheet.style((w, h)));
    let (pw, ph) = c.page();
    draw(&mut c, (0.0, 0.0, pw, ph));
    c.finish()
}

/// With `--style view` each histogram is the view's own of its knob, at
/// the thresholds now (its bars, their colours, the threshold's line, its
/// title, ticks and axis titles), and the pair the view's 2D panel of it;
/// `faba site-plot` at the same cut draws them alike. In either palette,
/// alone or not.
#[test]
fn the_view_style_draws_what_the_view_draws() {
    let t = table(M6A, 600);
    let n = cells(600);
    let f = SiteFilterArgs {
        site_min_cells: 10,
        ..SiteFilterArgs::default_values()
    };
    let mut p = picker(&t, Some(n.clone()), f.clone());
    let plan = PanelArgs {
        hist: vec!["pv".into(), "coverage".into(), "cells".into()],
        pairs: vec!["pv:cells".into()],
        style: PanelStyle::View,
        ..PanelArgs::default()
    }
    .plan()
    .unwrap();
    for palette in [Palette::Colour, Palette::Gray] {
        for alone in [false, true] {
            let sheet = sheet(palette, alone);
            p.set_scales(plan.scales, plan.bins_on(&sheet));
            let ours = Panels::of_view(&p.view().sites(), &p.filter, &plan, sheet)
                .unwrap()
                .figures();
            let theirs = Panels::new(&t, Some(&n), &Cut::of_qc(&f), &plan, sheet, true)
                .unwrap()
                .figures();
            assert_eq!(ours.len(), 4);
            for c in [
                Criterion::MaxPv,
                Criterion::MinCoverage,
                Criterion::MinCells,
            ] {
                focus_on(&mut p, c);
                let title = format!("{} · y: sites per bin", c.label());
                let view = drawn(sheet, |canvas, b| p.draw_hist(canvas, b, title));
                let stem = format!("m6a_hist_{}", Axis::of(c).key());
                let at = ours.iter().position(|g| g.stem == stem).unwrap();
                assert_eq!(ours[at].svg, view, "{stem} {palette:?} {alone}");
                // `faba site-plot` counts the cells over every cell (no
                // cell QC): its cells are not called kept.
                let kept = "kept cells with a converted read";
                let over_all = view.replace(kept, "cells with a converted read");
                assert_eq!(theirs[at].svg, over_all, "{stem} {palette:?} {alone}");
                let counted = c == Criterion::MinCells;
                assert_eq!(view.contains(kept), counted, "{stem}");
                assert!(!theirs[at].svg.contains(kept), "{stem}");
                assert_eq!(over_all.contains("cells with a converted read"), counted);
            }
            let g = crate::qc::pairs::grid(
                &p.view().sites(),
                (Axis::Pvalue, Axis::Cells),
                &f,
                crate::qc::pairs::BINS,
            );
            let title = if alone { "" } else { "cells by -log10 p" };
            let view = drawn(sheet, |canvas, b| g.draw(canvas, b, title));
            assert_eq!(ours[3].stem, "m6a_cells_vs_pv");
            assert_eq!(ours[3].svg, view, "pair {palette:?} {alone}");
            assert_eq!(theirs[3].svg, view, "pair {palette:?} {alone}");
        }
    }
}

/// `s` in `faba qc-view` saves the panel on screen in the style asked:
/// with `--style view`, the view's own histogram of the focused knob.
#[test]
fn s_saves_the_view_style_when_asked() {
    let t = table(M6A, 400);
    let mut p = picker(&t, Some(cells(400)), SiteFilterArgs::default_values());
    let sheet = sheet(Palette::Colour, false);
    let plan = PanelArgs {
        style: PanelStyle::View,
        ..PanelArgs::default()
    }
    .plan()
    .unwrap();
    p.set_scales(plan.scales, plan.bins_on(&sheet));
    focus_on(&mut p, Criterion::MinCoverage);
    let saved = p.screen_panel(&plan, sheet).unwrap();
    let title = format!("{} · y: sites per bin", Criterion::MinCoverage.label());
    assert_eq!(saved, drawn(sheet, |c, b| p.draw_hist(c, b, title)));
}

/// `--bins` sets how many bins the view style's histograms take, about
/// that many on every knob's scale, fewer asked fewer drawn; none given,
/// the view's own.
#[test]
fn bins_set_the_view_styles_histograms() {
    let t = table(M6A, 600);
    let n = cells(600);
    let p = picker(&t, Some(n.clone()), SiteFilterArgs::default_values());
    let s = p.view().sites();
    for &c in &s.criteria.to_vec() {
        for &scale in c.scales() {
            let own = Column::new(&s, c, scale).hist.counts;
            assert_eq!(Column::binned(&s, c, scale, None).hist.counts, own);
            let len = |k: usize| Column::binned(&s, c, scale, Some(k)).hist.counts.len();
            let (few, more) = (len(8), len(24));
            assert!(few <= 9 && more <= 25, "{c:?} {scale:?}: {few} {more}");
            assert!(few < more, "{c:?} {scale:?}: {few} {more}");
        }
    }
    let figures = |bins: Option<usize>| {
        let plan = PanelArgs {
            hist: vec!["pv".into(), "coverage".into(), "cells".into()],
            pairs: vec!["none".into()],
            style: PanelStyle::View,
            bins,
            ..PanelArgs::default()
        }
        .plan()
        .unwrap();
        let sheet = sheet(Palette::Colour, false);
        let cut = Cut::of_qc(&SiteFilterArgs::default_values());
        Panels::new(&t, Some(&n), &cut, &plan, sheet, true)
            .unwrap()
            .figures()
    };
    let (own, eight) = (figures(None), figures(Some(8)));
    assert_eq!(own.len(), 3);
    for (a, b) in own.iter().zip(&eight) {
        assert_eq!(a.stem, b.stem);
        assert_ne!(a.svg, b.svg, "{}: --bins 8 drawn as the view's own", a.stem);
    }
}
