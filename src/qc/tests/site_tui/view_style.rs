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
                // Its x title, the one thing apart, is placed to fit as long
                // as it is: the rest is the view's.
                let kept = "kept cells with a converted read";
                let all = "cells with a converted read";
                let untitled = |svg: &str| {
                    (svg.lines())
                        .filter(|l| !l.contains(all))
                        .collect::<Vec<_>>()
                        .join("\n")
                };
                assert_eq!(
                    untitled(&theirs[at].svg),
                    untitled(&view),
                    "{stem} {palette:?} {alone}"
                );
                let counted = c == Criterion::MinCells;
                assert_eq!(view.contains(kept), counted, "{stem}");
                assert!(!theirs[at].svg.contains(kept), "{stem}");
                assert_eq!(theirs[at].svg.contains(all), counted);
            }
            let g = crate::qc::pairs::grid(
                &p.view().sites(),
                (Axis::Pvalue, Axis::Cells),
                &f,
                crate::qc::pairs::BINS,
            );
            let title = if alone { "" } else { "cells by -log10 p" };
            let how = crate::qc::pairs::PairLook {
                arms: Default::default(),
                titles: Some(&plan.titles),
            };
            let none = (Default::default(), 0.0);
            let view = drawn(sheet, |canvas, b| g.draw_with(canvas, b, title, none, &how));
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

/// Made-up m6A sites with a control arm, spread as real ones are: most
/// near-uniform, a third with a strong signal; their cells per site.
fn made_up(n: usize) -> (SiteTable, Vec<usize>) {
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut cells = Vec::with_capacity(n);
    let rows: Vec<crate::qc::fixtures::Row> = (0..n)
        .map(|_| {
            let strong = next() < 0.35;
            let cov = (3.0 * (next() * 4.5).exp()) as u64 + 1;
            let ctl_cov = (2.0 * (next() * 4.2).exp()) as u64;
            let rate = if strong {
                0.1 + 0.5 * next()
            } else {
                0.02 * next()
            };
            let conv = ((cov as f64) * rate * (0.5 + next())) as u64 % (cov + 1);
            let ctl_conv = ((ctl_cov as f64) * 0.02 * next() * 2.0) as u64 % (ctl_cov + 1);
            let pv = match strong {
                true => 10f64.powf(-1.0 - 25.0 * next() * next()) as f32,
                false => next().max(1e-3) as f32,
            };
            let pv = if next() < 0.01 { 0.0 } else { pv };
            cells.push((conv as f64 * (0.3 + next())) as usize);
            (pv, cov, conv, ctl_cov, ctl_conv)
        })
        .collect();
    (
        crate::qc::fixtures::site_table(M6A, 50, &rows, |i| i),
        cells,
    )
}

/// Mock-ups of the view style, from made-up m6A sites with a control arm,
/// 2.4 by 1.8 inches at 7 points, in gray: histograms of pv, cells,
/// coverage and the edit ratio, two pairs with short marginals, one with
/// the control arm in front; and the view on screen at 160 by 44. Written
/// where `FABA_MOCKUPS_DIR` names, the screen where `FABA_SCREENSHOTS_DIR`
/// does.
#[test]
#[ignore = "writes mock-ups to look at; run by hand"]
fn view_style_mockups() {
    let Some(dir) = std::env::var_os("FABA_MOCKUPS_DIR") else {
        eprintln!("FABA_MOCKUPS_DIR is not set");
        return;
    };
    let dir = std::path::PathBuf::from(dir);
    std::fs::create_dir_all(&dir).unwrap();
    let (t, n) = made_up(6000);
    let cut = crate::qc::site_plot::Thresholds {
        site_max_pv: Some(0.05),
        site_min_cells: Some(10),
        site_min_coverage: Some(10),
        site_min_edit_ratio: Some(0.05),
        ..Default::default()
    }
    .cut();
    let sheet = sheet(Palette::Gray, true);
    let draw = |args: PanelArgs, suffix: &str| {
        let plan = args.plan().unwrap();
        let figures = Panels::new(&t, Some(&n), &cut, &plan, sheet, true)
            .unwrap()
            .figures();
        for f in figures {
            let prefix = dir.join(format!("{}{suffix}", f.stem));
            crate::figure::save_with_png(&f.svg, &prefix.to_string_lossy()).unwrap();
        }
    };
    let view = |hist: &[&str], pairs: &[&str]| PanelArgs {
        hist: hist.iter().map(|s| s.to_string()).collect(),
        pairs: pairs.iter().map(|s| s.to_string()).collect(),
        style: PanelStyle::View,
        marginals: vec!["both".into()],
        marginal_size: Some(0.12),
        key: true,
        ..PanelArgs::default()
    };
    draw(
        view(
            &["pv", "cells", "coverage", "edit_ratio"],
            &["pv:cells", "coverage:edit_ratio"],
        ),
        "",
    );
    let front = PanelArgs {
        front: Some(crate::qc::arms::Front::Mut),
        ..view(&["coverage"], &["coverage:edit_ratio"])
    };
    draw(front, "_front_mut");
    let mut p = picker(&t, Some(n.clone()), cut.filter.clone());
    focus_on(&mut p, Criterion::MinCoverage);
    let mut term = Terminal::new(TestBackend::new(160, 44)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    // Where `FABA_SCREENSHOTS_DIR` names (the same folder, set by hand).
    crate::tui::screenshot::save("qc_view_screen", term.backend().buffer());
    let save = |stem: &str, svg: &str| {
        let prefix = dir.join(stem);
        crate::figure::save_with_png(svg, &prefix.to_string_lossy()).unwrap();
    };
    save("qc_figure", &p.figure(Sheet::OWN).unwrap());
    save("qc_pairs_figure", &p.pairs_figure(0, Sheet::OWN));
    press(&mut p, KeyCode::Char('v'));
    term.draw(|f| p.render(f)).unwrap();
    crate::tui::screenshot::save("qc_view_screen_2d", term.backend().buffer());
}
