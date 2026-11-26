//! Every `--style view` panel is the view's own drawing of the same sites
//! at the same settings; the p-value's bars sit on log-spaced p under a
//! round cap; each count axis reads on the scale asked, ticked at round
//! counts where they sit; the key goes where asked; and each arm of the
//! p-value and the cells is drawn without changing any cut.

use super::*;
use crate::figure::style::{Palette, StyleArgs};
use crate::qc::arms::{self, CountScale, KeyAt, PAxis};
use crate::qc::site_plot::{Cut, PanelArgs, PanelStyle, Panels};

/// A table of `n` sites with each arm's cells, the control's its own.
fn with_arms(n: usize) -> SiteTable {
    let mut t = table(M6A, n);
    t.control_cells = Some((0..n).map(|i| i * 7 % 23).collect());
    t
}

/// A `w` by 1.8 inch panel at 7 points in `palette`, alone or not.
fn sheet(w: f64, palette: Palette, alone: bool) -> Sheet {
    Sheet {
        page: (Some(w * 72.0), Some(1.8 * 72.0)),
        look: StyleArgs {
            font_size: Some(7.0),
            panel_only: alone,
            palette,
        },
    }
}

/// The focused knob's histogram as the view draws it on `sheet`'s page.
fn view_drawn(p: &SitePicker, sheet: Sheet) -> String {
    let (w, h) = sheet.fill(|| unreachable!("both sides asked"));
    let mut c = Canvas::styled(w, h, &sheet.style((w, h)));
    let (pw, ph) = c.page();
    let title = format!("{} · y: sites per bin", p.criterion().label());
    p.draw_hist(&mut c, (0.0, 0.0, pw, ph), title);
    c.finish()
}

/// `args` in the view style, read.
fn plan(args: PanelArgs) -> crate::qc::site_plot::Plan {
    PanelArgs {
        style: PanelStyle::View,
        ..args
    }
    .plan()
    .unwrap()
}

/// For every knob, `faba site-plot`'s panel and `faba qc`'s (and so `qc
/// --plots`, `qc-view --plots` and the MCP tools, which draw through
/// them) are byte for byte the view's own drawing of it, at the same
/// thresholds, scales, arms and key, in either palette, alone or not. The
/// one thing apart: `faba site-plot` counts the cells over every cell (no
/// cell QC), so its cells' axis does not say "kept".
#[test]
fn every_knob_s_panel_is_the_view_s_drawing() {
    let t = with_arms(900);
    let n = cells(900);
    let f = SiteFilterArgs {
        site_min_cells: 10,
        ..SiteFilterArgs::default_values()
    };
    let mut p = picker(&t, Some(n.clone()), f.clone());
    let every: Vec<String> = Axis::ALL.iter().map(|a| a.key().to_string()).collect();
    for (palette, alone, key_at) in [
        (Palette::Gray, true, KeyAt::Right),
        (Palette::Gray, true, KeyAt::Inside),
        (Palette::Gray, true, KeyAt::Top),
        (Palette::Colour, false, KeyAt::Right),
        (Palette::Gray, false, KeyAt::Alone),
    ] {
        let plan = plan(PanelArgs {
            hist: every.clone(),
            pairs: vec!["none".into()],
            key: alone,
            key_at: Some(key_at),
            yscale: vec!["log1p".into(), "pv=sqrt".into()],
            ..PanelArgs::default()
        });
        let sheet = sheet(2.2, palette, alone);
        p.set_scales(plan.scales, plan.bins_on(&sheet));
        p.arm_look = plan.look;
        let ours = Panels::of_view(&p.view().sites(), &p.filter, &plan, sheet)
            .unwrap()
            .figures();
        let theirs = Panels::new(&t, Some(&n), &Cut::of_qc(&f), &plan, sheet, true)
            .unwrap()
            .figures();
        let knobs = p.view().criteria.clone();
        for q in Axis::ALL {
            // A quantity of two knobs draws the one on (of the edit ratio's
            // two, the lower).
            let on = |c: &&Criterion| knobs.contains(c) && !c.is_off(&p.filter);
            let first = (q.criteria().iter()).find(|c| knobs.contains(c));
            let Some(&c) = (q.criteria().iter()).find(on).or(first) else {
                continue;
            };
            focus_on(&mut p, c);
            let view = view_drawn(&p, sheet);
            let stem = format!("m6a_hist_{}", q.key());
            let of = |figs: &[Figure]| figs.iter().find(|g| g.stem == stem).unwrap().svg.clone();
            assert_eq!(of(&ours), view, "{stem} {palette:?} {alone} {key_at:?}");
            let untitled = |svg: &str| {
                (svg.lines())
                    .filter(|l| !l.contains("cells with a converted read"))
                    .collect::<Vec<_>>()
                    .join("\n")
            };
            assert_eq!(untitled(&of(&theirs)), untitled(&view), "{stem}");
            let keyed = key_at == KeyAt::Alone;
            let key_stem = format!("{stem}_key");
            assert_eq!(ours.iter().any(|g| g.stem == key_stem), keyed, "{key_stem}");
        }
    }
}

/// The p-value's last bar caps at a round p, the greatest of 1e-10,
/// 1e-20, 1e-30, 1e-50 ... the smallest p reaches, and holds every p at
/// or below it: no edge past the cap, the cap's tick "<= p" of it, and
/// every other tick a round p at the place that p sits.
#[test]
fn p_bins_cap_at_a_round_p_their_ticks_naming_their_place() {
    let mut t = with_arms(900);
    // Deep sites, so each arm's p reaches far below the cap: the
    // signal's on the cut's test, the control's on its own reads.
    for i in (0..t.len()).step_by(9) {
        let deep = 400 + i as u64;
        (t.control_coverage[i], t.control_converted[i]) = (deep, deep);
        t.pv[i] = 1e-36;
    }
    let mut p = picker(&t, Some(cells(900)), SiteFilterArgs::default_values());
    focus_on(&mut p, Criterion::MaxPv);
    let h = &p.arm_hist;
    let cap = h.clamped.0.expect("p-values past the cap");
    assert!([10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0].contains(&cap));
    assert_eq!(
        h.edges[h.edges.len() - 1],
        cap,
        "the last bar ends at the cap"
    );
    assert!(h.edges.iter().all(|&e| e <= cap), "{:?}", h.edges);
    let scale = p.scale();
    assert_eq!(scale, Scale::Log, "log-spaced p unless a scale is asked");
    for (at, said, end) in h.tick_choices(scale, PAxis::P) {
        let d = at.exp_m1();
        match end {
            true => assert_eq!(said, format!("<= {}", arms::fmt_p(10f64.powf(-cap)))),
            false => {
                assert_eq!(said, arms::fmt_p(10f64.powf(-d)), "at -log10 p {d}");
                // The p's leading figure: -d past its power of ten.
                let mantissa = 10f64.powf(-d - (-d + 1e-6).floor());
                assert!(
                    [1.0, 2.0, 5.0].iter().any(|m| (mantissa - m).abs() < 1e-3),
                    "{said} is a round p"
                );
            }
        }
    }
    // Every site is in a bar: those past the cap in the last.
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    assert_eq!(sum(&h.kept) + sum(&h.filtered), t.len());
}

/// The cap: the greatest round `-log10 p` the values reach, raised past
/// a threshold beyond it; none under the least.
#[test]
fn the_p_cap_is_round_and_past_every_threshold() {
    let cap = |hi: f64, cuts: &[f64]| arms::p_cap(hi, cuts.iter().copied());
    assert_eq!(cap(44.85, &[]), Some(30.0));
    assert_eq!(cap(9.0, &[]), None);
    assert_eq!(cap(300.0, &[2.0]), Some(300.0));
    assert_eq!(cap(44.85, &[35.0]), Some(50.0));
}

/// A count axis' ticks are round counts at the place each sits, on every
/// scale: 0, about the middle, and the least round count at or above the
/// most, where the axis ends.
#[test]
fn count_ticks_are_round_counts_where_they_sit() {
    let round = |v: f64| {
        let m = v / 10f64.powf(v.log10().floor());
        [1.0, 2.0, 2.5, 5.0].iter().any(|r| (m - r).abs() < 1e-9)
    };
    for scale in [CountScale::Linear, CountScale::Sqrt, CountScale::Log1p] {
        for most in [1.0, 7.0, 95_000.0, 139_087.0, 2.4e6] {
            let ticks = scale.ticks(most);
            assert_eq!(ticks[0], 0.0);
            assert!(ticks.windows(2).all(|w| w[0] < w[1]), "{scale:?} {ticks:?}");
            let top = ticks[ticks.len() - 1];
            assert!(top >= most, "{scale:?}: {top} under {most}");
            assert!(round(top), "{scale:?}: {top} is round");
            // The least such: no round count between the most and it.
            let lower = [0.5, 0.4, 0.8].map(|r| top * r);
            assert!(
                lower.iter().all(|&l| !round(l) || l < most),
                "{scale:?}: {top} for {most}"
            );
            for &v in &ticks {
                let back = scale.count(scale.of(v));
                assert!((back - v).abs() <= 1e-6 * v.max(1.0), "{scale:?} {v}");
                assert!(v < 1.0 || round(v), "{scale:?}: {v}");
            }
            // Deeper as it counts more.
            assert!((0..50).all(|n| scale.of(n as f64) < scale.of(n as f64 + 1.0)));
        }
    }
    // Sqrt's middle: a quarter of the top, rounded.
    assert_eq!(
        CountScale::Sqrt.ticks(139_087.0),
        vec![0.0, 50_000.0, 200_000.0]
    );
}

/// The count axis reaches over the tallest bar: its top tick, a round
/// count at or above it, is the highest gridline, and no bar runs past it.
#[test]
fn the_tallest_bar_stays_under_the_top_tick() {
    let t = with_arms(900);
    let mut p = picker(&t, Some(cells(900)), SiteFilterArgs::default_values());
    let sheet = sheet(2.2, Palette::Gray, true);
    let attr = |l: &str, k: &str| -> Option<f64> {
        l.split(&format!(" {k}=\""))
            .nth(1)?
            .split('"')
            .next()?
            .parse()
            .ok()
    };
    for c in [
        Criterion::MaxPv,
        Criterion::MinCoverage,
        Criterion::MinCells,
    ] {
        focus_on(&mut p, c);
        for scale in [CountScale::Linear, CountScale::Sqrt, CountScale::Log1p] {
            p.arm_look.y.set(Axis::of(c), scale);
            let svg = view_drawn(&p, sheet);
            let (page_w, _) = sheet.fill(|| unreachable!());
            let top_line = (svg.lines())
                .filter(|l| l.starts_with("<line") && !l.contains("dasharray"))
                .filter_map(|l| {
                    let (y1, y2) = (attr(l, "y1")?, attr(l, "y2")?);
                    (y1 == y2).then_some(y1)
                })
                .fold(f64::INFINITY, f64::min);
            let highest_bar = (svg.lines())
                .filter(|l| l.starts_with("<rect"))
                .filter(|l| attr(l, "width").is_some_and(|w| w < page_w - 1.0))
                .filter_map(|l| attr(l, "y"))
                .fold(f64::INFINITY, f64::min);
            assert!(
                highest_bar >= top_line - 0.015,
                "{c:?} {scale:?}: a bar at {highest_bar} over the top tick at {top_line}"
            );
            let top = scale.ticks(p.arm_hist.tallest() as f64);
            let said = crate::tui::ui::compact(top[top.len() - 1]);
            assert!(
                crate::figure::labels(&svg).iter().any(|l| l.1 == said),
                "{c:?} {scale:?}: {said}"
            );
        }
    }
}

/// Each label of the y axis is the count at its height: its baseline
/// where that count's bar would reach, on every scale; switching the
/// scale moves only y, the bars' places across and the threshold's line
/// as they were.
#[test]
fn the_y_scale_moves_only_y_and_its_labels_name_their_height() {
    let t = with_arms(900);
    let mut p = picker(&t, Some(cells(900)), SiteFilterArgs::default_values());
    focus_on(&mut p, Criterion::MinCoverage);
    let sheet = sheet(3.0, Palette::Gray, true);
    let across = |svg: &str| -> Vec<(String, String)> {
        (svg.lines())
            .filter(|l| l.starts_with("<rect") || l.contains("stroke-dasharray"))
            .map(|l| {
                let attr = |k: &str| {
                    let key = format!(" {k}=\"");
                    l.split(&key)
                        .nth(1)
                        .map_or("", |r| r.split('"').next().unwrap())
                        .to_string()
                };
                (attr("x"), attr("width"))
            })
            .collect()
    };
    let mut seen = Vec::new();
    for scale in [CountScale::Linear, CountScale::Sqrt, CountScale::Log1p] {
        p.arm_look.y.set(Axis::Coverage, scale);
        let svg = view_drawn(&p, sheet);
        seen.push(svg.clone());
        let most = p.arm_hist.tallest() as f64;
        let said: Vec<String> = crate::figure::labels(&svg)
            .into_iter()
            .map(|l| l.1)
            .collect();
        for v in scale.ticks(most) {
            let text = crate::tui::ui::compact(v);
            assert!(said.contains(&text), "{scale:?}: {text} in {said:?}");
        }
        assert_eq!(across(&svg), across(&seen[0]), "{scale:?}: only y moves");
    }
    assert!(seen[0] != seen[1] && seen[1] != seen[2]);
}

/// A threshold's words sit on their own side of its line, the low end's
/// left of it, the high end's right; a side with no room has none, never
/// the other side's place.
#[test]
fn a_threshold_s_words_keep_to_their_side() {
    let t = with_arms(600);
    let f = SiteFilterArgs {
        site_min_cells: 10,
        ..SiteFilterArgs::permissive()
    };
    let mut p = picker(&t, Some(cells(600)), f);
    focus_on(&mut p, Criterion::MinCells);
    let h = &p.arm_hist;
    let e = h.cuts[0].2;
    let width = |s: &str| s.chars().count() as f64;
    for line_at in [3.0, 40.0, 97.0] {
        let x_of = |a: f64| if a == h.at[e] { line_at } else { a };
        let placed = h.place_words(&x_of, (0.0, 100.0), width);
        for w in &placed {
            assert!(!w.text.contains('←') && !w.text.contains('→'), "{w:?}");
            match w.text.starts_with('<') {
                true => assert!(w.from + width(&w.text) <= line_at, "{w:?} at {line_at}"),
                false => assert!(w.from >= line_at, "{w:?} at {line_at}"),
            }
        }
        let n = placed.len();
        assert_eq!(
            n,
            match line_at as i64 {
                3 | 97 => 1,
                _ => 2,
            }
        );
    }
}

/// The plot's ends across: where its base line and gridlines start and
/// end.
fn plot_span(svg: &str) -> (f64, f64) {
    (svg.lines())
        .filter(|l| l.starts_with("<line"))
        .filter_map(|l| {
            let attr = |k: &str| -> Option<f64> {
                l.split(&format!(" {k}=\""))
                    .nth(1)?
                    .split('"')
                    .next()?
                    .parse()
                    .ok()
            };
            let (y1, y2) = (attr("y1")?, attr("y2")?);
            (y1 == y2).then_some((attr("x1")?, attr("x2")?))
        })
        .fold((f64::INFINITY, 0.0), |(a, b), (x1, x2)| {
            (a.min(x1), b.max(x2))
        })
}

/// Whether `svg` draws the key.
fn has_key(svg: &str) -> bool {
    crate::figure::labels(svg)
        .iter()
        .any(|l| l.1 == "retained (WT)" || l.1 == "WT kept")
}

/// `--key-at`: right narrows the plot by the key's width; top and inside
/// keep the plot as wide as with no key, the key above it or over its
/// corner; alone draws none in the panel and the key as a figure of its
/// own.
#[test]
fn the_key_goes_where_asked() {
    let t = with_arms(600);
    let mut p = picker(&t, Some(cells(600)), SiteFilterArgs::default_values());
    focus_on(&mut p, Criterion::MinCoverage);
    let sheet = sheet(4.0, Palette::Gray, true);
    let right_edge = |svg: &str| plot_span(svg).1;
    p.arm_look.key = false;
    let none = view_drawn(&p, sheet);
    assert!(!has_key(&none));
    p.arm_look.key = true;
    p.arm_look.key_at = KeyAt::Right;
    let right = view_drawn(&p, sheet);
    p.arm_look.key_at = KeyAt::Top;
    let top = view_drawn(&p, sheet);
    p.arm_look.key_at = KeyAt::Inside;
    let inside = view_drawn(&p, sheet);
    p.arm_look.key_at = KeyAt::Alone;
    let alone = view_drawn(&p, sheet);
    assert!(has_key(&right) && has_key(&top) && has_key(&inside) && !has_key(&alone));
    assert!(right_edge(&right) < right_edge(&none) - 20.0);
    assert_eq!(right_edge(&top), right_edge(&none));
    assert_eq!(right_edge(&inside), right_edge(&none));
    assert_eq!(alone, none, "alone, the panel is drawn as with no key");
    let plan = plan(PanelArgs {
        hist: vec!["coverage".into()],
        pairs: vec!["none".into()],
        key: true,
        key_at: Some(KeyAt::Alone),
        ..PanelArgs::default()
    });
    let figs = Panels::of_view(&p.view().sites(), &p.filter, &plan, sheet)
        .unwrap()
        .figures();
    let key = figs
        .iter()
        .find(|g| g.stem == "m6a_hist_coverage_key")
        .unwrap();
    for said in ["background (MUT)", "filtered out (WT)", "retained (WT)"] {
        assert!(
            crate::figure::labels(&key.svg).iter().any(|l| l.1 == said),
            "{said}"
        );
    }
}

/// A key at the right never leaves the plot under 60% of the panel's
/// width: on a narrow panel it is drawn smaller, else above the plot, as
/// `--key-at top` draws it; on a wide one beside the plot.
#[test]
fn a_key_at_the_right_leaves_the_plot_most_of_the_panel() {
    let t = with_arms(600);
    let mut p = picker(&t, Some(cells(600)), SiteFilterArgs::default_values());
    p.arm_look.key = true;
    for c in [
        Criterion::MinCells,
        Criterion::MaxPv,
        Criterion::MinCoverage,
    ] {
        focus_on(&mut p, c);
        for w in [1.55, 1.85, 2.2, 2.6, 3.0, 4.0, 6.0] {
            let sheet = sheet(w, Palette::Gray, true);
            p.arm_look.key_at = KeyAt::Right;
            let right = view_drawn(&p, sheet);
            assert!(has_key(&right), "{c:?} at {w} in");
            let (a, b) = plot_span(&right);
            assert!(
                b - a >= 0.6 * w * 72.0 - 0.5,
                "{c:?} at {w} in: the plot {} of {}",
                b - a,
                w * 72.0
            );
            p.arm_look.key_at = KeyAt::Top;
            let top = view_drawn(&p, sheet);
            match w >= 6.0 {
                true => assert!(plot_span(&right).1 < plot_span(&top).1 - 20.0, "{w}"),
                false if plot_span(&right) == plot_span(&top) => {
                    assert_eq!(right, top, "{c:?} at {w} in: right falls back to top")
                }
                false => {}
            }
        }
    }
}

/// Every name the key may list a series or its note by, long or short.
const KEY_NAMES: [&str; 9] = [
    "background (MUT)",
    "background (MUT, vs error)",
    "filtered out (WT)",
    "retained (WT)",
    "MUT",
    "WT filtered",
    "WT kept",
    "WT only",
    "no MUT arm",
];

/// The plot's top: where its highest upright line (the y axis, a
/// threshold's line) starts.
fn plot_top(svg: &str) -> f64 {
    (svg.lines())
        .filter(|l| l.starts_with("<line"))
        .filter_map(|l| {
            let attr = |k: &str| -> Option<f64> {
                l.split(&format!(" {k}=\""))
                    .nth(1)?
                    .split('"')
                    .next()?
                    .parse()
                    .ok()
            };
            let (x1, x2) = (attr("x1")?, attr("x2")?);
            (x1 == x2).then_some(attr("y1")?.min(attr("y2")?))
        })
        .fold(f64::INFINITY, f64::min)
}

/// The boxes of the key's names in `svg`.
fn key_boxes(svg: &str) -> Vec<(f64, f64, f64, f64)> {
    (crate::figure::text_boxes(svg).into_iter())
        .filter(|b| KEY_NAMES.contains(&b.1.as_str()))
        .map(|b| b.0)
        .collect()
}

/// The rectangles of `svg` short of the page's width: its bars and the
/// key's swatches.
fn rects(svg: &str, page_w: f64) -> Vec<(f64, f64, f64, f64)> {
    (svg.lines())
        .filter(|l| l.starts_with("<rect"))
        .filter_map(|l| {
            let attr = |k: &str| -> Option<f64> {
                l.split(&format!(" {k}=\""))
                    .nth(1)?
                    .split('"')
                    .next()?
                    .parse()
                    .ok()
            };
            Some((attr("x")?, attr("y")?, attr("width")?, attr("height")?))
        })
        .filter(|r| r.2 < page_w - 1.0)
        .collect()
}

/// Whether every one of the key's names and swatches in `svg` sits above
/// the plot: none over it.
fn key_off_the_plot(svg: &str, page_w: f64) -> bool {
    let top = plot_top(svg);
    let names = key_boxes(svg);
    !names.is_empty()
        && names.iter().all(|b| b.3 <= top + 0.01)
        && (rects(svg, page_w).iter()).all(|r| r.1 >= top - 0.01 || r.1 + r.3 <= top + 0.01)
}

/// The words either side of each threshold the focused knob draws.
fn both_words(p: &SitePicker) -> Vec<String> {
    (p.arm_hist.cuts.iter())
        .flat_map(|&(c, raw, _)| {
            let (lo, hi) = p.arm_hist.words(c, raw);
            [lo, hi]
        })
        .collect()
}

/// `--key-at top` on a narrow panel: the key in a row (or two) above the
/// plot, none of it over the plot; the plot as wide as with no key; the
/// words either side of each threshold drawn as with no key (on the
/// cells at 1.85 inches, both).
/// `--key-at right`, with no room beside the plot, draws it so too; and
/// `--key-at inside` still draws it over the plot, when asked.
#[test]
fn a_key_at_the_top_keeps_off_the_plot_and_its_words() {
    let t = with_arms(600);
    let mut p = picker(&t, Some(cells(600)), SiteFilterArgs::default_values());
    for c in [
        Criterion::MinCells,
        Criterion::MaxPv,
        Criterion::MinCoverage,
    ] {
        focus_on(&mut p, c);
        for w in [1.55, 1.85, 2.2] {
            let sheet = sheet(w, Palette::Gray, true);
            let page_w = w * 72.0;
            p.arm_look.key = false;
            let none = view_drawn(&p, sheet);
            p.arm_look.key = true;
            let said = |svg: &str| -> Vec<String> {
                (crate::figure::text_boxes(svg).into_iter())
                    .map(|b| b.1)
                    .collect()
            };
            for key_at in [KeyAt::Top, KeyAt::Right] {
                p.arm_look.key_at = key_at;
                let svg = view_drawn(&p, sheet);
                let at = format!("{c:?} {key_at:?} at {w} in");
                assert!(key_off_the_plot(&svg, page_w), "{at}: over the plot");
                assert_eq!(plot_span(&svg), plot_span(&none), "{at}: the plot narrowed");
                // Every word drawn with no key, and at 1.85 inches on
                // the cells both sides of its line.
                let (drawn, unkeyed) = (said(&svg), said(&none));
                let wanted = match (c, w) {
                    (Criterion::MinCells, 1.85) => both_words(&p),
                    _ => both_words(&p)
                        .into_iter()
                        .filter(|s| unkeyed.contains(s))
                        .collect(),
                };
                for word in wanted {
                    assert!(drawn.contains(&word), "{at}: {word:?} dropped");
                }
            }
            p.arm_look.key_at = KeyAt::Inside;
            let inside = view_drawn(&p, sheet);
            let (top, (a, b)) = (plot_top(&inside), plot_span(&inside));
            assert!(
                (key_boxes(&inside).iter()).any(|k| k.1 >= top && k.0 < b && k.2 > a),
                "{c:?} inside at {w} in: not over the plot"
            );
        }
    }
}

/// The p-value and the cells draw each arm: the signal arm on the p its
/// threshold reads (WT against MUT), the control arm on the p its own
/// reads give against sequencing error (the beta-binomial test the A-to-I
/// sites are called by); the cells counted on each arm's batches; and no
/// cut changes: the signal arm's kept are the sites every threshold
/// keeps, read off the same p-value as before.
#[test]
fn each_arm_s_p_and_cells_are_drawn_and_no_cut_changes() {
    use faba::hypothesis_tests::{betabinom_pvalue_greater, ERROR_RATE, OVERDISPERSION};
    let t = with_arms(700);
    let f = SiteFilterArgs::default_values();
    let n = cells(700);
    let mut p = picker(&t, Some(n.clone()), f.clone());
    let kept = kept_linear(&t, Some(&n), &f);
    let s = p.view().sites();
    for i in [0, 5, 97, 123, 699] {
        let want = |k: u64, n: u64| match betabinom_pvalue_greater(k, n, ERROR_RATE, OVERDISPERSION)
        {
            0.0 => f64::INFINITY,
            p => (-(p as f64).log10()).max(0.0),
        };
        assert_eq!(
            arms::wt(&s, Axis::Pvalue, i),
            Criterion::MaxPv.display(t.pv[i] as f64),
            "the p the cut reads"
        );
        assert_eq!(
            arms::mutant_hist(&s, Axis::Pvalue, i),
            want(t.control_converted[i], t.control_coverage[i])
        );
        assert_eq!(
            arms::mutant_hist(&s, Axis::Cells, i),
            t.control_cells.as_ref().unwrap()[i] as f64
        );
        assert_eq!(
            arms::wt(&s, Axis::Cells, i),
            n[i] as f64,
            "the cells the cut reads"
        );
    }
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    for c in [Criterion::MaxPv, Criterion::MinCells] {
        focus_on(&mut p, c);
        let h = &p.arm_hist;
        assert_eq!(sum(&h.kept), kept, "{c:?}: kept as the cut keeps");
        assert_eq!(sum(&h.kept) + sum(&h.filtered), t.len());
        assert_eq!(sum(h.mutant.as_ref().expect("each arm")), t.len());
    }
    // The 2D panels keep to the signal arm on the p-value and the cells.
    assert!(!arms::per_arm(Axis::Pvalue, &t) && !arms::per_arm(Axis::Cells, &t));
}

/// In gray the sparsest kept square reads clearly darker than the control
/// arm's light gray, its shade deeper as it holds more, on every count
/// scale, over a gradient that reads: from about half the ink up.
#[test]
fn gray_kept_squares_read_darker_than_the_background() {
    let t = with_arms(900);
    let v = SiteView::new(&t, Some(cells(900)), &SiteFilterArgs::permissive());
    let g = crate::qc::pairs::grid(
        &v.sites(),
        (Axis::EditRatio, Axis::Coverage),
        &SiteFilterArgs::permissive(),
        20,
    );
    let light = |rgb: (u8, u8, u8)| (rgb.0 as f64 + rgb.1 as f64 + rgb.2 as f64) / 3.0;
    let background = 0xc9 as f64;
    for scale in [CountScale::Linear, CountScale::Sqrt, CountScale::Log1p] {
        let mut shades: Vec<(u32, f64)> = Vec::new();
        for at in 0..g.squares.kept.len() {
            let (k, o) = (g.squares.kept[at], g.squares.other[at]);
            if let [Some(crate::qc::pairs::Tone::Fill(rgb)), _] =
                g.tones(at, Palette::Gray, (scale, None))
            {
                assert!(
                    background - light(rgb) >= 50.0,
                    "{scale:?}: {rgb:?} near the background"
                );
                shades.push((k + o, light(rgb)));
            }
        }
        shades.sort_by_key(|s| s.0);
        assert!(
            shades.windows(2).all(|w| w[0].1 >= w[1].1),
            "{scale:?}: darker as it holds more"
        );
        // The gradient on the default scale: the sparsest from not much
        // past half the ink, the busiest well past it.
        let ink = |l: f64| (255.0 - l) / (255.0 - light((0x26, 0x26, 0x26)));
        let (sparse, busy) = (ink(shades[0].1), ink(shades[shades.len() - 1].1));
        if scale == CountScale::Log1p {
            assert!(sparse <= 0.7, "{scale:?}: the sparsest {sparse} of the ink");
            assert!(busy - sparse >= 0.15, "{scale:?}: {sparse} to {busy}");
        }
    }
}

/// At a figure panel's width, a knob's axis takes at least four ticks
/// where they fit at the panel's type, each clear of the next by its
/// measured width, spread along the axis.
#[test]
fn a_narrow_panel_still_takes_four_ticks() {
    let n = 900;
    let mut t = with_arms(n);
    // Deep control sites, so the p-values reach far; cells in the
    // thousands on each arm.
    for i in (0..n).step_by(7) {
        let deep = 300 + i as u64;
        (t.control_coverage[i], t.control_converted[i]) = (deep, deep);
    }
    let wide = |k: usize| (0..n).map(|i| i * i * k % 6000).collect::<Vec<usize>>();
    t.control_cells = Some(wide(3));
    let mut p = picker(&t, Some(wide(7)), SiteFilterArgs::default_values());
    p.arm_look.key = true;
    for (c, w) in [
        (Criterion::MaxPv, 1.55),
        (Criterion::MinCells, 1.85),
        (Criterion::MinCoverage, 1.85),
        (Criterion::MinEditRatio, 2.2),
    ] {
        focus_on(&mut p, c);
        let sheet = sheet(w, Palette::Gray, true);
        let svg = view_drawn(&p, sheet);
        let choices: Vec<String> = (p.arm_hist.tick_choices(p.scale(), PAxis::P))
            .into_iter()
            .map(|t| t.1)
            .collect();
        let boxes: Vec<((f64, f64, f64, f64), String)> = crate::figure::text_boxes(&svg)
            .into_iter()
            .filter(|b| choices.contains(&b.1))
            .collect();
        // The tick row: the lowest of those.
        let row = boxes.iter().map(|b| b.0 .3).fold(0.0, f64::max);
        let mut ticks: Vec<(f64, String)> = (boxes.into_iter())
            .filter(|b| (b.0 .3 - row).abs() < 0.5)
            .map(|b| ((b.0 .0 + b.0 .2) / 2.0, b.1))
            .collect();
        ticks.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert!(ticks.len() >= 4, "{c:?} at {w} in: {ticks:?}");
        let size = Canvas::styled(1.0, 1.0, &sheet.style((1.0, 1.0))).label_size(8.0);
        for pair in ticks.windows(2) {
            let half = |s: &str| crate::figure::text_width(s, size) / 2.0;
            assert!(
                pair[1].0 - pair[0].0 >= half(&pair[0].1) + half(&pair[1].1),
                "{c:?} at {w} in: {pair:?}"
            );
        }
    }
}

/// The 2D panels' p-value axis ends at the histogram's cap: the round
/// cap of where the cut's p-values reach, its last squares holding every
/// p at or below it.
#[test]
fn the_2d_p_axis_ends_at_the_histogram_s_cap() {
    let mut t = with_arms(900);
    for i in (0..t.len()).step_by(5) {
        t.pv[i] = 10f32.powf(-(20.0 + (i % 17) as f32));
    }
    let f = SiteFilterArgs::default_values();
    let mut p = picker(&t, Some(cells(900)), f.clone());
    focus_on(&mut p, Criterion::MaxPv);
    let cap = p.arm_hist.clamped.0.expect("p-values past the cap");
    assert_eq!(cap, 30.0);
    let h = &p.arm_hist;
    assert_eq!(h.edges[h.edges.len() - 1], cap);
    let s = p.view().sites();
    for other in [Axis::Coverage, Axis::EditRatio] {
        let g = crate::qc::pairs::grid(&s, (Axis::Pvalue, other), &f, 20);
        let ruler = s.arm_ruler(Axis::Pvalue, &f);
        assert_eq!(g.xr.1, ruler.along(cap), "{other:?}");
        // Every site with both values is in a square.
        let held: u32 = g.squares.kept.iter().chain(&g.squares.other).sum();
        assert_eq!(held as usize, t.len(), "{other:?}");
    }
}

/// The p-value's panel names what each arm is drawn on, in its title and
/// in its key: the signal arm on the cut's own WT-vs-MUT test, the
/// background on the control arm's own test against error. Every signal
/// site right of the cut's line passes it and every one left fails it.
#[test]
fn the_p_panel_names_each_arm_s_test() {
    let t = with_arms(900);
    let f = SiteFilterArgs::default_values();
    let mut p = picker(&t, Some(cells(900)), f.clone());
    p.arm_look.key = true;
    // Over the corner, where the key keeps its long names.
    p.arm_look.key_at = KeyAt::Inside;
    focus_on(&mut p, Criterion::MaxPv);
    for (p_axis, title) in [
        (PAxis::P, "p: WT vs MUT (Fisher); MUT vs error"),
        (
            PAxis::Neglog10,
            "-log10 p: WT vs MUT (Fisher); MUT vs error",
        ),
    ] {
        p.arm_look.p_axis = p_axis;
        let svg = view_drawn(&p, sheet(3.0, Palette::Gray, true));
        let said: Vec<String> = crate::figure::labels(&svg)
            .into_iter()
            .map(|l| l.1)
            .collect();
        assert!(said.iter().any(|l| l == title), "{title} in {said:?}");
        assert!(said.iter().any(|l| l == "background (MUT, vs error)"));
    }
    // On a narrow panel the title breaks after its "; ", both arms named.
    p.arm_look.p_axis = PAxis::P;
    let narrow = view_drawn(&p, sheet(1.55, Palette::Gray, true));
    let said: Vec<String> = (crate::figure::labels(&narrow).into_iter())
        .map(|l| l.1)
        .collect();
    for line in ["p: WT vs MUT (Fisher);", "MUT vs error"] {
        assert!(said.iter().any(|l| l == line), "{line} in {said:?}");
    }
    // Another knob's key names its background as it was.
    focus_on(&mut p, Criterion::MinCoverage);
    let svg = view_drawn(&p, sheet(3.0, Palette::Gray, true));
    assert!(has_key(&svg));
    assert!(crate::figure::labels(&svg)
        .iter()
        .any(|l| l.1 == "background (MUT)"));
    // The line splits the signal arm's sites as the cut does.
    focus_on(&mut p, Criterion::MaxPv);
    let h = &p.arm_hist;
    let e = h.cuts[0].2;
    let (left, right): (usize, usize) = (h.kept[..e].iter().sum(), h.filtered[e..].iter().sum());
    assert_eq!(left, 0, "kept left of the line");
    let passing_but_filtered = (0..t.len())
        .filter(|&i| t.pv[i] <= f.site_max_pv && p.view().fails[i] != 0)
        .count();
    assert_eq!(right, passing_but_filtered);
}

/// Each 2D marginal's count is ticked at the panel's own type, on its own
/// count axis right next to it: the top one's beside it at the left, the
/// right one's over it, in the corner by the top one; each at the round
/// count its axis ends at, clear of every other label. The
/// marginals are on log1p unless asked, so one busy bar leaves the others
/// a depth that reads.
#[test]
fn a_marginal_s_count_is_ticked_beside_it_at_the_panel_s_type() {
    use crate::qc::pairs::{framed, grid, PairLook};
    assert_eq!(arms::ArmLook::default().marginal, CountScale::Log1p);
    let t = with_arms(900);
    let f = SiteFilterArgs::default_values();
    let v = SiteView::new(&t, Some(cells(900)), &f);
    let s = v.sites();
    let marginals = data_beans::figure::panel::Marginals {
        top: true,
        right: true,
    };
    for (x, y) in [
        (Axis::Pvalue, Axis::Coverage),
        (Axis::EditRatio, Axis::Coverage),
    ] {
        let g = grid(&s, (x, y), &f, 20);
        let sheet = sheet(2.6, Palette::Gray, true);
        let (w, h) = sheet.fill(|| unreachable!());
        let mut c = Canvas::styled(w, h, &sheet.style((w, h)));
        let (pw, ph) = c.page();
        let how = PairLook {
            arms: arms::ArmLook::default(),
            titles: None,
        };
        g.draw_with(&mut c, (0.0, 0.0, pw, ph), "", (marginals, 0.18), &how);
        let svg = c.finish();
        let lay = framed((0.0, 0.0, pw, ph), g.xr, g.yr, marginals, 0.18);
        let plot = &lay.plot;
        let font = |said: &str| -> String {
            let line = (svg.lines())
                .find(|l| l.starts_with("<text") && l.contains(&format!(">{said}</text>")))
                .unwrap_or_else(|| panic!("{said}"));
            line.split(" font-size=\"")
                .nth(1)
                .unwrap()
                .split('"')
                .next()
                .unwrap()
                .into()
        };
        let boxes = crate::figure::text_boxes(&svg);
        // An x tick's type: a label under the plot, inside it across.
        let x_tick = (boxes.iter())
            .find(|(b, _)| b.1 > plot.py + plot.ph && b.0 >= plot.px && b.2 <= plot.px + plot.pw)
            .expect("an x tick");
        let tick_font = font(&x_tick.1);
        let (cols, rows) = g.margins_by(crate::qc::arms::KeptBy::All);
        let (mcols, mrows) = g.mutant_margins();
        let most = |wt: &[(f64, f64)], m: &[f64]| {
            (wt.iter().map(|(k, o)| k + o))
                .chain(m.iter().copied())
                .fold(1.0, f64::max)
        };
        let scale = CountScale::Log1p;
        let end = |most: f64| crate::tui::ui::compact(*scale.ticks(most).last().unwrap());
        let (ty, th) = (lay.top.unwrap().1, lay.top.unwrap().3);
        let (rx, rw) = (lay.right.unwrap().0, lay.right.unwrap().2);
        let top_said = end(most(&cols, &mcols));
        let beside_top = (boxes.iter()).any(|(b, t)| {
            *t == top_said
                && b.2 <= plot.px
                && (b.1 + b.3) / 2.0 >= ty - 1.0
                && (b.1 + b.3) / 2.0 <= ty + th + 1.0
        });
        assert!(beside_top, "{x:?}: {top_said} beside the top marginal");
        let right_said = end(most(&rows, &mrows));
        let over_right = (boxes.iter()).any(|(b, t)| {
            *t == right_said
                && b.3 <= plot.py
                && b.1 >= ty + th - 2.0 * (b.3 - b.1)
                && (b.0 + b.2) / 2.0 >= rx - 1.0
                && (b.0 + b.2) / 2.0 <= rx + rw + 1.0
        });
        assert!(over_right, "{x:?}: {right_said} over the right marginal");
        // Clear of every other label.
        for said in [&top_said, &right_said] {
            let mine: Vec<_> = boxes.iter().filter(|(_, t)| t == said).collect();
            for (b, t) in boxes.iter().filter(|(_, t)| t != said) {
                for (m, _) in &mine {
                    let apart = b.2 <= m.0 || m.2 <= b.0 || b.3 <= m.1 || m.3 <= b.1;
                    assert!(apart, "{x:?}: {said} over {t}");
                }
            }
        }
        for said in [&top_said, &right_said] {
            assert_eq!(font(said), tick_font, "{x:?}: {said} at the panel's type");
        }
    }
}
