//! The view style's two arms, in the one drawing the view on screen and
//! every panel share: the series each bar holds add up to the site table,
//! each threshold's line sits on a cut between bars, the ticks name what
//! is at their place, the 2D squares and their marginals count the same
//! sites and arms, and `--front` changes the order drawn, nothing else.

use super::*;
use crate::figure::style::{Palette, StyleArgs};
use crate::qc::arms::{self, ArmLook, Front, PAxis};
use crate::qc::pairs::{self as qp, Axis};
use crate::qc::site_plot::{Cut, PanelArgs, PanelStyle, Panels};

/// A table with a control arm and cells, at thresholds that cut every
/// knob somewhere inside a bin: none sits on a bin's edge.
fn mid_bin() -> SiteFilterArgs {
    SiteFilterArgs {
        site_max_pv: 0.013,
        site_min_log_odds: 0.37,
        site_min_fold: 1.7,
        site_min_coverage: 23,
        site_min_converted: 7,
        site_min_edit_ratio: 0.083,
        site_max_edit_ratio: 0.91,
        site_min_cells: 10,
    }
}

/// The picker over a table of `n` sites with cells, at `f`.
fn arms_picker(t: &SiteTable, f: SiteFilterArgs) -> SitePicker<'_> {
    picker(t, Some(cells(t.len())), f)
}

/// Site `i`'s value of knob `c` on the signal arm and on the control arm
/// (NaN where it has none), as its histogram shows it.
fn values(p: &SitePicker, c: Criterion, i: usize) -> (f64, f64) {
    let s = p.view().sites();
    let q = Axis::of(c);
    (arms::wt(&s, q, i), arms::mutant_hist(&s, q, i))
}

/// The bar of `h` holding `v`, found by walking its edges (the end bars
/// take what lies past them, an infinity at the end it is drawn at).
fn walked(h: &ArmHist, v: f64) -> usize {
    let n = h.bars();
    (0..n)
        .find(|&j| h.edges[j] <= v && v < h.edges[j + 1])
        .unwrap_or(if v < h.edges[0] { 0 } else { n - 1 })
}

/// Every knob's bars hold every site of the table once on each arm it has
/// a value on: the signal arm's kept (by every threshold) and filtered
/// (by any) add up to its sites with a value, each in the bar whose edges
/// hold it; the control arm's every site's control value where the
/// quantity has one (the reads, the edit ratio), else none at all.
#[test]
fn the_series_add_up_to_the_site_table() {
    let t = table(M6A, 1200);
    for f in [SiteFilterArgs::default_values(), mid_bin()] {
        let mut p = arms_picker(&t, f.clone());
        for c in p.view().criteria.clone() {
            focus_on(&mut p, c);
            let h = &p.arm_hist;
            let n = h.bars();
            let (mut kept, mut filtered, mut mutant) = (vec![0; n], vec![0; n], vec![0; n]);
            let (lo, hi) = (h.edges[0], h.edges[n]);
            let place = |v: f64| match v {
                f64::INFINITY => Some(walked(h, h.clamped.0.unwrap())),
                f64::NEG_INFINITY => Some(walked(h, h.clamped.1.unwrap())),
                v if v.is_nan() => None,
                v => Some(walked(h, v)),
            };
            for i in 0..t.len() {
                let (wt, m) = values(&p, c, i);
                // Past a p-value's cap, in its last bar.
                let capped = h.clamped.0.is_some_and(|cap| wt > cap);
                assert!(
                    wt.is_nan() || wt.is_infinite() || capped || (lo..hi).contains(&wt),
                    "{c:?}"
                );
                if let Some(b) = place(wt) {
                    match p.view().fails[i] == 0 {
                        true => kept[b] += 1,
                        false => filtered[b] += 1,
                    }
                }
                if let Some(b) = place(m) {
                    mutant[b] += 1;
                }
            }
            assert_eq!(h.kept, kept, "{c:?}");
            assert_eq!(h.filtered, filtered, "{c:?}");
            let valued = (0..t.len())
                .filter(|&i| !values(&p, c, i).0.is_nan())
                .count();
            assert_eq!(
                kept.iter().sum::<usize>() + filtered.iter().sum::<usize>(),
                valued
            );
            let all_kept = p.view().fails.iter().filter(|&&m| m == 0).count();
            assert_eq!(kept.iter().sum::<usize>(), all_kept, "{c:?}");
            match arms::per_arm_hist(Axis::of(c), &t) {
                true => {
                    assert_eq!(h.mutant.as_ref(), Some(&mutant), "{c:?}");
                    let with = (0..t.len())
                        .filter(|&i| !values(&p, c, i).1.is_nan())
                        .count();
                    assert_eq!(mutant.iter().sum::<usize>(), with, "{c:?}");
                }
                false => assert!(h.mutant.is_none(), "{c:?}: no control arm to draw"),
            }
        }
    }
    // Every site's control reads are drawn, and its edit ratio wherever
    // its control arm has a read.
    let p = arms_picker(&t, mid_bin());
    let s = p.view().sites();
    for i in 0..t.len() {
        assert_eq!(
            arms::mutant(&s, Axis::Coverage, i),
            t.control_coverage[i] as f64
        );
        assert_eq!(
            arms::mutant(&s, Axis::Converted, i),
            t.control_converted[i] as f64
        );
        let ratio = t.control_converted[i] as f32 / t.control_coverage[i] as f32;
        assert_eq!(arms::mutant(&s, Axis::EditRatio, i), ratio as f64);
        assert_eq!(arms::wt(&s, Axis::Coverage, i), t.coverage[i] as f64);
    }
    // Without a control arm nothing is drawn for it.
    let atoi = table(ATOI, 300);
    let mut p = arms_picker(&atoi, mid_bin());
    focus_on(&mut p, Criterion::MinCoverage);
    assert!(p.arm_hist.mutant.is_none());
}

/// The bars are the knob's own bins, a threshold's bin split at it: each
/// edge not a threshold's is where one of the knob's bins starts.
#[test]
fn the_bars_are_the_knobs_bins() {
    let t = table(M6A, 1200);
    let mut p = arms_picker(&t, mid_bin());
    for c in p.view().criteria.clone() {
        focus_on(&mut p, c);
        let (h, bins) = (&p.arm_hist, p.column.hist.bins);
        let thresholds: Vec<f64> = (h.cuts.iter()).map(|cut| h.edges[cut.2]).collect();
        for &e in &h.edges[1..h.bars()] {
            if thresholds.contains(&e) {
                continue;
            }
            let below = if c.is_integer() {
                e - 1.0
            } else {
                e.next_down()
            };
            assert!(bins.key(e) > bins.key(below), "{c:?}: {e} starts no bin");
        }
    }
}

/// Each threshold's line sits on a cut between bars, and that cut is the
/// threshold's: every signal site in a bar on its failing side fails it,
/// every one on the other side passes, at thresholds inside a bin and on
/// a bin's edge alike (the coverage's adds both arms' reads, so on one
/// arm's axis its line marks its value and its words say so). No bar of
/// the drawing crosses the line.
#[test]
fn each_threshold_line_sits_at_the_cut_between_bars() {
    let t = table(M6A, 1200);
    // Thresholds at a site's own value, where the bound's side of it
    // matters: kept by an upper bound, by a lower one too.
    let at_values = SiteFilterArgs {
        site_max_pv: t.pv[5],
        site_min_edit_ratio: t.edit_ratio[7],
        site_max_edit_ratio: t.edit_ratio[11],
        site_min_converted: t.converted[13],
        ..mid_bin()
    };
    for f in [SiteFilterArgs::default_values(), mid_bin(), at_values] {
        let mut p = arms_picker(&t, f.clone());
        for c in p.view().criteria.clone() {
            focus_on(&mut p, c);
            let h = &p.arm_hist;
            let Some(&(_, raw, e)) = h.cuts.iter().find(|cut| cut.0 == c) else {
                assert!(c.is_off(&f), "{c:?}: on, but no line");
                continue;
            };
            assert_eq!(raw, c.get(&f));
            assert!(e > 0 && e < h.edges.len(), "{c:?}: the line at an end");
            if c == Criterion::MinCoverage {
                assert_eq!(h.edges[e], raw, "{c:?}");
                assert!(h.words(c, raw).0.contains("WT+MUT"));
                continue;
            }
            for i in 0..t.len() {
                let wt = values(&p, c, i).0;
                if wt.is_nan() {
                    continue;
                }
                let v = match wt {
                    f64::INFINITY => h.clamped.0.unwrap(),
                    f64::NEG_INFINITY => h.clamped.1.unwrap(),
                    v => v,
                };
                let left = walked(h, v) < e;
                let fails = c.fails(&f, &t, i, p.view().n_cells.as_ref().map(|n| n[i]));
                // The left side fails a lower bound and passes an upper one.
                let failing_side = if c.is_max() && !c.neg_log() {
                    !left
                } else {
                    left
                };
                assert_eq!(
                    fails, failing_side,
                    "{c:?} site {i}: {v} cut at {}",
                    h.edges[e]
                );
            }
            let svg = drawn_hist(&p, Sheet::OWN);
            let lines = threshold_xs(&svg);
            assert_eq!(lines.len(), h.cuts.len(), "{c:?}");
            for line in lines {
                // Each written to a hundredth: a bar ending at the line may
                // read a hundredth past it.
                for (x, w) in bar_rects(&svg) {
                    let clear = x + w <= line + 0.015 || x >= line - 0.015;
                    assert!(clear, "{c:?}: a bar at {x}+{w} crosses {line}");
                }
            }
        }
    }
}

/// The focused knob's histogram as the view draws it on `sheet`'s page.
fn drawn_hist(p: &SitePicker, sheet: Sheet) -> String {
    let (w, h) = sheet.fill(|| (390.0, 280.0));
    let mut c = Canvas::styled(w, h, &sheet.style((w, h)));
    let (pw, ph) = c.page();
    p.draw_hist(&mut c, (0.0, 0.0, pw, ph), "TITLE1".into());
    c.finish()
}

/// Where each threshold's dashed line is, across: the upright ones (the
/// key's swatch lies flat).
fn threshold_xs(svg: &str) -> Vec<f64> {
    let stroke = format!("stroke=\"{}\"", arms::THRESHOLD);
    let attr = |l: &str, k: &str| -> f64 {
        let v = l
            .split(&format!(" {k}=\""))
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap();
        v.parse().unwrap()
    };
    (svg.lines())
        .filter(|l| l.starts_with("<line") && l.contains(&stroke))
        .filter(|l| attr(l, "x1") == attr(l, "x2"))
        .map(|l| attr(l, "x1"))
        .collect()
}

/// Each bar's left edge and width: the rectangles over the plot, not the
/// page nor the key's swatches (7 across).
fn bar_rects(svg: &str) -> Vec<(f64, f64)> {
    let attr = |l: &str, k: &str| -> f64 {
        let v = l
            .split(&format!(" {k}=\""))
            .nth(1)
            .unwrap()
            .split('"')
            .next()
            .unwrap();
        v.parse().unwrap()
    };
    (svg.lines())
        .filter(|l| l.starts_with("<rect"))
        .map(|l| (attr(l, "x"), attr(l, "width")))
        .filter(|&(x, w)| x > 0.0 && (w - 7.0).abs() > 1e-9)
        .collect()
}

/// The words either side of a line name the threshold and its value, the
/// low end's first.
#[test]
fn the_words_name_both_sides() {
    let t = table(M6A, 600);
    let mut p = arms_picker(&t, mid_bin());
    for (c, low, high) in [
        (Criterion::MaxPv, "p > 0.013", "p <= 0.013"),
        (Criterion::MinCells, "< 10 cells", ">= 10 cells"),
        (Criterion::MinConverted, "< 7 reads", ">= 7 reads"),
        (Criterion::MaxEditRatio, "<= 0.91", "> 0.91"),
        (Criterion::MinCoverage, "WT+MUT < 23", "WT+MUT >= 23"),
    ] {
        focus_on(&mut p, c);
        assert_eq!(
            p.arm_hist.words(c, c.get(&p.filter)),
            (low.into(), high.into())
        );
        let labels: Vec<String> = (crate::figure::text_boxes(&drawn_hist(&p, Sheet::OWN)))
            .into_iter()
            .map(|l| l.1)
            .collect();
        // Each side's words where there is room for them on that side.
        let h = &p.arm_hist;
        let e = (h.cuts.iter()).find(|cut| cut.0 == c).unwrap().2;
        let (lo, hi) = (h.at[0], h.at[h.at.len() - 1]);
        let room = |to: f64| (to - h.at[e]).abs() / (hi - lo);
        for (said, side) in [(low, lo), (high, hi)] {
            if room(side) > 0.25 {
                assert!(
                    labels.iter().any(|l| l == said),
                    "{c:?}: {said} in {labels:?}"
                );
            }
        }
        assert!(!labels.iter().any(|l| l.contains('←') || l.contains('→')));
    }
}

/// The p-value's ticks name the p at their place along `-log10 p` (or
/// `-log10 p` itself), the bars where they were: where the p-values of 0
/// are drawn, "<= p" of the smallest other.
#[test]
fn p_ticks_name_the_p_at_their_place() {
    use crate::qc::scale::Scales;
    let mut t = table(M6A, 1200);
    // A control site whose own reads give a p of 0: drawn at the end.
    (t.control_coverage[0], t.control_converted[0]) = (5000, 5000);
    let mut p = arms_picker(&t, mid_bin());
    // On -log10 p as it is, as asked: a tick's place is its -log10 p.
    p.set_scales(Scales::parse(&["pv=linear"]).unwrap(), 50);
    focus_on(&mut p, Criterion::MaxPv);
    let h = &p.arm_hist;
    let hi = h.clamped.0.expect("p-values of 0 are drawn at the end");
    let scale = p.scale();
    assert_eq!(scale, Scale::Linear, "a tick's place is its -log10 p");
    let number = |s: &str| {
        s.trim_start_matches("<= ")
            .trim_start_matches(">= ")
            .parse::<f64>()
            .unwrap()
    };
    let p_ticks = h.tick_choices(scale, PAxis::P);
    let neg_ticks = h.tick_choices(scale, PAxis::Neglog10);
    assert_eq!(p_ticks.len(), neg_ticks.len());
    for ((at, said, end), (at2, neg, _)) in p_ticks.iter().zip(&neg_ticks) {
        assert_eq!(at, at2, "the ticks stay where they were");
        // The bins are linear on -log10 p: a tick's place is its value.
        let d = *at;
        let pv = 10f64.powf(-d);
        let read = number(said);
        assert!(
            (read - pv).abs() <= 0.06 * pv,
            "{said} at -log10 p {d} (p {pv})"
        );
        assert!(
            (number(neg) - d).abs() <= 0.05 * d.abs() + 0.05,
            "{neg} at {d}"
        );
        if *end {
            assert_eq!(d, hi);
            assert!(
                said.starts_with("<= ") && neg.starts_with(">= "),
                "{said} {neg}"
            );
        }
    }
    // The ticks drawn are among them, as each mode says.
    for (look, ticks) in [(PAxis::P, &p_ticks), (PAxis::Neglog10, &neg_ticks)] {
        p.arm_look.p_axis = look;
        let svg = drawn_hist(&p, Sheet::OWN);
        let labels: Vec<String> = crate::figure::labels(&svg)
            .into_iter()
            .map(|l| l.1)
            .collect();
        let drawn: Vec<&String> = ticks
            .iter()
            .map(|t| &t.1)
            .filter(|t| labels.contains(t))
            .collect();
        assert!(drawn.len() >= 3, "{look:?}: {drawn:?}");
    }
}

/// The cells on log bins are labelled at their bins' edges, the least
/// whole count each holds: a tick says the count its bar starts at.
#[test]
fn cell_ticks_are_the_true_bin_edges() {
    let t = table(M6A, 1200);
    let mut p = arms_picker(&t, mid_bin());
    focus_on(&mut p, Criterion::MinCells);
    assert_eq!(p.scale(), Scale::Log);
    let (h, bins) = (&p.arm_hist, p.column.hist.bins);
    let ticks = h.tick_choices(p.scale(), PAxis::P);
    assert_eq!(ticks[0].1, "0");
    for (at, said, _) in &ticks {
        let e: f64 = said.parse().unwrap();
        assert_eq!(e.fract(), 0.0, "{said}");
        let j = h
            .edges
            .iter()
            .position(|&x| x == e)
            .expect("a tick is an edge");
        assert_eq!(*at, h.at[j]);
        // The least whole count of its bin: the count before is in another.
        if e > 0.0 {
            assert!(bins.key(e) > bins.key(e - 1.0), "{said} starts no bin");
        }
        // And a bar starts there that holds it (past the last, none).
        if j < h.bars() {
            assert_eq!(walked(h, e), j);
        }
    }
}

/// The 2D squares and their marginals count the same sites and arms: the
/// signal arm's every site with both values, the control arm's every
/// site with both on it, where both axes have a value on each arm; each
/// marginal adds its squares up.
#[test]
fn the_squares_and_their_marginals_count_the_same_sites_and_arms() {
    let t = table(M6A, 1500);
    let f = mid_bin();
    let v = SiteView::new(&t, Some(cells(1500)), &f);
    let s = v.sites();
    let both = |a: Axis, b: Axis, value: fn(&crate::qc::pairs::Sites, Axis, usize) -> f64| {
        (0..t.len())
            .filter(|&i| !value(&s, a, i).is_nan() && !value(&s, b, i).is_nan())
            .count()
    };
    let g = qp::grid(&s, (Axis::Coverage, Axis::EditRatio), &f, qp::BINS);
    let wt: u32 = g.squares.kept.iter().chain(&g.squares.other).sum();
    assert_eq!(wt as usize, both(Axis::Coverage, Axis::EditRatio, arms::wt));
    let mutant = g
        .mutant
        .as_ref()
        .expect("both axes have a value on each arm");
    let m: u32 = mutant.iter().sum();
    assert_eq!(
        m as usize,
        both(Axis::Coverage, Axis::EditRatio, arms::mutant)
    );
    let (cols, rows) = g.margins_by(crate::qc::arms::KeptBy::All);
    let (mcols, mrows) = g.mutant_margins();
    for (wt_side, m_side) in [(&cols, &mcols), (&rows, &mrows)] {
        assert_eq!(wt_side.iter().map(|p| p.0 + p.1).sum::<f64>(), wt as f64);
        assert_eq!(m_side.iter().sum::<f64>(), m as f64);
    }
    // The histogram of the coverage counts the control arm's sites alike.
    let mut p = picker(&t, Some(cells(1500)), f.clone());
    focus_on(&mut p, Criterion::MinCoverage);
    let hist_mut: usize = p.arm_hist.mutant.as_ref().unwrap().iter().sum();
    assert_eq!(hist_mut, t.len());
    // A pair with an axis that compares the arms draws the signal alone.
    let g = qp::grid(&s, (Axis::Pvalue, Axis::Cells), &f, qp::BINS);
    assert!(g.mutant.is_none());
    assert_eq!(
        g.mutant_margins(),
        (vec![0.0; g.squares.bins], vec![0.0; g.squares.bins])
    );
    // Each control site sits in the square of its own two values, its
    // coverage the site's WT + MUT reads, as the coverage's cut reads it.
    let g = qp::grid(&s, (Axis::Coverage, Axis::Converted), &f, qp::BINS);
    let (rx, ry) = (
        s.arm_ruler(Axis::Coverage, &f),
        s.arm_ruler(Axis::Converted, &f),
    );
    // As the squares are counted, in single precision along the line.
    let square = |r: &crate::qc::scale::Ruler, v: f64| {
        let (lo, hi) = r.span();
        let k = (g.squares.bins as f64 / (hi - lo)) as f32;
        let last = g.squares.bins as f32 - 1.0;
        ((r.along(v) as f32 - lo as f32) * k)
            .floor()
            .clamp(0.0, last) as usize
    };
    let mut want = vec![0u32; g.squares.bins * g.squares.bins];
    for i in 0..t.len() {
        let x = (t.coverage[i] + t.control_coverage[i]) as f64;
        let y = t.control_converted[i] as f64;
        want[square(&ry, y) * g.squares.bins + square(&rx, x)] += 1;
    }
    assert_eq!(g.mutant.as_ref(), Some(&want));
    // And each marginal bar adds up its own column's, or row's, squares.
    let n = g.squares.bins;
    let (mcols, mrows) = g.mutant_margins();
    for k in 0..n {
        let col: u32 = (0..n).map(|j| want[j * n + k]).sum();
        let row: u32 = (0..n).map(|j| want[k * n + j]).sum();
        assert_eq!((mcols[k], mrows[k]), (col as f64, row as f64), "{k}");
    }
}

/// Every `<rect>` of a drawing, as written, in order, but the key's
/// swatches (7 by 7), which list the arm in front first.
fn rects(svg: &str) -> Vec<&str> {
    (svg.lines())
        .filter(|l| l.starts_with("<rect") && !l.contains(r#"width="7.00" height="7.00""#))
        .collect()
}

/// `--front` (`f`) changes which arm is drawn over the other, nothing
/// else: the same bars and squares, in another order, the key listing the
/// arm in front first.
#[test]
fn front_swaps_the_draw_order_only() {
    let t = table(M6A, 800);
    let mut p = arms_picker(&t, mid_bin());
    focus_on(&mut p, Criterion::MinCoverage);
    let hist = |p: &SitePicker| drawn_hist(p, Sheet::OWN);
    let wt = hist(&p);
    press(&mut p, KeyCode::Char('f'));
    assert_eq!(p.arm_look.front, Front::Mut);
    let mutant = hist(&p);
    let sorted = |svg: &str| {
        let mut r = rects(svg);
        r.sort_unstable();
        r.into_iter().map(String::from).collect::<Vec<_>>()
    };
    assert_eq!(sorted(&wt), sorted(&mutant), "the same bars");
    assert_ne!(rects(&wt), rects(&mutant), "in another order");
    // The control arm's bars come last with it in front, first without.
    let first_mut = |svg: &str| {
        let all = rects(svg);
        let mutant = format!("fill=\"{}\"", arms::MUT_FILL);
        let at = |is: &dyn Fn(&str) -> bool| all.iter().position(|r| is(r));
        (
            at(&|r| r.contains(&mutant)),
            at(&|r| r.contains(" stroke=")),
        )
    };
    let ((m1, f1), (m2, f2)) = (first_mut(&wt), first_mut(&mutant));
    assert!(m1 < f1 && m2 > f2, "{m1:?} {f1:?} / {m2:?} {f2:?}");
    // The squares and their marginals alike.
    let s = p.view().sites();
    let g = qp::grid(&s, (Axis::Coverage, Axis::EditRatio), &p.filter, qp::BINS);
    let pair = |front: Front| {
        let how = qp::PairLook {
            arms: ArmLook {
                front,
                ..Default::default()
            },
            titles: None,
        };
        crate::figure::svg(300.0, 300.0, |c| {
            g.draw_with(
                c,
                (0.0, 0.0, 300.0, 300.0),
                "",
                (both_marginals(), 0.15),
                &how,
            )
        })
    };
    let (a, b) = (pair(Front::Wt), pair(Front::Mut));
    assert_eq!(sorted(&a), sorted(&b));
    assert_ne!(rects(&a), rects(&b));
}

/// A box, `(x, y, w, h)`.
type Box4 = (f64, f64, f64, f64);

/// Marginals on both sides.
fn both_marginals() -> data_beans::figure::panel::Marginals {
    data_beans::figure::panel::Marginals {
        top: true,
        right: true,
    }
}

/// The view on screen and the panel draw the same series whichever arm is
/// in front and however the p-value reads: `s` in `faba qc-view` saves
/// what the view draws, and `faba site-plot --style view` draws it too.
#[test]
fn the_view_and_the_panel_show_the_same_series() {
    let t = table(M6A, 600);
    let n = cells(600);
    let f = mid_bin();
    let sheet = Sheet {
        page: (Some(2.4 * 72.0), Some(1.8 * 72.0)),
        look: StyleArgs {
            font_size: Some(7.0),
            panel_only: true,
            palette: Palette::Gray,
        },
    };
    for (front, p_axis) in [(Front::Wt, PAxis::P), (Front::Mut, PAxis::Neglog10)] {
        let args = PanelArgs {
            hist: vec!["pv".into(), "coverage".into()],
            pairs: vec!["coverage:edit_ratio".into()],
            style: PanelStyle::View,
            front: Some(front),
            p_axis: Some(p_axis),
            key: true,
            ..PanelArgs::default()
        };
        let plan = args.plan().unwrap();
        let mut p = picker(&t, Some(n.clone()), f.clone());
        p.view_only(plan.clone(), sheet);
        p.set_scales(plan.scales, plan.bins_on(&sheet));
        let panels = Panels::new(&t, Some(&n), &Cut::of_qc(&f), &plan, sheet, true)
            .unwrap()
            .figures();
        // What the view draws on screen, on the panel's page.
        let (w, h) = sheet.fill(|| unreachable!("both sides asked"));
        let on_screen = |draw: &dyn Fn(&mut Canvas, Box4)| {
            let mut c = Canvas::styled(w, h, &sheet.style((w, h)));
            let (pw, ph) = c.page();
            draw(&mut c, (0.0, 0.0, pw, ph));
            c.finish()
        };
        for c in [Criterion::MaxPv, Criterion::MinCoverage] {
            focus_on(&mut p, c);
            let saved = p.screen_panel(&plan, sheet).unwrap();
            let stem = format!("m6a_hist_{}", Axis::of(c).key());
            let panel = panels.iter().find(|g| g.stem == stem).unwrap();
            assert_eq!(saved, panel.svg, "{stem} {front:?} {p_axis:?}");
            let shown = on_screen(&|canvas, b| p.draw_hist(canvas, b, String::new()));
            assert_eq!(shown, panel.svg, "{stem} on screen {front:?} {p_axis:?}");
        }
        press(&mut p, KeyCode::Char('v'));
        p.views[0].pairs = vec![(Axis::Coverage, Axis::EditRatio, "chosen by hand".into())];
        let saved = p.screen_panel(&plan, sheet).unwrap();
        let panel = panels
            .iter()
            .find(|g| g.stem == "m6a_edit_ratio_vs_coverage")
            .unwrap();
        assert_eq!(saved, panel.svg, "pair {front:?} {p_axis:?}");
        let s = p.view().sites();
        let g = qp::grid(&s, (Axis::Coverage, Axis::EditRatio), &p.filter, qp::BINS);
        let how = qp::PairLook {
            arms: p.arm_look,
            titles: Some(&p.titles),
        };
        let marginals = (plan.marginals, plan.marginal_size);
        let shown = on_screen(&|c, b| g.draw_with(c, b, "", marginals, &how));
        assert_eq!(shown, panel.svg, "pair on screen {front:?} {p_axis:?}");
    }
    // On screen without images, the key counts the series the bars draw.
    let mut p = picker(&t, Some(n.clone()), f.clone());
    focus_on(&mut p, Criterion::MinCoverage);
    let text = screen_text(&mut p, 160, 44);
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    let h = &p.arm_hist;
    for said in [
        format!("retained (WT) {}", sum(&h.kept)),
        format!("filtered out (WT) {}", sum(&h.filtered)),
        format!("background (MUT) {}", sum(h.mutant.as_ref().unwrap())),
    ] {
        assert!(text.contains(&said), "{said}");
    }
}

/// The screen as text.
fn screen_text(p: &mut SitePicker, w: u16, h: u16) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    term.backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect()
}

/// The key goes with `--panel-only` unless `--key` keeps it; `--p-axis`,
/// `--front` and `--key` are the view style's, refused with the panel's.
#[test]
fn the_key_stays_alone_only_when_asked() {
    let t = table(M6A, 400);
    let mut p = arms_picker(&t, mid_bin());
    focus_on(&mut p, Criterion::MinCoverage);
    let alone = Sheet {
        look: StyleArgs {
            panel_only: true,
            ..Default::default()
        },
        ..Sheet::OWN
    };
    let has_key = |svg: &str| {
        crate::figure::labels(svg)
            .iter()
            .any(|l| l.1 == "retained (WT)")
    };
    assert!(has_key(&drawn_hist(&p, Sheet::OWN)));
    assert!(!has_key(&drawn_hist(&p, alone)));
    p.arm_look.key = true;
    assert!(has_key(&drawn_hist(&p, alone)));
    for args in [
        PanelArgs {
            key: true,
            ..PanelArgs::default()
        },
        PanelArgs {
            front: Some(Front::Mut),
            ..PanelArgs::default()
        },
        PanelArgs {
            p_axis: Some(PAxis::Neglog10),
            ..PanelArgs::default()
        },
    ] {
        let e = args.plan().unwrap_err().to_string();
        assert!(e.contains("--style view"), "{e}");
    }
    // A knob of a quantity that compares the arms says it draws WT alone.
    focus_on(&mut p, Criterion::MinFold);
    let labels = crate::figure::labels(&drawn_hist(&p, Sheet::OWN));
    assert!(labels.iter().any(|l| l.1 == "WT only"));
    assert!(!labels.iter().any(|l| l.1 == "background (MUT)"));
}

/// `f` puts the other arm in front and `u` reads the p-value the other
/// way, each said in the header; `t` types the x and then the y title of
/// the panel `s` saves, which it saves with them; Esc leaves them.
#[test]
fn f_u_and_t_restyle_the_panel_s_saves() {
    let t = table(M6A, 400);
    let panels = PanelArgs {
        style: PanelStyle::View,
        ..PanelArgs::default()
    };
    let sheet = Sheet {
        page: (Some(200.0), Some(150.0)),
        look: StyleArgs::default(),
    };
    let mut p = picker(&t, Some(cells(400)), SiteFilterArgs::default_values());
    p.view_only(panels.plan().unwrap(), sheet);
    let header = |p: &mut SitePicker| screen_text(p, 160, 44);
    assert!(header(&mut p).contains("WT in front · p axis p-value"));
    press(&mut p, KeyCode::Char('f'));
    press(&mut p, KeyCode::Char('u'));
    assert!(header(&mut p).contains("MUT in front · p axis -log10 p"));
    let (plan, sheet) = p.viewing.clone().unwrap();
    focus_on(&mut p, Criterion::MinCoverage);
    let typed = |p: &mut SitePicker, text: &str| {
        text.chars().for_each(|ch| press(p, KeyCode::Char(ch)));
        press(p, KeyCode::Enter);
    };
    press(&mut p, KeyCode::Char('t'));
    assert!(p.footer_line().to_string().contains("x title of coverage"));
    typed(&mut p, "TITLE1");
    assert!(p
        .footer_line()
        .to_string()
        .contains("y title of the counts"));
    typed(&mut p, "TITLE2");
    assert!(matches!(p.mode, Mode::Browse));
    let labels: Vec<String> = (crate::figure::labels(&p.screen_panel(&plan, sheet).unwrap()))
        .into_iter()
        .map(|l| l.1)
        .collect();
    for want in ["TITLE1", "TITLE2"] {
        assert!(labels.iter().any(|l| l == want), "{want}: {labels:?}");
    }
    // In 2D, the focused pair's axes; Esc keeps what was there.
    press(&mut p, KeyCode::Char('v'));
    let (x, y, _) = p.view().shown_pairs()[p.pair_at].clone();
    press(&mut p, KeyCode::Char('t'));
    assert!(p
        .footer_line()
        .to_string()
        .contains(&format!("x title of {}", x.short())));
    typed(&mut p, "TITLE3");
    assert!(p
        .footer_line()
        .to_string()
        .contains(&format!("y title of {}", y.short())));
    press(&mut p, KeyCode::Esc);
    assert!(matches!(p.mode, Mode::Browse));
    let pair = p.screen_panel(&plan, sheet).unwrap();
    assert!(!pair.contains("TITLE3"), "Esc leaves the titles");
    press(&mut p, KeyCode::Char('t'));
    typed(&mut p, "TITLE3");
    typed(&mut p, "TITLE4");
    let pair = p.screen_panel(&plan, sheet).unwrap();
    assert!(pair.contains("TITLE3") && pair.contains("TITLE4"));
}
