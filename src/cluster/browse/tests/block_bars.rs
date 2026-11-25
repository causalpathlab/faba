//! The modules grid's blocks as the browser's mirror bars: `--look bars`
//! draws each block's converted reads rising above a zero line and its
//! unconverted reads hanging below, each row on its own scale (as the
//! browser's rows are), the rows in
//! the browser's order (by shift range) with each row's range on its
//! right; `--block-bars` draws the same bars over another look's blocks.

use super::super::figure::{draw_with, Look, BLOCK_OFF, BLOCK_OFF_GRAY, BLOCK_UP, BLOCK_UP_GRAY};
use super::margins::{grid, grid_by, paper, rects, texts, Box_, Grid, F};
use super::{loaded, Modules};
use crate::figure::style::Palette;

/// Whether `lens` are in proportion to `values`, one scale for all.
fn proportional(lens: &[f64], values: &[f64]) -> bool {
    let scale: Vec<f64> = lens.iter().zip(values).map(|(l, v)| l / v).collect();
    lens.len() == values.len()
        && !scale.is_empty()
        && scale.iter().all(|s| (s - scale[0]).abs() < 0.02 * scale[0])
}

/// The bars in the grid's cells `at`, rows the gene modules `order`: each
/// block's converted reads up from its row's zero line (its centre) in
/// `up`, its unconverted reads down from it in `off`, centred in the
/// block's column, one scale across the row, the tallest side of its
/// tallest block as tall as every row's. The bars' indices in the SVG.
fn bars_are_the_blocks_reads(
    svg: &str,
    m: &Modules,
    order: &[usize],
    at: &Grid,
    (up, off): (&str, &str),
) -> Vec<usize> {
    let all = rects(svg);
    let inside = |r: &Box_| {
        let cx = r.x + r.w / 2.0;
        cx > at.left() && cx < at.right()
    };
    let (mut found, mut tallest) = (Vec::new(), Vec::new());
    for (i, &g) in order.iter().enumerate() {
        let (mut lens, mut values) = (Vec::new(), Vec::new());
        let zero = at.rows[i];
        for (j, r) in m.block(g).iter().enumerate() {
            let x = at.cols[j];
            for (fill, v, rising) in [(up, r.up(), true), (off, r.off(), false)] {
                let hit = (all.iter().enumerate()).find(|(_, b)| {
                    let end = if rising { b.y + b.h } else { b.y };
                    b.fill == fill
                        && inside(b)
                        && (b.x + b.w / 2.0 - x).abs() < 0.05
                        && (end - zero).abs() < 0.02
                });
                match (hit, v > 0.0) {
                    (Some((k, b)), true) => {
                        lens.push(b.h);
                        values.push(v);
                        found.push(k);
                    }
                    (None, false) => {}
                    (hit, _) => panic!("M{} column {j}: {v} reads, {hit:?}", g + 1),
                }
            }
        }
        if !lens.is_empty() {
            assert!(
                proportional(&lens, &values),
                "M{}: {lens:?} by {values:?}",
                g + 1
            );
            tallest.push(lens.iter().copied().fold(0.0, f64::max));
        }
    }
    // Each row's tallest the same height.
    assert!(
        tallest.iter().all(|t: &f64| (t - tallest[0]).abs() < 0.02),
        "{tallest:?}"
    );
    // Nothing else in the cells in the bars' fills.
    let drawn = (all.iter()).filter(|b| (b.fill == up || b.fill == off) && inside(b));
    assert_eq!(drawn.count(), found.len());
    found
}

/// Each row's shift range at its right, as the browser writes it.
fn ranges_are_on_the_right(svg: &str, m: &Modules, order: &[usize], at: &Grid) {
    let t = texts(svg);
    for (i, &g) in order.iter().enumerate() {
        let said = format!("Δ{:.2}", m.shift_range(g));
        let hit = t.iter().find(|t| t.s == said && t.x > at.right());
        let hit = hit.unwrap_or_else(|| panic!("M{}: no {said}", g + 1));
        assert!((hit.y - F / 2.5 - at.rows[i]).abs() < 0.5, "{hit:?}");
        assert_eq!(hit.size, F);
    }
}

fn fills(palette: Palette) -> (&'static str, &'static str) {
    match palette {
        Palette::Gray => (BLOCK_UP_GRAY, BLOCK_OFF_GRAY),
        _ => (BLOCK_UP, BLOCK_OFF),
    }
}

/// `--look bars`: the rows in the browser's order (by shift range, not
/// the other looks' diagonal), the columns in the panes' order, each
/// block a mirror bar of its reads on its row's scale, each row's range
/// on its right.
#[test]
fn the_bars_look_draws_each_blocks_reads_as_the_browser_does() {
    // The browser's order, and the other way round, so it is told apart
    // from the diagonal.
    let mut turned = loaded();
    turned.order.reverse();
    assert_ne!(turned.order, turned.diagonal);
    for m in [loaded(), turned, super::swapped()] {
        for palette in [Palette::Colour, Palette::Gray] {
            let svg = draw_with(&m, Look::Bars, &grid(|_| {}), None, "t", paper(palette)).unwrap();
            assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
            let at = grid_by(&svg, &m, &m.order);
            assert!(at.rows.windows(2).all(|w| w[0] < w[1]), "{:?}", at.rows);
            assert!(at.cols.windows(2).all(|w| w[0] < w[1]), "{:?}", at.cols);
            bars_are_the_blocks_reads(&svg, &m, &m.order, &at, fills(palette));
            ranges_are_on_the_right(&svg, &m, &m.order, &at);
            // The light bars outlined, so they read on white paper.
            let off = fills(palette).1;
            assert!(svg.lines().any(|l| l.starts_with("<rect")
                && l.contains(&format!("fill=\"{off}\""))
                && l.contains("stroke=")));
        }
    }
}

/// `--block-bars` over the Hinton look: its squares kept, sized by |z| as
/// without the bars, behind the same bars, in its own row order.
#[test]
fn block_bars_over_hinton_keep_the_square_behind() {
    let m = loaded();
    for palette in [Palette::Colour, Palette::Gray] {
        let g = grid(|g| g.block_bars = true);
        let svg = draw_with(&m, Look::Hinton, &g, None, "t", paper(palette)).unwrap();
        assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
        let at = grid_by(&svg, &m, &m.diagonal);
        let bars = bars_are_the_blocks_reads(&svg, &m, &m.diagonal, &at, fills(palette));
        ranges_are_on_the_right(&svg, &m, &m.diagonal, &at);
        let z: Vec<Vec<Option<f64>>> = (m.diagonal.iter())
            .map(|&g| super::super::standardized(m.block(g)))
            .collect();
        let top = z
            .iter()
            .flatten()
            .flatten()
            .fold(0.0f64, |a, v| a.max(v.abs()));
        let all = rects(&svg);
        let mut squares = 0;
        for (i, row) in z.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                let Some(v) = v.filter(|v| v.abs() > 1e-9) else {
                    continue;
                };
                let side = (at.cell - 2.0) * (v.abs() / top).sqrt();
                let (cx, cy) = (at.cols[j], at.rows[i]);
                let (k, _) = (all.iter().enumerate())
                    .find(|(_, b)| {
                        b.fill != "none"
                            && (b.w - side).abs() < 0.02
                            && (b.h - side).abs() < 0.02
                            && (b.x + b.w / 2.0 - cx).abs() < 0.05
                            && (b.y + b.h / 2.0 - cy).abs() < 0.05
                    })
                    .unwrap_or_else(|| panic!("no square at {i},{j}"));
                // Its block's bars drawn next, over it.
                assert!(bars.contains(&(k + 1)), "the square behind the bars");
                squares += 1;
            }
        }
        assert!(squares >= 2, "{squares}");
    }
}

/// Every label at the font size and on the page with every extra on, the
/// bars as a look and over another.
#[test]
fn the_bars_keep_every_label_on_the_page() {
    use super::super::figure::Side;
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = super::go_column::with_go(tmp.path());
    let table = super::super::figure::go_table(&m, Some(2), &ask).unwrap();
    for (look, block_bars) in [
        (Look::Bars, false),
        (Look::Hinton, true),
        (Look::Shift, true),
    ] {
        let g = grid(|g| {
            g.size_bars = vec![Side::Left, Side::Top];
            g.read_bars = vec![Side::Left, Side::Top];
            g.key = true;
            g.values = true;
            g.go_labels = Some(2);
            g.block_bars = block_bars;
        });
        for palette in [Palette::Gray, Palette::Colour] {
            let svg = draw_with(&m, look, &g, table.as_ref(), "t", paper(palette)).unwrap();
            assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
            for t in texts(&svg) {
                assert_eq!(t.size, F, "{t:?}");
                let w = t.s.chars().count() as f64 * F * 0.56;
                let (x0, x1) = match t.anchor.as_str() {
                    "start" => (t.x, t.x + w),
                    "end" => (t.x - w, t.x),
                    _ => (t.x - w / 2.0, t.x + w / 2.0),
                };
                assert!(x0 >= -0.5 && x1 <= 331.7, "{look:?}: {t:?}");
                assert!(t.y >= 0.7 * F && t.y <= 187.2, "{look:?}: {t:?}");
            }
        }
    }
}

/// `--look bars` and `--block-bars` parse, the switch only with a figure.
#[test]
fn the_bars_flags_parse() {
    let parse = |words: &[&str]| {
        let line = ["faba", "cluster-view", "out/"].iter().chain(words);
        let m =
            crate::matches(&mut crate::faba_command(), line.copied()).map_err(|e| e.to_string())?;
        crate::parse(&m).map_err(|e| e.to_string())?;
        Ok::<_, String>(())
    };
    assert!(parse(&["--figure", "f", "--look", "bars"]).is_ok());
    assert!(parse(&["--figure", "f", "--look", "hinton", "--block-bars"]).is_ok());
    assert!(parse(&["--block-bars"]).unwrap_err().contains("--figure"));
}
