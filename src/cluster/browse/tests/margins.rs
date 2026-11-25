//! The modules grid's margins and extras: the modules' sizes and reads as
//! bars beside the grid, in its order; GO terms on its right; the key kept
//! under `--panel-only`; the blocks' values; and the least size a refusal
//! names, which is itself drawn.

use super::super::figure::{
    draw, draw_with, go_table, GridArgs, Side, READ_OFF, READ_UP, SIZE_BAR,
};
use super::{loaded, Count};
use crate::figure::style::{Palette, Sheet, StyleArgs};

/// 4.6 by 2.6 inches, in points.
const PAPER: (Option<f64>, Option<f64>) = (Some(331.2), Some(187.2));
const F: f64 = 7.0;

/// The paper's panel: 7-point text, the plot alone, in `palette`.
fn paper(palette: Palette) -> Sheet {
    Sheet {
        page: PAPER,
        look: StyleArgs {
            font_size: Some(F),
            panel_only: true,
            palette,
        },
    }
}

/// `key="..."` of an element's line.
fn attr<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let at = line.find(&format!(" {key}=\""))? + key.len() + 3;
    line[at..].split('"').next()
}

fn num(line: &str, key: &str) -> f64 {
    attr(line, key).unwrap().parse().unwrap()
}

#[derive(Debug, Clone)]
struct Box_ {
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    fill: String,
}

fn rects(svg: &str) -> Vec<Box_> {
    (svg.lines())
        .filter(|l| l.starts_with("<rect"))
        .map(|l| Box_ {
            x: num(l, "x"),
            y: num(l, "y"),
            w: num(l, "width"),
            h: num(l, "height"),
            fill: attr(l, "fill").unwrap().to_string(),
        })
        .collect()
}

#[derive(Debug, Clone)]
struct Text {
    x: f64,
    y: f64,
    size: f64,
    anchor: String,
    s: String,
}

fn texts(svg: &str) -> Vec<Text> {
    (svg.lines())
        .filter(|l| l.starts_with("<text x="))
        .map(|l| Text {
            x: num(l, "x"),
            y: num(l, "y"),
            size: num(l, "font-size"),
            anchor: attr(l, "text-anchor").unwrap().to_string(),
            s: l[l.find('>').unwrap() + 1..l.rfind("</text>").unwrap()]
                .replace("&lt;", "<")
                .replace("&gt;", ">")
                .replace("&amp;", "&"),
        })
        .collect()
}

/// The grid's geometry off its labels: each column's centre (by the
/// columns' order), each row's centre (by the rows'), and the cell.
struct Grid {
    cols: Vec<f64>,
    rows: Vec<f64>,
    cell: f64,
}

fn grid_of(svg: &str, m: &super::Modules) -> Grid {
    let t = texts(svg);
    let cols: Vec<f64> = (m.columns.iter())
        .map(|&c| t.iter().find(|t| t.s == format!("C{}", c + 1)).unwrap().x)
        .collect();
    let rows: Vec<f64> = (m.diagonal.iter())
        .map(|&g| {
            let label = (t.iter().find(|t| t.s.starts_with(&format!("M{} (", g + 1)))).unwrap();
            label.y - label.size / 2.5
        })
        .collect();
    let cell = match (cols.len(), rows.len()) {
        (2.., _) => cols[1] - cols[0],
        (_, 2..) => rows[1] - rows[0],
        _ => panic!("one block"),
    };
    Grid { cols, rows, cell }
}

impl Grid {
    fn left(&self) -> f64 {
        self.cols[0] - self.cell / 2.0
    }
    fn right(&self) -> f64 {
        self.cols[self.cols.len() - 1] + self.cell / 2.0
    }
    fn top(&self) -> f64 {
        self.rows[0] - self.cell / 2.0
    }
    fn bottom(&self) -> f64 {
        self.rows[self.rows.len() - 1] + self.cell / 2.0
    }
}

/// Whether `lens` are in proportion to `values`, one scale for all.
fn proportional(lens: &[f64], values: &[f64]) -> bool {
    let scale = lens
        .iter()
        .zip(values)
        .map(|(l, v)| l / v)
        .collect::<Vec<_>>();
    lens.len() == values.len() && scale.iter().all(|s| (s - scale[0]).abs() < 0.02 * scale[0])
}

fn grid(f: impl FnOnce(&mut GridArgs)) -> GridArgs {
    let mut g = GridArgs::default();
    f(&mut g);
    g
}

/// The size bars: a gene module's on the left of its row, a cell module's
/// over (or under) its column, in the grid's own order, each as long as
/// its size against the largest; the figure at the size asked, under each
/// look. Both fixtures, the cell modules shown in either order.
#[test]
fn size_bars_follow_the_grids_order_at_lengths_by_the_fits_sizes() {
    for m in [loaded(), super::swapped()] {
        for look in [
            super::super::figure::Look::Shift,
            super::super::figure::Look::Hinton,
        ] {
            for side in [Side::Top, Side::Bottom] {
                let g = grid(|g| g.size_bars = vec![Side::Left, side]);
                let svg = draw_with(&m, look, &g, None, "t", paper(Palette::Colour)).unwrap();
                assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
                let at = grid_of(&svg, &m);
                let bars: Vec<Box_> = (rects(&svg).into_iter())
                    .filter(|r| r.fill == SIZE_BAR)
                    .collect();
                let (mut left, mut cells): (Vec<Box_>, Vec<Box_>) =
                    bars.into_iter().partition(|r| r.x + r.w <= at.left());
                left.sort_by(|a, b| a.y.total_cmp(&b.y));
                cells.sort_by(|a, b| a.x.total_cmp(&b.x));
                let genes: Vec<f64> = m.diagonal.iter().map(|&g| m.gene_sizes[g] as f64).collect();
                let sizes: Vec<f64> = m.columns.iter().map(|&c| m.cell_sizes[c]).collect();
                let lens: Vec<f64> = left.iter().map(|r| r.w).collect();
                assert!(proportional(&lens, &genes), "{lens:?} by {genes:?}");
                let lens: Vec<f64> = cells.iter().map(|r| r.h).collect();
                assert!(proportional(&lens, &sizes), "{lens:?} by {sizes:?}");
                // Each beside its own row or column, growing from the grid.
                for (r, y) in left.iter().zip(&at.rows) {
                    assert!((r.y + r.h / 2.0 - y).abs() < 0.05, "{r:?} at {y}");
                }
                let edge = left[0].x + left[0].w;
                assert!(left.iter().all(|r| (r.x + r.w - edge).abs() < 0.02));
                for (r, x) in cells.iter().zip(&at.cols) {
                    assert!((r.x + r.w / 2.0 - x).abs() < 0.05, "{r:?} at {x}");
                    match side {
                        Side::Top => assert!(r.y + r.h < at.top() - F),
                        _ => assert!(r.y > at.bottom()),
                    }
                }
                let base = |r: &Box_| if side == Side::Top { r.y + r.h } else { r.y };
                assert!(cells
                    .iter()
                    .all(|r| (base(r) - base(&cells[0])).abs() < 0.02));
                // The largest of each labelled, at the font size.
                let t = texts(&svg);
                let big = m.diagonal.iter().map(|&g| m.gene_sizes[g]).max().unwrap();
                assert!(t.iter().any(|t| t.s == format!("{big} genes")), "{svg}");
                let most = sizes.iter().copied().fold(0.0, f64::max);
                assert!(
                    t.iter()
                        .any(|t| t.s == format!("{} cells", m.cells_label(most))),
                    "{svg}"
                );
            }
        }
    }
}

/// `--marginal-size` is a bar's depth, the largest bar's length: that
/// fraction of the grid's shorter side.
#[test]
fn the_marginal_size_sets_the_bars_depth() {
    let m = loaded();
    let mut depths = Vec::new();
    for s in [0.1, 0.3] {
        let g = grid(|g| {
            g.size_bars = vec![Side::Left, Side::Top];
            g.marginal_size = s;
        });
        let svg = draw_with(
            &m,
            super::super::figure::Look::Shift,
            &g,
            None,
            "t",
            paper(Palette::Colour),
        )
        .unwrap();
        assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
        let at = grid_of(&svg, &m);
        let shorter = at.cell * m.columns.len().min(m.diagonal.len()) as f64;
        let bars: Vec<Box_> = (rects(&svg).into_iter())
            .filter(|r| r.fill == SIZE_BAR)
            .collect();
        let (left, top): (Vec<&Box_>, Vec<&Box_>) =
            bars.iter().partition(|r| r.x + r.w <= at.left());
        let left = left.iter().map(|r| r.w).fold(0.0, f64::max);
        let top = top.iter().map(|r| r.h).fold(0.0, f64::max);
        for depth in [left, top] {
            assert!(
                (depth - s * shorter).abs() < 0.05,
                "{depth} of {shorter} at {s}"
            );
        }
        depths.push(left);
    }
    assert!(depths[1] > depths[0] * 2.0, "{depths:?}");
}

/// The read bars: each module's reads, summed over its blocks in the
/// grid, as the browser's mirror bars draw them: converted reads on one
/// side of a zero line, unconverted on the other, one scale across a
/// side's modules, in the grid's order.
#[test]
fn read_bars_draw_each_modules_reads_on_one_scale() {
    for m in [loaded(), super::swapped()] {
        for side in [Side::Top, Side::Bottom] {
            let g = grid(|g| g.read_bars = vec![Side::Left, side]);
            let look = super::super::figure::Look::Hinton;
            let svg = draw_with(&m, look, &g, None, "t", paper(Palette::Colour)).unwrap();
            assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
            let at = grid_of(&svg, &m);
            let gene: Vec<Count> = (m.diagonal.iter())
                .map(|&g| m.block(g).iter().copied().sum())
                .collect();
            let cell: Vec<Count> = (0..m.columns.len())
                .map(|j| m.diagonal.iter().map(|&g| m.block(g)[j]).sum())
                .collect();
            let of = |fill: &str| -> (Vec<Box_>, Vec<Box_>) {
                let (mut l, mut c): (Vec<Box_>, Vec<Box_>) = (rects(&svg).into_iter())
                    .filter(|r| r.fill == fill)
                    .partition(|r| r.x + r.w <= at.left());
                l.sort_by(|a, b| a.y.total_cmp(&b.y));
                c.sort_by(|a, b| a.x.total_cmp(&b.x));
                (l, c)
            };
            let ((lu, cu), (ld, cd)) = (of(READ_UP), of(READ_OFF));
            // One scale across each side's modules, up and down alike.
            let lens: Vec<f64> = lu.iter().chain(&ld).map(|r| r.w).collect();
            let values: Vec<f64> = gene
                .iter()
                .map(|c| c.up())
                .chain(gene.iter().map(|c| c.off()))
                .collect();
            assert!(proportional(&lens, &values), "{lens:?} by {values:?}");
            let lens: Vec<f64> = cu.iter().chain(&cd).map(|r| r.h).collect();
            let values: Vec<f64> = cell
                .iter()
                .map(|c| c.up())
                .chain(cell.iter().map(|c| c.off()))
                .collect();
            assert!(proportional(&lens, &values), "{lens:?} by {values:?}");
            // Converted outward (left) and up, unconverted the other way,
            // from one zero line; each by its own row or column.
            let zero = lu[0].x + lu[0].w;
            for (u, d) in lu.iter().zip(&ld) {
                assert!((u.x + u.w - zero).abs() < 0.02 && (d.x - zero).abs() < 0.02);
            }
            for ((u, d), y) in lu.iter().zip(&ld).zip(&at.rows) {
                assert!((u.y + u.h / 2.0 - y).abs() < 0.05 && (d.y + d.h / 2.0 - y).abs() < 0.05);
            }
            let zero = cu[0].y + cu[0].h;
            for ((u, d), x) in cu.iter().zip(&cd).zip(&at.cols) {
                assert!((u.y + u.h - zero).abs() < 0.02 && (d.y - zero).abs() < 0.02);
                assert!((u.x + u.w / 2.0 - x).abs() < 0.05 && (d.x + d.w / 2.0 - x).abs() < 0.05);
            }
            match side {
                Side::Top => assert!(cd.iter().all(|r| r.y + r.h < at.top())),
                _ => assert!(cu.iter().all(|r| r.y > at.bottom())),
            }
        }
    }
    // In the gray palette, dark and light gray.
    let g = grid(|g| g.read_bars = vec![Side::Left]);
    let svg = draw_with(
        &loaded(),
        super::super::figure::Look::Shift,
        &g,
        None,
        "t",
        paper(Palette::Gray),
    )
    .unwrap();
    assert!(
        svg.contains("fill=\"#262626\"") && svg.contains("fill=\"#b3b3b3\""),
        "{svg}"
    );
}

/// `--key` keeps a key under `--panel-only`, under the grid, nothing
/// right of it (the side the GO terms take), every label at the font
/// size; for the Hinton look, squares at reference |z| too, each its
/// area as the grid's squares are.
#[test]
fn the_key_is_kept_under_the_panel_alone() {
    use super::super::figure::Look;
    let mut referenced = 0;
    for (m, look) in [Look::Shift, Look::Percent, Look::Z, Look::Hinton]
        .into_iter()
        .map(|look| (loaded(), look))
        .chain([(super::wide(), Look::Hinton)])
    {
        let none = draw_with(
            &m,
            look,
            &GridArgs::default(),
            None,
            "t",
            paper(Palette::Colour),
        )
        .unwrap();
        let keyed = draw_with(
            &m,
            look,
            &grid(|g| g.key = true),
            None,
            "t",
            paper(Palette::Colour),
        )
        .unwrap();
        assert_eq!(crate::figure::svg_size(&keyed), (331.2, 187.2));
        assert!(rects(&keyed).len() > rects(&none).len() + 10, "{look:?}");
        let at = grid_of(&keyed, &m);
        let t = texts(&keyed);
        assert!(t.iter().all(|t| t.size == F), "{look:?}");
        let new: Vec<&Text> = t
            .iter()
            .filter(|t| !texts(&none).iter().any(|n| n.s == t.s))
            .collect();
        assert!(!new.is_empty(), "{look:?}");
        assert!(new.iter().all(|t| t.y > at.bottom()), "{look:?}: {new:?}");
        assert!(
            rects(&keyed)
                .iter()
                .skip(1)
                .all(|r| r.x + r.w <= at.right() + 0.01),
            "{look:?}"
        );
        if look == Look::Hinton {
            let z: Vec<f64> = (m.diagonal.iter())
                .flat_map(|&g| super::super::standardized(m.block(g)))
                .flatten()
                .collect();
            let top = z.iter().fold(0.0f64, |a, v| a.max(v.abs()));
            let squares: Vec<Box_> = (rects(&keyed).into_iter())
                .filter(|r| r.y > at.bottom() && (r.w - r.h).abs() < 1e-9 && r.w > 1.0)
                .collect();
            assert!(
                !squares.is_empty() && t.iter().any(|t| t.s == "|z|"),
                "{keyed}"
            );
            for s in &squares {
                // Labelled under it with its |z|.
                let label = (t.iter())
                    .find(|t| (t.x - (s.x + s.w / 2.0)).abs() < 0.05 && t.y > s.y + s.h)
                    .unwrap_or_else(|| panic!("{s:?} unlabelled"));
                let v: f64 = label.s.parse().unwrap();
                let side = (at.cell - 2.0) * (v / top).sqrt();
                assert!((s.w - side).abs() < 0.02, "{s:?}: |z| {v} of {top}");
                referenced += 1;
            }
        }
    }
    // Squares at more than one |z|, so their areas are told apart.
    assert!(referenced >= 3, "{referenced}");
}

/// With every extra on, at 4.6 by 2.6 inches and 7 points under each
/// look and side: the figure at the size asked, every label at the font
/// size and inside the page.
#[test]
fn every_label_is_at_the_font_size_and_on_the_page() {
    use super::super::figure::Look;
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = super::go_column::with_go(tmp.path());
    let table = go_table(&m, Some(2), &ask).unwrap();
    for look in [Look::Shift, Look::Hinton] {
        for side in [Side::Top, Side::Bottom] {
            let g = grid(|g| {
                g.size_bars = vec![Side::Left, side];
                g.read_bars = vec![Side::Left, side];
                g.key = true;
                g.values = true;
                g.go_labels = Some(2);
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
                    assert!(x0 >= -0.5 && x1 <= 331.7, "{look:?} {side:?}: {t:?}");
                    assert!(t.y >= 0.7 * F && t.y <= 187.2, "{look:?} {side:?}: {t:?}");
                }
            }
        }
    }
}

/// `--go-labels N`: each gene module's row carries its top N terms on the
/// grid's right, in the table's order, each with its q; a module with
/// none says so; GO is read only with the flag.
#[test]
fn go_labels_put_each_modules_top_terms_on_its_row() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, mut ask) = super::go_column::with_go(tmp.path());
    let table = go_table(&m, Some(3), &ask).unwrap().unwrap();
    for n in [1, 2, 3] {
        let g = grid(|g| g.go_labels = Some(n));
        let look = super::super::figure::Look::Hinton;
        let svg = draw_with(&m, look, &g, Some(&table), "t", paper(Palette::Gray)).unwrap();
        assert_eq!(crate::figure::svg_size(&svg), (331.2, 187.2));
        let at = grid_of(&svg, &m);
        let right: Vec<Text> = (texts(&svg).into_iter())
            .filter(|t| t.x > at.right())
            .collect();
        for (&g, &y) in m.diagonal.iter().zip(&at.rows) {
            let mut row: Vec<&Text> = (right.iter())
                .filter(|t| (t.y - y).abs() < at.cell / 2.0 && t.anchor == "start")
                .collect();
            row.sort_by(|a, b| a.y.total_cmp(&b.y));
            let hits = &table.of(g).unwrap().hits;
            if hits.is_empty() {
                assert_eq!(row.len(), 1, "{row:?}");
                let said = row[0].s.trim_end_matches('…');
                assert!(
                    "none at q <= 0.05".starts_with(said) && said.len() > 3,
                    "{row:?}"
                );
                continue;
            }
            let want: Vec<&str> = hits.iter().take(n).map(|h| h.named()).collect();
            let got: Vec<&str> = row.iter().map(|t| t.s.trim_end_matches('…')).collect();
            assert_eq!(got.len(), want.len(), "M{}: {got:?}", g + 1);
            for (got, want) in got.iter().zip(&want) {
                assert!(want.starts_with(got) && got.len() > 3, "{got} of {want}");
            }
            let qs: Vec<&Text> = (right.iter())
                .filter(|t| (t.y - y).abs() < at.cell / 2.0 && t.anchor == "end")
                .collect();
            let said: Vec<String> = hits
                .iter()
                .take(n)
                .map(|h| format!("{:.1e}", h.q))
                .collect();
            let mut qs: Vec<(f64, String)> = qs.iter().map(|t| (t.y, t.s.clone())).collect();
            qs.sort_by(|a, b| a.0.total_cmp(&b.0));
            assert_eq!(qs.into_iter().map(|q| q.1).collect::<Vec<_>>(), said);
        }
    }
    // Not read without the flag: a source that is not there is no error.
    ask.gmt = Some(tmp.path().join("none.gmt").to_string_lossy().into());
    assert!(go_table(&m, None, &ask).unwrap().is_none());
    assert!(go_table(&m, Some(2), &ask).is_err());
}

/// `--values`: each block's value in its cell, at the font size, as the
/// look has it (z for the Hinton and z looks, the shift, the %), at the
/// cell's centre.
#[test]
fn values_are_printed_in_each_block() {
    use super::super::figure::Look;
    let m = loaded();
    for look in [Look::Hinton, Look::Z, Look::Shift, Look::Percent] {
        let svg = draw_with(
            &m,
            look,
            &grid(|g| g.values = true),
            None,
            "t",
            paper(Palette::Colour),
        )
        .unwrap();
        let at = grid_of(&svg, &m);
        let t = texts(&svg);
        for (i, &g) in m.diagonal.iter().enumerate() {
            let z = super::super::standardized(m.block(g));
            for (j, &c) in m.columns.iter().enumerate() {
                let want = match look {
                    Look::Hinton | Look::Z => z[j].map(|v| format!("{v:+.1}")),
                    Look::Shift => Some(format!("{:+.2}", m.shift(c, g))),
                    Look::Percent => m.block(g)[j].number(),
                };
                let Some(want) = want else { continue };
                let (x, y) = (at.cols[j], at.rows[i] + F / 2.5);
                let hit = t
                    .iter()
                    .find(|t| (t.x - x).abs() < 0.05 && (t.y - y).abs() < 0.05);
                let hit =
                    hit.unwrap_or_else(|| panic!("{look:?}: no value at M{} C{}", g + 1, c + 1));
                assert_eq!((hit.s.as_str(), hit.size), (want.as_str(), F), "{look:?}");
            }
        }
    }
    // Without it, the Hinton squares carry none.
    let svg = draw(
        &m,
        super::super::figure::Plot::Modules,
        Look::Hinton,
        "t",
        paper(Palette::Colour),
        &Default::default(),
    )
    .unwrap();
    assert!(!texts(&svg).iter().any(|t| t.s.starts_with('+')), "{svg}");
}

/// Every number in `s`.
fn numbers(s: &str) -> Vec<f64> {
    (s.split(|c: char| !(c.is_ascii_digit() || c == '.')))
        .filter_map(|w| w.trim_matches('.').parse().ok())
        .collect()
}

/// A figure refused as too small names the least it needs, in the sides
/// asked for; at that size it is drawn (never a size itself refused).
#[test]
fn a_refusal_names_a_size_that_is_drawn() {
    use super::super::figure::{Look, Plot};
    let batched = super::batches::batched();
    let wide = super::wide();
    // The small pages, and a wide grid's larger ones, where the text set
    // for the page grows with it (and needs more room).
    let small = [
        (Some(60.0), None),
        (Some(60.0), Some(50.0)),
        (None, Some(40.0)),
    ];
    let large = [
        (Some(200.0), None),
        (Some(200.0), Some(150.0)),
        (None, Some(30.0)),
    ];
    for (m, plot, pages) in [
        (&loaded(), Plot::Modules, small),
        (&batched, Plot::Batches, small),
        (&wide, Plot::Modules, large),
    ] {
        for font_size in [None, Some(7.0)] {
            for page in pages {
                let sheet = |page| Sheet {
                    page,
                    look: StyleArgs {
                        font_size,
                        ..Default::default()
                    },
                };
                let e =
                    draw(m, plot, Look::Shift, "t", sheet(page), &Default::default()).unwrap_err();
                assert!(e.contains("needs at least"), "{e}");
                let said: Vec<f64> = numbers(&e).into_iter().map(|v| v * 72.0).collect();
                let least = match page {
                    (Some(_), Some(_)) => (Some(said[0]), Some(said[1])),
                    (Some(_), None) => (Some(said[0]), None),
                    _ => (None, Some(said[0])),
                };
                let drawn = draw(m, plot, Look::Shift, "t", sheet(least), &Default::default());
                assert!(
                    drawn.is_ok(),
                    "{plot:?} {font_size:?} {page:?}: {e}, then {drawn:?}"
                );
                // And the least: a twentieth of an inch less on each side
                // asked for is refused.
                let less = |s: Option<f64>| s.map(|s| s - 3.6);
                let under = (less(least.0), less(least.1));
                let refused = draw(m, plot, Look::Shift, "t", sheet(under), &Default::default());
                assert!(
                    refused.is_err(),
                    "{plot:?} {page:?}: {e}, yet {under:?} drawn"
                );
            }
        }
    }
}

/// The flags parse, only with `--figure`: `--size-bars` bare is left and
/// top; top and bottom together are refused; `--marginal-size` is in
/// (0, 1].
#[test]
fn the_grids_flags_parse_with_a_figure() {
    let parse = |words: &[&str]| {
        let line = ["faba", "cluster-view", "out/"].iter().chain(words);
        let m =
            crate::matches(&mut crate::faba_command(), line.copied()).map_err(|e| e.to_string())?;
        let cli = crate::parse(&m).map_err(|e| e.to_string())?;
        match cli.commands {
            crate::Commands::ClusterView(a) => Ok::<_, String>(a.grid),
            _ => unreachable!(),
        }
    };
    let g = parse(&["--figure", "f", "--size-bars"]).unwrap();
    assert_eq!(g.size_bars, [Side::Left, Side::Top]);
    let g = parse(&[
        "--figure",
        "f",
        "--size-bars",
        "left,bottom",
        "--read-bars",
        "top",
    ])
    .unwrap();
    assert_eq!(
        (g.size_bars, g.read_bars),
        (vec![Side::Left, Side::Bottom], vec![Side::Top])
    );
    let g = parse(&[
        "--figure",
        "f",
        "--marginal-size",
        "0.25",
        "--key",
        "--values",
        "--go-labels",
        "2",
    ])
    .unwrap();
    assert_eq!(
        (g.marginal_size, g.key, g.values, g.go_labels),
        (0.25, true, true, Some(2))
    );
    assert_eq!(parse(&["--figure", "f"]).unwrap(), GridArgs::default());
    let both = parse(&["--figure", "f", "--size-bars", "top,bottom"]).unwrap();
    assert!(both
        .check()
        .unwrap_err()
        .to_string()
        .contains("top or bottom"));
    for bad in [
        &["--marginal-size", "0"][..],
        &["--marginal-size", "1.5"],
        &["--go-labels", "0"],
    ] {
        let with: Vec<&str> = bad.iter().copied().chain(["--figure", "f"]).collect();
        assert!(parse(&with).is_err(), "{bad:?}");
    }
    for alone in [
        &["--size-bars", "left"][..],
        &["--read-bars", "left"],
        &["--key"],
        &["--values"],
        &["--go-labels", "2"],
        &["--marginal-size", "0.2"],
    ] {
        assert!(parse(alone).unwrap_err().contains("--figure"), "{alone:?}");
    }
}
