//! The fit's figures: the modules grid under each look, laid out for a size
//! when one is given; the trace and the batches refused when there are none.

use super::super::figure::{draw, Look, Plot};
use super::loaded;
use crate::figure::style::{AxisTitles, Palette, Sheet, StyleArgs};

#[test]
fn the_modules_grid_is_drawn_under_each_look_and_for_a_size() {
    let m = loaded();
    for look in [Look::Shift, Look::Percent, Look::Z, Look::Hinton] {
        let svg = draw(
            &m,
            Plot::Modules,
            look,
            "a fit",
            Sheet::OWN,
            &AxisTitles::default(),
        )
        .unwrap();
        assert!(
            svg.contains("M1") && svg.contains("C1") && svg.contains("C2"),
            "{look:?}"
        );
    }
    let shift = draw(
        &m,
        Plot::Modules,
        Look::Shift,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default(),
    )
    .unwrap();
    assert!(shift.contains("shift (logit"));
    // 3.5 inches wide is 252 points, whatever the modules; given alone,
    // the height fits square cells, as tall as the grid's own is wide.
    let at = |page| {
        crate::figure::svg_size(
            &draw(
                &m,
                Plot::Modules,
                Look::Shift,
                "a fit",
                Sheet::of(page),
                &AxisTitles::default(),
            )
            .unwrap(),
        )
    };
    assert_eq!(at((Some(252.0), Some(216.0))), (252.0, 216.0));
    // One side alone: the other in the grid's own proportions.
    let own = crate::figure::svg_size(&shift);
    let (w, h) = at((Some(2.0 * own.0), None));
    assert!((w / h - own.0 / own.1).abs() < 1e-9, "{w}x{h} of {own:?}");
    // Too small for its cells: refused, saying the least.
    let e = draw(
        &m,
        Plot::Modules,
        Look::Shift,
        "a fit",
        Sheet::of((Some(60.0), Some(50.0))),
        &AxisTitles::default(),
    )
    .unwrap_err();
    assert!(e.contains("needs at least"), "{e}");
}

/// The Hinton figure draws each block as a true square, centred in its
/// cell, its area by its |z| against the largest: in gene module 1's row,
/// the block farthest from the row's mean the largest, the nearest the
/// least.
#[test]
fn the_hinton_figure_sizes_each_square_by_its_magnitude() {
    let m = super::wide();
    let svg = draw(
        &m,
        Plot::Modules,
        Look::Hinton,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default(),
    )
    .unwrap();
    let num = |rect: &str, key: &str| -> f64 {
        let at = rect.find(&format!("{key}=\"")).unwrap() + key.len() + 2;
        rect[at..].split('"').next().unwrap().parse().unwrap()
    };
    // Every square, by the row it centres on, left to right.
    let mut squares: Vec<(f64, f64, f64)> = (svg.lines())
        .filter(|l| l.starts_with("<rect") && !l.contains(r#"fill="none""#))
        .map(|r| (num(r, "x"), num(r, "y"), num(r, "width"), num(r, "height")))
        .filter(|&(.., w, h)| (w - h).abs() < 1e-9)
        .map(|(x, y, w, _)| ((y + w / 2.0).round(), x + w / 2.0, w))
        .collect();
    squares.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let first = squares[0].0;
    let row: Vec<f64> = (squares.iter())
        .filter(|s| s.0 == first)
        .map(|s| s.2)
        .collect();
    let g = m.diagonal[0];
    let z: Vec<f64> = super::super::standardized(m.block(g))
        .into_iter()
        .flatten()
        .collect();
    assert_eq!(row.len(), z.len(), "{svg}");
    let (far, near) = (0..z.len()).fold((0, 0), |(f, n), j| {
        (
            if z[j].abs() > z[f].abs() { j } else { f },
            if z[j].abs() < z[n].abs() { j } else { n },
        )
    });
    assert!(row[far] > row[near], "{row:?} by {z:?}");
    assert!(
        row.iter().all(|&w| w <= row[far] && w >= row[near]),
        "{row:?}"
    );
    assert!(svg.contains("square area by |z|"), "{svg}");
}

/// The Hinton figure colours each square as the browser's Hinton look
/// does: red above the row's mean and blue below, never paler than the
/// browser's floor, full red or blue at its top |z|; not the other looks'
/// pale palette.
#[test]
fn the_hinton_figure_colours_each_square_as_the_browser_does() {
    use super::super::view::{Ramp, Stretch};
    use ratatui::style::Color;
    let m = super::wide();
    let svg = draw(
        &m,
        Plot::Modules,
        Look::Hinton,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default(),
    )
    .unwrap();
    let stretch = Stretch::of(m.order.iter().map(|&g| m.block(g)), Ramp::Sqrt);
    let hex = |c: Color| match c {
        Color::Rgb(r, g, b) => format!("#{r:02x}{g:02x}{b:02x}"),
        c => panic!("{c:?}"),
    };
    let browser: Vec<String> = (m.diagonal.iter())
        .flat_map(|&g| super::super::standardized(m.block(g)))
        .flatten()
        .map(|z| hex(stretch.hinton(z).0))
        .collect();
    let squares: Vec<String> = (super::margins::rects(&svg).into_iter())
        .filter(|b| (b.w - b.h).abs() < 1e-9 && b.w > 0.0 && b.fill != "none")
        .map(|b| b.fill)
        .collect();
    assert_eq!(squares.len(), browser.len(), "{svg}");
    for fill in &squares {
        assert!(browser.contains(fill), "{fill} not the browser's");
    }
    assert!(squares.iter().any(|f| f == "#ff0000") && squares.iter().any(|f| f == "#0000ff"));
}

/// Under the Hinton look the colour scale beside the grid runs through the
/// squares' own colours: full red at its top, full blue at its bottom.
#[test]
fn the_hinton_scale_is_the_squares_colours() {
    let m = super::wide();
    let svg = draw(
        &m,
        Plot::Modules,
        Look::Hinton,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default(),
    )
    .unwrap();
    // The scale's steps: narrow, non-square rects, top to bottom.
    let steps: Vec<_> = (super::margins::rects(&svg).into_iter())
        .filter(|b| b.fill != "none" && b.h < b.w && b.h < 5.0)
        .collect();
    assert!(steps.len() >= 20, "{svg}");
    // Saturated as the squares are, not the pastel ramp: one channel full,
    // the other two near zero.
    let rgb = |h: &str| {
        let v = |a: usize| u8::from_str_radix(&h[a..a + 2], 16).unwrap();
        (v(1), v(3), v(5))
    };
    let (r, g, b) = rgb(&steps.first().unwrap().fill);
    assert!(r == 255 && g < 16 && b < 16, "top step full red");
    let (r, g, b) = rgb(&steps.last().unwrap().fill);
    assert!(b == 255 && r < 16 && g < 16, "bottom step full blue");
}

/// Under the Hinton look each square is outlined in gray90 at 0.1 pt;
/// the grid's cells are not.
#[test]
fn each_hinton_square_is_outlined_thin_in_gray90() {
    let m = super::wide();
    let svg = draw(
        &m,
        Plot::Modules,
        Look::Hinton,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default(),
    )
    .unwrap();
    assert!(!svg.contains(r#"fill="none""#), "no cell outlines\n{svg}");
    let squares: Vec<&str> = (svg.lines())
        .filter(|l| l.starts_with("<rect") && (l.contains("#ff") || l.contains("#00")))
        .filter(|l| {
            let w = super::margins::rects(l)[0].w;
            let h = super::margins::rects(l)[0].h;
            (w - h).abs() < 1e-9 && w > 0.0
        })
        .collect();
    assert!(!squares.is_empty(), "{svg}");
    for sq in &squares {
        assert!(sq.contains(r##"stroke="#e5e5e5""##), "{sq}");
        let w: f64 = (sq.split(r#"stroke-width=""#).nth(1).unwrap())
            .split('"')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert!((w - 0.1).abs() < 1e-9, "{sq}");
    }
}

/// `--fig-width` and `--fig-height` (0.19.0's `--width` and `--height`
/// too), and the look (`--font-size`, `--panel-only`, `--palette`), only
/// with `--figure`.
#[test]
fn a_size_goes_with_a_figure() {
    let parse = |words: &[&str]| {
        let line = ["faba", "cluster-view", "out/"].iter().chain(words);
        crate::matches(&mut crate::faba_command(), line.copied()).map_err(|e| e.to_string())
    };
    assert!(parse(&["--figure", "f", "--width", "3", "--height", "2"]).is_ok());
    assert!(parse(&["--figure", "f", "--fig-width", "3"]).is_ok());
    assert!(parse(&["--fig-width", "3"])
        .unwrap_err()
        .contains("--figure"));
    for look in [
        &["--font-size", "6"][..],
        &["--panel-only"],
        &["--palette", "gray"],
    ] {
        assert!(parse(look).unwrap_err().contains("--figure"), "{look:?}");
        let with: Vec<&str> = look.iter().copied().chain(["--figure", "f"]).collect();
        assert!(parse(&with).is_ok(), "{with:?}");
    }
}

/// Every `fill="#rrggbb"` and `stroke="#rrggbb"` of `svg`.
fn colours(svg: &str) -> Vec<&str> {
    ["fill=\"#", "stroke=\"#"]
        .iter()
        .flat_map(|at| {
            svg.match_indices(at)
                .map(move |(i, _)| &svg[i + at.len()..][..6])
        })
        .collect()
}

/// The panel look, on each figure: `panel_only` leaves the title and the
/// key out; `font_size` sets every label's size; `palette` gray draws in
/// grays alone; each at the size asked for.
#[test]
fn each_figure_takes_the_panels_look() {
    let traced = super::trace::load(true, true);
    let batched = super::batches::batched();
    let page = (Some(288.0), Some(216.0));
    for (m, plot) in [
        (&loaded(), Plot::Modules),
        (&traced, Plot::Trace),
        (&batched, Plot::Batches),
    ] {
        let drawn = |look: StyleArgs| {
            let sheet = Sheet { page, look };
            draw(
                m,
                plot,
                Look::Shift,
                "TITLE1",
                sheet,
                &AxisTitles::default(),
            )
            .unwrap()
        };
        let own = drawn(StyleArgs::default());
        assert!(own.contains("TITLE1"), "{plot:?}");
        assert_eq!(crate::figure::svg_size(&own), (288.0, 216.0), "{plot:?}");
        let alone = drawn(StyleArgs {
            panel_only: true,
            ..Default::default()
        });
        assert!(!alone.contains("TITLE1"), "{plot:?}: no title");
        assert_eq!(crate::figure::svg_size(&alone), (288.0, 216.0), "{plot:?}");
        let sized = drawn(StyleArgs {
            font_size: Some(6.5),
            ..Default::default()
        });
        assert!(sized.contains("font-size=\"6.5\""), "{plot:?}");
        assert!(!own.contains("font-size=\"6.5\""), "{plot:?}");
        let gray = drawn(StyleArgs {
            palette: Palette::Gray,
            ..Default::default()
        });
        let hues: Vec<&str> = (colours(&gray).into_iter())
            .filter(|c| !(c[..2] == c[2..4] && c[2..4] == c[4..]))
            .collect();
        assert!(hues.is_empty(), "{plot:?} in gray: {hues:?}");
        assert!(colours(&own).iter().any(|c| c[..2] != c[2..4]), "{plot:?}");
    }
    // The grid's key, the colour scale, is left out with the title.
    let m = loaded();
    let sheet = |panel_only| Sheet {
        page,
        look: StyleArgs {
            panel_only,
            ..Default::default()
        },
    };
    let key = |s: &str| s.contains("shift (logit");
    assert!(key(&draw(
        &m,
        Plot::Modules,
        Look::Shift,
        "t",
        sheet(false),
        &AxisTitles::default(),
    )
    .unwrap()));
    assert!(!key(&draw(
        &m,
        Plot::Modules,
        Look::Shift,
        "t",
        sheet(true),
        &AxisTitles::default(),
    )
    .unwrap()));
}

#[test]
fn a_figure_with_nothing_to_draw_says_why() {
    let m = loaded();
    let err = draw(
        &m,
        Plot::Batches,
        Look::Shift,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default(),
    )
    .unwrap_err();
    assert!(err.contains("one batch"), "{err}");
    assert!(draw(
        &m,
        Plot::Trace,
        Look::Shift,
        "a fit",
        Sheet::OWN,
        &AxisTitles::default()
    )
    .is_err());
}

/// The trace's and the batches' margins do not shrink with their text:
/// a side under 2 by 1.5 inches is refused before any fit is read, saying
/// so; the grid lays itself out and says what it needs as it is drawn.
#[test]
fn the_trace_and_the_batches_keep_their_floor() {
    use super::super::figure::page;
    use crate::figure::SizeArgs;
    for plot in [Plot::Trace, Plot::Batches] {
        for (w, h) in [(Some(1.9), None), (None, Some(1.4))] {
            let e = page(plot, &SizeArgs::new(w, h)).unwrap_err().to_string();
            assert!(e.contains("at least 2.0 by 1.5 inches"), "{plot:?}: {e}");
        }
        assert!(page(plot, &SizeArgs::new(Some(2.0), Some(1.5))).is_ok());
    }
    let small = SizeArgs::new(Some(1.2), Some(0.9));
    assert_eq!(page(Plot::Modules, &small).unwrap(), small.page().unwrap());
}

/// `--xlab` and `--ylab` title the trace's axes in place of its own; the
/// modules grid and the batches have no axis titles, are drawn as they
/// are, and say so in a note.
#[test]
fn the_axis_titles_asked_title_the_trace() {
    use super::super::figure::unused_titles;
    let traced = super::trace::load(true, true);
    let asked = AxisTitles {
        xlab: Some("TEXT1 (-log10)".into()),
        ylab: Some("TEXT2".into()),
    };
    let labels = |titles: &AxisTitles| -> Vec<String> {
        let svg = draw(&traced, Plot::Trace, Look::Shift, "t", Sheet::OWN, titles).unwrap();
        (crate::figure::labels(&svg).into_iter())
            .map(|l| l.1)
            .collect()
    };
    let own = labels(&AxisTitles::default());
    for t in ["k-means iteration, then epoch", "ELBO"] {
        assert!(own.iter().any(|l| l == t), "{t}: {own:?}");
    }
    let titled = labels(&asked);
    for t in ["TEXT1 (-log10)", "TEXT2"] {
        assert!(titled.iter().any(|l| l == t), "{t}: {titled:?}");
    }
    for t in ["k-means iteration, then epoch", "ELBO"] {
        assert!(!titled.iter().any(|l| l == t), "{t}: {titled:?}");
    }
    assert_eq!(unused_titles(Plot::Trace, &asked), None);
    let batched = super::batches::batched();
    for (m, plot) in [(&loaded(), Plot::Modules), (&batched, Plot::Batches)] {
        let at = |titles: &AxisTitles| draw(m, plot, Look::Shift, "t", Sheet::OWN, titles);
        assert_eq!(at(&asked), at(&AxisTitles::default()), "{plot:?}");
        let note = unused_titles(plot, &asked).unwrap();
        assert!(note.contains("--xlab") && note.contains("trace"), "{note}");
        assert_eq!(unused_titles(plot, &AxisTitles::default()), None);
    }
}
