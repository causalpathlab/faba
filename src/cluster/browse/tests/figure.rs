//! The fit's figures: the modules grid under each look, laid out for a size
//! when one is given; the trace and the batches refused when there are none.

use super::super::figure::{draw, Look, Plot};
use super::loaded;
use crate::figure::style::{Palette, Sheet, StyleArgs};

#[test]
fn the_modules_grid_is_drawn_under_each_look_and_for_a_size() {
    let m = loaded();
    for look in [Look::Shift, Look::Percent, Look::Z, Look::Hinton] {
        let svg = draw(&m, Plot::Modules, look, "a fit", Sheet::OWN).unwrap();
        assert!(
            svg.contains("M1") && svg.contains("C1") && svg.contains("C2"),
            "{look:?}"
        );
    }
    let shift = draw(&m, Plot::Modules, Look::Shift, "a fit", Sheet::OWN).unwrap();
    assert!(shift.contains("shift (logit"));
    // 3.5 inches wide is 252 points, whatever the modules; given alone,
    // the height fits square cells, as tall as the grid's own is wide.
    let at = |page| {
        crate::figure::svg_size(
            &draw(&m, Plot::Modules, Look::Shift, "a fit", Sheet::of(page)).unwrap(),
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
    let svg = draw(&m, Plot::Modules, Look::Hinton, "a fit", Sheet::OWN).unwrap();
    let num = |rect: &str, key: &str| -> f64 {
        let at = rect.find(&format!("{key}=\"")).unwrap() + key.len() + 2;
        rect[at..].split('"').next().unwrap().parse().unwrap()
    };
    // Every square, by the row it centres on, left to right.
    let mut squares: Vec<(f64, f64, f64)> = (svg.lines())
        .filter(|l| l.starts_with("<rect"))
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
            draw(m, plot, Look::Shift, "TITLE1", sheet).unwrap()
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
        sheet(false)
    )
    .unwrap()));
    assert!(!key(&draw(
        &m,
        Plot::Modules,
        Look::Shift,
        "t",
        sheet(true)
    )
    .unwrap()));
}

#[test]
fn a_figure_with_nothing_to_draw_says_why() {
    let m = loaded();
    let err = draw(&m, Plot::Batches, Look::Shift, "a fit", Sheet::OWN).unwrap_err();
    assert!(err.contains("one batch"), "{err}");
    assert!(draw(&m, Plot::Trace, Look::Shift, "a fit", Sheet::OWN).is_err());
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
