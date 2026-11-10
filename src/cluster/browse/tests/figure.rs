//! The fit's figures: the modules grid under each look, laid out for a size
//! when one is given; the trace and the batches refused when there are none.

use super::super::figure::{draw, Look, Plot};
use super::loaded;

/// The SVG's width attribute, in points.
fn svg_width(svg: &str) -> f64 {
    let at = svg.find("width=\"").unwrap() + 7;
    svg[at..at + svg[at..].find('"').unwrap()].parse().unwrap()
}

#[test]
fn the_modules_grid_is_drawn_under_each_look_and_for_a_size() {
    let m = loaded();
    for look in [Look::Shift, Look::Percent, Look::Z] {
        let svg = draw(&m, Plot::Modules, look, "a fit", None).unwrap();
        assert!(
            svg.contains("M1") && svg.contains("C1") && svg.contains("C2"),
            "{look:?}"
        );
    }
    let shift = draw(&m, Plot::Modules, Look::Shift, "a fit", None).unwrap();
    assert!(shift.contains("shift (logit"));
    // 3.5 inches wide is 252 points, whatever the modules.
    let sized = draw(&m, Plot::Modules, Look::Shift, "a fit", Some((3.5, 3.0))).unwrap();
    assert!(
        (svg_width(&sized) - 252.0).abs() < 1e-9,
        "{}",
        svg_width(&sized)
    );
}

#[test]
fn a_figure_with_nothing_to_draw_says_why() {
    let m = loaded();
    let err = draw(&m, Plot::Batches, Look::Shift, "a fit", None).unwrap_err();
    assert!(err.contains("one batch"), "{err}");
    assert!(draw(&m, Plot::Trace, Look::Shift, "a fit", None).is_err());
}
