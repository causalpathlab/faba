//! The fit's figures: the modules grid under each look, laid out for a size
//! when one is given; the trace and the batches refused when there are none.

use super::super::figure::{draw, Look, Plot};
use super::loaded;
use crate::figure::size::OWN;

#[test]
fn the_modules_grid_is_drawn_under_each_look_and_for_a_size() {
    let m = loaded();
    for look in [Look::Shift, Look::Percent, Look::Z] {
        let svg = draw(&m, Plot::Modules, look, "a fit", OWN).unwrap();
        assert!(
            svg.contains("M1") && svg.contains("C1") && svg.contains("C2"),
            "{look:?}"
        );
    }
    let shift = draw(&m, Plot::Modules, Look::Shift, "a fit", OWN).unwrap();
    assert!(shift.contains("shift (logit"));
    // 3.5 inches wide is 252 points, whatever the modules; given alone,
    // the height fits square cells, as tall as the grid's own is wide.
    let at = |page| {
        crate::figure::svg_size(&draw(&m, Plot::Modules, Look::Shift, "a fit", page).unwrap())
    };
    assert_eq!(at((Some(252.0), Some(216.0))), (252.0, 216.0));
    // Square cells: a wider page alone makes a taller grid.
    let (narrow, wide) = (at((Some(252.0), None)), at((Some(504.0), None)));
    assert_eq!((narrow.0, wide.0), (252.0, 504.0));
    assert!(wide.1 > narrow.1, "{narrow:?} {wide:?}");
}

/// `--fig-width` and `--fig-height` (0.19.0's `--width` and `--height`
/// too) only with `--figure`.
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
}

#[test]
fn a_figure_with_nothing_to_draw_says_why() {
    let m = loaded();
    let err = draw(&m, Plot::Batches, Look::Shift, "a fit", OWN).unwrap_err();
    assert!(err.contains("one batch"), "{err}");
    assert!(draw(&m, Plot::Trace, Look::Shift, "a fit", OWN).is_err());
}
