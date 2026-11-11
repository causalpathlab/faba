//! The fit's figures: the modules grid under each look, laid out for a size
//! when one is given; the trace and the batches refused when there are none.

use super::super::figure::{draw, Look, Plot};
use super::loaded;

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
    // 3.5 inches wide is 252 points, whatever the modules; given alone,
    // the height is the grid's own.
    let size = |w, h| crate::figure::SizeArgs {
        fig_width: w,
        fig_height: h,
    };
    let at = |s: crate::figure::SizeArgs| {
        let drawn = s.draw(super::super::figure::MIN, |size| {
            draw(&m, Plot::Modules, Look::Shift, "a fit", size).map_err(anyhow::Error::msg)
        });
        crate::figure::svg_size(&drawn.unwrap())
    };
    assert_eq!(at(size(Some(3.5), Some(3.0))), (252.0, 216.0));
    let own = crate::figure::svg_size(&shift);
    assert_eq!(at(size(Some(3.5), None)), (252.0, own.1));
    assert_eq!(at(size(None, Some(3.0))), (own.0, 216.0));
    let small = size(Some(1.0), None).draw(super::super::figure::MIN, |_| Ok(String::new()));
    assert!(small.unwrap_err().to_string().contains("at least 2 by 1.5"));
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
    let err = draw(&m, Plot::Batches, Look::Shift, "a fit", None).unwrap_err();
    assert!(err.contains("one batch"), "{err}");
    assert!(draw(&m, Plot::Trace, Look::Shift, "a fit", None).is_err());
}
