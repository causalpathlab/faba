//! A figure's size on the page, in inches, as every tool that writes one
//! takes it (`--fig-width`, `--fig-height`; MCP's `width` and `height`):
//! either or both, the other the figure's own. A figure given a size is laid
//! out for it, its text in points as it prints; one given none keeps the
//! size it has on its own.

use clap::Args;

/// Points per inch: a figure's SVG is drawn in points.
pub const PT: f64 = 72.0;

/// The smallest figure laid out, in inches: its labels need the room.
pub const MIN: (f64, f64) = (2.0, 1.5);

/// `--fig-width` and `--fig-height`.
#[derive(Args, Debug, Clone, Copy, Default, PartialEq, serde::Serialize)]
pub struct SizeArgs {
    #[arg(
        long = "fig-width",
        value_name = "INCHES",
        help = "The figure's width in inches (e.g. 3.5 for a column, 7 for a page); laid out for it, its text kept readable"
    )]
    pub fig_width: Option<f64>,

    #[arg(
        long = "fig-height",
        value_name = "INCHES",
        help = "The figure's height in inches; alone or with --fig-width, the other the figure's own"
    )]
    pub fig_height: Option<f64>,
}

impl SizeArgs {
    /// Whether a size is asked for.
    pub fn given(&self) -> bool {
        self.fig_width.is_some() || self.fig_height.is_some()
    }

    /// The sides given, in inches, each checked: a number, at least `min`
    /// inches. For a figure drawn in inches by a renderer of its own.
    pub fn inches(&self, min: (f64, f64)) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        for (side, v, least) in [
            ("--fig-width", self.fig_width, min.0),
            ("--fig-height", self.fig_height, min.1),
        ] {
            if let Some(v) = v {
                anyhow::ensure!(
                    v.is_finite() && v >= least,
                    "{side} {v}: a figure is at least {} by {} inches",
                    min.0,
                    min.1
                );
            }
        }
        Ok((self.fig_width, self.fig_height))
    }

    /// The size asked for in points: each side given (see
    /// [`SizeArgs::inches`]), else the figure's own (`natural`, in points);
    /// none when neither is given.
    pub fn points(
        &self,
        natural: impl FnOnce() -> (f64, f64),
        min: (f64, f64),
    ) -> anyhow::Result<Option<(f64, f64)>> {
        self.inches(min)?;
        if !self.given() {
            return Ok(None);
        }
        // The figure's own size, worked out only for a side not given.
        let own = match (self.fig_width, self.fig_height) {
            (Some(_), Some(_)) => (0.0, 0.0),
            _ => natural(),
        };
        Ok(Some((
            self.fig_width.map_or(own.0, |w| w * PT),
            self.fig_height.map_or(own.1, |h| h * PT),
        )))
    }

    /// The figure `draw` makes, at the size asked for (see
    /// [`SizeArgs::points`]): `draw(None)` its own size, `draw(Some(..))`
    /// laid out for one in points, a side not given taken from its own.
    pub fn draw(
        &self,
        min: (f64, f64),
        draw: impl Fn(Option<(f64, f64)>) -> anyhow::Result<String>,
    ) -> anyhow::Result<String> {
        // The figure's own size, when a side is missing: drawn once to learn
        // it (a drawing that fails is refused below, at the size asked for).
        let own = || draw(None).map_or((0.0, 0.0), |svg| super::svg_size(&svg));
        draw(self.points(own, min)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn size(w: Option<f64>, h: Option<f64>) -> SizeArgs {
        SizeArgs {
            fig_width: w,
            fig_height: h,
        }
    }

    /// Neither side: the figure's own; one side: the other its own; both:
    /// as given, in points; a side below the least or not a number refused.
    #[test]
    fn a_size_is_either_side_the_other_the_figures_own() {
        let own = || (500.0, 300.0);
        let min = (2.0, 1.5);
        assert_eq!(size(None, None).points(own, min).unwrap(), None);
        assert_eq!(
            size(Some(3.5), None).points(own, min).unwrap(),
            Some((252.0, 300.0))
        );
        assert_eq!(
            size(None, Some(2.0)).points(own, min).unwrap(),
            Some((500.0, 144.0))
        );
        assert_eq!(
            size(Some(7.0), Some(4.0)).points(own, min).unwrap(),
            Some((504.0, 288.0))
        );
        for bad in [
            size(Some(1.0), None),
            size(None, Some(0.0)),
            size(Some(f64::NAN), None),
        ] {
            let e = bad.points(own, min).unwrap_err().to_string();
            assert!(e.contains("at least 2 by 1.5 inches"), "{e}");
        }
    }

    /// `draw` lays the figure out once at its own size to learn it when a
    /// side is missing, then for the size asked for.
    #[test]
    fn draw_asks_its_own_size_only_for_a_missing_side() {
        let calls = std::cell::RefCell::new(Vec::new());
        let draw = |s: Option<(f64, f64)>| {
            calls.borrow_mut().push(s);
            let (w, h) = s.unwrap_or((400.0, 200.0));
            Ok(super::super::svg(w, h, |_| {}))
        };
        size(Some(3.0), None).draw((2.0, 1.5), draw).unwrap();
        assert_eq!(*calls.borrow(), [None, Some((216.0, 200.0))]);
        calls.borrow_mut().clear();
        size(None, None).draw((2.0, 1.5), draw).unwrap();
        assert_eq!(*calls.borrow(), [None]);
        calls.borrow_mut().clear();
        size(Some(3.0), Some(2.0)).draw((2.0, 1.5), draw).unwrap();
        assert_eq!(
            *calls.borrow(),
            [Some((216.0, 144.0))],
            "both given: drawn once"
        );
    }
}
