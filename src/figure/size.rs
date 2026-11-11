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

/// A page in points: each side asked for, none where the figure keeps its
/// own (see [`fill`]).
pub type Page = (Option<f64>, Option<f64>);

/// Neither side asked for: the figure's own size.
pub const OWN: Page = (None, None);

/// `page` with each side not asked for taken from `own` (the figure's own
/// size in points, worked out only when a side is missing).
pub fn fill(page: Page, own: impl FnOnce() -> (f64, f64)) -> (f64, f64) {
    match page {
        (Some(w), Some(h)) => (w, h),
        (w, h) => {
            let own = own();
            (w.unwrap_or(own.0), h.unwrap_or(own.1))
        }
    }
}

/// Refused when a side of `page` asked for is less than `least` points
/// (what the figure, `what`, needs for its labels); the sides not asked for
/// are the figure's own, which it lays out to fit.
pub fn at_least(page: Page, least: (f64, f64), what: &str) -> anyhow::Result<()> {
    let short = |side: Option<f64>, least: f64| side.is_some_and(|v| v < least);
    anyhow::ensure!(
        !short(page.0, least.0) && !short(page.1, least.1),
        "{what} is at least {:.1} by {:.1} inches: its labels need the room",
        least.0 / PT,
        least.1 / PT
    );
    Ok(())
}

impl SizeArgs {
    /// The sides given, in inches, each checked: a number, at least
    /// [`MIN`] inches.
    pub fn inches(&self) -> anyhow::Result<(Option<f64>, Option<f64>)> {
        for (side, v, least) in [
            ("--fig-width", self.fig_width, MIN.0),
            ("--fig-height", self.fig_height, MIN.1),
        ] {
            if let Some(v) = v {
                anyhow::ensure!(
                    v.is_finite() && v >= least,
                    "{side} {v}: a figure is at least {} by {} inches",
                    MIN.0,
                    MIN.1
                );
            }
        }
        Ok((self.fig_width, self.fig_height))
    }

    /// The page asked for, in points (see [`SizeArgs::inches`]): a figure
    /// lays itself out for it, filling a side not given with its own.
    pub fn page(&self) -> anyhow::Result<Page> {
        let (w, h) = self.inches()?;
        Ok((w.map(|w| w * PT), h.map(|h| h * PT)))
    }

    /// [`SizeArgs::page`], refused when a side asked for is less than
    /// `least` points (see [`at_least`]): checked before a command's work.
    pub fn page_at_least(&self, least: (f64, f64), what: &str) -> anyhow::Result<Page> {
        let page = self.page()?;
        at_least(page, least, what)?;
        Ok(page)
    }

    /// Of `w` and `h` inches, for a test.
    #[cfg(test)]
    pub fn new(w: Option<f64>, h: Option<f64>) -> Self {
        Self {
            fig_width: w,
            fig_height: h,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Neither side: none; one side or both: in points; a side below the
    /// least or not a number refused. A side not given is the figure's own.
    #[test]
    fn a_size_is_either_side_the_other_the_figures_own() {
        let page = |w, h| SizeArgs::new(w, h).page();
        assert_eq!(page(None, None).unwrap(), OWN);
        assert_eq!(page(Some(3.5), None).unwrap(), (Some(252.0), None));
        assert_eq!(
            page(Some(7.0), Some(4.0)).unwrap(),
            (Some(504.0), Some(288.0))
        );
        for (w, h) in [(Some(1.0), None), (None, Some(0.0)), (Some(f64::NAN), None)] {
            let e = page(w, h).unwrap_err().to_string();
            assert!(e.contains("at least 2 by 1.5 inches"), "{e}");
        }
        assert_eq!(fill((Some(252.0), None), || (500.0, 300.0)), (252.0, 300.0));
        let not_asked = || -> (f64, f64) { panic!("both sides given") };
        assert_eq!(fill((Some(1.0), Some(2.0)), not_asked), (1.0, 2.0));
        // A figure's own floor: only a side asked for is held to it.
        assert!(at_least((Some(300.0), None), (360.0, 200.0), "a profile").is_err());
        assert!(at_least((None, Some(250.0)), (360.0, 200.0), "a profile").is_ok());
        let e = SizeArgs::new(Some(4.0), None).page_at_least((360.0, 200.0), "a profile");
        assert!(e
            .unwrap_err()
            .to_string()
            .contains("at least 5.0 by 2.8 inches"));
    }
}
