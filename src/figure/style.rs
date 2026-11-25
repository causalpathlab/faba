//! How a figure's text and colours look on the page, as every tool that
//! writes a figure for a paper takes it: `--font-size` (points, as
//! printed), `--panel-only` (the plot, its axes and their titles, nothing
//! else) and `--palette`.

use super::size::{self, Page, SizeArgs};
use clap::{Args, ValueEnum};

/// The colours kept and removed sites are drawn in.
#[derive(ValueEnum, Debug, Clone, Copy, Default, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Palette {
    /// Kept in blue, removed in orange.
    #[default]
    Colour,
    /// Kept in dark gray, removed in light gray.
    Gray,
}

/// What each colour of a [`Palette`] stands for, and its name in a key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Roles {
    /// Kept and removed sites at full depth, in 2D plots (deeper as a
    /// square holds more) and in histograms alike.
    pub kept: &'static str,
    pub removed: &'static str,
    /// A histogram's bars, kept; what a threshold drops; the threshold's
    /// line.
    pub bar: &'static str,
    pub dropped: &'static str,
    pub line: &'static str,
    /// Barcodes under the near-empty floor, in the cell QC.
    pub near_empty: &'static str,
    /// The kept and removed colours' names, for a key; what a threshold
    /// drops is named as removed.
    pub kept_name: &'static str,
    pub removed_name: &'static str,
    /// The kept bars' name, and the palest bars' (every site).
    pub bar_name: &'static str,
    pub faint_name: &'static str,
}

/// [`Palette::Colour`]'s roles: kept in blue, removed in orange; the bars
/// dark, what a threshold drops in the accent of the screen.
const COLOUR: Roles = Roles {
    kept: crate::qc::density::KEPT,
    removed: crate::qc::density::REMOVED,
    bar: super::BAR,
    dropped: super::ACCENT,
    line: super::ACCENT,
    near_empty: "#8a8f8b",
    kept_name: "blue",
    removed_name: "orange",
    bar_name: "dark gray",
    faint_name: "light gray",
};

/// [`Palette::Gray`]'s roles: kept in dark gray, removed in light gray, the
/// threshold's line in ink.
const GRAY: Roles = Roles {
    kept: "#262626",
    removed: "#b3b3b3",
    bar: "#262626",
    dropped: "#b3b3b3",
    line: super::INK,
    near_empty: "#8c8c8c",
    kept_name: "dark gray",
    removed_name: "light gray",
    bar_name: "dark gray",
    faint_name: "palest",
};

impl Palette {
    /// What each of its colours stands for.
    pub fn roles(self) -> &'static Roles {
        match self {
            Palette::Colour => &COLOUR,
            Palette::Gray => &GRAY,
        }
    }
}

/// `--font-size`, `--panel-only` and `--palette`.
#[derive(Args, Debug, Clone, Copy, Default, PartialEq, serde::Serialize)]
pub struct StyleArgs {
    #[arg(
        long = "font-size",
        value_name = "PT",
        help = "Every label's size in points, as printed (default: fit to the figure)"
    )]
    pub font_size: Option<f64>,

    #[arg(
        long = "panel-only",
        default_value_t = false,
        help = "Draw the plot, its axes and their titles only: no title, file name, key or notes"
    )]
    pub panel_only: bool,

    #[arg(
        long = "palette",
        value_enum,
        default_value_t = Palette::Colour,
        help = "Kept and removed sites in blue and orange (colour) or dark and light gray (gray)"
    )]
    pub palette: Palette,
}

/// The text size a figure takes when none is asked for: 8 points on a
/// figure of its own size, else fit to the page asked for, from 5 points on
/// a small panel to 8 on a large one.
pub const OWN_FONT: f64 = 8.0;

/// A figure's look, worked out for its page.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Style {
    /// Labels' size in points: what the figures' 8-point labels are drawn
    /// at, the others in proportion.
    pub font: f64,
    /// Every label at `font`, as `--font-size` asks, not in proportion.
    pub uniform: bool,
    pub panel_only: bool,
    pub palette: Palette,
}

impl Default for Style {
    fn default() -> Self {
        Self {
            font: OWN_FONT,
            uniform: false,
            panel_only: false,
            palette: Palette::Colour,
        }
    }
}

impl Style {
    /// What the labels and the room around them are scaled by, from the
    /// size the figures were first laid out for ([`OWN_FONT`]).
    pub fn em(&self) -> f64 {
        self.font / OWN_FONT
    }

    /// This look for a figure that fits its text to its own boxes: the size
    /// asked for kept, else its labels left as the boxes fit them.
    pub fn for_boxes(self) -> Self {
        match self.uniform {
            true => self,
            false => Self {
                font: OWN_FONT,
                ..self
            },
        }
    }
}

impl StyleArgs {
    /// The size asked for, checked: a positive number of points.
    pub fn check(&self) -> anyhow::Result<()> {
        if let Some(pt) = self.font_size {
            anyhow::ensure!(
                pt.is_finite() && pt > 0.0,
                "--font-size {pt}: a size in points, above 0"
            );
        }
        Ok(())
    }

    /// The look on a page `w` by `h` points; `asked` whether the page was
    /// asked for (else the figure's own, at [`OWN_FONT`]).
    pub fn style(&self, (w, h): (f64, f64), asked: bool) -> Style {
        let font = match (self.font_size, asked) {
            (Some(pt), _) => pt,
            (None, true) => data_beans::figure::panel::fit_font(w, h),
            (None, false) => OWN_FONT,
        };
        Style {
            font,
            uniform: self.font_size.is_some(),
            panel_only: self.panel_only,
            palette: self.palette,
        }
    }
}

/// A figure's page and look, as asked for: what it lays itself out for.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Sheet {
    pub page: Page,
    pub look: StyleArgs,
}

impl Sheet {
    /// The figure's own size, in its own look.
    pub const OWN: Sheet = Sheet {
        page: size::OWN,
        look: StyleArgs {
            font_size: None,
            panel_only: false,
            palette: Palette::Colour,
        },
    };

    /// `page`, in the figure's own look.
    #[cfg(test)]
    pub fn of(page: Page) -> Self {
        Self { page, ..Self::OWN }
    }

    /// The page `size` asks for and the look `look` does, both checked.
    pub fn new(size: &SizeArgs, look: &StyleArgs) -> anyhow::Result<Self> {
        look.check()?;
        Ok(Self {
            page: size.page()?,
            look: *look,
        })
    }

    /// The figure's size in points: each side asked for, the others from
    /// `own` (see [`size::fill`]).
    pub fn fill(&self, own: impl FnOnce() -> (f64, f64)) -> (f64, f64) {
        size::fill(self.page, own)
    }

    /// The look on the figure's size `(w, h)`, its text fit to it when a
    /// side was asked for.
    pub fn style(&self, (w, h): (f64, f64)) -> Style {
        self.look.style((w, h), self.page != size::OWN)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The font is the one asked for, else fit to a page asked for, else
    /// the figures' own 8 points; a size asked for sets every label.
    #[test]
    fn the_font_is_asked_fit_or_the_figures_own() {
        let none = StyleArgs::default();
        assert_eq!(none.style((100.0, 80.0), false), Style::default());
        let fit = none.style((130.0, 100.0), true);
        assert!((fit.font - 100.0 / 18.0).abs() < 1e-12 && !fit.uniform);
        assert_eq!(none.style((1000.0, 800.0), true).font, OWN_FONT);
        assert_eq!(none.style((50.0, 40.0), true).font, 5.0);
        let asked = StyleArgs {
            font_size: Some(7.0),
            ..none
        };
        let s = asked.style((1000.0, 800.0), true);
        assert_eq!((s.font, s.uniform, s.em()), (7.0, true, 0.875));
        assert_eq!(s.for_boxes(), s);
        assert_eq!(fit.for_boxes().font, OWN_FONT);
        for pt in [0.0, -1.0, f64::NAN] {
            let bad = StyleArgs {
                font_size: Some(pt),
                ..none
            };
            let e = Sheet::new(&SizeArgs::default(), &bad).unwrap_err();
            assert!(e.to_string().contains("--font-size"), "{e}");
        }
    }

    /// A sheet's font follows its page: its own size at 8 points, a page
    /// asked for fit.
    #[test]
    fn a_sheet_fits_its_font_only_to_a_page_asked_for() {
        let own = Sheet::default();
        assert_eq!(own.style((100.0, 80.0)).font, OWN_FONT);
        let small = Sheet::of((Some(100.0), None));
        assert_eq!(small.style((100.0, 80.0)).font, 5.0);
        assert_eq!(small.fill(|| (1.0, 80.0)), (100.0, 80.0));
    }

    /// The look goes with the figure it lays out: on metagene and pwm only
    /// with --figure, on qc only with --plots, as the size does; refused
    /// rather than dropped.
    #[test]
    fn the_look_goes_with_a_figure() {
        let parse = |line: &[&str]| {
            let line = ["faba"].iter().chain(line);
            crate::matches(&mut crate::faba_command(), line.copied()).map_err(|e| e.to_string())
        };
        for (command, input, needs) in [
            ("metagene", "out/", "--figure"),
            ("pwm", "out/", "--figure"),
            ("qc", "in/", "--plots"),
        ] {
            for look in [
                &["--font-size", "6"][..],
                &["--panel-only"],
                &["--palette", "gray"],
            ] {
                let line: Vec<&str> = [command, input].iter().chain(look).copied().collect();
                let e = parse(&line).unwrap_err();
                assert!(e.contains(needs), "{line:?}: {e}");
                let with = match needs {
                    "--figure" => ["--figure", "f"].as_slice(),
                    _ => &["--plots", "--batch-process"],
                };
                let line: Vec<&str> = line.iter().copied().chain(with.iter().copied()).collect();
                assert!(parse(&line).is_ok(), "{line:?}");
            }
            assert!(parse(&[command, input]).is_ok(), "{command}");
        }
    }
}
