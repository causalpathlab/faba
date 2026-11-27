//! The browser's raw view (`r`): the reads under the fit in the modules
//! pane, drawn as the raw figure draws them (see `--raw`), as an image,
//! only in a terminal that shows images. The reads are read on a thread
//! of their own the first time it is shown.

use super::{Look, ModulesView};
use crate::cluster::browse::figure::{self, GridArgs, RawDraw};
use crate::cluster::browse::raw::{Plan, Raw, RawArgs};
use crate::figure::style::{Sheet, StyleArgs};
use crate::figure::term::PlotImage;
use ratatui::layout::Rect;
use ratatui::Frame;
use std::path::PathBuf;
use std::thread::JoinHandle;

/// The reads under the fit: being read, or read (or why not).
enum Load {
    Reading(JoinHandle<Result<Raw, String>>),
    Read(Result<Raw, String>),
}

/// The raw view's state: where the fit is and how to draw it, whether it
/// is shown, its reads, and its image.
#[derive(Default)]
pub(in crate::cluster::browse) struct RawView {
    source: Option<(PathBuf, RawArgs)>,
    on: bool,
    load: Option<Load>,
    image: PlotImage,
}

/// What the footer says while the reads are read, and where the terminal
/// shows no images.
const READING: &str = "reading the reads under the fit…";
const NO_IMAGES: &str =
    "the raw view is drawn as an image: it needs a terminal that shows them (kitty, iTerm2, sixel)";

/// Pixels a point a terminal image is drawn at: as many as its font's
/// cells have, a cell twelve points tall (see [`PlotImage::render`]).
fn per_point(picker: &ratatui_image::picker::Picker) -> f64 {
    f64::from(picker.font_size().height.max(1)) / 12.0
}

impl<'a> ModulesView<'a> {
    /// The raw view of the fit in `dir`, drawn as `args` asks (`r`), shown
    /// from the start with `--raw`.
    pub fn with_raw(mut self, dir: PathBuf, args: RawArgs) -> Self {
        let on = args.raw;
        self.raw.source = Some((dir, args));
        if on {
            self.show_raw();
        }
        self
    }

    /// Show the raw view, its reads read on a thread of their own the
    /// first time (or again after they failed); refused for a count fit.
    fn show_raw(&mut self) {
        let Some((dir, args)) = self.raw.source.clone() else {
            self.status = Some("no fit folder to read the reads from".into());
            return;
        };
        if !matches!(self.raw.load, Some(Load::Read(Err(_))) | None) {
            self.raw.on = true;
            return;
        }
        match Plan::of(self.m, args.unit) {
            Ok(plan) => {
                let read = move || plan.read(&dir, &args.data).map_err(|e| format!("{e:#}"));
                self.raw.load = Some(Load::Reading(std::thread::spawn(read)));
                self.raw.on = true;
                self.status = Some(READING.into());
            }
            Err(e) => self.status = Some(format!("{e:#}")),
        }
    }

    /// `r`: the raw view on or off; on only where the terminal shows
    /// images.
    pub(super) fn toggle_raw(&mut self) {
        self.raw.image.invalidate();
        match self.raw.on {
            true => self.raw.on = false,
            false if self.controls.images_on() => self.show_raw(),
            false => self.status = Some(NO_IMAGES.into()),
        }
    }

    /// The reads, once read.
    fn raw_read(&self) -> Option<&Raw> {
        match &self.raw.load {
            Some(Load::Read(Ok(raw))) => Some(raw),
            _ => None,
        }
    }

    /// The reads and how to draw them, while the raw view is drawn: what a
    /// saved screen draws in the grid's place.
    pub(super) fn raw_shown(&self) -> Option<(&Raw, &RawArgs)> {
        let raw = self.raw_read().filter(|_| self.raw_drawn())?;
        Some((raw, &self.raw.source.as_ref()?.1))
    }

    /// Whether the modules pane draws the raw view this frame.
    pub(super) fn raw_drawn(&self) -> bool {
        self.raw.on && self.controls.images_on()
    }

    /// Take the reads once read, saying why not when they failed (and
    /// leaving the raw view); true when that changes the screen.
    pub(super) fn take_raw(&mut self) -> bool {
        let Some(Load::Reading(h)) = &self.raw.load else {
            return false;
        };
        if !h.is_finished() {
            return false;
        }
        let Some(Load::Reading(h)) = self.raw.load.take() else {
            unreachable!("the load was being read");
        };
        let read = h.join().unwrap_or_else(|_| Err("the read failed".into()));
        if let Err(e) = &read {
            self.raw.on = false;
            self.status = Some(format!("raw view: {e}"));
        } else if self.status.as_deref() == Some(READING) {
            self.status = None;
        }
        self.raw.load = Some(Load::Read(read));
        self.raw.image.invalidate();
        true
    }

    /// The raw view in `area`, as an image; true when drawn (or still
    /// being read, the pane left blank), false where the grid is drawn
    /// instead: no images (said once), or no image made.
    pub(super) fn draw_raw(&mut self, frame: &mut Frame, area: Rect) -> bool {
        if !self.raw.on {
            return false;
        }
        let Some(picker) = self.controls.images().cloned() else {
            self.raw.on = false;
            self.status = Some(NO_IMAGES.into());
            return false;
        };
        // The rows as the grid's: by shift range under the bars, else
        // diagonal.
        let look = match self.looks[0] {
            Look::Bars => figure::Look::Bars,
            _ => figure::Look::Shift,
        };
        let grid = GridArgs {
            key: true,
            ..GridArgs::default()
        };
        let m = self.m;
        let RawView {
            source,
            load,
            image,
            ..
        } = &mut self.raw;
        let raw = match load {
            Some(Load::Read(Ok(raw))) => &*raw,
            Some(Load::Read(Err(_))) => return false,
            _ => return true,
        };
        let Some((_, args)) = source.as_ref() else {
            return false;
        };
        let how = RawDraw {
            per_point: per_point(&picker),
            ..RawDraw::of(raw, args)
        };
        image.render(frame, area, &picker, |w, h| {
            let sheet = Sheet {
                page: (Some(w), Some(h)),
                look: StyleArgs {
                    font_size: Some(9.0),
                    panel_only: true,
                    ..StyleArgs::default()
                },
            };
            figure::draw_raw(m, how, look, &grid, "", sheet).unwrap_or_else(|e| {
                crate::figure::svg(w, h, |c| {
                    c.text(
                        9.0,
                        18.0,
                        &e,
                        9.0,
                        crate::figure::Anchor::Start,
                        crate::figure::MUTED,
                    )
                })
            })
        })
    }
}
