//! Figures drawn as images inside the full-screen views, in terminals that
//! speak an image protocol (kitty, iTerm2, sixel). The image is the same SVG
//! the export writes, rasterised to the area's exact pixel size; elsewhere
//! the views keep their text plots, which beat half-block images.

use std::sync::OnceLock;

use image::{DynamicImage, RgbaImage};
use ratatui::layout::{Rect, Size};
use ratatui::Frame;
use ratatui_image::picker::{Picker, ProtocolType};
use ratatui_image::protocol::Protocol;
use ratatui_image::{Image, Resize};

use data_beans::interactive::ui::{run_screen, Screen};

use super::options;

/// Figure units per terminal cell height, so 9-unit text reads at about
/// three quarters of a text line.
const UNITS_PER_CELL: f32 = 12.0;

/// Whether the environment names a terminal with an image protocol. The
/// query is only sent then: a terminal that never answers it would lose the
/// next keystroke to the query's reader.
fn likely_images(var: impl Fn(&str) -> Option<String>) -> bool {
    let has = |k: &str| var(k).is_some_and(|v| !v.is_empty());
    let term = var("TERM").unwrap_or_default();
    let program = var("TERM_PROGRAM").unwrap_or_default();
    has("KITTY_WINDOW_ID")
        || has("WEZTERM_EXECUTABLE")
        || matches!(
            program.as_str(),
            "iTerm.app" | "WezTerm" | "ghostty" | "vscode"
        )
        || var("LC_TERMINAL").is_some_and(|v| v == "iTerm2")
        || ["kitty", "ghostty", "wezterm", "foot"]
            .iter()
            .any(|t| term.contains(t))
}

/// The terminal's image support, queried once. `None` without an image
/// protocol, in a terminal not known to have one (`FABA_TUI_IMAGES=1`
/// queries anyway), or when `FABA_TUI_IMAGES=0`. Call it before a view
/// opens: the query talks to the terminal.
pub fn picker() -> Option<Picker> {
    static PICKER: OnceLock<Option<Picker>> = OnceLock::new();
    PICKER
        .get_or_init(|| {
            let var = |k: &str| std::env::var(k).ok();
            match var("FABA_TUI_IMAGES").as_deref() {
                Some("0") => return None,
                Some("1") => {}
                _ if !likely_images(var) => return None,
                _ => {}
            }
            let p = Picker::from_query_stdio().ok()?;
            (p.protocol_type() != ProtocolType::Halfblocks).then_some(p)
        })
        .clone()
}

/// `tree` rendered at `scale` pixels per unit, on white.
pub(super) fn pixmap(tree: &usvg::Tree, scale: f32) -> anyhow::Result<resvg::tiny_skia::Pixmap> {
    let size = tree
        .size()
        .to_int_size()
        .scale_by(scale)
        .ok_or_else(|| anyhow::anyhow!("figure has no size"))?;
    let mut pixmap = resvg::tiny_skia::Pixmap::new(size.width(), size.height())
        .ok_or_else(|| anyhow::anyhow!("figure too large to rasterise"))?;
    pixmap.fill(resvg::tiny_skia::Color::WHITE);
    resvg::render(
        tree,
        resvg::tiny_skia::Transform::from_scale(scale, scale),
        &mut pixmap.as_mut(),
    );
    Ok(pixmap)
}

/// Run a full-screen view with logging paused: a log line printed while
/// the view owns the terminal scrolls it and tears the picture.
pub fn run(screen: &mut impl Screen) -> anyhow::Result<()> {
    let level = log::max_level();
    log::set_max_level(log::LevelFilter::Off);
    let result = run_screen(screen);
    log::set_max_level(level);
    result
}

/// Rasterise `svg` at `scale` pixels per unit, on white.
fn raster(svg: &str, scale: f32) -> anyhow::Result<DynamicImage> {
    let pixmap = pixmap(&usvg::Tree::from_str(svg, &options())?, scale)?;
    let (w, h) = (pixmap.width(), pixmap.height());
    let rgba = RgbaImage::from_raw(w, h, pixmap.take())
        .ok_or_else(|| anyhow::anyhow!("pixel buffer size mismatch"))?;
    Ok(DynamicImage::ImageRgba8(rgba))
}

/// Run a view when stdin and stdout are a terminal; otherwise say why not
/// and carry on without it.
pub fn when_terminal(view: impl FnOnce() -> anyhow::Result<()>) -> anyhow::Result<()> {
    if data_beans::interactive::tui_available() {
        view()
    } else {
        log::warn!("--interactive needs stdin and stdout on a terminal; skipping the view");
        Ok(())
    }
}

/// One plot area drawn as an image: rebuilt only when marked stale or when
/// the area changes size.
#[derive(Default)]
pub struct PlotImage {
    stale: bool,
    size: Option<(u16, u16)>,
    protocol: Option<Protocol>,
}

impl PlotImage {
    /// The view changed: rebuild on the next frame.
    pub fn invalidate(&mut self) {
        self.stale = true;
    }

    /// Draw into `area`. `draw(width, height)` returns the SVG, in figure
    /// units, for an area that size. False when no image could be made, so
    /// the caller can fall back to text.
    pub fn render(
        &mut self,
        frame: &mut Frame,
        area: Rect,
        picker: &Picker,
        draw: impl FnOnce(f64, f64) -> String,
    ) -> bool {
        if area.width == 0 || area.height == 0 {
            return true;
        }
        let size = (area.width, area.height);
        if self.stale || self.size != Some(size) || self.protocol.is_none() {
            let font = picker.font_size();
            let scale = font.height.max(1) as f32 / UNITS_PER_CELL;
            let (px, py) = (
                area.width as f32 * font.width as f32,
                area.height as f32 * font.height as f32,
            );
            let svg = draw((px / scale) as f64, (py / scale) as f64);
            self.protocol = raster(&svg, scale).ok().and_then(|img| {
                picker
                    .new_protocol(img, Size::new(area.width, area.height), Resize::Fit(None))
                    .ok()
            });
            self.size = Some(size);
            self.stale = false;
        }
        match &self.protocol {
            Some(p) => {
                frame.render_widget(Image::new(p), area);
                true
            }
            None => false,
        }
    }
}

#[cfg(test)]
#[path = "tests/term.rs"]
mod tests;
