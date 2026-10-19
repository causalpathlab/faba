//! Figures drawn as images inside the full-screen views, in terminals that
//! speak an image protocol (kitty, iTerm2, sixel). The image is the same SVG
//! the export writes, rasterised to the area's exact pixel size; elsewhere
//! the views keep their text plots, which beat half-block images.

use image::{DynamicImage, RgbaImage};
use legume_tui::graphics::Fallback;
use ratatui::layout::{Rect, Size};
use ratatui::Frame;
use ratatui_image::picker::{Picker, ProtocolType};
use ratatui_image::protocol::Protocol;
use ratatui_image::{Image, Resize};

use super::options;

/// Figure units per terminal cell height, so 9-unit text reads at about
/// three quarters of a text line.
const UNITS_PER_CELL: f32 = 12.0;

/// Whether the environment names a terminal with an image protocol. The
/// query is only sent then, so a terminal without one never waits for its
/// answer.
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

/// Whether to ask the terminal for its image protocol, by
/// `FABA_TUI_IMAGES` (`images`) and the terminal's environment (`var`):
/// never with `0`, always with `1`, else where the environment names a
/// terminal that has one.
fn asks(images: Option<&str>, var: impl Fn(&str) -> Option<String>) -> bool {
    match images {
        Some("0") => false,
        Some("1") => true,
        _ => likely_images(var),
    }
}

/// The image support of the terminal `session` is attached to now,
/// asked of it (`legume_tui::graphics`, with a timeout). `None` without
/// an image protocol (block characters, where the text plots do better),
/// in a terminal not known to have one (`FABA_TUI_IMAGES=1` asks anyway),
/// or when `FABA_TUI_IMAGES=0`. Call it with the screen up: the query
/// talks to the terminal, and the answer is the attached terminal's,
/// asked again after every attach (see [`Controls::refresh`]).
///
/// [`Controls::refresh`]: super::Controls::refresh
pub fn picker(session: &legume_tui::Session) -> Option<Picker> {
    let images = std::env::var("FABA_TUI_IMAGES").ok();
    if !asks(images.as_deref(), |k| session.terminal_env(k)) {
        return None;
    }
    let (p, _) = legume_tui::graphics::pick_protocol(session, None, Fallback::Probe);
    (p.protocol_type() != ProtocolType::Halfblocks).then_some(p)
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

/// Rasterise `svg` at `scale` pixels per unit, on white.
fn raster(svg: &str, scale: f32) -> anyhow::Result<DynamicImage> {
    let pixmap = pixmap(&usvg::Tree::from_str(svg, &options())?, scale)?;
    let (w, h) = (pixmap.width(), pixmap.height());
    let rgba = RgbaImage::from_raw(w, h, pixmap.take())
        .ok_or_else(|| anyhow::anyhow!("pixel buffer size mismatch"))?;
    Ok(DynamicImage::ImageRgba8(rgba))
}

/// `img` as a terminal image fitting `area`.
pub(super) fn protocol(picker: &Picker, img: DynamicImage, area: Rect) -> Option<Protocol> {
    picker
        .new_protocol(img, Size::new(area.width, area.height), Resize::Fit(None))
        .ok()
}

/// One plot area drawn as an image: rebuilt only when marked stale, when
/// the area changes size, or after a terminal attached (its picker may be
/// another, and kitty's images were deleted by the attach).
#[derive(Default)]
pub struct PlotImage {
    stale: bool,
    size: Option<(u16, u16)>,
    protocol: Option<Protocol>,
    /// The attach ([`legume_tui::attaches`]) the image was made at, and
    /// the picker's [`kind`].
    made: Option<(u64, PickerKind)>,
}

/// What an image made with a picker depends on: its protocol and cell size.
type PickerKind = (ProtocolType, u16, u16);

/// `picker`'s [`PickerKind`].
fn kind(picker: &Picker) -> PickerKind {
    let font = picker.font_size();
    (picker.protocol_type(), font.width, font.height)
}

impl PlotImage {
    /// The view changed: rebuild on the next frame.
    pub fn invalidate(&mut self) {
        self.stale = true;
    }

    /// Whether the next frame rebuilds it.
    #[cfg(test)]
    pub fn is_stale(&self) -> bool {
        self.stale
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
        self.render_at(frame, area, picker, legume_tui::attaches(), draw)
    }

    /// [`PlotImage::render`] after `attach` attaches: an image made at
    /// another, or with a picker of another protocol or cell size, is
    /// made again.
    fn render_at(
        &mut self,
        frame: &mut Frame,
        area: Rect,
        picker: &Picker,
        attach: u64,
        draw: impl FnOnce(f64, f64) -> String,
    ) -> bool {
        if area.width == 0 || area.height == 0 {
            return true;
        }
        let size = (area.width, area.height);
        let made = (attach, kind(picker));
        if self.made != Some(made) {
            self.invalidate();
            self.made = Some(made);
        }
        if self.stale || self.size != Some(size) || self.protocol.is_none() {
            let font = picker.font_size();
            let scale = font.height.max(1) as f32 / UNITS_PER_CELL;
            let (px, py) = (
                area.width as f32 * font.width as f32,
                area.height as f32 * font.height as f32,
            );
            let svg = draw((px / scale) as f64, (py / scale) as f64);
            self.protocol = raster(&svg, scale)
                .ok()
                .and_then(|img| protocol(picker, img, area));
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
