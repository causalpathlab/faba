//! A fit's figures, as `--figure` writes them (a PDF and a PNG) and MCP's
//! `cluster_view` returns them: the modules grid (gene modules against cell
//! modules, each block coloured by its fitted shift, a count fit's log rate,
//! its methylated %, UMIs per cell, or that standardized across the cell
//! modules, drawn as boxes or as a Hinton diagram's squares), the trace (each restart's ELBO by
//! step), and the cell modules against the batches. Each is laid out for a
//! size when one is given (inches), its text at a size that reads in print,
//! and drawn in the panel look asked for (see [`crate::figure::style`]):
//! the text at `--font-size`, the plot alone with `--panel-only` (no title,
//! key or notes), in grays with `--palette gray`. `--xlab` and `--ylab`
//! title the trace's axes; the modules grid and the batches have no axis
//! titles and leave them unused (see [`unused_titles`]).

use super::order::standardized;
use super::Modules;

mod genes;
mod grid;
pub mod screen;
mod sites;
use crate::figure::size::{self, Page};
use crate::figure::style::{AxisTitles, Palette, Sheet};
use crate::figure::{esc, mix, svg, Anchor, Canvas, FAINT, INK, MUTED};
#[cfg(test)]
pub(super) use grid::go_width;

/// Which figure.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Plot {
    /// Gene modules against cell modules.
    Modules,
    /// Each restart's ELBO by step.
    Trace,
    /// Cell modules against the batches: each one's share in each.
    Batches,
}

/// A figure's text size in points: 9 at its own size, less on a small
/// page (to 5, the least that reads in print), by the sides asked for.
fn font(page: Page) -> f64 {
    let side = [page.0, page.1]
        .into_iter()
        .flatten()
        .fold(f64::INFINITY, f64::min);
    match side.is_finite() {
        true => (side / 40.0).clamp(5.0, 9.0),
        false => 9.0,
    }
}

/// The least cell side a grid is drawn with, in points: its values need
/// the room.
const CELL_LEAST: f64 = 12.0;

/// Why a figure is not drawn: its page too small for its cells (the
/// least it needs, in points, on the sides asked for, at the font it was
/// laid out in), or another reason.
#[derive(Debug)]
enum Refused {
    Small(Page),
    Said(String),
}

impl From<String> for Refused {
    fn from(said: String) -> Self {
        Refused::Said(said)
    }
}

/// Square cells, `cols` by `rows` (in cells, a margin's share too),
/// beside `fixed` points across and down, `own` points each at the grid's
/// own size: a side not asked for follows the grid's own proportions, and
/// the cells are as large as the page allows. The cell's side and the
/// figure's size; refused when the cells would be too small to read
/// (under [`CELL_LEAST`]), with the least the sides asked for need.
fn cells(
    page: Page,
    own: f64,
    (fw, fh): (f64, f64),
    (cols, rows): (f64, f64),
) -> Result<(f64, (f64, f64)), Refused> {
    let (cols, rows) = (cols.max(1.0), rows.max(1.0));
    let (ow, oh) = (fw + own * cols, fh + own * rows);
    let (w, h) = match page {
        (None, None) => return Ok((own, (ow, oh))),
        (Some(w), None) => (w, w * oh / ow),
        (None, Some(h)) => (h * ow / oh, h),
        (Some(w), Some(h)) => (w, h),
    };
    let cell = ((w - fw) / cols).min((h - fh) / rows);
    // A size read back from its inches is not refused for a rounding.
    if cell < CELL_LEAST - 1e-9 {
        let (nw, nh) = (fw + CELL_LEAST * cols, fh + CELL_LEAST * rows);
        return Err(Refused::Small(match page {
            (Some(_), None) => (Some(nw.max(nh * ow / oh)), None),
            (None, Some(_)) => (None, Some(nh.max(nw * oh / ow))),
            _ => (Some(nw), Some(nh)),
        }));
    }
    Ok((cell, (w, h)))
}

/// A figure drawn on `sheet` by `one`; refused as too small, saying the
/// least page it is drawn on (see [`least`]).
fn drawn(sheet: Sheet, one: impl Fn(Sheet) -> Result<String, Refused>) -> Result<String, String> {
    match one(sheet) {
        Ok(drawn) => Ok(drawn),
        Err(Refused::Said(why)) => Err(why),
        Err(Refused::Small(need)) => Err(least(sheet, need, one)),
    }
}

/// The least page, on the sides `sheet` asks for, that `one` draws, from
/// the `need` it said: each side grown to what it needs, in twentieths
/// of an inch, and tried again (a larger page may set larger text, which
/// needs more), so the size said is one that is drawn.
fn least(sheet: Sheet, mut need: Page, one: impl Fn(Sheet) -> Result<String, Refused>) -> String {
    let step = size::PT / 20.0;
    let up = |v: f64| (v / step - 1e-6).ceil() * step;
    let mut page = sheet.page;
    for _ in 0..400 {
        let grown =
            |side: Option<f64>, need: Option<f64>| side.map(|s| up(s.max(need.unwrap_or(s))));
        let next = (grown(page.0, need.0), grown(page.1, need.1));
        // Nothing grew (a need rounded to the side): a step more each.
        page = match next == page {
            true => (page.0.map(|s| up(s + step)), page.1.map(|s| up(s + step))),
            false => next,
        };
        match one(Sheet { page, ..sheet }) {
            Err(Refused::Small(more)) => need = more,
            _ => break,
        }
    }
    let inches = |p: f64| format!("{:.2}", p / size::PT);
    let least = match page {
        (Some(w), Some(h)) => format!("{} by {} inches", inches(w), inches(h)),
        (Some(w), None) => format!("{} inches across", inches(w)),
        (_, h) => format!("{} inches tall", inches(h.unwrap_or_default())),
    };
    format!("the grid needs at least {least} for its labels and cells")
}

/// The width of `s` in `f`-point text, about.
fn width(s: &str, f: f64) -> f64 {
    s.chars().count() as f64 * f * 0.56
}

/// The figure's title at its top left, in `f + 1`-point text, shortened
/// from the front to fit a figure `w` points wide.
fn heading(c: &mut Canvas, title: &str, f: f64, w: f64) {
    let shown = title_text(title, f, w - f);
    c.text(f / 2.0, 1.6 * f, &esc(&shown), f + 1.0, Anchor::Start, INK);
}

/// The title in `f + 1`-point text, at most `room` points wide: shortened
/// from the front when it would not fit.
fn title_text(title: &str, f: f64, room: f64) -> String {
    let fits = (room / width("x", f + 1.0)) as usize;
    let n = title.chars().count();
    match n > fits && fits > 1 {
        true => std::iter::once('…')
            .chain(title.chars().skip(n + 1 - fits))
            .collect(),
        false => title.to_string(),
    }
}

/// A side of the modules grid a margin's bars go on: the gene modules'
/// on the left, by the rows; the cell modules' on top or at the bottom,
/// by the columns. The right is the GO terms'.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Side {
    Left,
    Top,
    Bottom,
}

/// What the modules grid draws besides its blocks (only with `--figure`
/// and `--plot modules`): the modules' sizes and reads on its margins,
/// GO terms on its right, the key under `--panel-only`, the blocks'
/// values.
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct GridArgs {
    #[arg(
        long = "size-bars",
        value_enum,
        value_name = "SIDES",
        value_delimiter = ',',
        num_args = 0..=1,
        default_missing_value = "left,top",
        requires = "drawn",
        help = "The modules grid's margins: each gene module's size in genes as a bar by its row (left), each cell module's in cells by its column (top or bottom); bare: left,top"
    )]
    pub size_bars: Vec<Side>,

    #[arg(
        long = "read-bars",
        value_enum,
        value_name = "SIDES",
        value_delimiter = ',',
        num_args = 0..=1,
        default_missing_value = "left,top",
        requires = "drawn",
        help = "The modules grid's margins: each module's reads summed over its blocks as the browser's mirror bars (converted out or up, unconverted the other way, one scale a side), sides as --size-bars"
    )]
    pub read_bars: Vec<Side>,

    #[arg(
        long = "marginal-size",
        value_name = "FRACTION",
        default_value_t = MARGINAL,
        value_parser = fraction,
        requires = "drawn",
        help = "The depth of the grid's margin bars, a fraction of its shorter side, in (0, 1]"
    )]
    pub marginal_size: f64,

    #[arg(
        long = "key",
        requires = "drawn",
        help = "Keep the grid's key under --panel-only, small, under the grid (for --look hinton, squares at reference |z| too)"
    )]
    pub key: bool,

    #[arg(
        long = "values",
        requires = "drawn",
        help = "Print each block's value in its cell at the font size (z for --look z and hinton, the shift, the %), left out where it does not fit"
    )]
    pub values: bool,

    #[arg(
        long = "go-labels",
        value_name = "N",
        value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..=10),
        requires = "drawn",
        help = "Label each gene module's row, right of the grid, with up to N of its GO terms (or --gmt sets), a line each with its q after it, by p and none mostly within a term already shown (--gaf, --obo, --gmt, --go-q, --go-min-weight); a row with none is left blank"
    )]
    pub go_labels: Option<usize>,

    #[arg(
        long = "go-none",
        requires = "go_shown",
        help = "With --go-labels (or --screen's --side go), say \"none at q <= Q\" on a row with no term, as the browser's GO column does, instead of leaving it blank"
    )]
    pub go_none: bool,

    #[arg(
        long = "block-bars",
        requires = "drawn",
        help = "Draw each block's reads as --look bars draws them (converted up, unconverted down, one scale) over the look's own block, lightened"
    )]
    pub block_bars: bool,

    /// Each row's shift range right of the grid under every look, as the
    /// browser writes it (the bars look has it anyway): the screen's.
    #[arg(skip)]
    pub ranges: bool,
}

/// The margin bars' depth when `--marginal-size` is not given.
const MARGINAL: f64 = 0.15;

impl Default for GridArgs {
    fn default() -> Self {
        Self {
            size_bars: Vec::new(),
            read_bars: Vec::new(),
            marginal_size: MARGINAL,
            key: false,
            values: false,
            go_labels: None,
            go_none: false,
            block_bars: false,
            ranges: false,
        }
    }
}

impl GridArgs {
    /// Refused: the cell modules' bars both on top and at the bottom.
    pub fn check(&self) -> anyhow::Result<()> {
        for (flag, sides) in [
            ("--size-bars", &self.size_bars),
            ("--read-bars", &self.read_bars),
        ] {
            anyhow::ensure!(
                !(sides.contains(&Side::Top) && sides.contains(&Side::Bottom)),
                "{flag}: the cell modules' bars go top or bottom, not both"
            );
        }
        Ok(())
    }
}

/// A fraction in (0, 1].
fn fraction(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("{s} is not a number"))?;
    match v > 0.0 && v <= 1.0 {
        true => Ok(v),
        false => Err(format!("{s} is not in (0, 1]")),
    }
}

/// The size bars' fill, in either palette.
pub const SIZE_BAR: &str = "#8c8c8c";
/// The read bars' fills: converted reads in the accent, unconverted in the
/// plain tone, as the browser's mirror bars are.
pub const READ_UP: &str = crate::figure::ACCENT;
pub const READ_OFF: &str = MUTED;
/// The same in grays: dark, and light.
const READ_UP_GRAY: &str = "#262626";
const READ_OFF_GRAY: &str = "#b3b3b3";
/// The unconverted side of the browser's mirror bars on paper (the pinned
/// genes' tracks, the sites'): a light gray in either palette, as the
/// terminal's plain light tone; and the least a side with reads is drawn,
/// an eighth of a terminal line (two columns deep), in font sizes, as
/// the browser's least block.
pub const MIRROR_OFF: &str = READ_OFF_GRAY;
pub const MIRROR_LEAST: f64 = 1.2 / 8.0;
/// The blocks' bars: converted reads in the accent, unconverted in a light
/// tone outlined so it reads on white paper; in grays, dark and light.
pub const BLOCK_UP: &str = crate::figure::ACCENT;
pub const BLOCK_OFF: &str = "#f2efe9";
pub const BLOCK_UP_GRAY: &str = READ_UP_GRAY;
pub const BLOCK_OFF_GRAY: &str = "#d9d9d9";
/// The light bars' outline, and its width in points.
const BLOCK_EDGE: (&str, f64) = (MUTED, 0.5);

/// The GO terms `--go-labels` puts on the grid, read as `ask` says; none
/// read without it.
pub fn go_table(
    m: &Modules,
    go_labels: Option<usize>,
    ask: &super::go::Ask,
) -> anyhow::Result<Option<super::go::Table>> {
    if go_labels.is_none() {
        return Ok(None);
    }
    let source = ask.source(m).map_err(anyhow::Error::msg)?;
    super::go::table(m, &source, ask.min_weight, ask.q).map(Some)
}

/// The figure `plot` of `m` under `look`, titled `title`, laid out for
/// the sheet's page (see [`crate::figure::SizeArgs`]) in its look, the
/// trace's axes titled as `titles` asks; why there is none (no trace, one
/// batch).
pub fn draw(
    m: &Modules,
    plot: Plot,
    look: Look,
    title: &str,
    sheet: Sheet,
    titles: &AxisTitles,
) -> Result<String, String> {
    match plot {
        Plot::Modules => draw_with(m, look, &GridArgs::default(), None, title, sheet),
        Plot::Trace => trace(m, title, sheet.page, Ink::of(&sheet), titles),
        Plot::Batches => drawn(sheet, |sheet| {
            batches(m, title, sheet.page, Ink::of(&sheet))
        }),
    }
}

/// The note `--xlab` or `--ylab` given to a figure with no axis titles
/// (the modules grid, the batches) is told: they are unused there.
pub fn unused_titles(plot: Plot, titles: &AxisTitles) -> Option<String> {
    let given = titles.xlab.is_some() || titles.ylab.is_some();
    let what = match plot {
        Plot::Trace => return None,
        Plot::Modules => "the modules grid",
        Plot::Batches => "the batches figure",
    };
    given.then(|| {
        format!("--xlab and --ylab title the trace's axes; {what} has none, so they are unused")
    })
}

/// The modules grid of `m` under `look` with `grid`'s extras, `go` the
/// terms for `--go-labels`.
pub fn draw_with(
    m: &Modules,
    look: Look,
    grid: &GridArgs,
    go: Option<&super::go::Table>,
    title: &str,
    sheet: Sheet,
) -> Result<String, String> {
    drawn(sheet, |sheet| {
        grid::heatmap(m, look, grid, go, title, sheet.page, Ink::of(&sheet))
    })
}

/// How a figure is drawn, from its sheet's look: its text size in points,
/// whether the plot is drawn alone, and whether in grays.
#[derive(Debug, Clone, Copy)]
struct Ink {
    f: f64,
    panel_only: bool,
    gray: bool,
}

/// The gray palette's ends: a value of 0 and the top value.
const GRAY_LOW: (f64, f64, f64) = (242.0, 242.0, 242.0);
const GRAY_HIGH: (f64, f64, f64) = (38.0, 38.0, 38.0);

impl Ink {
    /// The look `sheet` asks for: the text at its size, else fit to the
    /// page (see [`font`]).
    fn of(sheet: &Sheet) -> Self {
        Self {
            f: sheet.look.font_size.unwrap_or_else(|| font(sheet.page)),
            panel_only: sheet.look.panel_only,
            gray: sheet.look.palette == Palette::Gray,
        }
    }

    /// A value's fill, as [`fill`] gives it; in grays, darker as it is
    /// further from 0 either way.
    fn fill(&self, v: f64, top: f64, diverging: bool) -> (String, bool) {
        match self.gray {
            true => {
                let t = (v.abs() / top.max(1e-12)).min(1.0);
                (mix(GRAY_LOW, GRAY_HIGH, t), t > 0.6)
            }
            false => fill(v, top, diverging),
        }
    }

    /// The title, but for the plot alone.
    fn heading(&self, c: &mut Canvas, title: &str, w: f64) {
        if !self.panel_only {
            heading(c, title, self.f, w);
        }
    }

    /// `room` for a title, a key or a note: none for the plot alone.
    fn room(&self, room: f64) -> f64 {
        if self.panel_only {
            0.0
        } else {
            room
        }
    }
}

impl Plot {
    /// The figure it is, by its floor.
    pub fn figure(self) -> size::Figure {
        match self {
            Plot::Modules => size::Figure::Modules,
            Plot::Trace => size::Figure::Trace,
            Plot::Batches => size::Figure::Batches,
        }
    }
}

/// The page `asked` for `plot`: refused when a side is under its floor
/// (see [`size::Figure::least`]), the trace's and the batches' margins not
/// shrinking with their text; the grid lays itself out, saying the least
/// it needs as it is drawn.
pub fn page(plot: Plot, asked: &crate::figure::SizeArgs) -> anyhow::Result<Page> {
    plot.figure().page(asked)
}

/// What colours a block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Look {
    /// The block's fitted shift (a count fit's log rate), on a diverging
    /// scale.
    Shift,
    /// Its methylated % (a count fit's UMIs per cell), light to dark.
    Percent,
    /// That standardized across the cell modules, within each gene module.
    Z,
    /// That z as a Hinton diagram: a square per block, its area by |z|
    /// against the largest, coloured as the browser's Hinton look (red
    /// above the row's mean, blue below, deeper as |z| grows).
    Hinton,
    /// Its reads as the browser's mirror bars: converted rising above a
    /// zero line, unconverted hanging below, every block on one scale;
    /// the rows by shift range, as the browser orders them.
    Bars,
}

/// The diverging pair and its neutral midpoint, and the sequential ramp's
/// ends.
const LOW: (f64, f64, f64) = (42.0, 120.0, 214.0);
const MID: (f64, f64, f64) = (240.0, 239.0, 236.0);
const HIGH: (f64, f64, f64) = (227.0, 73.0, 72.0);
const DARK: (f64, f64, f64) = (26.0, 72.0, 140.0);

/// A value's fill: diverging around 0 to `top` either way, or sequential
/// from 0 to `top`; and whether it is dark enough for light ink.
fn fill(v: f64, top: f64, diverging: bool) -> (String, bool) {
    let t = (v.abs() / top.max(1e-12)).min(1.0);
    let far = match (diverging, v >= 0.0) {
        (false, _) => DARK,
        (true, true) => HIGH,
        (true, false) => LOW,
    };
    (mix(MID, far, t), t > 0.6)
}

/// A Hinton square's fill as the browser's Hinton look colours it (see
/// [`super::view::Stretch::hinton`]), `z_top` the |z| at full red or
/// blue; and whether it is dark enough for light ink.
fn hinton_fill(v: f64, z_top: f64) -> (String, bool) {
    let stretch = super::view::Stretch {
        lo: 0.0,
        hi: 1.0,
        z_top,
        deepest: 0.0,
    };
    let (colour, t) = stretch.hinton(v);
    (hex(colour), t > 0.5)
}

/// A terminal colour as SVG writes it.
fn hex(colour: ratatui::style::Color) -> String {
    match colour {
        ratatui::style::Color::Rgb(r, g, b) => format!("#{r:02x}{g:02x}{b:02x}"),
        _ => INK.to_string(),
    }
}

/// Each gene module's value in each cell module under `look`, rows the
/// gene modules with sites in the diagonal order (by shift range for the
/// bars, as the browser orders them), columns as the panes show them;
/// none for a block without reads.
fn values(m: &Modules, look: Look) -> (Vec<usize>, Vec<Vec<Option<f64>>>) {
    let rows = match look {
        Look::Bars => m.order.clone(),
        _ => m.diagonal.clone(),
    };
    let values = rows
        .iter()
        .map(|&g| match look {
            Look::Shift => (m.columns.iter()).map(|&c| Some(m.shift(c, g))).collect(),
            Look::Percent | Look::Bars => (m.block(g).iter())
                .map(|r| {
                    r.fraction()
                        .map(|f| if m.per_cell() { f } else { 100.0 * f })
                })
                .collect(),
            Look::Z | Look::Hinton => standardized(m.block(g)),
        })
        .collect();
    (rows, values)
}

/// A plot's frame: its y range ticked at the ends and the middle on the
/// left of `(x, y, w, h)`, its x range under it, and their labels.
fn axes(
    c: &mut Canvas,
    (x, y, w, h): (f64, f64, f64, f64),
    (x0, x1): (f64, f64),
    (y0, y1): (f64, f64),
    (xl, yl): (&str, &str),
    f: f64,
) {
    c.line(x, y + h, x + w, y + h, MUTED, 0.8);
    c.line(x, y, x, y + h, MUTED, 0.8);
    for (v, at) in [(y1, y), ((y0 + y1) / 2.0, y + h / 2.0), (y0, y + h)] {
        c.line(x - f * 0.3, at, x, at, MUTED, 0.8);
        c.text(
            x - f * 0.5,
            at + f * 0.35,
            &tick(v),
            f * 0.9,
            Anchor::End,
            MUTED,
        );
    }
    for k in 0..=4 {
        let (v, at) = (x0 + (x1 - x0) * k as f64 / 4.0, x + w * k as f64 / 4.0);
        c.line(at, y + h, at, y + h + f * 0.3, MUTED, 0.8);
        c.text(
            at,
            y + h + f * 1.3,
            &format!("{v:.0}"),
            f * 0.9,
            Anchor::Middle,
            MUTED,
        );
    }
    c.text(x + w / 2.0, y + h + f * 2.6, xl, f, Anchor::Middle, INK);
    c.rotated(x - TICK * f - f, y + h / 2.0, yl, f, Anchor::Middle, INK);
}

/// An axis value in a few characters: thousands and millions shortened.
fn tick(v: f64) -> String {
    match v.abs() {
        a if a >= 1e6 => format!("{:.3}M", v / 1e6),
        a if a >= 1e4 => format!("{:.1}k", v / 1e3),
        _ => format!("{v:.2}"),
    }
}

/// A tick label's width at most, in font sizes.
const TICK: f64 = 4.5;

/// The trace: each restart's ELBO over its steps, in a colour of its own,
/// the end of its k-means start marked.
fn trace(
    m: &Modules,
    title: &str,
    size: Page,
    ink: Ink,
    titles: &AxisTitles,
) -> Result<String, String> {
    let t = m.fit.trace.as_ref().map_err(Clone::clone)?;
    let points: Vec<(f64, f64)> = t
        .restarts
        .iter()
        .flatten()
        .map(|&(x, v)| (x as f64, v))
        .collect();
    if points.is_empty() {
        return Err("the trace holds no points".into());
    }
    let (x0, x1) = points
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.0), a.1.max(p.0)));
    let (y0, y1) = points
        .iter()
        .fold((f64::MAX, f64::MIN), |a, p| (a.0.min(p.1), a.1.max(p.1)));
    let f = ink.f;
    let (w, h) = size::fill(size, || (480.0, 300.0));
    // The plot alone: no title over it, no restarts named beside it.
    let (left, top, right, bottom) = (
        TICK * f + 2.5 * f,
        f + ink.room(2.0 * f),
        f + ink.room(width("restart 0", f)),
        3.5 * f,
    );
    let (pw, ph) = (w - left - right, h - top - bottom);
    let sx = |x: f64| left + pw * (x - x0) / (x1 - x0).max(1e-12);
    let sy = |v: f64| top + ph * (1.0 - (v - y0) / (y1 - y0).max(1e-12));
    let colours: &[&str] = match ink.gray {
        true => &["#1a1a1a", "#5c5c5c", "#8c8c8c", "#b3b3b3"],
        false => &[
            "#2a78d6", "#e34948", "#1f9e89", "#d89b23", "#7c5cc4", "#c2549d", "#5c7d8a", "#8a5c3b",
        ],
    };
    Ok(svg(w, h, |c: &mut Canvas| {
        ink.heading(c, title, w);
        let names = (titles.x("k-means iteration, then epoch"), titles.y("ELBO"));
        let names = (names.0.as_str(), names.1.as_str());
        axes(c, (left, top, pw, ph), (x0, x1), (y0, y1), names, f);
        if t.warmup as f64 > x0 && (t.warmup as f64) < x1 {
            let x = sx(t.warmup as f64);
            c.dashed(x, top, x, top + ph, FAINT, 0.8);
        }
        for (i, restart) in t.restarts.iter().enumerate() {
            let line: Vec<(f64, f64)> = restart
                .iter()
                .map(|&(x, v)| (sx(x as f64), sy(v)))
                .collect();
            let colour = colours[i % colours.len()];
            c.polyline(&line, colour, 1.5);
            if ink.panel_only {
                continue;
            }
            let name = format!("restart {i}");
            c.text(
                left + pw + f,
                top + 1.3 * f * (i as f64 + 1.0),
                &name,
                f,
                Anchor::Start,
                colour,
            );
        }
    }))
}

/// The cell modules against the batches, by condition: each module's share
/// in each batch (its cells by batch size), its evenness and label.
fn batches(m: &Modules, title: &str, size: Page, ink: Ink) -> Result<String, Refused> {
    let b = (m.batches.as_ref()).ok_or_else(|| String::from("the cells are of one batch"))?;
    let d = &b.design;
    let order = d.order();
    let f = ink.f;
    let said: Vec<String> = (m.columns.iter())
        .map(|&r| {
            let s = &b.spread[r];
            let even = s.evenness.map_or("-".into(), |e| format!("{e:.2}"));
            format!("{even}  {}", s.label.word())
        })
        .collect();
    let longest =
        |names: &mut dyn Iterator<Item = &String>| names.map(|s| width(s, f)).fold(0.0, f64::max);
    // Under the title, the batches' names upright, read from below; the
    // plot alone has no title over it and no conditions under it.
    let names = longest(&mut d.batches.iter()) * 1.1 + f;
    let (label, head) = (3.5 * f, ink.room(2.6 * f) + names);
    let (right, foot) = (longest(&mut said.iter()) + 2.0 * f, f + ink.room(2.0 * f));
    let (cols, n) = (order.len(), m.columns.len());
    let counts = (cols as f64, n as f64);
    let (cell, (w, h)) = cells(size, 40.0, (label + right, head + foot), counts)?;
    let font = (cell * 0.28).min(f);
    Ok(svg(w, h, |c: &mut Canvas| {
        ink.heading(c, title, w);
        for (j, &bt) in order.iter().enumerate() {
            let x = label + cell * (j as f64 + 0.5);
            c.rotated(
                x + f * 0.35,
                head - f * 0.6,
                &esc(&d.batches[bt]),
                f,
                Anchor::Start,
                INK,
            );
        }
        let right_x = label + cell * cols as f64 + f;
        c.text(
            right_x,
            head - f * 0.6,
            "even  label",
            f,
            Anchor::Start,
            MUTED,
        );
        for (i, (&role, said)) in m.columns.iter().zip(&said).enumerate() {
            let s = &b.spread[role];
            let y = head + cell * i as f64;
            c.text(
                label - f * 0.6,
                y + cell / 2.0 + f / 2.5,
                &format!("C{}", role + 1),
                f,
                Anchor::End,
                INK,
            );
            let top = s.per_batch.iter().copied().fold(0.0, f64::max);
            for (j, &bt) in order.iter().enumerate() {
                let p = s.per_batch[bt];
                let (colour, dark) = ink.fill(p, top, false);
                let x = label + cell * j as f64;
                c.rect(x + 1.0, y + 1.0, cell - 2.0, cell - 2.0, &colour);
                let text_ink = if dark { "#ffffff" } else { INK };
                c.text(
                    x + cell / 2.0,
                    y + cell / 2.0 + font / 2.5,
                    &format!("{:.0}%", 100.0 * p),
                    font,
                    Anchor::Middle,
                    text_ink,
                );
            }
            let colour = match s.label.is_batch() && !ink.gray {
                true => crate::figure::ACCENT,
                false => INK,
            };
            c.text(
                right_x,
                y + cell / 2.0 + f / 2.5,
                &esc(said),
                f,
                Anchor::Start,
                colour,
            );
        }
        if !ink.panel_only {
            let listing = d.listing().join(" · ");
            c.text(label, h - f, &esc(&listing), f * 0.9, Anchor::Start, MUTED);
        }
    }))
}
