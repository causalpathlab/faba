//! The pinned genes as the browser's genes pane draws them: side by side
//! in rows (packed as the browser packs them), each in its main gene
//! module, its gene model over it (the pileup figure's own drawing), then
//! a track for all fitted cells and one for each cell module, each with
//! its soft size on the left: the posterior-mean reads at each of the
//! gene's sites, methylated rising and unmethylated hanging, the "all"
//! track on a scale of its own and the cell modules' on one shared scale,
//! as the browser's are, on the scale asked for (linear, as the browser
//! opens; sqrt; log10(1 + reads)); the gene's symbol and module under it.
//! The model can go under the tracks instead (the name then over them),
//! the "all" track and the tracks' labels can be left out, and beside a
//! transposed grid the genes are one row, each cell module's track on
//! its grid row (see [`Aligned`]), every gene as wide as the widest and
//! drawn at its sites (`--genes-pad`).

use super::screen::{GenesLook, ModelAt, YScale};
use super::{width, Ink, BLOCK_UP, BLOCK_UP_GRAY, MIRROR_LEAST, MIRROR_OFF};
use crate::cluster::browse::data::bounds;
use crate::cluster::browse::view::{columns_on, pack, peak, segment_width, Pin};
use crate::cluster::browse::{Count, Modules};
use crate::figure::{Anchor, Canvas, FAINT, INK, MUTED};
use crate::site_analysis::miami::bin::BinEdges;
use crate::site_analysis::miami::genemodel::{draw_gene_model, GeneModel};
use crate::site_analysis::pileup::tui::Window;
use crate::site_analysis::pileup::{query_symbol, window_at_sites};
use crate::tui::ui::Scale;
use rustc_hash::FxHashMap;

/// A terminal column's width, in font sizes: the segments are as many
/// columns wide as the browser makes them, and packed in rows as it packs
/// them.
const COLUMN: f64 = 0.6;
/// The model's band and the names' line, in font sizes.
const MODEL: f64 = 1.6;
const NAME: f64 = 1.5;
/// A line of a name, in font sizes: a second line takes this more.
const LINE: f64 = 1.1;

/// The room the names take, in font sizes: one line, or two.
fn name_room(look: GenesLook) -> f64 {
    match look.two_line_names {
        true => NAME + LINE,
        false => NAME,
    }
}

/// How the genes are drawn: their bars' scale, the look asked for, and
/// for a transposed screen the grid's rows they share (see [`Aligned`]).
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct How {
    pub scale: YScale,
    pub look: GenesLook,
    pub aligned: Option<Aligned>,
}

impl How {
    /// As a screen's state and ask say, in packed rows.
    pub fn of(state: &super::screen::Placed, ask: &super::screen::Ask) -> Self {
        Self {
            scale: state.genes_yscale,
            look: ask.genes,
            aligned: None,
        }
    }
}

/// The grid's rows a transposed screen's genes share: the top of its
/// first cell module's row and the rows' pitch, in the canvas's points.
/// The genes are then one row of them, in the order given, each cell
/// module's track on its grid row (the "all" track a row over the first).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Aligned {
    pub top: f64,
    pub pitch: f64,
}

/// The width `pins` want in one row, at the browser's columns, with
/// their tracks' labels (unless `how` leaves them out).
pub(super) fn wanted(m: &Modules, pins: &[Pin], f: f64, how: How) -> f64 {
    let tracks = tracks(m, how.look);
    label_width(m, &tracks, f, how.look) + one_row(&beside(m, pins)) * COLUMN * f
}

/// How far over a transposed grid's first row the genes beside it reach:
/// the "all" track a row (`pitch`) over it, unless left out, and over
/// that the models (when there are any and they go on top) or the names.
pub(super) fn over(look: GenesLook, models: bool, f: f64, pitch: f64) -> f64 {
    let all = if look.no_all { 0.0 } else { pitch };
    let above = match (look.model_at, models) {
        (ModelAt::Top, true) => MODEL * f,
        (ModelAt::Top, false) => 0.0,
        (ModelAt::Bottom, _) => name_room(look) * f,
    };
    all + above
}

/// Each pin's segment's width, in the browser's columns.
fn widths(m: &Modules, pins: &[Pin]) -> Vec<u16> {
    (pins.iter())
        .map(|p| segment_width(m.genes(p.role)[p.gene].len()))
        .collect()
}

/// Each pin's segment's width beside a transposed grid, in the browser's
/// columns: every one the widest's, so each has room for its name.
pub(super) fn beside(m: &Modules, pins: &[Pin]) -> Vec<u16> {
    let widest = widths(m, pins).into_iter().max().unwrap_or(0);
    vec![widest; pins.len()]
}

/// The columns segments `widths` wide take in one row, a divider between
/// each two.
pub(super) fn one_row(widths: &[u16]) -> f64 {
    let dividers = widths.len().saturating_sub(1) as f64;
    widths.iter().map(|&w| f64::from(w)).sum::<f64>() + dividers
}

/// The tracks: "all" (unless left out), then each cell module, each with
/// its cells, and whether it is "all".
fn tracks(m: &Modules, look: GenesLook) -> Vec<(String, f64, bool)> {
    let fitted: f64 = m.cell_sizes.iter().sum();
    let all = (!look.no_all).then(|| ("all".to_string(), fitted, true));
    (all.into_iter())
        .chain((m.columns.iter()).map(|&c| (format!("C{}", c + 1), m.cell_sizes[c], false)))
        .collect()
}

/// The tracks' labels' width: none when they are left out.
fn label_width(m: &Modules, tracks: &[(String, f64, bool)], f: f64, look: GenesLook) -> f64 {
    match look.no_track_labels {
        true => 0.0,
        false => {
            (tracks.iter())
                .map(|(name, n, _)| width(&format!("{name}  {}", m.cells_label(*n)), f))
                .fold(0.0, f64::max)
                + f
        }
    }
}

/// Where a row of genes falls: its tracks' top and depth, its models'
/// band's top, and its names' baseline.
struct Band {
    tracks_top: f64,
    track_h: f64,
    model_y: f64,
    name_y: f64,
}

/// `pins` drawn in `(x, y, w, h)` on `c`, gene models from `models`, as
/// `how` says.
pub(super) fn draw(
    c: &mut Canvas,
    m: &Modules,
    pins: &[Pin],
    models: &FxHashMap<Box<str>, GeneModel>,
    (x, y, w, h): (f64, f64, f64, f64),
    ink: Ink,
    how: How,
) {
    if pins.is_empty() {
        return;
    }
    let f = ink.f;
    let tracks = tracks(m, how.look);
    let label_w = label_width(m, &tracks, f, how.look);
    // Beside a transposed grid every gene as wide as the widest, so each
    // has room for its name and its sites; else as the browser has them.
    let widths = match how.aligned {
        Some(_) => beside(m, pins),
        None => widths(m, pins),
    };
    let model_h = if models.is_empty() { 0.0 } else { MODEL * f };
    let top_model = how.look.model_at == ModelAt::Top;
    // The rows of genes, and the columns' width: packed as the browser
    // packs them, or one row narrowed to fit beside a transposed grid.
    let (rows, col) = match how.aligned {
        None => {
            let room = ((w - label_w) / (COLUMN * f)).max(1.0) as u16;
            (pack(&widths, room), COLUMN * f)
        }
        Some(_) => {
            let one = std::iter::once(0..pins.len()).collect();
            let cols = one_row(&widths).max(1.0);
            (one, (COLUMN * f).min((w - label_w) / cols))
        }
    };
    let row_h = h / rows.len().max(1) as f64;
    let band = |k: usize| -> Band {
        let n = tracks.len() as f64;
        let (tracks_top, track_h) = match how.aligned {
            Some(a) => (a.top - a.pitch * (n - m.columns.len() as f64), a.pitch),
            None => {
                let top = y + row_h * k as f64;
                let names = name_room(how.look) * f;
                let track_h = ((row_h - model_h - names) / n).max(1.0);
                let over = if top_model { model_h } else { names };
                (top + over, track_h)
            }
        };
        let bottom = tracks_top + track_h * n;
        let (model_y, name_y) = match top_model {
            true => (tracks_top - model_h, bottom + 1.1 * f),
            false => (bottom, tracks_top - 0.4 * f),
        };
        Band {
            tracks_top,
            track_h,
            model_y,
            name_y,
        }
    };
    let up = match ink.gray {
        true => BLOCK_UP_GRAY,
        false => BLOCK_UP,
    };
    let ink_of = |all: bool| match (all, ink.gray) {
        (true, false) => crate::figure::ACCENT,
        _ => INK,
    };
    for (k, row) in rows.iter().enumerate() {
        let b = band(k);
        if label_w > 0.0 {
            for (t, (name, n, all)) in tracks.iter().enumerate() {
                let mid = b.tracks_top + b.track_h * (t as f64 + 0.5) + f / 2.5;
                c.text(x, mid, name, f, Anchor::Start, ink_of(*all));
                let size = m.cells_label(*n);
                c.text(x + label_w - f, mid, &size, f, Anchor::End, MUTED);
            }
        }
        // Each gene as many columns wide as the browser makes it, a
        // divider a column, from the left: the rest of the row left blank.
        let mut at = x + label_w;
        let bottom = b.tracks_top + b.track_h * tracks.len() as f64;
        let (span_top, span_bottom) = match top_model {
            true => (b.model_y, bottom),
            false => (b.tracks_top, bottom + model_h),
        };
        for p in row.clone() {
            if p > row.start {
                let dx = at + col / 2.0;
                c.line(dx, span_top, dx, span_bottom, FAINT, 0.5);
                at += col;
            }
            let seg = Segment {
                x: at,
                col,
                cols: widths[p],
                x_scale: how.look.x_scale.into(),
            };
            let gene = m.gene_of(pins[p].role, m.genes(pins[p].role)[pins[p].gene].start);
            let model = models.get(gene);
            let extent = extent(m, pins[p], model, how);
            if let Some(model) = model {
                seg.model(c, model, extent, (b.model_y, model_h), f);
            }
            let reads = (b.tracks_top, b.track_h, how.scale, !how.look.no_all);
            seg.tracks(c, m, pins[p], reads, up, extent);
            // The name kept to its gene's width, lest it run into the next:
            // the symbol shortened, its module kept, a letter clear of the
            // next (a symbol's capitals set wider than the estimate); in
            // two lines, the symbol over its module, the whole width its.
            let module = format!("M{}", pins[p].role + 1);
            let wide = col * f64::from(widths[p]) - f;
            let lines = match how.look.two_line_names {
                true => {
                    let symbol = crate::figure::fit_text(&query_symbol(gene), wide, f);
                    vec![symbol, module]
                }
                false => {
                    let room = wide - crate::figure::text_width(&format!(" {module}"), f);
                    let symbol = crate::figure::fit_text(&query_symbol(gene), room, f);
                    vec![format!("{symbol} {module}")]
                }
            };
            // Under the tracks the first line at the names' baseline; over
            // them the last.
            let first = match top_model {
                true => b.name_y,
                false => b.name_y - LINE * f * (lines.len() - 1) as f64,
            };
            for (k, line) in lines.iter().enumerate() {
                let y = first + LINE * f * k as f64;
                c.text(at, y, line, f, Anchor::Start, INK);
            }
            at += col * f64::from(widths[p]);
        }
    }
}

/// The positions a pin's segment spans, as `how` asks: the whole gene
/// (its model's, else its sites') or at its sites (see
/// [`window_at_sites`], padded as asked) within the model; by default at
/// its sites beside a transposed grid, else the whole gene.
fn extent(m: &Modules, pin: Pin, model: Option<&GeneModel>, how: How) -> (i64, i64) {
    let span = m.genes(pin.role)[pin.gene].clone();
    let pos = m.members[pin.role][span]
        .iter()
        .map(|s| m.sites[s.site].pos);
    let sites = bounds(pos);
    let gene = model.map_or(sites, |g| (g.lo, g.hi - 1));
    let window = (how.look.window).unwrap_or(match how.aligned {
        Some(_) => Window::Sites,
        None => Window::Gene,
    });
    match window {
        Window::Gene => gene,
        Window::Sites => {
            let pad = how.look.pad.map(i64::from);
            window_at_sites(sites, pad, (gene.0.min(sites.0), gene.1.max(sites.1)), 1)
        }
    }
}

/// One pinned gene's place: its left edge, its columns' width and how
/// many.
struct Segment {
    x: f64,
    col: f64,
    cols: u16,
    /// The positions' scale across it, from its left end.
    x_scale: Scale,
}

impl Segment {
    /// The span its inner columns cover, as the browser maps positions
    /// to columns: from the first inner column's middle to the last's.
    fn span(&self) -> (f64, f64) {
        let last = self.cols.saturating_sub(2).max(1);
        let x0 = self.x + 1.5 * self.col;
        (x0, f64::from(last - 1) * self.col)
    }

    /// The gene's model across the segment, in the pileup figure's own
    /// drawing, its symbol left to the name under the tracks.
    fn model(
        &self,
        c: &mut Canvas,
        g: &GeneModel,
        (lo, hi): (i64, i64),
        (top, model_h): (f64, f64),
        f: f64,
    ) {
        let (x0, w) = self.span();
        let edges = BinEdges::new(lo, hi, usize::from(self.cols)).on(self.x_scale);
        let blank = GeneModel {
            symbol: "".into(),
            ..g.clone()
        };
        let band = (0.45 * model_h).min(0.8 * f);
        draw_gene_model(c, &blank, &edges, x0, w, top + model_h / 2.0, band);
    }

    /// The tracks from `top`, `track_h` deep each: all cells (when
    /// `all`), then each cell module, each column's reads of the gene's
    /// sites there as a mirror bar on `scale`.
    fn tracks(
        &self,
        c: &mut Canvas,
        m: &Modules,
        pin: Pin,
        (top, track_h, scale, all): (f64, f64, YScale, bool),
        up: &str,
        extent: (i64, i64),
    ) {
        let span = m.genes(pin.role)[pin.gene].clone();
        let members = &m.members[pin.role][span];
        let pos: Vec<i64> = members.iter().map(|s| m.sites[s.site].pos).collect();
        let on = BinEdges::new(extent.0, extent.1, 1).on(self.x_scale);
        let n = 1 + m.columns.len();
        let mut cols = vec![vec![Count::default(); n]; usize::from(self.cols)];
        for (member, col) in members.iter().zip(columns_on(&pos, &on, self.cols)) {
            for (t, &r) in member.reads.iter().enumerate() {
                cols[usize::from(col)][0].add(r);
                cols[usize::from(col)][t + 1].add(r);
            }
        }
        let all_peak = peak(cols.iter().map(|c| &c[0])).max(1e-12);
        let cell_peak = peak(cols.iter().flat_map(|c| &c[1..])).max(1e-12);
        // A column wide, as the browser's; a side with reads an eighth of
        // a terminal line deep at least.
        let bw = self.col;
        let least = MIRROR_LEAST * self.col / COLUMN;
        let reach = track_h / 2.0 * 0.9;
        let side = |reads: f64, peak: f64| match reads > 0.0 {
            true => (reach * scale.share(reads, peak)).max(least),
            false => 0.0,
        };
        let first = usize::from(!all);
        for t in first..n {
            let mid = top + track_h * ((t - first) as f64 + 0.5);
            let peak = if t == 0 { all_peak } else { cell_peak };
            let (x0, x1) = (
                self.x + self.col,
                self.x + self.col * f64::from(self.cols - 1),
            );
            c.line(x0, mid, x1, mid, FAINT, 0.3);
            for (j, col) in cols.iter().enumerate() {
                let r = col[t];
                if r.total <= 0.0 {
                    continue;
                }
                let bx = self.x + self.col * j as f64;
                let (u, d) = (side(r.up(), peak), side(r.off(), peak));
                c.rect(bx, mid - u, bw, u, up);
                c.rect(bx, mid, bw, d, MIRROR_OFF);
            }
        }
    }
}
