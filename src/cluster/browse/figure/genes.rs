//! The pinned genes as the browser's genes pane draws them: side by side
//! in rows (packed as the browser packs them), each in its main gene
//! module, its gene model over it (the pileup figure's own drawing), then
//! a track for all fitted cells and one for each cell module, each with
//! its soft size on the left: the posterior-mean reads at each of the
//! gene's sites, methylated rising and unmethylated hanging, the "all"
//! track on a scale of its own and the cell modules' on one shared scale,
//! as the browser's are; the gene's symbol and module under it.

use super::{width, Ink, BLOCK_UP, BLOCK_UP_GRAY, MIRROR_LEAST, MIRROR_OFF};
use crate::cluster::browse::data::bounds;
use crate::cluster::browse::view::{columns, pack, peak, segment_width, Pin};
use crate::cluster::browse::{Count, Modules};
use crate::figure::{esc, Anchor, Canvas, FAINT, INK, MUTED};
use crate::site_analysis::miami::bin::BinEdges;
use crate::site_analysis::miami::genemodel::{draw_gene_model, GeneModel};
use crate::site_analysis::pileup::query_symbol;
use rustc_hash::FxHashMap;

/// A terminal column's width, in font sizes: the segments are as many
/// columns wide as the browser makes them, and packed in rows as it packs
/// them.
const COLUMN: f64 = 0.6;
/// The model's band and the names' line, in font sizes.
const MODEL: f64 = 1.6;
const NAME: f64 = 1.5;

/// `pins` drawn in `(x, y, w, h)` on `c`, gene models from `models`.
pub(super) fn draw(
    c: &mut Canvas,
    m: &Modules,
    pins: &[Pin],
    models: &FxHashMap<Box<str>, GeneModel>,
    (x, y, w, h): (f64, f64, f64, f64),
    ink: Ink,
) {
    if pins.is_empty() {
        return;
    }
    let f = ink.f;
    let fitted: f64 = m.cell_sizes.iter().sum();
    let tracks: Vec<(String, f64, bool)> = std::iter::once(("all".to_string(), fitted, true))
        .chain((m.columns.iter()).map(|&c| (format!("C{}", c + 1), m.cell_sizes[c], false)))
        .collect();
    let label_w = (tracks.iter())
        .map(|(name, n, _)| width(&format!("{name}  {}", m.cells_label(*n)), f))
        .fold(0.0, f64::max)
        + f;
    let widths: Vec<u16> = (pins.iter())
        .map(|p| segment_width(m.genes(p.role)[p.gene].len()))
        .collect();
    let room = ((w - label_w) / (COLUMN * f)).max(1.0) as u16;
    let rows = pack(&widths, room);
    let row_h = h / rows.len().max(1) as f64;
    let model_h = if models.is_empty() { 0.0 } else { MODEL * f };
    let track_h = ((row_h - model_h - NAME * f) / tracks.len() as f64).max(1.0);
    let up = match ink.gray {
        true => BLOCK_UP_GRAY,
        false => BLOCK_UP,
    };
    let ink_of = |all: bool| match (all, ink.gray) {
        (true, false) => crate::figure::ACCENT,
        _ => INK,
    };
    for (k, row) in rows.iter().enumerate() {
        let top = y + row_h * k as f64;
        let tracks_top = top + model_h;
        for (t, (name, n, all)) in tracks.iter().enumerate() {
            let mid = tracks_top + track_h * (t as f64 + 0.5) + f / 2.5;
            c.text(x, mid, name, f, Anchor::Start, ink_of(*all));
            let size = m.cells_label(*n);
            c.text(x + label_w - f, mid, &size, f, Anchor::End, MUTED);
        }
        // Each gene as many columns wide as the browser makes it, a
        // divider a column, from the left: the rest of the row left blank.
        let col = COLUMN * f;
        let mut at = x + label_w;
        for p in row.clone() {
            if p > row.start {
                let dx = at + col / 2.0;
                let bottom = tracks_top + track_h * tracks.len() as f64;
                c.line(dx, top, dx, bottom, FAINT, 0.5);
                at += col;
            }
            let seg = Segment {
                x: at,
                top,
                col,
                cols: widths[p],
            };
            let gene = m.gene_of(pins[p].role, m.genes(pins[p].role)[pins[p].gene].start);
            if let Some(model) = models.get(gene) {
                seg.model(c, model, model_h, f);
            }
            seg.tracks(c, m, pins[p], (model_h, track_h), up, models.get(gene));
            let name = format!("{} M{}", query_symbol(gene), pins[p].role + 1);
            let ny = tracks_top + track_h * tracks.len() as f64 + 1.1 * f;
            c.text(at, ny, &esc(&name), f, Anchor::Start, INK);
            at += col * f64::from(widths[p]);
        }
    }
}

/// One pinned gene's place: its left edge, its row's top, its columns'
/// width and how many.
struct Segment {
    x: f64,
    top: f64,
    col: f64,
    cols: u16,
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
    fn model(&self, c: &mut Canvas, g: &GeneModel, model_h: f64, f: f64) {
        let (x0, w) = self.span();
        let edges = BinEdges::new(g.lo, g.hi - 1, usize::from(self.cols));
        let blank = GeneModel {
            symbol: "".into(),
            ..g.clone()
        };
        let band = (0.45 * model_h).min(0.8 * f);
        draw_gene_model(c, &blank, &edges, x0, w, self.top + model_h / 2.0, band);
    }

    /// The tracks: all cells, then each cell module, each column's reads
    /// of the gene's sites there as a mirror bar.
    fn tracks(
        &self,
        c: &mut Canvas,
        m: &Modules,
        pin: Pin,
        (model_h, track_h): (f64, f64),
        up: &str,
        model: Option<&GeneModel>,
    ) {
        let span = m.genes(pin.role)[pin.gene].clone();
        let members = &m.members[pin.role][span];
        let pos: Vec<i64> = members.iter().map(|s| m.sites[s.site].pos).collect();
        let extent = match model {
            Some(g) => (g.lo, g.hi - 1),
            None => bounds(pos.iter().copied()),
        };
        let n = 1 + m.columns.len();
        let mut cols = vec![vec![Count::default(); n]; usize::from(self.cols)];
        for (member, col) in members.iter().zip(columns(&pos, extent, self.cols)) {
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
            true => (reach * reads / peak).max(least),
            false => 0.0,
        };
        for t in 0..n {
            let mid = self.top + model_h + track_h * (t as f64 + 0.5);
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
