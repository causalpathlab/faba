//! The view saved as a figure: the window's rows drawn on a canvas, the
//! gene models on top as the figure draws them. `s` saves it at the view's
//! shape; `faba pileup --batch-process` draws the same rows on the page
//! asked for ([`PileupView::figure_on`]), and [`export`] many genes' at
//! once, each on its own page and side by side, as a stack's `s` does.

use super::{fmt_thousands, Bins, PileupView, Row, Show, GENE_BAND, GENE_LANE};
use crate::figure::size::PT;
use crate::figure::style::{Sheet, Style};
use crate::figure::{self, Anchor, Bars, Canvas, Panel, FAINT, INK, MUTED};
use crate::site_analysis::miami::genemodel::{draw_gene_body, draw_gene_model};
use crate::site_analysis::pileup::query_symbol;

/// Exon height and lane pitch of the gene models on a tight page.
const TIGHT_BAND: f64 = 9.0;
const TIGHT_LANE: f64 = 14.0;

/// The least a plot row's bars take, each side of a mirror's zero line,
/// laid out.
const LEAST_ARM: f64 = 10.0;

/// The least a plot takes across, laid out: room for its y labels and
/// some bars.
const LEAST_WIDE: f64 = 64.0;

/// Where a figure puts the key to its colours.
#[derive(clap::ValueEnum, Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum KeyAt {
    /// On the figure's page, under the plots (under the whole row of a
    /// combined figure).
    #[default]
    Inside,
    /// On the figure's page, beside the plots (beside the whole row).
    Right,
    /// On a page of its own, `{out}_key`: one key for a row of panels.
    Alone,
}

impl KeyAt {
    /// The next place `K` steps to.
    pub(super) fn next(self) -> Self {
        match self {
            KeyAt::Inside => KeyAt::Right,
            KeyAt::Right => KeyAt::Alone,
            KeyAt::Alone => KeyAt::Inside,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            KeyAt::Inside => "inside",
            KeyAt::Right => "right",
            KeyAt::Alone => "alone",
        }
    }
}

/// What a figure's rows take that never stretches, laid out: the head
/// (title and the line under it), the genes, the key under the plots and
/// the foot down; the key beside them across.
struct Fixed {
    head: f64,
    genes: f64,
    key: f64,
    foot: f64,
    right: f64,
}

impl Fixed {
    fn sum(&self) -> f64 {
        self.head + self.genes + self.key + self.foot
    }
}

/// A key: each colour (`"front"` for the palette's converted one) and
/// what it names.
type Key = [(&'static str, String); 2];

impl PileupView<'_> {
    /// The visible window as a figure at the view's shape, one panel per
    /// row, the genes on top.
    pub(super) fn figure(&self) -> String {
        let style = Style::default();
        let rows = self.figure_rows();
        let n_plots = Self::plot_count(&rows) as f64;
        let room = figure::tick_room(&self.tick_list(), true);
        let fixed = self.fixed(&rows, &style, None);
        // The shape the view has on screen, grown alike on both sides until
        // each plot has room for its axis; else a page of its own.
        let least = 100.0 + room;
        let page = self.controls.page(
            self.drawn,
            (360.0 + fixed.right, fixed.sum() + least * n_plots),
        );
        let size = page.unwrap_or((720.0 + fixed.right, fixed.sum() + (170.0 + room) * n_plots));
        self.draw_figure(size, &style)
    }

    /// What `s` saves: the figure, and its key on a page of its own
    /// (`_key`) when it goes alone.
    pub(super) fn saved(&self) -> Vec<(String, String)> {
        let mut files = vec![(String::new(), self.figure())];
        if self.key_at == KeyAt::Alone {
            files.extend(
                self.key_page(&Style::default(), None)
                    .map(|k| ("_key".into(), k)),
            );
        }
        files
    }

    /// The window as a figure on `sheet`'s page, a bar a column of the
    /// view (the bases per bar it shows): laid out from its font, tight to
    /// its labels; refused only when its labels cannot fit, naming the
    /// least size that is drawn.
    pub(crate) fn figure_on(&self, sheet: &Sheet) -> anyhow::Result<String> {
        let n_plots = Self::plot_count(&self.figure_rows()) as f64;
        let size = self.page_on(sheet);
        let (style, short) = self.least_on(sheet, size);
        if let Some(mut least) = short {
            // A font fit to the page grows with it: the least named is one
            // its own font fits.
            for _ in 0..8 {
                match self.least_on(sheet, least).1 {
                    Some(more) => least = more,
                    None => break,
                }
            }
            anyhow::bail!(
                "a pileup figure of {n_plots} plot(s) at {:.1} pt is at least {:.2} by {:.2} inches: its labels need the room",
                self.least_on(sheet, least).0.font,
                least.0 / PT,
                least.1 / PT
            );
        }
        Ok(self.draw_figure(size, &style))
    }

    /// The page `sheet` gives the figure, in points: as asked, else its
    /// own size (at 8-point text) grown to the font asked.
    fn page_on(&self, sheet: &Sheet) -> (f64, f64) {
        let rows = self.figure_rows();
        let n_plots = Self::plot_count(&rows) as f64;
        let style = Style {
            tight: true,
            tick_room: self.tick_room,
            ..sheet.style((360.0, 360.0))
        };
        let fixed = self.fixed(&rows, &style, Some(360.0));
        let h = fixed.sum() + n_plots * (self.least_plot() + 60.0);
        sheet.fill_own((360.0 + fixed.right, h))
    }

    /// The room the position labels take under the plots on `sheet`'s
    /// page, laid out: none when they lie flat. Each row's ticks are
    /// chosen as it draws them, from where its plot starts.
    fn tick_room_on(&self, sheet: &Sheet) -> f64 {
        let size = self.page_on(sheet);
        let style = Style {
            tight: true,
            ..sheet.style(size)
        };
        let mut c = Canvas::styled(size.0, size.1, &style);
        let rows = self.figure_rows();
        let fixed = self.fixed(&rows, &style, Some(c.page().0));
        let w = c.page().0 - fixed.right - 4.0;
        let bins = self.bins();
        // A row's depth, about: its labels are chosen for it.
        let row_h = c.page().1 / rows.len().max(1) as f64;
        let room = |p: &mut dyn Panel, c: &mut Canvas| {
            let (ticks, turned) = self.ticks_for(c, w, p.left(c, row_h));
            figure::tick_room(&ticks, turned)
        };
        (rows.into_iter())
            .filter_map(|row| self.with_panel(row, &bins, &mut c, None, String::new(), room))
            .fold(0.0, f64::max)
    }

    /// The look on a page `size` points, tight, and the least page (in
    /// points, each side at least as asked, rounded up to a hundredth of
    /// an inch) when its rows cannot fit `size`.
    fn least_on(&self, sheet: &Sheet, size: (f64, f64)) -> (Style, Option<(f64, f64)>) {
        let style = Style {
            tight: true,
            tick_room: self.tick_room,
            ..sheet.style(size)
        };
        let em = style.em();
        let (lw, lh) = (size.0 / em, size.1 / em);
        let rows = self.figure_rows();
        let n_plots = Self::plot_count(&rows) as f64;
        let fixed = self.fixed(&rows, &style, Some(lw));
        let least = (
            LEAST_WIDE + fixed.right,
            fixed.sum() + n_plots * self.least_plot(),
        );
        let fits = lw + 1e-9 >= least.0 && lh + 1e-9 >= least.1;
        let up = |v: f64, asked: f64| ((v * em / PT * 100.0).ceil() / 100.0 * PT).max(asked);
        let short = (!fits).then(|| (up(least.0, size.0), up(least.1, size.1)));
        (style, short)
    }

    /// The rows a figure draws: the genes on top, then the plots.
    fn figure_rows(&self) -> Vec<Row> {
        super::render::genes_first(self.rows())
    }

    /// The least a plot row takes on a tight page, laid out: its margins
    /// (untitled) and both arms, its ticks (in Mb) turned as the most they
    /// take.
    fn least_plot(&self) -> f64 {
        let room = figure::tick_room(&self.mb_ticks(3), true);
        8.0 + 2.0 * LEAST_ARM + 25.0 + room
    }

    /// What `rows` take that never stretches, in `style` on a page
    /// `width` across (laid out; the key on one line when none is given).
    fn fixed(&self, rows: &[Row], style: &Style, width: Option<f64>) -> Fixed {
        let head = if style.tight { 24.0 } else { 44.0 };
        let genes = if rows.contains(&Row::Genes) {
            let lanes = self.lanes() as f64;
            match style.tight {
                true => TIGHT_LANE * lanes,
                false => 24.0 + GENE_LANE * lanes,
            }
        } else {
            0.0
        };
        let key = self.drawn_key(style);
        let (key, right) = match (key, self.key_at, width) {
            (None, ..) | (_, KeyAt::Alone, _) => (0.0, 0.0),
            (Some(k), KeyAt::Right, _) => (0.0, key_width(&k)),
            (Some(_), KeyAt::Inside, None) => (KEY_LINE, 0.0),
            (Some(k), KeyAt::Inside, Some(w)) => (KEY_LINE * key_lines(&k, w - 8.0) as f64, 0.0),
        };
        Fixed {
            head,
            genes,
            key,
            foot: 6.0,
            right,
        }
    }

    /// The colours a key names: the converted reads in front of the total,
    /// when the plots draw both.
    pub(super) fn key(&self) -> Option<Key> {
        let both = self.show == Show::Both && self.tracks.iter().any(|t| t.behind().is_some());
        both.then(|| {
            [
                ("front", self.on.to_string()),
                (FAINT, self.off.to_string()),
            ]
        })
    }

    /// The key a figure in `style` draws on its own page: none on a
    /// panel alone unless asked for (`--key`), nor where it goes alone, nor
    /// on a panel of a combined figure, which draws it once.
    fn drawn_key(&self, style: &Style) -> Option<Key> {
        let asked = !style.panel_only || self.with_key;
        let here = self.key_at != KeyAt::Alone && !self.keyless;
        self.key().filter(|_| asked && here)
    }

    /// The key on a page of its own in `style`: its entries side by side
    /// when they fit `room` across (laid out), else one under another.
    fn key_page(&self, style: &Style, room: Option<f64>) -> Option<String> {
        let key = self.key()?;
        let em = style.em();
        let (w, lines) = match room {
            Some(r) => (r, key_lines(&key, r - 8.0)),
            None => (key_width(&key), key.len()),
        };
        let h = KEY_LINE * lines as f64 + 6.0;
        let mut c = Canvas::styled(w * em, h * em, style);
        let front = style.palette.roles().front;
        draw_key(&mut c, &key, front, (4.0, 0.0, w - 8.0));
        Some(c.finish())
    }

    /// The gene the figure is of, by its symbol: the one gene model drawn
    /// of its own, else the symbol in its name.
    pub(super) fn symbol(&self) -> String {
        match self.genes.as_slice() {
            [g] => g.symbol.to_string(),
            _ => query_symbol(&self.title).to_string(),
        }
    }

    /// The figure on a page `size` points in `style`: the symbol as title,
    /// the window under it, the genes, the plots sharing what is left, and
    /// a key when the colours need one, under or beside them.
    fn draw_figure(&self, size: (f64, f64), style: &Style) -> String {
        let (lo, hi) = self.window;
        let mut c = Canvas::styled(size.0, size.1, style);
        let (width, height) = c.page();
        let rows = self.figure_rows();
        let fixed = self.fixed(&rows, style, Some(width));
        let n_plots = Self::plot_count(&rows) as f64;
        let panel_h = (height - fixed.sum()) / n_plots;
        let across = width - fixed.right;
        let (x, w) = if style.tight {
            (2.0, across - 4.0)
        } else {
            (8.0, across - 16.0)
        };
        let bins = self.bins();
        // The plots first: the genes and the title line up with their bars.
        let mut y = fixed.head + fixed.genes;
        let mut plot = (x + 52.0, w - 64.0);
        for &row in rows.iter().filter(|&&r| r != Row::Genes) {
            let title = self.row_title(row);
            let (px, _, pw, _) = self.draw_row(row, &bins, &mut c, (x, y, w, panel_h), None, title);
            plot = (px, pw);
            y += panel_h;
        }
        if rows.contains(&Row::Genes) {
            let bbox = (x, fixed.head, w, fixed.genes);
            self.draw_genes(&mut c, bbox, plot, "genes".into());
        }
        let symbol = self.symbol();
        let window = format!(
            "{}:{}-{}, {}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.per_bar()
        );
        // The window, but for the bars' width alone on a panel, or where
        // the window does not fit across.
        let narrow = figure::text_width(&window, c.label_size(9.0)) > across - 4.0;
        let under = match style.panel_only || narrow {
            true => self.per_bar(),
            false => window,
        };
        if style.tight {
            let mid = plot.0 + plot.1 / 2.0;
            // Inside the row's box, a margin from the page's edges.
            let (at, size) = figure::place_bold_title(&c, &symbol, mid, (x, w), 12.0);
            c.bold_exact(at, 11.0, &symbol, size, Anchor::Middle, INK);
            let (at, u) = figure::place_title(&c, &under, mid, (x, w), 9.0);
            c.text(at, 21.0, &u, 9.0, Anchor::Middle, MUTED);
        } else {
            c.bold(16.0, 22.0, &symbol, 12.0, Anchor::Start, INK);
            c.text(16.0, 38.0, &under, 9.0, Anchor::Start, MUTED);
        }
        if let Some(key) = self.drawn_key(style) {
            let front = style.palette.roles().front;
            match self.key_at {
                KeyAt::Right => {
                    let top = fixed.head + fixed.genes;
                    draw_key(&mut c, &key, front, (across + 2.0, top, 0.0));
                }
                _ => {
                    let top = height - fixed.foot - fixed.key;
                    draw_key(&mut c, &key, front, (plot.0.min(x + 4.0), top, width - 8.0));
                }
            }
        }
        c.finish()
    }

    /// Draw one row (track, contrast, mirror or genes) in the box; the box
    /// its bars take, the box itself for the genes. Its position ticks are
    /// chosen for where its plot starts.
    pub(super) fn draw_row(
        &self,
        row: Row,
        bins: &Bins,
        c: &mut Canvas,
        bbox: (f64, f64, f64, f64),
        pointer: Option<usize>,
        title: String,
    ) -> (f64, f64, f64, f64) {
        let (x, _, w, _) = bbox;
        if row == Row::Genes {
            self.draw_genes(c, bbox, (x + 52.0, w - 64.0), title);
            return bbox;
        }
        let draw = |p: &mut dyn Panel, c: &mut Canvas| {
            let (ticks, rotated) = self.ticks_for(c, w, p.left(c, bbox.3));
            p.set_ticks(ticks, rotated);
            p.draw_in(c, bbox)
        };
        (self.with_panel(row, bins, c, pointer, title, draw)).unwrap_or(bbox)
    }

    /// `f` of row `row`'s panel (a track, contrast or mirror; none for the
    /// genes) on `c`, as yet without position ticks.
    fn with_panel<R>(
        &self,
        row: Row,
        bins: &Bins,
        c: &mut Canvas,
        pointer: Option<usize>,
        title: String,
        f: impl FnOnce(&mut dyn Panel, &mut Canvas) -> R,
    ) -> Option<R> {
        // A tight page names the chromosome and the ticks' unit.
        // On a sqrt or log x scale, the bases from the window's left end.
        let x_title = match (self.edges().is_linear(), c.style().tight) {
            (true, true) => format!("{} (Mb)", self.chr),
            (true, false) => format!("{} position", self.chr),
            (false, _) => format!(
                "bp from {} ({})",
                fmt_thousands(self.window.0),
                self.edges().x_scale.name()
            ),
        };
        let roles = c.style().palette.roles();
        let (bar, front) = (roles.bar, roles.front);
        match row {
            Row::Track(i) => {
                let t = &self.tracks[i];
                let (lo, hi) = self.window;
                let shown = &bins.shown[i];
                let marks = t.sites_in(lo, hi).into_iter();
                let colour = |_| if shown.accent { front } else { bar };
                let mut p = Bars {
                    rotated: false,
                    values: &shown.values,
                    front: shown.front.as_deref(),
                    // Converted reads in the palette's colour for bars in
                    // front of a total, as the mirrored rows draw them.
                    accent: &|_| false,
                    colour: &colour,
                    y_scale: self.y_scale,
                    y_max: self.row_max(row, bins),
                    ticks: Vec::new(),
                    pointer,
                    marks: marks.map(|s| bins.edges.col_of(s)).collect(),
                    dividers: Vec::new(),
                    title,
                    x_title,
                    y_title: self.titled(self.axis_title(i)),
                };
                Some(f(&mut p, c))
            }
            Row::Contrast(contrast) => {
                let values = self.contrast_values(contrast, bins);
                let measure = contrast.measure;
                let label = |v| measure.label(v);
                let y_title = match c.style().tight {
                    true => measure.short().into(),
                    false => measure.name(self.on),
                };
                let mut p = figure::Diverging {
                    rotated: false,
                    values: &values,
                    ticks: Vec::new(),
                    pointer,
                    title,
                    x_title,
                    y_title: self.titled(y_title),
                    label: &label,
                };
                Some(f(&mut p, c))
            }
            Row::Mirror(a, b) => {
                let half = |k: usize| {
                    let shown = &bins.shown[k];
                    // Converted alone draws all in the converted colour.
                    let front = match (&shown.front, shown.accent) {
                        (Some(f), _) => Some(f.as_slice()),
                        (None, true) => Some(shown.values.as_slice()),
                        (None, false) => None,
                    };
                    figure::Half {
                        values: &shown.values,
                        front,
                        colour: bar,
                        name: self.tracks[k].label.to_string(),
                    }
                };
                let mut p = figure::Mirror {
                    rotated: false,
                    up: half(a),
                    down: half(b),
                    y_scale: self.y_scale,
                    y_max: self.row_max(row, bins),
                    y_labels: Default::default(),
                    ticks: Vec::new(),
                    pointer,
                    title,
                    x_title,
                    y_title: self.titled(self.axis_title(a)),
                };
                Some(f(&mut p, c))
            }
            Row::Genes => None,
        }
    }

    /// A row's y title `s`, but none on a panel not first in its row of
    /// panels, whose plot takes that room.
    fn titled(&self, s: String) -> String {
        match self.y_untitled {
            true => String::new(),
            false => s,
        }
    }

    /// What track `k`'s bars count, on a figure's y axis: as its title
    /// says, but where a key names the colours.
    fn axis_title(&self, k: usize) -> String {
        let t = &self.tracks[k];
        let sites = self.show == Show::Sites && t.is_reads();
        if sites || t.behind().is_some() {
            self.show.axis(self.on, self.off)
        } else {
            t.name.to_string()
        }
    }

    /// About five labelled genomic ticks across the window (on a sqrt or
    /// log x scale, at round distances from its left end).
    fn tick_list(&self) -> Vec<(usize, String)> {
        match self.edges().is_linear() {
            true => self.ticks_every(self.columns / 5),
            false => self.distance_ticks(5),
        }
    }

    /// At most `n` ticks on a sqrt or log x scale: the left end (0) and
    /// round distances from it, evenly spread on that scale (as a y axis
    /// on it is ticked), the farthest kept.
    fn distance_ticks(&self, n: usize) -> Vec<(usize, String)> {
        let edges = self.edges();
        let top = figure::scaled(edges.x_scale, edges.span() as f64);
        let round = figure::round_y_ticks(edges.x_scale, top, figure::SQRT_SPREAD);
        let far = round.len().saturating_sub(n.saturating_sub(1));
        let at = |share: f64| ((share * self.columns as f64) as usize).min(self.columns.max(1) - 1);
        std::iter::once((0, "0".to_string()))
            .chain(
                round[far..]
                    .iter()
                    .map(|(share, l)| (at(*share), l.clone())),
            )
            .collect()
    }

    /// The ticks on `c` for a row `w` across whose plot starts `left` in,
    /// and whether they are turned: about five, turned; on a tight page as
    /// [`PileupView::tight_ticks`].
    fn ticks_for(&self, c: &Canvas, w: f64, left: f64) -> (Vec<(usize, String)>, bool) {
        match c.style().tight {
            true => self.tight_ticks(c.label_size(8.0), w, left),
            false => (self.tick_list(), true),
        }
    }

    /// A tight page's ticks for a row `w` across whose plot starts `left`
    /// in, its labels at `size` (laid out), and whether they are turned:
    /// three or else two in Mb (a round step's, or the first and last of
    /// three), flat where they fit across the row, else three turned (the
    /// chromosome under them names the unit).
    fn tight_ticks(&self, size: f64, w: f64, left: f64) -> (Vec<(usize, String)>, bool) {
        // Where a tick's label lies across the row, as it is drawn: moved
        // in at the row's ends; none where it cannot be.
        let bw = (w - left - figure::TIGHT_RIGHT) / self.columns.max(1) as f64;
        let span = |i: usize, l: &str| {
            let half = figure::tick_width(l, size) / 2.0;
            let at = figure::flat_label(left + (i as f64 + 0.5) * bw, half, (0.0, w))?;
            Some((at - half, at + half))
        };
        // Those whose labels lie inside the row, when two or more do and
        // none runs into the next.
        let flat = |ticks: Vec<(usize, String)>| {
            let kept: Vec<_> = (ticks.into_iter())
                .filter_map(|(i, l)| Some((span(i, &l)?, (i, l))))
                .collect();
            let apart = (kept.windows(2)).all(|p| p[1].0 .0 >= p[0].0 .1 + 0.25 * size);
            (kept.len() >= 2 && apart).then(|| kept.into_iter().map(|k| k.1).collect())
        };
        // Two: a round step's, else the first and last of three.
        let three = self.mb_ticks(3);
        let ends = match three.as_slice() {
            [first, .., last] => vec![first.clone(), last.clone()],
            _ => Vec::new(),
        };
        for ticks in [three.clone(), self.mb_ticks(2), ends] {
            if let Some(ticks) = flat(ticks) {
                return (ticks, false);
            }
        }
        (three, true)
    }

    /// At most `n` ticks (two at least where the window allows) at round
    /// positions, each labelled in Mb with just the decimals that tell
    /// them apart: a step of 1, 2 or 5 times a power of ten bases.
    pub(super) fn mb_ticks(&self, n: usize) -> Vec<(usize, String)> {
        if !self.edges().is_linear() {
            return self.distance_ticks(n);
        }
        let (lo, hi) = self.window;
        let n = n as i64;
        let count = |step: i64| hi.div_euclid(step) - (lo + step - 1).div_euclid(step) + 1;
        let steps = || (0..10u32).flat_map(|e| [1, 2, 5].map(|m| m * 10i64.pow(e)));
        let step = (steps().find(|&s| (2..=n).contains(&count(s))))
            .or_else(|| steps().find(|&s| (1..=n).contains(&count(s))));
        let Some(step) = step else {
            return Vec::new();
        };
        let decimals = (0..6u32)
            .find(|&d| step % 10i64.pow(6 - d) == 0)
            .unwrap_or(6) as usize;
        let edges = self.edges();
        let first = (lo + step - 1).div_euclid(step) * step;
        (0..count(step))
            .map(|k| first + k * step)
            .map(|p| (edges.col_of(p), format!("{:.*}", decimals, p as f64 / 1e6)))
            .collect()
    }

    /// A labelled genomic tick every `every` bars (at least one).
    fn ticks_every(&self, every: usize) -> Vec<(usize, String)> {
        (0..self.columns)
            .step_by(every.max(1))
            .map(|k| (k, fmt_thousands(self.bar_range(k).0)))
            .collect()
    }

    /// The genes row in the box, as the figure draws gene models, across
    /// the plots' bars `(left, width)`. A tight page draws them smaller,
    /// and the figure's own gene, which its title names, without a label.
    fn draw_genes(
        &self,
        c: &mut Canvas,
        (x, y, _w, _h): (f64, f64, f64, f64),
        (left, width): (f64, f64),
        title: String,
    ) {
        let style = *c.style();
        let (band, lane, first) = match style.tight {
            true => (TIGHT_BAND, TIGHT_LANE, TIGHT_LANE / 2.0),
            false => (GENE_BAND, GENE_LANE, 26.0),
        };
        // A tight page's genes go untitled, as its plots do.
        if !style.panel_only && !style.tight {
            c.bold(x + 52.0, y + 11.0, &title, 11.0, Anchor::Start, INK);
        }
        let edges = self.edges();
        let symbol = self.symbol();
        for (k, g) in self.gene_lanes() {
            let mid = y + first + lane * k as f64;
            if style.tight && *g.symbol == *symbol {
                draw_gene_body(c, g, &edges, (left, width), mid, band);
            } else {
                draw_gene_model(c, g, &edges, left, width, mid, band);
            }
        }
    }
}

/// Swatch side and the room after it, and between entries, in a key; a
/// key's line.
const SWATCH: f64 = 6.0;
const KEY_GAP: f64 = 8.0;
const KEY_LINE: f64 = 12.0;

/// How wide one key entry sets at 8 points on `c`'s labels.
fn entry_width(name: &str, size: f64) -> f64 {
    SWATCH + 3.0 + figure::text_width(name, size)
}

/// How wide a key of one entry a line is, laid out, with room either side.
fn key_width(key: &[(&'static str, String)]) -> f64 {
    let size = crate::figure::style::OWN_FONT;
    let widest = key
        .iter()
        .map(|(_, n)| entry_width(n, size))
        .fold(0.0, f64::max);
    widest + 8.0
}

/// Lines a key takes `room` across: one when its entries fit side by
/// side, else a line each.
fn key_lines(key: &[(&'static str, String)], room: f64) -> usize {
    let size = crate::figure::style::OWN_FONT;
    let all: f64 = key
        .iter()
        .map(|(_, n)| entry_width(n, size) + KEY_GAP)
        .sum();
    if all - KEY_GAP <= room {
        1
    } else {
        key.len()
    }
}

/// The key at `(x, top)`, `room` across: each colour's swatch and name,
/// side by side when they fit, else one under another. `"front"` is the
/// palette's converted colour.
fn draw_key(
    c: &mut Canvas,
    key: &[(&'static str, String)],
    front: &'static str,
    (x, top, room): (f64, f64, f64),
) {
    let size = c.label_size(8.0);
    let one_line = key_lines(key, room) == 1;
    let (mut at, mut y) = (x, top + 9.0);
    for (colour, name) in key {
        let fill = if *colour == "front" { front } else { colour };
        c.rect(at, y - SWATCH, SWATCH, SWATCH, fill);
        c.text(at + SWATCH + 3.0, y, name, 8.0, Anchor::Start, MUTED);
        if one_line {
            at += entry_width(name, size) + KEY_GAP;
        } else {
            y += KEY_LINE;
        }
    }
}

/// How an export lays out its combined figure: `cols` genes to a row,
/// named `name` (`row`, `grid`, `column`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Combine {
    pub cols: usize,
    pub name: &'static str,
}

/// The files exporting `views` writes, each a name to put after the
/// output prefix and its SVG: every gene's figure on `sheet`'s page
/// (`_{SYMBOL}`), with `combine` one figure of them all side by side, each
/// panel on that page and the key drawn once for all (`_{name}`), and the
/// key on a page of its own (`_key`) when it goes alone. Every gene is
/// drawn as the first is: its y scale, channels, difference row and key.
pub fn export(
    views: &mut [PileupView],
    sheet: &Sheet,
    combine: Option<Combine>,
) -> anyhow::Result<Vec<(String, String)>> {
    let Some(first) = views.first() else {
        anyhow::bail!("no gene to draw");
    };
    let settings = first.settings();
    for v in views.iter_mut() {
        v.apply(&settings);
    }
    share_tick_room(views, sheet);
    let mut out = Vec::new();
    let mut names: Vec<String> = Vec::new();
    for v in views.iter() {
        let mut name = format!("_{}", file_part(&v.symbol()));
        // Two genes of one symbol (their ids differ) keep both files.
        if names.contains(&name) {
            name = (2..)
                .map(|k| format!("{name}_{k}"))
                .find(|n| !names.contains(n))
                .expect("some suffix is free");
        }
        names.push(name.clone());
        out.push((name, v.figure_on(sheet)?));
    }
    let style = |svg: &str| {
        let size = figure::svg_size(svg);
        Style {
            tight: true,
            ..sheet.style(size)
        }
    };
    if let Some(Combine { cols, name }) = combine {
        // Only the first panel of each row titles its y axis; the others'
        // plots take that room.
        for (k, v) in views.iter_mut().enumerate() {
            v.y_untitled = k % cols != 0;
        }
        share_tick_room(views, sheet);
        let mut panels = Vec::with_capacity(views.len());
        for v in views.iter_mut() {
            v.keyless = true;
            let panel = v.figure_on(sheet);
            v.keyless = false;
            v.y_untitled = false;
            panels.push(panel?);
        }
        let first = &views[0];
        let grid = figure::tile_svgs(&panels, cols);
        let style = style(&panels[0]);
        let (gw, _) = figure::svg_size(&grid);
        let keyed = first.drawn_key(&style).is_some();
        let joined = match first.key_at {
            KeyAt::Inside if keyed => {
                let key = first.key_page(&style, Some(gw / style.em()));
                figure::tile_svgs(&[grid, key.unwrap_or_default()], 1)
            }
            KeyAt::Right if keyed => {
                let key = first.key_page(&style, None).unwrap_or_default();
                figure::beside_svgs(&grid, &key)
            }
            _ => grid,
        };
        out.push((format!("_{name}"), joined));
    }
    let first = &views[0];
    let style = style(&out[0].1);
    let asked = !style.panel_only || first.with_key;
    if first.key_at == KeyAt::Alone && asked {
        if let Some(key) = first.key_page(&style, None) {
            out.push(("_key".into(), key));
        }
    }
    for v in views.iter_mut() {
        v.tick_room = 0.0;
    }
    Ok(out)
}

/// One room for every view's position labels on `sheet`'s page, so their
/// plots line up side by side, in one figure or in files put together.
fn share_tick_room(views: &mut [PileupView], sheet: &Sheet) {
    let room = (views.iter())
        .map(|v| v.tick_room_on(sheet))
        .fold(0.0, f64::max);
    for v in views.iter_mut() {
        v.tick_room = room;
    }
}

/// `s` as a file name's part: letters, digits, `-` and `.` kept, the rest
/// `_`.
fn file_part(s: &str) -> String {
    s.chars()
        .map(
            |c| match c.is_ascii_alphanumeric() || c == '-' || c == '.' {
                true => c,
                false => '_',
            },
        )
        .collect()
}
