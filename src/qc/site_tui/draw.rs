//! Drawing: the threshold table, the three plots, and the saved figure.

use super::*;
use crate::site_analysis::miami::genemodel::draw_gene_model;

/// Each coding region's name at its middle bin, where a short UTR's label
/// does not run into the next one.
fn meta_ticks(regions: [usize; 3]) -> Vec<(usize, String)> {
    let mut at = 0;
    let mut out = Vec::new();
    for (r, &n) in regions.iter().enumerate() {
        if n > 0 {
            out.push((at + n / 2, REGION_NAMES[r].to_string()));
        }
        at += n;
    }
    out
}

/// The first bin of the CDS and of the 3'UTR: the region boundaries.
fn meta_marks(regions: [usize; 3]) -> Vec<usize> {
    vec![regions[0], regions[0] + regions[1]]
}

/// A colour key: a block in each bar style, then what it stands for.
fn key_line(items: &[(Style, String)]) -> Line<'static> {
    let mut spans = Vec::new();
    for (style, what) in items {
        spans.push(Span::styled("█ ", *style));
        spans.push(Span::styled(format!("{what}   "), DIM));
    }
    Line::from(spans)
}

/// The key of the gene and metagene plots, with the key that switches
/// what they add up.
fn kept_key(weight: Weight) -> Line<'static> {
    let mut line = key_line(&[
        (DIM, format!("all {}", weight.unit())),
        (PLAIN, "kept".into()),
    ]);
    line.push_span(Span::styled(format!("c: {}", weight.other().unit()), DIM));
    line
}

/// Bars of `all` sites behind the `kept` ones, in the box `(x, y, w, h)`.
#[allow(clippy::too_many_arguments)]
fn kept_bars(
    canvas: &mut Canvas,
    (all, kept): (&[usize], &[usize]),
    (x, y, w, h): (f64, f64, f64, f64),
    ticks: Vec<(usize, String)>,
    marks: Vec<usize>,
    title: String,
    x_title: &str,
    y_title: &str,
) {
    let f = |v: &[usize]| v.iter().map(|&n| n as f64).collect::<Vec<_>>();
    let (all, kept) = (f(all), f(kept));
    Bars {
        values: &all,
        front: Some(&kept),
        accent: &|_| false,
        y_scale: Scale::Linear,
        y_max: None,
        ticks,
        pointer: None,
        marks,
        title,
        x_title: x_title.into(),
        y_title: y_title.into(),
    }
    .draw(canvas, x, y, w, h);
}

/// A gene's profile in the box `(x, y, w, h)`, with its gene model in the
/// space Bars keeps for x labels (52 left, 12 right, 38 below).
fn draw_gene(canvas: &mut Canvas, p: &GeneProfile, (x, y, w, h): (f64, f64, f64, f64)) {
    let bars = (&p.all[..], &p.kept[..]);
    kept_bars(
        canvas,
        bars,
        (x, y, w, h),
        Vec::new(),
        Vec::new(),
        p.title(),
        "",
        p.unit,
    );
    if let Some(m) = &p.model {
        draw_gene_model(canvas, m, &p.edges, x + 52.0, w - 64.0, y + h - 24.0, 8.0);
    }
}

/// What the table's two counts mean, a line each for the screen and the
/// figure.
const LEGEND: [&str; 4] = [
    "filtered out: sites failing this threshold.",
    "  A site failing several thresholds counts in",
    "  each row, so the rows add up to more than the",
    "  sites removed.",
];

/// What the figure's bar colours mean; on screen each plot has its own key.
const FIGURE_KEY: [&str; 3] = [
    "Bars: light gray, all sites; dark gray, kept;",
    "  orange, filtered out only by the focused",
    "  threshold.",
];

/// The count column's header, a line each.
const COUNT_HEADER: (&str, &str) = ("sites", "filtered out");

/// Widths of the table's columns: marker and name, threshold, count.
const TABLE_WIDTHS: [usize; 3] = [17, 10, 14];

impl<'a> SitePicker<'a> {
    pub(super) fn render_gene_list(&self, frame: &mut Frame, area: Rect) {
        let block = panel(" genes: kept / all sites ".into(), false);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(genes) = &self.view().genes else {
            let why = "this site table has no gene column";
            frame.render_widget(Paragraph::new(Span::styled(why, DIM)), inner);
            return;
        };
        let list = self.gene_list();
        let at = self.list.at;
        let rows = inner.height.saturating_sub(1) as usize;
        let first = at.saturating_sub(rows.saturating_sub(1) / 2);
        let first = first.min(list.len().saturating_sub(rows));
        let width = inner.width as usize;
        let mut lines: Vec<Line> = list
            .iter()
            .enumerate()
            .skip(first)
            .take(rows)
            .map(|(j, &g)| {
                let g = g as usize;
                let (kept, all) = self.gene_kept(g);
                let counts = format!("{kept} / {all}");
                let name_w = width.saturating_sub(counts.len() + 3);
                let selected = j == at;
                let style = if selected { HIGHLIGHT } else { PLAIN };
                Line::from(vec![
                    Span::styled(if selected { "▸ " } else { "  " }, HIGHLIGHT),
                    Span::styled(format!("{:<name_w$.name_w$}", genes.symbol(g)), style),
                    Span::styled(format!(" {counts}"), DIM),
                ])
            })
            .collect();
        if list.is_empty() {
            lines.push(Line::from(Span::styled("  no gene matches", DIM)));
        }
        let find = match self.mode {
            Mode::Find => input_line("gene: ", &self.list.find, &[]),
            _ if !self.list.find.is_empty() => {
                input_line("gene: ", &self.list.find, &[("[ ]", "move"), ("/", "find")])
            }
            _ => {
                let mut hint = help_line(&[("[ ]", "move"), ("/", "find")]);
                hint.push_span(Span::styled(format!("  {} genes", list.len()), DIM));
                hint
            }
        };
        let [list_area, find_area] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(1)]).areas(inner);
        frame.render_widget(Paragraph::new(lines), list_area);
        frame.render_widget(find, find_area);
    }

    pub(super) fn render_gene(&self, frame: &mut Frame, area: Rect) {
        let n = area.width.saturating_sub(2 + GUTTER).max(10) as usize;
        let profile = self.gene_profile(n);
        let title = match &profile {
            Some(p) => format!(" {} ", p.title()),
            None => " gene ".into(),
        };
        let block = panel(title, true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(p) = profile else {
            let why = "no gene selected";
            frame.render_widget(Paragraph::new(Span::styled(why, DIM)), inner);
            return;
        };
        let [plot, track, key] = Layout::vertical([
            Constraint::Min(4),
            Constraint::Length(1),
            Constraint::Length(1),
        ])
        .areas(inner);
        frame.render_widget(Paragraph::new(kept_key(self.weight)), key);
        HistPlot {
            bins: Binning::with_width(Scale::Linear, 1.0),
            kmin: 0,
            counts: &p.all,
            style: &|_| PLAIN,
            subset: Some(&p.kept),
            y_scale: Scale::Linear,
            y_max: None,
            pointer: None,
            marks: Vec::new(),
            x_label: Some(&|_| None),
            tick_every: None,
        }
        .render(frame.buffer_mut(), plot);
        let model = match p.model {
            Some(_) => (0..p.edges.num_bins)
                .map(|b| {
                    if p.exonic(b) == Some(true) {
                        "▬"
                    } else {
                        "─"
                    }
                })
                .collect::<String>(),
            None => "(no gene model)".into(),
        };
        let line = Line::from(vec![
            Span::styled(format!("{:<6}", "exons"), DIM),
            Span::styled(model, DIM),
        ]);
        frame.render_widget(Paragraph::new(line), track);
    }

    /// The metagene's title: what the bars add up.
    fn meta_title(&self) -> String {
        format!("metagene · y: {} per bin", self.weight.unit())
    }

    /// The metagene in the box `bbox`: all sites behind, the kept in front.
    fn draw_meta(
        &self,
        canvas: &mut Canvas,
        m: &MetaCounts,
        bbox: (f64, f64, f64, f64),
        title: String,
    ) {
        let (ticks, marks) = (meta_ticks(m.regions), meta_marks(m.regions));
        let x_title = "metagene position (MetaPlotR scale)";
        kept_bars(
            canvas,
            (m.all, m.kept),
            bbox,
            ticks,
            marks,
            title,
            x_title,
            self.weight.unit(),
        );
    }

    pub(super) fn render_meta(&mut self, frame: &mut Frame, area: Rect) {
        let mut image = std::mem::take(&mut self.meta_plot);
        self.draw_meta_panel(frame, area, &mut image);
        self.meta_plot = image;
    }

    /// The metagene panel; `image` is taken out of the picker because the
    /// counts drawn into it borrow the picker.
    fn draw_meta_panel(&self, frame: &mut Frame, area: Rect, image: &mut PlotImage) {
        let block = panel(format!(" {} ", self.meta_title()), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(m) = self.meta_counts() else {
            frame.render_widget(
                Paragraph::new(Line::from(Span::styled(self.meta_status(), DIM))),
                inner,
            );
            return;
        };
        let [stats, plot, key] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(4),
            Constraint::Length(1),
        ])
        .areas(inner);
        frame.render_widget(Paragraph::new(kept_key(self.weight)), key);
        let region_sum = |c: &[usize], r: usize| {
            let start: usize = m.regions[..r].iter().sum();
            c.get(start..start + m.regions[r])
                .map_or(0, |bins| bins.iter().sum::<usize>())
        };
        let dim = |t: String| Span::styled(t, DIM);
        let mut line = vec![dim("kept / all  ".into())];
        for (r, name) in REGION_NAMES[..3].iter().enumerate() {
            line.push(dim(format!("{name} ")));
            line.push(Span::raw(format!(
                "{}/{}   ",
                region_sum(m.kept, r),
                region_sum(m.all, r)
            )));
        }
        line.push(dim(format!("off-transcript {}", m.unassigned)));
        frame.render_widget(Paragraph::new(Line::from(line)), stats);

        let drawn = self.controls.images().is_some_and(|picker| {
            image.render(frame, plot, picker, |w, h| {
                figure::svg(w, h, |c| {
                    self.draw_meta(c, &m, (0.0, 0.0, w, h), String::new())
                })
            })
        });
        if drawn {
            return;
        }
        // Stretch the bins over the whole chart, one column each: HistPlot
        // gives a bin a whole number of columns, which leaves the rest of a
        // panel empty.
        let n = m.all.len().max(1);
        let cols = (plot.width.saturating_sub(GUTTER) as usize).max(n);
        let bin_of = |x: usize| x * n / cols;
        let stretch = |v: &[usize]| (0..cols).map(|x| v[bin_of(x)]).collect::<Vec<_>>();
        let (all, kept) = (stretch(m.all), stretch(m.kept));
        let ticks: Vec<(usize, String)> = meta_ticks(m.regions)
            .into_iter()
            .map(|(b, t)| ((2 * b + 1) * cols / (2 * n), t))
            .collect();
        let label = |k: i32| {
            ticks
                .iter()
                .find(|(i, _)| *i as i32 == k)
                .map(|(_, t)| t.clone())
        };
        HistPlot {
            bins: Binning::with_width(Scale::Linear, 1.0),
            kmin: 0,
            counts: &all,
            style: &|_| PLAIN,
            subset: Some(&kept),
            y_scale: Scale::Linear,
            y_max: None,
            pointer: None,
            marks: Vec::new(),
            x_label: Some(&label),
            tick_every: Some(1),
        }
        .render(frame.buffer_mut(), plot);
    }

    pub(super) fn criteria_lines(&self) -> Vec<Line<'static>> {
        let dim = |t: String| Span::styled(t, DIM);
        let [name_w, value_w, count_w] = TABLE_WIDTHS;
        let header = |first: &str, second: &str, count: &str| {
            Line::from(dim(format!(
                "{first:<name_w$}{second:>value_w$} │{count:>w$}",
                w = count_w - 1
            )))
        };
        let rule = format!(
            "{}┼{}",
            "─".repeat(name_w + value_w + 1),
            "─".repeat(count_w)
        );
        let mut lines = vec![
            header("", "", COUNT_HEADER.0),
            header("  threshold", "value", COUNT_HEADER.1),
            Line::from(dim(rule)),
        ];
        for (j, &c) in self.view().criteria.iter().enumerate() {
            let focused = j == self.focus;
            let (value, value_style) = match c.shown(&self.filter) {
                None => ("off".to_string(), DIM),
                Some(v) => (v, if focused { HIGHLIGHT } else { PLAIN }),
            };
            lines.push(Line::from(vec![
                Span::styled(if focused { "▸ " } else { "  " }, HIGHLIGHT),
                Span::styled(
                    format!("{:<w$}", c.label(), w = name_w - 2),
                    if focused { HIGHLIGHT } else { PLAIN },
                ),
                Span::styled(format!("{value:>value_w$}"), value_style),
                dim(" │".into()),
                Span::styled(
                    format!("{:>w$}", self.tally.fail[c as usize], w = count_w - 1),
                    ACCENTED,
                ),
            ]));
        }
        let n = self.view().table.len();
        lines.push(Line::from(""));
        lines.push(Line::from(vec![
            dim("  kept ".into()),
            Span::styled(format!("{}", self.tally.kept), HIGHLIGHT),
            dim(format!(" / {n} sites ({:.2}%)", pct(self.tally.kept, n))),
        ]));
        lines.push(Line::from(vec![
            dim("  in ".into()),
            Span::raw(format!("{}", self.tally.genes)),
            dim(format!(" / {} genes", self.view().table.n_genes)),
        ]));
        lines.push(Line::from(""));
        lines.extend(LEGEND.iter().map(|l| Line::from(dim(format!("  {l}")))));
        lines
    }

    /// The current view as a figure: the threshold table beside the
    /// focused column's histogram.
    pub(super) fn figure(&self) -> String {
        let c = self.criterion();
        let view = self.view();
        let n = view.table.len();
        let meta = self.meta_counts();
        let gene = self.gene_profile(80);
        let height = 340.0
            + if gene.is_some() { 220.0 } else { 0.0 }
            + if meta.is_some() { 220.0 } else { 0.0 };
        let mut canvas = Canvas::new(720.0, height);
        canvas.bold(16.0, 22.0, &self.title, 12.0, Anchor::Start, INK);
        canvas.text(
            16.0,
            38.0,
            &format!(
                "{}: kept {} of {} sites ({:.2}%) in {} of {} genes",
                view.table.modality,
                self.tally.kept,
                n,
                pct(self.tally.kept, n),
                self.tally.genes,
                view.table.n_genes
            ),
            9.0,
            Anchor::Start,
            MUTED,
        );

        let (x0, mut y) = (16.0, 60.0);
        // Right edges of the value and count columns, and the rule between.
        let cols = [x0 + 150.0, x0 + 230.0];
        let sep = cols[0] + 8.0;
        canvas.text(x0, y + 10.0, "threshold", 8.0, Anchor::Start, MUTED);
        canvas.text(cols[0], y + 10.0, "value", 8.0, Anchor::End, MUTED);
        canvas.text(cols[1], y, COUNT_HEADER.0, 8.0, Anchor::End, MUTED);
        canvas.text(cols[1], y + 10.0, COUNT_HEADER.1, 8.0, Anchor::End, MUTED);
        let top = y - 9.0;
        y += 15.0;
        canvas.line(x0, y, cols[1] + 4.0, y, MUTED, 0.6);
        let rows_end = y + 16.0 * view.criteria.len() as f64 + 5.0;
        canvas.line(sep, top, sep, rows_end, MUTED, 0.6);
        y -= 5.0;
        for &k in &view.criteria {
            y += 16.0;
            let value = k.shown(&self.filter).unwrap_or_else(|| "off".into());
            let colour = if k == c { figure::ACCENT } else { INK };
            canvas.text(x0, y, k.label(), 9.0, Anchor::Start, colour);
            canvas.text(cols[0], y, &value, 9.0, Anchor::End, colour);
            let n = self.tally.fail[k as usize].to_string();
            canvas.text(cols[1], y, &n, 9.0, Anchor::End, INK);
        }
        y += 8.0;
        for line in LEGEND.iter().chain(&FIGURE_KEY) {
            y += 10.0;
            canvas.text(x0, y, line, 7.5, Anchor::Start, MUTED);
        }

        let title = format!("{} · y: sites per bin", c.label());
        self.draw_hist(&mut canvas, (320.0, 50.0, 390.0, 280.0), title);
        let mut below = 340.0;
        if let Some(p) = &gene {
            draw_gene(&mut canvas, p, (320.0, below, 390.0, 210.0));
            below += 220.0;
        }
        if let Some(m) = &meta {
            let title = format!("{} {}", view.table.modality, self.meta_title());
            self.draw_meta(&mut canvas, m, (320.0, below, 390.0, 210.0), title);
        }
        canvas.finish()
    }

    /// The focused column's histogram in the box `(x, y, w, h)`, as the
    /// export and the in-terminal image draw it.
    pub(super) fn draw_hist(
        &self,
        canvas: &mut Canvas,
        (x, y, w, h): (f64, f64, f64, f64),
        title: String,
    ) {
        let c = self.criterion();
        // Histogram ticks: the smallest value in about six bins.
        let col = &self.column;
        let (nb, lowest) = (col.lowest.len(), &col.lowest);
        let every = (nb / 6).max(1);
        let ticks = (0..nb)
            .step_by(every)
            .filter_map(|b| (b..nb).find(|&j| lowest[j].is_finite()))
            .map(|b| (b, c.fmt_short(lowest[b])))
            .collect::<Vec<_>>();
        let key = self.pointer_key();
        let pointer = key.map(|k| (k - col.hist.kmin) as usize);
        let dropped = |b: usize| c.drops(col.hist.kmin + b as i32, key);
        let values: Vec<f64> = col.hist.counts.iter().map(|&v| v as f64).collect();
        let front: Vec<f64> = self.tally.subset.iter().map(|&v| v as f64).collect();
        Bars {
            values: &values,
            front: Some(&front),
            accent: &dropped,
            y_scale: self.y_scale,
            y_max: None,
            ticks,
            pointer,
            marks: Vec::new(),
            title,
            x_title: format!("{} ({} bins)", c.axis(), self.scale().name()),
            y_title: "sites".into(),
        }
        .draw(canvas, x, y, w, h);
    }

    pub(super) fn render_hist(&mut self, frame: &mut Frame, area: Rect) {
        let c = self.criterion();
        let block = panel(format!(" {} · y: sites per bin ", c.axis()), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let front: usize = self.tally.subset.iter().sum();
        let n = self.view().table.len();
        let key_text = key_line(&[
            (DIM, format!("all sites {n}")),
            (PLAIN, format!("pass the other thresholds {front}")),
            (
                ACCENTED,
                format!("filtered out only by this {}", self.tally.only),
            ),
        ]);
        let key_rows = if key_text.width() > inner.width as usize {
            2
        } else {
            1
        };
        let [stats, plot, key] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(5),
            Constraint::Length(key_rows),
        ])
        .areas(inner);

        let col = &self.column;
        let dim = |t: String| Span::styled(t, DIM);
        let mut first = vec![
            dim("min ".into()),
            Span::raw(c.fmt_short(col.lo)),
            dim("   max ".into()),
            Span::raw(c.fmt_short(col.hi)),
        ];
        if col.n_inf > 0 {
            first.push(dim(format!(
                "   {} infinite, drawn at the end bin",
                col.n_inf
            )));
        }
        frame.render_widget(Paragraph::new(Line::from(first)), stats);
        frame.render_widget(
            Paragraph::new(key_text).wrap(ratatui::widgets::Wrap { trim: true }),
            key,
        );

        let mut image = std::mem::take(&mut self.plot);
        let drawn = self.controls.images().is_some_and(|picker| {
            image.render(frame, plot, picker, |w, h| {
                figure::svg(w, h, |c| self.draw_hist(c, (0.0, 0.0, w, h), String::new()))
            })
        });
        self.plot = image;
        if drawn {
            return;
        }
        let col = &self.column;
        let pointer = self.pointer_key();
        let style = |k: i32| if c.drops(k, pointer) { ACCENTED } else { PLAIN };
        HistPlot {
            bins: col.hist.bins,
            kmin: col.hist.kmin,
            counts: &col.hist.counts,
            style: &style,
            subset: Some(&self.tally.subset),
            y_scale: self.y_scale,
            y_max: None,
            pointer,
            marks: Vec::new(),
            x_label: None,
            tick_every: None,
        }
        .render(frame.buffer_mut(), plot);
    }
}
