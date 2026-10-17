//! Drawing a frame: the rows fitted to the screen, each as a picture or
//! text, the readout, and the footer's keys.

use super::*;

impl PileupView<'_> {
    /// Lines the genes row takes: a line per lane, and its frame.
    pub(super) fn genes_height(&self) -> u16 {
        self.lanes() as u16 + 2
    }

    /// Fit `all` the rows to `area` less `footer` lines of keys, as a frame
    /// there draws them, the genes row (if one of them) `genes` lines.
    pub(super) fn fit_in(&mut self, all: &[Row], area: Rect, footer: usize, genes: u16) {
        let genes = if all.contains(&Row::Genes) { genes } else { 0 };
        self.fit_rows(all, area.height.saturating_sub(2 + footer as u16 + genes));
    }

    /// Draw in `area`: the header, the rows and the readout, and with
    /// `footer` the keys and the save pop-up (a stack draws those once).
    pub(super) fn render_in(&mut self, frame: &mut Frame, area: Rect, footer: bool) {
        self.take_genes();
        let all = self.rows();
        let genes_height = self.genes_height();
        // The fit leaves the track scrolled to as it was, so only the last
        // one counts.
        let help = if footer {
            fitted_help(self, area.width, Self::footer_line, |v, lines| {
                v.fit_in(&all, area, lines, genes_height);
                v.rows_hidden()
            })
        } else {
            Vec::new()
        };
        self.fit_in(&all, area, help.len(), genes_height);
        let rows = self.shown_rows(&all);
        let plot_rows = Self::plot_count(&rows);
        let n_plots = plot_rows as u32;
        let mut layout = vec![Constraint::Length(1)];
        layout.extend(rows.iter().map(|&r| match r {
            Row::Genes => Constraint::Length(genes_height),
            _ => Constraint::Ratio(1, n_plots),
        }));
        layout.extend([Constraint::Length(1), Constraint::Length(help.len() as u16)]);
        let areas = Layout::vertical(layout).split(area);
        let (top, readout, footer_area) =
            (areas[0], areas[areas.len() - 2], areas[areas.len() - 1]);

        let inner_width = areas[1].width.saturating_sub(2 + GUTTER);
        self.columns = (inner_width as usize).max(1);
        // Inside the panel's border, past the y axis.
        (self.drawn, self.chart_x) = (area, areas[1].x + 1 + GUTTER);
        self.chart_rows = (areas[1].y, areas[plot_rows].bottom());
        self.clamp_window();
        let (lo, hi) = self.window;
        let extra = format!(
            "{}:{}-{}  {} bp/bar · {} · y {}{}{}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.bin_width(),
            self.show.label(self.on, self.off),
            self.y_scale.name(),
            match self.y_lock {
                Some(max) => format!(" ≤ {max:.0} locked"),
                None => String::new(),
            },
            self.controls.tag()
        );
        frame.render_widget(header(self.badge, &self.title, &extra), top);

        let every = TICK_SPACING.max(self.columns / 5);
        let ticks: Vec<Option<String>> = (0..self.columns)
            .map(|k| {
                (k.is_multiple_of(every) && k + 10 < self.columns)
                    .then(|| fmt_thousands(self.bar_range(k).0))
            })
            .collect();
        let bins = self.bins();
        let pointer = self.cursor_col();
        let mut plots = std::mem::take(&mut self.plots);
        plots.resize_with(rows.len(), PlotImage::default);
        for (k, (&row, &area)) in rows.iter().zip(&areas[1..areas.len() - 2]).enumerate() {
            let title = self.marked_title(k, &rows, self.row_title(row));
            let block = panel(format!(" {title} "), true);
            let inner = block.inner(area);
            frame.render_widget(block, area);
            let drawn = self.controls.images().is_some_and(|picker| {
                plots[k].render(frame, inner, picker, |w, h| {
                    let bbox = (0.0, 0.0, w, h);
                    crate::figure::svg(w, h, |c| {
                        self.draw_row(row, &bins, c, bbox, Some(pointer), String::new())
                    })
                })
            });
            if drawn {
                continue;
            }
            self.render_row_text(frame.buffer_mut(), row, inner, &bins, pointer, &ticks);
        }
        self.plots = plots;

        frame.render_widget(Paragraph::new(self.readout(&bins)), readout);
        if footer {
            frame.render_widget(Paragraph::new(help), footer_area);
            self.controls.render_save(frame, frame.area());
        }
    }

    /// The genes row as text, a lane per line: heavy exons on a thin arrowed intron line.
    pub(super) fn render_genes(
        &self,
        buf: &mut ratatui::buffer::Buffer,
        area: ratatui::layout::Rect,
    ) {
        let chart_x = area.x + GUTTER;
        if area.width <= GUTTER || area.height == 0 {
            return;
        }
        let right = area.right();
        let edges = self.edges();
        let col = |pos: i64| (chart_x + edges.col_of(pos) as u16).min(right - 1);
        for (lane, g) in self.gene_lanes() {
            let y = area.y + lane as u16;
            if y >= area.bottom() {
                break;
            }
            let (c0, c1) = draw_gene_line(buf, y, col, g, self.window, PLAIN);
            // The symbol after the gene if it fits, else before it.
            let label = g.symbol.as_ref();
            let len = label.chars().count() as u16;
            let lx = if c1 + 2 + len <= right {
                c1 + 2
            } else {
                c0.saturating_sub(len + 1).max(area.x)
            };
            buf.set_string(lx, y, label, DIM);
        }
    }

    /// The footer: the search being typed, a status, the save prompt, or
    /// the keys, `[/]` among them when rows are `hidden`.
    pub(super) fn footer_line(&self, hidden: bool) -> Line<'static> {
        if let Some(buf) = self.search.text() {
            input_line(
                "search gene or chr:start-end: ",
                buf,
                &[("Enter", "go"), ("Esc", "back")],
            )
        } else if let Some(msg) = &self.status {
            status_line(msg)
        } else {
            self.controls.footer().unwrap_or_else(|| self.help(hidden))
        }
    }

    pub(super) fn help(&self, hidden: bool) -> Line<'static> {
        let mut keys = vec![
            ("←/→", "bar"),
            ("PgUp/PgDn", "half a view"),
            ("n/N", "next/prev site"),
            ("+/-", "zoom"),
            ("0", "whole"),
            ("/", "search"),
            ("y", "scale"),
            (
                "Y",
                if self.y_lock.is_some() {
                    "unlock y"
                } else {
                    "lock y"
                },
            ),
            ("g/Esc", "gene list"),
        ];
        if self.tracks.iter().any(|t| t.is_reads()) {
            keys.push(("c", self.show.next().short()));
        }
        if self.contrast.is_some() {
            keys.push(("d", "difference/fold"));
        }
        if self.mirrorable() {
            keys.push(("v", if self.mirror { "split" } else { "mirror" }));
        }
        if hidden {
            keys.push(("[/]", "tracks"));
        }
        if self.switchable {
            keys.push(("m", "m6A/A-to-I"));
        }
        self.controls.help_keys(&mut keys);
        keys.push(("q", "quit"));
        help_line(&keys)
    }

    /// Row `row` as text in `inner`, the cursor at column `pointer` and
    /// the x axis labelled `ticks` (a label or none per column).
    fn render_row_text(
        &self,
        buf: &mut ratatui::buffer::Buffer,
        row: Row,
        inner: Rect,
        bins: &Bins,
        pointer: usize,
        ticks: &[Option<String>],
    ) {
        let (lo, hi) = self.window;
        let column_label = |k: usize| ticks.get(k).cloned().flatten();
        let label = |k: i32| column_label(k as usize);
        match row {
            Row::Track(i) => {
                let marks = self.tracks[i].sites_in(lo, hi).into_iter();
                let shown = &bins.shown[i];
                HistPlot {
                    bins: Binning::with_width(Scale::Linear, 1.0),
                    kmin: 0,
                    counts: &shown.values,
                    style: &|_| {
                        if shown.accent {
                            ACCENTED
                        } else {
                            PLAIN
                        }
                    },
                    subset: shown.front.as_deref(),
                    y_scale: self.y_scale,
                    y_max: self.row_max(row, bins),
                    pointer: Some(pointer as i32),
                    marks: marks
                        .map(|s| (bins.edges.col_of(s) as i32, "+", DIM))
                        .collect(),
                    x_label: Some(&label),
                    tick_every: Some(1),
                }
                .render(buf, inner);
            }
            Row::Contrast(contrast) => {
                let values = self.contrast_values(contrast, bins);
                let (up, down) = crate::figure::split_signed(&values);
                let max = up.iter().chain(&down).fold(0.0, |m: f64, &v| m.max(v));
                let side = |counts, style| MirrorSide {
                    counts,
                    subset: None,
                    style,
                    name: "",
                };
                let measure = contrast.measure;
                MirrorPlot {
                    up: side(&up, ACCENTED),
                    down: side(&down, PLAIN),
                    y_scale: Scale::Linear,
                    y_max: Some(max),
                    y_labels: Some([max, 0.0, -max].map(|v| measure.label(v))),
                    pointer: Some(pointer),
                    x_label: Some(&column_label),
                }
                .render(buf, inner);
            }
            Row::Mirror(a, b) => {
                let side = |k: usize| {
                    let shown = &bins.shown[k];
                    MirrorSide {
                        counts: &shown.values,
                        subset: shown.front.as_deref(),
                        style: if shown.accent { ACCENTED } else { PLAIN },
                        name: self.tracks[k].label,
                    }
                };
                MirrorPlot {
                    up: side(a),
                    down: side(b),
                    y_scale: self.y_scale,
                    y_max: self.row_max(row, bins),
                    y_labels: None,
                    pointer: Some(pointer),
                    x_label: Some(&column_label),
                }
                .render(buf, inner);
            }
            Row::Genes => self.render_genes(buf, inner),
        }
    }

    /// The line under the rows: where the cursor is, and each row's value
    /// in its bar.
    pub(super) fn readout(&self, bins: &Bins) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let col = self.cursor_col();
        let (start, stop) = self.bar_range(col);
        let mut spans = vec![
            Span::styled(
                format!(" {}:{}", self.chr, fmt_thousands(self.cursor)),
                HIGHLIGHT,
            ),
            dim(format!(
                "   bar {}-{}",
                fmt_thousands(start),
                fmt_thousands(stop)
            )),
        ];
        if let Some(c) = self.contrast {
            let (a, b) = (self.tracks[c.a].label, self.tracks[c.b].label);
            let value = self.contrast_values(c, bins).get(col).copied().flatten();
            spans.push(dim(format!("   {a} vs {b} ")));
            spans.push(Span::styled(
                value.map_or("-".into(), |v| c.measure.label(v)),
                HIGHLIGHT,
            ));
        }
        for (t, shown) in self.tracks.iter().zip(&bins.shown) {
            spans.push(dim(format!("   {} ", t.label)));
            let v = shown.values.get(col).copied().unwrap_or(0.0);
            spans.push(Span::raw(
                match shown.front.as_ref().and_then(|f| f.get(col)) {
                    Some(f) => format!("{f:.0}/{v:.0}"),
                    None => compact(v),
                },
            ));
        }
        let here = self
            .sites
            .iter()
            .filter(|&&s| s >= start && s < stop.max(start + 1))
            .count();
        spans.push(dim(format!("   {here} site(s) in bar")));
        Line::from(spans)
    }
}

/// The footer `line` gives for `view`, wrapped to `width`, naming `[/]`
/// (`line(view, true)`) only when `hidden`, fitting the rows over that
/// many lines of keys, says rows are left out.
pub(super) fn fitted_help<T>(
    view: &mut T,
    width: u16,
    line: impl Fn(&T, bool) -> Line<'static>,
    mut hidden: impl FnMut(&mut T, usize) -> bool,
) -> Vec<Line<'static>> {
    let help = crate::tui::wrap_help(line(view, false), width);
    if hidden(view, help.len()) {
        crate::tui::wrap_help(line(view, true), width)
    } else {
        help
    }
}
