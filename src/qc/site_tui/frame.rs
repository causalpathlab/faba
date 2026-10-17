//! The frame: the layout of the screen, and its footer.

use super::*;

impl SitePicker<'_> {
    /// The footer: a prompt, a hint, or the keys of the mode.
    pub(super) fn footer_line(&self) -> Line<'static> {
        // A report's outcome shows while browsing, past a save's outcome
        // but not its prompt, until the next key.
        if let Some(said) = self
            .status
            .as_ref()
            .filter(|_| self.browsing() && !self.enter_hint)
        {
            return figure::status_line(said);
        }
        match (&self.mode, self.controls.footer()) {
            (_, Some(line)) => line,
            (Mode::Edit(buf), None) => input_line(
                &format!("{} {}: ", self.criterion().label(), self.criterion().flag()),
                buf,
                &[("Enter", "set"), ("Esc", "back")],
            ),
            (Mode::Confirm, None) => {
                help_line(&[(APPLY_KEYS, "apply and write"), ("Esc/n", "back")])
            }
            (Mode::Writing, None) => help_line(&[("", "writing; the view closes when done")]),
            (Mode::Leave(picked), None) => help_line(&[
                match picked {
                    Picked::PrintOnly(_) => ("p", "print the flags and leave"),
                    _ => ("q", "leave without writing"),
                },
                ("Esc/n", "back"),
            ]),
            (Mode::Gff(_), None) => help_line(&[
                ("type", "find"),
                ("↑/↓", "move"),
                ("Enter", "open / read"),
                ("←", "up"),
                ("Esc", "clear find / back"),
            ]),
            (Mode::Find, None) => input_line(
                "gene: ",
                &self.list.find,
                &[("↑/↓", "move"), ("Enter", "focus"), ("Esc", "cancel")],
            ),
            (Mode::Browse, None) => {
                let apply = (APPLY_KEYS, "apply");
                if self.enter_hint {
                    return Line::styled(format!(" {}", enter_hint("applies")), DIM);
                }
                let mut keys = vec![("Tab", "panel")];
                keys.extend(match self.panel {
                    Panel::Thresholds => [("↑/↓", "knob"), ("←/→", "bin")],
                    Panel::Genes => [("↑/↓", "gene"), ("/", "find")],
                });
                // What matters most first.
                keys.extend([apply, ("p", "print flags"), ("q", "cancel")]);
                self.controls.help_keys(&mut keys);
                keys.extend([
                    ("m/M", "modality"),
                    ("[ ]", "gene"),
                    ("-/+", "±1"),
                    ("0-9/e", "type"),
                    ("o", "off"),
                    ("r", "reset"),
                    ("x/y", "scale"),
                    ("a", "annotation"),
                    (REPORT_KEY, "report"),
                ]);
                help_line(&keys)
            }
        }
    }

    /// Draw the whole screen: header, tabs, panels, footer, and the
    /// pop-up of the mode over them.
    pub(super) fn render_frame(&mut self, frame: &mut Frame) {
        self.drawn = frame.area();
        // The footer takes as many lines as its keys need at this width.
        let help = crate::tui::wrap_help(self.footer_line(), frame.area().width);
        let [top, tabs, body, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(help.len() as u16),
        ])
        .areas(frame.area());

        let scales = format!(
            "x {} · y {}{}",
            self.scale().name(),
            self.y_scale.name(),
            self.controls.tag()
        );
        frame.render_widget(header("qc", &self.title, &scales), top);

        self.hits.clear();
        let names: Vec<&str> = self.views.iter().map(|v| &*v.table.modality).collect();
        let button = matches!(self.mode, Mode::Browse).then_some(("✓ apply", Hit::Apply));
        let tabs_line = tab_bar(
            tabs,
            &names,
            self.modality,
            button,
            &self.hits,
            Hit::Modality,
        );
        // Keep the mark clear of the header's text and the modality tabs.
        let reserve = (self.title.chars().count() + scales.chars().count() + 12)
            .max(tabs_line.width()) as u16;
        frame.render_widget(tabs_line, tabs);
        let corner = Rect::new(top.x, top.y, top.width, 2);
        crate::figure::logo::draw_mini_logo(frame.buffer_mut(), corner, reserve);

        let [left, right] =
            Layout::horizontal([Constraint::Length(56), Constraint::Fill(1)]).areas(body);
        let table = self.criteria_lines();
        let [table_area, genes] = Layout::vertical([
            Constraint::Length(table.len() as u16 + 2),
            Constraint::Fill(1),
        ])
        .areas(left);
        let block = panel(" site thresholds ".into(), self.panel == Panel::Thresholds);
        let inner = block.inner(table_area);
        frame.render_widget(block, table_area);
        frame.render_widget(Paragraph::new(table), inner);
        self.hits.add(table_area, Hit::Thresholds(None));
        // The knobs come under the table's three header lines.
        let knobs = Rect::new(
            inner.x,
            inner.y + 3,
            inner.width,
            inner.height.saturating_sub(3),
        );
        let n = self.view().criteria.len();
        self.hits.rows(knobs, 0, n, |j| Hit::Thresholds(Some(j)));
        self.hits.add(genes, Hit::Genes(None));
        self.render_gene_list(frame, genes);
        let [hist, gene, meta] = Layout::vertical([
            Constraint::Fill(3),
            Constraint::Fill(2),
            Constraint::Fill(2),
        ])
        .areas(right);
        self.render_hist(frame, hist);
        self.render_gene(frame, gene);
        self.render_meta(frame, meta);

        frame.render_widget(Paragraph::new(help), footer);
        match self.mode {
            Mode::Confirm => self.render_confirm(frame, body),
            Mode::Writing => self.render_writing(frame, body),
            Mode::Gff(ref browser) => render_gff(frame, body, browser),
            Mode::Leave(ref picked) => render_leave(frame, body, picked),
            _ => {}
        }
        self.controls.render_save(frame, frame.area());
    }
}
