//! The frame: the layout of the screen, and its footer.

use super::*;

impl SitePicker<'_> {
    /// The footer: a prompt, a hint, or the keys of the mode.
    pub(super) fn footer_line(&self) -> Line<'static> {
        if let Some(keys) = self.footer_pairs() {
            return help_line(&keys);
        }
        if let Some(said) = self.said() {
            return figure::status_line(said);
        }
        if let Some(line) = self.controls.footer() {
            return line;
        }
        match &self.mode {
            Mode::Edit(buf) => input_line(
                &format!("{} {}: ", self.criterion().label(), self.criterion().flag()),
                buf,
                &[("Enter", "set"), ("Esc", "back")],
            ),
            Mode::Title { y, buf, .. } => {
                let (x_axis, y_axis) = self.titled();
                let (what, next) = match (y, y_axis) {
                    (false, _) => (format!("x title of {}", x_axis.short()), "next"),
                    (true, Some(a)) => (format!("y title of {}", a.short()), "set"),
                    (true, None) => ("y title of the counts".into(), "set"),
                };
                input_line(
                    &format!("{what} (empty: its own): "),
                    buf,
                    &[("Enter", next), ("Esc", "back")],
                )
            }
            Mode::Writing => help_line(&[("", "writing; the view closes when done")]),
            Mode::Find => input_line(
                "gene: ",
                &self.list.find,
                &[("↑/↓", "move"), ("Enter", "focus"), ("Esc", "cancel")],
            ),
            // Browsing, past the Enter that did nothing.
            _ => Line::styled(format!(" {}", enter_hint("applies")), DIM),
        }
    }

    /// A report's outcome, which shows while browsing, past a save's
    /// outcome but not its prompt, until the next key.
    fn said(&self) -> Option<&String> {
        self.status
            .as_ref()
            .filter(|_| self.browsing() && !self.enter_hint)
    }

    /// The keys the footer names, where it names keys: none while a
    /// report's or a save's outcome, a prompt or the Enter hint shows.
    pub(super) fn footer_pairs(&self) -> Option<Vec<(&'static str, &'static str)>> {
        let hint = self.enter_hint && matches!(self.mode, Mode::Browse);
        if self.said().is_some() || self.controls.footer().is_some() || hint {
            return None;
        }
        self.mode_keys()
    }

    /// The keys of a mode that is no prompt, as the footer names them.
    pub(super) fn mode_keys(&self) -> Option<Vec<(&'static str, &'static str)>> {
        Some(match &self.mode {
            Mode::Confirm => vec![(APPLY_KEYS, "apply and write"), ("Esc/n", "back")],
            Mode::Leave(picked) => vec![
                match picked {
                    Picked::PrintOnly(_) => ("p", "print the flags and leave"),
                    _ => ("q", "leave without writing"),
                },
                ("Esc/n", "back"),
            ],
            Mode::Gff(_) => vec![
                ("type", "find"),
                ("↑/↓", "move"),
                ("Enter", "open / read"),
                ("←", "up"),
                ("Esc", "clear find / back"),
            ],
            Mode::Browse => self.browse_footer(),
            Mode::Edit(_) | Mode::Title { .. } | Mode::Find | Mode::Writing => return None,
        })
    }

    /// The footer's keys while browsing: those that act on the panel
    /// focused, as [`SitePicker::browse_help`] has them, put short.
    fn browse_footer(&self) -> Vec<(&'static str, &'static str)> {
        let mut keys = vec![("Tab/Shift+Tab", "panel")];
        match self.panel {
            Panel::Thresholds => keys.extend([("↑/↓", "knob"), ("←/→", "bin")]),
            Panel::Genes => keys.extend([("↑/↓", "gene"), ("PgUp/PgDn", "page"), ("/", "find")]),
            Panel::Pairs => {
                // The arrows step from panel to panel, where there are two.
                if self.view().pairs.len() > 1 {
                    keys.push(("arrows", "panel"));
                }
                keys.push(("x/y X/Y", "axis"));
            }
        }
        // What matters most first.
        match self.read_only() {
            true => keys.push(("q", "quit")),
            false => keys.extend([(APPLY_KEYS, "apply"), ("p", "print flags"), ("q", "cancel")]),
        }
        self.controls.help_keys(&mut keys);
        if self.views.len() > 1 {
            keys.push(("m/M", "modality"));
        }
        keys.extend([
            ("[ ]", "gene"),
            ("-/+", "±1"),
            ("0-9/e", "type"),
            ("o", "off"),
            ("r", "reset"),
        ]);
        // On the panels `x` and `y` step the axes instead; `x` steps the
        // knob's scale only where it has more than one.
        if self.panel != Panel::Pairs {
            keys.push(match self.criterion().scales().len() > 1 {
                true => ("x/y", "scale"),
                false => ("y", "scale"),
            });
        }
        keys.push(("v", if self.two_d { "histograms" } else { "2D" }));
        keys.extend([("f", "front"), ("u", "p axis"), ("t", "titles")]);
        // In 2D the weight's plots are hidden: `w` changes nothing seen.
        match self.two_d {
            true => keys.push(("A", "auto pairs")),
            false => keys.push(("w", "weight")),
        }
        keys.extend([("a", "annotation"), (REPORT_KEY, "report"), ("?", "keys")]);
        keys
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
            "x {} · y {} · {} in front · p axis {}{}",
            self.scale().name(),
            self.y_scale.name(),
            self.arm_look.front.name(),
            self.arm_look.p_axis.name(),
            self.controls.tag()
        );
        let badge = if self.read_only() { "qc-view" } else { "qc" };
        frame.render_widget(header(badge, &self.title, &scales), top);

        self.hits.clear();
        let names: Vec<&str> = self.views.iter().map(|v| &*v.table.modality).collect();
        let button = (matches!(self.mode, Mode::Browse) && !self.read_only())
            .then_some(("✓ apply", Hit::Apply));
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
        if self.two_d {
            self.render_pairs(frame, right);
        } else {
            let [hist, gene, meta] = Layout::vertical([
                Constraint::Fill(3),
                Constraint::Fill(2),
                Constraint::Fill(2),
            ])
            .areas(right);
            self.render_hist(frame, hist);
            self.render_gene(frame, gene);
            self.render_meta(frame, meta);
        }

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
