//! Keys: the report's chord, the key that leaves, and what each mode
//! does with a key ([`Screen::handle_key`] hands it on by mode).

use super::*;

/// The report's chord, the one every view writes its page with.
pub(super) use crate::tui::{is_report, REPORT_KEY};

/// The key that asks to leave as `picked` says, and confirms it.
pub(super) fn leave_key(picked: &Picked) -> char {
    match picked {
        Picked::PrintOnly(_) => 'p',
        _ => 'q',
    }
}

impl SitePicker<'_> {
    /// Browsing, with no save prompt open: where the report's chord acts
    /// and its outcome shows.
    pub(super) fn browsing(&self) -> bool {
        matches!(self.mode, Mode::Browse) && !self.controls.prompting()
    }

    /// A key while typing a threshold.
    pub(super) fn edit_key(&mut self, key: KeyEvent) {
        let Mode::Edit(buf) = &mut self.mode else {
            return;
        };
        match key.code {
            // The apply key sets the value, as Enter does.
            _ if is_apply(&key) || key.code == KeyCode::Enter => {
                let typed = buf.parse::<f64>().ok();
                self.mode = Mode::Browse;
                if let Some(v) = typed {
                    self.set_threshold(v);
                }
            }
            KeyCode::Char(ch)
                if (ch.is_ascii_digit() || "-+.eEinf".contains(ch)) && buf.len() < 24 =>
            {
                buf.push(ch)
            }
            KeyCode::Backspace => {
                buf.pop();
            }
            KeyCode::Esc => self.mode = Mode::Browse,
            _ => {}
        }
    }

    /// A key while typing into the gene list's filter. Letters type, so
    /// only the arrow keys move through the matches.
    pub(super) fn find_key(&mut self, key: KeyEvent) {
        match key.code {
            _ if is_apply(&key) => self.end_find(true),
            KeyCode::Char(ch) if self.list.find.len() < 32 => self.set_find(|f| f.push(ch)),
            KeyCode::Backspace => self.set_find(|f| {
                f.pop();
            }),
            KeyCode::Up => self.step_gene(-1),
            KeyCode::Down => self.step_gene(1),
            KeyCode::PageUp => self.step_gene(-10),
            KeyCode::PageDown => self.step_gene(10),
            KeyCode::Enter => self.end_find(true),
            KeyCode::Esc => self.end_find(false),
            _ => {}
        }
    }

    /// A key while asking before the thresholds are applied.
    pub(super) fn confirm_key(&mut self, key: KeyEvent) {
        match key.code {
            _ if is_apply(&key) => self.apply(),
            // `q` is quit everywhere, so never back here.
            KeyCode::Esc | KeyCode::Char('n') => self.mode = Mode::Browse,
            _ => {}
        }
    }

    /// A key while asking before leaving without writing.
    pub(super) fn leave_key_press(&mut self, key: KeyEvent) {
        let Mode::Leave(picked) = &self.mode else {
            return;
        };
        match key.code {
            KeyCode::Char(ch) if ch == leave_key(picked) => {
                let picked = picked.clone();
                self.decide(picked);
            }
            KeyCode::Esc | KeyCode::Char('n') => self.mode = Mode::Browse,
            _ => {}
        }
    }

    /// Every key while browsing, for the `?` pop-up.
    pub(super) fn browse_help(&self) -> Vec<(&'static str, &'static str)> {
        let mut keys = vec![
            ("Tab/Shift+Tab", "next/previous panel"),
            ("↑/↓ k/j", "knob; gene in the genes; panel in 2D"),
            ("←/→ h/l", "bin; panel in 2D"),
            ("PgUp/PgDn", "ten genes up/down, in the genes"),
            ("-/+ , .", "threshold down/up a step (= as +)"),
            ("0-9 e", "type a threshold"),
            ("o", "threshold off"),
            ("r", "reset the threshold"),
            ("x/y", "x/y scale; in 2D the panel's axes"),
            ("X/Y", "in 2D, the panel's axes back"),
            ("v", "2D panels or histograms"),
        ];
        if self.two_d {
            keys.push(("A", "auto pairs"));
        }
        keys.extend([
            ("m/M", "next/previous modality"),
            ("[ ]", "previous/next gene"),
            ("/", "find a gene"),
            ("w", "weight: coverage, converted reads or sites"),
            ("a", "read an annotation"),
        ]);
        self.controls.help_keys(&mut keys);
        keys.push((REPORT_KEY, "write the HTML report"));
        // Read only, nothing is applied and nothing is lost on leaving.
        match self.read_only() {
            true => keys.push(("q", "quit")),
            false => keys.extend([
                (APPLY_KEYS, "apply and write"),
                ("p", "print the flags and leave"),
                ("q", "leave (asks when changed)"),
            ]),
        }
        keys.push(("?", "these keys"));
        keys
    }

    /// A key while browsing.
    pub(super) fn browse_key(&mut self, key: KeyEvent) {
        let writes = !self.read_only();
        match key.code {
            _ if is_apply(&key) && writes => self.mode = Mode::Confirm,
            KeyCode::Enter => self.enter_hint = writes,
            KeyCode::Tab => {
                self.panel = match (self.panel, self.two_d) {
                    (Panel::Thresholds, _) => Panel::Genes,
                    (Panel::Genes, true) => Panel::Pairs,
                    _ => Panel::Thresholds,
                }
            }
            KeyCode::BackTab => {
                self.panel = match (self.panel, self.two_d) {
                    (Panel::Thresholds, true) => Panel::Pairs,
                    (Panel::Thresholds, false) | (Panel::Pairs, _) => Panel::Genes,
                    (Panel::Genes, _) => Panel::Thresholds,
                }
            }
            // `d` as it was before `v` cycled every view's look.
            KeyCode::Char('v' | 'd') => self.toggle_two_d(),
            KeyCode::Char('A') if self.two_d => self.auto_pairs(),
            KeyCode::Left | KeyCode::Char('h') if self.panel == Panel::Pairs => self.step_pair(-1),
            KeyCode::Right | KeyCode::Char('l') if self.panel == Panel::Pairs => self.step_pair(1),
            KeyCode::Up | KeyCode::Char('k') if self.panel == Panel::Pairs => self.step_pair(-2),
            KeyCode::Down | KeyCode::Char('j') if self.panel == Panel::Pairs => self.step_pair(2),
            KeyCode::Char('x') if self.panel == Panel::Pairs => self.step_axis(true, 1),
            KeyCode::Char('X') if self.panel == Panel::Pairs => self.step_axis(true, -1),
            KeyCode::Char('y') if self.panel == Panel::Pairs => self.step_axis(false, 1),
            KeyCode::Char('Y') if self.panel == Panel::Pairs => self.step_axis(false, -1),
            KeyCode::Char('m') => self.switch_modality(1),
            KeyCode::Char('M') => self.switch_modality(-1),
            KeyCode::Char('[') => self.step_gene(-1),
            KeyCode::Char(']') => self.step_gene(1),
            KeyCode::Char('/') => self.start_find(),
            KeyCode::Char('w') => self.switch_weight(),
            KeyCode::Char('a') => self.open_gff_browser(),
            KeyCode::Up | KeyCode::Char('k') if self.panel == Panel::Genes => self.step_gene(-1),
            KeyCode::Down | KeyCode::Char('j') if self.panel == Panel::Genes => self.step_gene(1),
            KeyCode::PageUp if self.panel == Panel::Genes => self.step_gene(-10),
            KeyCode::PageDown if self.panel == Panel::Genes => self.step_gene(10),
            KeyCode::Up | KeyCode::Char('k') => self.move_focus(-1),
            KeyCode::Down | KeyCode::Char('j') => self.move_focus(1),
            KeyCode::Left | KeyCode::Char('h') => self.step_bin(-1),
            KeyCode::Right | KeyCode::Char('l') => self.step_bin(1),
            KeyCode::Char('-' | ',') => self.nudge(-1),
            KeyCode::Char('+' | '=' | '.') => self.nudge(1),
            KeyCode::Char('o') => self.set_off(),
            KeyCode::Char('r') => self.reset(),
            KeyCode::Char('x') => self.cycle_x_scale(),
            KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
            KeyCode::Char('e') => self.mode = Mode::Edit(String::new()),
            KeyCode::Char(ch) if ch.is_ascii_digit() => self.mode = Mode::Edit(ch.to_string()),
            // Esc backs out of every pop-up, so it never leaves the view.
            // Read only, there are no flags to print and nothing to lose.
            KeyCode::Char('p') if writes => self.leave(Picked::PrintOnly(self.filter.clone())),
            KeyCode::Char('q') if writes => self.leave(Picked::Cancelled),
            KeyCode::Char('q') => self.decide(Picked::Cancelled),
            _ => {}
        }
    }
}
