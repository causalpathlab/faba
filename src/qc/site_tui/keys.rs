//! Keys: the report's chord, the key that leaves, and what each mode
//! does with a key ([`Screen::handle_key`] hands it on by mode).

use super::*;

/// The key that writes the HTML report, as the footer names it.
pub(super) const REPORT_KEY: &str = "Ctrl+E";

/// Whether `key` is the report's: Ctrl+E and no other modifier, but the
/// Shift a terminal may add for a capital E.
pub(super) fn is_report(key: &KeyEvent) -> bool {
    key.modifiers.difference(KeyModifiers::SHIFT) == KeyModifiers::CONTROL
        && matches!(key.code, KeyCode::Char('e' | 'E'))
}

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
            KeyCode::Esc | KeyCode::Char('n' | 'q') => self.mode = Mode::Browse,
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

    /// A key while browsing.
    pub(super) fn browse_key(&mut self, key: KeyEvent) {
        match key.code {
            _ if is_apply(&key) => self.mode = Mode::Confirm,
            KeyCode::Enter => self.enter_hint = true,
            KeyCode::Tab | KeyCode::BackTab => {
                self.panel = match self.panel {
                    Panel::Thresholds => Panel::Genes,
                    Panel::Genes => Panel::Thresholds,
                }
            }
            KeyCode::Char('m') => self.switch_modality(1),
            KeyCode::Char('M') => self.switch_modality(-1),
            KeyCode::Char('[') => self.step_gene(-1),
            KeyCode::Char(']') => self.step_gene(1),
            KeyCode::Char('/') => self.start_find(),
            KeyCode::Char('c') => self.switch_weight(),
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
            KeyCode::Char('p') => self.leave(Picked::PrintOnly(self.filter.clone())),
            KeyCode::Char('q') => self.leave(Picked::Cancelled),
            _ => {}
        }
    }
}
