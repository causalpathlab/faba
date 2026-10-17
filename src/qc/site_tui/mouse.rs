//! The mouse: what a click lands on, and what it does there.

use super::*;

/// What a click lands on, as the last draw laid it out.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(super) enum Hit {
    Modality(usize),
    /// The thresholds panel, or one of its knobs.
    Thresholds(Option<usize>),
    /// The gene list, or one of its genes (by place in the list).
    Genes(Option<usize>),
    /// The buttons: apply (as the apply key), and in the confirmation
    /// write and back.
    Apply,
    Write,
    Back,
}

impl SitePicker<'_> {
    /// A click focuses and selects what is under it; the wheel then moves
    /// as ↑/↓ would. Only while browsing, with no prompt open.
    pub(super) fn click(&mut self, m: MouseEvent) {
        let Some(hit) = self.hits.at(m.column, m.row) else {
            return;
        };
        let clicked = matches!(m.kind, MouseEventKind::Down(_));
        // A button does what its key does, the same way.
        let browsing = matches!(self.mode, Mode::Browse) && self.controls.footer().is_none();
        let confirming = matches!(self.mode, Mode::Confirm);
        let button = match hit {
            _ if !clicked => None,
            Hit::Apply if browsing => Some(apply_key()),
            Hit::Write if confirming => Some(apply_key()),
            Hit::Back if confirming => Some(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE)),
            _ => None,
        };
        if let Some(key) = button {
            return self.handle_key(key);
        }
        if !browsing {
            return;
        }
        self.enter_hint = false;
        // The plots are drawn again only when what they show changes.
        let shown = (self.modality, self.focus, self.list.at);
        match hit {
            Hit::Modality(i) if clicked => {
                self.switch_modality(i as isize - self.modality as isize);
            }
            Hit::Thresholds(j) => {
                self.panel = Panel::Thresholds;
                if let Some(j) = j.filter(|_| clicked) {
                    self.move_focus(j as isize - self.focus as isize);
                }
            }
            Hit::Genes(i) => {
                self.panel = Panel::Genes;
                if let Some(i) = i.filter(|_| clicked) {
                    self.step_gene(i as isize - self.list.at as isize);
                }
            }
            _ => return,
        }
        if shown != (self.modality, self.focus, self.list.at) {
            self.plot.invalidate();
        }
        if let Some(key) = wheel_key(&m) {
            self.handle_key(key);
        }
    }
}
