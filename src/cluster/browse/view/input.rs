//! Keys and the mouse: presses, drags and the wheel over what the last
//! frame drew, and every key the panes take.

use super::{Look, ModulesView, Pane, Sort, Spot, PANES};
use crate::site_analysis::pileup::query_symbol;
use ratatui::crossterm::event::{
    KeyCode, KeyEvent, KeyModifiers, MouseButton, MouseEvent, MouseEventKind,
};

/// Modifiers that turn a press into a drag. Terminals mostly keep Shift
/// with the mouse for selecting text, so Alt is the one that arrives.
const DRAG: KeyModifiers = KeyModifiers::ALT
    .union(KeyModifiers::SHIFT)
    .union(KeyModifiers::CONTROL);

/// What a drag holds until it is let go: a site row (dropped on the strip,
/// its gene is pinned there) or a pinned gene (dropped along the strip it
/// moves there, off it it is unpinned).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Held {
    Site(usize),
    Pin(usize),
}

impl ModulesView<'_> {
    /// A press: focus what it lands on. With a drag modifier a site row or
    /// a pinned gene is held until it is let go; a plain press only selects.
    fn press(&mut self, spot: Spot, drag: bool) {
        self.held = None;
        let held = match spot {
            Spot::Module(g) => {
                self.pane = Pane::Modules;
                self.select_module(g);
                None
            }
            Spot::Site(i) => {
                self.pane = Pane::Sites;
                self.select_site(i);
                Some(Held::Site(i))
            }
            Spot::Gene(k) | Spot::Strip(k, _) => {
                let site = match spot {
                    Spot::Strip(_, i) => i,
                    _ => self.pin_site,
                };
                self.pane = Pane::Genes;
                self.focus_pin(k, Some(site));
                Some(Held::Pin(k))
            }
            Spot::Pane(p) => {
                self.pane = p;
                None
            }
        };
        if !drag {
            return;
        }
        self.held = held;
        let say = |gene: &str, what: &str| format!("holding {}: {what}", query_symbol(gene));
        self.status = match held {
            Some(Held::Site(i)) => self.role().map(|g| {
                let gene = self.m.gene_of(g, i);
                say(gene, "let go on the genes pane to pin it there")
            }),
            Some(Held::Pin(k)) => self.pins.get(k).map(|p| {
                let gene = self.m.gene_of(p.role, self.m.genes(p.role)[p.gene].start);
                say(
                    gene,
                    "let go along the genes to move it, off them to unpin it",
                )
            }),
            None => None,
        };
    }

    /// A release: what a drag held is dropped where it is let go.
    fn release(&mut self, spot: Option<Spot>) -> bool {
        let Some(held) = self.held.take() else {
            return false;
        };
        self.status = None;
        // Where along the strip it lands: a pinned gene's place, or the end.
        let slot = match spot {
            Some(Spot::Gene(k) | Spot::Strip(k, _)) => Some(k),
            Some(Spot::Pane(Pane::Genes)) => Some(self.pins.len()),
            _ => None,
        };
        match (held, slot) {
            (Held::Site(i), Some(at)) => {
                self.pin_site_gene(i, at);
                self.pane = Pane::Genes;
            }
            (Held::Pin(k), Some(at)) => self.move_pin(k, at),
            (Held::Pin(k), None) => self.unpin(k),
            (Held::Site(_), None) => {}
        }
        true
    }

    /// The wheel over `spot`: the sites list scrolls and carries the focus
    /// with it; elsewhere it steps.
    fn wheel(&mut self, spot: Spot, up: bool) {
        let by: isize = if up { -1 } else { 1 };
        match spot {
            Spot::Site(_) | Spot::Pane(Pane::Sites) => {
                self.pane = Pane::Sites;
                let last = self.list.len().saturating_sub(self.rows);
                self.from = crate::tui::moved(self.from, 3 * by, last);
                let at = self.listed_at().unwrap_or(0);
                let end = (self.from + self.rows - 1).min(self.list.len().saturating_sub(1));
                if let Some(&i) = self.list.get(at.clamp(self.from, end)) {
                    self.site = i;
                }
            }
            Spot::Gene(_) | Spot::Strip(..) | Spot::Pane(Pane::Genes) => {
                self.pane = Pane::Genes;
                self.step_strip(by);
            }
            Spot::Module(_) | Spot::Pane(Pane::Modules) => {
                self.pane = Pane::Modules;
                self.step(by);
            }
        }
    }

    /// A mouse event, for `Screen::mouse`: whether it did anything.
    pub(super) fn on_mouse(&mut self, event: MouseEvent) -> bool {
        if self.find.is_some() || self.tracing {
            return false;
        }
        let spot = self.hits.at(event.column, event.row);
        let drag = event.modifiers.intersects(DRAG);
        // Sideways, or Shift with the wheel: the cell modules scroll.
        let shift = event.modifiers.contains(KeyModifiers::SHIFT);
        let sideways = match event.kind {
            MouseEventKind::ScrollLeft => Some(-1),
            MouseEventKind::ScrollRight => Some(1),
            MouseEventKind::ScrollUp if shift => Some(-1),
            MouseEventKind::ScrollDown if shift => Some(1),
            _ => None,
        };
        if let Some(by) = sideways {
            self.scroll_columns(by);
            return true;
        }
        match (event.kind, spot) {
            (MouseEventKind::ScrollUp, Some(s)) => self.wheel(s, true),
            (MouseEventKind::ScrollDown, Some(s)) => self.wheel(s, false),
            (MouseEventKind::Down(MouseButton::Left), Some(s)) => self.press(s, drag),
            (MouseEventKind::Up(MouseButton::Left), spot) => return self.release(spot),
            _ => return false,
        }
        true
    }

    /// A key, for `Screen::handle_key`.
    pub(super) fn on_key(&mut self, key: KeyEvent) {
        self.status = None;
        if self.find.is_some() {
            return self.find_key(key);
        }
        if self.tracing {
            return self.trace_key(key);
        }
        // The grid in focus, a place in `looks`.
        let grid = (self.pane != Pane::Genes).then(|| self.pane.index());
        let (sites, genes) = (self.pane == Pane::Sites, self.pane == Pane::Genes);
        match key.code {
            KeyCode::Char('q') => self.done = true,
            KeyCode::Char('e') => self.export(),
            KeyCode::Char('t') => self.open_trace(),
            KeyCode::Char('/') => self.open_find(),
            KeyCode::Tab => self.pane = self.pane.cycled(1),
            KeyCode::BackTab => self.pane = self.pane.cycled(PANES.len() - 1),
            KeyCode::Enter => match self.pane {
                Pane::Modules => self.pane = Pane::Sites,
                Pane::Sites => self.pane = Pane::Genes,
                Pane::Genes => {
                    if let Some(p) = self.pins.get(self.pin) {
                        let first = self.m.genes(p.role)[p.gene].start;
                        self.open = Some(self.m.gene_of(p.role, first).to_string());
                        self.done = true;
                    }
                }
            },
            // Back up a level; nothing above the modules.
            KeyCode::Esc => match self.pane {
                Pane::Genes => self.pane = Pane::Sites,
                Pane::Sites => self.pane = Pane::Modules,
                Pane::Modules => {}
            },
            // Pin, or select when pinned already; never unpin.
            KeyCode::Char(' ') if sites => self.pin_site_gene(self.site, usize::MAX),
            KeyCode::Delete | KeyCode::Backspace if sites => {
                let p = self.pin_of(self.site);
                if let Some(k) = self.pins.iter().position(|&q| Some(q) == p) {
                    self.unpin(k);
                }
            }
            KeyCode::Delete | KeyCode::Backspace if genes => self.unpin(self.pin),
            KeyCode::Char('X') => {
                self.pins.clear();
                (self.pin, self.row_from) = (0, 0);
            }
            KeyCode::Char('v') => {
                if let Some(k) = grid {
                    self.looks[k] = self.looks[k].next();
                }
            }
            KeyCode::Char('y') if grid.is_some_and(|k| self.looks[k] == Look::Percent) => {
                self.ramp = self.ramp.next()
            }
            KeyCode::Char('o') if sites => self.next_sort(),
            KeyCode::Char('O') if sites && self.sort != Sort::Module => {
                self.resort(self.sort, !self.descending)
            }
            KeyCode::Char('[') => self.scroll_columns(-1),
            KeyCode::Char(']') => self.scroll_columns(1),
            KeyCode::Down | KeyCode::Char('j') => self.step(1),
            KeyCode::Up | KeyCode::Char('k') => self.step(-1),
            KeyCode::Right | KeyCode::Char('l' | 'n') if genes => self.step_strip(1),
            KeyCode::Left | KeyCode::Char('h' | 'N') if genes => self.step_strip(-1),
            KeyCode::PageDown if sites => self.step_listed(self.rows as isize),
            KeyCode::PageUp if sites => self.step_listed(-(self.rows as isize)),
            KeyCode::Home if sites => self.step_listed(isize::MIN / 2),
            KeyCode::End if sites => self.step_listed(isize::MAX / 2),
            _ => {}
        }
    }
}
