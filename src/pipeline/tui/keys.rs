//! What each key does, screen by screen.

use std::path::PathBuf;

use ratatui::crossterm::event::{KeyCode, KeyEvent};

use super::form::Kind;
use super::inputs::{Browser, InputsFocus, Row};
use super::steps::Step;
use super::{script, App, Page, Target};
use crate::figure::{Edit, LineInput};
use crate::tui::is_go;

/// Longest text a line input takes.
const MAX_TYPED: usize = 4096;

impl App {
    /// Route `key`: an open line input first, then the file pop-up, the
    /// preview, the keys every screen shares, and the screen's own.
    pub(super) fn key(&mut self, key: KeyEvent) {
        self.note = None;
        if self.editing.is_some() {
            return self.key_editing(key);
        }
        if self.picking.is_some() {
            return self.key_picking(key);
        }
        if self.preview {
            return self.key_preview(key);
        }
        if is_go(&key) {
            self.preview = true;
            return;
        }
        match key.code {
            KeyCode::Tab => return self.turn(1),
            KeyCode::BackTab => return self.turn(-1),
            KeyCode::Char(c @ '1'..='4') => {
                let page = self.pages().get(c as usize - '1' as usize).copied();
                if let Some(p) = page {
                    self.page = p;
                }
                return;
            }
            KeyCode::Char('q') => return self.leave(),
            _ => {}
        }
        match self.page {
            Page::Inputs => self.key_inputs(key),
            Page::Steps => self.key_steps(key),
            Page::Flags => self.key_flags(key),
            Page::Run => self.key_run(key),
        }
    }

    /// The screens a tab reaches: Run once there is a run.
    fn pages(&self) -> Vec<Page> {
        let mut v = vec![Page::Inputs, Page::Steps, Page::Flags];
        if self.job.is_some() {
            v.push(Page::Run);
        }
        v
    }

    fn turn(&mut self, d: isize) {
        let pages = self.pages();
        let n = pages.len() as isize;
        let at = pages.iter().position(|p| *p == self.page).unwrap_or(0) as isize;
        if let Some(p) = pages.get((at + d).rem_euclid(n) as usize) {
            self.page = *p;
        }
    }

    /// Open a line input for `target`, holding `initial`.
    fn edit(&mut self, target: Target, initial: &str) {
        let mut line = LineInput::new(MAX_TYPED);
        line.open(initial);
        self.editing = Some((target, line));
    }

    fn key_editing(&mut self, key: KeyEvent) {
        let Some((_, line)) = self.editing.as_mut() else {
            return;
        };
        match line.handle(key) {
            Edit::Typing => {}
            Edit::Cancelled => self.editing = None,
            Edit::Submitted(t) => {
                let Some((target, _)) = self.editing.take() else {
                    return;
                };
                match target {
                    Target::Output => self.inputs.output = t,
                    Target::DepthKb => self.steps.depth_kb = t,
                    Target::Flag(i) => {
                        if let Some(f) = self.form.fields.get_mut(i) {
                            f.set_text(&t);
                        }
                    }
                    Target::Find => {
                        self.find = t;
                        self.flags_at = 0;
                    }
                }
            }
        }
    }

    fn key_picking(&mut self, key: KeyEvent) {
        let Some((row, browser)) = self.picking.as_mut() else {
            return;
        };
        let row = *row;
        match key.code {
            KeyCode::Up | KeyCode::Char('k') => browser.step(-1),
            KeyCode::Down | KeyCode::Char('j') => browser.step(1),
            KeyCode::PageUp => browser.step(-10),
            KeyCode::PageDown => browser.step(10),
            KeyCode::Left | KeyCode::Char('h') => browser.up(),
            KeyCode::Enter | KeyCode::Right => {
                if let Some(file) = browser.enter() {
                    self.inputs.set(row, Some(&file));
                    self.picking = None;
                }
            }
            KeyCode::Esc => self.picking = None,
            _ => {}
        }
    }

    fn key_preview(&mut self, key: KeyEvent) {
        if is_go(&key) || key.code == KeyCode::Char('y') {
            if let Err(e) = self.start() {
                self.note = Some(format!("cannot start: {e}"));
            }
            return;
        }
        match key.code {
            KeyCode::Char('c') => {
                let text = script::command_lines(&self.argv()).join(" \\\n  ");
                self.copy(&text);
            }
            KeyCode::Esc => self.preview = false,
            _ => {}
        }
    }

    fn key_inputs(&mut self, key: KeyEvent) {
        match self.inputs.focus {
            InputsFocus::Bams => {
                let bams = &mut self.inputs.bams;
                match key.code {
                    KeyCode::Up | KeyCode::Char('k') => bams.step(-1),
                    KeyCode::Down | KeyCode::Char('j') => bams.step(1),
                    KeyCode::PageUp => bams.step(-10),
                    KeyCode::PageDown => bams.step(10),
                    KeyCode::Enter | KeyCode::Right => {
                        // Folders open; a BAM is picked with Space.
                        if bams.entries.get(bams.at).is_some_and(|e| e.dir) {
                            bams.enter();
                        }
                    }
                    KeyCode::Left | KeyCode::Backspace => bams.up(),
                    KeyCode::Char(' ') => self.inputs.toggle(),
                    KeyCode::Char('b') => self.inputs.flip(),
                    KeyCode::Char('l') => self.inputs.focus = InputsFocus::Rows,
                    _ => {}
                }
            }
            InputsFocus::Rows => match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    self.inputs.row = self.inputs.row.saturating_sub(1);
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    self.inputs.row = (self.inputs.row + 1).min(Row::ALL.len() - 1);
                }
                KeyCode::Enter => self.open_row(),
                KeyCode::Char('h') | KeyCode::Esc => self.inputs.focus = InputsFocus::Bams,
                _ => {}
            },
        }
    }

    /// Enter on an inputs row: a file pop-up, or a line to type in.
    fn open_row(&mut self) {
        let Some(&row) = Row::ALL.get(self.inputs.row) else {
            return;
        };
        match row {
            Row::File(f) => {
                let start: PathBuf = f
                    .get(&self.inputs)
                    .and_then(|p| p.parent())
                    .filter(|p| p.is_dir())
                    .map_or_else(|| self.inputs.bams.cwd.clone(), |p| p.to_path_buf());
                self.picking = Some((f, Browser::new(start, f.ext())));
            }
            Row::Output => self.edit(Target::Output, &self.inputs.output()),
            Row::Threads => {
                let i = self
                    .form
                    .fields
                    .iter()
                    .position(|f| f.long == "max-threads");
                if let Some(i) = i {
                    let value = self.form.fields[i].value.clone();
                    self.edit(Target::Flag(i), &value);
                }
            }
        }
    }

    fn key_steps(&mut self, key: KeyEvent) {
        let n = Step::ALL.len();
        match key.code {
            KeyCode::Up | KeyCode::Char('k') => self.steps.at = self.steps.at.saturating_sub(1),
            KeyCode::Down | KeyCode::Char('j') => self.steps.at = (self.steps.at + 1).min(n - 1),
            KeyCode::Char(' ') => {
                let has_bg = !self.inputs.bg().is_empty();
                self.note = self.steps.toggle(has_bg).map(String::from);
                let depth = Step::ALL.get(self.steps.at) == Some(&Step::Depth);
                if depth && self.steps.is_on(Step::Depth) {
                    let kb = self.steps.depth_kb.clone();
                    self.edit(Target::DepthKb, &kb);
                }
            }
            KeyCode::Enter if Step::ALL.get(self.steps.at) == Some(&Step::Depth) => {
                let kb = self.steps.depth_kb.clone();
                self.edit(Target::DepthKb, &kb);
            }
            _ => {}
        }
    }

    fn key_flags(&mut self, key: KeyEvent) {
        let n = self.visible_flags().len();
        let last = n.saturating_sub(1);
        let at = self.flags_at.min(last);
        let row = self.flag_row();
        match key.code {
            KeyCode::Up | KeyCode::Char('k') => self.flags_at = at.saturating_sub(1),
            KeyCode::Down | KeyCode::Char('j') => self.flags_at = (at + 1).min(last),
            KeyCode::PageUp => self.flags_at = at.saturating_sub(10),
            KeyCode::PageDown => self.flags_at = (at + 10).min(last),
            KeyCode::Home => self.flags_at = 0,
            KeyCode::End => self.flags_at = last,
            KeyCode::Char('/') => {
                let find = self.find.clone();
                self.edit(Target::Find, &find);
            }
            KeyCode::Char('a') => {
                self.advanced ^= true;
                self.flags_at = 0;
            }
            KeyCode::Char('R') => self.form.fields.iter_mut().for_each(|f| f.reset()),
            _ => {
                let Some(i) = row else { return };
                let Some(f) = self.form.fields.get_mut(i) else {
                    return;
                };
                match (key.code, &f.kind) {
                    (KeyCode::Char(' '), Kind::Flag { .. }) => f.toggle(),
                    (KeyCode::Left, Kind::Choice(_)) => f.cycle(-1),
                    (KeyCode::Right | KeyCode::Char(' '), Kind::Choice(_)) => f.cycle(1),
                    (KeyCode::Char('r'), _) => f.reset(),
                    (KeyCode::Enter, Kind::Flag { .. }) => f.toggle(),
                    (KeyCode::Enter, _) => {
                        let value = f.value.clone();
                        self.edit(Target::Flag(i), &value);
                    }
                    _ => {}
                }
            }
        }
    }

    fn key_run(&mut self, key: KeyEvent) {
        let Some(job) = self.job.as_mut() else {
            return;
        };
        // As far up as the Run screen can show: the top line in its top row.
        let lines = job.log.lock().map_or(0, |l| l.lines.len());
        let top = lines.saturating_sub(job.rows.get());
        let scroll = job.scroll.min(top);
        match key.code {
            KeyCode::Up | KeyCode::Char('k') => job.scroll = (scroll + 1).min(top),
            KeyCode::Down | KeyCode::Char('j') => job.scroll = scroll.saturating_sub(1),
            KeyCode::PageUp => job.scroll = (scroll + 10).min(top),
            KeyCode::PageDown => job.scroll = scroll.saturating_sub(10),
            KeyCode::End => job.scroll = 0,
            KeyCode::Char('s') if job.running() => {
                if job.asking || job.stopper.is_stopped() {
                    job.asking = false;
                    job.stopper.stop();
                } else {
                    job.asking = true;
                }
            }
            KeyCode::Esc => job.asking = false,
            _ => {}
        }
    }
}
