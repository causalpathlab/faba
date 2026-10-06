//! What each key does, screen by screen.

use std::path::PathBuf;

use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent, MouseEventKind};

use super::fetch::{bam_chr_named, Catalogue, Fetch};
use super::form::Kind;
use super::inputs::{InputsFocus, Row};
use super::steps::Step;
use super::{script, App, Hit, Page, Target};
use crate::figure::{Edit, LineInput};
use crate::tui::browser::Nav;
use crate::tui::{apply_key, enter_hint, is_apply, moved, output_problem, typed_path, wheel_key};

/// Longest text a line input takes.
const MAX_TYPED: usize = 4096;

/// A move down a list (up when negative) for the keys every list shares:
/// ↑/k, ↓/j, PgUp and PgDn.
fn nav(code: KeyCode) -> Option<isize> {
    match code {
        KeyCode::Up | KeyCode::Char('k') => Some(-1),
        KeyCode::Down | KeyCode::Char('j') => Some(1),
        KeyCode::PageUp => Some(-10),
        KeyCode::PageDown => Some(10),
        _ => None,
    }
}

impl App {
    /// Take `key`, then check the command it leaves.
    pub(super) fn key(&mut self, key: KeyEvent) {
        self.route(key);
        self.refresh();
    }

    /// Route `key`: an open line input first, then the file pop-up, the
    /// preview, the keys every screen shares, and the screen's own.
    fn route(&mut self, key: KeyEvent) {
        self.note = None;
        let leaving = std::mem::take(&mut self.leaving);
        if self.editing.is_some() {
            return self.key_editing(key);
        }
        if self.ended_shown {
            // The run's end, told: Enter or Esc closes it, onto the log.
            if matches!(key.code, KeyCode::Enter | KeyCode::Esc) {
                self.ended_shown = false;
            }
            return;
        }
        if self.picking.is_some() {
            return self.key_picking(key);
        }
        if self.catalogue.is_some() {
            return self.key_catalogue(key);
        }
        if self.fetch_shown && self.fetch.is_some() {
            return self.key_fetch(key);
        }
        if self.preview {
            return self.key_preview(key);
        }
        if is_apply(&key) {
            if !self.ask_output(key) {
                self.preview = true;
            }
            return;
        }
        // The BAM list takes its keys first: letters type there, so they
        // are no commands. Enter on a BAM picks nothing; Space does.
        if self.page == Page::Inputs && self.inputs.focus == InputsFocus::Bams {
            match self.inputs.bams.key(key) {
                Nav::Moved => return,
                Nav::Picked(_) => return self.enter_hint(),
                Nav::Ignored => {}
            }
        }
        match key.code {
            KeyCode::Char('c' | 'p') if self.ask_output(key) => return,
            KeyCode::Char('c') => return self.copy_command(),
            KeyCode::Char('p') => return self.print_and_leave(),
            KeyCode::Tab => return self.turn(1),
            KeyCode::BackTab => return self.turn(-1),
            KeyCode::Char(c @ '1'..='4') => {
                let page = self.pages().get(c as usize - '1' as usize).copied();
                if let Some(p) = page {
                    self.page = p;
                }
                return;
            }
            // A set-up run is lost on leaving: the first `q` only asks.
            KeyCode::Char('q')
                if !leaving && self.job.is_none() && !self.inputs.picked.is_empty() =>
            {
                self.leaving = true;
                self.note = Some("q again leaves, losing this setup; p prints its command".into());
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

    /// A click or a turn of the wheel: focus and select what is under it;
    /// the wheel then moves as ↑/↓ would. Nothing while a pop-up is up.
    pub(super) fn click(&mut self, m: MouseEvent) {
        let Some(hit) = self.hits.at(m.column, m.row) else {
            return;
        };
        let clicked = matches!(m.kind, MouseEventKind::Down(_));
        let popped = self.popped();
        // A button does what its key does, the same way.
        let esc = KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE);
        let button = match hit {
            _ if !clicked => None,
            Hit::Preview if !popped => Some(apply_key()),
            Hit::Start if self.preview => Some(apply_key()),
            Hit::Cancel if self.preview => Some(esc),
            Hit::Ok if self.ended_shown => Some(esc),
            Hit::Stop if !popped => {
                self.page = Page::Run;
                Some(KeyEvent::new(KeyCode::Char('s'), KeyModifiers::NONE))
            }
            _ => None,
        };
        if let Some(key) = button {
            return self.key(key);
        }
        if popped {
            return;
        }
        // Focus and selection only: nothing the command depends on changes.
        self.note = None;
        match hit {
            Hit::Tab(p) if clicked && self.pages().contains(&p) => self.page = p,
            Hit::Bams(i) => {
                self.inputs.focus = InputsFocus::Bams;
                if let Some(i) = i.filter(|_| clicked) {
                    self.inputs.bams.list.at = i;
                }
            }
            Hit::Rows(i) => {
                self.inputs.focus = InputsFocus::Rows;
                if let Some(i) = i.filter(|_| clicked) {
                    self.inputs.row = i;
                }
            }
            Hit::Step(i) if clicked => self.steps.at = i,
            Hit::Flag(k) if clicked => self.flags_at = k,
            _ => {}
        }
        if let Some(key) = wheel_key(&m) {
            self.key(key);
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

    /// Tab: the BAM list, then the inputs rows, then the next screen;
    /// Shift+Tab back.
    fn turn(&mut self, d: isize) {
        if self.page == Page::Inputs {
            let to = match (self.inputs.focus, d > 0) {
                (InputsFocus::Bams, true) => Some(InputsFocus::Rows),
                (InputsFocus::Rows, false) => Some(InputsFocus::Bams),
                _ => None,
            };
            if let Some(focus) = to {
                self.inputs.focus = focus;
                return;
            }
        }
        let pages = self.pages();
        let n = pages.len() as isize;
        let at = pages.iter().position(|p| *p == self.page).unwrap_or(0) as isize;
        if let Some(p) = pages.get((at + d).rem_euclid(n) as usize) {
            self.page = *p;
        }
        if self.page == Page::Inputs {
            self.inputs.focus = if d > 0 {
                InputsFocus::Bams
            } else {
                InputsFocus::Rows
            };
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
            Edit::Cancelled => {
                self.editing = None;
                self.asked_output = None;
            }
            Edit::Submitted(t) => {
                let Some((target, _)) = self.editing.take() else {
                    return;
                };
                match target {
                    Target::Output if t.is_empty() => {
                        self.inputs.output = t;
                        self.asked_output = None;
                    }
                    Target::Output => {
                        let out = typed_path(&t).to_string_lossy().into_owned();
                        // A folder with files in it is refused at once.
                        if let Some(why) = output_problem(&out) {
                            self.edit(Target::Output, &t);
                            self.note = Some(why);
                            return;
                        }
                        self.inputs.output = out;
                        // Asked for by a key: now go on with what it asked.
                        if let Some(asked) = self.asked_output.take() {
                            self.route(asked);
                        }
                    }
                    Target::DepthKb => self.steps.depth_kb = t,
                    Target::Flag(i) => {
                        if let Some(f) = self.form.fields.get_mut(i) {
                            f.set_text(&t);
                        }
                    }
                    Target::RefDir => {
                        if let Some(pick) = self.ref_pick.take().filter(|_| !t.is_empty()) {
                            self.fetch = Some(Fetch::start(pick, typed_path(&t)));
                            self.fetch_shown = true;
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
        match browser.key(key) {
            Nav::Moved => {}
            Nav::Picked(file) => {
                self.inputs.set(row, Some(&file));
                self.picking = None;
            }
            Nav::Ignored if key.code == KeyCode::Esc => self.picking = None,
            Nav::Ignored => {}
        }
    }

    /// Enter did nothing: say what the apply key does here, and what to do
    /// on a terminal that cannot send it.
    pub(super) fn enter_hint(&mut self) {
        let does = if self.preview {
            "starts the run"
        } else {
            "previews the run"
        };
        self.note = Some(format!(
            "{}; or c copies the faba run --batch-process command, p prints it and leaves",
            enter_hint(does)
        ));
    }

    /// `d` on the inputs rows: the references the FTP sites offer, unless
    /// one is coming already.
    fn offer_references(&mut self) {
        if self.fetch.is_some() {
            self.fetch_shown = true;
            return;
        }
        let bam_chr = self.inputs.fg().first().and_then(|b| bam_chr_named(b));
        self.catalogue = Some(Catalogue::new(bam_chr));
    }

    /// The download's pop-up: Esc hides it, the download going on; `s`
    /// stops it, keeping what came for the next try.
    fn key_fetch(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Esc => self.fetch_shown = false,
            KeyCode::Char('s') => {
                if let Some(f) = &self.fetch {
                    f.stop();
                }
            }
            _ => {}
        }
    }

    fn key_catalogue(&mut self, key: KeyEvent) {
        let Some(c) = self.catalogue.as_mut() else {
            return;
        };
        if c.key(&key) {
            return;
        }
        match key.code {
            KeyCode::Enter | KeyCode::Right => {
                if let Some(pick) = c.enter() {
                    self.catalogue = None;
                    let base = crate::tui::home().unwrap_or_else(|| self.inputs.bams.cwd.clone());
                    let dir = base.join("faba_refs").join(pick.dir_name());
                    self.ref_pick = Some(pick);
                    self.edit(Target::RefDir, &dir.to_string_lossy());
                }
            }
            KeyCode::Esc | KeyCode::Left if !c.back() => self.catalogue = None,
            _ => {}
        }
    }

    /// No output folder named yet: open the output row's line, holding a
    /// suggestion, and take `key` again once one is named. The run never
    /// picks a folder by itself.
    fn ask_output(&mut self, key: KeyEvent) -> bool {
        if !self.inputs.output().is_empty() {
            return false;
        }
        self.page = Page::Inputs;
        self.inputs.focus = InputsFocus::Rows;
        self.inputs.row = Row::ALL.iter().position(|r| *r == Row::Output).unwrap_or(0);
        self.asked_output = Some(key);
        let suggested = self.inputs.suggested_output();
        self.edit(Target::Output, &suggested);
        true
    }

    /// Copy the exact command the run would start.
    fn copy_command(&mut self) {
        let text = script::command_lines(&self.argv_out()).join(" \\\n  ");
        self.copy(&text);
    }

    fn key_preview(&mut self, key: KeyEvent) {
        // Only the apply key starts: a stray key cannot start a run.
        if is_apply(&key) {
            if let Err(e) = self.start() {
                self.note = Some(format!("cannot start: {e}"));
            }
            return;
        }
        match key.code {
            KeyCode::Char('c') => self.copy_command(),
            KeyCode::Char('p') => self.print_and_leave(),
            KeyCode::Enter => self.enter_hint(),
            KeyCode::Esc => self.preview = false,
            _ => {}
        }
    }

    fn key_inputs(&mut self, key: KeyEvent) {
        match self.inputs.focus {
            // The rest of the list's keys went to the browser in `route`.
            InputsFocus::Bams => {
                if key.code == KeyCode::Char(' ') {
                    self.inputs.toggle();
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
                KeyCode::Char('d') => self.offer_references(),
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
                self.picking = Some((f, f.browser(start)));
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
            KeyCode::Enter => self.enter_hint(),
            _ => {}
        }
    }

    fn key_flags(&mut self, key: KeyEvent) {
        let n = self.visible_flags().len();
        let last = n.saturating_sub(1);
        let at = self.flags_at.min(last);
        let row = self.flag_row();
        if let Some(d) = nav(key.code) {
            self.flags_at = moved(at, d, last);
            return;
        }
        match key.code {
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
                let Some(i) = row else {
                    if key.code == KeyCode::Enter {
                        self.enter_hint();
                    }
                    return;
                };
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
        if let Some(d) = nav(key.code) {
            // Scrolled up is lines from the end: a move up adds to it.
            job.scroll = moved(scroll, -d, top);
            return;
        }
        match key.code {
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
            KeyCode::Enter => self.enter_hint(),
            _ => {}
        }
    }
}
