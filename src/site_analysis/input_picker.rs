//! The pop-up `pileup` and `metagene` open when the command line names no
//! input: browse, mark the files to read together, or choose a whole faba
//! output folder.

use std::path::PathBuf;

use data_beans::interactive::ui::{header, help_line, Screen, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::qc::layout::looks_like_faba_dir;
use crate::qc::path_tui::faba_browser;
use crate::tui::browser::{Browser, Nav};
use crate::tui::{is_apply, popup_frame, run_view, tilde, View, APPLY_KEYS};

/// The choice, independent of the terminal so it can be tested.
pub struct InputPicker {
    /// Marking files, in the order marked, from any folder.
    browser: Browser,
    /// The command asking, for the header.
    command: &'static str,
    /// What the listed files are, for the title.
    what: &'static str,
    /// Why the last choice was refused.
    pub error: Option<String>,
    /// The files or the folder chosen; `Some(None)` when cancelled.
    pub decision: Option<Option<Vec<PathBuf>>>,
}

impl InputPicker {
    /// Listing `cwd`'s folders and the files `keep` accepts by name.
    pub fn new(
        cwd: PathBuf,
        command: &'static str,
        what: &'static str,
        keep: fn(&str) -> bool,
    ) -> Self {
        let mut browser = faba_browser(cwd.clone(), keep).marking();
        browser.open(cwd);
        Self {
            browser,
            command,
            what,
            error: None,
            decision: None,
        }
    }

    /// Take the whole of `dir`, if it is a faba output folder.
    fn choose_dir(&mut self, dir: PathBuf) {
        if looks_like_faba_dir(&dir) {
            self.decision = Some(Some(vec![dir]));
        } else {
            self.error = Some(format!("{} is not a faba output folder", tilde(&dir)));
        }
    }
}

impl View for InputPicker {
    /// The apply key takes the marked files, or the folder shown.
    fn takes_apply(&self) -> bool {
        true
    }
}

impl Screen for InputPicker {
    fn reports_chords(&self) -> bool {
        true
    }

    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn interrupt(&mut self) {
        self.decision = Some(None);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.error = None;
        if is_apply(&key) {
            let marked = self.browser.marked.take().unwrap_or_default();
            if marked.is_empty() {
                self.browser.marked = Some(marked);
                return self.choose_dir(self.browser.cwd.clone());
            }
            self.decision = Some(Some(marked));
            return;
        }
        match self.browser.key(key) {
            Nav::Picked(file) => return self.browser.toggle_mark(file),
            Nav::Moved => return,
            Nav::Ignored => {}
        }
        match key.code {
            KeyCode::Char(' ') => {
                let Some(e) = self.browser.highlighted() else {
                    return;
                };
                let (path, dir, up) = (e.path.clone(), e.dir, e.name == "..");
                if !dir {
                    self.browser.toggle_mark(path);
                } else if !up {
                    self.choose_dir(path);
                }
            }
            KeyCode::Esc => self.decision = Some(None),
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        frame.render_widget(header(self.command, "choose the input", ""), area);
        let w = area.width.saturating_sub(4).clamp(20, 90);
        let rows = area.height.saturating_sub(8).clamp(3, 24);
        let n_marked = self.browser.marked.as_ref().map_or(0, Vec::len);
        let marked = match n_marked {
            0 => String::new(),
            n => format!("· {n} marked "),
        };
        let title = format!(
            " {}, or a faba output folder · {} {}{marked}",
            self.what,
            tilde(&self.browser.cwd),
            self.browser.list.find.tag()
        );
        let inner = popup_frame(frame, area, w, rows + 4, title);
        let mut lines = self.browser.lines(
            rows as usize,
            inner.width as usize,
            "faba output",
            "no subfolder or file to read",
        );
        lines.push(match &self.error {
            Some(why) => Line::from(Span::styled(format!(" {why}"), HIGHLIGHT)),
            None => Line::from(""),
        });
        let apply = if n_marked == 0 {
            "choose this folder"
        } else {
            "read the marked"
        };
        lines.push(help_line(&[
            ("type", "find"),
            ("↑/↓", "move"),
            ("Enter/→", "open"),
            ("←", "up"),
            ("Space", "mark file / choose folder"),
            (APPLY_KEYS, apply),
            ("Esc", "clear find / cancel"),
        ]));
        frame.render_widget(Paragraph::new(lines), inner);
    }
}

/// The input when the command line named none: the files marked in a
/// browser (those `keep` accepts by name), or one faba output folder.
/// `Ok(None)` when the user cancels; an error when there is no browsing
/// (`batch`, or no terminal), saying what to name.
pub fn ask_inputs(
    command: &'static str,
    batch: bool,
    what: &'static str,
    keep: fn(&str) -> bool,
) -> anyhow::Result<Option<Vec<Box<str>>>> {
    anyhow::ensure!(
        !batch && data_beans::interactive::tui_available(),
        "name {what}, or a faba output directory"
    );
    let mut picker = InputPicker::new(std::env::current_dir()?, command, what, keep);
    run_view(&mut picker)?;
    Ok(picker.decision.flatten().map(|paths| {
        paths
            .iter()
            .map(|p| p.to_string_lossy().into_owned().into_boxed_str())
            .collect()
    }))
}

#[cfg(test)]
#[path = "tests/input_picker.rs"]
mod tests;
