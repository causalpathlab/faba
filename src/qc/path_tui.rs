//! Pop-up asking `faba qc` for the input and output directories when the
//! command line left them out: browse to a faba output directory, then name
//! a new one to write to.

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::{header, help_line, input_line, Screen, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;

use super::browser::{Browser, Listing, Nav};
use super::layout::looks_like_faba_dir;
use crate::figure::{Edit, LineInput};
use crate::tui::{as_apply, is_stray, next_free, output_problem, popup_frame};

enum Step {
    /// Browsing for the input directory.
    Input,
    /// Naming the output directory for the chosen `input`.
    Output { input: PathBuf, line: LineInput },
}

/// The two questions, independent of the terminal so they can be tested.
struct PathPicker {
    step: Step,
    /// Subdirectories only, tagged when they hold faba matrices or site
    /// tables; listed when the browser first shows.
    browser: Browser,
    /// Directories already checked for faba files.
    faba: FxHashMap<PathBuf, bool>,
    /// Why the last choice was refused.
    error: Option<String>,
    /// The output as given on the command line, if it was.
    given_output: Option<String>,
    decision: Option<Option<(String, String)>>,
}

/// Whether `dir` holds faba matrices or site tables, checked once.
fn cached_faba(cache: &mut FxHashMap<PathBuf, bool>, dir: &Path) -> bool {
    *cache
        .entry(dir.to_path_buf())
        .or_insert_with(|| looks_like_faba_dir(dir))
}

/// A directory next to `input` named after it, that does not exist yet.
fn suggest_output(input: &Path) -> String {
    let base = input.to_string_lossy();
    // Joined onto an empty base, the stem is the path itself.
    next_free(Path::new(""), &format!("{}_qc", base.trim_end_matches('/')))
        .to_string_lossy()
        .into_owned()
}

impl PathPicker {
    fn new(cwd: PathBuf, input: Option<PathBuf>, given_output: Option<String>) -> Self {
        let mut p = Self {
            step: Step::Input,
            browser: Browser::new(cwd.clone()),
            faba: FxHashMap::default(),
            error: None,
            given_output,
            decision: None,
        };
        // Skip what the command line already answered.
        if let Some(dir) = input {
            p.choose(dir);
        }
        if p.decision.is_none() && matches!(p.step, Step::Input) {
            p.open(cwd);
        }
        p
    }

    /// `f` on the browser, listing subdirectories only and tagging faba
    /// output directories.
    fn browse<R>(&mut self, f: impl FnOnce(&mut Browser, &mut Listing) -> R) -> R {
        let faba = &mut self.faba;
        let mut tag = |d: &Path| cached_faba(faba, d);
        let mut listing = Listing {
            keep: &|_| false,
            tag: &mut tag,
        };
        f(&mut self.browser, &mut listing)
    }

    fn open(&mut self, dir: PathBuf) {
        self.browse(|b, l| b.open(dir, l));
    }

    fn is_faba(&mut self, dir: &Path) -> bool {
        cached_faba(&mut self.faba, dir)
    }

    /// Take `dir` as the input, if it is a faba output directory.
    fn choose(&mut self, dir: PathBuf) {
        if !self.is_faba(&dir) {
            self.error = Some(format!(
                "{} has no faba matrices or site tables",
                dir.display()
            ));
            return;
        }
        let out = self
            .given_output
            .clone()
            .unwrap_or_else(|| suggest_output(&dir));
        if self.given_output.is_some() && output_problem(&out).is_none() {
            return self.finish(&dir, out);
        }
        self.error = self.given_output.as_deref().and_then(output_problem);
        let mut line = LineInput::new(4096);
        line.open(&out);
        self.step = Step::Output { input: dir, line };
    }

    fn finish(&mut self, input: &Path, out: String) {
        self.decision = Some(Some((input.to_string_lossy().into_owned(), out)));
    }

    fn browse_key(&mut self, key: KeyEvent) {
        if is_stray(&key) {
            return;
        }
        if !matches!(self.browse(|b, l| b.key(key, l)), Nav::Ignored) {
            return;
        }
        match key.code {
            KeyCode::Char(' ') => {
                if let Some((dir, _)) = self.browser.highlighted() {
                    self.choose(dir);
                }
            }
            KeyCode::Char('.') => self.choose(self.browser.cwd.clone()),
            KeyCode::Esc | KeyCode::Char('q') => self.decision = Some(None),
            _ => {}
        }
    }
}

impl Screen for PathPicker {
    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn interrupt(&mut self) {
        self.decision = Some(None);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        let key = as_apply(key);
        self.error = None;
        let Step::Output { input, line } = &mut self.step else {
            return self.browse_key(key);
        };
        match line.handle(key) {
            Edit::Typing => {}
            Edit::Submitted(out) => match output_problem(&out) {
                None => {
                    let input = input.clone();
                    self.finish(&input, out);
                }
                Some(why) => {
                    line.open(&out);
                    self.error = Some(why);
                }
            },
            Edit::Cancelled => {
                self.step = Step::Input;
                if self.browser.entries.is_empty() {
                    self.open(self.browser.cwd.clone());
                }
            }
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        frame.render_widget(header("qc", "choose the directories", ""), area);
        let w = area.width.saturating_sub(4).clamp(20, 90);
        let (title, help, rows) = match &self.step {
            Step::Input => (
                format!(
                    " input: a faba output directory · {} ",
                    self.browser.cwd.display()
                ),
                help_line(&[
                    ("↑/↓", "move"),
                    ("Enter/→", "open"),
                    ("←", "up"),
                    ("Space", "choose highlighted"),
                    (".", "choose this one"),
                    ("Esc", "cancel"),
                ]),
                area.height.saturating_sub(8).clamp(3, 24),
            ),
            Step::Output { .. } => (
                " output: a new or empty directory ".to_string(),
                help_line(&[("Enter", "next"), ("Esc", "back")]),
                3,
            ),
        };
        let inner = popup_frame(frame, area, w, rows + 4, title);
        let mut lines = match &self.step {
            Step::Input => self.browser.lines(
                rows as usize,
                inner.width as usize,
                "faba output",
                "no subdirectory",
            ),
            Step::Output { input, line } => vec![
                Line::from(vec![
                    Span::styled(" input  ", DIM),
                    Span::raw(input.display().to_string()),
                ]),
                Line::from(""),
                input_line("output ", line.text().unwrap_or_default(), &[]),
            ],
        };
        lines.push(match &self.error {
            Some(why) => Line::from(Span::styled(format!(" {why}"), HIGHLIGHT)),
            None => Line::from(""),
        });
        lines.push(help);
        frame.render_widget(Paragraph::new(lines), inner);
    }
}

/// Ask for whichever of `input` and `output` is missing. `None` when the
/// user cancels.
pub fn ask_paths(
    input: Option<&str>,
    output: Option<&str>,
) -> anyhow::Result<Option<(String, String)>> {
    let cwd = std::env::current_dir()?;
    let mut picker = PathPicker::new(cwd, input.map(PathBuf::from), output.map(String::from));
    if !picker.done() {
        data_beans::interactive::ui::run_screen(&mut picker)?;
    }
    Ok(picker.decision.flatten())
}

#[cfg(test)]
#[path = "tests/path_tui.rs"]
mod tests;
