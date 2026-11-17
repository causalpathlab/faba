//! Pop-up asking `faba qc` for the input and output directories when the
//! command line left them out: browse to a faba output directory, then name
//! a new one to write to.

use std::borrow::Cow;
use std::path::{Path, PathBuf};

use crate::tui::ui::{header, help_line, input_line, Screen, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;

use super::layout::looks_like_faba_dir;
use crate::figure::{Edit, LineInput};
use crate::tui::browser::{self, Browser, Outcome, Simple};
use crate::tui::{is_apply, next_free, output_problem, popup_frame, run_view, View, APPLY_KEYS};

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
    /// Why the last choice was refused.
    error: Option<String>,
    /// The output as given on the command line, if it was.
    given_output: Option<String>,
    decision: Option<Option<(String, String)>>,
}

/// What a browser of subdirectories and the files `keep` accepts wants (a
/// popup `title` for `what`): the folders that hold faba matrices or site
/// tables tagged, each checked once.
pub fn faba_want(
    title: &'static str,
    what: &'static str,
    keep: impl Fn(&str) -> bool + Send + Sync + 'static,
) -> Simple {
    let mut seen: FxHashMap<PathBuf, bool> = FxHashMap::default();
    let tag = move |path: &Path, is_dir: bool| {
        if !is_dir {
            return None;
        }
        let faba = *seen
            .entry(path.to_path_buf())
            .or_insert_with(|| looks_like_faba_dir(path));
        faba.then_some(Cow::Borrowed("faba output"))
    };
    browser::want(title, what, keep, tag)
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
            browser: browser::open(cwd, faba_want("input", "subdirectory", |_| false), None),
            error: None,
            given_output,
            decision: None,
        };
        // Skip what the command line already answered.
        if let Some(dir) = input {
            p.choose(dir);
        }
        p
    }

    /// Take `dir` as the input, if it is a faba output directory.
    fn choose(&mut self, dir: PathBuf) {
        if !looks_like_faba_dir(&dir) {
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
        // Letters narrow the listing, so the one being shown is chosen with
        // the apply key.
        if is_apply(&key) {
            return self.choose(self.browser.dir.clone());
        }
        match self.browser.key(key) {
            Outcome::Cancelled => self.decision = Some(None),
            Outcome::Ignored if key.code == KeyCode::Char(' ') => {
                if let Some(dir) = self.browser.current_path() {
                    self.choose(dir);
                }
            }
            _ => {}
        }
    }
}

impl View for PathPicker {
    /// The apply key chooses the folder shown.
    fn takes_apply(&self) -> bool {
        true
    }
}

impl Screen for PathPicker {
    fn reports_chords(&self) -> bool {
        true
    }

    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn tick(&mut self) -> bool {
        browser::tick(&mut self.browser)
    }

    fn interrupt(&mut self) {
        self.decision = Some(None);
    }

    fn handle_key(&mut self, key: KeyEvent) {
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
            Edit::Cancelled => self.step = Step::Input,
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let area = frame.area();
        frame.render_widget(header("qc", "choose the directories", ""), area);
        let w = area.width.saturating_sub(4).clamp(20, 90);
        let (title, help, rows) = match &self.step {
            Step::Input => (
                format!(
                    " input: a faba output directory · {}",
                    browser::place(&self.browser)
                ),
                help_line(&[
                    ("type", "find"),
                    ("↑/↓", "move"),
                    ("Enter/→", "open"),
                    ("←", "up"),
                    ("Space", "choose highlighted"),
                    (APPLY_KEYS, "choose this one"),
                    ("Esc", "clear find / cancel"),
                ]),
                area.height.saturating_sub(8).clamp(3, 24),
            ),
            Step::Output { .. } => (
                " output: a new or empty directory ".to_string(),
                help_line(&[("Enter", "next"), ("Esc", "back")]),
                3,
            ),
        };
        // The border takes two columns; the help wraps within the rest.
        let help = crate::tui::wrap_help(help, w.saturating_sub(2));
        let inner = popup_frame(frame, area, w, rows + 3 + help.len() as u16, title);
        let mut lines = match &self.step {
            Step::Input => browser::lines(
                &self.browser,
                rows as usize,
                inner.width as usize,
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
        lines.extend(help);
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
        run_view(&mut picker)?;
    }
    Ok(picker.decision.flatten())
}

#[cfg(test)]
#[path = "tests/path_tui.rs"]
mod tests;
