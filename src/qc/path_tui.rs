//! Pop-up asking `faba qc` for the input and output directories when the
//! command line left them out: browse to a faba output directory, then name
//! a new one to write to.

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::{header, help_line, input_line, Screen, DIM, HIGHLIGHT, PLAIN};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;

use super::layout::looks_like_faba_dir;
use crate::figure::{Edit, LineInput};
use crate::tui::{first_visible, popup_frame};

/// One subdirectory in the browser.
struct Entry {
    name: String,
    /// It holds faba matrices or site tables.
    faba: bool,
}

enum Step {
    /// Browsing for the input directory.
    Input,
    /// Naming the output directory for the chosen `input`.
    Output { input: PathBuf, line: LineInput },
}

/// The two questions, independent of the terminal so they can be tested.
struct PathPicker {
    step: Step,
    cwd: PathBuf,
    /// `cwd`'s subdirectories; listed when the browser first shows.
    entries: Vec<Entry>,
    at: usize,
    /// Directories already checked for faba files.
    faba: FxHashMap<PathBuf, bool>,
    /// Why the last choice was refused.
    error: Option<String>,
    /// The output as given on the command line, if it was.
    given_output: Option<String>,
    decision: Option<Option<(String, String)>>,
}

/// A directory next to `input` named after it, that does not exist yet.
fn suggest_output(input: &Path) -> String {
    let base = input.to_string_lossy();
    let base = base.trim_end_matches('/');
    let mut out = format!("{base}_qc");
    let mut n = 2;
    while Path::new(&out).exists() {
        out = format!("{base}_qc{n}");
        n += 1;
    }
    out
}

/// Why `out` cannot take the filtered fileset, if it cannot.
fn output_problem(out: &str) -> Option<String> {
    let path = Path::new(out);
    if out.is_empty() {
        return Some("name a directory".into());
    }
    if path.is_file() {
        return Some(format!("{out} is a file"));
    }
    let non_empty = std::fs::read_dir(path).is_ok_and(|mut d| d.next().is_some());
    non_empty.then(|| format!("{out} already contains files; choose an empty one"))
}

impl PathPicker {
    fn new(cwd: PathBuf, input: Option<PathBuf>, given_output: Option<String>) -> Self {
        let mut p = Self {
            step: Step::Input,
            cwd: cwd.clone(),
            entries: Vec::new(),
            at: 0,
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

    fn is_faba(&mut self, dir: &Path) -> bool {
        *self
            .faba
            .entry(dir.to_path_buf())
            .or_insert_with(|| looks_like_faba_dir(dir))
    }

    /// Show `dir`'s visible subdirectories, sorted, with `..` first unless at
    /// the root; going up lands on the directory just left.
    fn open(&mut self, dir: PathBuf) {
        let mut names: Vec<String> = std::fs::read_dir(&dir)
            .map(|rd| {
                rd.flatten()
                    .filter(|e| e.file_type().is_ok_and(|t| t.is_dir()) || e.path().is_dir())
                    .map(|e| e.file_name().to_string_lossy().into_owned())
                    // A zarr store is a directory, but it is one matrix.
                    .filter(|n| !n.starts_with('.') && !n.ends_with(".zarr"))
                    .collect()
            })
            .unwrap_or_default();
        names.sort();
        if dir.parent().is_some() {
            names.insert(0, "..".into());
        }
        let entries = names
            .into_iter()
            .map(|name| Entry {
                faba: name != ".." && self.is_faba(&dir.join(&name)),
                name,
            })
            .collect();
        let left = (Some(dir.as_path()) == self.cwd.parent())
            .then(|| self.cwd.file_name())
            .flatten()
            .map(|n| n.to_string_lossy().into_owned());
        self.entries = entries;
        self.at = left
            .and_then(|l| self.entries.iter().position(|e| e.name == l))
            .unwrap_or(0);
        self.cwd = dir;
    }

    fn highlighted(&self) -> Option<PathBuf> {
        let e = self.entries.get(self.at)?;
        Some(match e.name.as_str() {
            ".." => self.cwd.parent()?.to_path_buf(),
            name => self.cwd.join(name),
        })
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

    fn step_at(&mut self, delta: isize) {
        let n = self.entries.len() as isize;
        if n > 0 {
            self.at = (self.at as isize + delta).clamp(0, n - 1) as usize;
        }
    }

    fn browse_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Up | KeyCode::Char('k') => self.step_at(-1),
            KeyCode::Down | KeyCode::Char('j') => self.step_at(1),
            KeyCode::PageUp => self.step_at(-10),
            KeyCode::PageDown => self.step_at(10),
            KeyCode::Enter | KeyCode::Right | KeyCode::Char('l') => {
                if let Some(dir) = self.highlighted() {
                    self.open(dir);
                }
            }
            KeyCode::Left | KeyCode::Backspace | KeyCode::Char('h') => {
                if let Some(up) = self.cwd.parent().map(Path::to_path_buf) {
                    self.open(up);
                }
            }
            KeyCode::Char(' ') => {
                if let Some(dir) = self.highlighted() {
                    self.choose(dir);
                }
            }
            KeyCode::Char('.') => self.choose(self.cwd.clone()),
            KeyCode::Esc | KeyCode::Char('q') => self.decision = Some(None),
            _ => {}
        }
    }

    fn lines(&self, rows: usize, width: usize) -> Vec<Line<'static>> {
        let first = first_visible(self.at, self.entries.len(), rows);
        let mut lines: Vec<Line> = self
            .entries
            .iter()
            .enumerate()
            .skip(first)
            .take(rows)
            .map(|(i, e)| {
                let selected = i == self.at;
                let tag = if e.faba { "  faba output" } else { "" };
                let name_w = width.saturating_sub(tag.len() + 3);
                Line::from(vec![
                    Span::styled(if selected { "▸ " } else { "  " }, HIGHLIGHT),
                    Span::styled(
                        format!("{:<name_w$.name_w$}", format!("{}/", e.name)),
                        if selected { HIGHLIGHT } else { PLAIN },
                    ),
                    Span::styled(tag, DIM),
                ])
            })
            .collect();
        if self.entries.is_empty() {
            lines.push(Line::from(Span::styled("  no subdirectory", DIM)));
        }
        lines
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
                if self.entries.is_empty() {
                    self.open(self.cwd.clone());
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
                format!(" input: a faba output directory · {} ", self.cwd.display()),
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
            Step::Input => self.lines(rows as usize, inner.width as usize),
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
