//! The pop-up `pileup` and `metagene` open when the command line names no
//! input: browse, mark the files to read together, or choose a whole faba
//! output folder. It opens where the last choice was made.

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::{header, help_line, Screen, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::qc::layout::looks_like_faba_dir;
use crate::qc::path_tui::faba_browser;
use crate::tui::browser::{Browser, Nav};
use crate::tui::{is_apply, popup_frame, run_view, tilde, View, APPLY_KEYS};

/// What was chosen: the marked files, or one folder.
#[derive(Debug, PartialEq)]
pub struct Chosen {
    pub paths: Vec<PathBuf>,
    /// Marked files are to be one track each rather than one together.
    pub separate: bool,
}

/// The choice, independent of the terminal so it can be tested.
pub struct InputPicker {
    /// Marking files, in the order marked, from any folder.
    browser: Browser,
    /// The command asking, for the header.
    command: &'static str,
    /// What the listed files are, for the title.
    what: &'static str,
    /// Whether marked files read as one track each (Tab); `None` when the
    /// command reads them only together.
    separate: Option<bool>,
    /// Why the last choice was refused.
    pub error: Option<String>,
    /// The choice; `Some(None)` when cancelled.
    pub decision: Option<Option<Chosen>>,
}

impl InputPicker {
    /// Listing `cwd`'s folders and the files `keep` accepts by name, with
    /// Tab to read marked files apart when `separable`.
    pub fn new(
        cwd: PathBuf,
        command: &'static str,
        what: &'static str,
        keep: fn(&str) -> bool,
        separable: bool,
    ) -> Self {
        let mut browser = faba_browser(cwd.clone(), keep).marking();
        browser.open(cwd);
        Self {
            browser,
            command,
            what,
            separate: separable.then_some(false),
            error: None,
            decision: None,
        }
    }

    /// Open at `last`'s folder with `last` highlighted.
    fn start_at(&mut self, last: &Path) {
        let (Some(dir), Some(name)) = (last.parent(), last.file_name()) else {
            return;
        };
        if dir.is_dir() {
            self.browser.open(dir.to_path_buf());
            self.browser.list.select(&name.to_string_lossy());
        }
    }

    fn finish(&mut self, paths: Vec<PathBuf>) {
        let separate = self.separate.unwrap_or(false);
        self.decision = Some(Some(Chosen { paths, separate }));
    }

    /// Take the whole of `dir`, if it is a faba output folder.
    fn choose_dir(&mut self, dir: PathBuf) {
        if looks_like_faba_dir(&dir) {
            self.finish(vec![dir]);
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
            return self.finish(marked);
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
            KeyCode::Tab => {
                if let Some(s) = &mut self.separate {
                    *s = !*s;
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
        let marked = match (n_marked, self.separate) {
            (0, _) => String::new(),
            (n, Some(true)) => format!("· {n} marked, a track each "),
            (n, Some(false)) => format!("· {n} marked, one track "),
            (n, None) => format!("· {n} marked "),
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
        let mut keys = vec![
            ("type", "find"),
            ("↑/↓", "move"),
            ("Enter/→", "open"),
            ("←", "up"),
            ("Space", "mark file / choose folder"),
        ];
        match self.separate {
            Some(true) => keys.push(("Tab", "one track")),
            Some(false) => keys.push(("Tab", "a track each")),
            None => {}
        }
        keys.extend([(APPLY_KEYS, apply), ("Esc", "clear find / cancel")]);
        lines.push(help_line(&keys));
        frame.render_widget(Paragraph::new(lines), inner);
    }
}

/// Where the last choice is kept, so the picker opens there next time.
fn last_choice_file() -> Option<PathBuf> {
    Some(crate::tui::home()?.join(".cache/faba/last_input"))
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
    separable: bool,
) -> anyhow::Result<Option<(Vec<Box<str>>, bool)>> {
    anyhow::ensure!(
        !batch && data_beans::interactive::tui_available(),
        "name {what}, or a faba output directory"
    );
    let mut picker = InputPicker::new(std::env::current_dir()?, command, what, keep, separable);
    let memo = last_choice_file();
    if let Some(last) = memo.as_ref().and_then(|m| std::fs::read_to_string(m).ok()) {
        picker.start_at(Path::new(last.trim()));
    }
    run_view(&mut picker)?;
    let Some(chosen) = picker.decision.flatten() else {
        return Ok(None);
    };
    // Remembering is a convenience: a failure to write is no failure.
    if let (Some(memo), Some(first)) = (memo, chosen.paths.first()) {
        if let Some(dir) = memo.parent() {
            let _ = std::fs::create_dir_all(dir);
        }
        let _ = std::fs::write(memo, first.to_string_lossy().as_bytes());
    }
    let paths = chosen
        .paths
        .iter()
        .map(|p| p.to_string_lossy().into_owned().into_boxed_str())
        .collect();
    Ok(Some((paths, chosen.separate)))
}

#[cfg(test)]
#[path = "tests/input_picker.rs"]
mod tests;
