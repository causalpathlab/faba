//! The pop-up `pileup` and `metagene` open when the command line names no
//! input: browse, mark the files to read together, or choose a whole faba
//! output folder. It opens where the last choice was made from the same
//! working folder.

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
    pub paths: Vec<Box<str>>,
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
    /// Listing the folders and the files `keep` accepts by name, at
    /// `last`'s folder with it highlighted when that is still there, else
    /// at `cwd`; with Tab to read marked files apart when `separable`.
    pub fn new(
        cwd: PathBuf,
        last: Option<&Path>,
        command: &'static str,
        what: &'static str,
        keep: fn(&str) -> bool,
        separable: bool,
    ) -> Self {
        let at = last.and_then(|l| Some((l.parent().filter(|d| d.is_dir())?, l.file_name()?)));
        let start = at.map_or(cwd, |(dir, _)| dir.to_path_buf());
        let mut browser = faba_browser(start.clone(), keep).marking();
        browser.open(start);
        if let Some((_, name)) = at {
            browser.list.select(&name.to_string_lossy());
        }
        Self {
            browser,
            command,
            what,
            separate: separable.then_some(false),
            error: None,
            decision: None,
        }
    }

    fn finish(&mut self, paths: Vec<PathBuf>) {
        let paths = paths
            .iter()
            .map(|p| p.to_string_lossy().into_owned().into_boxed_str())
            .collect();
        let separate = self.separate.unwrap_or(false);
        self.decision = Some(Some(Chosen { paths, separate }));
    }

    /// The marked count and, when Tab applies, how they read now and what
    /// Tab makes of them.
    fn track_mode(&self) -> Option<(&'static str, &'static str)> {
        self.separate.map(|apart| match apart {
            true => ("a track each", "one track"),
            false => ("one track", "a track each"),
        })
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
        let marked = match (n_marked, self.track_mode()) {
            (0, _) => String::new(),
            (n, Some((now, _))) => format!("· {n} marked, {now} "),
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
        if let Some((_, tab)) = self.track_mode() {
            keys.push(("Tab", tab));
        }
        keys.extend([(APPLY_KEYS, apply), ("Esc", "clear find / cancel")]);
        lines.push(help_line(&keys));
        frame.render_widget(Paragraph::new(lines), inner);
    }
}

/// Where the last choices are kept, one line per working folder
/// (`folder<TAB>choice`), so the picker opens there next time.
fn memo_file() -> Option<PathBuf> {
    let cache = std::env::var_os("XDG_CACHE_HOME")
        .map(PathBuf::from)
        .filter(|p| p.is_absolute())
        .or_else(|| Some(crate::tui::home()?.join(".cache")))?;
    Some(cache.join("faba/last_input"))
}

/// The last choice made from `cwd`.
fn recall(memo: &Path, cwd: &Path) -> Option<PathBuf> {
    let text = std::fs::read_to_string(memo).ok()?;
    let cwd = cwd.to_string_lossy();
    text.lines()
        .find_map(|l| l.split_once('\t').filter(|(from, _)| *from == cwd))
        .map(|(_, last)| PathBuf::from(last))
}

/// Keep `choice` as the last made from `cwd`, the most recent folders
/// first and at most [`MEMO_LINES`] of them. A convenience: a failure to
/// write is no failure.
fn remember(memo: &Path, cwd: &Path, choice: &str) {
    let cwd = cwd.to_string_lossy();
    let old = std::fs::read_to_string(memo).unwrap_or_default();
    let others = old
        .lines()
        .filter(|l| l.split_once('\t').is_some_and(|(from, _)| from != cwd));
    let lines: Vec<String> = std::iter::once(format!("{cwd}\t{choice}"))
        .chain(others.map(String::from))
        .take(MEMO_LINES)
        .collect();
    if let Some(dir) = memo.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let _ = std::fs::write(memo, lines.join("\n") + "\n");
}

/// Working folders whose last choice is kept.
const MEMO_LINES: usize = 64;

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
) -> anyhow::Result<Option<Chosen>> {
    anyhow::ensure!(
        !batch && data_beans::interactive::tui_available(),
        "name {what}, or a faba output directory"
    );
    let cwd = std::env::current_dir()?;
    let memo = memo_file();
    let last = memo.as_deref().and_then(|m| recall(m, &cwd));
    let mut picker = InputPicker::new(cwd.clone(), last.as_deref(), command, what, keep, separable);
    run_view(&mut picker)?;
    let chosen = picker.decision.flatten();
    if let (Some(memo), Some(first)) = (&memo, chosen.as_ref().and_then(|c| c.paths.first())) {
        remember(memo, &cwd, first);
    }
    Ok(chosen)
}

#[cfg(test)]
#[path = "tests/input_picker.rs"]
mod tests;
