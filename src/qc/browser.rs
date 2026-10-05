//! A directory browser for the `faba qc` pop-ups: one directory's
//! subdirectories, then the files a caller asks for, with a cursor.

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::{DIM, HIGHLIGHT, PLAIN};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};

use super::widgets::first_visible;

/// One listed name.
pub struct Entry {
    pub name: String,
    pub dir: bool,
    /// Marked by the caller's `tag`, e.g. a faba output directory.
    pub tagged: bool,
}

/// What to list and mark: files `keep` accepts (directories are always
/// listed), and the entries `tag` marks.
pub struct Listing<'a> {
    pub keep: &'a dyn Fn(&str) -> bool,
    pub tag: &'a mut dyn FnMut(&Path) -> bool,
}

/// What a key did in the browser.
pub enum Nav {
    /// It moved or opened a directory.
    Moved,
    /// Enter on a file.
    Picked(PathBuf),
    /// Not a browser key; the caller's to handle.
    Ignored,
}

/// The browser's state, independent of the terminal so it can be tested.
pub struct Browser {
    pub cwd: PathBuf,
    /// `cwd`'s listing: `..` first unless at the root, then subdirectories,
    /// then files, each sorted.
    pub entries: Vec<Entry>,
    pub at: usize,
}

impl Browser {
    /// A browser at `cwd` that has not listed it yet.
    pub fn new(cwd: PathBuf) -> Self {
        Self {
            cwd,
            entries: Vec::new(),
            at: 0,
        }
    }

    /// Show `dir`'s visible entries; going up lands on the directory just
    /// left.
    pub fn open(&mut self, dir: PathBuf, listing: &mut Listing) {
        let (mut dirs, mut files): (Vec<String>, Vec<String>) = (Vec::new(), Vec::new());
        for e in std::fs::read_dir(&dir).into_iter().flatten().flatten() {
            let name = e.file_name().to_string_lossy().into_owned();
            if name.starts_with('.') {
                continue;
            }
            // A zarr store is a directory, but it is one matrix.
            let is_dir = e.file_type().is_ok_and(|t| t.is_dir()) || e.path().is_dir();
            if is_dir && !name.ends_with(".zarr") {
                dirs.push(name);
            } else if !is_dir && (listing.keep)(&name) {
                files.push(name);
            }
        }
        dirs.sort();
        files.sort();
        if dir.parent().is_some() {
            dirs.insert(0, "..".into());
        }
        let mut entries: Vec<Entry> = dirs
            .into_iter()
            .map(|name| Entry {
                tagged: name != ".." && (listing.tag)(&dir.join(&name)),
                name,
                dir: true,
            })
            .collect();
        entries.extend(files.into_iter().map(|name| Entry {
            tagged: (listing.tag)(&dir.join(&name)),
            name,
            dir: false,
        }));
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

    /// The path under the cursor, and whether it is a directory.
    pub fn highlighted(&self) -> Option<(PathBuf, bool)> {
        let e = self.entries.get(self.at)?;
        let path = match e.name.as_str() {
            ".." => self.cwd.parent()?.to_path_buf(),
            name => self.cwd.join(name),
        };
        Some((path, e.dir))
    }

    fn step(&mut self, delta: isize) {
        let n = self.entries.len() as isize;
        if n > 0 {
            self.at = (self.at as isize + delta).clamp(0, n - 1) as usize;
        }
    }

    /// Move, open a directory or go up; Enter on a file picks it.
    pub fn key(&mut self, key: KeyEvent, listing: &mut Listing) -> Nav {
        match key.code {
            KeyCode::Up | KeyCode::Char('k') => self.step(-1),
            KeyCode::Down | KeyCode::Char('j') => self.step(1),
            KeyCode::PageUp => self.step(-10),
            KeyCode::PageDown => self.step(10),
            KeyCode::Enter | KeyCode::Right | KeyCode::Char('l') => match self.highlighted() {
                Some((dir, true)) => self.open(dir, listing),
                Some((file, false)) if key.code == KeyCode::Enter => return Nav::Picked(file),
                _ => {}
            },
            KeyCode::Left | KeyCode::Backspace | KeyCode::Char('h') => {
                if let Some(up) = self.cwd.parent().map(Path::to_path_buf) {
                    self.open(up, listing);
                }
            }
            _ => return Nav::Ignored,
        }
        Nav::Moved
    }

    /// `rows` lines of the listing around the cursor, `width` wide; tagged
    /// entries carry `tag`, and an empty listing says `empty`.
    pub fn lines(&self, rows: usize, width: usize, tag: &str, empty: &str) -> Vec<Line<'static>> {
        let first = first_visible(self.at, self.entries.len(), rows);
        let mut lines: Vec<Line> = self
            .entries
            .iter()
            .enumerate()
            .skip(first)
            .take(rows)
            .map(|(i, e)| {
                let selected = i == self.at;
                let tag = if e.tagged {
                    format!("  {tag}")
                } else {
                    String::new()
                };
                let name_w = width.saturating_sub(tag.len() + 3);
                let name = if e.dir {
                    format!("{}/", e.name)
                } else {
                    e.name.clone()
                };
                Line::from(vec![
                    Span::styled(if selected { "▸ " } else { "  " }, HIGHLIGHT),
                    Span::styled(
                        format!("{name:<name_w$.name_w$}"),
                        if selected { HIGHLIGHT } else { PLAIN },
                    ),
                    Span::styled(tag, DIM),
                ])
            })
            .collect();
        if self.entries.is_empty() {
            lines.push(Line::from(Span::styled(format!("  {empty}"), DIM)));
        }
        lines
    }
}

#[cfg(test)]
#[path = "tests/browser.rs"]
mod tests;
