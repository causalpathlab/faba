//! The file browser every view's pop-ups and panes share: one folder's
//! subfolders, then the files it keeps, with a cursor, narrowed to what has
//! been typed, and typed paths walked (see [`Find`]).

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::{DIM, HIGHLIGHT, PLAIN};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};

use super::{marker, window, FindList};

/// One spelling for a path: absolute, with `.` and `..` resolved by name.
/// Symlinks are kept as given, so a linked file keeps the name it was given.
pub fn normalize(path: &Path) -> PathBuf {
    let abs = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    let mut out = PathBuf::new();
    for c in abs.components() {
        match c {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                out.pop();
            }
            c => out.push(c),
        }
    }
    out
}

/// Whether `name` is an annotation: a GTF or GFF, gzipped or not, any case.
pub fn is_annotation(name: &str) -> bool {
    let name = name.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    [".gtf", ".gff", ".gff3"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

/// One listed name.
#[derive(Debug)]
pub struct Entry {
    pub name: String,
    /// The folder's own path, or the file's; `..` is the parent.
    pub path: PathBuf,
    pub dir: bool,
    /// Marked by the browser's tagger: a faba output folder, a BAM with
    /// its index.
    pub tagged: bool,
}

impl AsRef<str> for Entry {
    fn as_ref(&self) -> &str {
        &self.name
    }
}

/// What a key did in the browser.
pub enum Nav {
    /// It moved, opened a folder or narrowed the listing.
    Moved,
    /// Enter on a file.
    Picked(PathBuf),
    /// Not a browser key; the caller's to handle.
    Ignored,
}

/// The browser's state, independent of the terminal so it can be tested.
pub struct Browser {
    pub cwd: PathBuf,
    /// The listing: `..` first unless at the root, then subfolders, then
    /// files, each sorted, narrowed by what has been typed.
    pub list: FindList<Entry>,
    /// The files marked, in the order marked, from any folder; `None` when
    /// this browser does not mark.
    pub marked: Option<Vec<PathBuf>>,
    keep: Box<dyn Fn(&str) -> bool>,
    tag: Box<dyn FnMut(&Path) -> bool>,
}

impl Browser {
    /// A browser at `cwd` that has not listed it yet: files `keep` accepts
    /// by name (folders are always listed), the entries `tag` marks.
    pub fn new(
        cwd: PathBuf,
        keep: impl Fn(&str) -> bool + 'static,
        tag: impl FnMut(&Path) -> bool + 'static,
    ) -> Self {
        Self {
            cwd,
            list: FindList::default(),
            marked: None,
            keep: Box::new(keep),
            tag: Box::new(tag),
        }
    }

    /// Marking files, with none marked yet.
    pub fn marking(mut self) -> Self {
        self.marked = Some(Vec::new());
        self
    }

    /// Mark `file`, or unmark it if it was; nothing when not marking.
    pub fn toggle_mark(&mut self, file: PathBuf) {
        let Some(marked) = &mut self.marked else {
            return;
        };
        match marked.iter().position(|m| *m == file) {
            Some(i) => {
                marked.remove(i);
            }
            None => marked.push(file),
        }
    }

    /// Listed at its folder.
    pub fn opened(mut self) -> Self {
        self.open(self.cwd.clone());
        self
    }

    /// Show `dir`'s visible entries, with nothing typed; going up lands on
    /// the folder just left.
    pub fn open(&mut self, dir: PathBuf) {
        let dir = normalize(&dir);
        let (mut dirs, mut files): (Vec<String>, Vec<String>) = (Vec::new(), Vec::new());
        for e in std::fs::read_dir(&dir).into_iter().flatten().flatten() {
            let name = e.file_name().to_string_lossy().into_owned();
            if name.starts_with('.') {
                continue;
            }
            // A zarr store is a folder, but it is one matrix.
            let is_dir = e.file_type().is_ok_and(|t| t.is_dir()) || e.path().is_dir();
            if is_dir && !name.ends_with(".zarr") {
                dirs.push(name);
            } else if (self.keep)(&name) {
                files.push(name);
            }
        }
        dirs.sort();
        files.sort();
        let mut all: Vec<Entry> = dir
            .parent()
            .map(|p| Entry {
                name: "..".into(),
                path: p.to_path_buf(),
                dir: true,
                tagged: false,
            })
            .into_iter()
            .collect();
        for (names, is_dir) in [(dirs, true), (files, false)] {
            for name in names {
                let path = dir.join(&name);
                all.push(Entry {
                    tagged: (self.tag)(&path),
                    name,
                    path,
                    dir: is_dir,
                });
            }
        }
        let left = (Some(dir.as_path()) == self.cwd.parent())
            .then(|| self.cwd.file_name())
            .flatten()
            .map(|n| n.to_string_lossy().into_owned());
        self.list.set(all, 0);
        if let Some(l) = left {
            self.list.select(&l);
        }
        self.cwd = dir;
    }

    pub fn up(&mut self) {
        if let Some(p) = self.cwd.parent().map(Path::to_path_buf) {
            self.open(p);
        }
    }

    /// The entry under the cursor.
    pub fn highlighted(&self) -> Option<&Entry> {
        self.list.highlighted()
    }

    /// Walk a typed path, narrow by typed text, move, open a folder or go
    /// up; Enter on a file picks it. Letters type, so only arrows move.
    pub fn key(&mut self, key: KeyEvent) -> Nav {
        let into = self
            .highlighted()
            .filter(|e| e.dir && !self.list.find.text.is_empty())
            .map(|e| e.path.clone());
        if let Some(to) = self.list.find.jump(&key, &self.cwd, into) {
            self.open(to);
            return Nav::Moved;
        }
        if self.list.key(&key) {
            return Nav::Moved;
        }
        match key.code {
            KeyCode::Enter | KeyCode::Right => match self.highlighted() {
                Some(e) if e.dir => self.open(e.path.clone()),
                Some(e) if key.code == KeyCode::Enter => return Nav::Picked(e.path.clone()),
                _ => {}
            },
            KeyCode::Left | KeyCode::Backspace => self.up(),
            _ => return Nav::Ignored,
        }
        Nav::Moved
    }

    /// `rows` lines of the listing around the cursor, `width` wide; marked
    /// entries carry `tag`, and an empty listing says `empty`.
    /// With marking on, each file has a `[x]`/`[ ]` box.
    pub fn lines(&self, rows: usize, width: usize, tag: &str, empty: &str) -> Vec<Line<'static>> {
        if self.list.is_empty() {
            let why = self.list.empty_note(empty);
            return vec![Line::from(Span::styled(format!("  {why}"), DIM))];
        }
        let lines = self
            .list
            .shown()
            .enumerate()
            .map(|(i, e)| {
                let selected = i == self.list.at;
                let tag = if e.tagged {
                    format!("  {tag}")
                } else {
                    String::new()
                };
                let check = match &self.marked {
                    None => "",
                    Some(_) if e.dir => "    ",
                    Some(m) if m.contains(&e.path) => "[x] ",
                    Some(_) => "[ ] ",
                };
                let name_w = width.saturating_sub(tag.len() + check.len() + 3);
                let name = if e.dir {
                    format!("{}/", e.name)
                } else {
                    e.name.clone()
                };
                let style = match (selected, e.dir) {
                    (true, _) => HIGHLIGHT,
                    (false, true) => DIM,
                    (false, false) => PLAIN,
                };
                Line::from(vec![
                    marker(selected),
                    Span::styled(check, HIGHLIGHT),
                    Span::styled(format!("{name:<name_w$.name_w$}"), style),
                    Span::styled(tag, DIM),
                ])
            })
            .collect();
        window(lines, self.list.at, rows)
    }
}

#[cfg(test)]
#[path = "tests/browser.rs"]
mod tests;
