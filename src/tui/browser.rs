//! The file browser every view's pop-ups and panes share: one folder's
//! subfolders, then the files it keeps, with a cursor, narrowed to what has
//! been typed, and typed paths walked: legume-tui's
//! [`legume_tui::browse::Browser`], noting folders as well as files, with
//! the files marked kept here.

use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};

use ratatui::text::Line;

pub use legume_tui::browse::{normalize, Entry, Nav};

/// Whether `name` is an annotation: a GTF or GFF, gzipped or not, any case.
pub fn is_annotation(name: &str) -> bool {
    let name = name.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    [".gtf", ".gff", ".gff3"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

/// Whether `name` is a genome FASTA, gzipped or not, any case.
pub fn is_fasta(name: &str) -> bool {
    let name = name.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    [".fa", ".fasta", ".fna"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

/// The browser's state, independent of the terminal so it can be tested:
/// legume-tui's (its folder, listing, keys), and the files marked.
pub struct Browser {
    shared: legume_tui::browse::Browser,
    /// The files marked, in the order marked, from any folder; `None` when
    /// this browser does not mark.
    pub marked: Option<Vec<PathBuf>>,
}

impl Browser {
    /// A browser at `cwd` that has not listed it yet: files `keep` accepts
    /// by name (folders are always listed), each file with what `note` says
    /// of it (a BAM with its index).
    pub fn new(
        cwd: PathBuf,
        keep: impl Fn(&str) -> bool + 'static,
        note: impl FnMut(&Path) -> Option<String> + 'static,
    ) -> Self {
        Self {
            shared: legume_tui::browse::Browser::new(cwd, keep, note),
            marked: None,
        }
    }

    /// Noting folders as well as files (see
    /// [`legume_tui::browse::Browser::noting_folders`]).
    pub fn noting_folders(mut self) -> Self {
        self.shared = self.shared.noting_folders();
        self
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
        self.shared = self.shared.opened();
        self
    }

    /// `rows` lines of the listing around the cursor, `width` wide, each
    /// entry's note after it; an empty listing says `empty`. With marking
    /// on, each file has a `[x]`/`[ ]` box.
    pub fn lines(&self, rows: usize, width: usize, empty: &str) -> Vec<Line<'static>> {
        let mark_width = if self.marked.is_some() { 3 } else { 0 };
        let mark = |e: &Entry| {
            let on = self.marked.as_ref()?.contains(&e.path);
            Some(if on { "[x]" } else { "[ ]" }.to_string())
        };
        self.shared
            .lines(rows, width, &super::LOOK, mark_width, mark, empty)
    }
}

impl Deref for Browser {
    type Target = legume_tui::browse::Browser;

    fn deref(&self) -> &Self::Target {
        &self.shared
    }
}

impl DerefMut for Browser {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.shared
    }
}
