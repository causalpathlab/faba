//! The file browser every view's pop-ups and panes share: legume-tui's
//! [`legume_tui::browse::Browser`] for a [`Simple`] want (one folder's
//! subfolders, then the files it keeps, narrowed to what has been typed,
//! typed paths walked, listed on a worker thread), with what faba adds:
//! which names are annotations or genomes, the listing's wait, and the
//! rows drawn with what the listing says under them.

use std::borrow::Cow;
use std::path::Path;
use std::time::Duration;

use ratatui::text::{Line, Span};

pub use legume_tui::browse::{normalize, Entry, Header, Outcome, Simple, NOTHING_MATCHES};

/// faba's browser: files kept by name, each with what a noter says of it.
pub type Browser = legume_tui::browse::Browser<Simple>;

/// How long opening a folder waits for its listing before showing it is
/// still listing: legume-tui's wait, or, in faba's tests, long enough for
/// any listing to finish, so no test races the worker thread.
#[cfg(not(test))]
const GRACE: Duration = legume_tui::browse::LIST_GRACE;
#[cfg(test)]
const GRACE: Duration = Duration::from_secs(60);

/// A browser at `dir` for `want`, listed (within [`GRACE`]), the cursor on
/// `select` when that is listed.
pub fn open(dir: impl AsRef<Path>, want: Simple, select: Option<&str>) -> Browser {
    Browser::open_within(dir, want, select, GRACE)
}

/// A [`Simple`] want titled `title` for `what`, keeping the files `keep`
/// takes by name, each with what `note` says of it (told whether it speaks
/// of a folder).
pub fn want(
    title: &'static str,
    what: &'static str,
    keep: impl Fn(&str) -> bool + Send + Sync + 'static,
    note: impl FnMut(&Path, bool) -> Option<Cow<'static, str>> + Send + 'static,
) -> Simple {
    Simple::new(Header::new(title, what), keep, note)
}

/// [`want`] with nothing said of the files.
pub fn plain(
    title: &'static str,
    what: &'static str,
    keep: impl Fn(&str) -> bool + Send + Sync + 'static,
) -> Simple {
    Simple::keeping(Header::new(title, what), keep)
}

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

/// Take what the folder's listing found; whether to redraw: something
/// came, or it is still listing (its spinner and seconds turn).
pub fn tick(b: &mut Browser) -> bool {
    b.poll() || b.listing().is_some()
}

/// Where a pop-up's title says the browser is: its folder, `~` for home,
/// then what was typed.
pub fn place(b: &Browser) -> String {
    format!("{} {}", super::tilde(&b.dir), b.tag())
}

/// What the listing says under its rows, if anything: that the folder is
/// still being listed, that nothing matches what was typed, or `empty`
/// when it holds nothing but `..`.
pub fn note(b: &Browser, empty: Option<&str>) -> Option<Cow<'static, str>> {
    if b.listing().is_some() {
        return b.status();
    }
    if b.shown().iter().any(|e| !matches!(e, Entry::Up)) {
        return None;
    }
    if !b.filter.is_empty() {
        return Some(Cow::Borrowed(NOTHING_MATCHES));
    }
    empty.map(|e| Cow::Owned(e.to_string()))
}

/// `rows` lines of the listing around the cursor, `width` wide, in faba's
/// look ([`super::LOOK`]); the last is what [`note`] says, when it says
/// anything.
pub fn lines(b: &Browser, rows: usize, width: usize, empty: &str) -> Vec<Line<'static>> {
    let note = note(b, Some(empty));
    let room = rows.saturating_sub(usize::from(note.is_some())).max(1);
    let mut out = b.rows(room, width, &super::LOOK);
    if let Some(n) = note {
        out.push(Line::from(Span::styled(format!("  {n}"), super::LOOK.hint)));
    }
    out
}
