//! The figures the views have saved: a log of each PDF with a small
//! thumbnail, kept in `.faba-view/` where faba runs, so the list lasts from
//! one session to the next. The save prompt shows it beside the name.

use image::RgbaImage;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

/// Width of a thumbnail on disk, in pixels.
const THUMB_WIDTH: u32 = 240;
/// Saves remembered; older ones go, with their thumbnails.
const KEEP: usize = 200;

/// One saved figure.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Entry {
    /// The PDF, as an absolute path.
    pub path: PathBuf,
    /// The view it came from.
    pub what: String,
    /// When, in seconds since the Unix epoch.
    pub when: u64,
    /// The thumbnail's file name in `thumbs/`.
    thumb: String,
}

impl Entry {
    /// The file's name, without its folder.
    pub fn name(&self) -> String {
        self.path
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned())
    }
}

/// Every save remembered in one directory, newest first.
pub struct Gallery {
    dir: PathBuf,
    entries: Vec<Entry>,
}

impl Gallery {
    /// The saves logged in `dir`, less those whose file or thumbnail is
    /// gone. An unreadable log starts an empty gallery.
    pub fn open(dir: &Path) -> Self {
        let mut g = Self {
            dir: dir.into(),
            entries: std::fs::read(dir.join("saved.json"))
                .ok()
                .and_then(|b| serde_json::from_slice(&b).ok())
                .unwrap_or_default(),
        };
        let thumbs = g.thumbs();
        g.entries
            .retain(|e| e.path.exists() && thumbs.join(&e.thumb).exists());
        g
    }

    fn thumbs(&self) -> PathBuf {
        self.dir.join("thumbs")
    }

    /// The gallery of the directory faba runs in; under test, one apart.
    pub fn here() -> Self {
        if cfg!(test) {
            return Self::open(
                &std::env::temp_dir().join(format!("faba-view-{}", std::process::id())),
            );
        }
        Self::open(Path::new(".faba-view"))
    }

    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// Where `e`'s thumbnail is.
    pub fn thumb(&self, e: &Entry) -> PathBuf {
        self.thumbs().join(&e.thumb)
    }

    /// Remember that `path` was just saved from `what`, showing `picture`:
    /// at the top, replacing an earlier save of the same file.
    pub fn add(&mut self, path: &Path, what: &str, picture: &RgbaImage) -> anyhow::Result<()> {
        let path = path.canonicalize()?;
        if let Some(i) = self.entries.iter().position(|e| e.path == path) {
            let old = self.entries.remove(i);
            let _ = std::fs::remove_file(self.thumb(&old));
        }
        let now = SystemTime::now().duration_since(UNIX_EPOCH)?;
        let thumbs = self.thumbs();
        std::fs::create_dir_all(&thumbs)?;
        let thumb = format!("{}.png", now.as_nanos());
        thumbnail(picture).save(thumbs.join(&thumb))?;
        self.entries.insert(
            0,
            Entry {
                path,
                what: what.into(),
                when: now.as_secs(),
                thumb,
            },
        );
        for old in self.entries.split_off(self.entries.len().min(KEEP)) {
            let _ = std::fs::remove_file(self.thumb(&old));
        }
        std::fs::write(
            self.dir.join("saved.json"),
            serde_json::to_string_pretty(&self.entries)?,
        )?;
        Ok(())
    }
}

/// `picture` shrunk to [`THUMB_WIDTH`] wide at its aspect.
fn thumbnail(picture: &RgbaImage) -> RgbaImage {
    let (w, h) = picture.dimensions();
    let th = (THUMB_WIDTH as f32 * h as f32 / w.max(1) as f32)
        .round()
        .max(1.0) as u32;
    image::imageops::resize(
        picture,
        THUMB_WIDTH,
        th,
        image::imageops::FilterType::Triangle,
    )
}

/// How long ago `when` was, at `now` (both seconds since the epoch).
pub fn ago(when: u64, now: u64) -> String {
    match now.saturating_sub(when) {
        s if s < 60 => "just now".into(),
        s if s < 3600 => format!("{} min ago", s / 60),
        s if s < 86_400 => format!("{} h ago", s / 3600),
        s => format!("{} d ago", s / 86_400),
    }
}

#[cfg(test)]
#[path = "tests/gallery.rs"]
mod tests;
