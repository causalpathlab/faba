//! Self-contained HTML reports of the full-screen views, for sharing or
//! posting (a GitHub page, a lab wiki): one file, nothing loaded from
//! anywhere, its data embedded as JSON.
//!
//! Every view's report is a page fragment of its own (its panes, styles and
//! drawing code) inside one shared shell (`report/shell.html`): the theme,
//! the window with the faba mark, the key footer, the legend, the canvas
//! helpers, and the data. Each full-screen view names its report key in
//! its footer: `e` in the modules browser, `Ctrl+E` in the QC site view
//! (where `e` types a value), a chord it takes through
//! [`crate::tui::View::takes_chord`]. A view that can run without a
//! terminal also writes it by a `--html FILE` flag. Either way the page is
//! written on a thread of its own ([`Writing`]).

use crate::figure::logo::{BEAN, TONES};
use ratatui::style::Color;
use serde::Serialize;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;

/// The page around every report; see the module docs.
const SHELL: &str = include_str!("report/shell.html");

/// The faba mark, as the views draw it in their top-right corner: the
/// terminal mark's beans, each in its own tone, the name under the last.
fn mark() -> String {
    let hex = |tone: Color| match tone {
        Color::Rgb(r, g, b) => format!("#{r:02x}{g:02x}{b:02x}"),
        _ => unreachable!("the mark's tones are RGB"),
    };
    let row = |glyphs: &str| {
        TONES
            .iter()
            .map(|&tone| format!(r#"<span style="color:{}">{glyphs}</span>"#, hex(tone)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    format!(
        r#"<span class="mark" aria-label="faba"><i>{}</i><i>{} <b style="color:{}">faba</b></i></span>"#,
        row(BEAN[0]),
        row(BEAN[1]),
        hex(TONES[0])
    )
}

/// A writer passing on JSON with every `<` as `\u003c`, the same string to
/// a JSON reader, so nothing in the data can end or reopen the script
/// element holding it (`</script>`, `<!--`), however the writes split it.
/// JSON has `<` only inside strings, so this is all it need look at.
struct Escaped<W>(W);

impl<W: Write> Write for Escaped<W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        for part in buf.split_inclusive(|&b| b == b'<') {
            match part.split_last() {
                Some((b'<', before)) => {
                    self.0.write_all(before)?;
                    self.0.write_all(b"\\u003c")?;
                }
                _ => self.0.write_all(part)?,
            }
        }
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.0.flush()
    }
}

/// Write the page titled `title`, with `fragment` (a view's markup, styles
/// and code, `{{MARK}}` where the mark goes) and `data` embedded for it to
/// read from `#data`, to `out` as it goes.
fn stream(
    mut out: impl Write,
    title: &str,
    fragment: &str,
    data: &impl Serialize,
) -> anyhow::Result<()> {
    // Each slot filled once in the shell itself, so nothing filled in
    // (a gene named `{{BODY}}`, say) is taken for a slot.
    let (head, rest) = SHELL
        .split_once("{{DATA}}")
        .expect("the shell has a data slot");
    let (middle, tail) = rest
        .split_once("{{BODY}}")
        .expect("the shell has a body slot");
    out.write_all(
        head.replace("{{TITLE}}", &crate::figure::esc(title))
            .as_bytes(),
    )?;
    serde_json::to_writer(Escaped(&mut out), data)?;
    for part in [middle, &fragment.replace("{{MARK}}", &mark()), tail] {
        out.write_all(part.as_bytes())?;
    }
    out.flush()?;
    Ok(())
}

/// Write the page (see [`stream`]) to `path`, put in place only once whole
/// (see [`crate::common::write_replacing`]).
pub fn write(
    path: &Path,
    title: &str,
    fragment: &str,
    data: &impl Serialize,
) -> anyhow::Result<()> {
    crate::common::write_replacing(path, |w| stream(w, title, fragment, data))
        .map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))
}

/// To a tenth, which is all a page shows of a read count.
pub(crate) fn tenth(x: f64) -> f64 {
    (10.0 * x).round() / 10.0
}

/// To a thousandth, enough for two decimals or a fraction of a site.
pub(crate) fn thousandth(x: f64) -> f64 {
    (1000.0 * x).round() / 1000.0
}

/// A report written on a thread of its own, one at a time: asked again
/// while one is under way, one more is written once it is done.
#[derive(Default)]
pub struct Writing {
    /// The file being written, and the thread writing it, until polled.
    under_way: Option<(PathBuf, JoinHandle<anyhow::Result<()>>)>,
    /// Asked again meanwhile.
    again: bool,
}

/// A write [`Writing::poll`] found done.
pub struct Written {
    /// Where it went, or why not.
    pub said: String,
    /// It was asked for again meanwhile: write it anew now.
    pub again: bool,
}

impl Written {
    /// What to say: where it went, or why not, and when it was asked for
    /// again, what `restart` (writing it anew) says after it.
    pub fn then(self, restart: impl FnOnce() -> String) -> String {
        if self.again {
            format!("{}; {}", self.said, restart())
        } else {
            self.said
        }
    }
}

impl Writing {
    /// While a write is under way, or done but not yet polled, queue one
    /// more and say so; none when [`Self::spawn`] may start one.
    pub fn busy(&mut self) -> Option<String> {
        let (path, _) = self.under_way.as_ref()?;
        self.again = true;
        Some(format!("writing {}, then again", crate::tui::tilde(path)))
    }

    /// Run `work`, which writes `path`, on a thread of its own, once
    /// [`Self::busy`] says none is under way: what to say.
    pub fn spawn<W>(&mut self, path: PathBuf, work: W) -> String
    where
        W: FnOnce() -> anyhow::Result<()> + Send + 'static,
    {
        debug_assert!(self.under_way.is_none(), "one write at a time");
        let said = format!("writing {}", crate::tui::tilde(&path));
        self.under_way = Some((path, std::thread::spawn(work)));
        said
    }

    /// The write once it is done, however long ago, taken so it is said
    /// once; none while it is under way or when there is none.
    pub fn poll(&mut self) -> Option<Written> {
        if !self.under_way.as_ref()?.1.is_finished() {
            return None;
        }
        let said = match self.join()? {
            Ok(msg) | Err(msg) => msg,
        };
        let again = std::mem::take(&mut self.again);
        Some(Written { said, again })
    }

    /// The file of the write still running, if any.
    pub fn running(&self) -> Option<&Path> {
        let (path, thread) = self.under_way.as_ref()?;
        (!thread.is_finished()).then_some(path.as_path())
    }

    /// Wait for the write, if any: where it went, or why not.
    fn join(&mut self) -> Option<Result<String, String>> {
        let (path, thread) = self.under_way.take()?;
        let at = crate::tui::tilde(&path);
        Some(match thread.join() {
            Ok(Ok(())) => Ok(format!("wrote {at}")),
            Ok(Err(e)) => Err(format!("{e:#}")),
            Err(_) => Err(format!("{at}: writing it failed")),
        })
    }
}

impl Drop for Writing {
    /// A write under way is finished, not cut short with the view; the
    /// wait and a failure are said, not hidden.
    fn drop(&mut self) {
        if let Some(path) = self.running() {
            log::warn!("finishing {}", crate::tui::tilde(path));
        }
        match self.join() {
            Some(Ok(msg)) => log::info!("{msg}"),
            Some(Err(msg)) => log::warn!("{msg}"),
            None => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The page as [`stream`] writes it.
    fn streamed(title: &str, fragment: &str, data: &impl Serialize) -> String {
        let mut out = Vec::new();
        stream(&mut out, title, fragment, data).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn a_page_embeds_its_data_and_fragment_safely() {
        let data = serde_json::json!({"name": "</script><b>x</b>", "slot": "{{BODY}}"});
        let (title, fragment) = ("a <b> & c {{DATA}}", "<div>{{MARK}}</div>");
        let html = streamed(title, fragment, &data);
        assert!(html.starts_with("<!doctype html>"));
        // The data read back from its script element as it went in.
        let at = html
            .find(r#"<script id="data" type="application/json">"#)
            .expect("a data element");
        let json = &html[at..];
        let json = &json[json.find('>').unwrap() + 1..json.find("</script>").unwrap()];
        let back: serde_json::Value = serde_json::from_str(json).unwrap();
        assert_eq!(back, data);
        assert!(html.contains("<title>a &lt;b&gt; &amp; c {{DATA}}</title>"));
        assert!(!html.contains("</script><b>"), "data cannot end its script");
        assert!(html.contains(r#"\u003c/script>\u003cb>x\u003c/b>"#));
        assert!(
            html.contains(r#""slot":"{{BODY}}""#),
            "data is left as it is"
        );
        assert!(
            html.contains("faba</b></i></span></div>"),
            "the mark in place"
        );
        // Every slot filled once: the only braces left are the title's and the data's own.
        assert_eq!(html.matches("{{").count(), 2);
        // Nothing from the network.
        assert!(!html.contains("http://") && !html.contains("https://"));
    }

    #[test]
    fn no_markup_in_the_data_can_end_its_script() {
        let text = "a</script><!--<script>x</b>";
        let data = serde_json::json!({ "name": text });
        let html = streamed("t", "{{MARK}}", &data);
        assert!(
            !html.contains("<!--<script>") && !html.contains("x</b>"),
            "{html}"
        );
        // The same string to a JSON reader, every split of the writes alike.
        let json = serde_json::to_string(&data).unwrap();
        for at in 0..=json.len() {
            let mut out = Escaped(Vec::new());
            out.write_all(&json.as_bytes()[..at]).unwrap();
            out.write_all(&json.as_bytes()[at..]).unwrap();
            let back: serde_json::Value = serde_json::from_slice(&out.0).unwrap();
            assert_eq!(back, data, "at {at}");
        }
    }

    /// Wait for `w`'s write to finish, and what it said.
    fn finished(w: &mut Writing) -> Written {
        let mut done = None;
        crate::test_util::wait_for(|| {
            done = w.poll();
            done.is_some()
        });
        done.expect("waited for")
    }

    /// Ask `w` for a write of `path` by `work`, as a view does: queued
    /// while one is under way. What it says.
    fn start<W>(w: &mut Writing, path: &Path, work: W) -> String
    where
        W: FnOnce() -> anyhow::Result<()> + Send + 'static,
    {
        w.busy()
            .unwrap_or_else(|| w.spawn(path.to_path_buf(), work))
    }

    #[test]
    fn a_write_says_where_it_went_or_why_not() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.html");
        let mut w = Writing::default();
        assert!(w.poll().is_none(), "nothing under way");
        let to = path.clone();
        let said = start(&mut w, &path, move || Ok(std::fs::write(&to, "x")?));
        assert_eq!(said, format!("writing {}", crate::tui::tilde(&path)));
        let done = finished(&mut w);
        assert_eq!(done.said, format!("wrote {}", crate::tui::tilde(&path)));
        assert!(!done.again);
        assert!(w.poll().is_none(), "said once");
        start(&mut w, &path, || anyhow::bail!("no room"));
        assert_eq!(finished(&mut w).said, "no room");
        start(&mut w, &path, || panic!("a bug"));
        assert_eq!(
            finished(&mut w).said,
            format!("{}: writing it failed", crate::tui::tilde(&path))
        );
    }

    #[test]
    fn a_press_while_writing_queues_one_more_and_a_late_poll_still_says() {
        let path = PathBuf::from("r.html");
        let (go, wait) = std::sync::mpsc::channel::<()>();
        let mut w = Writing::default();
        start(&mut w, &path, move || {
            wait.recv().ok();
            Ok(())
        });
        assert_eq!(w.running(), Some(path.as_path()));
        assert_eq!(w.busy().as_deref(), Some("writing r.html, then again"));
        go.send(()).unwrap();
        // Polled long after it finished, it still says so, and once.
        crate::test_util::wait_for(|| w.running().is_none());
        assert_eq!(
            w.busy().as_deref(),
            Some("writing r.html, then again"),
            "done but not polled is still busy"
        );
        let done = finished(&mut w);
        assert_eq!(done.said, "wrote r.html");
        assert!(done.again, "the queued one");
        start(&mut w, &path, || Ok(()));
        assert!(!finished(&mut w).again, "queued once only");
    }

    #[test]
    fn a_rerun_says_the_first_writes_outcome_before_its_start() {
        let done = |again| Written {
            said: "wrote r.html".into(),
            again,
        };
        assert_eq!(
            done(true).then(|| "writing r2.html".into()),
            "wrote r.html; writing r2.html"
        );
        assert_eq!(
            done(false).then(|| unreachable!("not asked again")),
            "wrote r.html"
        );
    }

    #[test]
    fn dropping_a_write_under_way_finishes_it() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.html");
        let mut w = Writing::default();
        let to = path.clone();
        start(&mut w, &path, move || {
            std::thread::sleep(std::time::Duration::from_millis(50));
            Ok(std::fs::write(&to, "x")?)
        });
        drop(w);
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "x");
    }

    #[test]
    fn a_failed_write_keeps_the_page_that_was_there() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("r.html");
        write(&path, "t", "{{MARK}}", &serde_json::json!({ "a": 1 })).unwrap();
        let before = std::fs::read_to_string(&path).unwrap();
        // A map keyed by pairs cannot be JSON: the write fails part way.
        let bad: std::collections::BTreeMap<(u8, u8), u8> = [((1, 2), 3)].into();
        assert!(write(&path, "t", "{{MARK}}", &bad).is_err());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), before);
        let left: Vec<_> = std::fs::read_dir(dir.path()).unwrap().collect();
        assert_eq!(left.len(), 1, "nothing left beside it");
    }
}
