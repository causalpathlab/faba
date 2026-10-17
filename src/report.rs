//! Self-contained HTML reports of the full-screen views, for sharing or
//! posting (a GitHub page, a lab wiki): one file, nothing loaded from
//! anywhere, its data embedded as JSON.
//!
//! Every view's report is a page fragment of its own (its panes, styles and
//! drawing code) inside one shared shell (`report/shell.html`): the theme,
//! the window with the faba mark, the key footer, the legend, the canvas
//! helpers, and the data. A view offers it the same way everywhere: a
//! `--html FILE` flag, and the `e` key in its full-screen view.

use crate::figure::logo::{BEAN, TONES};
use ratatui::style::Color;
use serde::Serialize;
use std::io::Write;
use std::path::Path;

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
