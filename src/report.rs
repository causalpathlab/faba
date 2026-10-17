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

/// The page titled `title`, with `fragment` (a view's markup, styles and
/// code, `{{MARK}}` where the mark goes) and `data` embedded for it to read
/// from `#data`.
pub fn page(title: &str, fragment: &str, data: &impl Serialize) -> anyhow::Result<String> {
    // JSON in a script element ends at the first `</`, whatever it says.
    let data = serde_json::to_string(data)?.replace("</", "<\\/");
    // Each slot filled once in the shell itself, so nothing filled in
    // (a gene named `{{BODY}}`, say) is taken for a slot.
    let (head, rest) = SHELL
        .split_once("{{DATA}}")
        .expect("the shell has a data slot");
    let (middle, tail) = rest
        .split_once("{{BODY}}")
        .expect("the shell has a body slot");
    Ok([
        head.replace("{{TITLE}}", &crate::figure::esc(title)),
        data,
        middle.to_string(),
        fragment.replace("{{MARK}}", &mark()),
        tail.to_string(),
    ]
    .concat())
}

/// Write `page` to `path`.
pub fn write(path: &Path, page: &str) -> anyhow::Result<()> {
    std::fs::write(path, page).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_page_embeds_its_data_and_fragment_safely() {
        let data = serde_json::json!({"name": "</script><b>x</b>", "slot": "{{BODY}}"});
        let html = page("a <b> & c {{DATA}}", "<div>{{MARK}}</div>", &data).unwrap();
        assert!(html.starts_with("<!doctype html>"));
        assert!(html.contains("<title>a &lt;b&gt; &amp; c {{DATA}}</title>"));
        assert!(!html.contains("</script><b>"), "data cannot end its script");
        assert!(html.contains(r#"<\/script><b>x<\/b>"#));
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
}
