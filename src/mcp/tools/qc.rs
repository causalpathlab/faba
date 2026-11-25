//! `qc`, `qc_view` and `site_plot`: cutting a directory's cells and
//! sites, reading back what a cut kept and dropped, and drawing its site
//! statistics.

use crate::mcp::args::Arguments;
use crate::mcp::run::{
    absolute, command, exists, fresh, named, pngs, qc_images, run_out, tsv_rows,
};
use anyhow::Result;
use legume_mcp::{Call, Content};
use serde_json::json;
use std::path::{Path, PathBuf};

pub(in crate::mcp) fn qc(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let input = a.path("input")?;
    let mut argv = vec!["qc".into(), input.clone()];
    argv.extend(a.options("qc")?);
    // A page that exists is refused by `faba qc` before it reads anything.
    let html = a.path_opt("html")?;
    if let Some(h) = &html {
        argv.push(format!("--html={h}"));
    }
    let plots = a.flag("plots")?;
    if plots {
        argv.push("--plots".into());
    }
    if a.flag("dry_run")? {
        argv.push("--dry-run".into());
        // The figures, without cutting: under `out`, new, never a fileset.
        let out = match plots {
            true => {
                let out = a.out_or(cx, (&named(&input, "qc_plots"), ""), exists)?;
                argv.push(format!("--output={out}"));
                Some(out)
            }
            false => None,
        };
        let printed = cx.run(command(&argv)?)?;
        // The counts, then after a blank line Auto's pairs: two tables.
        let (counts, pairs) = printed.split_once("\n\n").unwrap_or((&printed, ""));
        let mut content = vec![Content::Json(json!({
            "dry_run": tsv_rows(counts),
            "pairs": tsv_rows(pairs),
            "html": html,
            "plots": out.as_ref().map(|o| format!("{o}/qc_plots")),
            "note": "no fileset written, only `plots` and `html` when asked; `removes` counts every site failing a threshold, `only_this` those failing it alone; `pairs` are those Auto picks for the 2D panels, `split` the sites one of their thresholds keeps and the other drops",
        }))];
        if let Some(out) = &out {
            content.extend(qc_images(out)?);
        }
        return Ok(content);
    }
    argv.push("--batch-process".into());
    let out = a.out_or(cx, (&named(&input, "qc"), ""), exists)?;
    let what = format!("qc of {input}");
    let v = run_out(a, cx, argv, &out, &what)?;
    let mut content = vec![Content::Json(v)];
    content.extend(qc_images(&out)?);
    Ok(content)
}

/// `faba site-plot`'s command line for `a`, short of `-o`.
fn site_plot_argv(a: &Arguments) -> Result<Vec<String>> {
    let mut argv = vec!["site-plot".into(), a.path("input")?];
    argv.extend(a.options("site_plot")?);
    Ok(argv)
}

/// Most images attached to one `site_plot` call.
const SITE_PLOT_IMAGES: usize = 12;

pub(in crate::mcp) fn site_plot(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let input = a.path("input")?;
    let argv = site_plot_argv(a)?;
    let out = a.out_or(cx, (&named(&input, "site_plot"), ""), exists)?;
    let what = format!("site panels of {input}");
    let v = run_out(a, cx, argv, &out, &what)?;
    let mut content = vec![Content::Json(v)];
    for p in pngs(Path::new(&out)).iter().take(SITE_PLOT_IMAGES) {
        content.push(Content::Png(std::fs::read(p)?));
    }
    Ok(content)
}

/// `faba qc-view`'s command line for `a`, the figures going under `out`
/// when `plots` asks for them: always `--batch-process`, so the table is
/// printed and no view opens.
fn qc_view_argv(a: &Arguments, out: Option<&str>) -> Result<Vec<String>> {
    let mut argv = vec!["qc-view".into(), a.path("input")?];
    argv.extend(a.options("qc_view")?);
    if let Some(h) = a.path_opt("html")? {
        argv.push(format!("--html={h}"));
    }
    if let Some(out) = out {
        argv.extend(["--plots".into(), format!("--output={out}")]);
    }
    argv.push("--batch-process".into());
    Ok(argv)
}

pub(in crate::mcp) fn qc_view(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let input = a.path("input")?;
    let html = a.path_opt("html")?;
    let out = match a.flag("plots")? {
        true => Some(a.out_or(cx, (&named(&input, "qc_view"), ""), exists)?),
        false => None,
    };
    // The panel `s` saves, as a state arranges the view.
    let figure = a.path_opt("figure")?;
    let figures: Vec<String> = match &figure {
        Some(prefix) => {
            let formats: Vec<String> = a.strings("format")?;
            let formats = match formats.is_empty() {
                true => vec!["pdf".to_string(), "png".to_string()],
                false => formats,
            };
            formats.iter().map(|f| format!("{prefix}.{f}")).collect()
        }
        None => Vec::new(),
    };
    // Neither the page, the figures' folder nor a figure is ever replaced.
    let made: Vec<PathBuf> = (html.iter().chain(&out).chain(&figures))
        .map(PathBuf::from)
        .collect();
    for p in html.iter().chain(&out).chain(&figures) {
        fresh(p)?;
    }
    let mut argv = qc_view_argv(a, out.as_deref())?;
    if let Some(prefix) = &figure {
        argv.insert(2, format!("--figure={prefix}"));
    }
    let cmd = command(&argv)?;
    let printed = match made.is_empty() {
        true => cx.run_reading(cmd)?,
        false => cx.run_making(cmd, &made, &format!("the qc view of {input}"))?,
    };
    // The counts, then after a blank line Auto's pairs: two tables.
    let (counts, pairs) = printed.split_once("\n\n").unwrap_or((&printed, ""));
    let mut content = vec![Content::Json(json!({
        "table": tsv_rows(counts),
        "pairs": tsv_rows(pairs),
        "html": html.as_deref().map(absolute),
        "plots": out.as_ref().map(|o| absolute(&format!("{o}/qc_plots"))),
        "figure": figures.iter().map(|f| absolute(f)).collect::<Vec<_>>(),
        "note": "nothing cut or written but `plots`, `html` and `figure` when asked; the counts are at the thresholds the qc run applied (none for a directory no qc run wrote): `removes` counts every site failing a threshold, `only_this` those failing it alone",
    }))];
    if let Some(png) = figures.iter().find(|f| f.ends_with(".png")) {
        content.push(Content::Png(std::fs::read(png)?));
    }
    if let Some(out) = &out {
        content.extend(qc_images(out)?);
    }
    Ok(content)
}

#[cfg(test)]
mod tests {
    use super::{qc_view_argv, site_plot_argv};
    use crate::mcp::args::{typed_flags, with_typed, Arguments};
    use crate::mcp::schema::tools;
    use serde_json::json;
    /// `site_plot` offers every panel option, the thresholds and the look,
    /// and runs `faba site-plot` with them as flags; `qc` takes the same
    /// panel options.
    #[test]
    fn site_plot_takes_the_panels_and_builds_its_command() {
        let spec = tools()
            .into_iter()
            .find(|t| t.name == "site_plot")
            .expect("site_plot");
        let props = &spec.input_schema["properties"];
        for key in [
            "input",
            "modality",
            "hist",
            "pairs",
            "marginals",
            "style",
            "marginal_size",
            "scale",
            "bins",
            "site_max_pv",
            "site_min_log_odds",
            "site_min_coverage",
            "site_min_converted",
            "site_min_edit_ratio",
            "site_max_edit_ratio",
            "site_min_fold",
            "site_min_cells",
            "width",
            "height",
            "font_size",
            "panel_only",
            "palette",
            "out",
            "options",
        ] {
            assert!(props.get(key).is_some(), "site_plot: no `{key}`");
        }
        assert_eq!(spec.input_schema["required"], json!(["input"]));
        let args = json!({
            "input": "runs/a", "modality": "m6a", "hist": ["pv", "cells"],
            "pairs": ["pv:cells", "fold:pv"], "marginals": ["both"],
            "style": "view", "marginal_size": 0.1, "scale": ["pv=log"], "bins": 30, "site_max_pv": 0.05, "site_min_cells": 10,
            "width": 1.8, "height": 1.4, "font_size": 6, "panel_only": true, "palette": "gray",
        });
        let argv = site_plot_argv(&Arguments(with_typed("site_plot", args).unwrap())).unwrap();
        let has = |w: &[&str]| argv.windows(w.len()).any(|x| x == w);
        assert_eq!(argv[..2], ["site-plot", "runs/a"]);
        for w in [
            &["--modality", "m6a"][..],
            &["--hist", "pv", "--hist", "cells"],
            &["--pairs", "pv:cells", "--pairs", "fold:pv"],
            &["--marginals", "both"],
            &["--style", "view", "--marginal-size", "0.1"],
            &["--scale", "pv=log"],
            &["--bins", "30"],
            &["--site-max-pv", "0.05"],
            &["--site-min-cells", "10"],
            &["--fig-width", "1.8", "--fig-height", "1.4"],
            &["--font-size", "6", "--panel-only", "--palette", "gray"],
        ] {
            assert!(has(w), "{w:?} in {argv:?}");
        }
        assert!(
            !argv.iter().any(|a| a == "-o" || a == "--batch-process"),
            "{argv:?}"
        );
        let qc = typed_flags(
            "qc",
            &Arguments(json!({ "pairs": ["pv:cells"], "bins": 20 })),
        );
        assert_eq!(qc.unwrap(), ["--pairs", "pv:cells", "--bins", "20"]);
        // The style and the marginals' size reach every tool that draws the
        // site panels.
        for tool in ["site_plot", "qc", "qc_view"] {
            let a = json!({ "style": "view", "marginal_size": 0.25 });
            let flags = typed_flags(tool, &Arguments(a)).unwrap();
            let want = ["--style", "view", "--marginal-size", "0.25"];
            assert!(flags.windows(4).any(|w| w == want), "{tool}: {flags:?}");
            // And the arms: which in front, how the p-value reads, the key.
            let a = json!({ "style": "view", "front": "mut", "p_axis": "neglog10", "key": true });
            let flags = typed_flags(tool, &Arguments(a)).unwrap();
            for want in [
                &["--front", "mut"][..],
                &["--p-axis", "neglog10"],
                &["--key"],
            ] {
                let has = flags.windows(want.len()).any(|w| w == want);
                assert!(has, "{tool}: {want:?} in {flags:?}");
            }
        }
    }

    /// `qc_view` offers the page, the figures and every panel option and
    /// look, and runs `faba qc-view --batch-process` with them as flags,
    /// the figures under `out` only when `plots` asks; its view's session
    /// and its output are not offered.
    #[test]
    fn qc_view_takes_the_panels_and_builds_its_command() {
        let spec = tools()
            .into_iter()
            .find(|t| t.name == "qc_view")
            .expect("qc_view");
        let props = &spec.input_schema["properties"];
        for key in [
            "input",
            "html",
            "plots",
            "out",
            "gff",
            "genes",
            "hist",
            "pairs",
            "marginals",
            "scale",
            "bins",
            "width",
            "height",
            "font_size",
            "panel_only",
            "palette",
            "options",
        ] {
            assert!(props.get(key).is_some(), "qc_view: no `{key}`");
        }
        for key in ["output", "batch_process", "session", "no_detach"] {
            assert!(props.get(key).is_none(), "qc_view: `{key}` offered");
        }
        assert_eq!(spec.input_schema["required"], json!(["input"]));
        let args = json!({
            "input": "runs/qc", "html": "page.html", "plots": true, "hist": ["pv"],
            "marginals": ["top"], "width": 1.8, "height": 1.4, "panel_only": true, "gff": "g.gtf",
        });
        let a = Arguments(with_typed("qc_view", args).unwrap());
        let argv = qc_view_argv(&a, Some("figs")).unwrap();
        let has = |w: &[&str]| argv.windows(w.len()).any(|x| x == w);
        assert_eq!(argv[..2], ["qc-view", "runs/qc"]);
        assert_eq!(argv.last().map(String::as_str), Some("--batch-process"));
        for w in [
            &["--hist", "pv"][..],
            &["--marginals", "top"],
            &["--fig-width", "1.8", "--fig-height", "1.4"],
            &["--panel-only"],
            &["--gff", "g.gtf"],
            &["--html=page.html"],
            &["--plots", "--output=figs"],
        ] {
            assert!(has(w), "{w:?} in {argv:?}");
        }
        let bare = qc_view_argv(&Arguments(json!({ "input": "runs/qc" })), None).unwrap();
        assert_eq!(bare, ["qc-view", "runs/qc", "--batch-process"]);
        // What the command line asks for parses as `faba qc-view` takes it.
        let line = std::iter::once("faba".to_string()).chain(argv);
        assert!(crate::matches(&mut crate::faba_command(), line).is_ok());
        let e = Arguments(json!({ "options": ["--batch-process"] })).options("qc_view");
        assert!(e.unwrap_err().to_string().contains("adds it itself"));
    }
}
