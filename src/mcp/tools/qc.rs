//! `qc` and `site_plot`: cutting a directory's cells and sites, and
//! drawing its site statistics.

use crate::mcp::args::{push_opt, Arguments};
use crate::mcp::run::{command, exists, named, pngs, qc_images, run_out, tsv_rows};
use crate::mcp::schema::SITE_THRESHOLDS;
use anyhow::Result;
use legume_mcp::{Call, Content};
use serde_json::json;
use std::path::Path;

pub(in crate::mcp) fn qc(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let input = a.path("input")?;
    let mut argv = vec!["qc".into(), input.clone()];
    argv.extend(a.numbers(&SITE_THRESHOLDS)?);
    argv.extend(a.numbers(&["row_nnz_cutoff", "column_nnz_cutoff", "qc_mads"])?);
    if a.flag("auto_cutoff")? {
        argv.push("--auto-cutoff".into());
    }
    if a.flag("no_cell_qc")? {
        argv.push("--no-cell-qc".into());
    }
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
    push_opt(&mut argv, "--modality", a.str("modality"));
    argv.extend(a.numbers(&SITE_THRESHOLDS)?);
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

#[cfg(test)]
mod tests {
    use super::site_plot_argv;
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
            "pairs": ["pv:cells", "fold:pv"], "marginals": ["top", "right"],
            "scale": ["pv=log"], "bins": 30, "site_max_pv": 0.05, "site_min_cells": 10,
            "width": 1.8, "height": 1.4, "font_size": 6, "panel_only": true, "palette": "gray",
        });
        let argv = site_plot_argv(&Arguments(with_typed("site_plot", args).unwrap())).unwrap();
        let has = |w: &[&str]| argv.windows(w.len()).any(|x| x == w);
        assert_eq!(argv[..2], ["site-plot", "runs/a"]);
        for w in [
            &["--modality", "m6a"][..],
            &["--hist", "pv", "--hist", "cells"],
            &["--pairs", "pv:cells", "--pairs", "fold:pv"],
            &["--marginals", "top", "--marginals", "right"],
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
    }
}
