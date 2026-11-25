//! `pileup`, `metagene` and `pwm`: a directory's sites read back as
//! tables and figures.

use crate::mcp::args::{comma, no_comma, push_flag, Arguments};
use crate::mcp::run::{
    absolute, command, exists, figure_prefix, fresh, named, outcome, run_text, svg_png,
    text_outcome, with_figure,
};
use anyhow::{ensure, Result};
use legume_mcp::{Call, Content};
use serde_json::json;
use std::path::{Path, PathBuf};

/// `faba pileup`'s command line for a call, and whether it draws the
/// figure.
pub(in crate::mcp) fn pileup_argv(a: &Arguments) -> Result<(Vec<String>, bool)> {
    let data = a.paths_req("data")?;
    let genes = a.strings("genes")?;
    let regions = a.strings("regions")?;
    ensure!(
        !genes.is_empty() || !regions.is_empty() || a.str("gene_file").is_some(),
        "give `genes`, `gene_file` or `regions`"
    );
    let bam = a.paths("bam")?;
    let gtf = a.path_opt("gtf")?;
    let figure = a.flag("figure")? || gtf.is_some() || !bam.is_empty();
    let mut argv = vec!["pileup".into()];
    argv.extend(data.iter().cloned());
    if !genes.is_empty() {
        push_flag(&mut argv, "-q", &comma("genes", &genes)?);
    }
    if !regions.is_empty() {
        push_flag(&mut argv, "--regions", &comma("regions", &regions)?);
    }
    argv.extend(a.options("pileup")?);
    argv.push("--batch-process".into());
    if figure {
        // Neither a gene model nor a depth track makes the command draw:
        // `--svg` asks for the figure.
        if gtf.is_none() && bam.is_empty() {
            argv.push("--svg".into());
        }
        if let Some(gtf) = &gtf {
            push_flag(&mut argv, "--gtf", gtf);
        }
        if !bam.is_empty() {
            argv.push("--bam".into());
            argv.extend(bam);
        }
    }
    Ok((argv, figure))
}

pub(in crate::mcp) fn pileup(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.paths_req("data")?;
    let (mut argv, figure) = pileup_argv(a)?;
    if figure && several(a)? {
        return pileup_genes(a, cx, argv, &data);
    }
    if figure {
        // `--out` is the prefix; the file written is `{prefix}.miami.svg`.
        let prefix = a.out_or(cx, (&named(&data[0], "pileup"), ""), |p| {
            Path::new(&format!("{p}.miami.svg")).exists()
        })?;
        let svg = format!("{prefix}.miami.svg");
        fresh(&svg)?;
        push_flag(&mut argv, "--out", &prefix);
        let what = format!("Miami figure of {}", data.join(", "));
        let printed = cx.run_making(command(&argv)?, &[PathBuf::from(&svg)], &what)?;
        let mut content = vec![Content::Json(outcome(&svg, &printed))];
        match svg_png(&svg) {
            Ok(png) => content.push(Content::Png(png)),
            Err(e) => log::warn!("{svg}: not drawn as an image: {e:#}"),
        }
        return Ok(content);
    }
    let tsv = a.out_or(cx, (&named(&data[0], "pileup"), ".tsv"), exists)?;
    let what = format!("pileup of {}", data.join(", "));
    Ok(vec![Content::Json(run_text(cx, argv, &tsv, &what)?)])
}

/// Whether a figure call draws several genes' figures: more than one in
/// `genes`, a `gene_file`, or a `row` or `grid` of them.
fn several(a: &Arguments) -> Result<bool> {
    let many = a.strings("genes")?.len() > 1 || a.str("gene_file").is_some();
    Ok(many || a.flag("row")? || a.str("grid").is_some())
}

/// Most figures returned as images by one call of several genes.
const GENE_IMAGES: usize = 8;

/// Several genes' figures: `{out}_{SYMBOL}.svg` each, `{out}_row.svg` or
/// `{out}_grid.svg` of them all, `{out}_key.svg` for a key alone; listed,
/// and the combined one (else the first few) returned as images.
fn pileup_genes(
    a: &Arguments,
    cx: &Call,
    mut argv: Vec<String>,
    data: &[String],
) -> Result<Vec<Content>> {
    let combined = match (a.flag("row")?, a.str("grid")) {
        (true, _) => Some("row"),
        (false, Some(_)) => Some("grid"),
        _ => None,
    };
    let prefix = a.out_or(cx, (&named(&data[0], "pileup"), ""), |p| {
        written_under(p).next().is_some()
    })?;
    ensure!(
        written_under(&prefix).next().is_none(),
        "figures under {prefix}_ exist; not replacing them: choose another `out`"
    );
    push_flag(&mut argv, "--out", &prefix);
    let outputs: Vec<PathBuf> = combined
        .map(|c| PathBuf::from(format!("{prefix}_{c}.svg")))
        .into_iter()
        .collect();
    let what = format!("pileup figures of {}", data.join(", "));
    let printed = cx.run_making(command(&argv)?, &outputs, &what)?;
    let mut svgs: Vec<String> = written_under(&prefix).collect();
    svgs.sort();
    let mut content = vec![Content::Json(json!({
        "written": svgs.iter().map(|p| absolute(p)).collect::<Vec<_>>(),
        "stdout": printed.trim(),
    }))];
    // The figure of them all first, then each gene's.
    let first = combined.map(|c| format!("{prefix}_{c}.svg"));
    let shown = first
        .iter()
        .chain(svgs.iter().filter(|p| Some(*p) != first.as_ref()));
    for svg in shown.take(GENE_IMAGES) {
        match svg_png(svg) {
            Ok(png) => content.push(Content::Png(png)),
            Err(e) => log::warn!("{svg}: not drawn as an image: {e:#}"),
        }
    }
    Ok(content)
}

/// The SVGs written under `prefix` as several genes' figures are:
/// `{prefix}_*.svg`.
fn written_under(prefix: &str) -> impl Iterator<Item = String> {
    let path = Path::new(prefix);
    let dir = match path.parent() {
        Some(d) if !d.as_os_str().is_empty() => d.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let stem = format!(
        "{}_",
        path.file_name()
            .map(|f| f.to_string_lossy().to_string())
            .unwrap_or_default()
    );
    let files = std::fs::read_dir(&dir).into_iter().flatten().flatten();
    files.filter_map(move |e| {
        let name = e.file_name().to_string_lossy().to_string();
        let ours = name.starts_with(&stem) && name.ends_with(".svg");
        ours.then(|| dir.join(name).to_string_lossy().to_string())
    })
}

pub(in crate::mcp) fn metagene(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.paths_req("data")?;
    let mut argv = vec!["metagene".into()];
    argv.extend(data.iter().cloned());
    if let Some(gff) = a.path_opt("gff")? {
        push_flag(&mut argv, "-g", &gff);
    }
    argv.extend(a.options("metagene")?);
    argv.push("--batch-process".into());
    argv.push("--print".into());
    let tsv = a.out_or(cx, (&named(&data[0], "metagene"), ".tsv"), exists)?;
    fresh(&tsv)?;
    push_flag(&mut argv, "-o", &tsv);
    let mut outputs = vec![PathBuf::from(&tsv)];
    if let Some(dist) = a.path_opt("dist_measures")? {
        fresh(&dist)?;
        push_flag(&mut argv, "--dist-measures", &dist);
        outputs.push(PathBuf::from(&dist));
    }
    let png = figure_arg(a, &tsv, &mut argv, &mut outputs)?;
    let what = format!("metagene of {}", data.join(", "));
    let printed = cx.run_making(command(&argv)?, &outputs, &what)?;
    let mut v = text_outcome(&tsv, &printed);
    if let (Some(dist), Some(m)) = (a.str("dist_measures"), v.as_object_mut()) {
        m.insert("dist_measures".into(), json!(absolute(dist)));
    }
    with_figure(v, png.as_deref())
}

/// `--figure` beside `table` unless the call says `figure: false`; the PNG
/// it writes.
fn figure_arg(
    a: &Arguments,
    table: &str,
    argv: &mut Vec<String>,
    outputs: &mut Vec<PathBuf>,
) -> Result<Option<String>> {
    if !a.parse::<bool>("figure")?.unwrap_or(true) {
        return Ok(None);
    }
    let (prefix, png) = figure_prefix(table, outputs)?;
    push_flag(argv, "--figure", &prefix);
    Ok(Some(png))
}

pub(in crate::mcp) fn pwm(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.paths_req("data")?;
    no_comma("data", &data)?;
    let mut argv = vec!["pwm".into()];
    argv.extend(data.iter().cloned());
    if let Some(genome) = a.path_opt("genome")? {
        push_flag(&mut argv, "-f", &genome);
    }
    argv.extend(a.options("pwm")?);
    argv.push("--batch-process".into());
    let tsv = a.out_or(cx, (&named(&data[0], "pwm"), ".tsv"), exists)?;
    fresh(&tsv)?;
    push_flag(&mut argv, "-o", &tsv);
    let mut outputs = vec![PathBuf::from(&tsv)];
    let png = figure_arg(a, &tsv, &mut argv, &mut outputs)?;
    let what = format!("the PWM around sites of {}", data.join(", "));
    let printed = cx.run_making(command(&argv)?, &outputs, &what)?;
    with_figure(text_outcome(&tsv, &printed), png.as_deref())
}

#[cfg(test)]
mod tests {
    use crate::mcp::args::Arguments;

    use super::pileup_argv;
    use serde_json::json;
    /// `figure` alone draws: the command gets `--svg`, which puts it in
    /// figure mode (a GTF or BAMs would, so it is not added then).
    #[test]
    fn a_pileup_figure_asked_alone_is_drawn() {
        let call = |args| pileup_argv(&Arguments(args)).unwrap();
        let (argv, figure) = call(json!({ "data": ["s.zarr"], "genes": ["G1"], "figure": true }));
        assert!(figure && argv.contains(&"--svg".to_string()), "{argv:?}");
        let (argv, figure) = call(json!({ "data": ["s.zarr"], "genes": ["G1"], "gtf": "g.gtf" }));
        assert!(figure && !argv.contains(&"--svg".to_string()), "{argv:?}");
        let (argv, _) =
            call(json!({ "data": ["s.zarr"], "genes": ["G1"], "figure": true, "bam": ["a.bam"] }));
        assert!(!argv.contains(&"--svg".to_string()), "{argv:?}");
        let (argv, figure) = call(json!({ "data": ["s.zarr"], "genes": ["G1"] }));
        assert!(!figure && !argv.contains(&"--svg".to_string()), "{argv:?}");
    }
}
