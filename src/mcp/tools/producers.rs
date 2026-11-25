//! The producers: each writes one modality, or `run` the pipeline, into
//! a new directory.

use crate::mcp::args::{comma, push_flag, Arguments};
use crate::mcp::run::{exists, named, note, run_out};
use anyhow::{ensure, Result};
use legume_mcp::Call;
use serde_json::Value;

pub(in crate::mcp) fn count(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path("gff")?;
    note(cx, a, &bam);
    note(cx, a, std::slice::from_ref(&gff));
    let out = a.out_or(cx, (&named(&bam[0], "count"), ""), exists)?;
    let mut argv = vec!["count".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "-g", &gff);
    argv.extend(a.options("count")?);
    let what = format!("gene counts of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn depth(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    ensure!(
        a.0.get("resolution_kb").is_some_and(Value::is_number),
        "`resolution_kb` is required"
    );
    note(cx, a, &bam);
    let out = a.out_or(cx, (&named(&bam[0], "depth"), ""), exists)?;
    let mut argv = vec!["depth".into(), comma("bam", &bam)?];
    if let Some(cells) = a.path_opt("valid_cells")? {
        push_flag(&mut argv, "--valid-cells", &cells);
    }
    argv.extend(a.options("depth")?);
    let what = format!("read depth of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn atoi(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path("gff")?;
    let genome = a.path("genome")?;
    note(cx, a, &bam);
    note(cx, a, &[gff.clone(), genome.clone()]);
    let out = a.out_or(cx, (&named(&bam[0], "atoi"), ""), exists)?;
    let mut argv = vec!["atoi".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "-g", &gff);
    push_flag(&mut argv, "-f", &genome);
    argv.extend(a.options("atoi")?);
    let what = format!("A-to-I sites of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn dartseq(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let control = a.paths_req("control")?;
    let gff = a.path("gff")?;
    let genome = a.path("genome")?;
    note(cx, a, &bam);
    note(cx, a, &control);
    note(cx, a, &[gff.clone(), genome.clone()]);
    let out = a.out_or(cx, (&named(&bam[0], "dartseq"), ""), exists)?;
    let mut argv = vec!["dartseq".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "--control-bam", &comma("control", &control)?);
    push_flag(&mut argv, "-g", &gff);
    push_flag(&mut argv, "-f", &genome);
    argv.extend(a.options("dartseq")?);
    let what = format!("DART-seq m6A sites of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn apa(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path_opt("gff")?;
    let utr_bed = a.path_opt("utr_bed")?;
    let opts = a.options("apa")?;
    ensure!(
        gff.is_some() || utr_bed.is_some(),
        "give `gff` or `utr_bed`"
    );
    note(cx, a, &bam);
    if let Some(gff) = &gff {
        note(cx, a, std::slice::from_ref(gff));
    }
    let out = a.out_or(cx, (&named(&bam[0], "apa"), ""), exists)?;
    let mut argv = vec!["apa".into(), comma("bam", &bam)?];
    if let Some(gff) = &gff {
        push_flag(&mut argv, "-g", gff);
    }
    if let Some(bed) = &utr_bed {
        push_flag(&mut argv, "--utr-bed", bed);
    }
    argv.extend(opts);
    let what = format!("alternative polyadenylation of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn snp(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let genome = a.path("genome")?;
    note(cx, a, &bam);
    note(cx, a, std::slice::from_ref(&genome));
    let out = a.out_or(cx, (&named(&bam[0], "snp"), ""), exists)?;
    let mut argv = vec!["snp".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "-f", &genome);
    if let Some(gff) = a.path_opt("gff")? {
        note(cx, a, std::slice::from_ref(&gff));
        push_flag(&mut argv, "-g", &gff);
    }
    if let Some(known) = a.path_opt("known_snps")? {
        push_flag(&mut argv, "--known-snps", &known);
    }
    argv.extend(a.options("snp")?);
    let what = format!("SNP genotypes of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}

pub(in crate::mcp) fn run(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path("gff")?;
    let genome = a.path("genome")?;
    let control = a.paths("control")?;
    note(cx, a, &bam);
    note(cx, a, &control);
    note(cx, a, &[gff.clone(), genome.clone()]);
    let out = a.out_or(cx, (&named(&bam[0], "run"), ""), exists)?;
    let mut argv = vec!["run".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "-g", &gff);
    push_flag(&mut argv, "-f", &genome);
    if !control.is_empty() {
        push_flag(&mut argv, "--control-bam", &comma("control", &control)?);
    }
    if let Some(known) = a.path_opt("known_snps")? {
        push_flag(&mut argv, "--known-snps", &known);
    }
    argv.extend(a.options("run")?);
    argv.push("--batch-process".into());
    let what = format!("the pipeline on {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
}
