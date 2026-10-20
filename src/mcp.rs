//! `faba mcp`: faba's commands as Model Context Protocol tools, served over
//! stdio with `legume-mcp`, so an agent can run a profiling step and read
//! what it wrote.
//!
//! Every tool runs as a child `faba` process, as from the shell: stdout
//! carries the protocol, a child's log lines come back as progress, and a
//! cancelled call kills its child. A command that opens a view is run with
//! `--batch-process`, so the server never waits on a terminal.
//!
//! Writes take an `out` that must not exist yet: no tool replaces a
//! directory or a file. The project's records are legume-mcp's, laid out
//! after Noble (2009), "A quick guide to organizing computational biology
//! projects", PLoS Comput Biol 5(7):e1000424: each writing step goes into
//! the day's replayable `results/<date>/runall.sh` with its outputs
//! (skipped on replay when they exist) and a line in `results/INDEX.tsv`;
//! an output the agent does not name goes in that day's folder. `source`
//! notes an input's origin in `data/README.md`.

use anyhow::{bail, ensure, Context, Result};
use clap::Args;
use legume_mcp::{Call, Content, Handler, ServerInfo, ToolSpec};
use serde_json::{json, Map, Value};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Args, Debug)]
pub struct McpArgs {}

const INSTRUCTIONS: &str = "\
faba turns coordinate-sorted BAM alignments into per-cell feature matrices: \
gene counts, read depth, A-to-I editing, DART-seq m6A, alternative \
polyadenylation and SNP genotypes. Row names are \
{gene_key}/{modality}/{detail}.

Producers write every candidate and apply no p-value cutoff. qc is where \
cells, features and sites are cut, and it writes a new directory.

1. docs lists the method write-ups. With `query` it prints the sections \
named (a topic such as profiling, a command such as qc, or a section number).
2. count, atoi, dartseq, apa, depth and snp each write one modality. run \
writes the pipeline (SNP, then count, A-to-I, m6A, APA) into one directory. \
dartseq needs control BAMs; without them, run skips m6A.
3. qc filters a directory into a new one, with the site thresholds given. \
It does not edit the input.
4. pileup, metagene and pwm read a directory's sites. pileup returns a \
binned table for a gene or region (and, with `figure`, a Miami SVG). \
modules prints the gene and cell modules of a directory that already holds \
cluster's block tables.
5. cluster fits cells by their own site reads (link, blocks or crp) and \
writes its tables to `out`.

Every write goes to `out`, which must not exist. Leave `out` out and the \
output goes in today's results folder, named after the input and the tool. \
Paths resolve against the server's working directory. Further flags go in \
`options`, as `faba <command> --help` lists them; `options` cannot name \
the output or open a view. Results give the paths written.";

fn schema(properties: Value, required: &[&str]) -> Value {
    json!({ "type": "object", "properties": properties, "required": required })
}

fn tools() -> Vec<ToolSpec> {
    let bam = || {
        json!({
            "type": "array", "items": { "type": "string" }, "minItems": 1,
            "description": "Coordinate-sorted BAM files, each with a .bai index"
        })
    };
    let gff = || json!({ "type": "string", "description": "Gene annotation (GFF/GTF)" });
    let genome = || {
        json!({
            "type": "string",
            "description": "Reference FASTA (.fa or .fasta); indexed on use if it has no .fai"
        })
    };
    let out = || {
        json!({
            "type": "string",
            "description": "Output directory; must not exist. Omitted, a new directory in today's results folder"
        })
    };
    let file_out = || {
        json!({
            "type": "string",
            "description": "Output file; must not exist. Omitted, a new file in today's results folder"
        })
    };
    let options = || {
        json!({
            "type": "array", "items": { "type": "string" },
            "description": "Further flags, as `faba <command> --help` lists them. A short option stands apart from its value (`-g X`)"
        })
    };
    let source = || {
        json!({
            "type": "string",
            "description": "Where the inputs came from (a URL or accession, and its version), kept in data/README.md"
        })
    };
    let modality = || json!({ "type": "string", "enum": ["m6a", "atoi"] });
    let dir = || {
        json!({
            "type": "string",
            "description": "A faba output directory (from a producer, `run` or `qc`)"
        })
    };
    vec![
        ToolSpec::new(
            "docs",
            "The method write-ups compiled into faba. Without `query`, the list of topics and sections; with it, those sections (a topic, a command named in a heading, or a section number).",
            schema(
                json!({ "query": { "type": "array", "items": { "type": "string" } } }),
                &[],
            ),
        ),
        ToolSpec::new(
            "count",
            "Count reads per gene (`faba count`): one sparse matrix per BAM, rows `{gene}/count/{spliced|unspliced}`. Writes a new directory.",
            schema(
                json!({ "bam": bam(), "gff": gff(), "out": out(), "source": source(), "options": options() }),
                &["bam", "gff"],
            ),
        ),
        ToolSpec::new(
            "depth",
            "Bin the genome and count read coverage per cell (`faba depth`). `resolution_kb` is the bin size in kilobases. Writes a new directory.",
            schema(
                json!({
                    "bam": bam(),
                    "resolution_kb": { "type": "number", "exclusiveMinimum": 0 },
                    "valid_cells": { "type": "string", "description": "A `faba count` directory whose {batch}_cells.tsv.gz names the cells to keep" },
                    "out": out(), "source": source(), "options": options(),
                }),
                &["bam", "resolution_kb"],
            ),
        ),
        ToolSpec::new(
            "atoi",
            "Call A-to-I editing sites (`faba atoi`). Every putative site is written with its p-value; `qc` is where a cutoff is applied. Writes a new directory.",
            schema(
                json!({ "bam": bam(), "gff": gff(), "genome": genome(), "out": out(), "source": source(), "options": options() }),
                &["bam", "gff", "genome"],
            ),
        ),
        ToolSpec::new(
            "dartseq",
            "Call DART-seq m6A sites by a signal-versus-control contrast on C-to-T (`faba dartseq`). `control` is the catalytically dead arm; without it a genomic C/T cannot be told from an m6A. Writes a new directory.",
            schema(
                json!({
                    "bam": bam(), "control": bam(), "gff": gff(), "genome": genome(),
                    "out": out(), "source": source(), "options": options(),
                }),
                &["bam", "control", "gff", "genome"],
            ),
        ),
        ToolSpec::new(
            "apa",
            "Quantify alternative polyadenylation (`faba apa`). `method` is mixture (SCAPE, the default) or simple. Needs `gff`, or `--utr-bed` in `options`. Writes a new directory.",
            schema(
                json!({
                    "bam": bam(), "gff": gff(),
                    "method": { "type": "string", "enum": ["mixture", "simple"] },
                    "out": out(), "source": source(), "options": options(),
                }),
                &["bam"],
            ),
        ),
        ToolSpec::new(
            "snp",
            "Discover and genotype SNPs (`faba snp`). With `gff`, also writes a per-cell allele-frequency matrix; `bulk` writes the call set only. Writes a new directory.",
            schema(
                json!({
                    "bam": bam(), "genome": genome(), "gff": gff(),
                    "known_snps": { "type": "string", "description": "Known sites (VCF/BCF/Parquet) to force-call" },
                    "bulk": { "type": "boolean" },
                    "skip_discovery": { "type": "boolean", "description": "Genotype `known_snps` only" },
                    "out": out(), "source": source(), "options": options(),
                }),
                &["bam", "genome"],
            ),
        ),
        ToolSpec::new(
            "run",
            "The whole pipeline (`faba run --batch-process`): SNP, gene counts, A-to-I, m6A, APA, into one new directory. m6A runs only when `control` is given. No view.",
            schema(
                json!({
                    "bam": bam(), "gff": gff(), "genome": genome(), "control": bam(),
                    "known_snps": { "type": "string" },
                    "out": out(), "source": source(), "options": options(),
                }),
                &["bam", "gff", "genome"],
            ),
        ),
        ToolSpec::new(
            "qc",
            "Filter a faba directory into a new one (`faba qc --batch-process`). Cells are cut on `{batch}_count` and the same cells are kept in every matrix; sites are cut on the thresholds given (marginal p-values, no FDR). Omit a threshold to keep faba's default.",
            schema(qc_properties(&dir(), &out(), &options()), &["input"]),
        ),
        ToolSpec::new(
            "cluster",
            "Cluster cells by their site reads (`faba cluster`). `data` is a faba directory or site matrices. `model` is link (the default), blocks or crp. Tables go to `out`, a new directory, not back into `data`.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "minItems": 1 },
                    "modality": modality(),
                    "model": { "type": "string", "enum": ["link", "blocks", "crp"] },
                    "cell_roles": { "type": "integer", "minimum": 1 },
                    "gene_roles": { "type": "integer", "minimum": 1 },
                    "chains": { "type": "integer", "minimum": 1 },
                    "sweeps": { "type": "integer", "minimum": 1 },
                    "burn_in": { "type": "integer", "minimum": 0 },
                    "fit": { "type": "string", "description": "Fit only matrices whose path matches this pattern (`*`, `?`); the rest are assigned after" },
                    "out": out(), "options": options(),
                }),
                &["data"],
            ),
        ),
        ToolSpec::new(
            "modules",
            "Print the gene and cell modules `cluster` wrote (`faba modules --batch-process`). `data` must already hold the block tables; this does not fit. With `html`, write one self-contained page to a new file instead.",
            schema(
                json!({
                    "data": dir(),
                    "modality": modality(),
                    "gtf": { "type": "string", "description": "Gene annotation for the models drawn in the HTML page" },
                    "html": { "type": "string", "description": "Write the page here; must not exist" },
                    "options": options(),
                }),
                &["data"],
            ),
        ),
        ToolSpec::new(
            "pileup",
            "A gene's or region's sites (`faba pileup --batch-process`). Returns a binned TSV. The ASCII histogram is the call's progress. With `figure` (and `gtf` or `bam`), write a Miami SVG instead; `out` is then the figure prefix (`{out}.miami.svg`).",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "minItems": 1, "description": "Site matrices, or a faba output directory" },
                    "genes": { "type": "array", "items": { "type": "string" }, "description": "Gene symbols or Ensembl ids" },
                    "regions": { "type": "array", "items": { "type": "string" }, "description": "Regions as chr:lb-ub" },
                    "modality": modality(),
                    "figure": { "type": "boolean", "description": "Write the Miami SVG instead of the TSV" },
                    "gtf": { "type": "string", "description": "Gene annotation for the figure's gene model" },
                    "bam": { "type": "array", "items": { "type": "string" }, "description": "BAMs for the figure's read-depth track" },
                    "out": file_out(),
                    "options": options(),
                }),
                &["data"],
            ),
        ),
        ToolSpec::new(
            "metagene",
            "Where sites sit along 5'UTR, CDS and 3'UTR (`faba metagene --batch-process`), as a TSV. Bin heights are raw counts. The ASCII histogram is the call's progress.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "minItems": 1, "description": "Site parquets, or faba output directories" },
                    "gff": gff(),
                    "modality": modality(),
                    "bins": { "type": "integer", "minimum": 1, "description": "Bins across the metagene (default 200)" },
                    "dist_measures": { "type": "string", "description": "Also write MetaPlotR's per-site distance table here; must not exist" },
                    "out": file_out(), "options": options(),
                }),
                &["data"],
            ),
        ),
        ToolSpec::new(
            "pwm",
            "Base frequencies in a window around the sites (`faba pwm --batch-process`), as a TSV. `source` here is reads or reference, not where the file came from.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "minItems": 1, "description": "Site parquets or faba output directories; BAMs too when source is reads" },
                    "genome": genome(),
                    "modality": modality(),
                    "source": { "type": "string", "enum": ["reference", "reads"] },
                    "window": { "type": "integer", "minimum": 0, "description": "Half-window in bp (default 10)" },
                    "out": file_out(), "options": options(),
                }),
                &["data"],
            ),
        ),
    ]
}

/// `qc`'s arguments: the directory in, the directory out, and the site cuts.
fn qc_properties(dir: &Value, out: &Value, options: &Value) -> Value {
    json!({
        "input": dir,
        "out": out,
        "site_max_pv": { "type": "number", "description": "Keep a site with p-value at most this (marginal, no FDR)" },
        "site_min_log_odds": { "type": "number", "description": "m6A: keep a site whose signal-versus-control log odds ratio is at least this" },
        "site_min_coverage": { "type": "integer", "minimum": 0 },
        "site_min_converted": { "type": "integer", "minimum": 0 },
        "site_min_edit_ratio": { "type": "number" },
        "site_max_edit_ratio": { "type": "number" },
        "site_min_fold": { "type": "number", "description": "m6A: keep a site whose signal edit rate is at least this times the control's" },
        "site_min_cells": { "type": "integer", "minimum": 0, "description": "Keep a site seen as converted in at least this many kept cells" },
        "row_nnz_cutoff": { "type": "integer", "minimum": 0 },
        "column_nnz_cutoff": { "type": "integer", "minimum": 0 },
        "qc_mads": { "type": "number", "exclusiveMinimum": 0, "description": "MAD band for the cell outlier drops" },
        "auto_cutoff": { "type": "boolean" },
        "no_cell_qc": { "type": "boolean" },
        "options": options,
    })
}

pub fn run_mcp(_args: &McpArgs) -> Result<()> {
    let info = ServerInfo {
        name: "faba".into(),
        version: env!("CARGO_PKG_VERSION").into(),
        instructions: Some(INSTRUCTIONS.into()),
    };
    legume_mcp::serve(info, Faba)
}

struct Faba;

impl Handler for Faba {
    fn tools(&self) -> Vec<ToolSpec> {
        tools()
    }

    fn call(&mut self, name: &str, args: Value, cx: &Call) -> Result<Vec<Content>> {
        let a = Arguments(args);
        let v = match name {
            "docs" => return docs(&a, cx),
            "count" => count(&a, cx)?,
            "depth" => depth(&a, cx)?,
            "atoi" => atoi(&a, cx)?,
            "dartseq" => dartseq(&a, cx)?,
            "apa" => apa(&a, cx)?,
            "snp" => snp(&a, cx)?,
            "run" => run(&a, cx)?,
            "qc" => qc(&a, cx)?,
            "cluster" => cluster(&a, cx)?,
            "modules" => return modules(&a, cx),
            "pileup" => return pileup(&a, cx),
            "metagene" => metagene(&a, cx)?,
            "pwm" => pwm(&a, cx)?,
            _ => bail!("no tool {name}"),
        };
        Ok(vec![Content::Json(v)])
    }
}

/// A tool's arguments.
struct Arguments(Value);

impl Arguments {
    fn str(&self, key: &str) -> Option<&str> {
        self.0.get(key).and_then(Value::as_str)
    }

    fn req(&self, key: &str) -> Result<String> {
        self.str(key)
            .map(str::to_string)
            .with_context(|| format!("`{key}` is required"))
    }

    /// A path, which goes on the command line where a leading `-` would
    /// read as an option.
    fn path(&self, key: &str) -> Result<String> {
        let p = self.req(key)?;
        ensure!(!p.starts_with('-'), "`{key}` {p}: write it as ./{p}");
        Ok(p)
    }

    fn path_opt(&self, key: &str) -> Result<Option<String>> {
        match self.str(key) {
            None => Ok(None),
            Some(_) => self.path(key).map(Some),
        }
    }

    fn parse<T: serde::de::DeserializeOwned>(&self, key: &str) -> Result<Option<T>> {
        match self.0.get(key) {
            None | Some(Value::Null) => Ok(None),
            Some(v) => serde_json::from_value(v.clone())
                .map(Some)
                .with_context(|| format!("reading `{key}`")),
        }
    }

    fn flag(&self, key: &str) -> Result<bool> {
        Ok(self.parse(key)?.unwrap_or(false))
    }

    fn strings(&self, key: &str) -> Result<Vec<String>> {
        Ok(self.parse(key)?.unwrap_or_default())
    }

    fn paths(&self, key: &str) -> Result<Vec<String>> {
        let v = self.strings(key)?;
        if let Some(p) = v.iter().find(|p| p.starts_with('-')) {
            bail!("`{key}` {p}: write it as ./{p}");
        }
        Ok(v)
    }

    fn paths_req(&self, key: &str) -> Result<Vec<String>> {
        let v = self.paths(key)?;
        ensure!(!v.is_empty(), "`{key}` is required");
        Ok(v)
    }

    /// The numbers given among `keys`, as `--key-name value`. A value that
    /// starts like a flag is `--key-name=-1`, so clap keeps the number.
    fn numbers(&self, keys: &[&str]) -> Result<Vec<String>> {
        let mut v = vec![];
        for k in keys {
            match self.0.get(*k) {
                None | Some(Value::Null) => {}
                Some(n @ Value::Number(_)) => {
                    v.extend(crate::script::emit(
                        &format!("--{}", k.replace('_', "-")),
                        crate::script::Many::One,
                        vec![n.to_string()],
                    ));
                }
                Some(other) => bail!("`{k}` is a number, not {other}"),
            }
        }
        Ok(v)
    }

    /// `out`, or, when the agent names none, `{name}{ext}` in today's
    /// results folder, `-2`, `-3`… added to `name` (before `ext`) while
    /// `taken` says a name is in use.
    fn out_or(
        &self,
        cx: &Call,
        (name, ext): (&str, &str),
        taken: impl Fn(&str) -> bool,
    ) -> Result<String> {
        if self.0.get("out").is_some_and(|v| !v.is_null()) {
            return self.path("out");
        }
        let dir = cx
            .results_dir()
            .context("`out` is required: this server keeps no results folder")?;
        (1..10_000)
            .map(|k| match k {
                1 => dir.join(format!("{name}{ext}")),
                k => dir.join(format!("{name}-{k}{ext}")),
            })
            .map(|p| p.to_string_lossy().into_owned())
            .find(|p| !taken(p))
            .context("no free name in the results folder")
    }

    /// `options`, refusing those that would name the output, open a view,
    /// or hide a short option inside a cluster (`-rX`). `also` are flags
    /// this one tool reserves for its own arguments.
    fn options(&self, also: &[&str]) -> Result<Vec<String>> {
        let v = self.strings("options")?;
        const TAKEN: &[&str] = &[
            "-o",
            "--output",
            "--out",
            "--batch-process",
            "--batch-mode",
            "--session",
            "--no-detach",
            "--html",
            "--dist-measures",
            "--refit",
            "-I",
            "--interactive",
            "--",
        ];
        for o in &v {
            let name = o.split('=').next().unwrap_or(o);
            let cluster = o.starts_with('-')
                && !o.starts_with("--")
                && o.len() > 2
                && o.chars().nth(1).is_some_and(|c| c.is_ascii_alphabetic());
            ensure!(
                !cluster,
                "{o}: give a short option apart from its value (`-g X`)"
            );
            ensure!(
                !TAKEN.contains(&name) && !also.contains(&name),
                "{o} is not taken in `options`: outputs go to `out`, and no view opens here"
            );
        }
        Ok(v)
    }
}

/// `faba -v` with `argv`. `-v` so the child's log lines come back as
/// progress. The binary is this server's own, so a `cargo install` that
/// replaces it while the server runs still reaches the running image.
fn command(argv: &[String]) -> Result<Command> {
    let mut args = Vec::with_capacity(argv.len() + 1);
    args.push("-v".to_string());
    args.extend(argv.iter().cloned());
    legume_mcp::own_command(args)
}

fn absolute(p: &str) -> String {
    std::path::absolute(p)
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| p.to_string())
}

/// A file or directory a tool is about to write: refused if it exists.
fn fresh(path: &str) -> Result<()> {
    ensure!(
        !Path::new(path).exists(),
        "{path} exists; not replacing it: choose another `out`"
    );
    Ok(())
}

/// Where an input's name comes from: the file name, without a suffix, so a
/// chain of unnamed outputs can be read (`sample.count`).
fn stem_of(path: &str) -> String {
    let name = Path::new(path)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "faba".into());
    name.split('.').next().unwrap_or("faba").to_string()
}

fn named(input: &str, tool: &str) -> String {
    format!("{}.{}", stem_of(input), tool)
}

/// `source`, when the agent gave one, against each input.
fn note(cx: &Call, a: &Arguments, paths: &[String]) {
    let Some(source) = a.str("source") else {
        return;
    };
    for p in paths {
        cx.provenance(Path::new(p), source);
    }
}

fn push_flag(argv: &mut Vec<String>, flag: &str, value: &str) {
    argv.extend([flag.to_string(), value.to_string()]);
}

fn push_opt(argv: &mut Vec<String>, flag: &str, value: Option<&str>) {
    if let Some(value) = value {
        push_flag(argv, flag, value);
    }
}

/// What a write returns: the path, the files now in a directory, and
/// anything the command printed.
fn outcome(path: &str, printed: &str) -> Value {
    let mut v = Map::new();
    v.insert("written".into(), json!(absolute(path)));
    let files = entries(path);
    if !files.is_empty() {
        v.insert("files".into(), json!(files));
    }
    let printed = printed.trim();
    if !printed.is_empty() {
        v.insert("stdout".into(), json!(printed));
    }
    Value::Object(v)
}

/// [`outcome`] for a text file, with its first lines, so the table can be
/// read without a second call.
fn text_outcome(path: &str, printed: &str) -> Value {
    let mut v = outcome(path, printed);
    let Ok(text) = std::fs::read_to_string(path) else {
        return v;
    };
    let lines: Vec<&str> = text.lines().collect();
    let head = lines
        .iter()
        .take(40)
        .copied()
        .collect::<Vec<_>>()
        .join("\n");
    let Some(m) = v.as_object_mut() else {
        return v;
    };
    m.insert("lines".into(), json!(lines.len()));
    m.insert("head".into(), json!(head));
    v
}

/// The files a write left at `path`: the file itself, or a directory's
/// immediate children. A `.zarr` tree is one entry, not every chunk.
fn entries(path: &str) -> Vec<String> {
    let p = Path::new(path);
    if p.is_file() {
        return vec![absolute(path)];
    }
    let mut v: Vec<String> = std::fs::read_dir(p)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| absolute(&e.path().to_string_lossy()))
        .collect();
    v.sort();
    v
}

/// Run `argv` so that it writes `out` (made of `what`) with `-o`.
fn run_out(cx: &Call, mut argv: Vec<String>, out: &str, what: &str) -> Result<Value> {
    fresh(out)?;
    argv.extend(["-o".into(), out.into()]);
    let printed = cx.run_making(command(&argv)?, &[PathBuf::from(out)], what)?;
    Ok(outcome(out, &printed))
}

/// [`run_out`] for a text file, returned with its first lines.
fn run_text(cx: &Call, mut argv: Vec<String>, out: &str, what: &str) -> Result<Value> {
    fresh(out)?;
    argv.extend(["-o".into(), out.into()]);
    let printed = cx.run_making(command(&argv)?, &[PathBuf::from(out)], what)?;
    Ok(text_outcome(out, &printed))
}

/// One comma-separated argument. Clap splits that argument on commas, so a
/// value that already contains one would be read as several.
fn comma(key: &str, paths: &[String]) -> Result<String> {
    no_comma(key, paths)?;
    Ok(paths.join(","))
}

fn no_comma(key: &str, paths: &[String]) -> Result<()> {
    if let Some(p) = paths.iter().find(|p| p.contains(',')) {
        bail!("`{key}` {p}: a comma here is read as another item; rename it");
    }
    Ok(())
}

fn docs(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let mut argv = vec!["docs".to_string()];
    argv.extend(a.strings("query")?);
    let text = cx.run_reading(command(&argv)?)?;
    Ok(vec![Content::Text(text)])
}

fn count(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path("gff")?;
    note(cx, a, &bam);
    note(cx, a, std::slice::from_ref(&gff));
    let out = a.out_or(cx, (&named(&bam[0], "count"), ""), exists)?;
    let mut argv = vec!["count".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "-g", &gff);
    argv.extend(a.options(&[])?);
    let what = format!("gene counts of {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn depth(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    ensure!(
        a.0.get("resolution_kb").is_some_and(Value::is_number),
        "`resolution_kb` is required"
    );
    note(cx, a, &bam);
    let out = a.out_or(cx, (&named(&bam[0], "depth"), ""), exists)?;
    let mut argv = vec!["depth".into(), comma("bam", &bam)?];
    argv.extend(a.numbers(&["resolution_kb"])?);
    if let Some(cells) = a.path_opt("valid_cells")? {
        push_flag(&mut argv, "--valid-cells", &cells);
    }
    argv.extend(a.options(&[])?);
    let what = format!("read depth of {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn atoi(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path("gff")?;
    let genome = a.path("genome")?;
    note(cx, a, &bam);
    note(cx, a, &[gff.clone(), genome.clone()]);
    let out = a.out_or(cx, (&named(&bam[0], "atoi"), ""), exists)?;
    let mut argv = vec!["atoi".into(), comma("bam", &bam)?];
    push_flag(&mut argv, "-g", &gff);
    push_flag(&mut argv, "-f", &genome);
    argv.extend(a.options(&[])?);
    let what = format!("A-to-I sites of {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn dartseq(a: &Arguments, cx: &Call) -> Result<Value> {
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
    argv.extend(a.options(&[])?);
    let what = format!("DART-seq m6A sites of {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn apa(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path_opt("gff")?;
    let opts = a.options(&[])?;
    ensure!(
        gff.is_some()
            || opts
                .iter()
                .any(|o| o == "--utr-bed" || o.starts_with("--utr-bed=")),
        "`gff` is required, or pass `--utr-bed` in `options`"
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
    if let Some(method) = a.str("method") {
        ensure!(
            method == "mixture" || method == "simple",
            "`method` is mixture or simple, not {method}"
        );
        push_flag(&mut argv, "--method", method);
    }
    argv.extend(opts);
    let what = format!("alternative polyadenylation of {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn snp(a: &Arguments, cx: &Call) -> Result<Value> {
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
    if a.flag("bulk")? {
        argv.push("--bulk".into());
    }
    if a.flag("skip_discovery")? {
        argv.push("--skip-discovery".into());
    }
    argv.extend(a.options(&[])?);
    let what = format!("SNP genotypes of {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn run(a: &Arguments, cx: &Call) -> Result<Value> {
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
    argv.extend(a.options(&[])?);
    argv.push("--batch-process".into());
    let what = format!("the pipeline on {}", bam.join(", "));
    run_out(cx, argv, &out, &what)
}

fn qc(a: &Arguments, cx: &Call) -> Result<Value> {
    let input = a.path("input")?;
    let out = a.out_or(cx, (&named(&input, "qc"), ""), exists)?;
    let mut argv = vec!["qc".into(), input.clone()];
    argv.extend(a.numbers(&[
        "site_max_pv",
        "site_min_log_odds",
        "site_min_coverage",
        "site_min_converted",
        "site_min_edit_ratio",
        "site_max_edit_ratio",
        "site_min_fold",
        "site_min_cells",
        "row_nnz_cutoff",
        "column_nnz_cutoff",
        "qc_mads",
    ])?);
    if a.flag("auto_cutoff")? {
        argv.push("--auto-cutoff".into());
    }
    if a.flag("no_cell_qc")? {
        argv.push("--no-cell-qc".into());
    }
    argv.extend(a.options(&[])?);
    argv.push("--batch-process".into());
    let what = format!("qc of {input}");
    run_out(cx, argv, &out, &what)
}

fn cluster(a: &Arguments, cx: &Call) -> Result<Value> {
    let data = a.paths_req("data")?;
    let out = a.out_or(cx, (&named(&data[0], "cluster"), ""), exists)?;
    let mut argv = vec!["cluster".into()];
    argv.extend(data.iter().cloned());
    push_opt(&mut argv, "--modality", a.str("modality"));
    if let Some(model) = a.str("model") {
        ensure!(
            matches!(model, "link" | "blocks" | "crp"),
            "`model` is link, blocks or crp, not {model}"
        );
        push_flag(&mut argv, "--model", model);
    }
    push_opt(&mut argv, "--fit", a.str("fit"));
    argv.extend(a.numbers(&["cell_roles", "gene_roles", "chains", "sweeps", "burn_in"])?);
    argv.extend(a.options(&[])?);
    let what = format!("clusters of {}", data.join(", "));
    run_out(cx, argv, &out, &what)
}

/// Flags that would make `modules` fit again, over the tables it was asked
/// to read.
const MODULES_TAKEN: &[&str] = &[
    "--model",
    "--cell-roles",
    "--communities",
    "--gene-roles",
    "--shift",
    "--global",
    "--profile-strength",
    "--prior-strength",
    "--chains",
    "--sweeps",
    "--burn-in",
];

fn modules(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.path("data")?;
    block_tables(&data, a.str("modality"))?;
    let mut argv = vec!["modules".into(), data.clone()];
    push_opt(&mut argv, "--modality", a.str("modality"));
    if let Some(gtf) = a.path_opt("gtf")? {
        push_flag(&mut argv, "--gtf", &gtf);
    }
    argv.extend(a.options(MODULES_TAKEN)?);
    argv.push("--batch-process".into());
    if let Some(html) = a.path_opt("html")? {
        fresh(&html)?;
        push_flag(&mut argv, "--html", &html);
        let what = format!("modules of {data} as a page");
        cx.run_making(command(&argv)?, &[PathBuf::from(&html)], &what)?;
        return Ok(vec![Content::Json(json!({ "written": absolute(&html) }))]);
    }
    let text = cx.run_reading(command(&argv)?)?;
    Ok(vec![Content::Text(text)])
}

/// `dir` holds the block tables `cluster` writes, so `modules` has nothing
/// to fit and nothing to replace.
fn block_tables(dir: &str, modality: Option<&str>) -> Result<()> {
    let has = |m: &str| {
        Path::new(dir)
            .join(format!("{m}_block_sites.parquet"))
            .is_file()
    };
    match modality {
        Some(m) => ensure!(
            has(m),
            "no {m} block tables in {dir}: run cluster into a new directory first"
        ),
        None => ensure!(
            has("m6a") || has("atoi"),
            "no block tables in {dir}: run cluster into a new directory first"
        ),
    }
    Ok(())
}

/// Flags that switch `pileup` into the Miami figure. They are arguments of
/// their own, because the figure writes an SVG and does not write the TSV.
const PILEUP_TAKEN: &[&str] = &["--gtf", "--bam", "--svg", "--png", "--format", "--no-pdf"];

fn pileup(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.paths_req("data")?;
    let genes = a.strings("genes")?;
    let regions = a.strings("regions")?;
    ensure!(
        !genes.is_empty() || !regions.is_empty(),
        "give `genes` or `regions`"
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
    push_opt(&mut argv, "--modality", a.str("modality"));
    argv.extend(a.options(PILEUP_TAKEN)?);
    argv.push("--batch-process".into());
    if figure {
        if let Some(gtf) = &gtf {
            push_flag(&mut argv, "--gtf", gtf);
        }
        if !bam.is_empty() {
            argv.push("--bam".into());
            argv.extend(bam);
        }
        // `--out` is the prefix; the file written is `{prefix}.miami.svg`.
        let prefix = a.out_or(cx, (&named(&data[0], "pileup"), ""), |p| {
            Path::new(&format!("{p}.miami.svg")).exists()
        })?;
        let svg = format!("{prefix}.miami.svg");
        fresh(&svg)?;
        push_flag(&mut argv, "--out", &prefix);
        let what = format!("Miami figure of {}", data.join(", "));
        let printed = cx.run_making(command(&argv)?, &[PathBuf::from(&svg)], &what)?;
        return Ok(vec![Content::Json(outcome(&svg, &printed))]);
    }
    let tsv = a.out_or(cx, (&named(&data[0], "pileup"), ".tsv"), exists)?;
    let what = format!("pileup of {}", data.join(", "));
    Ok(vec![Content::Json(run_text(cx, argv, &tsv, &what)?)])
}

fn metagene(a: &Arguments, cx: &Call) -> Result<Value> {
    let data = a.paths_req("data")?;
    let mut argv = vec!["metagene".into()];
    argv.extend(data.iter().cloned());
    if let Some(gff) = a.path_opt("gff")? {
        push_flag(&mut argv, "-g", &gff);
    }
    push_opt(&mut argv, "--modality", a.str("modality"));
    argv.extend(a.numbers(&["bins"])?);
    argv.extend(a.options(&[])?);
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
    let what = format!("metagene of {}", data.join(", "));
    let printed = cx.run_making(command(&argv)?, &outputs, &what)?;
    let mut v = text_outcome(&tsv, &printed);
    if let (Some(dist), Some(m)) = (a.str("dist_measures"), v.as_object_mut()) {
        m.insert("dist_measures".into(), json!(absolute(dist)));
    }
    Ok(v)
}

fn pwm(a: &Arguments, cx: &Call) -> Result<Value> {
    let data = a.paths_req("data")?;
    no_comma("data", &data)?;
    let mut argv = vec!["pwm".into()];
    argv.extend(data.iter().cloned());
    if let Some(genome) = a.path_opt("genome")? {
        push_flag(&mut argv, "-f", &genome);
    }
    push_opt(&mut argv, "--modality", a.str("modality"));
    if let Some(source) = a.str("source") {
        ensure!(
            source == "reference" || source == "reads",
            "`source` is reference or reads, not {source}"
        );
        push_flag(&mut argv, "--source", source);
    }
    argv.extend(a.numbers(&["window"])?);
    argv.extend(a.options(&[])?);
    argv.push("--batch-process".into());
    let tsv = a.out_or(cx, (&named(&data[0], "pwm"), ".tsv"), exists)?;
    let what = format!("the PWM around sites of {}", data.join(", "));
    run_text(cx, argv, &tsv, &what)
}

fn exists(path: &str) -> bool {
    Path::new(path).exists()
}
