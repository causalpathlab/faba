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
named (a topic such as profiling, a command such as qc, or a section \
number).
2. count, atoi, dartseq, apa, depth and snp each write one modality. run \
writes the pipeline (SNP, then count, A-to-I, m6A, APA) into one directory. \
dartseq needs control BAMs; without them, run skips m6A.
3. qc filters a directory into a new one, with the site thresholds given. It \
does not edit the input. With `dry_run` it writes nothing and returns what \
the cell QC keeps and what each threshold removes: try thresholds that way \
first.
4. pileup, metagene and pwm read a directory's sites. metagene and pwm \
return their table and its figure as an image. pileup returns a binned table \
for a gene or region (and, with `figure`, a Miami SVG, as an image too). \
cluster_view shows a cluster fit level by level, from its run record \
({m}_blockmodel.run.json) or the directory holding its tables: the modules, \
then a gene module, a gene, the cell modules, the trace and verdict, the \
settings.
5. cluster fits cells by their own site reads (a blockmodel of cell and \
gene modules, by variational inference) and writes its tables to `out`. In two rounds: first `modality` count, cell \
modules by gene counts (it writes count_blockmodel.run.json); then, with \
`within` naming that record, m6a or atoi within each of its cell modules. \
The second round writes {m}_within.run.json, which cluster_view opens: the \
list of its fits, and with `group` the fit within one cell module.

The producers (count, depth, atoi, dartseq, apa, snp, run, qc, cluster) take \
`background`: on large BAMs set it, and the call returns a job id at once.

Every write goes to `out`, which must not exist. Leave `out` out and the \
output goes in today's results folder, named after the input and the tool. \
Paths resolve against the server's working directory. Common flags have \
their own arguments; further flags go in `options`, as `faba <command> \
--help` lists them; `options` cannot name the output or open a view. Results \
give the paths written and what the run recorded (`results`: a fit's ELBO \
verdict, its cluster sizes).";

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
    let fit_modality = || {
        json!({
            "type": "string", "enum": ["m6a", "atoi", "count"],
            "description": "m6a or atoi: site reads. count: gene counts, the first of two rounds"
        })
    };
    let dir = || {
        json!({
            "type": "string",
            "description": "A faba output directory (from a producer, `run` or `qc`)"
        })
    };
    let mut specs = vec![
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
            "Quantify alternative polyadenylation (`faba apa`). `method` is mixture (SCAPE, the default) or simple. Needs `gff` or `utr_bed`. Writes a new directory.",
            schema(
                json!({
                    "bam": bam(), "gff": gff(),
                    "utr_bed": { "type": "string", "description": "3'-UTR intervals as BED, each taken as one block; in place of `gff`" },
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
            "Filter a faba directory into a new one (`faba qc --batch-process`). Cells are cut on `{batch}_count` and the same cells are kept in every matrix; sites are cut on the thresholds given (marginal p-values, no FDR). Omit a threshold to keep faba's default. Try thresholds first with `dry_run`: it returns what the cell QC keeps per batch and what each threshold removes per modality, writing nothing. With `plots`, the view's figures go to `qc_plots/` and come back as images.",
            schema(qc_properties(&dir(), &out(), &options()), &["input"]),
        ),
        ToolSpec::new(
            "cluster",
            "Cluster cells by their site reads (`faba cluster`). `data` is a faba directory or site matrices. The model is a blockmodel of cell modules and gene modules, a methylation shift per pair and each cell's own level, fitted by variational inference from k-means starts (`chains` restarts, the best by ELBO kept; `svi_epochs`, `batch_size`, `epochs`, `tol`). The MCMC fits of faba before 0.19 are gone: `model` link or crp, `sweeps`, `burn_in`, `shift`, `profile_strength` and `holdout` are refused. Tables go to `out`, a new directory, not back into `data`. With `modality` count, cluster cells by their gene counts instead (count matrices), fitted alike: its verdict is whether the best restart's ELBO settled. With `within`, the run record of such a count fit (count_blockmodel.run.json), fit the sites within each of its cell modules, one fit a module, into `out`/within/{m}/{group}/, with `out`/{m}_within.run.json naming them all and {m}_within_cells.parquet giving each cell its nested label (C3.2); `data` may then be left out, and the matrices are those beside the count fit's inputs.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "description": "A faba directory, or site (or, with `modality` count, count) matrices; required unless `within` is given" },
                    "modality": fit_modality(),
                    "cell_roles": { "type": "integer", "minimum": 1 },
                    "gene_roles": { "type": "integer", "minimum": 1 },
                    "chains": { "type": "integer", "minimum": 1 },
                    "svi_epochs": { "type": "integer", "minimum": 0, "description": "Stochastic epochs before coordinate ascent" },
                    "batch_size": { "type": "integer", "minimum": 1, "description": "Cells per minibatch" },
                    "epochs": { "type": "integer", "minimum": 1, "description": "Epochs at most" },
                    "tol": { "type": "number", "minimum": 0, "description": "Converged when an epoch moves the ELBO by less than this per count (per read for a site fit)" },
                    "fit": { "type": "string", "description": "Fit only matrices whose path matches this pattern (`*`, `?`); the rest are assigned after" },
                    "out": out(), "options": options(),
                }),
                &[],
            ),
        ),
        ToolSpec::new(
            "cluster_view",
            "View a `cluster` fit (`faba cluster-view --batch-process`), level by level as the browser does: by default the gene modules against the cell modules (level 1). `module` prints one gene module's genes and sites with their values in each cell module (level 2); `gene` one gene's sites in each cell module over every gene module (levels 3 and 4); `cells` the cell modules' sizes; `go` the GO terms (or `gmt` gene sets) over-represented in each gene module, against the fit's own annotated genes (files default to lupin's cache); `batches` how each cell module spreads over the batches (shared, one condition's, or batch-specific: uneven over a condition's replicates), the conditions read off the batch names less a replicate token or given by `batch_groups`; `trace` each restart's ELBO by epoch (an older MCMC fit's: each chain's log joint), the verdict and its advice; with `figure`, a figure instead (`plot`: the modules grid coloured by `look`, the trace, or the batches; `width` and `height` in inches lay it out for a page), returned as an image; `settings` the settings the fit ran with. Several may be given; rows are tab-separated. Values are % converted for a site fit, UMIs per cell for a count fit. `data` is the fit's run record or the directory holding its tables; this never fits. A within record ({m}_within.run.json, from `cluster` with `within`) prints the count fit's modules, each group's fit (cells, verdict; why one was skipped) and the modules within each fitted group; with `group`, every argument applies to the fit within that cell module, and without it the levels asked for are the count fit's. With `html`, write one self-contained page to a new file instead; a within record without `group` writes a page a fit, `html` for the count fit and {stem}.{group}.html beside it for each group, linked to each other.",
            schema(
                json!({
                    "data": { "type": "string", "description": "A `cluster` fit: its run record `{m}_blockmodel.run.json`, the output directory holding its block tables, or a within record `{m}_within.run.json`" },
                    "modality": fit_modality(),
                    "group": { "type": "string", "description": "Of a within record: the cell module of the count fit whose fit to view, as C3 or 3" },
                    "module": { "type": "string", "description": "A gene module to print, as M3 or 3" },
                    "gene": { "type": "string", "description": "A gene to print, by key or symbol" },
                    "cells": { "type": "boolean", "description": "Print the cell modules and their sizes" },
                    "batches": { "type": "boolean", "description": "Print how each cell module spreads over the batches: shared, one condition's, or batch-specific" },
                    "go": { "type": "boolean", "description": "Print the GO terms over-represented in each gene module (hypergeometric against the fit's annotated genes, BH over every module and term)" },
                    "gaf": { "type": "string", "description": "GO annotations for `go` (default: lupin's cached GO Consortium file for the genes' species)" },
                    "obo": { "type": "string", "description": "The GO graph for `gaf` (default: lupin's cached go-basic.obo)" },
                    "gmt": { "type": "string", "description": "Gene sets to test for `go` instead of GO, a GMT file" },
                    "go_min_weight": { "type": "number", "minimum": 0, "maximum": 1, "description": "A gene counts in its main gene module for `go` at this weight or more (default 0.5)" },
                    "batch_groups": { "type": "string", "description": "Which condition each batch is a replicate of, as rep1_mut=mut,rep2_mut=mut; the rest are read off their names" },
                    "trace": { "type": "boolean", "description": "Print each restart's trace, the verdict and its advice" },
                    "settings": { "type": "boolean", "description": "Print the settings the fit ran with" },
                    "gtf": { "type": "string", "description": "Gene annotation for the models drawn in the HTML page" },
                    "html": { "type": "string", "description": "Write the page here; must not exist" },
                    "figure": { "type": "string", "description": "Write a figure as {figure}.pdf and {figure}.png (neither may exist) and return it as an image" },
                    "plot": { "type": "string", "enum": ["modules", "trace", "batches"], "description": "Which figure: the modules grid (gene modules against cell modules; default), the fit's trace, or the cell modules against the batches" },
                    "look": { "type": "string", "enum": ["shift", "percent", "z"], "description": "What colours the modules grid's blocks: the fitted shift (a count fit's log rate; default), the methylated % (UMIs per cell), or that standardized across the cell modules" },
                    "width": { "type": "number", "exclusiveMinimum": 0, "description": "The figure's width in inches, with `height`: laid out for that size, its text kept readable (e.g. 3.5 for a column, 7 for a page)" },
                    "height": { "type": "number", "exclusiveMinimum": 0, "description": "The figure's height in inches, with `width`" },
                    "options": options(),
                }),
                &["data"],
            ),
        ),
        ToolSpec::new(
            "pileup",
            "A gene's or region's sites (`faba pileup --batch-process`). Returns a binned TSV. The ASCII histogram is the call's progress. With `figure` (or `gtf` or `bam`), write a Miami SVG instead, returned as an image too; `out` is then the figure prefix (`{out}.miami.svg`). A `gtf` draws the gene model (else the annotation recorded with the data, if any), `bam` the read-depth track.",
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
            "Where sites sit along 5'UTR, CDS and 3'UTR (`faba metagene --batch-process`), as a TSV and a figure. Bin heights are raw counts. The ASCII histogram is the call's progress.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "minItems": 1, "description": "Site parquets, or faba output directories" },
                    "gff": gff(),
                    "modality": modality(),
                    "bins": { "type": "integer", "minimum": 1, "description": "Bins across the metagene (default 200)" },
                    "dist_measures": { "type": "string", "description": "Also write MetaPlotR's per-site distance table here; must not exist" },
                    "figure": { "type": "boolean", "description": "Also draw the figure beside the table (`{out stem}.pdf` and `.png`) and return it as an image; default true" },
                    "out": file_out(), "options": options(),
                }),
                &["data"],
            ),
        ),
        ToolSpec::new(
            "pwm",
            "Base frequencies in a window around the sites (`faba pwm --batch-process`), as a TSV and a sequence logo. `source` here is reads or reference, not where the file came from.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "minItems": 1, "description": "Site parquets or faba output directories; BAMs too when source is reads" },
                    "genome": genome(),
                    "modality": modality(),
                    "source": { "type": "string", "enum": ["reference", "reads"] },
                    "window": { "type": "integer", "minimum": 0, "description": "Half-window in bp (default 10)" },
                    "figure": { "type": "boolean", "description": "Also draw the figure beside the table (`{out stem}.pdf` and `.png`) and return it as an image; default true" },
                    "out": file_out(), "options": options(),
                }),
                &["data"],
            ),
        ),
    ];
    for s in specs.iter_mut() {
        for &(_, key, kind, about) in TYPED.iter().filter(|t| t.0.contains(&s.name.as_str())) {
            s.input_schema["properties"][key] = typed_schema(kind, about);
        }
    }
    // A producer on large BAMs outlasts a client's wait on one call: it can
    // run as a job instead.
    for s in specs
        .iter_mut()
        .filter(|s| PRODUCERS.contains(&s.name.as_str()))
    {
        s.input_schema["properties"]["background"] = json!({
            "type": "boolean",
            "description": "Run as a background job: the call returns its id at once; follow it with job_wait or job_status, and use `written` only once it is done"
        });
    }
    specs
}

/// How a typed argument becomes flags: `--{key with dashes}`, then its value.
#[derive(Clone, Copy, Debug)]
enum Kind {
    Text,
    /// A file: no leading `-`, which would read as a flag.
    Path,
    /// The bare flag when true.
    Switch,
    Number,
    /// The flag once per value.
    List,
    Choice(&'static [&'static str]),
}

/// The producers that read cell barcodes from BAM tags.
const TAGGED: &[&str] = &["count", "depth", "atoi", "dartseq", "apa", "snp", "run"];
/// The producers that also count genes and UMIs.
const TAGGED_GENES: &[&str] = &["count", "atoi", "dartseq", "apa", "snp", "run"];

/// Flags common enough to have their own argument: the tools taking each,
/// its key, kind and description. Given ones go ahead of `options`, which
/// still checks them.
const TYPED: &[(&[&str], &str, Kind, &str)] = &[
    (
        TAGGED,
        "cell_barcode_tag",
        Kind::Text,
        "BAM tag holding the cell barcode",
    ),
    (
        TAGGED_GENES,
        "gene_barcode_tag",
        Kind::Text,
        "BAM tag holding the gene",
    ),
    (
        TAGGED_GENES,
        "umi_tag",
        Kind::Text,
        "BAM tag holding the UMI",
    ),
    (
        TAGGED_GENES,
        "gene_type",
        Kind::Text,
        "Keep genes of this biotype only (protein_coding, lncRNA, pseudogene)",
    ),
    (
        &["atoi", "dartseq", "apa", "run"],
        "mixture",
        Kind::Switch,
        "Also fit the per-gene mixture matrices (a slow EM)",
    ),
    (&["run"], "skip_snp", Kind::Switch, "Skip the SNP step"),
    (&["run"], "skip_count", Kind::Switch, "Skip the gene counts"),
    (&["run"], "skip_atoi", Kind::Switch, "Skip A-to-I"),
    (&["run"], "skip_apa", Kind::Switch, "Skip APA"),
    (
        &["run"],
        "depth_resolution_kb",
        Kind::Number,
        "Read-depth bin width in kb",
    ),
    (
        &["cluster"],
        "seed",
        Kind::Number,
        "Random seed; restart i uses seed + i (default 1)",
    ),
    (
        &["cluster"],
        "within",
        Kind::Path,
        "A count fit's run record (count_blockmodel.run.json): fit the sites within each of its cell modules",
    ),
    (
        &["cluster"],
        "within_min_cells",
        Kind::Number,
        "With `within`: a group with fewer cells to fit is skipped (default 50)",
    ),
    (
        &["pileup"],
        "cell_membership",
        Kind::Path,
        "Cell barcode to cell type (TSV/CSV/Parquet): one panel per type",
    ),
    (
        &["pileup"],
        "top_modality",
        Kind::List,
        "Restrict the top track to these modalities, e.g. m6A",
    ),
    (
        &["pileup"],
        "track",
        Kind::List,
        "Labelled tracks as LABEL=PATTERN over the input paths (`*`, `?`)",
    ),
    (
        &["metagene"],
        "isoforms",
        Kind::Choice(&["weighted", "all"]),
        "A site in k isoforms adds 1/k to each (weighted, the default) or 1 (all)",
    ),
    (
        &["metagene"],
        "include_non_coding",
        Kind::Switch,
        "Also profile non-coding genes, on their own track",
    ),
];

/// The JSON schema of a typed argument.
fn typed_schema(kind: Kind, about: &str) -> Value {
    match kind {
        Kind::Text | Kind::Path => json!({ "type": "string", "description": about }),
        Kind::Switch => json!({ "type": "boolean", "description": about }),
        Kind::Number => json!({ "type": "number", "description": about }),
        Kind::List => {
            json!({ "type": "array", "items": { "type": "string" }, "description": about })
        }
        Kind::Choice(of) => json!({ "type": "string", "enum": of, "description": about }),
    }
}

/// `tool`'s typed arguments given in `a`, as flags.
fn typed_flags(tool: &str, a: &Arguments) -> Result<Vec<String>> {
    let mut argv = Vec::new();
    for &(_, key, kind, _) in TYPED.iter().filter(|t| t.0.contains(&tool)) {
        let Some(v) = a.0.get(key).filter(|v| !v.is_null()) else {
            continue;
        };
        let flag = format!("--{}", key.replace('_', "-"));
        match (kind, v) {
            (Kind::Switch, Value::Bool(on)) => {
                if *on {
                    argv.push(flag);
                }
            }
            (Kind::Number, Value::Number(n)) => push_flag(&mut argv, &flag, &n.to_string()),
            (Kind::Text, Value::String(s)) => push_flag(&mut argv, &flag, s),
            (Kind::Path, Value::String(_)) => push_flag(&mut argv, &flag, &a.path(key)?),
            (Kind::Choice(of), Value::String(s)) if of.contains(&s.as_str()) => {
                push_flag(&mut argv, &flag, s)
            }
            (Kind::List, Value::Array(items)) => {
                for x in items {
                    let x = x
                        .as_str()
                        .with_context(|| format!("`{key}` holds strings, not {x}"))?;
                    push_flag(&mut argv, &flag, x);
                }
            }
            (Kind::Choice(of), v) => bail!("`{key}` is one of {}, not {v}", of.join(", ")),
            (kind, v) => bail!(
                "`{key}` is a {}, not {v}",
                typed_schema(kind, "")["type"].as_str().unwrap_or("value")
            ),
        }
    }
    Ok(argv)
}

/// `args` with `tool`'s typed arguments put ahead of its `options`.
fn with_typed(tool: &str, args: Value) -> Result<Value> {
    let a = Arguments(args);
    let typed = typed_flags(tool, &a)?;
    let Arguments(mut args) = a;
    if typed.is_empty() {
        return Ok(args);
    }
    let given: Vec<String> = Arguments(args.clone()).strings("options")?;
    if let Some(m) = args.as_object_mut() {
        m.insert("options".into(), json!([typed, given].concat()));
    }
    Ok(args)
}

/// The tools that write a directory with a long run of faba's, which may run
/// as a job ([`run_out`]).
const PRODUCERS: [&str; 9] = [
    "count", "depth", "atoi", "dartseq", "apa", "snp", "run", "qc", "cluster",
];

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
        "dry_run": { "type": "boolean", "description": "Count what the thresholds remove, per batch and modality; write nothing (`out` is ignored)" },
        "plots": { "type": "boolean", "description": "Also draw each modality's and knob's figure under `qc_plots/` and return them as images" },
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
        let a = Arguments(with_typed(name, args)?);
        let v = match name {
            "docs" => return docs(&a, cx),
            "count" => count(&a, cx)?,
            "depth" => depth(&a, cx)?,
            "atoi" => atoi(&a, cx)?,
            "dartseq" => dartseq(&a, cx)?,
            "apa" => apa(&a, cx)?,
            "snp" => snp(&a, cx)?,
            "run" => run(&a, cx)?,
            "qc" => return qc(&a, cx),
            "cluster" => cluster(&a, cx)?,
            "cluster_view" | "modules" => return cluster_view(&a, cx),
            "pileup" => return pileup(&a, cx),
            "metagene" => return metagene(&a, cx),
            "pwm" => return pwm(&a, cx),
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
            "--figure",
            "--plots",
            "--dry-run",
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
                !TAKEN.contains(&name),
                "{o} is not taken in `options`: outputs go to `out`, and no view opens here"
            );
            ensure!(
                !also.contains(&name),
                "{o} is not taken in `options`: this tool has its own argument for it"
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
    if let Some(results) = recorded_results(path) {
        v.insert("results".into(), results);
    }
    Value::Object(v)
}

/// Largest run-record `results` returned inline, in bytes of JSON.
const RESULTS_INLINE: usize = 16 * 1024;

/// The `results` of every run record in the directory `dir`, by job: what a
/// run reports beyond its files (a fit's verdict, its cluster sizes), which
/// the child logs to stderr and the call would otherwise only show as
/// progress.
fn recorded_results(dir: &str) -> Option<Value> {
    let mut m = Map::new();
    for e in std::fs::read_dir(dir).ok()?.flatten() {
        let name = e.file_name().to_string_lossy().into_owned();
        let Some(job) = name.strip_suffix(".run.json") else {
            continue;
        };
        let Some(results) = std::fs::read_to_string(e.path())
            .ok()
            .and_then(|t| serde_json::from_str::<Value>(&t).ok())
            .and_then(|mut v| v.get_mut("results").map(Value::take))
            .filter(|r| !r.is_null())
        else {
            continue;
        };
        let size = results.to_string().len();
        let results = if size <= RESULTS_INLINE {
            results
        } else {
            json!(format!(
                "{size} bytes: read `results` in {}",
                absolute(&e.path().to_string_lossy())
            ))
        };
        m.insert(job.into(), results);
    }
    (!m.is_empty()).then_some(Value::Object(m))
}

/// The prefix `--figure` takes beside the table `table`, after checking the
/// PDF and PNG it writes are new; both are added to `outputs`.
fn figure_prefix(table: &str, outputs: &mut Vec<PathBuf>) -> Result<(String, String)> {
    let prefix = table.strip_suffix(".tsv").unwrap_or(table);
    // Exactly the files a figure saved under `prefix` writes.
    let [pdf, png] = crate::figure::saved_paths(prefix);
    fresh(&pdf)?;
    fresh(&png)?;
    outputs.extend([PathBuf::from(&pdf), PathBuf::from(&png)]);
    Ok((prefix.to_string(), png))
}

/// `v`, with the figure `png` (and its PDF beside it) named in it and
/// attached as an image.
fn with_figure(mut v: Value, png: Option<&str>) -> Result<Vec<Content>> {
    let Some(png) = png else {
        return Ok(vec![Content::Json(v)]);
    };
    let pdf = format!("{}.pdf", png.strip_suffix(".png").unwrap_or(png));
    if let Some(m) = v.as_object_mut() {
        m.insert("figure".into(), json!([absolute(&pdf), absolute(png)]));
    }
    let image = std::fs::read(png).with_context(|| format!("reading {png}"))?;
    Ok(vec![Content::Json(v), Content::Png(image)])
}

/// Widest image returned for an SVG, in pixels.
const IMAGE_WIDTH: f64 = 1600.0;

/// The SVG at `path` as a PNG at most [`IMAGE_WIDTH`] pixels wide.
fn svg_png(path: &str) -> Result<Vec<u8>> {
    let svg = std::fs::read_to_string(path)?;
    let (w, _) = crate::figure::svg_size(&svg);
    let scale = (IMAGE_WIDTH / w.max(1.0)).min(2.0) as f32;
    crate::figure::png(&svg, scale)
}

/// Most `qc_plots/` images attached to one `qc` call.
const QC_IMAGES: usize = 8;

/// A TSV with a header line as one object per row.
fn tsv_rows(text: &str) -> Value {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return json!([]);
    };
    let keys: Vec<&str> = header.split('\t').collect();
    let rows: Vec<Value> = lines
        .filter(|l| !l.is_empty())
        .map(|l| {
            let m: Map<String, Value> = keys
                .iter()
                .zip(l.split('\t'))
                .map(|(k, x)| {
                    let x = x.parse::<i64>().map_or_else(|_| json!(x), |n| json!(n));
                    (k.to_string(), x)
                })
                .collect();
            Value::Object(m)
        })
        .collect();
    json!(rows)
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

/// Run `argv` so that it writes `out` (made of `what`) with `-o`; with
/// `background`, as a job whose id comes back at once.
fn run_out(
    a: &Arguments,
    cx: &Call,
    mut argv: Vec<String>,
    out: &str,
    what: &str,
) -> Result<Value> {
    fresh(out)?;
    argv.extend(["-o".into(), out.into()]);
    if a.flag("background")? {
        let job = cx.spawn_job(command(&argv)?, &[PathBuf::from(out)], what)?;
        return Ok(json!({
            "job": job,
            "written": absolute(out),
            "note": "running as a background job: job_wait or job_status until it is done, then use `written`",
        }));
    }
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
    run_out(a, cx, argv, &out, &what)
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
    run_out(a, cx, argv, &out, &what)
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
    run_out(a, cx, argv, &out, &what)
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
    run_out(a, cx, argv, &out, &what)
}

fn apa(a: &Arguments, cx: &Call) -> Result<Value> {
    let bam = a.paths_req("bam")?;
    let gff = a.path_opt("gff")?;
    let utr_bed = a.path_opt("utr_bed")?;
    let opts = a.options(&["-u", "--utr-bed"])?;
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
    if let Some(method) = a.str("method") {
        ensure!(
            method == "mixture" || method == "simple",
            "`method` is mixture or simple, not {method}"
        );
        push_flag(&mut argv, "--method", method);
    }
    argv.extend(opts);
    let what = format!("alternative polyadenylation of {}", bam.join(", "));
    run_out(a, cx, argv, &out, &what)
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
    run_out(a, cx, argv, &out, &what)
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
    run_out(a, cx, argv, &out, &what)
}

fn qc(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let input = a.path("input")?;
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
    if a.flag("dry_run")? {
        argv.push("--dry-run".into());
        let printed = cx.run(command(&argv)?)?;
        return Ok(vec![Content::Json(json!({
            "dry_run": tsv_rows(&printed),
            "note": "nothing written; `removes` counts every site failing a threshold, `only_this` those failing it alone",
        }))]);
    }
    let plots = a.flag("plots")?;
    if plots {
        argv.push("--plots".into());
    }
    argv.push("--batch-process".into());
    let out = a.out_or(cx, (&named(&input, "qc"), ""), exists)?;
    let what = format!("qc of {input}");
    let v = run_out(a, cx, argv, &out, &what)?;
    let mut images: Vec<PathBuf> = std::fs::read_dir(Path::new(&out).join("qc_plots"))
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|x| x == "png"))
        .collect();
    images.sort();
    let mut content = vec![Content::Json(v)];
    for p in images.iter().take(QC_IMAGES) {
        content.push(Content::Png(std::fs::read(p)?));
    }
    Ok(content)
}

/// `cluster`'s arguments of the MCMC fits faba had before 0.19, refused as
/// `faba cluster` refuses their flags, rather than dropped: none is in the
/// schema (`model` blocks, the model left, is taken).
fn refuse_removed(a: &Arguments) -> Result<()> {
    let mut removed: Vec<&str> = (a.0.as_object().into_iter().flatten())
        .filter(|(_, v)| !v.is_null())
        .filter(|(k, v)| {
            let value = v.as_str().map_or_else(|| v.to_string(), str::to_string);
            crate::cluster::refused(&format!("--{}", k.replace('_', "-")), &value)
        })
        .map(|(k, _)| k.as_str())
        .collect();
    removed.sort_unstable();
    ensure!(
        removed.is_empty(),
        "`{}` {}",
        removed.join("`, `"),
        crate::cluster::REMOVED_SAID
    );
    Ok(())
}

fn cluster(a: &Arguments, cx: &Call) -> Result<Value> {
    // Before anything is resolved (`out` may make today's folder).
    refuse_removed(a)?;
    let data = a.paths("data")?;
    // The second round reads the matrices beside the count fit's when none
    // are named; `--within` itself is among the typed options.
    let within = a.path_opt("within")?;
    let (from, tool, what) = match (&within, data.first()) {
        (Some(record), first) => (
            first.unwrap_or(record),
            "within",
            format!("clusters within the cell modules of {record}"),
        ),
        (None, Some(first)) => (first, "cluster", format!("clusters of {}", data.join(", "))),
        (None, None) => bail!("`data` is required (or `within`, a count fit's run record)"),
    };
    let out = a.out_or(cx, (&named(from, tool), ""), exists)?;
    let mut argv = vec!["cluster".into()];
    argv.extend(data.iter().cloned());
    push_opt(&mut argv, "--modality", a.str("modality"));
    push_opt(&mut argv, "--fit", a.str("fit"));
    argv.extend(a.numbers(&[
        "cell_roles",
        "gene_roles",
        "chains",
        "svi_epochs",
        "batch_size",
        "epochs",
        "tol",
    ])?);
    argv.extend(a.options(&[])?);
    run_out(a, cx, argv, &out, &what)
}

fn cluster_view(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.path("data")?;
    let html = a.path_opt("html")?;
    // What cannot be opened as it is is refused before the command runs;
    // the pages are those `html` writes, a within record's one a fit.
    let pages = crate::cluster::opens(
        &[data.as_str().into()],
        a.str("modality"),
        a.str("group"),
        html.as_deref().map(Path::new),
    )?;
    let mut argv = vec!["cluster-view".into(), data.clone()];
    push_opt(&mut argv, "--modality", a.str("modality"));
    push_opt(&mut argv, "--group", a.str("group"));
    push_opt(&mut argv, "--batch-groups", a.str("batch_groups"));
    for file in ["gaf", "obo", "gmt"] {
        if let Some(path) = a.path_opt(file)? {
            push_flag(&mut argv, &format!("--{file}"), &path);
        }
    }
    argv.extend(a.numbers(&["go_min_weight"])?);
    if let Some(gtf) = a.path_opt("gtf")? {
        push_flag(&mut argv, "--gtf", &gtf);
    }
    // A fitting option is refused by the command itself, naming `faba cluster`.
    argv.extend(a.options(&[])?);
    argv.push("--batch-process".into());
    push_opt(&mut argv, "--module", a.str("module"));
    push_opt(&mut argv, "--gene", a.str("gene"));
    for level in ["cells", "batches", "go", "trace", "settings"] {
        if a.flag(level)? {
            argv.push(format!("--{level}"));
        }
    }
    if let Some(figure) = a.path_opt("figure")? {
        let mut outputs = Vec::new();
        // The prefix checked is the one written under.
        let (prefix, png) = figure_prefix(&figure, &mut outputs)?;
        push_flag(&mut argv, "--figure", &prefix);
        push_opt(&mut argv, "--plot", a.str("plot"));
        push_opt(&mut argv, "--look", a.str("look"));
        argv.extend(a.numbers(&["width", "height"])?);
        let what = format!("a figure of {data}");
        cx.run_making(command(&argv)?, &outputs, &what)?;
        return with_figure(json!({}), Some(&png));
    }
    if let Some(html) = html {
        for page in &pages {
            fresh(&page.to_string_lossy())?;
        }
        push_flag(&mut argv, "--html", &html);
        let what = format!("modules of {data} as a page");
        cx.run_making(command(&argv)?, &pages, &what)?;
        let written: Vec<String> = pages
            .iter()
            .map(|p| absolute(&p.to_string_lossy()))
            .collect();
        return Ok(vec![Content::Json(match &written[..] {
            [page] => json!({ "written": page }),
            pages => json!({ "written": pages[0], "pages": pages }),
        })]);
    }
    let text = cx.run_reading(command(&argv)?)?;
    Ok(vec![Content::Text(text)])
}

/// Flags that switch `pileup` into the Miami figure. They are arguments of
/// their own, because the figure writes an SVG and does not write the TSV.
const PILEUP_TAKEN: &[&str] = &["--gtf", "--bam", "--svg", "--png", "--format", "--no-pdf"];

/// `faba pileup`'s command line for a call, and whether it draws the
/// figure.
fn pileup_argv(a: &Arguments) -> Result<(Vec<String>, bool)> {
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

fn pileup(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let data = a.paths_req("data")?;
    let (mut argv, figure) = pileup_argv(a)?;
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

fn metagene(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
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

fn pwm(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
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
    fresh(&tsv)?;
    push_flag(&mut argv, "-o", &tsv);
    let mut outputs = vec![PathBuf::from(&tsv)];
    let png = figure_arg(a, &tsv, &mut argv, &mut outputs)?;
    let what = format!("the PWM around sites of {}", data.join(", "));
    let printed = cx.run_making(command(&argv)?, &outputs, &what)?;
    with_figure(text_outcome(&tsv, &printed), png.as_deref())
}

fn exists(path: &str) -> bool {
    Path::new(path).exists()
}

#[cfg(test)]
mod tests {
    use super::{
        pileup_argv, recorded_results, tools, tsv_rows, typed_flags, with_typed, Arguments, TYPED,
    };
    use clap::CommandFactory;
    use serde_json::json;
    use std::path::{Path, PathBuf};

    /// Every typed argument names a flag its command has, and its tool's
    /// schema offers it.
    #[test]
    fn typed_arguments_are_flags_of_their_commands() {
        let cli = crate::Cli::command();
        let specs = tools();
        for &(tools, key, _, _) in TYPED {
            let long = key.replace('_', "-");
            for tool in tools {
                // A tool is named after its command, `_` for `-`.
                let sub = cli
                    .find_subcommand(tool.replace('_', "-"))
                    .unwrap_or_else(|| panic!("no command {tool}"));
                assert!(
                    sub.get_arguments()
                        .any(|a| a.get_long() == Some(long.as_str())),
                    "`faba {tool}` has no --{long}"
                );
                let spec = specs.iter().find(|s| s.name == *tool).expect("a tool");
                assert!(
                    spec.input_schema["properties"].get(key).is_some(),
                    "{tool}: no `{key}` in its schema"
                );
            }
        }
    }

    #[test]
    fn typed_arguments_become_flags_ahead_of_options() {
        let args = json!({
            "data": ["d"], "seed": 7, "within_min_cells": 20,
            "options": ["--epochs", "10"],
        });
        let a = Arguments(with_typed("cluster", args).unwrap());
        let argv = a.options(&[]).unwrap();
        assert_eq!(
            argv,
            ["--seed", "7", "--within-min-cells", "20", "--epochs", "10"]
        );
        let flags = |tool, args| typed_flags(tool, &Arguments(args));
        let e = flags("metagene", json!({ "isoforms": "sometimes" })).unwrap_err();
        assert!(e.to_string().contains("one of weighted, all"), "{e}");
        assert!(flags("cluster", json!({ "seed": "7" })).is_err());
        assert!(flags("pileup", json!({ "cell_membership": "-m.tsv" })).is_err());
        assert_eq!(
            flags("pileup", json!({ "track": ["wt=*_wt_*", "mut=*_mut_*"] })).unwrap(),
            ["--track", "wt=*_wt_*", "--track", "mut=*_mut_*"]
        );
        assert_eq!(
            flags(
                "run",
                json!({ "skip_apa": true, "skip_snp": false, "mixture": true })
            )
            .unwrap(),
            ["--mixture", "--skip-apa"]
        );
        // Not this tool's: left alone.
        assert!(flags("depth", json!({ "gene_type": "lncRNA" }))
            .unwrap()
            .is_empty());
    }

    #[test]
    fn a_dry_run_table_reads_as_rows() {
        let rows = tsv_rows(
            "level	name	removes	only_this
cells	b1	3	
sites	m6A	10	2
",
        );
        assert_eq!(
            rows,
            json!([
                { "level": "cells", "name": "b1", "removes": 3, "only_this": "" },
                { "level": "sites", "name": "m6A", "removes": 10, "only_this": 2 },
            ])
        );
    }

    /// A run's results come back with its files, by job; a record too
    /// large to quote is pointed to.
    #[test]
    fn run_record_results_come_back() {
        let dir = tempfile::tempdir().unwrap();
        let d = dir.path();
        std::fs::write(
            d.join("cluster.run.json"),
            json!({ "job": "cluster", "results": { "blockmodel": { "rhat_log_joint": 1.01 } } })
                .to_string(),
        )
        .unwrap();
        std::fs::write(d.join("qc.run.json"), json!({ "job": "qc" }).to_string()).unwrap();
        let big = "x".repeat(super::RESULTS_INLINE + 1);
        std::fs::write(
            d.join("big.run.json"),
            json!({ "results": big }).to_string(),
        )
        .unwrap();
        let r = recorded_results(&d.to_string_lossy()).unwrap();
        assert_eq!(r["cluster"]["blockmodel"]["rhat_log_joint"], json!(1.01));
        assert!(r.get("qc").is_none());
        assert!(r["big"].as_str().unwrap().contains("big.run.json"));
    }

    /// Every blockmodel option, under any of its flags, passes through
    /// `options` to `faba cluster-view`, which refuses it as clap does,
    /// saying to run `faba cluster` with it (a viewer never fits), or, for
    /// the MCMC fits' flags, that they are gone.
    #[test]
    fn cluster_view_refuses_every_option_that_would_fit() {
        let removed = [
            &["--model", "link"][..],
            &["--prior-strength", "2"],
            &["--shift", "on"],
            &["--global", "on"],
            &["--profile-strength", "inf"],
            &["--community-gamma", "3"],
            &["--sweeps", "50"],
            &["--burn-in", "10"],
        ];
        for option in [
            &["--cell-roles", "4"][..],
            &["--communities", "4"],
            &["--gene-roles", "3"],
            &["--membership-alpha=2"],
            &["--gene-alpha", "1"],
            &["--shift-sd", "1"],
            &["--cell-sd", "1"],
            &["--chains", "2"],
            &["--epochs", "5"],
        ]
        .into_iter()
        .chain(removed)
        {
            let a = Arguments(json!({ "options": option }));
            let passed = a.options(&[]).unwrap();
            let line = ["faba", "cluster-view", "out/", "--batch-process"];
            let all = line.iter().copied().map(String::from).chain(passed);
            let err = crate::matches(&mut crate::faba_command(), all)
                .unwrap_err()
                .to_string();
            let said = match removed.contains(&option) {
                true => err.contains("removed in faba 0.19"),
                false => err.contains("never fits") && err.contains("`faba cluster`"),
            };
            assert!(said, "{option:?}: {err}");
        }
    }

    /// `cluster_view` opens a fit's record or its tables, never anything
    /// it would have to fit first: another command's record or site
    /// matrices without tables beside them, or a directory without tables,
    /// are refused.
    #[test]
    fn cluster_view_opens_only_what_needs_no_fit() {
        let tmp = tempfile::tempdir().unwrap();
        let d = tmp.path();
        let at = |name: &str| d.join(name).to_string_lossy().into_owned();
        let opens = |data: String, group| crate::cluster::opens(&[data.into()], None, group, None);
        std::fs::write(at("pipeline.run.json"), r#"{"job":"run"}"#).unwrap();
        std::fs::write(at("x.zarr"), "").unwrap();
        for refused in [at("pipeline.run.json"), at("x.zarr"), at("")] {
            let err = opens(refused.clone(), None).unwrap_err().to_string();
            assert!(err.contains("run `faba cluster` first"), "{refused}: {err}");
        }
        std::fs::write(at("m6a_block_sites.parquet"), "").unwrap();
        assert!(opens(at(""), None).unwrap().is_empty());
        let fit = r#"{"job":"cluster","results":{"model":{},"tables":{"modality":"m6a"}}}"#;
        std::fs::write(at("m6a_blockmodel.run.json"), fit).unwrap();
        assert!(opens(at("m6a_blockmodel.run.json"), None)
            .unwrap()
            .is_empty());
        // One fit has no groups to name.
        let err = opens(at("m6a_blockmodel.run.json"), Some("C1")).unwrap_err();
        assert!(err.to_string().contains("within record"), "{err}");
    }

    /// A within record opened whole names its fitted groups, whose pages
    /// `html` writes beside the count fit's; with `group` it is one fit,
    /// refused when that group has none.
    #[test]
    fn cluster_view_opens_a_within_record_whole_or_by_group() {
        let tmp = tempfile::tempdir().unwrap();
        let d = tmp.path();
        let sub = d.join("within/m6a/C1");
        std::fs::create_dir_all(&sub).unwrap();
        let fit = r#"{"job":"cluster","results":{"model":{},"tables":{"modality":"m6a"}}}"#;
        std::fs::write(sub.join("m6a_blockmodel.run.json"), fit).unwrap();
        std::fs::write(sub.join("m6a_block_sites.parquet"), b"").unwrap();
        let record = json!({ "job": "cluster", "results": { "within": {
            "modality": "m6a", "count_record": d.join("count_blockmodel.run.json"),
            "groups": [
                { "group": "C1", "cells": 80, "status": "fitted", "dir": "within/m6a/C1",
                  "record": "m6a_blockmodel.run.json" },
                { "group": "C2", "cells": 9, "status": "skipped", "why": "too few" },
            ],
        } } });
        let path = d.join("m6a_within.run.json").to_string_lossy().into_owned();
        std::fs::write(&path, record.to_string()).unwrap();
        let html = Path::new("out/fit.html");
        let opens = |group| crate::cluster::opens(&[path.as_str().into()], None, group, Some(html));
        assert_eq!(
            opens(None).unwrap(),
            [PathBuf::from("out/fit.html"), "out/fit.C1.html".into()]
        );
        assert_eq!(opens(Some("c1")).unwrap(), [html]);
        let err = opens(Some("C2")).unwrap_err().to_string();
        assert!(err.contains("C2 has no fit to open"), "{err}");
        let err = opens(Some("C9")).unwrap_err().to_string();
        assert!(err.contains("no group C9"), "{err}");
    }

    /// The second round needs no `data`: its matrices are those beside the
    /// count fit's, and its output is named after the count fit's record.
    #[test]
    fn cluster_within_a_count_fit_needs_no_data() {
        let args = json!({ "within": "out/count_blockmodel.run.json", "within_min_cells": 20 });
        let a = Arguments(with_typed("cluster", args).unwrap());
        assert_eq!(
            a.options(&[]).unwrap(),
            [
                "--within",
                "out/count_blockmodel.run.json",
                "--within-min-cells",
                "20"
            ]
        );
        assert!(tools().iter().any(|t| {
            t.name == "cluster" && t.input_schema["properties"]["within"].is_object()
        }));
    }

    /// The MCMC fits' arguments are refused by name, `model` blocks taken;
    /// the fit's own are not touched.
    #[test]
    fn cluster_refuses_the_mcmc_arguments() {
        let refused = |args| super::refuse_removed(&Arguments(args)).map_err(|e| e.to_string());
        let err = refused(json!({ "data": ["d"], "sweeps": 100, "model": "link", "thin": 2 }))
            .unwrap_err();
        assert!(err.starts_with("`model`, `sweeps`, `thin` "), "{err}");
        assert!(err.contains("removed in faba 0.19"), "{err}");
        assert!(refused(json!({ "data": ["d"], "model": "blocks", "epochs": 5 })).is_ok());
        assert!(refused(json!({ "profile_strength": 5 })).is_err());
        assert!(
            refused(json!({ "sweeps": null, "holdout": null })).is_ok(),
            "null: not given"
        );
    }

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
