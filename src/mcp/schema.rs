//! The tools' JSON schemas: the hand-written arguments, and one made
//! from every other flag of a tool's command, read from its clap
//! definition.

use legume_mcp::ToolSpec;
use serde_json::{json, Map, Value};
use std::collections::BTreeMap;
use std::sync::OnceLock;

fn schema(properties: Value, required: &[&str]) -> Value {
    json!({ "type": "object", "properties": properties, "required": required })
}

/// Every tool: its hand-written arguments and one for every other flag of
/// its command (see [`generate`]).
pub(super) fn tools() -> Vec<ToolSpec> {
    all_tools().to_vec()
}

/// [`tools`], made once.
fn all_tools() -> &'static [ToolSpec] {
    static MADE: OnceLock<Vec<ToolSpec>> = OnceLock::new();
    MADE.get_or_init(|| {
        let mut specs = hand_written();
        for s in specs.iter_mut() {
            for g in generated(&s.name) {
                s.input_schema["properties"][g.key.as_str()] = g.schema.clone();
            }
        }
        specs
    })
}

/// The arguments `tool` takes, by key; none for a tool there is not.
pub(super) fn properties(tool: &str) -> Option<&'static Map<String, Value>> {
    (all_tools().iter())
        .find(|s| s.name == tool)
        .and_then(|s| s.input_schema["properties"].as_object())
}

/// The tools as written here: their inputs, outputs and common flags.
fn hand_written() -> Vec<ToolSpec> {
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
            "description": "Anything else for the command, as `faba <command> --help` lists it; every flag is an argument of its own already, and one given both ways is refused. A short option stands apart from its value (`-g X`)"
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
            "Filter a faba directory into a new one (`faba qc --batch-process`). Cells are cut on `{batch}_count` and the same cells are kept in every matrix; sites are cut on the thresholds given (marginal p-values, no FDR). Omit a threshold to keep faba's default. Try thresholds first with `dry_run`: it returns what the cell QC keeps per batch and what each threshold removes per modality, writing nothing. The dry run's `pairs` are the pairs of thresholds Auto picks to show two at a time: those splitting the most sites (pass one, fail the other) on axes that agree least. With `plots`, the view's figures go to `qc_plots/`: each threshold's histogram, each modality's thresholds two at a time ({m}_pairs: kept sites against removed, the box where both pass; not drawn with `panel_only`, whose pairs are the panels) and each batch's cell QC ({batch}_cells); and the site panels `site_plot` draws go to `qc_plots/panels/`. Both come back as images, at most 8: the view's figures first, then the panels, each in name order, half the room each and the rest to whichever has more. With `html`, a page to share goes to that file (dry run too); a file that exists is refused, never replaced: ask the user for another path.",
            schema(qc_properties(&dir(), &out(), &options()), &["input"]),
        ),
        ToolSpec::new(
            "site_plot",
            "Draw a faba directory's site statistics as paper panels (`faba site-plot`), no QC: for one `modality`, a histogram of each `hist` quantity (pv, fold, log_odds, coverage, converted, edit_ratio, cells; default every one that applies) and each of `pairs` two at a time (\"x:y\", e.g. \"pv:cells\"; default the pairs `qc` picks), with `marginals` (top, right) on the pair's own axes. Each quantity has one scale, range and set of bins in every figure (`scale`: \"pv=log\", \"fold=linear\"; defaults: -log10 p, edit ratio and log odds linear, fold and counts log). cells counts the cells with a converted read at a site, over every cell. Thresholds are optional: given ones are drawn as dashed lines and the sites they remove drawn apart; none given, every site is kept. Writes only figures to `out`, a new directory: {modality}_hist_{q} and {modality}_{y}_vs_{x}, each .pdf (labels as text) and .png; returns their paths and the PNGs as images. `qc` with `plots` writes the same panels under qc_plots/panels/.",
            schema(site_plot_properties(&dir(), &out(), &options(), &modality()), &["input"]),
        ),
        ToolSpec::new(
            "qc_view",
            "Read back what a `qc` run kept and dropped (`faba qc-view --batch-process`), cutting nothing: `input` is the qc output directory, whose qc.run.json gives the thresholds it applied; its kept sites and its dropped ones (with their reason) are read together, each site's cells counted as qc counted them (off the input it read, when that is still there; else the dropped sites read 0 cells, as it warns). Any other faba directory opens too, with no thresholds: every site kept. Returns the table `qc`'s `dry_run` returns, at the thresholds applied: `table`, what the cell QC kept per batch (off its reports) and what each threshold removes per modality (`removes` every site failing it, `only_this` those failing it alone), and `pairs`, the pairs Auto picks to show two at a time. With `plots`, the figures `qc` draws with `plots` go to `out`/qc_plots/ (the site panels to qc_plots/panels/) and come back as images, as `qc`'s do; with `html`, the page to share goes to that file. Neither may exist: an existing file or directory is refused, never replaced: ask the user for another path.",
            schema(qc_view_properties(&dir(), &options()), &["input"]),
        ),
        ToolSpec::new(
            "cluster",
            "Cluster cells by their site reads (`faba cluster`). `data` is a faba directory or site matrices. The model is a blockmodel of cell modules and gene modules, a methylation shift per pair and each cell's own level, fitted by variational inference from k-means starts (`chains` restarts, the best by ELBO kept; `svi_epochs`, `batch_size`, `epochs`, `tol`). The MCMC fits of faba before 0.19 are gone: `model` link or crp, `sweeps`, `burn_in`, `shift`, `profile_strength` and `holdout` are refused. Tables go to `out`, a new directory, not back into `data`. With `modality` count, cluster cells by their gene counts instead (count matrices), fitted alike: its verdict is whether the best restart's ELBO settled. With `within`, the run record of such a count fit (count_blockmodel.run.json), fit the sites within each of its cell modules, one fit a module, into `out`/within/{m}/{group}/, with `out`/{m}_within.run.json naming them all and {m}_within_cells.parquet giving each cell its nested label (C3.2); `data` may then be left out, and the matrices are those beside the count fit's inputs. The sites' batches pair with the count fit's by name (the same, or one named with its folder: `b1` and `runs/a/b1`); when that is ambiguous the run is refused, listing each choice: ask the user which is meant, then run again with `batch_map`. While it fits, a progress line every 15 s says the restart and epoch, the ELBO per read (or count) and its last move against `tol`, the share of steps done and about how long is left (an upper bound: a restart that settles ends early); with `background`, job_status shows it as a `progress` stage (steps done of all, seconds taken and left) and the latest line in `last_lines`.",
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
                    "batch_map": {
                        "description": "With `within`: which of the count fit's batches each of the sites' batches is, taken before any pairing by name. An object {\"SITES\": \"COUNT\"} or a list of \"SITES=COUNT\"; an empty COUNT leaves that batch unmatched. Give it when the run is refused as ambiguous: the refusal lists each choice",
                        "anyOf": [
                            { "type": "object", "additionalProperties": { "type": "string" } },
                            { "type": "array", "items": { "type": "string" } },
                        ],
                    },
                    "out": out(), "options": options(),
                }),
                &[],
            ),
        ),
        ToolSpec::new(
            "cluster_view",
            "View a `cluster` fit (`faba cluster-view --batch-process`), level by level as the browser does: by default the gene modules against the cell modules (level 1). `module` prints one gene module's genes and sites with their values in each cell module (level 2); `gene` one gene's sites in each cell module over every gene module (levels 3 and 4); `cells` the cell modules' sizes; `go` the GO terms (or `gmt` gene sets) over-represented in each gene module, against the fit's own annotated genes (files default to lupin's cache); `batches` how each cell module spreads over the batches (shared, one condition's, or batch-specific: uneven over a condition's replicates), the conditions read off the batch names less a replicate token or given by `batch_groups`; `trace` each restart's ELBO by epoch, the verdict and its advice; with `figure`, a figure instead (`plot`: the modules grid coloured by `look`, the trace, or the batches; `width` and `height` in inches lay it out for a page), returned as an image; `settings` the settings the fit ran with. Several may be given; rows are tab-separated. Values are % converted for a site fit, UMIs per cell for a count fit. `data` is the fit's run record or the directory holding its tables; this never fits, and the record of an MCMC fit of faba before 0.19 is refused. A within record ({m}_within.run.json, from `cluster` with `within`) prints the count fit's modules, each group's fit (cells, verdict; why one was skipped) and the modules within each fitted group; with `group`, every argument applies to the fit within that cell module, and without it the levels asked for are the count fit's. With `html`, write one self-contained page to a new file instead; a within record without `group` writes a page a fit, `html` for the count fit and {stem}.{group}.html beside it for each group, linked to each other.",
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
                    "go_q": { "type": "number", "exclusiveMinimum": 0, "maximum": 1, "description": "A term (or `gmt` set) is listed for `go` at this BH q or less (default 0.05)" },
                    "batch_groups": { "type": "string", "description": "Which condition each batch is a replicate of, as rep1_mut=mut,rep2_mut=mut; the rest are read off their names" },
                    "trace": { "type": "boolean", "description": "Print each restart's trace, the verdict and its advice" },
                    "settings": { "type": "boolean", "description": "Print the settings the fit ran with" },
                    "gtf": { "type": "string", "description": "Gene annotation for the models drawn in the HTML page" },
                    "html": { "type": "string", "description": "Write the page here; must not exist" },
                    "figure": { "type": "string", "description": "Write a figure as {figure}.pdf and {figure}.png (neither may exist) and return it as an image" },
                    "plot": { "type": "string", "enum": ["modules", "trace", "batches"], "description": "Which figure: the modules grid (gene modules against cell modules; default), the fit's trace, or the cell modules against the batches" },
                    "look": { "type": "string", "enum": ["shift", "percent", "z"], "description": "What colours the modules grid's blocks: the fitted shift (a count fit's log rate; default), the methylated % (UMIs per cell), or that standardized across the cell modules" },
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

/// The producers that read cell barcodes from BAM tags.
const TAGGED: &[&str] = &["count", "depth", "atoi", "dartseq", "apa", "snp", "run"];

/// The producers that also count genes and UMIs.
const TAGGED_GENES: &[&str] = &["count", "atoi", "dartseq", "apa", "snp", "run"];

/// Flags whose arguments say more than their help: the tools taking each,
/// its key and description. Each is an argument as every other flag
/// is (see [`generate`]); this only words it.
const TYPED: &[(&[&str], &str, &str)] = &[
    (
        TAGGED,
        "cell_barcode_tag",
        "BAM tag holding the cell barcode",
    ),
    (
        TAGGED_GENES,
        "gene_barcode_tag",
        "BAM tag holding the gene",
    ),
    (
        TAGGED_GENES,
        "umi_tag",
        "BAM tag holding the UMI",
    ),
    (
        TAGGED_GENES,
        "gene_type",
        "Keep genes of this biotype only (protein_coding, lncRNA, pseudogene)",
    ),
    (
        &["atoi", "dartseq", "apa", "run"],
        "mixture",
        "Also fit the per-gene mixture matrices (a slow EM)",
    ),
    (&["run"], "skip_snp", "Skip the SNP step"),
    (
        &["run"],
        "skip_count",
        "Skip the gene counts",
    ),
    (&["run"], "skip_atoi", "Skip A-to-I"),
    (&["run"], "skip_apa", "Skip APA"),
    (
        &["run"],
        "depth_resolution_kb",
        "Read-depth bin width in kb",
    ),
    (
        &["cluster"],
        "seed",
        "Random seed; restart i uses seed + i (default 1)",
    ),
    (
        &["cluster"],
        "within",
        "A count fit's run record (count_blockmodel.run.json): fit the sites within each of its cell modules",
    ),
    (
        &["cluster"],
        "within_min_cells",
        "With `within`: a group with fewer cells to fit is skipped (default 50)",
    ),
    (
        &["pileup"],
        "cell_membership",
        "Cell barcode to cell type (TSV/CSV/Parquet): one panel per type",
    ),
    (
        &["pileup"],
        "top_modality",
        "Restrict the top track to these modalities, e.g. m6A",
    ),
    (
        &["pileup"],
        "track",
        "Labelled tracks as LABEL=PATTERN over the input paths (`*`, `?`)",
    ),
    (
        &["metagene"],
        "isoforms",
        "A site in k isoforms adds 1/k to each (weighted, the default) or 1 (all)",
    ),
    (
        &["metagene"],
        "include_non_coding",
        "Also profile non-coding genes, on their own track",
    ),
    (
        SITE_PANELS,
        "hist",
        "Site panels: quantities to draw a histogram of (pv, fold, log_odds, coverage, converted, edit_ratio, cells); default every one that applies, [\"none\"] for none",
    ),
    (
        SITE_PANELS,
        "pairs",
        "Site panels: pairs to draw two at a time, each \"x:y\" (e.g. \"pv:cells\"); default the pairs qc picks, [\"none\"] for none",
    ),
    (
        SITE_PANELS,
        "marginals",
        "Site panels: histograms beside each pair on its own axes, \"top\" (x) and/or \"right\" (y)",
    ),
    (
        SITE_PANELS,
        "scale",
        "A quantity's scale in every figure of it, each \"q=linear\" or \"q=log\" (e.g. \"pv=log\"); log odds is linear only",
    ),
    (
        SITE_PANELS,
        "bins",
        "Site panels: bins along each quantity's axis (default from the figure's width, 20 to 60)",
    ),
    (
        PANELS,
        "font_size",
        "Every label's size in points, as printed (default: fit to the figure)",
    ),
    (
        PANELS,
        "panel_only",
        "Draw the plot, its axes and their titles only: no title, file name, key or notes",
    ),
    (
        PANELS,
        "palette",
        "Kept and removed sites in blue and orange (colour) or dark and light gray (gray)",
    ),
];

/// The tools whose site figures take a paper panel's look.
const PANELS: &[&str] = &["metagene", "pwm", "qc", "qc_view", "site_plot"];

/// The tools drawing the site panels, which take their choice of panels.
const SITE_PANELS: &[&str] = &["qc", "qc_view", "site_plot"];

/// The flags no tool offers as an argument, each for the tools named
/// (none: every tool), with why. Every other flag of a tool's command is
/// an argument of the tool. The MCMC flags `faba cluster` refuses since
/// 0.19 are hidden, so never offered, and refused by [`known_keys`](super::args::known_keys).
const EXCLUDED: &[(&[&str], &str, &str)] = &[
    (
        &[],
        "output",
        "the tool's `out`: it names the output, never one that exists",
    ),
    (
        &[],
        "batch-process",
        "the tool adds it itself: no view opens under the server",
    ),
    (
        &[],
        "no-detach",
        "opens the view in this terminal: no view opens under the server",
    ),
    (
        &[],
        "session",
        "names a view's session: no view opens under the server",
    ),
    (
        &[],
        "verbose",
        "the server runs every command with it: its log lines are the call's progress",
    ),
    (&[], "help", "prints the help and runs nothing"),
    (
        &["metagene"],
        "print",
        "the tool adds it itself: the histogram is the call's progress",
    ),
    (
        &["pileup"],
        "svg",
        "the figure is the Miami SVG the tool returns: `figure`, `gtf` and `bam` ask for it",
    ),
    (
        &["pileup"],
        "png",
        "the figure is the Miami SVG the tool returns: `figure`, `gtf` and `bam` ask for it",
    ),
    (
        &["pileup"],
        "no-pdf",
        "the figure is the Miami SVG the tool returns: `figure`, `gtf` and `bam` ask for it",
    ),
    (
        &["pileup"],
        "format",
        "the figure is the Miami SVG the tool returns: `figure`, `gtf` and `bam` ask for it",
    ),
];

/// Why `tool` does not offer `--{long}`, if it does not.
fn excluded(tool: &str, long: &str) -> Option<&'static str> {
    (EXCLUDED.iter())
        .find(|(only, l, _)| *l == long && (only.is_empty() || only.contains(&tool)))
        .map(|e| e.2)
}

/// What one value of a flag is, by its clap definition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum Ty {
    Integer { unsigned: bool },
    Number,
    Text,
    Choice(Vec<String>),
}

/// The type of `a`'s values: its possible values, else what its value
/// parser makes.
pub(super) fn value_type(a: &clap::Arg) -> Ty {
    use std::any::TypeId;
    let choices: Vec<String> = (a.get_possible_values().iter())
        .filter(|v| !v.is_hide_set())
        .map(|v| v.get_name().to_string())
        .collect();
    if !choices.is_empty() {
        return Ty::Choice(choices);
    }
    let made = a.get_value_parser().type_id();
    let any = |ids: &[TypeId]| ids.iter().any(|id| made == *id);
    let unsigned = [
        TypeId::of::<u8>(),
        TypeId::of::<u16>(),
        TypeId::of::<u32>(),
        TypeId::of::<u64>(),
        TypeId::of::<usize>(),
    ];
    let signed = [
        TypeId::of::<i8>(),
        TypeId::of::<i16>(),
        TypeId::of::<i32>(),
        TypeId::of::<i64>(),
        TypeId::of::<isize>(),
    ];
    if any(&unsigned) {
        Ty::Integer { unsigned: true }
    } else if any(&signed) {
        Ty::Integer { unsigned: false }
    } else if any(&[TypeId::of::<f32>(), TypeId::of::<f64>()]) {
        Ty::Number
    } else {
        Ty::Text
    }
}

/// How a flag is given.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum Take {
    /// The bare flag, when true.
    Switch,
    /// The bare flag, that many times.
    Count,
    /// The flag and a value, or (unless [`Many::One`]) several.
    Value(Ty),
}

/// A tool argument made from a flag of its command.
#[derive(Debug)]
pub(super) struct Generated {
    pub(super) key: String,
    /// `--{long}`.
    pub(super) flag: String,
    /// Every way to write the flag in `options`: `--{long}`, its aliases,
    /// its short name.
    pub(super) names: Vec<String>,
    pub(super) take: Take,
    pub(super) many: legume_tui::form::Many,
    /// A file: no leading `-`, which would read as a flag.
    pub(super) path: bool,
    /// The property offered.
    schema: Value,
}

/// What a tool's command offers through the tool, read from its flags.
#[derive(Debug, Default)]
struct Flags {
    /// An argument for every flag offered, in the command's order.
    generated: Vec<Generated>,
    /// Every name of a flag not offered, but the hidden ones, with why:
    /// [`EXCLUDED`], or the tool has its own argument for it.
    refused: Vec<(String, String)>,
    /// The keys of the command's hidden flags.
    hidden: Vec<String>,
}

/// A JSON schema type's name, as an error says what was wanted.
pub(super) fn wanted(schema: &Value) -> String {
    match &schema["type"] {
        Value::String(t) if t == "array" => format!("a list of {}", wanted(&schema["items"])),
        Value::String(t) if t == "integer" => "an integer".into(),
        Value::String(t) => format!("a {t}"),
        _ => "a value".into(),
    }
}

impl Generated {
    /// The property for the flag: its type, and `about` as description.
    pub(super) fn property(take: &Take, many: legume_tui::form::Many, about: &str) -> Value {
        let mut one = match take {
            Take::Switch => json!({ "type": "boolean" }),
            Take::Count => json!({ "type": "integer", "minimum": 0 }),
            Take::Value(Ty::Integer { unsigned: true }) => {
                json!({ "type": "integer", "minimum": 0 })
            }
            Take::Value(Ty::Integer { unsigned: false }) => json!({ "type": "integer" }),
            Take::Value(Ty::Number) => json!({ "type": "number" }),
            Take::Value(Ty::Text) => json!({ "type": "string" }),
            Take::Value(Ty::Choice(of)) => json!({ "type": "string", "enum": of }),
        };
        let several = matches!(take, Take::Value(_)) && many != legume_tui::form::Many::One;
        if several {
            one = json!({ "type": "array", "items": one });
        }
        one["description"] = json!(about);
        one
    }
}

/// Every tool's [`Flags`], by tool, made once.
fn all_flags() -> &'static BTreeMap<String, Flags> {
    static MADE: OnceLock<BTreeMap<String, Flags>> = OnceLock::new();
    MADE.get_or_init(generate)
}

/// `tool`'s arguments made from its command's flags, in the command's
/// order.
pub(super) fn generated(tool: &str) -> &'static [Generated] {
    all_flags().get(tool).map_or(&[], |f| &f.generated)
}

/// Why `tool` refuses `name` (`--output`, `-o`) in `options`: a flag of
/// its command it does not offer, but a hidden one.
pub(super) fn refused(tool: &str, name: &str) -> Option<&'static str> {
    let flags = all_flags().get(tool)?;
    (flags.refused.iter())
        .find(|(n, _)| n == name)
        .map(|(_, why)| why.as_str())
}

/// Whether `key` names a hidden flag of `tool`'s command.
pub(super) fn hidden(tool: &str, key: &str) -> bool {
    all_flags()
        .get(tool)
        .is_some_and(|f| f.hidden.iter().any(|h| h == key))
}

/// A flag offered under another key than its own by the tools named.
struct Renamed {
    tools: &'static [&'static str],
    long: &'static str,
    key: &'static str,
    /// What a tool's property says beyond the flag's type and help.
    patch: fn(&str) -> Value,
}

/// The flags offered under another key than their own.
const RENAMED: &[Renamed] = &[
    // `height` is the figure's.
    Renamed {
        tools: &["pileup"],
        long: "height",
        key: "ascii_height",
        patch: |_| json!({}),
    },
    Renamed {
        tools: &SIZED,
        long: "fig-width",
        key: "width",
        patch: |tool| {
            json!({ "minimum": least_size(tool).0,
                "description": "The figure's width in inches, laid out for it, its text kept readable (e.g. 3.5 for a column, 7 for a page); alone or with `height`, the other the figure's own" })
        },
    },
    Renamed {
        tools: &SIZED,
        long: "fig-height",
        key: "height",
        patch: |tool| {
            json!({ "minimum": least_size(tool).1,
                "description": "The figure's height in inches; alone or with `width`, its text fit to it. A qc figure with its table needs 10 inches across unless `panel_only` (it says so)" })
        },
    },
];

/// For each tool named after a command (`_` for `-`), an argument for
/// every flag of the command its hand-written schema does not have
/// (under the flag's key, as [`RENAMED`] may give it, or an alias's),
/// but those [`EXCLUDED`] and the hidden ones; those not offered, by
/// every name, are refused in `options`.
fn generate() -> BTreeMap<String, Flags> {
    use clap::{ArgAction, ValueHint};
    use legume_tui::form::{many_of, Many};
    let cli = crate::faba_command();
    let mut made = BTreeMap::new();
    for spec in hand_written() {
        let tool = spec.name.as_str();
        let Some(sub) = cli.find_subcommand(tool.replace('_', "-")) else {
            continue;
        };
        let flags: &mut Flags = made.entry(tool.to_string()).or_default();
        let props = &spec.input_schema["properties"];
        for a in sub.get_arguments() {
            let Some(long) = a.get_long() else {
                continue;
            };
            if a.is_hide_set() {
                flags.hidden.push(long.replace('-', "_"));
                continue;
            }
            let aliases: Vec<&str> = a.get_all_aliases().into_iter().flatten().collect();
            let mut names = vec![format!("--{long}")];
            names.extend(aliases.iter().map(|l| format!("--{l}")));
            names.extend(
                (a.get_short_and_visible_aliases().into_iter().flatten()).map(|c| format!("-{c}")),
            );
            let renamed = (RENAMED.iter()).find(|r| r.tools.contains(&tool) && r.long == long);
            let key = renamed.map_or_else(|| long.replace('-', "_"), |r| r.key.to_string());
            let own = (std::iter::once(key.clone()))
                .chain(aliases.iter().map(|l| l.replace('-', "_")))
                .find(|k| props.get(k).is_some());
            let why = match (excluded(tool, long), own) {
                (Some(why), _) => Some(why.to_string()),
                (None, Some(own)) => {
                    Some(format!("this tool has its own argument for it, `{own}`"))
                }
                (None, None) => None,
            };
            if let Some(why) = why {
                (flags.refused).extend(names.into_iter().map(|n| (n, why.clone())));
                continue;
            }
            let take = match a.get_action() {
                ArgAction::SetTrue | ArgAction::SetFalse => Take::Switch,
                ArgAction::Count => Take::Count,
                ArgAction::Set | ArgAction::Append => Take::Value(value_type(a)),
                _ => continue,
            };
            // Repeated before joined: `--f a --f b` is each value whole.
            let many = match a.get_action() {
                ArgAction::Append => Many::Repeated,
                _ => many_of(a),
            };
            let typed = TYPED.iter().find(|t| t.0.contains(&tool) && t.1 == key);
            let path = take == Take::Value(Ty::Text)
                && matches!(
                    a.get_value_hint(),
                    ValueHint::FilePath | ValueHint::DirPath | ValueHint::AnyPath
                );
            let about = match typed {
                Some(t) => t.2.to_string(),
                None => help_of(a, &take),
            };
            let mut schema = Generated::property(&take, many, &about);
            if let (Some(r), Some(m)) = (renamed, schema.as_object_mut()) {
                m.extend((r.patch)(tool).as_object().cloned().unwrap_or_default());
            }
            flags.generated.push(Generated {
                flag: format!("--{long}"),
                schema,
                key,
                names,
                take,
                many,
                path,
            });
        }
    }
    made
}

/// `a`'s help, with the value clap takes when it is not given.
fn help_of(a: &clap::Arg, take: &Take) -> String {
    let help = (a.get_help().or_else(|| a.get_long_help()))
        .map(ToString::to_string)
        .unwrap_or_default();
    let default: Vec<String> = (a.get_default_values().iter())
        .map(|v| v.to_string_lossy().into_owned())
        .filter(|v| !v.is_empty())
        .collect();
    match take {
        Take::Value(_) if !default.is_empty() && !help.contains("default") => {
            format!("{help} (default {})", default.join(","))
        }
        _ => help,
    }
}

/// The tools that write a directory with a long run of faba's, which may run
/// as a job ([`run_out`](super::run::run_out)).
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
        "plots": { "type": "boolean", "description": "Also draw each modality's and knob's figure, each modality's thresholds two at a time ({m}_pairs) and each batch's cell QC ({batch}_cells) under `qc_plots/`, and the site panels `site_plot` draws, at these thresholds, under `qc_plots/panels/` (`hist`, `pairs`, `marginals`, `scale`, `bins` choose them); both come back as images (see the tool's description)" },
        "html": { "type": "string", "description": "Also write the view's page here, to share (with `dry_run` too): the thresholds recounted live, their pairs two at a time and each batch's cell QC, carrying each modality's strongest 50,000 sites by p-value. Must not exist: an existing file is refused, never replaced" },
        "options": options,
    })
}

/// `qc_view`'s arguments: the directory in, the page and the figures out.
fn qc_view_properties(dir: &Value, options: &Value) -> Value {
    json!({
        "input": dir,
        "html": { "type": "string", "description": "Also write the view's page here, to share, at the thresholds applied (each modality's strongest 50,000 sites by p-value). Must not exist: an existing file is refused, never replaced" },
        "plots": { "type": "boolean", "description": "Also draw the figures `qc` draws with `plots`, at the thresholds applied, under `out`/qc_plots/ (`hist`, `pairs`, `marginals`, `scale`, `bins` choose the site panels under qc_plots/panels/); they come back as images" },
        "out": { "type": "string", "description": "With `plots`: the directory the figures go to; must not exist. Omitted, a new directory in today's results folder" },
        "options": options,
    })
}

/// The site thresholds `qc` and `site_plot` take, as numbers.
pub(super) const SITE_THRESHOLDS: [&str; 8] = [
    "site_max_pv",
    "site_min_log_odds",
    "site_min_coverage",
    "site_min_converted",
    "site_min_edit_ratio",
    "site_max_edit_ratio",
    "site_min_fold",
    "site_min_cells",
];

/// `site_plot`'s arguments: the directory in, the modality, the optional
/// thresholds (as `qc`'s), the directory out.
fn site_plot_properties(dir: &Value, out: &Value, options: &Value, modality: &Value) -> Value {
    let qc = qc_properties(dir, out, options);
    let mut p = json!({
        "input": dir,
        "modality": modality,
        "out": { "type": "string", "description": "Directory the figures go to; must not exist. Omitted, a new directory in today's results folder" },
        "options": options,
    });
    for key in SITE_THRESHOLDS {
        let mut t = qc[key].clone();
        t["description"] = json!(format!(
            "Optional; drawn as a dashed line, the sites it removes drawn apart. {}",
            t["description"].as_str().unwrap_or("As qc's")
        ));
        p[key] = t;
    }
    p
}

/// The tools that draw a figure, each taking its size as `width` and
/// `height` in inches, passed as `--fig-width` and `--fig-height` (see
/// [`with_typed`](super::args::with_typed)).
const SIZED: [&str; 7] = [
    "cluster_view",
    "pileup",
    "metagene",
    "pwm",
    "qc",
    "qc_view",
    "site_plot",
];

/// The least `width` and `height` a [`SIZED`] tool's figure takes, in
/// inches: pileup's and cluster_view's margins do not shrink with their
/// text (see [`crate::figure::size::FIXED_MIN`]).
fn least_size(tool: &str) -> (f64, f64) {
    match tool {
        "pileup" | "cluster_view" => crate::figure::size::FIXED_MIN,
        _ => crate::figure::size::MIN,
    }
}

#[cfg(test)]
/// The tools that run a `faba` command: named after it, `_` for `-`.
pub(super) fn commands() -> Vec<(String, clap::Command)> {
    let cli = crate::faba_command();
    tools()
        .iter()
        .filter_map(|t| {
            let sub = cli.find_subcommand(t.name.replace('_', "-"))?;
            Some((t.name.clone(), sub.clone()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::commands;
    use crate::mcp::args::{typed_flags, with_typed, Arguments};
    use clap::CommandFactory;

    use super::{excluded, generated, tools, EXCLUDED, SIZED, TYPED};
    use crate::mcp::tools::sites::pileup_argv;
    use serde_json::json;
    /// Every typed argument names a flag its command has, and its tool's
    /// schema offers it.
    #[test]
    fn typed_arguments_are_flags_of_their_commands() {
        let cli = crate::Cli::command();
        let specs = tools();
        for &(tools, key, _) in TYPED {
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

    /// Every flag of every tool's command is an argument made from it, a
    /// hand-written one (by its long name or an alias), or is in
    /// [`EXCLUDED`]: nothing reaches the command only through `options`.
    /// No two arguments of a tool share a key.
    #[test]
    fn every_flag_is_an_argument_or_excluded() {
        let specs = tools();
        let mut gaps = Vec::new();
        for (tool, sub) in commands() {
            let spec = specs.iter().find(|s| s.name == tool).unwrap();
            let props = spec.input_schema["properties"].as_object().unwrap();
            let made = generated(&tool);
            for a in sub.get_arguments().filter(|a| !a.is_hide_set()) {
                let Some(long) = a.get_long() else { continue };
                let offered = made.iter().any(|g| g.flag == format!("--{long}"));
                let key = |l: &str| l.replace('-', "_");
                let written = std::iter::once(long)
                    .chain(a.get_all_aliases().into_iter().flatten())
                    .any(|l| props.contains_key(&key(l)) && !made.iter().any(|g| g.key == key(l)));
                if !(offered || written || excluded(&tool, long).is_some()) {
                    gaps.push(format!("{tool} --{long}"));
                }
            }
            let mut keys: Vec<&str> = made.iter().map(|g| g.key.as_str()).collect();
            keys.sort_unstable();
            keys.dedup();
            assert_eq!(keys.len(), made.len(), "{tool}: two arguments share a key");
        }
        assert!(gaps.is_empty(), "no argument for:\n{}", gaps.join("\n"));
    }

    /// Each exclusion names a flag one of its tools' commands has, says
    /// why, and is offered by none of them.
    #[test]
    fn every_exclusion_is_a_flag_left_out() {
        let specs = tools();
        let commands = commands();
        for &(only, long, why) in EXCLUDED {
            assert!(!why.is_empty(), "--{long}: no reason");
            let of = |tool: &str| only.is_empty() || only.contains(&tool);
            let mut found = false;
            for (tool, sub) in commands.iter().filter(|(t, _)| of(t)) {
                found |= sub.get_arguments().any(|a| a.get_long() == Some(long));
                let spec = specs.iter().find(|s| &s.name == tool).unwrap();
                // A property by its name, if any, is made from another flag.
                let key = long.replace('-', "_");
                let made = generated(tool);
                assert!(
                    !made.iter().any(|g| g.flag == format!("--{long}"))
                        && (spec.input_schema["properties"].get(&key).is_none()
                            || made.iter().any(|g| g.key == key)),
                    "{tool} offers --{long}, excluded"
                );
            }
            assert!(found, "--{long}: no tool's command has it");
        }
    }

    /// Every figure tool takes `width` and `height` in inches, as the
    /// command's `--fig-width` and `--fig-height`.
    #[test]
    fn every_figure_tool_takes_its_size_alike() {
        for tool in tools().iter().filter(|t| SIZED.contains(&t.name.as_str())) {
            for key in ["width", "height"] {
                assert!(
                    tool.input_schema["properties"][key].is_object(),
                    "{} {key}",
                    tool.name
                );
            }
        }
        // The schema and the flags from one list: each tool's argv has them.
        for tool in SIZED {
            let a = Arguments(json!({ "width": 3.5, "height": 2 }));
            assert_eq!(
                typed_flags(tool, &a).unwrap(),
                ["--fig-width", "3.5", "--fig-height", "2"],
                "{tool}"
            );
        }
        let args = json!({ "data": ["s.zarr"], "genes": ["G1"], "figure": true, "width": 6 });
        let (argv, _) = pileup_argv(&Arguments(with_typed("pileup", args).unwrap())).unwrap();
        assert!(
            argv.windows(2).any(|w| w == ["--fig-width", "6"]),
            "{argv:?}"
        );
        assert!(typed_flags("pwm", &Arguments(json!({ "width": "wide" }))).is_err());
        // The schema states the least a figure takes: pileup's and
        // cluster_view's, whose margins do not shrink, more.
        let specs = tools();
        let least = |tool: &str, key: &str| {
            let t = specs.iter().find(|t| t.name == tool).unwrap();
            t.input_schema["properties"][key]["minimum"].clone()
        };
        let (min, fixed) = (crate::figure::size::MIN, crate::figure::size::FIXED_MIN);
        for tool in ["pwm", "metagene", "qc", "site_plot"] {
            assert_eq!(least(tool, "width"), json!(min.0), "{tool}");
            assert_eq!(least(tool, "height"), json!(min.1), "{tool}");
        }
        for tool in ["pileup", "cluster_view"] {
            assert_eq!(least(tool, "width"), json!(fixed.0), "{tool}");
            assert_eq!(least(tool, "height"), json!(fixed.1), "{tool}");
        }
    }
}
