//! The tools' JSON schemas: the hand-written arguments, and one made
//! from every other flag of a tool's command, read from its clap
//! definition.

use crate::figure::size::Figure;
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
                    "out": out(), "source": source(), "options": options(),
                }),
                &["bam"],
            ),
        ),
        ToolSpec::new(
            "snp",
            "Discover and genotype SNPs (`faba snp`). With `gff`, also writes a per-cell allele-frequency matrix; `bulk` writes the call set only; `skip_discovery` genotypes the `known_snps` only. Writes a new directory.",
            schema(
                json!({
                    "bam": bam(), "genome": genome(), "gff": gff(),
                    "known_snps": { "type": "string", "description": "Known sites (VCF/BCF/Parquet) to force-call" },
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
            "Filter a faba directory into a new one (`faba qc --batch-process`). Cells are cut on `{batch}_count` and the same cells are kept in every matrix; sites are cut on the thresholds given (marginal p-values, no FDR). Omit a threshold to keep faba's default. Try thresholds first with `dry_run`: it returns what the cell QC keeps per batch and what each threshold removes per modality, writing nothing. The dry run's `pairs` are the pairs of thresholds Auto picks to show two at a time: those splitting the most sites (pass one, fail the other) on axes that agree least. With `plots`, the view's figures go to `qc_plots/`: each threshold's histogram, each modality's thresholds two at a time ({m}_pairs: kept sites against removed, the box where both pass; not drawn with `panel_only`, whose pairs are the panels) and each batch's cell QC ({batch}_cells); and the site panels `site_plot` draws go to `qc_plots/panels/`, chosen by `hist`, `pairs`, `marginals`, `style`, `marginal_size`, `p_axis`, `front`, `key`, `scale` and `bins` as `site_plot`'s are. Both come back as images, at most 8: the view's figures first, then the panels, each in name order, half the room each and the rest to whichever has more. With `html`, a page to share goes to that file (dry run too); a file that exists is refused, never replaced: ask the user for another path.",
            schema(qc_properties(&dir(), &out(), &options()), &["input"]),
        ),
        ToolSpec::new(
            "site_plot",
            "Draw a faba directory's site statistics as paper panels (`faba site-plot`), no QC: for one `modality`, a histogram of each `hist` quantity (pv, fold, log_odds, coverage, converted, edit_ratio, cells; default every one that applies) and each of `pairs` two at a time (\"x:y\", e.g. \"pv:cells\"; default the pairs `qc` picks), with `marginals` (top, right, both) on the pair's own axes. `style` \"view\" draws each panel as the qc view draws it (a histogram as its knob's, a pair as its 2D panel): the control arm (MUT) light gray behind the signal arm (WT), WT's filtered sites white and outlined over its kept dark, MUT only where the quantity has a value on each arm (coverage, converted, edit_ratio), each threshold dashed orange-red with the words for either side; its marginals `marginal_size` of the plot's shorter side deep, `front` \"mut\" the control arm in front, `p_axis` \"neglog10\" the p-value's ticks as -log10 p (else p), `key` the histogram's key kept with `panel_only`. Each quantity has one scale, range and set of bins in every figure (`scale`: \"pv=log\", \"fold=linear\"; defaults: -log10 p, edit ratio and log odds linear, fold and counts log). cells counts the cells with a converted read at a site, over every cell. Thresholds are optional: given ones are drawn as dashed lines and the sites they remove drawn apart; none given, every site is kept. Writes only figures to `out`, a new directory: {modality}_hist_{q} and {modality}_{y}_vs_{x}, each .pdf (labels as text) and .png; returns their paths and the PNGs as images. `qc` with `plots` writes the same panels under qc_plots/panels/.",
            schema(
                json!({
                    "input": dir(),
                    "out": { "type": "string", "description": "Directory the figures go to; must not exist. Omitted, a new directory in today's results folder" },
                    "options": options(),
                }),
                &["input"],
            ),
        ),
        ToolSpec::new(
            "qc_view",
            "Read back what a `qc` run kept and dropped (`faba qc-view --batch-process`), cutting nothing: `input` is the qc output directory, whose qc.run.json gives the thresholds it applied; its kept sites and its dropped ones (with their reason) are read together, each site's cells counted as qc counted them (off the input it read, when that is still there; else the dropped sites read 0 cells, as it warns). Any other faba directory opens too, with no thresholds: every site kept. Returns the table `qc`'s `dry_run` returns, at the thresholds applied: `table`, what the cell QC kept per batch (off its reports) and what each threshold removes per modality (`removes` every site failing it, `only_this` those failing it alone), and `pairs`, the pairs Auto picks to show two at a time. With `plots`, the figures `qc` draws with `plots` go to `out`/qc_plots/ (the site panels to qc_plots/panels/, chosen by `hist`, `pairs`, `marginals`, `style`, `marginal_size`, `p_axis`, `front`, `key`, `scale` and `bins` as `site_plot`'s are) and come back as images, as `qc`'s do; with `html`, the page to share goes to that file. Neither may exist: an existing file or directory is refused, never replaced: ask the user for another path.",
            schema(qc_view_properties(&dir(), &options()), &["input"]),
        ),
        ToolSpec::new(
            "cluster",
            "Cluster cells by their site reads (`faba cluster`). `data` is a faba directory or site matrices. The model is a blockmodel of cell modules and gene modules, a methylation shift per pair and each cell's own level, fitted by variational inference from k-means starts (`chains` restarts, the best by ELBO kept; `svi_epochs`, `batch_size`, `epochs`, `tol`). The MCMC fits of faba before 0.19 are gone: `model` link or crp, `sweeps`, `burn_in`, `shift`, `profile_strength` and `holdout` are refused. Tables go to `out`, a new directory, not back into `data`. With `modality` count, cluster cells by their gene counts instead (count matrices), fitted alike: its verdict is whether the best restart's ELBO settled. With `within`, the run record of such a count fit (count_blockmodel.run.json), fit the sites within each of its cell modules, one fit a module, into `out`/within/{m}/{group}/, with `out`/{m}_within.run.json naming them all and {m}_within_cells.parquet giving each cell its nested label (C3.2); `data` may then be left out, and the matrices are those beside the count fit's inputs. `within` may instead be a cell table (.parquet, such as lupin's {out}.cell_annotation.parquet, or a .tsv/.csv of cell and group, header optional): its groups, taken as given, are fitted as the count fit's cell modules are, each named by its value made a file stem, into the same files, and cluster_view's `group` opens one by that name; `within_column` names the group column (default label, else cluster) and `within_cell_column` the cell column (default cell, barcode, cell_barcode or obs_names); a cell {barcode}@{batch} is that barcode of that batch less a matrix kind (b1_count is b1), its batch paired as a count fit's is (`batch_map` too), and a bare barcode is looked for in each site batch nothing pairs with; `data` defaults to the table's directory; a table whose cells match none of the sites' is refused, naming a cell of each. The sites' batches pair with the count fit's by name (the same, or one named with its folder: `b1` and `runs/a/b1`); when that is ambiguous the run is refused, listing each choice: ask the user which is meant, then run again with `batch_map`. While it fits, a progress line every 15 s says the restart and epoch, the ELBO per read (or count) and its last move against `tol`, the share of steps done and about how long is left (an upper bound: a restart that settles ends early); with `background`, job_status shows it as a `progress` stage (steps done of all, seconds taken and left) and the latest line in `last_lines`.",
            schema(
                json!({
                    "data": { "type": "array", "items": { "type": "string" }, "description": "A faba directory, or site (or, with `modality` count, count) matrices; required unless `within` is given" },
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
            "View a `cluster` fit (`faba cluster-view --batch-process --json`), level by level as the browser does; every level comes back as JSON, an object a level, its rows objects with named fields (a level argument's own text says it prints, as the command's help does: here it is that JSON). With no level named, the tool returns the `overview` (a within record opened whole, with no `group`: its listing instead): the verdict (`verdict` ELBO ok or check, `said`), how many gene and cell modules, each cell module's size (`cell_sizes`) and each gene module's genes and top genes by weight (`gene_sizes`); GO is never read for it. `grid` gives the gene modules against the cell modules (level 1: `modules`, a row a gene module, its `values` a value and reads (cells) in each cell module); `module` one gene module's genes and sites with their values in each cell module (level 2); `gene` one gene's sites in each cell module over every gene module (levels 3 and 4); `cells` the cell modules' sizes; `go` the GO terms (or `gmt` gene sets) over-represented in each gene module, against the fit's own annotated genes (files default to lupin's cache), read only when asked: per module its `terms` (id, name, namespace, k, n, K, N, odds, p, q, genes) or `why_none`, with the `go_q` a term needs and the `min_weight` a gene needs; `batches` how each cell module spreads over the batches (shared, one condition's, or batch-specific: uneven over a condition's replicates), the conditions read off the batch names less a replicate token or given by `batch_groups`; `trace` each restart's ELBO by epoch, the verdict and its advice; `settings` the settings the fit ran with. Several may be given. Values are % converted for a site fit, UMIs per cell for a count fit (`unit`). With `figure` or `plots`, figures instead, each returned as an image with its PDF and PNG paths (`figures`): `plots` several in one call (modules, trace, batches; `plot` one), each `width` by `height` inches, its text at `font_size` points, `panel_only` the plot alone (no title, key or notes), `palette` gray in grays, `xlab` and `ylab` the trace's axis titles (the modules grid and the batches have none and leave them unused); the modules grid coloured by `look`, with `size_bars` and `read_bars` on its margins (sides: left for the gene modules, top or bottom for the cell modules; `marginal_size` their depth), `go_labels` (up to N of each gene module's GO terms on its right, a line each with its q after it, by p and none mostly within a term already shown, read as `go` reads them; a module with none left blank unless `go_none`), `key` (a small key kept under `panel_only`), `values` (each block's value in it) and `block_bars` (each block's reads as the bars look draws them, over the look's own block). `data` is the fit's run record or the directory holding its tables; this never fits, and the record of an MCMC fit of faba before 0.19 is refused. A within record ({m}_within.run.json, from `cluster` with `within`) gives the count fit's modules, each group's fit (cells, status; why one was skipped) and the modules within each fitted group; with `group`, every argument applies to the fit within that cell module, and without it the levels asked for are the count fit's. With `html`, write one self-contained page to a new file instead; the page shows every level, so no level argument is taken beside `html` (leave them out, or drop `html` to get them as JSON); a within record without `group` writes a page a fit, `html` for the count fit and {stem}.{group}.html beside it for each group, linked to each other.",
            schema(
                json!({
                    "data": { "type": "string", "description": "A `cluster` fit: its run record `{m}_blockmodel.run.json`, the output directory holding its block tables, or a within record `{m}_within.run.json`" },
                    "html": { "type": "string", "description": "Write the page here; must not exist" },
                    "figure": { "type": "string", "description": "Write a figure as {figure}.pdf and {figure}.png, or with several `plots` each as {figure}.{plot}.pdf and .png (none may exist), and return each as an image. Omitted with `plots`, in today's results folder" },
                    "plot": { "type": "string", "enum": ["modules", "trace", "batches"], "description": "Which figure: the modules grid (gene modules against cell modules; default), the fit's trace, or the cell modules against the batches" },
                    "plots": { "type": "array", "items": { "type": "string", "enum": ["modules", "trace", "batches"] }, "minItems": 1, "uniqueItems": true, "description": "Several figures in one call, as `plot` names them, each drawn at `width` by `height` in the look asked for and returned as an image" },
                    "look": { "type": "string", "enum": ["shift", "percent", "z", "hinton", "bars"], "description": "What colours the modules grid's blocks: the fitted shift (a count fit's log rate; default), the methylated % (UMIs per cell), that standardized across the cell modules, that as a Hinton diagram (a square per block, its area by |z|), or the block's reads as the browser's mirror bars (converted up, unconverted down, one scale; rows by shift range, each row's range on its right)" },
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
                    "figure": { "type": "boolean", "description": "Also draw the figure beside the table (`{out stem}.pdf` and `.png`) and return it as an image; default true" },
                    "out": file_out(), "options": options(),
                }),
                &["data"],
            ),
        ),
    ];
    specs.extend(view_tools());
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

/// The tools that put a view in front of the user and read back what
/// they arranged in it (see [`crate::tui::live`]).
fn view_tools() -> Vec<ToolSpec> {
    let session = || json!({ "type": "string", "description": "The view's session, as view_open returned it (or its label, or the end of its id)" });
    let views: Vec<&str> = super::tools::views::VIEWS.iter().map(|v| v.0).collect();
    vec![
        ToolSpec::new(
            "view_open",
            "Use when the user wants to arrange a figure by eye (thresholds, scales, the panel, the look) rather than describe it in words. Starts the view as a detached session and returns at once with `session` and `attach`, the exact command for the user to run in their own terminal, of the form `faba attach faba-qc-view-1` (faba's full path; `LEGUME_TUI_DIR=...` in front when sessions live elsewhere). Give them that line verbatim. They arrange the view and tell you (Ctrl+Alt+D detaches and keeps it; q ends it, and then nothing can be read); then view_state reads what they chose and view_export draws the figure. qc_view and cluster_view answer view_state and view_export; pileup, metagene and pwm open for looking only. `args` are the arguments of the tool `view` names (`input` for qc_view, `data` for the others), without those that write files instead (`html`, `plots`, `out`, `figure`, `screen`). `open_window` (macOS) also tries a Terminal window for a few seconds; if macOS blocks it, `window` is false, `window_note` says why, and the session is still there to attach to.",
            schema(
                json!({
                    "view": { "type": "string", "enum": views, "description": "Which view" },
                    "args": { "type": "object", "description": "The arguments of the tool `view` names, as that tool's schema gives them" },
                    "session": { "type": "string", "description": "A label for the session (letters, digits, `-`, `_`, `.`), unique among the live ones; its id is then faba-{screen}-{label}. Default: the next number" },
                    "open_window": { "type": "boolean", "description": "Also open a Terminal window attached to it (macOS only); default false: give the user `attach`" },
                }),
                &["view"],
            ),
        ),
        ToolSpec::new(
            "view_list",
            "faba's live view sessions: each one's `session`, label, whether a terminal is attached, its command and status, whether it `answers` view_state and view_export, and `attach`, the command that shows it.",
            schema(json!({}), &[]),
        )
        .read_only(),
        ToolSpec::new(
            "view_state",
            "Read what the user arranged in a live view (from view_open), as JSON; with `set`, change it first (keys left out stay; the user sees it change) and return the state after. qc_view: `modality`, `panel`, `focus` (the knob), `thresholds` (null: off), `x_scale` per knob, `y_scale` per quantity, `marginal_scale`, `count_scale`, `two_d`, `pairs`, `pair_at`, `front`, `p_axis`, `key`, `key_at`, `titles`, `weight`, `scales`, `bins`, `style`, `marginals`, `marginal_size`; read only: `saves` (what `s` saves now), `gene`, `flags`. cluster_view: `preset`, `cell_order` (read only), `module`, `side` (sites, go), `go_terms`, `pins`. `command` is the view's command line.",
            schema(
                json!({
                    "session": session(),
                    "set": { "type": "object", "description": "A state, or part of one, as view_state returns it, to arrange the view by" },
                }),
                &["session"],
            ),
        ),
        ToolSpec::new(
            "view_export",
            "Draw the figure the live view's own `s` key saves, from what is on screen now, to `out` (formats default pdf and png); `width`, `height` (inches), `font_size`, `panel_only`, `palette` override the view's own. Writes the state beside it as {out}.view.json and returns `replay`, the batch command that draws the same figure again without the view, plus the paths, the state and the image. Never replaces a file.",
            schema(
                json!({
                    "session": session(),
                    "out": { "type": "string", "description": "The figure's prefix: {out}.pdf, {out}.png, {out}.svg as `formats` ask, and {out}.view.json; none may exist. Omitted, in today's results folder" },
                    "formats": { "type": "array", "items": { "type": "string", "enum": ["svg", "pdf", "png"] }, "minItems": 1, "uniqueItems": true, "description": "What to write (default pdf and png)" },
                    "width": { "type": "number", "minimum": crate::figure::size::MIN.0, "description": "The figure's width in inches; omitted, the view's own" },
                    "height": { "type": "number", "minimum": crate::figure::size::MIN.1, "description": "The figure's height in inches; omitted, the view's own" },
                    "font_size": { "type": "number", "exclusiveMinimum": 0, "description": "Every label's size in points; omitted, the view's (fit to the figure)" },
                    "panel_only": { "type": "boolean", "description": "The plot, its axes and their titles only" },
                    "palette": { "type": "string", "enum": ["colour", "gray"], "description": "Kept and removed sites in colour or in grays" },
                }),
                &["session"],
            ),
        ),
    ]
}

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
        &["cluster_view"],
        "json",
        "the tool adds it itself: every level comes back as JSON",
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
}

/// The flags offered under another key than their own.
const RENAMED: &[Renamed] = &[
    // `height` is the figure's.
    Renamed {
        tools: &["pileup"],
        long: "height",
        key: "ascii_height",
    },
    Renamed {
        tools: &SIZED,
        long: "fig-width",
        key: "width",
    },
    Renamed {
        tools: &SIZED,
        long: "fig-height",
        key: "height",
    },
];

/// What the arguments made from flags of the tools named say beyond each
/// flag's type and help: a bound or the values it takes, which its value
/// parser does not tell, or wording of the tool's own arguments.
struct Patched {
    tools: &'static [&'static str],
    keys: &'static [&'static str],
    /// Merged over a tool's property.
    patch: fn(&str) -> Value,
}

/// The arguments made from flags that say more than their flags.
const PATCHED: &[Patched] = &[
    Patched {
        tools: &SIZED,
        keys: &["width"],
        patch: |tool| {
            json!({ "minimum": least_size(tool).0,
                "description": format!("The figure's width in inches, laid out for it, its text kept readable (e.g. 3.5 for a column, 7 for a page); alone or with `height`, the other the figure's own{}", floors(tool)) })
        },
    },
    Patched {
        tools: &SIZED,
        keys: &["height"],
        patch: |tool| {
            json!({ "minimum": least_size(tool).1,
                "description": format!("The figure's height in inches; alone or with `width`, its text fit to it. A qc figure with its table needs 10 inches across unless `panel_only` (it says so){}", floors(tool)) })
        },
    },
    Patched {
        tools: &["depth"],
        keys: &["resolution_kb"],
        patch: |_| json!({ "exclusiveMinimum": 0 }),
    },
    Patched {
        tools: &["qc"],
        keys: &["qc_mads"],
        patch: |_| json!({ "exclusiveMinimum": 0 }),
    },
    Patched {
        tools: &["cluster"],
        keys: &["cell_roles", "gene_roles", "chains", "batch_size", "epochs"],
        patch: |_| json!({ "minimum": 1 }),
    },
    Patched {
        tools: &["cluster"],
        keys: &["tol"],
        patch: |_| json!({ "minimum": 0 }),
    },
    Patched {
        tools: &["cluster", "cluster_view"],
        keys: &["modality"],
        patch: |_| json!({ "enum": ["m6a", "atoi", "count"] }),
    },
    Patched {
        tools: &["pileup", "metagene", "pwm", "site_plot"],
        keys: &["modality"],
        patch: |_| json!({ "enum": ["m6a", "atoi"] }),
    },
    Patched {
        tools: &["cluster_view"],
        keys: &["go_min_weight"],
        patch: |_| json!({ "minimum": 0, "maximum": 1 }),
    },
    Patched {
        tools: &["cluster_view"],
        keys: &["go_q"],
        patch: |_| json!({ "exclusiveMinimum": 0, "maximum": 1 }),
    },
    Patched {
        tools: &["cluster_view"],
        keys: &["marginal_size"],
        patch: |_| json!({ "exclusiveMinimum": 0, "maximum": 1 }),
    },
    Patched {
        tools: &["cluster_view"],
        keys: &["go_labels"],
        patch: |_| json!({ "minimum": 1, "maximum": 10 }),
    },
    Patched {
        tools: &["metagene"],
        keys: &["bins"],
        patch: |_| json!({ "minimum": 1 }),
    },
    Patched {
        tools: &["pwm"],
        keys: &["window"],
        patch: |_| json!({ "minimum": 0 }),
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
            let path = take == Take::Value(Ty::Text)
                && matches!(
                    a.get_value_hint(),
                    ValueHint::FilePath | ValueHint::DirPath | ValueHint::AnyPath
                );
            let mut schema = Generated::property(&take, many, &help_of(a, &take));
            let patches = (PATCHED.iter())
                .filter(|p| p.tools.contains(&tool) && p.keys.contains(&key.as_str()));
            for p in patches {
                if let (Value::Object(m), Value::Object(more)) = (&mut schema, (p.patch)(tool)) {
                    m.extend(more);
                }
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

/// `qc`'s arguments: the directory in, the directory out, and the dry run,
/// figures and page it plans; the site cuts are made from their flags.
fn qc_properties(dir: &Value, out: &Value, options: &Value) -> Value {
    json!({
        "input": dir,
        "out": out,
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
        "figure": { "type": "string", "description": "Also write the panel the view's `s` saves, as `view_state` arranges the view (a state view_export kept, {prefix}.view.json), to {figure}.pdf and .png (`format`: svg, pdf, png), returned as an image: what view_export drew, without the view. None may exist" },
        "options": options,
    })
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

/// The figures a [`SIZED`] tool draws.
fn figures(tool: &str) -> &'static [Figure] {
    match tool {
        "cluster_view" => &[Figure::Modules, Figure::Trace, Figure::Batches],
        "pileup" => &[Figure::Pileup],
        "metagene" => &[Figure::Metagene],
        "pwm" => &[Figure::Logo],
        "qc" | "qc_view" => &[Figure::Qc, Figure::SitePanels],
        "site_plot" => &[Figure::SitePanels],
        _ => &[],
    }
}

/// The least `width` and `height` a [`SIZED`] tool takes, in inches: the
/// least of its figures' (see [`Figure::least`]).
fn least_size(tool: &str) -> (f64, f64) {
    (figures(tool).iter().map(|f| f.least()))
        .reduce(|a, b| (a.0.min(b.0), a.1.min(b.1)))
        .unwrap_or(crate::figure::size::MIN)
}

/// What a tool's size says of its figures that need more than
/// [`least_size`]: each by name, with its own least.
fn floors(tool: &str) -> String {
    let least = least_size(tool);
    let more: Vec<String> = (figures(tool).iter())
        .filter(|f| f.least() != least)
        .map(|f| {
            format!(
                "{} at least {} by {} inches",
                f.name(),
                f.least().0,
                f.least().1
            )
        })
        .collect();
    match more.is_empty() {
        true => String::new(),
        false => format!("; {}", more.join(", ")),
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

    use super::{
        excluded, generated, least_size, properties, tools, Figure, EXCLUDED, PATCHED, SIZED,
    };
    use crate::mcp::tools::sites::pileup_argv;
    use serde_json::json;
    /// Each tool's arguments, a line each: its key, its JSON type and the
    /// bounds and choices it states.
    fn shapes() -> String {
        let mut lines = Vec::new();
        for spec in tools() {
            let props = spec.input_schema["properties"].as_object().unwrap();
            for (key, p) in props {
                let mut line = format!("{} {key}", spec.name);
                let mut part = |name: &str, v: &serde_json::Value| {
                    if !v.is_null() {
                        line.push_str(&format!(" {name}={v}"));
                    }
                };
                part("type", &p["type"]);
                part("items", &p["items"]["type"]);
                part("anyOf", &json!(p["anyOf"].as_array().map(|a| a.len())));
                let mut choices: Vec<_> = (p["enum"].as_array().into_iter().flatten())
                    .chain(p["items"]["enum"].as_array().into_iter().flatten())
                    .filter_map(|c| c.as_str())
                    .collect();
                choices.sort_unstable();
                part("enum", &json!(Some(choices).filter(|c| !c.is_empty())));
                // A bound as a number: 1 and 1.0 are one bound to an agent.
                for bound in ["minimum", "exclusiveMinimum", "maximum", "minItems"] {
                    part(bound, &json!(p[bound].as_f64()));
                }
                lines.push(line);
            }
        }
        lines.sort_unstable();
        lines.join("\n") + "\n"
    }

    /// Every tool's arguments keep their keys, JSON types, bounds and
    /// choices, as `properties.snap` lists them: an argument moved between
    /// the hand-written ones and those made from flags is the same to an
    /// agent.
    #[test]
    fn every_tools_arguments_keep_their_keys_and_types() {
        let (got, want) = (shapes(), include_str!("properties.snap"));
        let diff: Vec<String> = (want.lines().filter(|l| !got.lines().any(|g| g == *l)))
            .map(|l| format!("- {l}"))
            .chain(
                (got.lines().filter(|l| !want.lines().any(|w| w == *l))).map(|l| format!("+ {l}")),
            )
            .collect();
        assert!(diff.is_empty(), "arguments changed:\n{}", diff.join("\n"));
    }

    /// A hand-written argument under a flag's key is one the tool does
    /// more with than pass on: a file it checks or notes, a figure or page
    /// it plans, a value it joins or checks against another. Every other
    /// flag's argument is made from the flag, its help and [`PATCHED`].
    #[test]
    fn a_hand_written_argument_for_a_flag_is_one_the_tool_composes() {
        let composed: &[(&str, &[&str])] = &[
            // Noted as inputs, or read before the command runs.
            ("count", &["gff"]),
            ("depth", &["valid_cells"]),
            ("atoi", &["gff", "genome"]),
            ("dartseq", &["gff", "genome"]),
            ("apa", &["gff", "utr_bed"]),
            ("snp", &["gff", "genome", "known_snps"]),
            ("run", &["gff", "genome", "known_snps"]),
            ("metagene", &["gff"]),
            ("pwm", &["genome"]),
            // A run without the cut, the figures and the page it plans.
            ("qc", &["dry_run", "plots", "html"]),
            ("qc_view", &["plots", "html", "figure"]),
            // The figures, each planned and checked new.
            ("cluster", &["batch_map"]),
            ("cluster_view", &["html", "figure", "plot", "look"]),
            ("pileup", &["genes", "regions", "gtf", "bam", "out"]),
            ("metagene", &["dist_measures", "figure"]),
            ("pwm", &["figure"]),
        ];
        let made = |tool: &str, key: &str| generated(tool).iter().any(|g| g.key == key);
        let mut extra = Vec::new();
        for (tool, sub) in commands() {
            let props = properties(&tool).unwrap();
            for a in sub.get_arguments().filter(|a| !a.is_hide_set()) {
                let Some(long) = a.get_long() else { continue };
                let key = long.replace('-', "_");
                let kept =
                    (composed.iter()).any(|(t, keys)| *t == tool && keys.contains(&key.as_str()));
                if props.contains_key(&key) && !made(&tool, &key) && !kept {
                    extra.push(format!("{tool} {key}"));
                }
            }
        }
        assert!(
            extra.is_empty(),
            "written by hand, but only passed on:\n{}",
            extra.join("\n")
        );
    }

    /// Each figure tool's `width` and `height` state the least of its
    /// figures' floors, and its command refuses a side under a figure's
    /// floor and takes one at it: one floor a figure (see
    /// [`Figure::least`]). cluster_view's least is its modules grid's; its
    /// description names the trace's and the batches', which are more.
    #[test]
    fn the_schema_states_the_floor_each_command_keeps() {
        use crate::figure::SizeArgs;
        let specs = tools();
        for tool in SIZED {
            let spec = specs.iter().find(|s| s.name == tool).unwrap();
            let props = &spec.input_schema["properties"];
            let (w, h) = (
                props["width"]["minimum"].as_f64(),
                props["height"]["minimum"].as_f64(),
            );
            let drawn = super::figures(tool);
            assert!(!drawn.is_empty(), "{tool} draws no figure");
            let least = |side: fn((f64, f64)) -> f64| {
                drawn
                    .iter()
                    .map(|f| side(f.least()))
                    .fold(f64::INFINITY, f64::min)
            };
            assert_eq!(
                (w, h),
                (Some(least(|l| l.0)), Some(least(|l| l.1))),
                "{tool}"
            );
            for f in drawn {
                let (fw, fh) = f.least();
                for under in [
                    SizeArgs::new(Some(fw - 0.01), None),
                    SizeArgs::new(None, Some(fh - 0.01)),
                ] {
                    let e = f.page(&under).unwrap_err().to_string();
                    assert!(e.contains("at least"), "{tool} {f:?}: {e}");
                }
                assert!(
                    f.page(&SizeArgs::new(Some(fw), Some(fh))).is_ok(),
                    "{tool} {f:?}"
                );
                if f.least() != (w.unwrap(), h.unwrap()) {
                    let said = props["width"]["description"].as_str().unwrap();
                    let floor = format!("{} at least {fw} by {fh} inches", f.name());
                    assert!(said.contains(&floor), "{tool}: {said}");
                }
            }
        }
        // Through the commands, which check the size before reading
        // anything: a side under a figure's floor is refused saying so, one
        // at it is not refused for its size.
        let run = |line: Vec<String>| {
            let m = crate::faba_command().try_get_matches_from(line).unwrap();
            match crate::parse(&m).unwrap().commands {
                crate::Commands::Pileup(ref args) => crate::run_pileup(args),
                crate::Commands::ClusterView(ref args) => crate::run_modules(args),
                _ => unreachable!(),
            }
        };
        let view = ["cluster-view", "/no/such", "--figure", "f", "--plot"];
        let lines = [
            ("pileup", Figure::Pileup, vec!["pileup", "/no/such.zarr"]),
            (
                "cluster_view",
                Figure::Modules,
                [&view[..], &["modules"]].concat(),
            ),
            (
                "cluster_view",
                Figure::Trace,
                [&view[..], &["trace"]].concat(),
            ),
            (
                "cluster_view",
                Figure::Batches,
                [&view[..], &["batches"]].concat(),
            ),
        ];
        let sized = |line: &[&str], w: f64, h: f64| {
            let size = [
                "--fig-width".to_string(),
                w.to_string(),
                "--fig-height".into(),
                h.to_string(),
            ];
            let all = std::iter::once("faba")
                .chain(line.iter().copied())
                .map(String::from);
            let e = run(all.chain(size).collect()).unwrap_err().to_string();
            !e.contains("inches")
        };
        for tool in ["pileup", "cluster_view"] {
            let (w, h) = least_size(tool);
            let of = || lines.iter().filter(|l| l.0 == tool);
            // The schema's least is one a figure of the tool takes, and
            // none takes less.
            assert!(
                of().any(|(_, _, line)| sized(line, w, h)),
                "{tool}: {w} by {h}"
            );
            for (_, f, line) in of() {
                assert!(
                    !sized(line, w - 0.01, h) && !sized(line, w, h - 0.01),
                    "{tool} {f:?}"
                );
                let (fw, fh) = f.least();
                assert!(sized(line, fw, fh), "{f:?}");
                assert!(
                    !sized(line, fw - 0.01, fh) && !sized(line, fw, fh - 0.01),
                    "{f:?}"
                );
            }
        }
    }

    /// Every patch is of an argument its tools make from a flag, and is in
    /// the schema they offer; a choice made from a flag is one its command
    /// takes.
    #[test]
    fn every_patch_is_of_an_argument_made_from_a_flag() {
        for p in PATCHED {
            for (tool, key) in (p.tools.iter()).flat_map(|t| p.keys.iter().map(move |k| (*t, *k))) {
                assert!(
                    generated(tool).iter().any(|g| g.key == key),
                    "{tool}: `{key}` is not made from a flag"
                );
                let offered = &properties(tool).unwrap()[key];
                for (k, v) in (p.patch)(tool).as_object().unwrap() {
                    assert_eq!(&offered[k], v, "{tool} {key} {k}");
                }
            }
        }
        let apa = |method| typed_flags("apa", &Arguments(json!({ "method": method })));
        assert_eq!(apa("simple").unwrap(), ["--method", "simple"]);
        let e = apa("fast").unwrap_err().to_string();
        assert!(e.contains("one of simple, mixture"), "{e}");
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
        // The schema states the least a figure takes: pileup's, whose
        // margins do not shrink, more; cluster_view's, its modules grid's.
        let specs = tools();
        let least = |tool: &str, key: &str| {
            let t = specs.iter().find(|t| t.name == tool).unwrap();
            t.input_schema["properties"][key]["minimum"].clone()
        };
        let (min, fixed) = (crate::figure::size::MIN, Figure::Pileup.least());
        for tool in ["pwm", "metagene", "qc", "site_plot", "cluster_view"] {
            assert_eq!(least(tool, "width"), json!(min.0), "{tool}");
            assert_eq!(least(tool, "height"), json!(min.1), "{tool}");
        }
        assert_eq!(least("pileup", "width"), json!(fixed.0));
        assert_eq!(least("pileup", "height"), json!(fixed.1));
    }
}
