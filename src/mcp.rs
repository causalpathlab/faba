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

mod args;
mod run;
mod schema;
mod tools;

use crate::mcp::args::{known_keys, with_typed, Arguments};
use crate::mcp::run::Binary;
use crate::mcp::schema::tools;
use crate::mcp::tools::cluster::{cluster, cluster_view};
use crate::mcp::tools::docs::docs;
use crate::mcp::tools::producers::{apa, atoi, count, dartseq, depth, run, snp};
use crate::mcp::tools::qc::{qc, qc_view, site_plot};
use crate::mcp::tools::sites::{metagene, pileup, pwm};
use crate::mcp::tools::views::{view_export, view_list, view_open, view_state};
use anyhow::{bail, Result};
use clap::Args;
use legume_mcp::{Call, Content, Handler, ServerInfo, ToolSpec};
use serde_json::Value;

#[derive(Args, Debug)]
pub struct McpArgs {}

const INSTRUCTIONS: &str = "\
faba turns coordinate-sorted BAM alignments into per-cell feature matrices: \
gene counts, read depth, A-to-I editing, DART-seq m6A, alternative \
polyadenylation and SNP genotypes. Row names are \
{gene_key}/{modality}/{detail}.

The steps, in order:
1. docs lists the method write-ups; with `query`, it prints the sections \
named (a topic such as profiling, a command such as qc, or a section number).
2. Producers write every candidate, with no p-value cutoff, each into a new \
directory: count, atoi, dartseq, apa, depth and snp one modality each; run \
the pipeline (SNP, then count, A-to-I, m6A, APA) at once. dartseq needs \
control BAMs; without them, run skips m6A.
3. qc cuts cells, features and sites into a new directory; the input is never \
changed. Try thresholds first with its `dry_run`, which writes nothing and \
returns what the cell QC keeps and what each threshold removes. qc_view reads \
back a cut already made: the same counts at the thresholds it applied.
4. cluster fits cells by their own reads (a blockmodel of cell and gene \
modules, by variational inference), in two rounds: first with `modality` \
count (gene counts; it writes count_blockmodel.run.json), then with `within` \
naming that record and `modality` m6a or atoi, a fit within each cell \
module ({m}_within.run.json). cluster_view reads a fit back: JSON levels \
(an overview by default), figures, the browser's screen or a page.
5. pileup (one gene's or region's sites), metagene (where sites sit over all \
genes) and pwm (bases around the sites) read a directory's sites: a table, \
and a figure as an image.

Which figure tool:
- site_plot: paper panels (histograms, pairs, marginals) of a directory's \
site statistics at any thresholds, no cut.
- qc with `plots`: the same panels and the qc view's figures, at a cut being \
made (or a `dry_run`); qc_view with `plots`: those of a cut already made.
- qc_view with `figure`: the one panel a saved view state arranges, with no \
live view.
- cluster_view with `figure` or `plots`: paper figures of a fit (the modules \
grid, the trace, the batches); with `screen`: the browser's screen as drawn \
(the grid, the sites or GO beside it, the pinned genes), and with \
`each_fit`, `each_group` or `each_module` one per fit, group or module.
- view_export: what the user arranged by eye in a live view.

When the user wants to arrange a figure by eye: view_open (qc_view or \
cluster_view; pileup, metagene and pwm open for looking only) returns an \
`attach` line such as `faba attach faba-qc-view-1`; give it to the user to \
run in their terminal. They detach with Ctrl+] when done (q ends the view, \
and then nothing can be read). view_state reads what they chose (and with \
`set` changes it), view_export draws that figure and returns `replay`, the \
batch command that draws it again. view_list lists the live views.

The producers (count, depth, atoi, dartseq, apa, snp, run, qc, cluster) take \
`background`: on large BAMs set it, and the call returns a `job` id at once; \
follow it with job_wait or job_status.

Every write goes to `out`, which must not exist: no tool replaces a file or \
a directory (when one exists, ask the user for another path). Leave `out` \
out and the output goes in today's results folder, named after the input \
and the tool. Paths resolve against the server's working directory. Every \
flag of a tool's command is an argument of the tool, typed, named after the \
flag (`--min-coverage` is `min_coverage`), but those naming the output or \
opening a view; `options` remains for anything else, as `faba <command> \
--help` lists it, and a flag given both ways is refused. A write returns \
the paths `written` and what the run recorded (`results`: a fit's ELBO \
verdict, its cluster sizes). docs gives this server's `version`. A call is \
refused when the faba it would run is another version than the server's, \
or is the server's own file rebuilt since it started: the server then \
needs a restart or a reconnect.";

pub fn run_mcp(_args: &McpArgs) -> Result<()> {
    let info = ServerInfo {
        name: "faba".into(),
        version: env!("CARGO_PKG_VERSION").into(),
        instructions: Some(INSTRUCTIONS.into()),
    };
    legume_mcp::serve(
        info,
        Faba {
            binary: Binary::this_server(),
        },
    )
}

struct Faba {
    binary: Binary,
}

impl Handler for Faba {
    fn tools(&self) -> Vec<ToolSpec> {
        tools()
    }

    fn call(&mut self, name: &str, args: Value, cx: &Call) -> Result<Vec<Content>> {
        // `modules` is `cluster_view`'s older name.
        let name = match name {
            "modules" => "cluster_view",
            name => name,
        };
        // Before anything is resolved (`out` may make today's folder).
        known_keys(name, &args)?;
        // Every tool runs faba: refused when it is not this server's.
        self.binary.check()?;
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
            "qc_view" => return qc_view(&a, cx),
            "site_plot" => return site_plot(&a, cx),
            "cluster" => cluster(&a, cx)?,
            "cluster_view" => return cluster_view(&a, cx),
            "pileup" => return pileup(&a, cx),
            "metagene" => return metagene(&a, cx),
            "pwm" => return pwm(&a, cx),
            "view_open" => return view_open(&a, cx),
            "view_list" => return view_list(&a, cx),
            "view_state" => return view_state(&a, cx),
            "view_export" => return view_export(&a, cx),
            _ => bail!("no tool {name}"),
        };
        Ok(vec![Content::Json(v)])
    }
}
