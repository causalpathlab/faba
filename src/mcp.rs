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
first. qc_view reads back a qc output: the same counts at the thresholds it \
applied, its dropped sites with its kept, and with `plots` or `html` its \
figures and page.
4. pileup, metagene and pwm read a directory's sites. metagene and pwm \
return their table and its figure as an image. site_plot draws a \
directory's site statistics as paper panels (histograms, pairs, \
marginals), thresholds optional, and returns them as images. pileup returns a binned table \
for a gene or region (and, with `figure`, a Miami SVG, as an image too). \
cluster_view shows a cluster fit level by level, from its run record \
({m}_blockmodel.run.json) or the directory holding its tables, each level \
as JSON rows: by default an overview (the ELBO verdict, the modules' sizes, \
each gene module's top genes), then the modules grid, a gene module, a gene, \
the cell modules, the batches, the trace, the settings, and GO terms only \
with `go`; with `plots`, several figures as images.
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
Paths resolve against the server's working directory. Every flag of a \
tool's command is an argument of the tool, typed, named after the flag \
(`--min-coverage` is `min_coverage`), but those naming the output or \
opening a view; `options` remains for anything else, as `faba <command> \
--help` lists it, and a flag given both ways is refused. Results \
give the paths written and what the run recorded (`results`: a fit's ELBO \
verdict, its cluster sizes). \
docs gives this server's `version`. A call is refused when the faba it \
would run is another version than the server's, or is the server's own \
file rebuilt since it started: the server then needs a restart or a \
reconnect.";

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
            _ => bail!("no tool {name}"),
        };
        Ok(vec![Content::Json(v)])
    }
}
