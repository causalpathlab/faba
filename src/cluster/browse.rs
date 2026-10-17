//! `faba modules`: the blockmodel's gene modules and cell modules, browsed.
//!
//! Reads the tables `faba cluster --model blocks` writes (running it first
//! when they are missing, see [`data`]) and opens a view of four linked
//! levels, each a margin of the next over one table of reads by gene module,
//! site and cell module (see [`view`]): every gene module against every
//! cell module, the sites of the module in focus, the genes pinned side by
//! side, and one gene's full pileup from the site matrices (see [`gene`]).
//! `--html` (or `e`) writes the same as one page (see [`report`]).

use super::{ClusterArgs, ModelKind};
use crate::site_analysis::input_picker::{ask_reference, Reference};
use crate::site_analysis::miami::genemodel::load_gene_models_where;
use clap::Args;
use std::path::Path;

mod data;
mod gene;
mod order;
mod report;
mod view;

pub use data::{gene_models, Count, Modules};
#[cfg(test)]
pub use order::diagonal;
pub use order::standardized;
pub use view::ModulesView;

#[derive(Args, Debug, Clone)]
pub struct ModulesArgs {
    #[command(flatten)]
    cluster: ClusterArgs,

    #[arg(
        long = "gtf",
        value_name = "FILE",
        help = "Gene annotation GTF/GFF for the gene models over the pileups",
        long_help = "Gene annotation GTF/GFF for the gene models over the pileups (exons,\n\
                     introns, strand). Without it, the GFF recorded next to the tables is\n\
                     used; when there is none on this machine, the browser asks for one, or\n\
                     downloads one, as `faba pileup` does."
    )]
    gtf: Option<Box<str>>,

    #[arg(
        long = "refit",
        default_value_t = false,
        help = "Cluster again with these options even when block tables exist",
        long_help = "Cluster again (`faba cluster --model blocks` with these options) even when\n\
                     the output directory has block tables. Without it, existing tables are\n\
                     opened as they are and the clustering options are not used."
    )]
    refit: bool,

    #[arg(
        long = "html",
        value_name = "FILE",
        help = "Write the modules as one self-contained HTML page, to share or post, instead of browsing",
        long_help = "Write the modules as one self-contained HTML page instead of browsing: the\n\
                     browser's levels and keys, its data embedded, nothing loaded from anywhere,\n\
                     so it can be posted (a GitHub page) or sent as it is. Gene models from\n\
                     --gtf or the recorded GFF, when there is one. In the browser, `e` writes it."
    )]
    html: Option<Box<str>>,

    #[arg(
        long = "batch-process",
        visible_alias = "batch-mode",
        default_value_t = false,
        help = "Print the modules instead of browsing them full screen"
    )]
    batch_process: bool,
}

impl ModulesArgs {
    /// Notes which clustering options `matches` were given, for the
    /// command `faba cluster` saves.
    pub fn record_given(&mut self, matches: &clap::ArgMatches) {
        self.cluster.record_given(matches);
    }

    /// The clustering options, with the blockmodel these tables come from.
    pub(super) fn blocks(&self) -> ClusterArgs {
        let mut cluster = self.cluster.clone();
        cluster.model = ModelKind::Blocks;
        cluster
    }
}

/// The modality whose block tables `dir` holds: the one asked for, else
/// m6a, else atoi; none when there are none.
fn tables(dir: &Path, asked: Option<&str>) -> Option<String> {
    let has = |m: &str| dir.join(format!("{m}_block_sites.parquet")).is_file();
    match asked {
        Some(m) => Some(m.to_ascii_lowercase()).filter(|m| has(m)),
        None => crate::qc::layout::SITE_MODALITIES
            .iter()
            .find(|m| has(m))
            .map(|m| m.to_string()),
    }
}

/// The modules as text: per gene module, its genes and each cell module's
/// methylated % and reads.
fn print(m: &Modules) {
    let mut head = format!("{:<8}{:>7}{:>8}", "module", "genes", "range");
    for &c in &m.columns {
        head.push_str(&format!("{:>14}", format!("C{} %/reads", c + 1)));
    }
    println!("{head}");
    for &g in &m.order {
        let mut line = format!(
            "{:<8}{:>7}{:>8.2}",
            format!("M{}", g + 1),
            m.gene_sizes[g],
            m.shift_range(g)
        );
        for r in m.block(g) {
            line.push_str(&format!("{:>14}", format!("{}/{:.0}", r.pct(), r.total)));
        }
        println!("{line}");
    }
}

pub fn run_modules(args: &ModulesArgs) -> anyhow::Result<()> {
    // No inputs named: browse for them, as `faba pileup` does.
    let batch = args.batch_process || args.html.is_some();
    let Some(cluster) = super::with_inputs(&args.cluster, "modules", batch)? else {
        return Ok(());
    };
    let args = &ModulesArgs {
        cluster,
        ..args.clone()
    };
    let dir = super::output_dir(&args.cluster);
    let modality = ensure_tables(args, &dir)?;
    // Shared with the thread that writes the report on `e`.
    let modules = std::sync::Arc::new(Modules::load(&dir, &modality)?);
    let gtf = crate::run_record::explicit_or_recorded(
        args.gtf.as_deref(),
        &dir.to_string_lossy(),
        "gff",
        "gene models",
    );
    if let Some(path) = &args.html {
        return write_html(&modules, gtf.as_deref(), path);
    }
    if crate::tui::batch(args.batch_process, "modules") {
        print(&modules);
        return Ok(());
    }
    browse(args, &dir, &modules, gtf)
}

/// The modality of the block tables in `dir`: those there, unless asked to
/// cluster again, else those `faba cluster --model blocks` writes now.
fn ensure_tables(args: &ModulesArgs, dir: &Path) -> anyhow::Result<String> {
    let asked = args.cluster.modality.as_deref();
    if let Some(m) = tables(dir, asked).filter(|_| !args.refit) {
        log::info!(
            "opening the {m} block tables in {} as they are (--refit clusters again)",
            dir.display()
        );
        return Ok(m);
    }
    log::info!(
        "running `faba cluster --model blocks` in {} first",
        dir.display()
    );
    // The tables of the modality just fitted, not older ones beside them.
    super::run_cluster(&args.blocks())?
        .and_then(|m| tables(dir, Some(&m)))
        .ok_or_else(|| anyhow::anyhow!("`faba cluster` wrote no block tables in {}", dir.display()))
}

/// `--html`: the modules as one page at `path`, gene models from `gtf`.
fn write_html(modules: &Modules, gtf: Option<&str>, path: &str) -> anyhow::Result<()> {
    let all = match gtf {
        Some(gtf) => load_gene_models_where(gtf, |_| true)?,
        None => Vec::new(),
    };
    report::write(modules, &gene_models(modules, &all), Path::new(path))?;
    log::info!("wrote {path}");
    Ok(())
}

/// The browser, until it is quit: Enter on a gene opens its full pileup,
/// and leaving that comes back to the browser as it was.
fn browse(
    args: &ModulesArgs,
    dir: &Path,
    modules: &std::sync::Arc<Modules>,
    gtf: Option<Box<str>>,
) -> anyhow::Result<()> {
    // The matrices are indexed while the browser opens and the annotation
    // is asked for: those of the modality of the tables opened, whatever
    // the default.
    let cluster = super::ClusterArgs {
        modality: Some(modules.modality.as_str().into()),
        ..args.cluster.clone()
    };
    let matrices = gene::Matrices::start(super::inputs(&cluster)?, modules);
    let gtf = match gtf {
        Some(gtf) => Some(gtf),
        None => ask_reference("modules", Reference::Annotation, None)?
            .map(|p| p.to_string_lossy().into()),
    };
    let models = crate::site_analysis::pileup::start_gene_models(gtf.as_deref());
    let title = format!("{} · {}", dir.display(), modules.modality);
    let mut view = ModulesView::new(&title, modules)
        .with_genes(models.clone())
        .with_report(
            dir.join(format!("{}_modules.html", modules.modality)),
            modules.clone(),
        );
    let mut settings = Default::default();
    loop {
        crate::tui::run_view(&mut view)?;
        let Some(key) = view.take_open() else {
            return Ok(());
        };
        // A gene that cannot be read leaves the browser as it was, saying why.
        let at = (key.as_str(), view.first_column());
        match gene::show(&matrices, modules, at, models.as_ref(), &mut settings) {
            Ok(true) => return Ok(()),
            Ok(false) => {}
            Err(e) => view.set_status(format!("{e:#}")),
        }
    }
}

#[cfg(test)]
mod tests;
