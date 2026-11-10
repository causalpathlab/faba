//! `faba cluster-view`: the blockmodel's gene modules and cell modules, browsed.
//!
//! Reads the tables `faba cluster`'s blockmodel writes (see [`data`]), as
//! they are: it never fits, and has none of the options that would (see
//! [`fitting_hint`]). It opens a view of four linked levels, each a
//! margin of the next over the posterior-mean reads by gene
//! role, site and cell role (see [`view`]): every gene module against every
//! cell module, the sites of the module in focus, the genes pinned side by
//! side (each in its main gene module), and one gene's full pileup over
//! every gene role (see [`gene`]).
//! `--html` (or `e`) writes the same as one page (see [`report`]).
//! A within record (`faba cluster --within`) opens as the count fit and
//! the fit within each of its cell modules (see [`within`]).

use super::command::{is_count, COUNT_MODALITY};
use crate::site_analysis::input_picker::{ask_reference, Reference};
use crate::site_analysis::miami::genemodel::load_gene_models_where;
use crate::site_analysis::pileup::{start_gene_models, SharedModels};
use anyhow::Context;
use clap::Args;
use std::path::{Path, PathBuf};

mod data;
mod gene;
mod order;
mod report;
mod text;
mod view;
mod within;

pub(super) use data::read_record;
pub use data::{gene_models, Count, Modules};
#[cfg(test)]
pub use order::diagonal;
pub use order::standardized;
pub use view::ModulesView;

#[derive(Args, Debug, Clone)]
pub struct ModulesArgs {
    #[arg(
        num_args = 0..,
        value_name = "RECORD|DIR",
        value_hint = clap::ValueHint::AnyPath,
        help = "A fit's run record, or the output directory (or site matrices) whose block tables to open; browse for one when left out",
        long_help = "A `faba cluster` fit's run record, {modality}_blockmodel.run.json, or\n\
                     the output directory (or site matrices) whose block tables to open, as\n\
                     they are. A within record, {modality}_within.run.json (`faba cluster\n\
                     --within`), opens the count fit and the fit within each of its cell\n\
                     modules. Left out, a browser asks, as `faba pileup`'s does."
    )]
    data_files: Vec<Box<str>>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi|count",
        help = "Which fit's tables to open from a directory (default: m6a, else atoi); count: the count fit's"
    )]
    modality: Option<Box<str>>,

    #[arg(
        long = "group",
        value_name = "C",
        help = "Of a within record: open the fit within this cell module of the count fit, as C3 or 3",
        long_help = "Of a within record ({modality}_within.run.json, which `faba cluster\n\
                     --within` writes): open the fit within this cell module of the count\n\
                     fit, as C3 or 3. Every other option\n\
                     then applies to that fit. Without it, the record opens as a list of\n\
                     its fits, the count fit first; the levels printed are the count\n\
                     fit's, and --html writes a page a fit, FILE for the count fit and\n\
                     {FILE's stem}.{group}.html beside it, each linking to the others."
    )]
    group: Option<String>,

    #[arg(
        short = 'o',
        long = "output",
        value_name = "DIR",
        value_hint = clap::ValueHint::DirPath,
        help = "The directory the tables are in (default: the input's directory)"
    )]
    output: Option<Box<str>>,

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

    /// The levels printed under --batch-process instead of the modules.
    #[command(flatten)]
    shown: text::Shown,

    /// `--no-detach` and `--session`: how the view runs as a session.
    #[command(flatten)]
    pub session: legume_tui::cli::SessionArgs,
}

impl ModulesArgs {
    /// The session the browser runs in; none under --batch-process or
    /// --html. Matrices named start a new one; with none, live sessions
    /// are offered first.
    pub fn view_session(&self) -> Option<legume_tui::Options> {
        let target = !self.data_files.is_empty();
        let batch = self.batch_process || self.html.is_some();
        (!batch).then(|| self.session.options(target))
    }
}

/// What `faba cluster-view` says when it is given an option it does not
/// have (`viewing`: the command line named it): `error` with one more tip
/// when that option is one of `faba cluster`'s (the model's, the sampler's,
/// the reads' filters), as the viewer opens fits as they are and never
/// fits; any other error as it is.
pub fn fitting_hint(mut error: clap::Error, viewing: bool) -> clap::Error {
    use clap::error::{ContextKind, ContextValue, ErrorKind};
    let Some(ContextValue::String(given)) = error.get(ContextKind::InvalidArg) else {
        return error;
    };
    let flag = given.split('=').next().unwrap_or(given);
    let of_cluster = super::command::definition().get_arguments().any(|a| {
        let longs = a
            .get_long()
            .into_iter()
            .chain(a.get_all_aliases().into_iter().flatten());
        let shorts = a
            .get_short()
            .into_iter()
            .chain(a.get_all_short_aliases().into_iter().flatten());
        longs
            .map(|l| format!("--{l}"))
            .chain(shorts.map(|c| format!("-{c}")))
            .any(|f| f == flag)
    });
    if error.kind() != ErrorKind::UnknownArgument || !viewing || !of_cluster {
        return error;
    }
    let mut tips = match error.get(ContextKind::Suggested) {
        Some(ContextValue::StyledStrs(tips)) => tips.clone(),
        _ => Vec::new(),
    };
    tips.push(FITTING_HINT.into());
    error.insert(ContextKind::Suggested, ContextValue::StyledStrs(tips));
    error
}

/// [`fitting_hint`]'s tip.
const FITTING_HINT: &str = "`faba cluster-view` opens fits as they are and never fits: run \
     `faba cluster` (or `faba cluster --modality count`) with that option, then open its run \
     record here";

/// What the browser offers to open: a fit's run record, or site matrices
/// (or their output directory) whose fit's tables are beside them.
fn is_view_input(name: &str) -> bool {
    crate::run_record::is_run_record(name) || crate::qc::layout::is_site_matrix_name(name)
}

/// A run record: its file name and contents.
pub(crate) type Record = (String, serde_json::Value);

/// A blockmodel fit as its run record names it: the record's directory,
/// its file name and contents, and the fit's modality.
#[derive(Debug)]
struct FitRecord {
    pub dir: PathBuf,
    pub record: Record,
    pub modality: String,
}

/// The fit `inputs` name when one is a `faba cluster` run record (see
/// [`read_fit_record`]); none for other inputs, another command's record
/// among them (a faba output directory's, taken as before). A fit's record
/// with other inputs is refused, as one view shows one fit.
fn fit_record(inputs: &[Box<str>]) -> anyhow::Result<Option<FitRecord>> {
    let mut fits = Vec::new();
    for p in inputs {
        if let Some(fit) = read_fit_record(Path::new(p.as_ref()))? {
            fits.push(fit);
        }
    }
    anyhow::ensure!(
        fits.is_empty() || inputs.len() == 1,
        "name one fit's run record alone, not with other inputs"
    );
    Ok(fits.pop())
}

/// The run record at `path`: its directory, its file name and its
/// contents; none when the file there is not named as one.
fn record_at(path: &Path) -> anyhow::Result<Option<(PathBuf, String, serde_json::Value)>> {
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or_default();
    if !crate::run_record::is_run_record(name) || !path.is_file() {
        return Ok(None);
    }
    let dir = crate::run_record::output_dir_of(&path.to_string_lossy())
        .with_context(|| format!("{} is not a run record", path.display()))?;
    let record = data::read_record(&dir, name).map_err(anyhow::Error::msg)?;
    Ok(Some((dir, name.to_string(), record)))
}

/// The blockmodel fit the file at `path` describes when it is a `faba
/// cluster` run record, from its contents: none for any other file
/// (another command's record among them), refused for a crp fit's record,
/// which has no modules. The modality is the one its `tables` say; an older
/// record's, the one its name begins with, else the one whose tables are
/// beside it. Each table it names is looked for beside it.
fn read_fit_record(path: &Path) -> anyhow::Result<Option<FitRecord>> {
    let Some((dir, name, record)) = record_at(path)? else {
        return Ok(None);
    };
    let name = name.as_str();
    if record["job"].as_str() != Some("cluster") {
        return Ok(None);
    }
    let results = &record["results"];
    let blockmodel = !results["model"].is_null() || !results["blockmodel"].is_null();
    anyhow::ensure!(
        blockmodel || results["posterior"].is_null(),
        "{} is the record of a `faba cluster --model crp` fit, which has no modules: \
         name a blockmodel fit's (--model link or blocks)",
        path.display()
    );
    anyhow::ensure!(blockmodel, "{} records no blockmodel fit", path.display());
    let entry = &results["tables"];
    let modality = entry["modality"]
        .as_str()
        .map(str::to_string)
        .or_else(|| {
            crate::qc::layout::SITE_MODALITIES
                .iter()
                .find(|m| name.starts_with(&format!("{m}_")))
                .map(|m| m.to_string())
        })
        .or_else(|| tables(&dir, None))
        .with_context(|| format!("{name} names no modality, and no tables are beside it"))?;
    let sides = ["cells", "genes"].map(|side| entry[side].as_array());
    let listed = sides
        .into_iter()
        .flatten()
        .flatten()
        .chain(Some(&entry["trace"]))
        .filter_map(|v| v.as_str());
    for table in listed {
        // Older fits wrote fewer tables: those still read approximately.
        if !dir.join(table).is_file() {
            log::warn!("{table}, named in {name}, is not beside it");
        }
    }
    anyhow::ensure!(
        tables(&dir, Some(&modality)).is_some(),
        "no {modality} block tables beside {}",
        path.display()
    );
    Ok(Some(FitRecord {
        dir,
        record: (name.to_string(), record),
        modality,
    }))
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

/// A fit to open: its tables' directory, its modality, and its run record
/// when it was named by one.
pub(crate) type Opened = (PathBuf, String, Option<Record>);

/// What `--group` is told without a within record.
const NO_GROUP: &str = "--group names a fit within a count fit's cell module: name a within \
     record ({modality}_within.run.json), which `faba cluster --within` writes";

pub fn run_modules(args: &ModulesArgs) -> anyhow::Result<()> {
    // No inputs named: browse for them, as `faba pileup` does.
    let batch = args.batch_process || args.html.is_some();
    let offered = (
        "a fit's run record or site matrices",
        is_view_input as fn(&str) -> bool,
    );
    let Some(inputs) = super::with_inputs(&args.data_files, "cluster-view", batch, offered)? else {
        return Ok(());
    };
    let (modality, output) = (args.modality.as_deref(), args.output.as_deref());
    let recorded = |dir: &Path| {
        let near = dir.to_string_lossy();
        crate::run_record::explicit_or_recorded(args.gtf.as_deref(), &near, "gff", "gene models")
    };
    match target(&inputs, modality, output, args.group.as_deref())? {
        Target::Fit(opened, near) => show(args, opened, recorded(&near)),
        Target::Nested(nested) => nested.show(args, recorded(nested.dir())),
    }
}

/// What `faba cluster-view` opens.
enum Target {
    /// One fit, and the directory a GFF is recorded in for it: its own,
    /// or the within record's, which serves every fit under it.
    Fit(Opened, PathBuf),
    /// A within record, whole.
    Nested(within::Nested),
}

/// What `inputs` open: a within record whole, the fit within its `group`,
/// or one fit as [`resolve`] opens it, which has no group to name.
fn target(
    inputs: &[Box<str>],
    modality: Option<&str>,
    output: Option<&str>,
    group: Option<&str>,
) -> anyhow::Result<Target> {
    match within::Nested::find(inputs, modality, output)? {
        Ok(nested) => match group {
            Some(group) => {
                let opened = nested.open(Some(nested.group(group)?))?;
                Ok(Target::Fit(opened, nested.dir().to_path_buf()))
            }
            None => Ok(Target::Nested(nested)),
        },
        Err(fit) => {
            anyhow::ensure!(group.is_none(), NO_GROUP);
            let opened = opened(fit, inputs, modality, output)?;
            let near = opened.0.clone();
            Ok(Target::Fit(opened, near))
        }
    }
}

/// One fit, shown as `args` ask: its page written, its levels printed, or
/// browsed; gene models from `gtf`.
fn show(
    args: &ModulesArgs,
    (dir, modality, record): Opened,
    gtf: Option<Box<str>>,
) -> anyhow::Result<()> {
    // Shared with the thread that writes the report on `e`.
    let modules = std::sync::Arc::new(Modules::load(&dir, &modality, record)?);
    if let Some(path) = &args.html {
        return write_html(&modules, gtf.as_deref(), path);
    }
    // A level asked for by name is printed as text, terminal or not.
    if args.shown.any() || crate::tui::batch(args.batch_process, "cluster-view") {
        text::printed(&modules, &args.shown)?
            .iter()
            .for_each(|line| println!("{line}"));
        return Ok(());
    }
    let models = start_gene_models(or_asked(gtf)?.as_deref());
    browse(&dir, &modules, models)
}

/// What `faba cluster-view` would open from `data_files`, refused as the
/// command refuses it, with nothing read but run records: the pages
/// `--html html` would write, one for one fit, and for a within record
/// opened whole one more for each fitted group; none without `html`.
pub(crate) fn opens(
    data_files: &[Box<str>],
    modality: Option<&str>,
    group: Option<&str>,
    html: Option<&Path>,
) -> anyhow::Result<Vec<PathBuf>> {
    let target = target(data_files, modality, None, group)?;
    let Some(html) = html else {
        return Ok(Vec::new());
    };
    Ok(match target {
        Target::Fit(..) => vec![html.to_path_buf()],
        Target::Nested(nested) => nested.pages(html).into_iter().map(|(_, p)| p).collect(),
    })
}

/// What `data_files` open, as they are: the directory, the modality and
/// the run record (its file name and contents) of a fit's run record
/// named alone, else of the block tables in `-o` (`output`), else in the
/// first input's directory. With no tables there, the fit `faba cluster`
/// would write is to be run first: nothing is ever fitted here.
pub(crate) fn resolve(
    data_files: &[Box<str>],
    modality: Option<&str>,
    output: Option<&str>,
) -> anyhow::Result<Opened> {
    opened(fit_record(data_files)?, data_files, modality, output)
}

/// [`resolve`], the fit's run record `data_files` name already read (see
/// [`fit_record`]).
fn opened(
    fit: Option<FitRecord>,
    data_files: &[Box<str>],
    modality: Option<&str>,
    output: Option<&str>,
) -> anyhow::Result<Opened> {
    let Some(fit) = fit else {
        let dir = super::output_dir(output, data_files);
        let modality = open_tables(&dir, modality)?;
        return Ok((dir, modality, None));
    };
    // A fit's run record: its tables, as they are.
    anyhow::ensure!(
        output.is_none(),
        "a run record opens its fit as it is (no -o): to fit again, run `faba cluster`"
    );
    if let Some(m) = modality {
        anyhow::ensure!(
            m.eq_ignore_ascii_case(&fit.modality),
            "the run record is of a {} fit, not --modality {m}",
            fit.modality
        );
    }
    Ok((fit.dir, fit.modality, Some(fit.record)))
}

/// The modality of the block tables in `dir` to open: those of `modality`,
/// else as [`tables`] chooses; with none there, the fit `faba cluster`
/// would write is to be run first.
fn open_tables(dir: &Path, modality: Option<&str>) -> anyhow::Result<String> {
    let count = is_count(modality);
    let asked = if count {
        Some(COUNT_MODALITY)
    } else {
        modality
    };
    tables(dir, asked).with_context(|| {
        let (what, fits) = match asked {
            Some(_) if count => ("count".to_string(), "faba cluster --modality count"),
            Some(m) => (format!("{} block", m.to_ascii_lowercase()), "faba cluster"),
            None => ("block".to_string(), "faba cluster"),
        };
        format!("no {what} tables in {}: run `{fits}` first", dir.display())
    })
}

/// The browser's header: the directory, the modality, and the model the
/// tables come from (its preset, custom settings, or a count fit).
fn title(dir: &Path, modules: &Modules) -> String {
    let model = match &modules.run {
        Ok(r) => r.label(),
        Err(why) if why == data::COUNT_FIT => "count blockmodel".to_string(),
        Err(_) => "settings not recorded".to_string(),
    };
    format!("{} · {} · {model}", dir.display(), modules.modality)
}

/// Every gene model of `gtf`; none without one.
fn all_models(
    gtf: Option<&str>,
) -> anyhow::Result<Vec<crate::site_analysis::miami::genemodel::GeneModel>> {
    match gtf {
        Some(gtf) => load_gene_models_where(gtf, |_| true),
        None => Ok(Vec::new()),
    }
}

/// `--html`: the modules as one page at `path`, gene models from `gtf`.
fn write_html(modules: &Modules, gtf: Option<&str>, path: &str) -> anyhow::Result<()> {
    let all = all_models(gtf)?;
    report::write(modules, &gene_models(modules, &all), Path::new(path))?;
    log::info!("wrote {path}");
    Ok(())
}

/// `gtf`, else the gene annotation the browser asks for (or downloads), as
/// `faba pileup` does; none when none is chosen.
fn or_asked(gtf: Option<Box<str>>) -> anyhow::Result<Option<Box<str>>> {
    match gtf {
        Some(gtf) => Ok(Some(gtf)),
        None => Ok(ask_reference("cluster-view", Reference::Annotation, None)?
            .map(|p| p.to_string_lossy().into())),
    }
}

/// The browser, until it is quit: Enter on a gene opens its full pileup,
/// and leaving that comes back to the browser as it was. `models` are the
/// gene models, as [`start_gene_models`] reads them.
fn browse(
    dir: &Path,
    modules: &std::sync::Arc<Modules>,
    models: Option<SharedModels>,
) -> anyhow::Result<()> {
    let title = title(dir, modules);
    let mut view = ModulesView::new(&title, modules)
        .with_genes(models.clone())
        .with_report(
            dir.join(format!("{}_cluster_view.html", modules.modality)),
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
        match gene::show(modules, at, models.as_ref(), &mut settings) {
            Ok(true) => return Ok(()),
            Ok(false) => {}
            Err(e) => view.set_status(format!("{e:#}")),
        }
    }
}

#[cfg(test)]
mod tests;
