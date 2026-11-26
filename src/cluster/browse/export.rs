//! `--screen`: the browser's screen written as a figure (see
//! [`super::figure::screen`]), of one fit or of many in one call: every
//! fit under a sweep's folder (`--each-fit`), every group of a within
//! record (`--each-group`), every gene module (`--each-module`), each file
//! named for what it shows, `{fit}_{group}_{module}_screen.pdf`.

use super::figure::screen::{self, Beside, Preset, Region, State};
use super::{gene_models, go, Modules, ModulesArgs, Opened, Target};
use crate::site_analysis::miami::genemodel::GeneModel;
use anyhow::Context;
use rustc_hash::FxHashMap;
use std::path::{Path, PathBuf};

/// A screen file's kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum Format {
    Pdf,
    Svg,
    Png,
}

impl Format {
    fn ext(self) -> &'static str {
        match self {
            Format::Pdf => "pdf",
            Format::Svg => "svg",
            Format::Png => "png",
        }
    }
}

/// The screen's options: where it goes, what it shows, and how many.
#[derive(clap::Args, Debug, Clone, Default, PartialEq)]
pub struct ScreenArgs {
    #[arg(
        long = "screen",
        value_name = "PREFIX|DIR",
        value_hint = clap::ValueHint::AnyPath,
        conflicts_with_all = ["html", "figure"],
        help = "Write the browser's screen (the modules grid, the sites or GO beside it, the pinned genes) as PREFIX.pdf and PREFIX.png, laid out and drawn as the browser draws it, instead of browsing; with --each-fit, --each-group or --each-module, DIR holds a file each, {fit}_{group}_{module}_screen.pdf"
    )]
    pub screen: Option<Box<str>>,

    #[arg(
        long = "region",
        value_enum,
        default_value_t = Region::All,
        requires = "screen",
        help = "Which of the screen --screen writes: all of it, the modules grid, what is beside it (--side), or the pinned genes"
    )]
    pub region: Region,

    #[arg(
        long = "side",
        value_enum,
        requires = "screen",
        help = "Beside the modules grid on --screen: the sites of the gene module --module names (default: the grid's first row), each gene module's top GO terms on its row with the q beside each (--go-labels of them, default 3, fewer where a row cannot hold them; --gaf, --obo, --gmt, --go-q; a row with none blank, or saying so with --go-none), or nothing"
    )]
    pub side: Option<Beside>,

    #[arg(
        long = "pin",
        value_name = "GENE,...",
        value_delimiter = ',',
        requires = "screen",
        help = "The genes pinned under the grid on --screen, by key or symbol, in order, each in its main gene module; with none pinned (no --pin or --pin-top) the screen has no genes region, the grid and what is beside it taking the page"
    )]
    pub pin: Vec<String>,

    #[arg(
        long = "pin-top",
        value_name = "N",
        conflicts_with = "pin",
        requires = "screen",
        value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..),
        help = "Pin each gene module's first N genes (the most sites first) on --screen, the modules in the grid's order; with --each-module, that module's"
    )]
    pub pin_top: Option<usize>,

    #[arg(
        long = "preset",
        value_enum,
        requires = "screen",
        help = "The grid on --screen as the browser draws it: diagonal (rows and columns in diagonal order, a Hinton square of z a block, red above the row's mean and blue below; default) or blocks (rows by shift range, each block's reads as mirror bars)"
    )]
    pub preset: Option<Preset>,

    #[arg(
        long = "genes-yscale",
        value_enum,
        requires = "screen",
        help = "The scale of the pinned genes' mirror bars on --screen (each side the posterior-mean reads at a site, summed over a track's cells): linear (default, as the browser opens), sqrt, or log10 (log10(1 + reads); log1p names it too), as faba pileup's; the \"all\" track on its own peak, the cell modules on one shared; in the browser, y on the genes pane steps them"
    )]
    pub genes_yscale: Option<screen::YScale>,

    /// How the pinned genes are drawn on --screen.
    #[command(flatten)]
    pub genes: screen::GenesLook,

    #[arg(
        long = "cell-order",
        value_name = "C,...",
        value_delimiter = ',',
        requires = "screen",
        help = "The cell modules first on --screen, in this order (C6,C2 or 6,2), the rest after them in diagonal order"
    )]
    pub cell_order: Vec<String>,

    #[arg(
        long = "state",
        value_name = "FILE",
        value_hint = clap::ValueHint::FilePath,
        requires = "screen",
        help = "A saved view's state (the JSON the browser's s writes beside its figure): its preset, cell order, module, side and pins, each overridden by the option given"
    )]
    pub state: Option<Box<str>>,

    #[arg(
        long = "format",
        value_enum,
        value_delimiter = ',',
        default_values = ["pdf", "png"],
        requires = "screen",
        help = "The files --screen writes: pdf (vector), svg, png"
    )]
    pub format: Vec<Format>,

    #[arg(
        long = "each-fit",
        requires = "screen",
        help = "A screen for every fit's run record under the folder named (a sweep's), each named by its folder there ({fit}); with --each-group, every group of every within record under it instead"
    )]
    pub each_fit: bool,

    #[arg(
        long = "each-group",
        requires = "screen",
        help = "A screen for every fitted group of a within record ({group})"
    )]
    pub each_group: bool,

    #[arg(
        long = "each-module",
        requires = "screen",
        help = "A screen for every gene module with sites ({module}, as M3), its sites beside the grid"
    )]
    pub each_module: bool,
}

impl ScreenArgs {
    /// Whether the screens go into a directory, a file each.
    fn many(&self) -> bool {
        self.each_fit || self.each_group || self.each_module
    }
}

/// A fit to draw: what names it in a file name, and how it is opened.
struct Fit {
    fit: Option<String>,
    group: Option<String>,
    opened: Opened,
    /// Where a GFF is recorded for it.
    near: PathBuf,
}

/// A screen to write: its fit (a place in the fits), its state, and its
/// files.
struct Job {
    fit: usize,
    state: State,
    paths: Vec<PathBuf>,
}

/// `--screen`: every screen `args` ask for written; the paths printed.
/// Every file is checked new before the first is drawn.
pub(super) fn run(args: &ModulesArgs, inputs: &[Box<str>]) -> anyhow::Result<()> {
    let s = &args.screen;
    let out = s.screen.as_deref().expect("--screen");
    let sheet = crate::figure::style::Sheet {
        page: super::figure::page(super::figure::Plot::Modules, &args.size)?,
        look: args.style,
    };
    args.style.check()?;
    args.grid.check()?;
    let base = match &s.state {
        Some(path) => read_state(Path::new(path.as_ref()))?,
        None => State::default(),
    };
    let fits = fits(args, inputs)?;
    anyhow::ensure!(!fits.is_empty(), "no fit to draw a screen of");
    let mut loaded = Vec::new();
    let mut jobs = Vec::new();
    for (k, fit) in fits.iter().enumerate() {
        let (dir, modality, record) = fit.opened.clone();
        let mut m = Modules::load(&dir, &modality, record)?;
        m.regroup(&args.batch_groups()?)?;
        let state = asked_state(args, &base, &m)?;
        state.order(&mut m)?;
        let modules: Vec<Option<usize>> = match s.each_module {
            true => rows_of(&m, &state).into_iter().map(Some).collect(),
            false => vec![None],
        };
        for module in modules {
            let mut state = state.clone();
            if let Some(g) = module {
                state.module = Some(format!("M{}", g + 1));
                state.side = match state.side {
                    Beside::None => Beside::Sites,
                    side => side,
                };
                if let Some(n) = s.pin_top {
                    state.pins = screen::top_genes(&m, &[g], n);
                }
            }
            let parts = [
                fit.fit.clone(),
                fit.group.clone(),
                module.map(|g| format!("M{}", g + 1)),
            ];
            let paths = paths(out, s, &parts);
            jobs.push(Job {
                fit: k,
                state,
                paths,
            });
        }
        loaded.push(m);
    }
    for path in jobs.iter().flat_map(|j| &j.paths) {
        anyhow::ensure!(
            !path.exists(),
            "{} exists: name a new --screen",
            path.display()
        );
    }
    let mut drawn = Vec::with_capacity(jobs.len());
    let mut models_of: FxHashMap<Option<Box<str>>, Vec<GeneModel>> = FxHashMap::default();
    for job in &jobs {
        let (fit, m) = (&fits[job.fit], &loaded[job.fit]);
        let placed = job.state.place(m)?;
        let gtf = crate::run_record::explicit_or_recorded(
            args.gtf.as_deref(),
            &fit.near.to_string_lossy(),
            "gff",
            "gene models",
        );
        // Gene models only where genes are pinned, each annotation read once.
        let models = match placed.pins.is_empty() {
            true => FxHashMap::default(),
            false => {
                if !models_of.contains_key(&gtf) {
                    let all = super::all_models(gtf.as_deref())?;
                    models_of.insert(gtf.clone(), all);
                }
                gene_models(m, &models_of[&gtf])
            }
        };
        let go_table = match placed.side {
            Beside::Go => {
                let source = args.shown.go_ask().source(m).map_err(anyhow::Error::msg)?;
                let ask = args.shown.go_ask();
                Some(go::table(m, &source, ask.min_weight, ask.q)?)
            }
            _ => None,
        };
        let title = super::figure_title(&fit.opened.0, m, fit.group.as_deref());
        let ask = screen::Ask {
            go: go_table.as_ref(),
            models: &models,
            grid: &args.grid,
            region: s.region,
            genes: s.genes,
            title: &title,
            sheet,
        };
        let svg = screen::draw(m, &placed, &ask)
            .map_err(anyhow::Error::msg)
            .with_context(|| format!("{} not drawn", job.paths[0].display()))?;
        drawn.push(svg);
    }
    // Every screen drawn before any is written: one refused writes none.
    if s.many() {
        std::fs::create_dir_all(out).with_context(|| format!("--screen {out}"))?;
    }
    for (job, svg) in jobs.iter().zip(&drawn) {
        for path in &job.paths {
            write(svg, path)?;
            println!("{}", path.display());
        }
    }
    Ok(())
}

/// `svg` written as `path`'s kind says.
fn write(svg: &str, path: &Path) -> anyhow::Result<()> {
    let p = path.to_string_lossy();
    match path.extension().and_then(|e| e.to_str()) {
        Some("svg") => std::fs::write(path, svg)?,
        Some("png") => std::fs::write(path, crate::figure::png(svg, 2.0)?)?,
        _ => {
            crate::figure::save(svg, &p)?;
        }
    }
    Ok(())
}

/// The files a screen is written as: `out` with each format's extension,
/// or in the directory `out` named by `parts`, those given, then
/// `screen`.
fn paths(out: &str, s: &ScreenArgs, parts: &[Option<String>]) -> Vec<PathBuf> {
    let stem = match s.many() {
        true => {
            let named: Vec<&str> = parts.iter().flatten().map(String::as_str).collect();
            Path::new(out).join(format!("{}_screen", named.join("_")))
        }
        false => {
            let bare = [".pdf", ".svg", ".png"]
                .iter()
                .find_map(|e| out.strip_suffix(e))
                .unwrap_or(out);
            PathBuf::from(bare)
        }
    };
    (s.format.iter())
        .map(|f| PathBuf::from(format!("{}.{}", stem.display(), f.ext())))
        .collect()
}

/// The grid's rows, in the order the state draws them: the modules
/// `--each-module` draws a screen of.
fn rows_of(m: &Modules, state: &State) -> Vec<usize> {
    match state.preset {
        Preset::Blocks => m.order.clone(),
        Preset::Diagonal => m.diagonal.clone(),
    }
}

/// The state `args` ask for of `m`: `base` (a saved view's, or the
/// browser's own as it opens) with each option given in its place.
fn asked_state(args: &ModulesArgs, base: &State, m: &Modules) -> anyhow::Result<State> {
    let s = &args.screen;
    let mut state = base.clone();
    if let Some(preset) = s.preset {
        state.preset = preset;
    }
    if !s.cell_order.is_empty() {
        state.cell_order = s.cell_order.clone();
    }
    if let Some(module) = &args.shown.module {
        state.module = Some(module.clone());
    }
    if let Some(side) = s.side {
        state.side = side;
    }
    if let Some(n) = args.grid.go_labels {
        state.go_terms = Some(n);
    }
    if !s.pin.is_empty() {
        state.pins = s.pin.clone();
    }
    if let Some(scale) = s.genes_yscale {
        state.genes_yscale = scale;
    }
    if let Some(n) = s.pin_top {
        state.pins = screen::top_genes(m, &rows_of(m, &state), n);
    }
    // Refused here, before any file is written, for a name `m` lacks.
    state.place(m)?;
    Ok(state)
}

/// A saved view's state: the browser's `s`'s, or a live view's export's
/// (its `view` and `command` said beside the state, see
/// [`crate::tui::live::snapshot`]).
fn read_state(path: &Path) -> anyhow::Result<State> {
    let text = std::fs::read_to_string(path).with_context(|| format!("{}", path.display()))?;
    let mut v: serde_json::Value = serde_json::from_str(&text)
        .with_context(|| format!("{} is no saved view state", path.display()))?;
    crate::tui::live::check_view(&v, "cluster_view")?;
    if let Some(m) = v.as_object_mut() {
        m.remove("view");
        m.remove("command");
    }
    serde_json::from_value(v).with_context(|| format!("{} is no saved view state", path.display()))
}

/// The fits `args` draw screens of: every fit under `--each-fit`, every
/// group of a within record with `--each-group`, else the one `inputs`
/// open (a within record's count fit without `--group`).
fn fits(args: &ModulesArgs, inputs: &[Box<str>]) -> anyhow::Result<Vec<Fit>> {
    let s = &args.screen;
    let modality = args.modality.as_deref();
    if !s.each_fit {
        return of_target(args, inputs, None);
    }
    let [root] = inputs else {
        anyhow::bail!("--each-fit takes one folder, the fits under it: name it alone");
    };
    let root = Path::new(root.as_ref());
    let records = records_under(root, s.each_group)?;
    anyhow::ensure!(
        !records.is_empty(),
        "no {} under {}",
        if s.each_group {
            "within record"
        } else {
            "fit's run record"
        },
        root.display()
    );
    let mut out = Vec::new();
    for record in records {
        let dir = record.parent().unwrap_or(root);
        let name = match dir.strip_prefix(root) {
            Ok(rel) if !rel.as_os_str().is_empty() => (rel.components())
                .map(|c| c.as_os_str().to_string_lossy().into_owned())
                .collect::<Vec<_>>()
                .join("-"),
            _ => root
                .file_name()
                .map_or("fit".into(), |n| n.to_string_lossy().into_owned()),
        };
        if let (Some(m), false) = (modality, s.each_group) {
            // A sweep's folders may hold fits of several modalities.
            let fit = super::read_fit_record(&record)?;
            if fit.is_some_and(|f| !f.modality.eq_ignore_ascii_case(m)) {
                continue;
            }
        }
        let input = [Box::<str>::from(record.to_string_lossy())];
        out.extend(of_target(args, &input, Some(name))?);
    }
    Ok(out)
}

/// The fits `inputs` open as `args` ask, `fit` naming them in file names.
fn of_target(
    args: &ModulesArgs,
    inputs: &[Box<str>],
    fit: Option<String>,
) -> anyhow::Result<Vec<Fit>> {
    let s = &args.screen;
    let (modality, output) = (args.modality.as_deref(), args.output.as_deref());
    match super::target(inputs, modality, output, args.group.as_deref())? {
        Target::Fit(opened, near) => {
            anyhow::ensure!(
                !s.each_group,
                "--each-group draws the groups of a within record: name one"
            );
            let group = args.group.clone();
            Ok(vec![Fit {
                fit,
                group: group.filter(|_| s.many()),
                opened,
                near,
            }])
        }
        Target::Nested(nested) => {
            let near = nested.dir().to_path_buf();
            if !s.each_group {
                return Ok(vec![Fit {
                    fit,
                    group: None,
                    opened: nested.open(None)?,
                    near,
                }]);
            }
            (nested.fitted())
                .map(|g| {
                    Ok(Fit {
                        fit: fit.clone(),
                        group: Some(g.group.clone()),
                        opened: nested.open(Some(g))?,
                        near: near.clone(),
                    })
                })
                .collect()
        }
    }
}

/// Every run record under `root`, sorted: fits' (`faba cluster`'s
/// blockmodel records), or within records with `within`; none under a
/// within record's own fits, which `--each-group` draws.
fn records_under(root: &Path, within: bool) -> anyhow::Result<Vec<PathBuf>> {
    let mut found = Vec::new();
    let mut dirs = vec![root.to_path_buf()];
    while let Some(dir) = dirs.pop() {
        let entries = std::fs::read_dir(&dir).with_context(|| format!("{}", dir.display()))?;
        for entry in entries {
            let path = entry?.path();
            let name = path.file_name().map(|n| n.to_string_lossy().into_owned());
            let name = name.unwrap_or_default();
            if path.is_dir() && !crate::qc::layout::is_site_matrix_name(&name) {
                dirs.push(path);
                continue;
            }
            let wanted = match within {
                true => name.ends_with("_within.run.json"),
                false => name.ends_with("_blockmodel.run.json"),
            };
            if wanted {
                found.push(path);
            }
        }
    }
    found.sort();
    Ok(found)
}
