//! `faba cluster`'s options: what they are, which settings are refused, and
//! the input matrices and output directory they name.

use super::{read, unified};
use clap::Args;
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;

#[derive(Args, Debug, Clone, serde::Serialize)]
pub struct ClusterArgs {
    #[arg(
        num_args = 0..,
        value_name = "MATRIX|DIR",
        value_hint = clap::ValueHint::AnyPath,
        help = "Site matrices (zarr or h5), or a faba output directory (or a run record in it); browse for one when left out",
        long_help = "Site matrices (zarr or h5): `{batch}_m6a_site` or `{batch}_atoi_site`,\n\
                     or a faba output directory (or a run record in it), whose site matrices\n\
                     of --modality are taken, one batch each. Left out, a browser asks, as\n\
                     `faba pileup`'s does: choose an output directory, or mark site matrices."
    )]
    pub data_files: Vec<Box<str>>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi|count",
        help = "Which sites to take from an output directory (default: m6a, else atoi); count: cluster the gene count matrices instead",
        long_help = "Which sites to take from an output directory (default: m6a, else atoi).\n\
                     count: a Poisson blockmodel of the `{batch}_count` matrices' spliced counts\n\
                     over the --count-genes most variable genes (Ball, Karrer and Newman\n\
                     2011; Schein et al. 2016), fitted by variational inference (Gopalan,\n\
                     Hofman and Blei 2015) from a k-means start: --cell-roles, --gene-roles\n\
                     (default 8), --svi-epochs, --batch-size, --epochs, --tol and --chains\n\
                     (restarts) apply. Its tables are count_*.parquet."
    )]
    pub(super) modality: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = 2000,
        help = "With --modality count: the most variable genes to fit"
    )]
    pub(super) count_genes: usize,

    #[arg(
        short = 'o',
        long = "output",
        value_name = "DIR",
        value_hint = clap::ValueHint::DirPath,
        help = "Where to write the tables and the run record (default: the input's directory)"
    )]
    pub(super) output: Option<Box<str>>,

    /// A pattern over the inputs' paths, so the saved command makes it
    /// absolute as a path, unless it starts with a wildcard.
    #[arg(
        long = "fit",
        value_name = "PATTERN",
        value_hint = clap::ValueHint::AnyPath,
        help = "Fit only the cells of matrices whose path matches PATTERN (`*`, `?`); assign the rest",
        long_help = "Fit only the cells of matrices whose path matches PATTERN (`*` any run,\n\
                     `?` one character), e.g. `--fit '*_wt_*'`. The other matrices' cells are\n\
                     assigned to the fitted clusters afterwards, as a control, from their reads\n\
                     as read: batch weights are not applied to them, so a control keeps its own\n\
                     level, and its `offset` measures that level against the fitted clusters.\n\
                     Default: fit all."
    )]
    pub(super) fit: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = false,
        help = "Do not weight reads to take the fitted batches out",
        long_help = "Do not weight reads to take the fitted batches out. By default, with more\n\
                     than one fitted batch, each read is weighted by the ratio of its site's\n\
                     pooled fraction to its batch's there (importance weights, capped and\n\
                     self-normalised per batch and site), so clusters cannot be made of\n\
                     batches. That also removes a site's difference between batches that\n\
                     differ in their cells."
    )]
    pub(super) no_batch_correction: bool,

    #[arg(
        long = "batch-process",
        visible_alias = "batch-mode",
        default_value_t = false,
        help = "Fit with a progress bar on stderr rather than the live view of the chains, even on a terminal"
    )]
    pub(super) batch_process: bool,

    #[arg(long, default_value_t = 50, help = "Drop cells covering fewer sites")]
    pub(super) min_cell_sites: usize,

    #[arg(
        long,
        default_value_t = 20,
        help = "Drop sites covered by fewer fitted cells"
    )]
    pub(super) min_site_cells: usize,

    #[arg(
        long = "profile-strength",
        visible_alias = "prior-strength",
        value_name = "S|inf",
        help = "Profile strength: pseudo-reads pulling a cluster's or block's site fraction toward its centre; `inf` fixes it there (default: 5, `inf` with --model blocks)",
        long_help = "Profile strength S: pseudo-reads pulling a cluster's or block's site fraction\n\
                     toward its centre, the site's pooled fraction (shifted by the block's\n\
                     shift under --shift on). `inf` fixes a block's fraction at the centre, so\n\
                     its reads are Binomial, and gives each cell its own level on top (prior sd\n\
                     --cell-sd); it needs --shift on, or no role reaches the reads. A finite\n\
                     strength integrates each block's fraction out, so no cell level can enter.\n\
                     Default: 5, `inf` with --model blocks."
    )]
    pub(super) profile_strength: Option<Strength>,

    #[arg(
        long,
        value_enum,
        default_value_t = ModelKind::Link,
        help = "link: each (cell, gene) link in a cell role (a community), cells with memberships; blocks: cell roles and gene roles (gene modules) with a methylation shift per pair; crp: each cell in one cluster",
        long_help = "link and blocks are presets of one blockmodel over (cell, gene) links, each\n\
                     link in a cell role and a gene role:\n\
                     - link: --gene-roles 1 --shift off --profile-strength 5 --global hdp\n\
                     - blocks: --gene-roles 8 --shift on --profile-strength inf --global hdp\n\
                       (with each cell's own level, --cell-sd)\n\
                     Each of those switches given overrides its preset.\n\
                     crp: each cell in one cluster, the number of clusters learned."
    )]
    pub(super) model: ModelKind,

    #[arg(
        long,
        visible_alias = "communities",
        default_value_t = 8,
        help = "Cell roles (link communities), with --model link or blocks"
    )]
    pub(super) cell_roles: usize,

    #[arg(
        long,
        help = "Gene roles (gene modules), with --model link or blocks; more than 1 needs --shift on (default: 1 with link, 8 with blocks)"
    )]
    pub(super) gene_roles: Option<usize>,

    #[arg(
        long,
        value_enum,
        value_name = "on|off",
        help = "A methylation shift per (cell role, gene role) block, on the logit scale (default: off with link, on with blocks)"
    )]
    pub(super) shift: Option<Switch>,

    #[arg(
        long,
        value_enum,
        help = "The cell roles' and gene roles' global weights: hdp, drawn with the labels so popular roles draw links and unneeded ones empty; symmetric, equal and fixed (default: hdp)"
    )]
    pub(super) global: Option<Global>,

    #[arg(
        long,
        default_value_t = 0.5,
        help = "With --model link or blocks: concentration of each cell's membership (smaller: cells condense onto fewer roles)"
    )]
    pub(super) membership_alpha: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "With --global hdp: concentration of the cell roles' and gene roles' global weights"
    )]
    pub(super) community_gamma: f64,

    #[arg(
        long,
        default_value_t = 0.5,
        help = "With --model link or blocks: concentration of each gene's role weights"
    )]
    pub(super) gene_alpha: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "With --shift on: prior sd of each (cell role, gene role) block's methylation shift, on the logit scale"
    )]
    pub(super) shift_sd: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "At --profile-strength inf (--model blocks): prior sd of each cell's own methylation level, on the logit scale"
    )]
    pub(super) cell_sd: f64,

    #[arg(
        long,
        default_value_t = 0.0,
        help = "Share of (cell, site) entries held out of the fit and scored at every kept draw, to compare models; the written tables then hold only the fitted share of the reads"
    )]
    pub(super) holdout: f64,

    #[arg(
        long,
        default_value_t = 10,
        help = "Smallest cluster reported: smaller ones' cells go to their most probable larger cluster"
    )]
    pub(super) min_cluster_cells: usize,

    #[arg(
        long,
        default_value_t = 4,
        help = "Independent chains, run in parallel; for a fit by variational inference \
                (--modality count, --model blocks), restarts, each from k-means seeds of its own, \
                the starts in parallel, then the fits one after another, the best kept"
    )]
    pub(super) chains: usize,

    #[arg(
        long,
        value_name = "N",
        default_value_t = crate::cluster::count::vi::Options::default().svi_epochs,
        help = "A fit by variational inference (--modality count, --model blocks): stochastic (minibatch) epochs before coordinate ascent"
    )]
    pub(super) svi_epochs: usize,

    #[arg(
        long,
        value_name = "N",
        default_value_t = crate::cluster::count::vi::Options::default().batch_size,
        help = "A fit by variational inference: cells per minibatch; every cell in one batch is plain coordinate ascent"
    )]
    pub(super) batch_size: usize,

    #[arg(
        long,
        value_name = "N",
        default_value_t = crate::cluster::count::vi::Options::default().epochs,
        help = "A fit by variational inference: epochs at most, stochastic and coordinate ascent together"
    )]
    pub(super) epochs: usize,

    #[arg(
        long,
        default_value_t = crate::cluster::count::vi::Options::default().tol,
        help = "A fit by variational inference (--modality count, --model blocks): converged when an epoch of coordinate ascent moves the ELBO by less than this per count (per read for --model blocks)"
    )]
    pub(super) tol: f64,

    #[arg(
        long,
        default_value_t = 150,
        help = "Sweeps per chain, burn-in included"
    )]
    pub(super) sweeps: usize,

    #[arg(
        long,
        default_value_t = 50,
        help = "Sweeps discarded at the start of each chain"
    )]
    pub(super) burn_in: usize,

    #[arg(
        long,
        default_value_t = 2,
        help = "Keep every this-many-th sweep after the burn-in"
    )]
    pub(super) thin: usize,

    #[arg(long, default_value_t = 10, help = "Split-merge proposals per sweep")]
    pub(super) split_merge: usize,

    #[arg(
        long,
        default_value_t = 2,
        help = "Restricted Gibbs scans before each split-merge proposal"
    )]
    pub(super) restricted_scans: usize,

    #[arg(
        long,
        default_value_t = 20,
        help = "Chains start from k-means into a number of clusters drawn from 1 up to this"
    )]
    pub(super) start_clusters: usize,

    #[arg(
        long,
        default_value_t = 1,
        help = "Random seed (chain i uses seed + i)"
    )]
    pub(super) seed: u64,

    #[arg(
        long = "within",
        value_name = "RECORD",
        value_hint = clap::ValueHint::FilePath,
        help = "Fit within each cell module of a count fit, named by its run record (count_blockmodel.run.json)",
        long_help = "The second round of two: fit the sites within each cell module of a count\n\
                     fit, named by its run record (`count_blockmodel.run.json`, as `faba cluster\n\
                     --modality count` writes it). Each cell goes to its main count module; each\n\
                     module of at least --within-min-cells cells is fitted on its own, one after\n\
                     the other, into `within/{modality}/{group}/` under the output directory\n\
                     (default: the count fit's), with `{modality}_within.run.json` naming them all and\n\
                     `{modality}_within_cells.parquet` giving each cell its nested label (`C3.2`).\n\
                     The site matrices default to those beside the count fit's inputs; --model\n\
                     defaults to blocks. Open the result with `faba cluster-view`."
    )]
    pub(super) within: Option<Box<str>>,

    #[arg(
        long = "within-min-cells",
        requires = "within",
        value_name = "N",
        default_value_t = 50,
        help = "With --within: a group with fewer cells to fit is skipped"
    )]
    pub(super) within_min_cells: usize,

    /// The options given, for the command saved beside the run.
    #[arg(skip)]
    #[serde(skip)]
    pub(super) given: super::command::Given,
}

/// Which model `faba cluster` fits: the mixture, or one of the
/// blockmodel's two presets.
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelKind {
    Crp,
    Link,
    Blocks,
}

impl ModelKind {
    /// Its name as `--model` takes it, the one spelling of it everywhere.
    pub fn name(self) -> String {
        use clap::ValueEnum;
        self.to_possible_value()
            .expect("no model is skipped")
            .get_name()
            .to_string()
    }
}

impl serde::Serialize for ModelKind {
    /// By its name, so the saved command says it as `--model` takes it.
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(&self.name())
    }
}

/// The blockmodel's block shift, on or off.
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub(super) enum Switch {
    On,
    Off,
}

/// The blockmodel's global weights over the cell roles.
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub(super) enum Global {
    Hdp,
    Symmetric,
}

/// A profile strength: a number, or `inf`, saved as written.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Strength(pub(super) f64);

impl std::str::FromStr for Strength {
    type Err = std::num::ParseFloatError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.parse().map(Strength)
    }
}

impl serde::Serialize for Strength {
    /// A number, or `"inf"` (see [`unified::strength`]), so the saved
    /// command says it.
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        unified::strength::serialize(&self.0, s)
    }
}

/// The profile strength of the mixture and of the link preset.
const DEFAULT_STRENGTH: f64 = 5.0;

/// The mixture's profile strength.
pub(super) fn crp_strength(args: &ClusterArgs) -> f64 {
    args.profile_strength.map_or(DEFAULT_STRENGTH, |s| s.0)
}

/// The link preset's value, or with `blocks` the blocks preset's.
fn pick<T>(blocks: bool, link: T, blocks_preset: T) -> T {
    if blocks {
        blocks_preset
    } else {
        link
    }
}

/// The blockmodel `args` fit: their preset (`--model link` or `blocks`),
/// each switch given overriding it.
pub(super) fn blockmodel(args: &ClusterArgs) -> unified::Settings {
    let blocks = args.model == ModelKind::Blocks;
    let shift = args.shift.unwrap_or(pick(blocks, Switch::Off, Switch::On));
    unified::Settings {
        cell_roles: args.cell_roles,
        gene_roles: args.gene_roles.unwrap_or(pick(blocks, 1, 8)),
        alpha: args.membership_alpha,
        gene_alpha: args.gene_alpha,
        weights: match args.global.unwrap_or(Global::Hdp) {
            Global::Hdp => unified::Weights::Hdp(args.community_gamma),
            Global::Symmetric => unified::Weights::Symmetric,
        },
        shift_sd: (shift == Switch::On).then_some(args.shift_sd),
        strength: args
            .profile_strength
            .map_or(pick(blocks, DEFAULT_STRENGTH, f64::INFINITY), |s| s.0),
        cell_sd: args.cell_sd,
    }
}

/// The preset `args` fit, when their switches leave it as it is: given
/// or not, each the preset's.
fn preset(args: &ClusterArgs) -> Option<ModelKind> {
    let preset = ClusterArgs {
        gene_roles: None,
        shift: None,
        global: None,
        profile_strength: None,
        ..args.clone()
    };
    (args.model != ModelKind::Crp && blockmodel(&preset) == blockmodel(args)).then_some(args.model)
}

/// The blockmodel `args` fit, as a run's record says it.
pub(super) fn recorded(args: &ClusterArgs) -> unified::Recorded {
    unified::Recorded {
        preset: preset(args).map(ModelKind::name),
        settings: blockmodel(args),
    }
}

/// The matrices to read: each with its batch and whether it is fitted; and
/// the modality asked for.
pub(super) fn inputs(args: &ClusterArgs) -> anyhow::Result<Vec<read::Input>> {
    let mut out = Vec::new();
    let fits = |path: &str| {
        args.fit
            .as_deref()
            .is_none_or(|p| crate::common::wildcard(p, path))
    };
    for f in &args.data_files {
        let folder: Box<str> = crate::run_record::default_output_dir(f)
            .to_string_lossy()
            .into();
        match crate::run_record::output_dir_of(f) {
            Some(dir) => {
                let matrices = if args.is_count() {
                    count_matrices(&dir)?
                } else {
                    crate::site_analysis::output_dir::output_site_matrices(
                        &dir,
                        args.modality.as_deref(),
                    )?
                    .matrices
                };
                for (batch, path) in matrices {
                    out.push(read::Input {
                        fit: fits(&path),
                        path,
                        folder: folder.clone(),
                        batch,
                    });
                }
            }
            None => {
                let name = crate::qc::layout::file_name(f);
                out.push(read::Input {
                    path: f.clone(),
                    folder,
                    batch: crate::qc::layout::matrix_batch(&name).unwrap_or(name),
                    fit: fits(f),
                });
            }
        }
    }
    anyhow::ensure!(
        out.iter().any(|i| i.fit),
        "no input matrix matches --fit {}",
        args.fit.as_deref().unwrap_or("")
    );
    qualify_batches(&mut out);
    Ok(out)
}

/// The `{batch}_count` matrices of output directory `dir`, as (batch, path).
fn count_matrices(dir: &std::path::Path) -> anyhow::Result<Vec<(Box<str>, Box<str>)>> {
    let layout = crate::qc::layout::scan_input_dir(&dir.to_string_lossy())?;
    let found: Vec<_> = (layout.matrices.iter())
        .filter(|m| &*m.kind == super::command::COUNT_MODALITY)
        .map(|m| (m.batch.clone(), m.path.clone()))
        .collect();
    anyhow::ensure!(!found.is_empty(), "{} has no count matrices", dir.display());
    Ok(found)
}

/// A batch name met in several folders is a different batch in each: those
/// are named `{folder}/{batch}`.
pub(super) fn qualify_batches(inputs: &mut [read::Input]) {
    let mut folders: FxHashMap<Box<str>, FxHashSet<Box<str>>> = FxHashMap::default();
    for i in inputs.iter() {
        folders
            .entry(i.batch.clone())
            .or_default()
            .insert(i.folder.clone());
    }
    for i in inputs.iter_mut() {
        if folders[&i.batch].len() > 1 {
            i.batch = format!("{}/{}", i.folder, i.batch).into();
        }
    }
}

/// `-o` (`output`), else the first input directory, else the first
/// matrix's directory.
pub(super) fn output_dir(output: Option<&str>, data_files: &[Box<str>]) -> PathBuf {
    match output {
        Some(o) => PathBuf::from(o),
        None => crate::run_record::default_output_dir(&data_files[0]),
    }
}

/// Whether `args` fit by variational inference: the count fit, and the
/// site blockmodel at infinite strength (the blocks preset's).
pub(super) fn is_vi(args: &ClusterArgs) -> bool {
    args.is_count() || (args.model == ModelKind::Blocks && blockmodel(args).strength.is_infinite())
}

/// Settings that would make the fit meaningless, refused up front.
pub(super) fn check(args: &ClusterArgs) -> anyhow::Result<()> {
    anyhow::ensure!(args.chains >= 1, "--chains must be at least 1");
    let positive = |x: f64| x > 0.0 && x.is_finite();
    if is_vi(args) {
        // A fit by variational inference: no draws to keep.
        for flag in ["sweeps", "burn_in", "thin"] {
            anyhow::ensure!(
                !args.given(flag),
                "--{} sets a sampler's draws; this fit is by variational inference \
                 (--epochs, --batch-size, --tol)",
                flag.replace('_', "-")
            );
        }
        anyhow::ensure!(
            args.holdout == 0.0,
            "--holdout is not for a fit by variational inference"
        );
        // Its role weights' global weights are fixed: nothing to draw.
        for flag in ["global", "community_gamma"] {
            anyhow::ensure!(
                !args.given(flag),
                "--{} sets the global weights a sampler draws; this fit by variational \
                 inference keeps them uniform",
                flag.replace('_', "-")
            );
        }
        anyhow::ensure!(
            args.epochs >= 1 && args.batch_size >= 1,
            "--epochs and --batch-size must be at least 1"
        );
        anyhow::ensure!(
            args.tol >= 0.0 && args.tol.is_finite(),
            "--tol must be 0 or more"
        );
    } else {
        for flag in ["svi_epochs", "batch_size", "epochs", "tol"] {
            anyhow::ensure!(
                !args.given(flag),
                "--{} applies to a fit by variational inference only (--modality count, --model blocks)",
                flag.replace('_', "-")
            );
        }
        // Split R-hat needs at least two draws in each half of every chain.
        anyhow::ensure!(
            args.sweeps.saturating_sub(args.burn_in) / args.thin.max(1) >= 4,
            "--sweeps less --burn-in must leave at least 4 draws per chain (after --thin)"
        );
    }
    anyhow::ensure!(
        (0.0..1.0).contains(&args.holdout),
        "--holdout must be at least 0 and less than 1"
    );
    if args.within.is_some() {
        anyhow::ensure!(
            !args.is_count(),
            "--within fits sites within a count fit's cell modules: give --modality m6a or atoi, not count"
        );
        anyhow::ensure!(
            args.model != ModelKind::Crp,
            "--within fits the blockmodel (--model link or blocks), not --model crp"
        );
        anyhow::ensure!(
            args.within_min_cells >= 2,
            "--within-min-cells must be at least 2"
        );
    }
    if args.is_count() {
        anyhow::ensure!(
            args.cell_roles >= 1 && args.gene_roles.is_none_or(|b| b >= 1),
            "--cell-roles and --gene-roles must be at least 1"
        );
        anyhow::ensure!(args.count_genes >= 1, "--count-genes must be at least 1");
        // What only the site fits do, refused rather than ignored.
        anyhow::ensure!(
            args.holdout == 0.0,
            "--holdout does not apply to --modality count"
        );
        anyhow::ensure!(
            !args.given("model"),
            "--model chooses a site fit's model; --modality count fits the count blockmodel"
        );
        anyhow::ensure!(
            !args.given("community_gamma"),
            "--community-gamma does not apply to --modality count (its role weights have a fixed shape)"
        );
        return Ok(());
    }
    if args.model == ModelKind::Crp {
        anyhow::ensure!(
            positive(crp_strength(args)),
            "--profile-strength must be positive and finite with --model crp"
        );
        return Ok(());
    }
    let s = blockmodel(args);
    anyhow::ensure!(s.strength > 0.0, "--profile-strength must be positive");
    s.check()?;
    anyhow::ensure!(positive(args.gene_alpha), "--gene-alpha must be positive");
    anyhow::ensure!(positive(args.shift_sd), "--shift-sd must be positive");
    anyhow::ensure!(positive(args.cell_sd), "--cell-sd must be positive");
    anyhow::ensure!(
        positive(args.membership_alpha),
        "--membership-alpha must be positive"
    );
    anyhow::ensure!(
        positive(args.community_gamma),
        "--community-gamma must be positive"
    );
    Ok(())
}

/// `data_files`, or when none are named the inputs asked for in a browser:
/// the files `keep` accepts, called `what`; none when none are chosen.
pub(super) fn with_inputs(
    data_files: &[Box<str>],
    command: &'static str,
    batch: bool,
    (what, keep): (&'static str, fn(&str) -> bool),
) -> anyhow::Result<Option<Vec<Box<str>>>> {
    if !data_files.is_empty() {
        return Ok(Some(data_files.to_vec()));
    }
    let ask = crate::site_analysis::input_picker::ask_inputs(command, batch, what, keep, false);
    Ok(ask?.map(|chosen| chosen.paths))
}
