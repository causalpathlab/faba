//! `faba cluster`'s options: what they are, which settings are refused, and
//! the input matrices and output directory they name.

use super::{read, settings};
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
        help = "Fit with a progress bar on stderr rather than the live view of the restarts, even on a terminal"
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
        long,
        visible_alias = "communities",
        default_value_t = 8,
        help = "Cell roles (cell modules)"
    )]
    pub(super) cell_roles: usize,

    #[arg(long, help = "Gene roles (gene modules) (default: 8)")]
    pub(super) gene_roles: Option<usize>,

    #[arg(
        long,
        default_value_t = 0.5,
        help = "Concentration of each cell's role weights (smaller: cells condense onto fewer roles)"
    )]
    pub(super) membership_alpha: f64,

    #[arg(
        long,
        default_value_t = 0.5,
        help = "Concentration of each gene's role weights"
    )]
    pub(super) gene_alpha: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "A site fit: prior sd of each (cell role, gene role) block's methylation shift, on the logit scale"
    )]
    pub(super) shift_sd: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "A site fit: prior sd of each cell's own methylation level, on the logit scale"
    )]
    pub(super) cell_sd: f64,

    #[arg(
        long,
        default_value_t = 4,
        help = "Restarts, each from k-means seeds of its own: the starts in parallel, then \
                the fits one after another, the best by ELBO kept"
    )]
    pub(super) chains: usize,

    #[arg(
        long,
        value_name = "N",
        default_value_t = crate::cluster::vi::Options::default().svi_epochs,
        help = "Stochastic (minibatch) epochs before coordinate ascent"
    )]
    pub(super) svi_epochs: usize,

    #[arg(
        long,
        value_name = "N",
        default_value_t = crate::cluster::vi::Options::default().batch_size,
        help = "Cells per minibatch; every cell in one batch is plain coordinate ascent"
    )]
    pub(super) batch_size: usize,

    #[arg(
        long,
        value_name = "N",
        default_value_t = crate::cluster::vi::Options::default().epochs,
        help = "Epochs at most, stochastic and coordinate ascent together"
    )]
    pub(super) epochs: usize,

    #[arg(
        long,
        default_value_t = crate::cluster::vi::Options::default().tol,
        help = "Converged when an epoch of coordinate ascent moves the ELBO by less than this per count (per read for a site fit)"
    )]
    pub(super) tol: f64,

    #[arg(
        long,
        default_value_t = 1,
        help = "Random seed (restart i uses seed + i)"
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
                     The site matrices default to those beside the count fit's inputs. Open the\n\
                     result with `faba cluster-view`."
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

    /// The MCMC fits' flags (faba before 0.19), refused when given (see
    /// [`refuse_removed`]); hidden, and never recorded.
    #[command(flatten)]
    #[serde(skip)]
    pub(super) removed: Removed,

    /// The options given, for the command saved beside the run.
    #[arg(skip)]
    #[serde(skip)]
    pub(super) given: super::command::Given,
}

/// The removed flags (see [`Removed`]), each once: its field, its name on
/// the command line, and the old name it also answered to.
macro_rules! removed {
    ($($field:ident = $flag:literal $(| $alias:literal)?),* $(,)?) => {
        /// The flags of the MCMC fits faba had before 0.19 (`--model link`,
        /// `crp`, the sampler's draws, the held-out reads), kept only to be
        /// refused by name; `--model blocks`, the only model left, is taken
        /// as it is.
        #[derive(Args, Debug, Clone, Default)]
        pub struct Removed {
            $(
                #[arg(long = $flag, hide = true $(, alias = $alias)?)]
                $field: Option<String>,
            )*
        }

        impl Removed {
            /// Each removed flag given, whatever its value.
            fn given_any(&self) -> Vec<&'static str> {
                let mut given = Vec::new();
                $(
                    if self.$field.is_some() {
                        given.push(concat!("--", $flag));
                    }
                )*
                given
            }

            /// Each removed flag given that is refused (see [`refused`]),
            /// by its name on the command line.
            fn given(&self) -> Vec<&'static str> {
                let mut given = Vec::new();
                $(
                    let flag = concat!("--", $flag);
                    if self.$field.as_deref().is_some_and(|v| refused(flag, v)) {
                        given.push(flag);
                    }
                )*
                given
            }
        }

        /// The removed flags' names on the command line, old names too.
        const REMOVED: &[&str] = &[$(concat!("--", $flag) $(, concat!("--", $alias))?),*];
    };
}

removed! {
    model = "model",
    profile_strength = "profile-strength" | "prior-strength",
    shift = "shift",
    global = "global",
    community_gamma = "community-gamma",
    holdout = "holdout",
    min_cluster_cells = "min-cluster-cells",
    sweeps = "sweeps",
    burn_in = "burn-in",
    thin = "thin",
    split_merge = "split-merge",
    restricted_scans = "restricted-scans",
    start_clusters = "start-clusters",
}

/// Whether `flag` (as `--name`) is one of the MCMC fits' flags faba no
/// longer has.
pub(crate) fn is_removed(flag: &str) -> bool {
    REMOVED.contains(&flag)
}

/// The removed flags' values that ask for the model left (the blocks
/// preset as 0.18 spelled it out), taken as they are.
const KEPT: [(&str, &str); 4] = [
    ("--model", "blocks"),
    ("--shift", "on"),
    ("--profile-strength", "inf"),
    ("--prior-strength", "inf"),
];

/// Whether `flag` given `value` is refused: a removed flag, but for the
/// values that ask for the model left (see [`KEPT`]).
pub(crate) fn refused(flag: &str, value: &str) -> bool {
    is_removed(flag) && !KEPT.contains(&(flag, value))
}

/// What a removed flag is told: every fit is by variational inference now.
pub(crate) const REMOVED_SAID: &str = "was for the MCMC fits (--model link and crp, \
    the sampler's sweeps and held-out reads), removed in faba 0.19: every fit is by \
    variational inference now, the site blockmodel (the old --model blocks) or the count \
    blockmodel (--modality count), set by --cell-roles, --gene-roles, --epochs, --tol \
    and --chains (restarts)";

/// The removed flags `args` give, refused before anything is read or
/// written: at the start of `faba cluster`, `--within` included.
pub(super) fn refuse_removed(args: &ClusterArgs) -> anyhow::Result<()> {
    let given = args.removed.given();
    anyhow::ensure!(given.is_empty(), "{} {REMOVED_SAID}", given.join(", "));
    Ok(())
}

/// The site blockmodel `args` fit: cell roles and gene roles, a shift per
/// block and each cell's own level, the role weights' global weights
/// uniform (as the variational fit keeps them).
pub(super) fn blockmodel(args: &ClusterArgs) -> settings::Settings {
    settings::Settings {
        cell_roles: args.cell_roles,
        gene_roles: args.gene_roles.unwrap_or(GENE_ROLES),
        alpha: args.membership_alpha,
        gene_alpha: args.gene_alpha,
        weights: settings::Weights::Symmetric,
        shift_sd: Some(args.shift_sd),
        strength: f64::INFINITY,
        cell_sd: args.cell_sd,
    }
}

/// Gene roles unless `--gene-roles` says.
pub(super) const GENE_ROLES: usize = 8;

/// The blockmodel `args` fit, as a run's record says it: the blocks preset
/// (as faba before 0.19 named it) unless `--gene-roles` changes it.
pub(super) fn recorded(args: &ClusterArgs) -> settings::Recorded {
    settings::Recorded {
        preset: args
            .gene_roles
            .is_none_or(|b| b == GENE_ROLES)
            .then(|| "blocks".to_string()),
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

/// Settings that would make the fit meaningless, refused up front.
pub(super) fn check(args: &ClusterArgs) -> anyhow::Result<()> {
    anyhow::ensure!(args.chains >= 1, "--chains must be at least 1");
    let positive = |x: f64| x > 0.0 && x.is_finite();
    anyhow::ensure!(
        args.epochs >= 1 && args.batch_size >= 1,
        "--epochs and --batch-size must be at least 1"
    );
    anyhow::ensure!(
        args.tol >= 0.0 && args.tol.is_finite(),
        "--tol must be 0 or more"
    );
    if args.within.is_some() {
        anyhow::ensure!(
            !args.is_count(),
            "--within fits sites within a count fit's cell modules: give --modality m6a or atoi, not count"
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
        anyhow::ensure!(
            args.removed.given_any().is_empty(),
            "{} chose a site fit's model; --modality count fits the count blockmodel",
            args.removed.given_any().join(", ")
        );
        return Ok(());
    }
    blockmodel(args).check()?;
    anyhow::ensure!(positive(args.gene_alpha), "--gene-alpha must be positive");
    anyhow::ensure!(positive(args.shift_sd), "--shift-sd must be positive");
    anyhow::ensure!(positive(args.cell_sd), "--cell-sd must be positive");
    anyhow::ensure!(
        positive(args.membership_alpha),
        "--membership-alpha must be positive"
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
