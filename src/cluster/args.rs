//! `faba cluster`'s options: what they are, which settings are refused, and
//! the input matrices and output directory they name.

use super::read;
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
        value_name = "m6a|atoi",
        help = "Which sites to take from an output directory (default: m6a, else atoi)"
    )]
    pub(super) modality: Option<Box<str>>,

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
        default_value_t = 5.0,
        help = "Prior strength: pseudo-reads pulling a cluster's site fraction toward the site's pooled fraction"
    )]
    pub(super) prior_strength: f64,

    #[arg(
        long,
        value_enum,
        default_value = "link",
        help = "link: each (cell, gene) link in a community, cells with memberships; blocks: cell roles and gene roles (gene modules) with a methylation shift per pair; crp: each cell in one cluster"
    )]
    pub(super) model: ModelKind,

    #[arg(long, default_value_t = 8, help = "Communities, with --model link")]
    pub(super) communities: usize,

    #[arg(
        long,
        default_value_t = 0.5,
        help = "With --model link or blocks: concentration of each cell's membership (smaller: cells condense onto fewer communities or roles)"
    )]
    pub(super) membership_alpha: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "With --model link: concentration of the communities' global weights"
    )]
    pub(super) community_gamma: f64,

    #[arg(long, default_value_t = 8, help = "Cell roles, with --model blocks")]
    pub(super) cell_roles: usize,

    #[arg(
        long,
        default_value_t = 8,
        help = "Gene roles (gene modules), with --model blocks"
    )]
    pub(super) gene_roles: usize,

    #[arg(
        long,
        default_value_t = 0.5,
        help = "With --model blocks: concentration of each gene's role weights"
    )]
    pub(super) gene_alpha: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "With --model blocks: prior sd of each (cell role, gene role) block's methylation shift, on the logit scale"
    )]
    pub(super) shift_sd: f64,

    #[arg(
        long,
        default_value_t = 1.0,
        help = "With --model blocks: prior sd of each cell's own methylation level, on the logit scale"
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
        help = "Independent chains, run in parallel"
    )]
    pub(super) chains: usize,

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

    /// The options given, for the command saved beside the run.
    #[arg(skip)]
    #[serde(skip)]
    pub(super) given: super::command::Given,
}

/// Which mixture `faba cluster` fits.
#[derive(clap::ValueEnum, Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum ModelKind {
    Crp,
    Link,
    Blocks,
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
                let sites = crate::site_analysis::output_dir::output_site_matrices(
                    &dir,
                    args.modality.as_deref(),
                )?;
                for (batch, path) in sites.matrices {
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

/// `-o`, else the first input directory, else the first matrix's directory.
pub(super) fn output_dir(args: &ClusterArgs) -> PathBuf {
    if let Some(o) = &args.output {
        return PathBuf::from(o.as_ref());
    }
    crate::run_record::default_output_dir(&args.data_files[0])
}

/// Settings that would make the fit meaningless, refused up front.
pub(super) fn check(args: &ClusterArgs) -> anyhow::Result<()> {
    anyhow::ensure!(
        args.prior_strength > 0.0 && args.prior_strength.is_finite(),
        "--prior-strength must be positive"
    );
    anyhow::ensure!(args.chains >= 1, "--chains must be at least 1");
    // Split R-hat needs at least two draws in each half of every chain.
    anyhow::ensure!(
        args.sweeps.saturating_sub(args.burn_in) / args.thin.max(1) >= 4,
        "--sweeps less --burn-in must leave at least 4 draws per chain (after --thin)"
    );
    // Labels are stored as u16.
    anyhow::ensure!(
        (1..u16::MAX as usize).contains(&args.communities),
        "--communities must be between 1 and {}",
        u16::MAX - 1
    );
    anyhow::ensure!(
        (1..u16::MAX as usize).contains(&args.cell_roles)
            && (1..u16::MAX as usize).contains(&args.gene_roles),
        "--cell-roles and --gene-roles must be between 1 and {}",
        u16::MAX - 1
    );
    let positive = |x: f64| x > 0.0 && x.is_finite();
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
    anyhow::ensure!(
        (0.0..1.0).contains(&args.holdout),
        "--holdout must be at least 0 and less than 1"
    );
    Ok(())
}

/// `args` with its inputs, asked for in a browser when none are named.
pub(super) fn with_inputs(
    args: &ClusterArgs,
    command: &'static str,
    batch: bool,
) -> anyhow::Result<Option<ClusterArgs>> {
    let mut args = args.clone();
    if args.data_files.is_empty() {
        let ask = crate::site_analysis::input_picker::ask_inputs(
            command,
            batch,
            "site matrices",
            crate::qc::layout::is_site_matrix_name,
            false,
        );
        let Some(chosen) = ask? else { return Ok(None) };
        args.data_files = chosen.paths;
    }
    Ok(Some(args))
}
