//! The reads under a fit, cell by cell: every fitted cell's converted and
//! all reads at each gene (or site) of the gene modules, read again from
//! the matrices the fit's run record names (or `--raw-data`), for the raw
//! view (`--raw`): the modules grid's layout, each block the heatmap of its
//! cells against its genes.

use super::data::{fitted_rows, run};
use super::Modules;
use crate::cluster::args::ClusterArgs;
use crate::cluster::read;
use crate::common::{column, read_record_batch};
use crate::qc::arms::Rgb;
use anyhow::Context;
use arrow::array::{Array, Float32Array, Int32Array, StringArray};
use rustc_hash::FxHashMap;
use std::path::Path;

/// What a raw view's column is: a gene (its sites' reads summed) or a site.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
pub enum RawUnit {
    #[default]
    Gene,
    Site,
}

impl RawUnit {
    /// What one is called.
    pub fn noun(self) -> &'static str {
        match self {
            RawUnit::Gene => "gene",
            RawUnit::Site => "site",
        }
    }
}

/// What colours a raw view's pixel, of the reads of the cells and the
/// features it pools.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
pub enum RawValue {
    /// Converted reads of all reads.
    #[default]
    Fraction,
    /// Converted reads per cell and feature.
    Methylated,
    /// log(1 + reads per cell and feature).
    Log1p,
}

/// How a raw view's block orders its cells and features, which are
/// exchangeable within it.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
pub enum RawSort {
    /// Each by its dots in the block (its pairs at `--raw-threshold` or
    /// above), most first.
    #[default]
    Value,
    /// Each by its reads in the block, most first.
    Reads,
    /// Cells by membership, strongest first; features in their module's
    /// order; a cell on one row across its band.
    None,
}

/// How the raw view's modules share the plot.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, clap::ValueEnum, serde::Serialize)]
pub enum RawBands {
    /// Every module as wide (and as tall) as the grid's blocks.
    Equal,
    /// Each module as its cells (or features) are many, within
    /// `--raw-band-min` and `--raw-band-max` of an equal band.
    #[default]
    Size,
}

/// A band sized by its module's size is no less than the first and no
/// more than the second of an equal band, unless `--raw-band-min` and
/// `--raw-band-max` say.
const BAND_CAP: (f64, f64) = (0.5, 2.0);

/// A dot's side in points, unless `--raw-dot` says.
const DOT: f64 = 1.0;

/// A pair is a dot at this beta or above, unless `--raw-threshold` says,
/// with this many converted reads, unless `--raw-min-methylated` says.
const THRESHOLD: f64 = 0.5;
const MIN_METHYLATED: u32 = 2;

/// A pixel under this many reads is blank, unless `--raw-min-reads` says.
const MIN_READS: u32 = 3;

/// The raw view's options.
#[derive(clap::Args, Debug, Clone, PartialEq)]
pub struct RawArgs {
    #[arg(
        long = "raw",
        help = "Draw the reads under the fit in place of the modules grid (--figure, --screen; in the browser, r on an image terminal): the grid's layout, each block its cells against its genes (or sites, --raw-unit), rasterised: solid dots where a cell's reads at a gene are methylated (--raw-threshold), sorted by them in each block (--raw-sort), or shades (--raw-ramp); its modules' names and lines drawn over it"
    )]
    pub raw: bool,

    #[arg(
        long = "raw-unit",
        value_enum,
        default_value_t,
        help = "The raw view's columns: genes (their sites' reads summed) or sites, each in its gene's main gene module"
    )]
    pub unit: RawUnit,

    #[arg(
        long = "raw-value",
        value_enum,
        default_value_t,
        help = "What colours the raw view: the methylated fraction (converted of all reads), the converted reads per cell and feature, or log(1 + reads) per cell and feature; each pixel pools the cells and features it covers"
    )]
    pub value: RawValue,

    #[arg(
        long = "raw-min-reads",
        value_name = "N",
        default_value_t = MIN_READS,
        value_parser = clap::value_parser!(u32).range(1..),
        help = "Leave a raw view's pixel blank (white) under N reads, so that a cell with a read or two is not drawn at 0 or 100%"
    )]
    pub min_reads: u32,

    #[arg(
        long = "raw-threshold",
        default_value_t = THRESHOLD,
        value_parser = beta,
        value_name = "BETA",
        help = "A cell and feature pair is a dot where its beta (converted of all reads) is at least BETA, as 0.5 for half methylated, with --raw-min-methylated converted reads; a drawn dot is solid where a pair under it is one, else white: no shades, which small dots do not show"
    )]
    pub threshold: f64,

    #[arg(
        long = "raw-min-methylated",
        value_name = "N",
        default_value_t = MIN_METHYLATED,
        help = "A cell and feature pair is a dot only with N converted reads at least, besides its beta at --raw-threshold"
    )]
    pub min_methylated: u32,

    #[arg(
        long = "raw-dot",
        value_name = "POINTS",
        default_value_t = DOT,
        value_parser = positive,
        help = "A raw view's dot, its side in points, the same across and down: each dot pools the cells and features under it, fewer and larger dots as it grows"
    )]
    pub dot: f64,

    #[arg(
        long = "raw-ramp",
        help = "Colour a raw view's pixels on a ramp from light to dark by their value, its top the 99th percentile (or --grid-cap's), instead of dots at --raw-threshold"
    )]
    pub ramp: bool,

    #[arg(
        long = "raw-colour",
        visible_alias = "raw-color",
        value_name = "COLOUR",
        value_parser = dot_colour,
        help = "A raw view's dots' colour, #rrggbb, #rgb or a name (white, black, gray, red, blue, green, orange, purple, yellow), e.g. black or #d62728; with --raw-ramp, the ramp's dark end, from a tint of it (default: dark blue, or near-black with --palette gray)"
    )]
    pub colour: Option<Rgb>,

    #[arg(
        long = "raw-sort",
        value_enum,
        default_value_t,
        help = "How each raw view block orders its cells and features, exchangeable within it: by their dots in the block (cell and feature pairs at --raw-threshold or above, with --raw-min-reads), most first, by their reads there, or none (cells by membership, features in their module's order, a cell then on one row across its band)"
    )]
    pub sort: RawSort,

    #[arg(
        long = "raw-bands",
        value_enum,
        default_value_t,
        help = "How the raw view's modules share the plot: size (each module as its cells or features are many, the data's own proportions, within --raw-band-min and --raw-band-max of an equal band), or equal (the grid's blocks); beside pinned genes on a screen they are equal, the genes' tracks on equal rows"
    )]
    pub bands: RawBands,

    #[arg(
        long = "raw-band-min",
        value_name = "TIMES",
        default_value_t = BAND_CAP.0,
        value_parser = positive,
        help = "With --raw-bands size, each band no less than TIMES an equal band (at most 1), so that a small module stays readable"
    )]
    pub band_min: f64,

    #[arg(
        long = "raw-band-max",
        value_name = "TIMES",
        default_value_t = BAND_CAP.1,
        value_parser = positive,
        help = "With --raw-bands size, each band no more than TIMES an equal band (at least 1), so that a large module does not crowd the rest"
    )]
    pub band_max: f64,

    #[arg(
        long = "raw-data",
        value_name = "MATRIX|DIR",
        num_args = 1..,
        help = "The matrices the fit was fitted on, as `faba cluster` was given them, when they are no longer where its run record says"
    )]
    pub data: Vec<Box<str>>,
}

impl RawArgs {
    /// The least and the most a band sized by its module is, of an equal
    /// band: never more than one, never less.
    pub fn band_cap(&self) -> (f64, f64) {
        (self.band_min.min(1.0), self.band_max.max(1.0))
    }
}

impl Default for RawArgs {
    fn default() -> Self {
        Self {
            raw: false,
            unit: RawUnit::Gene,
            value: RawValue::Fraction,
            min_reads: MIN_READS,
            threshold: THRESHOLD,
            min_methylated: MIN_METHYLATED,
            dot: DOT,
            ramp: false,
            colour: None,
            sort: RawSort::Value,
            bands: RawBands::Size,
            band_min: BAND_CAP.0,
            band_max: BAND_CAP.1,
            data: Vec::new(),
        }
    }
}

/// A colour, as `--raw-colour` names it.
fn dot_colour(s: &str) -> Result<Rgb, String> {
    crate::qc::arms::colour(s).map_err(|e| e.to_string())
}

/// A beta, 0 to 1.
fn beta(s: &str) -> Result<f64, String> {
    ranged(
        s,
        |v| (0.0..=1.0).contains(&v),
        "a beta from 0 to 1, as 0.5 for half methylated",
    )
}

/// A number above 0.
fn positive(s: &str) -> Result<f64, String> {
    ranged(s, |v| v > 0.0 && v.is_finite(), "a number above 0")
}

/// `s` as a number `ok` takes, else what to give, `said`.
fn ranged(s: &str, ok: impl Fn(f64) -> bool, said: &str) -> Result<f64, String> {
    match s.trim().parse::<f64>() {
        Ok(v) if ok(v) => Ok(v),
        Ok(_) => Err(format!("give {said}")),
        Err(e) => Err(format!("{s}: {e}")),
    }
}

/// One cell's reads: `(feature, converted, all)` by feature.
type CellReads = Vec<(u32, u32, u32)>;

/// The reads under a fit, cell by cell, its features numbered gene role by
/// gene role.
#[derive(Debug)]
pub struct Raw {
    pub unit: RawUnit,
    /// Per cell role (not column), its fitted cells' reads, the strongest
    /// member first.
    cells: Vec<Vec<CellReads>>,
    /// Per gene role, its first feature; one more, past the last.
    starts: Vec<u32>,
}

/// A pooled pixel: converted reads, all reads, and the (cell, feature)
/// pairs it covers.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Pooled {
    pub on: f64,
    pub total: f64,
    pub pairs: f64,
    pub dots: f64,
}

impl Pooled {
    /// Whether it is drawn as a dot: a pair under it is one.
    pub fn dot(self) -> bool {
        self.dots > 0.0
    }

    /// One pair's reads, and whether it is a dot as `shown` says.
    fn pair(on: u32, total: u32, shown: &RawArgs) -> Self {
        let mut p = Pooled {
            on: on as f64,
            total: total as f64,
            pairs: 1.0,
            ..Pooled::default()
        };
        // With the least reads, a dot at the threshold's beta and the
        // least methylated reads.
        if let Some(beta) = p.value(RawValue::Fraction, shown.min_reads) {
            let dot = beta >= shown.threshold && on >= shown.min_methylated;
            p.dots = f64::from(u8::from(dot));
        }
        p
    }

    /// `other`'s reads and dots added.
    fn add(&mut self, other: Pooled) {
        self.on += other.on;
        self.total += other.total;
        self.dots += other.dots;
    }

    /// Its value under `value`; none under `min_reads` reads (one at
    /// least).
    pub fn value(self, value: RawValue, min_reads: u32) -> Option<f64> {
        if self.total < min_reads.max(1) as f64 {
            return None;
        }
        Some(match value {
            RawValue::Fraction => self.on / self.total,
            RawValue::Methylated => self.on / self.pairs.max(1.0),
            RawValue::Log1p => (self.total / self.pairs.max(1.0)).ln_1p(),
        })
    }
}

/// One block's pixels, `rows` by `cols`, row by row.
#[derive(Debug, Clone, PartialEq)]
pub struct Tile {
    pub rows: usize,
    pub cols: usize,
    pub pixels: Vec<Pooled>,
}

impl Raw {
    /// Cell role `c`'s fitted cells.
    pub fn cells(&self, c: usize) -> usize {
        self.cells[c].len()
    }

    /// Gene role `g`'s features.
    pub fn features(&self, g: usize) -> usize {
        (self.starts[g + 1] - self.starts[g]) as usize
    }

    /// The tile of cell role `c` against gene role `g`, at most `rows`
    /// pixels for its cells and `cols` for its features (none fewer than
    /// one, none more than it has): each pixel pooling the cells and the
    /// features it covers, both in the order `shown` sorts them by.
    pub fn tile(&self, c: usize, g: usize, (rows, cols): (usize, usize), shown: &RawArgs) -> Tile {
        let (n, k) = (self.cells[c].len(), self.features(g));
        let (rows, cols) = (rows.clamp(1, n.max(1)), cols.clamp(1, k.max(1)));
        let (from, to) = (self.starts[g], self.starts[g + 1]);
        // Each cell's reads in the block: role `g`'s, one run, as a
        // cell's reads are by feature.
        let block: Vec<&[(u32, u32, u32)]> = (self.cells[c].iter())
            .map(|r| &r[r.partition_point(|e| e.0 < from)..r.partition_point(|e| e.0 < to)])
            .collect();
        let (cell_at, feature_at) = self.places(&block, (from, k), shown);
        let mut pixels = vec![Pooled::default(); rows * cols];
        // Each pixel's pairs: the cells of its row by the features of its
        // column, as an even split puts them.
        let share = |i: usize, of: usize, px: usize| (i * px / of.max(1)).min(px - 1);
        let mut row_cells = vec![0usize; rows];
        (0..n).for_each(|i| row_cells[share(i, n, rows)] += 1);
        let mut col_features = vec![0usize; cols];
        (0..k).for_each(|j| col_features[share(j, k, cols)] += 1);
        for (r, &nc) in row_cells.iter().enumerate() {
            for (q, &nf) in col_features.iter().enumerate() {
                pixels[r * cols + q].pairs = (nc * nf) as f64;
            }
        }
        for (i, reads) in block.iter().enumerate() {
            let r = share(cell_at[i], n, rows);
            for &(f, on, total) in reads.iter() {
                let q = share(feature_at[(f - from) as usize], k, cols);
                pixels[r * cols + q].add(Pooled::pair(on, total, shown));
            }
        }
        Tile { rows, cols, pixels }
    }

    /// Each of a block's cells' and features' place in it, as `shown`
    /// sorts them: as they are (`none`), by their dots in the block
    /// (`value`), or by their reads there (`reads`), most first, a tie as
    /// they are.
    fn places(
        &self,
        block: &[&[(u32, u32, u32)]],
        (from, k): (u32, usize),
        shown: &RawArgs,
    ) -> (Vec<usize>, Vec<usize>) {
        let n = block.len();
        if shown.sort == RawSort::None {
            return ((0..n).collect(), (0..k).collect());
        }
        let mut cells = vec![Pooled::default(); n];
        let mut features = vec![Pooled::default(); k];
        for (i, reads) in block.iter().enumerate() {
            for &(f, on, total) in reads.iter() {
                let p = Pooled::pair(on, total, shown);
                cells[i].add(p);
                features[(f - from) as usize].add(p);
            }
        }
        let key = |p: &Pooled| match shown.sort {
            RawSort::Reads => p.total,
            _ => p.dots,
        };
        // Most first; a tie as they are.
        let places = |sums: &[Pooled]| super::order::rank(sums.iter().map(|p| -key(p)).collect());
        (places(&cells), places(&features))
    }
}

/// A site: its gene, chromosome and position.
type SiteKey = (Box<str>, Box<str>, i64);

/// What reading the reads under a fit needs of its modules, owned, so
/// that they are read on a thread of their own: its modality, its sites'
/// features by `unit`, numbered gene role by gene role in each role's
/// order (genes with the most sites first, a gene's sites by position),
/// each role's first feature, and its cell roles.
#[derive(Debug, Clone)]
pub struct Plan {
    unit: RawUnit,
    modality: String,
    /// Each site's feature, by gene, chromosome and position.
    by_site: Vec<(SiteKey, u32)>,
    starts: Vec<u32>,
    cell_roles: usize,
}

impl Plan {
    /// Refused for a count fit, whose cells' counts are not reads at sites.
    pub fn of(m: &Modules, unit: RawUnit) -> anyhow::Result<Self> {
        anyhow::ensure!(
            !m.per_cell(),
            "the raw view draws a site fit's reads: a count fit has none"
        );
        let mut by_site = Vec::new();
        let mut starts = vec![0u32];
        let mut next = 0u32;
        for g in 0..m.gene_roles {
            for gene in m.genes(g) {
                for i in gene.clone() {
                    if i == gene.start || unit == RawUnit::Site {
                        next += 1;
                    }
                    let s = &m.sites[m.members[g][i].site];
                    by_site.push(((s.gene.clone(), s.chr.clone(), s.pos), next - 1));
                }
            }
            starts.push(next);
        }
        Ok(Self {
            unit,
            modality: m.modality.clone(),
            by_site,
            starts,
            cell_roles: m.cell_roles,
        })
    }

    /// The reads of the fit in `dir`: from `data` when given, else the
    /// matrices its run record names.
    pub fn read(self, dir: &Path, data: &[Box<str>]) -> anyhow::Result<Raw> {
        let inputs = fitted_inputs(dir, &self.modality, data)?;
        let read = read::read_inputs(&inputs, &crate::cluster::progress::Progress::new(1))?;
        // Each matrix site's feature, when it is one of the modules'.
        let index: FxHashMap<(&str, &str, i64), usize> = (read.sites.iter().enumerate())
            .map(|(i, s)| ((&*s.gene, &*s.chr, s.pos), i))
            .collect();
        let mut feature_of: Vec<Option<u32>> = vec![None; read.sites.len()];
        for ((gene, chr, pos), k) in &self.by_site {
            if let Some(&i) = index.get(&(&**gene, &**chr, *pos)) {
                feature_of[i] = Some(*k);
            }
        }
        let path = dir.join(format!("{}_cell_roles.parquet", self.modality));
        let table = read_record_batch(&path).with_context(|| path.display().to_string())?;
        let fitted = fitted_cells(dir, &self.modality, &table, self.cell_roles)?;
        anyhow::ensure!(
            !fitted.is_empty(),
            "{} has no fitted cell with a cell module: no cells to draw",
            path.display()
        );
        let rows: FxHashMap<(&str, &str), usize> = (read.cells.iter().enumerate())
            .map(|(i, c)| ((&*c.batch, &*c.barcode), i))
            .collect();
        let mut batches: Vec<&str> = read.cells.iter().map(|c| &*c.batch).collect();
        batches.sort_unstable();
        batches.dedup();
        let mut cells: Vec<Vec<(f32, CellReads)>> = vec![Vec::new(); self.cell_roles];
        let mut missed = 0usize;
        for f in &fitted {
            let batch = matched_batch(f.batch, &batches);
            let Some(&i) = batch.and_then(|b| rows.get(&(b, f.barcode))) else {
                missed += 1;
                continue;
            };
            let (sites, on, off) = read.reads.of(i);
            let mut reads: CellReads = (0..sites.len())
                .filter_map(|e| Some((feature_of[sites[e] as usize]?, on[e], on[e] + off[e])))
                .collect();
            // A gene's sites summed: by feature, its entries merged.
            reads.sort_unstable_by_key(|e| e.0);
            reads.dedup_by(|e, kept| {
                let same = e.0 == kept.0;
                if same {
                    (kept.1, kept.2) = (kept.1 + e.1, kept.2 + e.2);
                }
                same
            });
            cells[f.role].push((f.weight, reads));
        }
        anyhow::ensure!(
            missed < fitted.len(),
            "none of the fit's cells is in {}: name the matrices it was fitted on with --raw-data",
            inputs
                .iter()
                .map(|i| &*i.path)
                .collect::<Vec<_>>()
                .join(", ")
        );
        if missed > 0 {
            log::warn!(
                "{missed} of the fit's cells are not in its matrices: left out of the raw view"
            );
        }
        let cells = (cells.into_iter())
            .map(|mut role| {
                role.sort_by(|a, b| b.0.total_cmp(&a.0));
                role.into_iter().map(|(_, r)| r).collect()
            })
            .collect();
        Ok(Raw {
            unit: self.unit,
            cells,
            starts: self.starts,
        })
    }
}

/// A fitted cell of the cell table: its batch and barcode, main role and
/// that role's weight.
struct Fitted<'a> {
    batch: &'a str,
    barcode: &'a str,
    role: usize,
    weight: f32,
}

/// The fitted cells of `table` (`{modality}_cell_roles` in `dir`), as the
/// grid's sizes count them, each with its main role of `k`.
fn fitted_cells<'a>(
    dir: &Path,
    modality: &str,
    table: &'a arrow::record_batch::RecordBatch,
    k: usize,
) -> anyhow::Result<Vec<Fitted<'a>>> {
    let batch: &StringArray = column(table, "batch")?;
    let barcode: &StringArray = column(table, "barcode")?;
    let role: &Int32Array = column(table, "role")?;
    let weight: &Float32Array = column(table, "weight")?;
    Ok((fitted_rows(dir, modality, table, k)?.into_iter())
        .filter(|&i| role.is_valid(i))
        .filter_map(|i| {
            let r = usize::try_from(role.value(i)).ok()?.checked_sub(1)?;
            (r < k).then(|| Fitted {
                batch: batch.value(i),
                barcode: barcode.value(i),
                role: r,
                weight: weight.value(i),
            })
        })
        .collect())
}

/// The batch of the matrices read that the cell table's `batch` is: the
/// same name; else the one it qualifies (`faba cluster` names a batch two
/// folders have `{folder}/{batch}` with the folder as it was given, which
/// the record holds made absolute); else the one named as its last part.
/// None when two would be.
fn matched_batch<'a>(batch: &str, read: &[&'a str]) -> Option<&'a str> {
    if let Some(b) = read.iter().find(|b| **b == batch) {
        return Some(b);
    }
    let found: Vec<&str> = (read.iter().copied())
        .filter(|b| crate::cluster::within::qualifies(b, batch))
        .collect();
    let found = match found.is_empty() {
        false => found,
        true => {
            let plain = batch.rsplit('/').next()?;
            (read.iter().copied())
                .filter(|b| b.rsplit('/').next() == Some(plain))
                .collect()
        }
    };
    match found.as_slice() {
        [b] => Some(b),
        _ => None,
    }
}

/// The matrices `faba cluster` read for the fit in `dir` (of `modality`):
/// `data` as it would take them when given; else as its run record lists
/// them, each with its batch (`results.inputs`), or in a record before
/// that, as it was given them (`inputs.matrices`, made absolute).
fn fitted_inputs(
    dir: &Path,
    modality: &str,
    data: &[Box<str>],
) -> anyhow::Result<Vec<read::Input>> {
    let (name, record) = run::record(dir, modality).map_err(anyhow::Error::msg)?;
    let listed = |at: &str| -> Vec<&serde_json::Value> {
        (record
            .pointer(at)
            .and_then(|v| v.as_array())
            .into_iter()
            .flatten())
        .collect()
    };
    let inputs = match (data.is_empty(), listed("/results/inputs")) {
        (true, read) if !read.is_empty() => (read.into_iter())
            .filter_map(|i| {
                Some(read::Input {
                    path: there(i)?.into(),
                    folder: "".into(),
                    batch: i["batch"].as_str()?.into(),
                    fit: true,
                })
            })
            .collect(),
        (given, _) => {
            let data_files = match given {
                true => (listed("/inputs/matrices").into_iter())
                    .filter_map(|m| there(m).map(Into::into))
                    .collect(),
                false => data.to_vec(),
            };
            anyhow::ensure!(
                !data_files.is_empty(),
                "{name} lists no matrices: name them with --raw-data"
            );
            let mut args = cluster_args(&["--modality".into(), modality.into()])?;
            args.data_files = data_files;
            crate::cluster::args::inputs(&args)?
        }
    };
    let missing: Vec<&str> = (inputs.iter())
        .map(|i| &*i.path)
        .filter(|p| !Path::new(p).exists())
        .collect();
    anyhow::ensure!(
        missing.is_empty(),
        "the fit's matrices are not where its run record says ({}): name them with --raw-data",
        missing.join(", ")
    );
    Ok(inputs)
}

/// A recorded input's path: as given, else the file it resolved to when
/// only that is there still.
fn there(entry: &serde_json::Value) -> Option<&str> {
    let path = entry["path"].as_str()?;
    let resolved = entry["resolved"].as_str();
    match Path::new(path).exists() {
        true => Some(path),
        false => resolved.filter(|r| Path::new(r).exists()).or(Some(path)),
    }
}

/// `faba cluster`'s options from `words` after it.
fn cluster_args(words: &[String]) -> anyhow::Result<ClusterArgs> {
    #[derive(clap::Parser)]
    struct Again {
        #[command(flatten)]
        args: ClusterArgs,
    }
    let words = std::iter::once("cluster").chain(words.iter().map(String::as_str));
    Ok(<Again as clap::Parser>::try_parse_from(words)?.args)
}

#[cfg(test)]
#[path = "tests/raw.rs"]
mod tests;
