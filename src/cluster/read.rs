//! The site matrices' reads, kept per cell over its covered sites, for clustering.
//!
//! Nothing here is dense over cells x sites: each cell holds only the sites it
//! covers, as `(site, converted, unconverted)` entries in one flat array.

use crate::site_analysis::site_io::parse_row_channel;
use data_beans::hdf5_io::resolve_backend_file;
use data_beans::sparse_io::open_sparse_matrix;
use log::info;
use rustc_hash::FxHashMap;

/// Cells read from a matrix at a time.
const COLUMN_BLOCK: usize = 4096;

/// One input matrix: its path, its batch, and whether its cells are fitted
/// (the rest are assigned to the fitted clusters afterwards).
pub struct Input {
    pub path: Box<str>,
    /// The folder it came from: its output directory, else the one it sits in.
    pub folder: Box<str>,
    pub batch: Box<str>,
    pub fit: bool,
}

/// One site, shared by every matrix that has it.
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub gene: Box<str>,
    pub chr: Box<str>,
    pub pos: i64,
}

#[derive(Debug, Clone)]
pub struct Cell {
    pub barcode: Box<str>,
    pub batch: Box<str>,
    pub fit: bool,
}

/// Each cell's reads at the sites it covers: cell `c`'s entries are
/// `ptr[c]..ptr[c + 1]` of `site`, `on` (converted) and `off` (unconverted).
#[derive(Debug, Default, Clone)]
pub struct Reads {
    pub ptr: Vec<usize>,
    pub site: Vec<u32>,
    pub on: Vec<u32>,
    pub off: Vec<u32>,
}

impl Reads {
    pub fn cells(&self) -> usize {
        self.ptr.len().saturating_sub(1)
    }

    /// Cell `c`'s `(sites, converted, unconverted)`.
    pub fn of(&self, c: usize) -> (&[u32], &[u32], &[u32]) {
        let r = self.ptr[c]..self.ptr[c + 1];
        (&self.site[r.clone()], &self.on[r.clone()], &self.off[r])
    }

    /// From each cell's `(site, converted, unconverted)` entries, in any order.
    pub fn from_cells(cells: Vec<Vec<(u32, u32, u32)>>) -> Self {
        let mut reads = Reads {
            ptr: vec![0],
            ..Default::default()
        };
        for mut entries in cells {
            entries.sort_unstable_by_key(|e| e.0);
            let start = reads.site.len();
            for (s, on, off) in entries {
                // A site met twice (its two channels) adds up.
                if reads.site.len() > start && reads.site.last() == Some(&s) {
                    *reads.on.last_mut().unwrap() += on;
                    *reads.off.last_mut().unwrap() += off;
                    continue;
                }
                reads.site.push(s);
                reads.on.push(on);
                reads.off.push(off);
            }
            reads.ptr.push(reads.site.len());
        }
        reads
    }

    /// The reads of cells `cells`, in that order.
    pub fn subset(&self, cells: &[usize]) -> Self {
        let mut out = Reads {
            ptr: Vec::with_capacity(cells.len() + 1),
            ..Default::default()
        };
        out.ptr.push(0);
        for &c in cells {
            let r = self.ptr[c]..self.ptr[c + 1];
            out.site.extend_from_slice(&self.site[r.clone()]);
            out.on.extend_from_slice(&self.on[r.clone()]);
            out.off.extend_from_slice(&self.off[r]);
            out.ptr.push(out.site.len());
        }
        out
    }
}

/// Reads as a fit sees them: each entry's converted and unconverted counts,
/// weighted (importance weights that take batches out) or as read.
pub struct Weighted<'a> {
    pub reads: &'a Reads,
    pub on: Vec<f64>,
    pub off: Vec<f64>,
}

impl<'a> Weighted<'a> {
    /// `weights`: each entry's converted and unconverted reads' weights.
    pub fn new(reads: &'a Reads, weights: Option<(&[f64], &[f64])>) -> Self {
        let scale = |n: &[u32], w: Option<&[f64]>| -> Vec<f64> {
            match w {
                Some(w) => n.iter().zip(w).map(|(&n, &w)| n as f64 * w).collect(),
                None => n.iter().map(|&n| n as f64).collect(),
            }
        };
        Weighted {
            reads,
            on: scale(&reads.on, weights.map(|w| w.0)),
            off: scale(&reads.off, weights.map(|w| w.1)),
        }
    }

    pub fn cells(&self) -> usize {
        self.reads.cells()
    }

    /// Cell `c`'s `(sites, converted, unconverted)`.
    pub fn of(&self, c: usize) -> (&[u32], &[f64], &[f64]) {
        let r = self.reads.ptr[c]..self.reads.ptr[c + 1];
        (
            &self.reads.site[r.clone()],
            &self.on[r.clone()],
            &self.off[r],
        )
    }
}

/// Sites, cells and their reads.
pub struct SiteData {
    pub modality: Box<str>,
    pub sites: Vec<Site>,
    pub cells: Vec<Cell>,
    pub reads: Reads,
}

/// Each site's raw reads over every cell read, converted and not, and how
/// many cells have any there.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct SiteStats {
    pub on: Vec<u64>,
    pub off: Vec<u64>,
    pub cells: Vec<u64>,
}

impl SiteStats {
    /// `sites` sites with nothing counted.
    pub fn zeros(sites: usize) -> Self {
        Self {
            on: vec![0; sites],
            off: vec![0; sites],
            cells: vec![0; sites],
        }
    }

    pub fn of(reads: &Reads, sites: usize) -> Self {
        let mut out = Self::zeros(sites);
        for e in 0..reads.site.len() {
            let s = reads.site[e] as usize;
            out.on[s] += reads.on[e] as u64;
            out.off[s] += reads.off[e] as u64;
            out.cells[s] += u64::from(reads.on[e] + reads.off[e] > 0);
        }
        out
    }
}

/// Read every input's site rows (both channels) of the first modality met;
/// [`super::progress::Cancelled`] before the next input, or the next block
/// of an input's columns, once `progress` is.
pub fn read_inputs(
    inputs: &[Input],
    progress: &super::progress::Progress,
) -> anyhow::Result<SiteData> {
    let mut modality: Option<Box<str>> = None;
    let mut sites: Vec<Site> = Vec::new();
    let mut index: FxHashMap<(Box<str>, Box<str>, i64), u32> = FxHashMap::default();
    let mut cells: Vec<Cell> = Vec::new();
    let mut per_cell: Vec<Vec<(u32, u32, u32)>> = Vec::new();
    for input in inputs {
        let mut visit = SiteVisit {
            input,
            rows: Vec::new(),
            other: None,
            first: cells.len(),
            modality: &mut modality,
            sites: &mut sites,
            index: &mut index,
            cells: &mut cells,
            per_cell: &mut per_cell,
        };
        walk_matrix(&input.path, progress, &mut visit)?;
    }
    let modality = modality.ok_or_else(|| anyhow::anyhow!("no site rows in the inputs"))?;
    Ok(SiteData {
        modality,
        sites,
        cells,
        reads: Reads::from_cells(per_cell),
    })
}

/// What [`walk_matrix`] gives a matrix's reader: its row names, then its
/// column names (cells), then each entry above 0, rounded to a whole
/// count, as (row, column, count).
pub(super) trait Visit {
    fn rows(&mut self, names: &[Box<str>]) -> anyhow::Result<()>;
    fn columns(&mut self, barcodes: &[Box<str>]);
    fn entry(&mut self, row: usize, col: usize, n: u32);
}

/// One matrix at `path` read through to `visit`, columns a block at a
/// time, preloaded when they fit the budget.
pub(super) fn walk_matrix(
    path: &str,
    progress: &super::progress::Progress,
    visit: &mut impl Visit,
) -> anyhow::Result<()> {
    progress.check()?;
    let (backend, resolved) = resolve_backend_file(path, None)?;
    let mut data = open_sparse_matrix(&resolved, &backend)?;
    // Every column is read once: preload them when they fit the budget.
    let nnz = data.num_non_zeros().unwrap_or(0);
    let _reservation = data_beans::sparse_io::reserve_preload(nnz, path);
    if _reservation.is_some() {
        data.preload_columns()?;
    }
    visit.rows(&data.row_names()?)?;
    let barcodes = data.column_names()?;
    visit.columns(&barcodes);
    let ncol = barcodes.len();
    for start in (0..ncol).step_by(COLUMN_BLOCK) {
        progress.check()?;
        let block: Vec<usize> = (start..(start + COLUMN_BLOCK).min(ncol)).collect();
        let (_, _, triplets) = data.read_triplets_by_columns(block)?;
        for (row, col, value) in triplets {
            let n = value.round().max(0.0) as u32;
            if n > 0 {
                visit.entry(row as usize, start + col as usize, n);
            }
        }
    }
    info!("{path}: {ncol} cells");
    Ok(())
}

/// A site matrix's reader: its rows matched to the sites (and modality)
/// of the inputs before it, its cells and their entries appended.
struct SiteVisit<'v> {
    input: &'v Input,
    /// Each row: its site, and whether it is the converted channel.
    rows: Vec<Option<(u32, bool)>>,
    other: Option<Box<str>>,
    first: usize,
    modality: &'v mut Option<Box<str>>,
    sites: &'v mut Vec<Site>,
    index: &'v mut FxHashMap<(Box<str>, Box<str>, i64), u32>,
    cells: &'v mut Vec<Cell>,
    per_cell: &'v mut Vec<Vec<(u32, u32, u32)>>,
}

impl Visit for SiteVisit<'_> {
    fn rows(&mut self, names: &[Box<str>]) -> anyhow::Result<()> {
        for name in names {
            let Some((gene, m, chr, pos, channel)) = parse_row_channel(name) else {
                self.rows.push(None);
                continue;
            };
            let (Some(converted), false) = (channel, chr.is_empty()) else {
                self.rows.push(None);
                continue;
            };
            if self.modality.get_or_insert_with(|| m.into()).as_ref() != m {
                self.other.get_or_insert_with(|| m.into());
                self.rows.push(None);
                continue;
            }
            let key = (Box::<str>::from(gene), Box::<str>::from(chr), pos);
            let next = self.sites.len() as u32;
            let sites = &mut *self.sites;
            let s = *self.index.entry(key).or_insert_with(|| {
                sites.push(Site {
                    gene: gene.into(),
                    chr: chr.into(),
                    pos,
                });
                next
            });
            self.rows.push(Some((s, converted)));
        }
        // Inputs are clustered together only if they hold the same sites.
        if self.rows.iter().all(Option::is_none) {
            if let (Some(m), Some(o)) = (&*self.modality, &self.other) {
                anyhow::bail!(
                    "{} holds {o} sites, not {m} like the inputs before it",
                    self.input.path
                );
            }
        }
        Ok(())
    }

    fn columns(&mut self, barcodes: &[Box<str>]) {
        self.cells.extend(barcodes.iter().map(|b| Cell {
            barcode: b.clone(),
            batch: self.input.batch.clone(),
            fit: self.input.fit,
        }));
        self.per_cell.resize_with(self.cells.len(), Vec::new);
    }

    fn entry(&mut self, row: usize, col: usize, n: u32) {
        let Some(Some((s, converted))) = self.rows.get(row).copied() else {
            return;
        };
        // The two channels at a site merge in `Reads::from_cells`.
        let entry = if converted { (s, n, 0) } else { (s, 0, n) };
        self.per_cell[self.first + col].push(entry);
    }
}

/// Drop cells covering fewer than `min_cell_sites` sites, then sites covered by
/// fewer than `min_site_cells` of the remaining fitted cells, then cells left
/// with fewer than `min_cell_sites` of the kept sites; renumber the sites of
/// `data` to the kept ones. Also the kept cells' indices in `data`.
pub fn filter(
    data: &SiteData,
    min_cell_sites: usize,
    min_site_cells: usize,
) -> (SiteData, Vec<usize>) {
    let n = data.reads.cells();
    let kept_cells: Vec<usize> = (0..n)
        .filter(|&c| data.reads.ptr[c + 1] - data.reads.ptr[c] >= min_cell_sites)
        .collect();
    let mut site_cells = vec![0usize; data.sites.len()];
    for &c in kept_cells.iter().filter(|&&c| data.cells[c].fit) {
        for &s in data.reads.of(c).0 {
            site_cells[s as usize] += 1;
        }
    }
    let kept_sites: Vec<usize> = (0..data.sites.len())
        .filter(|&s| site_cells[s] >= min_site_cells)
        .collect();
    let mut renumber = vec![u32::MAX; data.sites.len()];
    for (i, &s) in kept_sites.iter().enumerate() {
        renumber[s] = i as u32;
    }
    let entries = |c: usize| -> Vec<(u32, u32, u32)> {
        let (s, on, off) = data.reads.of(c);
        (0..s.len())
            .filter(|&i| renumber[s[i] as usize] != u32::MAX)
            .map(|i| (renumber[s[i] as usize], on[i], off[i]))
            .collect()
    };
    // A cell whose sites went with the site filter has too little left.
    let (kept_cells, per_cell): (Vec<usize>, Vec<_>) = kept_cells
        .into_iter()
        .map(|c| (c, entries(c)))
        .filter(|(_, e)| e.len() >= min_cell_sites)
        .unzip();
    let out = SiteData {
        modality: data.modality.clone(),
        sites: kept_sites.iter().map(|&s| data.sites[s].clone()).collect(),
        cells: kept_cells.iter().map(|&c| data.cells[c].clone()).collect(),
        reads: Reads::from_cells(per_cell),
    };
    (out, kept_cells)
}
