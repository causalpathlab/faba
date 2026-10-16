//! The site matrices' reads, kept per cell over its covered sites, for clustering.
//!
//! Nothing here is dense over cells x sites: each cell holds only the sites it
//! covers, as `(site, converted, unconverted)` entries in one flat array.

use crate::site_analysis::pileup::parse_row_channel;
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
}

/// Sites, cells and their reads.
pub struct SiteData {
    pub modality: Box<str>,
    pub sites: Vec<Site>,
    pub cells: Vec<Cell>,
    pub reads: Reads,
}

/// Read every input's site rows (both channels) of the first modality met.
pub fn read_inputs(inputs: &[Input]) -> anyhow::Result<SiteData> {
    let mut modality: Option<Box<str>> = None;
    let mut sites: Vec<Site> = Vec::new();
    let mut index: FxHashMap<(Box<str>, Box<str>, i64), u32> = FxHashMap::default();
    let mut cells: Vec<Cell> = Vec::new();
    let mut per_cell: Vec<Vec<(u32, u32, u32)>> = Vec::new();
    for input in inputs {
        let (backend, resolved) = resolve_backend_file(&input.path, None)?;
        let mut data = open_sparse_matrix(&resolved, &backend)?;
        // Every column is read once: preload them when they fit the budget.
        let nnz = data.num_non_zeros().unwrap_or(0);
        let _reservation = data_beans::sparse_io::reserve_preload(nnz, &input.path);
        if _reservation.is_some() {
            data.preload_columns()?;
        }
        // Each row: its site, and whether it is the converted channel.
        let mut rows: Vec<Option<(u32, bool)>> = Vec::new();
        for name in data.row_names()? {
            let Some((gene, m, chr, pos, channel)) = parse_row_channel(&name) else {
                rows.push(None);
                continue;
            };
            let (Some(converted), false) = (channel, chr.is_empty()) else {
                rows.push(None);
                continue;
            };
            if modality.get_or_insert_with(|| m.into()).as_ref() != m {
                rows.push(None);
                continue;
            }
            let key = (Box::<str>::from(gene), Box::<str>::from(chr), pos);
            let next = sites.len() as u32;
            let s = *index.entry(key).or_insert_with(|| {
                sites.push(Site {
                    gene: gene.into(),
                    chr: chr.into(),
                    pos,
                });
                next
            });
            rows.push(Some((s, converted)));
        }
        let barcodes = data.column_names()?;
        let first = cells.len();
        cells.extend(barcodes.iter().map(|b| Cell {
            barcode: b.clone(),
            batch: input.batch.clone(),
            fit: input.fit,
        }));
        per_cell.resize_with(cells.len(), Vec::new);
        let ncol = barcodes.len();
        for start in (0..ncol).step_by(COLUMN_BLOCK) {
            let block: Vec<usize> = (start..(start + COLUMN_BLOCK).min(ncol)).collect();
            let (_, _, triplets) = data.read_triplets_by_columns(block)?;
            for (row, col, value) in triplets {
                let Some(Some((s, converted))) = rows.get(row as usize).copied() else {
                    continue;
                };
                let n = value.round().max(0.0) as u32;
                if n == 0 {
                    continue;
                }
                // The two channels at a site merge in `Reads::from_cells`.
                let entry = if converted { (s, n, 0) } else { (s, 0, n) };
                per_cell[first + start + col as usize].push(entry);
            }
        }
        info!("{}: {} cells", input.path, ncol);
    }
    let modality = modality.ok_or_else(|| anyhow::anyhow!("no site rows in the inputs"))?;
    Ok(SiteData {
        modality,
        sites,
        cells,
        reads: Reads::from_cells(per_cell),
    })
}

/// Who survived [`filter`]: kept cells, as indices into the input.
pub struct Kept {
    pub cells: Vec<usize>,
}

/// Drop cells covering fewer than `min_cell_sites` sites, then sites covered by
/// fewer than `min_site_cells` of the remaining fitted cells; renumber the
/// sites of `data` to the kept ones.
pub fn filter(data: &SiteData, min_cell_sites: usize, min_site_cells: usize) -> (SiteData, Kept) {
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
    let per_cell = kept_cells
        .iter()
        .map(|&c| {
            let (s, on, off) = data.reads.of(c);
            (0..s.len())
                .filter(|&i| renumber[s[i] as usize] != u32::MAX)
                .map(|i| (renumber[s[i] as usize], on[i], off[i]))
                .collect()
        })
        .collect();
    let out = SiteData {
        modality: data.modality.clone(),
        sites: kept_sites.iter().map(|&s| data.sites[s].clone()).collect(),
        cells: kept_cells.iter().map(|&c| data.cells[c].clone()).collect(),
        reads: Reads::from_cells(per_cell),
    };
    (out, Kept { cells: kept_cells })
}
