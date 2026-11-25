//! `faba cluster --within TABLE`: the groups of a cell table, taken as
//! given, in place of a count fit's cell modules.
//!
//! The table is a Parquet file with a cell column and a group column (as
//! lupin's annotation export writes, `{out}.cell_annotation.parquet`, a
//! `cell` column and per-cell `cluster` and `label`), or a TSV or CSV of
//! cell and group, the header optional. `--within-column` names the group
//! column (default: `label`, else `cluster`, else the other of two) and
//! `--within-cell-column` the cell's (default: the first of
//! [`CELL_COLUMNS`] there).
//!
//! A cell's name is matched to the site matrices' cells as a count fit's
//! cells are (see [`Names::read`]): `{barcode}@{tag}` is barcode
//! `{barcode}` of the batch `{tag}` names, less a matrix kind (`b1_count`
//! is batch `b1`, as cells are named when several matrices are read
//! together), and that batch pairs with the sites' as a count fit's
//! batches do. A bare barcode is looked for in every site batch nothing
//! pairs with, and left out where it is in more than one.

use super::{ByBarcode, ByBatch, Group};
use crate::common::read_record_batch;
use anyhow::Context;
use arrow::array::{Array, StringArray};
use arrow::datatypes::DataType;
use rustc_hash::FxHashMap;
use std::path::Path;

/// The cell columns looked for when `--within-cell-column` names none, in
/// this order, case aside.
pub(in crate::cluster) const CELL_COLUMNS: [&str; 4] =
    ["cell", "barcode", "cell_barcode", "obs_names"];

/// The group columns looked for when `--within-column` names none.
const GROUP_COLUMNS: [&str; 2] = ["label", "cluster"];

/// Whether `--within` names a cell table rather than a count fit: a file
/// of a table's extension.
pub(in crate::cluster) fn is_table(named: &str) -> bool {
    let path = Path::new(named);
    let ext = path
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase());
    path.is_file() && matches!(ext.as_deref(), Some("parquet" | "tsv" | "csv" | "txt"))
}

/// A cell table's groups, its cells under them.
#[derive(Debug)]
pub(in crate::cluster) struct Names {
    /// Each group's name, as its directory is named: [`Group`] `g` is
    /// `names[g]`, in order.
    pub names: Vec<String>,
    /// The tagged cells, by batch then barcode.
    pub cells: ByBatch,
    /// The bare barcodes.
    pub loose: ByBarcode,
    /// The cells named in each group.
    pub named: Vec<usize>,
    /// The columns read, the cell's and the group's (none for a headless
    /// TSV).
    pub columns: Option<(String, String)>,
}

impl Names {
    /// The table at `path`, its cells in `group` (`--within-column`) by
    /// `cell` (`--within-cell-column`). A row without a group is left out.
    pub(in crate::cluster) fn read(
        path: &Path,
        group: Option<&str>,
        cell: Option<&str>,
    ) -> anyhow::Result<Self> {
        let parquet = (path.extension()).is_some_and(|e| e.eq_ignore_ascii_case("parquet"));
        let (rows, columns) = match parquet {
            true => parquet_rows(path, group, cell)?,
            false => text_rows(path, group, cell)?,
        };
        Self::of(rows, columns)
    }

    /// `rows` (each a cell's name and group) grouped.
    fn of(rows: Vec<(String, String)>, columns: Option<(String, String)>) -> anyhow::Result<Self> {
        let mut by_name: FxHashMap<String, String> = FxHashMap::default();
        let mut names: Vec<String> = Vec::new();
        for (_, given) in &rows {
            let name = sanitised(given);
            match by_name.get(&name) {
                Some(other) if other != given => anyhow::bail!(
                    "the groups `{other}` and `{given}` are both named {name} for a file name: \
                     rename one in the table"
                ),
                Some(_) => {}
                None => {
                    by_name.insert(name.clone(), given.clone());
                    names.push(name);
                }
            }
        }
        anyhow::ensure!(!names.is_empty(), "the table names no cell with a group");
        names.sort_by(|a, b| natural(a).cmp(&natural(b)));
        let index: FxHashMap<&str, Group> = (names.iter().enumerate())
            .map(|(g, n)| (n.as_str(), Group(g)))
            .collect();
        let mut out = Self {
            named: vec![0; names.len()],
            names: Vec::new(),
            cells: FxHashMap::default(),
            loose: FxHashMap::default(),
            columns,
        };
        for (cell, given) in &rows {
            let g = index[sanitised(given).as_str()];
            let (map, barcode) = match cell.rsplit_once('@') {
                Some((barcode, tag)) => {
                    let batch = crate::qc::layout::split_stem(tag).0;
                    (out.cells.entry(batch).or_default(), barcode)
                }
                None => (&mut out.loose, cell.as_str()),
            };
            if let Some(was) = map.insert(barcode.into(), g) {
                anyhow::ensure!(was == g, "the cell {cell} is named twice, in two groups");
            } else {
                out.named[g.0] += 1;
            }
        }
        out.names = names;
        Ok(out)
    }
}

/// A group's value as a file stem: letters, digits, `.`, `_` and `-` kept,
/// anything else `_`; never starting with `.`, never empty.
pub(in crate::cluster) fn sanitised(value: &str) -> String {
    let mut out: String = (value.trim().chars())
        .map(
            |c| match c.is_ascii_alphanumeric() || matches!(c, '.' | '_' | '-') {
                true => c,
                false => '_',
            },
        )
        .collect();
    if out.starts_with('.') || out.is_empty() {
        out.insert(0, '_');
    }
    out
}

/// A name's sort key: by number when it is one, else after the numbers by
/// text.
fn natural(name: &str) -> (i64, &str) {
    (name.parse::<i64>().unwrap_or(i64::MAX), name)
}

/// The column the table names `wanted` (case aside), else the first of
/// `usual` there; refused, naming the columns found, when none fits.
fn choose<'a>(
    found: &'a [String],
    wanted: Option<&str>,
    usual: &[&str],
    flag: &str,
) -> anyhow::Result<&'a str> {
    let at = |name: &str| found.iter().find(|f| f.eq_ignore_ascii_case(name));
    let hit = match wanted {
        Some(w) => at(w),
        None => usual.iter().find_map(|u| at(u)),
    };
    hit.map(String::as_str).with_context(|| {
        let looked = match wanted {
            Some(w) => format!("no column {w}"),
            None => format!("none of {}", usual.join(", ")),
        };
        format!(
            "the table has {looked}: name one with {flag} (its columns: {})",
            found.join(", ")
        )
    })
}

/// The group column: `wanted`, else `label`, else `cluster`, else the one
/// column besides `cell` in a table of two.
fn group_column<'a>(
    found: &'a [String],
    wanted: Option<&str>,
    cell: &str,
) -> anyhow::Result<&'a str> {
    if wanted.is_none() && found.len() == 2 {
        if let Some(other) = found.iter().find(|f| *f != cell) {
            if !GROUP_COLUMNS
                .iter()
                .any(|g| found.iter().any(|f| f.eq_ignore_ascii_case(g)))
            {
                return Ok(other);
            }
        }
    }
    choose(found, wanted, &GROUP_COLUMNS, "--within-column")
}

type Rows = (Vec<(String, String)>, Option<(String, String)>);

/// A Parquet table's rows: each cell's name and group, any type read as
/// text; a null or empty group left out.
fn parquet_rows(path: &Path, group: Option<&str>, cell: Option<&str>) -> anyhow::Result<Rows> {
    let t = read_record_batch(path).with_context(|| format!("reading {}", path.display()))?;
    let found: Vec<String> = (t.schema().fields().iter())
        .map(|f| f.name().clone())
        .collect();
    let cell = choose(&found, cell, &CELL_COLUMNS, "--within-cell-column")?.to_string();
    let group = group_column(&found, group, &cell)?.to_string();
    let text = |name: &str| -> anyhow::Result<StringArray> {
        let column = t.column_by_name(name).expect("a column found");
        let cast = arrow::compute::cast(column, &DataType::Utf8)
            .with_context(|| format!("the column {name} cannot be read as text"))?;
        Ok(cast
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("cast")
            .clone())
    };
    let (cells, groups) = (text(&cell)?, text(&group)?);
    let rows = (0..t.num_rows())
        .filter(|&i| cells.is_valid(i) && groups.is_valid(i) && !groups.value(i).trim().is_empty())
        .map(|i| (cells.value(i).to_string(), groups.value(i).to_string()))
        .collect();
    Ok((rows, Some((cell, group))))
}

/// A TSV's (or CSV's: commas where the first line has no tab) rows. The
/// first line is a header when it names a cell column (`cell` or the
/// one `--within-cell-column` names); without one the cell is the first
/// column and the group the second, and naming a column is refused.
fn text_rows(path: &Path, group: Option<&str>, cell: Option<&str>) -> anyhow::Result<Rows> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let first = lines.next().context("the table is empty")?;
    let sep = if first.contains('\t') { '\t' } else { ',' };
    let split = |line: &str| -> Vec<String> {
        (line.split(sep))
            .map(|f| f.trim().trim_matches('"').to_string())
            .collect()
    };
    let head = split(first);
    let names_cell = |f: &String| match cell {
        Some(c) => f.eq_ignore_ascii_case(c),
        None => CELL_COLUMNS.iter().any(|u| f.eq_ignore_ascii_case(u)),
    };
    let header = head.iter().any(names_cell);
    let (at_cell, at_group, columns) = if header {
        let c = choose(&head, cell, &CELL_COLUMNS, "--within-cell-column")?;
        let g = group_column(&head, group, c)?;
        let at = |n: &str| head.iter().position(|h| h == n).expect("a column found");
        (at(c), at(g), Some((c.to_string(), g.to_string())))
    } else {
        anyhow::ensure!(
            group.is_none() && cell.is_none(),
            "{} has no header naming its columns (its first line: {first}): its cell is the first \
             column and its group the second, so --within-column and --within-cell-column are \
             not taken",
            path.display()
        );
        anyhow::ensure!(
            head.len() >= 2,
            "{}: write a cell and its group a line",
            path.display()
        );
        (0, 1, None)
    };
    let body = (!header).then_some(first).into_iter().chain(lines);
    let mut rows = Vec::new();
    for (n, line) in body.enumerate() {
        let f = split(line);
        let (Some(c), Some(g)) = (f.get(at_cell), f.get(at_group)) else {
            anyhow::bail!(
                "{} line {}: too few columns: {line}",
                path.display(),
                n + 1 + header as usize
            );
        };
        if !c.is_empty() && !g.is_empty() {
            rows.push((c.clone(), g.clone()));
        }
    }
    Ok((rows, columns))
}

#[cfg(test)]
mod tests;
