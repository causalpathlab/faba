//! `faba cluster --within`: the site blockmodel fitted within each cell
//! module of a count fit, the second round of two.
//!
//! Counts and site reads carry unequal weight (every gene counts in every
//! cell, a site has reads in few), so they are not fitted jointly. The
//! count fit (`faba cluster --modality count`) groups the cells first: each
//! in its main cell module. Each group with at least
//! `--within-min-cells` cells to fit is then fitted on its own site reads,
//! one after the other, into `within/{modality}/{group}/` under the output
//! directory: the tables and run record `faba cluster` writes for any site
//! fit. `{modality}_within.run.json` names them all (see [`GroupFit`]), and
//! `{modality}_within_cells.parquet` gives each cell its nested label
//! (`C3.2`: count module 3, site module 2 within it). `faba cluster-view`
//! opens the record as the count fit with a fit under each of its cell
//! modules.

use super::command::COUNT_MODALITY;
use super::{args, progress, read, write, ClusterArgs};
use crate::common::{column, read_record_batch};
use crate::run_record::RunRecord;
use anyhow::Context;
use arrow::array::{Array, Int32Array, StringArray};
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

/// The count fit's run record, as `faba cluster --modality count` names it.
pub(super) const COUNT_RECORD: &str = super::run::count::RECORD;

/// The within record of modality `m` (lowercased), which names each
/// group's fit.
pub(super) fn record_name(m: &str) -> String {
    format!("{m}_within.run.json")
}

/// A group of the count fit's cells: those of one cell module, numbered
/// from 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(super) struct Group(pub(super) usize);

impl Group {
    /// `C3` for the third cell module, as the count fit's view names it.
    pub(super) fn name(self) -> String {
        format!("C{}", self.0 + 1)
    }

    /// The group of a count fit's cell, from its main `role` (1..). None
    /// for a role that is no module's.
    fn of(role: i32) -> Option<Self> {
        Some(Group(usize::try_from(role.checked_sub(1)?).ok()?))
    }

    /// Its cell module as the count fit's table numbers it (1..).
    fn role(self) -> i32 {
        self.0 as i32 + 1
    }
}

/// What became of one group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub(in crate::cluster) enum Status {
    /// Fitted and written, in its directory.
    Fitted,
    /// Too few cells to fit.
    Skipped,
    /// Its fit failed (too few cells or sites past the filters, mostly).
    Failed,
    /// Cancelled while it was fitted: nothing written for it.
    Cancelled,
    /// The run was stopped before its turn.
    NotRun,
}

impl Status {
    /// In words.
    pub(in crate::cluster) fn said(self) -> &'static str {
        match self {
            Status::Fitted => "fitted",
            Status::Skipped => "skipped",
            Status::Failed => "failed",
            Status::Cancelled => "cancelled",
            Status::NotRun => "not run",
        }
    }
}

/// One group's fit, as the within record lists it under
/// `results.within.groups`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub(in crate::cluster) struct GroupFit {
    /// `C3`.
    pub group: String,
    /// Its cells among those read.
    pub cells: usize,
    pub status: Status,
    /// Why it was skipped or failed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub why: Option<String>,
    /// Its fit's directory, from the within record's (`within/m6a/C3`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dir: Option<String>,
    /// Its fit's run record, in that directory.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub record: Option<String>,
    /// Its chains' split R-hat.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rhat: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub converged: Option<bool>,
    /// Its chains were stopped early: its tables are from fewer draws.
    #[serde(default)]
    pub stopped_early: bool,
}

impl GroupFit {
    fn new(group: Group, cells: usize, status: Status) -> Self {
        Self {
            group: group.name(),
            cells,
            status,
            why: None,
            dir: None,
            record: None,
            rhat: None,
            converged: None,
            stopped_early: false,
        }
    }

    /// What became of it, in words: its status, its chains' R-hat and
    /// verdict, and why it was skipped or failed.
    pub(in crate::cluster) fn what(&self) -> String {
        let mut what = self.status.said().to_string();
        if let Some(rhat) = self.rhat {
            what.push_str(&format!(", R-hat {rhat:.2}"));
        }
        match self.converged {
            Some(true) => what.push_str(", converged"),
            Some(false) => what.push_str(", not converged"),
            None => {}
        }
        if self.stopped_early {
            what.push_str(", stopped early");
        }
        if let Some(why) = &self.why {
            what.push_str(&format!(": {why}"));
        }
        what
    }

    /// One line of the summary: the group, its cells, what became of it.
    pub(in crate::cluster) fn line(&self) -> String {
        format!("{}\t{} cells\t{}", self.group, self.cells, self.what())
    }
}

/// By barcode: a cell's main count module (1..) and its group.
type ByBarcode = FxHashMap<Box<str>, Group>;

/// By batch, then barcode.
type ByBatch = FxHashMap<Box<str>, ByBarcode>;

/// Whether batch `long` is batch `short` named with its folder
/// (`runs/a/b1` for `b1`, or for `a/b1`).
fn qualifies(long: &str, short: &str) -> bool {
    (long.strip_suffix(short)).is_some_and(|folder| folder.ends_with('/'))
}

/// A count fit's cells, each with its main cell module and its group.
#[derive(Debug)]
pub(super) struct CountCells {
    /// The count fit's run record, an absolute path.
    record: PathBuf,
    /// Where it is.
    dir: PathBuf,
    /// Its contents.
    json: serde_json::Value,
    /// Batch, then barcode: the cell's main module (1..) and its group.
    cells: ByBatch,
}

/// The site cells against the count fit's.
#[derive(Debug, Default, PartialEq)]
struct Grouped {
    /// Each site cell's group, in the order read; none for a cell the
    /// count fit has no group for.
    of: Vec<Option<Group>>,
    /// The cells of each group, as indices in that order.
    groups: BTreeMap<Group, Vec<usize>>,
}

impl Grouped {
    /// How many cells the count fit has a group for.
    fn matched(&self) -> usize {
        self.groups.values().map(Vec::len).sum()
    }
}

impl CountCells {
    /// The count fit `--within` names (its run record, or the directory it
    /// is in), its cells grouped. Refused when it is not a count fit's
    /// record.
    fn read(named: &str) -> anyhow::Result<Self> {
        let path = match Path::new(named) {
            dir if dir.is_dir() => dir.join(COUNT_RECORD),
            file => file.to_path_buf(),
        };
        const FIRST: &str = "run `faba cluster --modality count` first, and name its run record";
        anyhow::ensure!(
            path.is_file(),
            "--within {named}: no {COUNT_RECORD} there: {FIRST}"
        );
        let name = path.file_name().unwrap_or_default().to_string_lossy();
        let dir = crate::run_record::default_output_dir(&path.to_string_lossy());
        let json = super::browse::read_record(&dir, &name).map_err(anyhow::Error::msg)?;
        let modality = json.pointer("/results/tables/modality");
        anyhow::ensure!(
            json["job"].as_str() == Some("cluster")
                && modality.and_then(|m| m.as_str()) == Some(COUNT_MODALITY),
            "--within {named} is not a count fit's run record: {FIRST}"
        );
        // The cell table is the one the record names.
        let table = json
            .pointer("/results/tables/cells/0")
            .and_then(|t| t.as_str());
        let table = dir.join(table.with_context(|| format!("{named} names no cell table"))?);
        let cells = Self::cells(&table)
            .with_context(|| format!("reading the count fit's cells, {}", table.display()))?;
        Ok(Self {
            record: crate::run_record::absolute(&path.to_string_lossy()).into(),
            dir,
            json,
            cells,
        })
    }

    /// The cells of the count fit's cell table at `path`, those it kept.
    fn cells(path: &Path) -> anyhow::Result<ByBatch> {
        let t = read_record_batch(path)?;
        let barcode = column::<StringArray>(&t, "barcode")?;
        let batch = column::<StringArray>(&t, "batch")?;
        let role = column::<Int32Array>(&t, "role")?;
        if t.schema().column_with_name("level").is_some() {
            log::warn!(
                "{}: a count fit with a background (faba before 0.18.13); its background \
                 cells now go to their main modules: fit the counts again for modules \
                 without one",
                path.display()
            );
        }
        let mut out: ByBatch = FxHashMap::default();
        for i in 0..t.num_rows() {
            // A cell the count fit did not keep has no role, and no group.
            if !role.is_valid(i) {
                continue;
            }
            if let Some(group) = Group::of(role.value(i)) {
                out.entry(batch.value(i).into())
                    .or_default()
                    .insert(barcode.value(i).into(), group);
            }
        }
        anyhow::ensure!(!out.is_empty(), "the count fit kept no cells");
        Ok(out)
    }

    /// Where the count fit's matrices are: the output directories whose
    /// site matrices are fitted when none are named. One that is no longer
    /// there (the output was moved) is taken to be the record's own.
    fn folders(&self) -> Vec<Box<str>> {
        let matrices = self.json.pointer("/inputs/matrices");
        let mut out: Vec<Box<str>> = Vec::new();
        for path in (matrices.and_then(|m| m.as_array()).into_iter().flatten())
            .filter_map(|m| m["path"].as_str())
        {
            let recorded = crate::run_record::default_output_dir(path);
            let there = Path::new(path).exists();
            let folder = if there { &recorded } else { &self.dir };
            let folder: Box<str> = folder.to_string_lossy().into();
            if !out.contains(&folder) {
                out.push(folder);
            }
        }
        out
    }

    /// The count fit's batch of each of the sites' batches (those of
    /// `cells`). A batch met in several folders is named with its folder,
    /// as typed in each round (see [`args::qualify_batches`]): the same
    /// name, or else the one batch that is the other named with its folder,
    /// when no other batch of the sites is that one too. A count fit of one
    /// batch and sites of one are the same batch, whatever its names.
    fn batches<'a>(&'a self, cells: &'a [read::Cell]) -> FxHashMap<&'a str, &'a ByBarcode> {
        let mut ours: Vec<&str> = cells.iter().map(|c| &*c.batch).collect();
        ours.sort_unstable();
        ours.dedup();
        if let ([ours], 1) = (ours.as_slice(), self.cells.len()) {
            return self.cells.values().map(|theirs| (*ours, theirs)).collect();
        }
        // Each of the sites' batches with the count fit's batch it names.
        let alike: Vec<(&str, &str)> = (ours.iter())
            .filter_map(|&ours| {
                if let Some((theirs, _)) = self.cells.get_key_value(ours) {
                    return Some((ours, &**theirs));
                }
                let mut alike = (self.cells.keys())
                    .filter(|theirs| qualifies(ours, theirs) || qualifies(theirs, ours));
                let theirs = alike.next()?;
                alike.next().is_none().then_some((ours, &**theirs))
            })
            .collect();
        let mut named: FxHashMap<&str, usize> = FxHashMap::default();
        for (_, theirs) in &alike {
            *named.entry(theirs).or_default() += 1;
        }
        (alike.into_iter())
            .filter(|(_, theirs)| named[theirs] == 1)
            .map(|(ours, theirs)| (ours, &self.cells[theirs]))
            .collect()
    }

    /// `cells` (the sites') by their group in the count fit. Refused when
    /// none is among the count fit's cells: the two are of other cells, or
    /// name them otherwise.
    fn group(&self, cells: &[read::Cell]) -> anyhow::Result<Grouped> {
        let batches = self.batches(cells);
        let mut out = Grouped::default();
        for (i, cell) in cells.iter().enumerate() {
            let of = (batches.get(&*cell.batch)).and_then(|b| b.get(&cell.barcode).copied());
            if let Some(group) = of {
                out.groups.entry(group).or_default().push(i);
            }
            out.of.push(of);
        }
        if out.groups.is_empty() {
            let theirs = (self.cells.iter())
                .find_map(|(batch, cells)| Some(format!("{batch} {}", cells.keys().next()?)));
            let ours = cells.first().map(|c| format!("{} {}", c.batch, c.barcode));
            anyhow::bail!(
                "none of the {} cells with sites is among the count fit's (batch and barcode: \
                 a site cell is {}, a count cell {})",
                cells.len(),
                ours.unwrap_or_default(),
                theirs.unwrap_or_default()
            );
        }
        Ok(out)
    }
}

/// `faba cluster --within`: each group of the count fit's cells fitted on
/// its own site reads, one after the other (see the module). The modality
/// fitted; none when it was cancelled before any group was written.
pub(super) fn run(args: &ClusterArgs) -> anyhow::Result<Option<String>> {
    let named = args.within.as_deref().expect("--within");
    let mut args = args.clone();
    args::check(&args)?;
    let count = CountCells::read(named)?;
    if args.data_files.is_empty() {
        args.data_files = count.folders();
        anyhow::ensure!(
            !args.data_files.is_empty(),
            "the count fit's record names no inputs: name the site matrices, or their directory"
        );
    }
    let args = &args;
    let dir = match args.output.as_deref() {
        Some(o) => PathBuf::from(o),
        None => count.dir.clone(),
    };
    // The within record, and each group's, name the same inputs.
    let start = |dir: &Path| {
        RunRecord::start("cluster", &dir.to_string_lossy())
            .inputs("matrices", &args.data_files)
            .input("count_record", Some(&count.record.to_string_lossy()))
            .options(args)
    };
    let record = start(&dir);
    // A failure before any group is fitted leaves its record; a cancel
    // leaves nothing.
    let failed = |record: RunRecord, error: anyhow::Error| {
        if !error.is::<progress::Cancelled>() {
            record.finish(&Err::<(), _>(anyhow::anyhow!("{error:#}")));
        }
        super::unless_cancelled(args, error)
    };
    let title = format!("reading the sites to fit within {}", count.record.display());
    let (_, all) = super::watched(args, &title, false, |p| {
        read::read_inputs(&args::inputs(args)?, p)
    });
    let all = match all {
        Ok(all) => all,
        Err(e) => return failed(record, e),
    };
    let m = all.modality.to_ascii_lowercase();
    let name = record_name(&m);
    let mut record = record.file_name(&name);
    let grouped = match count.group(&all.cells) {
        Ok(grouped) => grouped,
        Err(e) => return failed(record, e),
    };
    let matched = grouped.matched();
    log::info!(
        "{m}: {} cells read, {} in the count fit's {} groups",
        all.cells.len(),
        matched,
        grouped.groups.len()
    );

    let mut fits: Vec<GroupFit> = Vec::new();
    let mut site_roles: Vec<Option<i32>> = vec![None; all.cells.len()];
    // Cancelled or stopped early: the groups left are not run.
    let mut halted = false;
    for (&group, cells) in &grouped.groups {
        let fitted = cells.iter().filter(|&&c| all.cells[c].fit).count();
        if halted {
            fits.push(GroupFit::new(group, cells.len(), Status::NotRun));
            continue;
        }
        if fitted < args.within_min_cells {
            let mut skipped = GroupFit::new(group, cells.len(), Status::Skipped);
            skipped.why = Some(format!(
                "{fitted} cells to fit, fewer than --within-min-cells {}",
                args.within_min_cells
            ));
            fits.push(skipped);
            continue;
        }
        let rel = format!("within/{m}/{}", group.name());
        let sub = dir.join(&rel);
        let data = read::SiteData {
            modality: all.modality.clone(),
            sites: all.sites.clone(),
            cells: cells.iter().map(|&c| all.cells[c].clone()).collect(),
            reads: all.reads.subset(cells),
        };
        let about = serde_json::json!({
            "group": group.name(),
            "cells": cells.len(),
            "count_record": count.record,
            "within_record": crate::run_record::absolute(&dir.join(&name).to_string_lossy()),
        });
        let title = format!(
            "fitting the site blockmodel within {} ({} of {} groups, {} cells) in {}",
            group.name(),
            fits.len() + 1,
            grouped.groups.len(),
            cells.len(),
            sub.display()
        );
        log::info!("within {}: {} cells", group.name(), cells.len());
        let (progress, out) = super::watched(args, &title, false, |p| {
            let mut record = start(&sub);
            record.set_result("within_group", &about);
            super::recorded(
                record,
                |record| super::site_fit(args, &sub, record, p, data),
                |_| (),
            )
        });
        let mut fit = GroupFit::new(group, cells.len(), Status::Fitted);
        match out {
            Ok(written) => {
                progress::spark::tell(&progress);
                let table = sub.join(format!("{m}_cell_roles.parquet"));
                match main_roles(&table, cells.len()) {
                    Ok(roles) => (cells.iter().zip(roles)).for_each(|(&c, r)| site_roles[c] = r),
                    Err(e) => log::warn!("{}: {e:#}", table.display()),
                }
                fit.dir = Some(rel);
                fit.record = Some(written.record);
                fit.converged = progress.verdict().map(|v| v.converged());
                fit.stopped_early = progress.stopped_early().is_some();
            }
            Err(e) if e.is::<progress::Cancelled>() => fit.status = Status::Cancelled,
            Err(e) => {
                log::warn!("within {}: {e:#}", group.name());
                fit.status = Status::Failed;
                fit.why = Some(format!("{e:#}"));
            }
        }
        // A cancel that came too late for this group stops the rest.
        halted = progress.cancelled()
            || progress.stopped_early().is_some()
            || progress.interrupted()
            || progress.asked_late();
        fits.push(fit);
    }

    let fitted = fits.iter().filter(|f| f.status == Status::Fitted).count();
    // Cancelled with nothing written; a group that failed is still said.
    if fitted == 0 && halted && fits.iter().all(|f| f.status != Status::Failed) {
        return super::unless_cancelled(args, progress::Cancelled.into());
    }
    let table = format!("{m}_within_cells.parquet");
    let rows: Vec<write::NestedRow> = (all.cells.iter().zip(&grouped.of).zip(&site_roles))
        .map(|((cell, of), &site_role)| write::NestedRow {
            cell,
            count_role: of.map(|g| g.role()),
            group: of.map(|g| g.name()),
            site_role,
        })
        .collect();
    record.set_result(
        "within",
        &serde_json::json!({
            "modality": m,
            "count_record": count.record,
            "min_cells": args.within_min_cells,
            "cells": {
                "read": all.cells.len(),
                "matched": matched,
                "unmatched": all.cells.len() - matched,
            },
            "model": args::recorded(args),
            "table": table,
            "groups": fits,
            "halted": halted,
        }),
    );
    let outcome = std::fs::create_dir_all(&dir)
        .map_err(anyhow::Error::from)
        .and_then(|()| write::nested_cells(&rows, &dir.join(&table)))
        .and_then(|()| {
            anyhow::ensure!(
                fitted > 0,
                "no group could be fitted: {}",
                fits.iter()
                    .map(|f| f.line().replace('\t', " "))
                    .collect::<Vec<_>>()
                    .join("; ")
            );
            Ok(())
        });
    super::save_command(args, &dir, &name);
    record.finish(&outcome);
    outcome?;
    let mut lines = vec![format!(
        "{m} within the count fit's cell modules, {fitted} of {} groups fitted:",
        fits.len()
    )];
    lines.extend(fits.iter().map(|f| format!(" {}", f.line())));
    lines.push(format!(
        "open them with: faba cluster-view {}",
        dir.join(&name).display()
    ));
    progress::spark::say(&lines, log::Level::Info);
    Ok(Some(m))
}

/// Each cell's main role (1.., none when not kept) in the cell table at
/// `path`, which is to have `cells` rows: a group's cells in the order
/// they were fitted.
fn main_roles(path: &Path, cells: usize) -> anyhow::Result<Vec<Option<i32>>> {
    let t = read_record_batch(path)?;
    let role = column::<Int32Array>(&t, "role")?;
    anyhow::ensure!(role.len() == cells, "{} rows for {cells} cells", role.len());
    Ok(role.iter().collect())
}

/// What a within record (its contents `json`) lists: its modality, the
/// count fit's run record, and each group's fit; none when `json` is not
/// one.
pub(in crate::cluster) fn listed(
    json: &serde_json::Value,
) -> Option<(String, PathBuf, Vec<GroupFit>)> {
    let within = json.pointer("/results/within")?;
    let groups = serde_json::from_value(within["groups"].clone()).ok()?;
    Some((
        within["modality"].as_str()?.to_string(),
        PathBuf::from(within["count_record"].as_str()?),
        groups,
    ))
}

#[cfg(test)]
mod tests;
