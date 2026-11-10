//! How each cell module of a fit spreads over the batches, read off its cell
//! table after the fit: shared by the replicates of each condition, specific
//! to one condition (and even over its replicates), or held by a few batches.
//!
//! Batches are replicate x condition. A module's cells are counted in each
//! batch by their shares there and divided by the batch's cells, so a small
//! batch weighs as much as a large one. Within each condition of two or more
//! replicates, how evenly the module spreads over them is its normalized
//! entropy (1: evenly; 0: one replicate alone); the module's evenness is the
//! mean over conditions, each by the module's share in it. Uneven, the
//! module is batch-specific: a technical batch effect, as a condition's
//! replicates should agree. Even but mostly in one condition, it is
//! condition-specific: biology the design asks about, not a batch effect.

use std::collections::BTreeMap;

/// Under this evenness a module is batch-specific.
pub const EVEN: f64 = 0.8;

/// At least this share in one condition makes a module that condition's.
pub const CONDITION: f64 = 0.8;

/// The share of a module the batches named as holding a batch-specific
/// one hold together.
const MOST: f64 = 0.8;

/// A module with fewer cells' worth of shares than this is too small to
/// tell.
pub const FEWEST: f64 = 10.0;

/// The batches as replicates of conditions.
#[derive(Debug, Clone, PartialEq)]
pub struct Design {
    pub batches: Vec<String>,
    pub conditions: Vec<String>,
    /// Each batch's condition, by index into `conditions`.
    pub condition_of: Vec<usize>,
}

impl Design {
    /// The design of `batches`: their conditions from `groups` (batch to
    /// condition) where given; else each name less its replicate token (see
    /// [`condition_name`]). A batch is a condition of its own when nothing
    /// says otherwise; a module then has no replicates there to be judged
    /// on.
    pub fn of(batches: &[String], groups: &BTreeMap<String, String>) -> Self {
        let named: Vec<String> = (batches.iter())
            .map(|b| groups.get(b).cloned().unwrap_or_else(|| condition_name(b)))
            .collect();
        let mut conditions = named.clone();
        conditions.sort();
        conditions.dedup();
        let condition_of = (named.iter())
            .map(|c| conditions.binary_search(c).expect("listed"))
            .collect();
        Self {
            batches: batches.to_vec(),
            conditions,
            condition_of,
        }
    }

    /// Whether some condition has two replicates or more.
    pub fn replicated(&self) -> bool {
        (0..self.conditions.len()).any(|c| self.replicates(c).nth(1).is_some())
    }

    /// Each condition's batches.
    pub fn replicates(&self, c: usize) -> impl Iterator<Item = usize> + '_ {
        (0..self.batches.len()).filter(move |&b| self.condition_of[b] == c)
    }

    /// The batches by condition, then by name: the order they are shown in.
    pub fn order(&self) -> Vec<usize> {
        let mut order: Vec<usize> = (0..self.batches.len()).collect();
        order.sort_by(|&a, &b| {
            (self.condition_of[a], &self.batches[a]).cmp(&(self.condition_of[b], &self.batches[b]))
        });
        order
    }

    /// Each condition and its batches, a line each: `mut: rep1_mut, rep2_mut`.
    pub fn listing(&self) -> Vec<String> {
        (self.conditions.iter().enumerate())
            .map(|(c, name)| {
                let reps: Vec<&str> = self
                    .replicates(c)
                    .map(|b| self.batches[b].as_str())
                    .collect();
                format!("{name}: {}", reps.join(", "))
            })
            .collect()
    }
}

/// A batch's condition by its name: past a folder (`runA/mut`, a batch
/// named alike in two folders, each folder a replicate), less its
/// replicate token (`rep1_mut`, `mut-r2`, `wt_replicate3`: `rep`, `r` or
/// `replicate` then digits, split at `_`, `-` and `.`); the name itself
/// when it has none, or nothing else.
fn condition_name(batch: &str) -> String {
    let batch = batch.rsplit('/').next().unwrap_or(batch);
    let tokens: Vec<&str> = batch.split(['_', '-', '.']).collect();
    let is_replicate = |t: &str| {
        let t = t.to_ascii_lowercase();
        let digits = ["replicate", "rep", "r"]
            .iter()
            .find_map(|p| t.strip_prefix(p).filter(|d| !d.is_empty()));
        digits.is_some_and(|d| d.bytes().all(|b| b.is_ascii_digit()))
    };
    let rest: Vec<&str> = tokens
        .iter()
        .copied()
        .filter(|t| !is_replicate(t))
        .collect();
    match rest.len() {
        n if n == tokens.len() || n == 0 => batch.to_string(),
        _ => rest.join("_"),
    }
}

/// What a module is, by how it spreads.
#[derive(Debug, Clone, PartialEq)]
pub enum Label {
    /// Even over every condition's replicates, in more than one condition.
    Shared,
    /// Even over the replicates of the one condition it is mostly in.
    Condition(String),
    /// Uneven: the batches that hold most of it.
    Batch(Vec<String>),
    /// Mostly where no condition has two replicates to tell.
    Unknown,
    /// Too few cells to tell (see [`FEWEST`]).
    Few,
}

impl Label {
    pub fn word(&self) -> String {
        match self {
            Label::Shared => "shared".into(),
            Label::Condition(c) => format!("{c} only"),
            Label::Batch(b) => format!("batch: {}", b.join(", ")),
            Label::Unknown => "no replicates to tell".into(),
            Label::Few => "too few cells".into(),
        }
    }

    /// Whether the module is a batch effect.
    pub fn is_batch(&self) -> bool {
        matches!(self, Label::Batch(_))
    }
}

/// How one cell module spreads over the batches.
#[derive(Debug, Clone, PartialEq)]
pub struct Spread {
    /// Its share in each batch, each batch's cells weighing alike; summing
    /// to one.
    pub per_batch: Vec<f64>,
    /// Its evenness over each condition's replicates (see the module); none
    /// when no condition has two.
    pub evenness: Option<f64>,
    pub label: Label,
}

/// Each of the `modules` cell modules' spread, from the cells in each
/// batch (of `design.batches`) and the modules' shares summed there
/// (modules x batches).
pub fn spread(design: &Design, cells: &[f64], units: &[f64], modules: usize) -> Vec<Spread> {
    let nb = cells.len();
    (0..modules)
        .map(|k| {
            let rate: Vec<f64> = (0..nb)
                .map(|b| units[k * nb + b] / cells[b].max(1.0))
                .collect();
            let held: f64 = (0..nb).map(|b| units[k * nb + b]).sum();
            let mut s = module(design, &rate);
            if held < FEWEST {
                s.label = Label::Few;
            }
            s
        })
        .collect()
}

/// One module's spread from its rate in each batch (its cells there by
/// their shares, over the batch's cells).
fn module(design: &Design, rate: &[f64]) -> Spread {
    let total: f64 = rate.iter().sum::<f64>().max(1e-300);
    let per_batch: Vec<f64> = rate.iter().map(|r| r / total).collect();
    // Each condition's share: its replicates' mean rate, so a condition
    // of more replicates weighs no more.
    let mean: Vec<f64> = (0..design.conditions.len())
        .map(|c| {
            let (n, sum) = design
                .replicates(c)
                .fold((0, 0.0), |(n, s), b| (n + 1, s + rate[b]));
            sum / f64::from(n.max(1))
        })
        .collect();
    let all: f64 = mean.iter().sum::<f64>().max(1e-300);
    let per_condition: Vec<f64> = mean.iter().map(|m| m / all).collect();
    // Each condition of two replicates or more: its evenness, by its share.
    let (mut even, mut weight) = (0.0, 0.0);
    for (c, &share) in per_condition.iter().enumerate() {
        let reps: Vec<f64> = design.replicates(c).map(|b| rate[b]).collect();
        let sum: f64 = reps.iter().sum();
        if reps.len() < 2 || sum <= 0.0 {
            continue;
        }
        let entropy: f64 = (reps.iter())
            .filter(|&&r| r > 0.0)
            .map(|&r| -(r / sum) * (r / sum).ln())
            .sum();
        even += share * entropy / (reps.len() as f64).ln();
        weight += share;
    }
    // Judged by its replicates only when most of it is where there are any.
    let evenness = (weight >= 0.5).then(|| even / weight);
    let main = (0..per_condition.len())
        .max_by(|&a, &b| per_condition[a].total_cmp(&per_condition[b]))
        .unwrap_or(0);
    let one = design.conditions.len() > 1 && per_condition[main] >= CONDITION;
    let label = match evenness {
        Some(e) if e < EVEN => Label::Batch(holding(design, &per_batch)),
        Some(_) if one => Label::Condition(design.conditions[main].clone()),
        Some(_) => Label::Shared,
        None => Label::Unknown,
    };
    Spread {
        per_batch,
        evenness,
        label,
    }
}

/// The fewest batches that hold [`MOST`] of a module, largest first.
fn holding(design: &Design, per_batch: &[f64]) -> Vec<String> {
    let mut order: Vec<usize> = (0..per_batch.len()).collect();
    order.sort_by(|&a, &b| per_batch[b].total_cmp(&per_batch[a]));
    let mut held = 0.0;
    order
        .into_iter()
        .take_while(|&b| {
            let more = held < MOST;
            held += per_batch[b];
            more
        })
        .map(|b| design.batches[b].clone())
        .collect()
}

#[cfg(test)]
mod tests;
