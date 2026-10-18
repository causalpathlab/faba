//! The fit's convergence as `faba cluster` left it beside the tables: the
//! split R-hat and the schedule from the run record, and each chain's log
//! joint from `{m}_trace.parquet`.

use crate::cluster::progress::advice::Verdict;
use crate::common::{column, read_record_batch};
use arrow::array::{Float64Array, Int32Array, Int64Array};
use serde_json::Value;
use std::path::Path;

/// Each chain's traced log joint, by sweep, and the burn-in's sweeps.
#[derive(Debug, Clone, PartialEq)]
pub struct Trace {
    pub chains: Vec<Vec<(usize, f64)>>,
    pub burn_in: usize,
}

/// What the fit's run record and trace say of its convergence.
#[derive(Debug)]
pub struct Fit {
    /// Sweeps each chain ran, 0 when not recorded.
    pub sweeps: usize,
    /// The trace, or in one line why there is none.
    pub trace: Result<Trace, String>,
    /// The verdict on the fit (see [`Fit::verdict`]).
    verdict: Option<Verdict>,
}

impl Fit {
    /// The fit of `modality`'s tables in `dir`, from `record` (the run
    /// record's name and contents, or why there is none).
    pub(super) fn read(
        dir: &Path,
        modality: &str,
        record: &Result<(String, Value), String>,
    ) -> Self {
        let at = |key: &str| record.as_ref().ok()?.1.pointer(key);
        let number = |key: &str| at(key).and_then(Value::as_u64).map(|n| n as usize);
        // Not finite when the record has none (too few draws); none
        // without a blockmodel record.
        let rhat = at("/results/blockmodel").map(|b| {
            b.get("rhat_log_joint")
                .and_then(Value::as_f64)
                .unwrap_or(f64::NAN)
        });
        let sweeps = number("/results/trace/sweeps")
            .or_else(|| number("/options/sweeps"))
            .unwrap_or(0);
        let ran = number("/results/blockmodel/chains")
            .or_else(|| number("/options/chains"))
            .unwrap_or(0);
        let name = format!("{modality}_trace.parquet");
        let trace = match (dir.join(&name).is_file(), number("/results/trace/burn_in")) {
            (false, _) => Err(format!(
                "no {name} beside the tables: older fits kept no trace (rerun with --refit)"
            )),
            (true, None) => Err(match record {
                Ok((r, _)) => format!("{r} gives no burn-in for {name}"),
                Err(why) => why.clone(),
            }),
            (true, Some(burn_in)) => chains(&dir.join(&name))
                .map_err(|e| format!("{name} cannot be read: {e:#}"))
                .and_then(|chains| match chains.iter().all(Vec::is_empty) {
                    true => Err(format!("{name} holds no points")),
                    false => Ok(Trace { chains, burn_in }),
                }),
        };
        // What pooling the chains found short of convergence: none in an
        // older record.
        let unsettled: Vec<String> = at("/results/pooling/unsettled")
            .and_then(Value::as_array)
            .map(|v| {
                v.iter()
                    .filter_map(Value::as_str)
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default();
        let verdict = rhat.map(|rhat| {
            match &trace {
                Ok(t) => Verdict::of(&t.chains, t.burn_in, sweeps, rhat),
                Err(_) => Verdict::of(&vec![Vec::new(); ran], 0, sweeps, rhat),
            }
            .with_unsettled(unsettled)
        });
        Self {
            sweeps,
            trace,
            verdict,
        }
    }

    /// The verdict on the fit (see [`Verdict::of`]): its trace's and its
    /// R-hat's, without a trace the R-hat's alone over the chains the
    /// record gives, and on what the record says pooling found unsettled;
    /// none without an R-hat.
    pub fn verdict(&self) -> Option<&Verdict> {
        self.verdict.as_ref()
    }

    /// The sweeps a plot of the trace spans: those run, or the last traced
    /// when later (a record that gives none).
    pub fn span(&self) -> usize {
        let last = self.trace.iter().flat_map(|t| t.chains.iter().flatten());
        self.sweeps.max(last.map(|p| p.0).max().unwrap_or(0))
    }
}

/// Each chain's points in the trace table at `path`, by chain and then in
/// the table's order.
fn chains(path: &Path) -> anyhow::Result<Vec<Vec<(usize, f64)>>> {
    let t = read_record_batch(path)?;
    let chain: &Int32Array = column(&t, "chain")?;
    let sweep: &Int64Array = column(&t, "sweep")?;
    let joint: &Float64Array = column(&t, "log_joint")?;
    let mut out: Vec<Vec<(usize, f64)>> = Vec::new();
    for i in 0..t.num_rows() {
        let c = usize::try_from(chain.value(i))?;
        if out.len() <= c {
            out.resize(c + 1, Vec::new());
        }
        out[c].push((usize::try_from(sweep.value(i))?, joint.value(i)));
    }
    Ok(out)
}
