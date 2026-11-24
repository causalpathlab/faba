//! The fit's convergence as `faba cluster` left it beside the tables: the
//! verdict and the schedule from the run record, and each restart's ELBO
//! from `{m}_trace.parquet`.

use crate::cluster::progress::advice::Verdict;
use crate::cluster::vi::Unit;
use crate::common::{column, read_record_batch};
use arrow::array::{Float64Array, Int32Array, Int64Array};
use serde_json::Value;
use std::path::Path;

/// Each restart's traced ELBO by step, its k-means start's `warmup` steps
/// first.
#[derive(Debug, Clone, PartialEq)]
pub struct Trace {
    pub restarts: Vec<Vec<(usize, f64)>>,
    pub warmup: usize,
}

/// What the fit's run record and trace say of its convergence.
#[derive(Debug)]
pub struct Fit {
    /// Steps each restart was scheduled, 0 when not recorded.
    pub steps: usize,
    /// The trace, or in one line why there is none.
    pub trace: Result<Trace, String>,
    /// The verdict on the fit (see [`Fit::verdict`]).
    verdict: Option<Verdict>,
}

/// What a `faba cluster` run record says of its fit by variational
/// inference: under `fit`, or (0.18) under the fit's kind, `count` or
/// `blocks`. None in a record of the MCMC fits of faba before 0.19.
pub(in crate::cluster) fn vi_fit(record: &Value) -> Option<&Value> {
    ["/results/fit", "/results/count/vi", "/results/blocks/vi"]
        .into_iter()
        .find_map(|at| record.pointer(at))
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
        let steps = number("/results/trace/sweeps").unwrap_or(0);
        let name = format!("{modality}_trace.parquet");
        let trace = match (dir.join(&name).is_file(), number("/results/trace/burn_in")) {
            (false, _) => Err(format!("no {name} beside the tables")),
            (true, None) => Err(match record {
                Ok((r, _)) => format!("{r} gives no k-means steps for {name}"),
                Err(why) => why.clone(),
            }),
            (true, Some(warmup)) => restarts(&dir.join(&name))
                .map_err(|e| format!("{name} cannot be read: {e:#}"))
                .and_then(|restarts| match restarts.iter().all(Vec::is_empty) {
                    true => Err(format!("{name} holds no points")),
                    false => Ok(Trace { restarts, warmup }),
                }),
        };
        let verdict = (record.as_ref().ok())
            .and_then(|(_, r)| vi_fit(r))
            .map(|vi| {
                let ran = (vi["restarts"].as_array().into_iter().flatten())
                    .filter(|r| r.get("elbo").is_some())
                    .count();
                let number = |k: &str| vi[k].as_u64().unwrap_or(0) as usize;
                let spread = match vi["unit"].as_str() {
                    Some(u) => (
                        vi["spread_per_unit"].as_f64().unwrap_or(0.0),
                        Unit::named(u),
                    ),
                    // 0.18: the unit in the field's name.
                    None => match vi["spread_per_count"].as_f64() {
                        Some(x) => (x, Unit::Count),
                        None => (vi["spread_per_read"].as_f64().unwrap_or(0.0), Unit::Read),
                    },
                };
                Verdict::vi(
                    vi["converged"].as_bool().unwrap_or(false),
                    vi["stopped_early"].as_bool().unwrap_or(false),
                    (number("epochs"), number("asked")),
                    ran,
                    spread,
                )
            });
        Self {
            steps,
            trace,
            verdict,
        }
    }

    /// The verdict on the fit, as its record gives it; none without a
    /// record.
    pub fn verdict(&self) -> Option<&Verdict> {
        self.verdict.as_ref()
    }

    /// The steps a plot of the trace spans: those scheduled, or the last
    /// traced when later (a record that gives none).
    pub fn span(&self) -> usize {
        let last = self.trace.iter().flat_map(|t| t.restarts.iter().flatten());
        self.steps.max(last.map(|p| p.0).max().unwrap_or(0))
    }
}

/// Each restart's points in the trace table at `path`, by restart and
/// then in the table's order (its columns named as the MCMC fits named
/// them: `chain`, `sweep`, `log_joint`).
fn restarts(path: &Path) -> anyhow::Result<Vec<Vec<(usize, f64)>>> {
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
