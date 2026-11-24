//! The model the tables come from, as `faba cluster`'s blockmodel recorded
//! it beside them (`{m}_blockmodel.run.json`, see [`crate::run_record`]).

use crate::cluster::record_name;
use crate::cluster::settings::Recorded;
use std::path::Path;

/// The run record of the run that wrote `modality`'s tables in `dir`, its
/// file name and contents, or in one line why there is none.
pub(in crate::cluster::browse) fn record(
    dir: &Path,
    modality: &str,
) -> Result<(String, serde_json::Value), String> {
    let name = record_name(modality);
    if !dir.join(&name).is_file() {
        return Err(format!("no run record ({name}) beside the tables"));
    }
    read(dir, &name).map(|record| (name, record))
}

/// Refused when `record` (the run record `name`) is of one of the MCMC fits
/// of faba before 0.19, which faba no longer reads: what only they
/// recorded (a split R-hat, a crp fit's posterior, the chains' pooling, the
/// sampler's sweeps) and no fit by variational inference.
pub(in crate::cluster) fn refuse_mcmc(
    name: &str,
    record: &serde_json::Value,
) -> anyhow::Result<()> {
    let mcmc = [
        "/results/blockmodel",
        "/results/count/rhat_log_likelihood",
        "/results/posterior",
        "/results/pooling",
        "/options/sweeps",
    ]
    .iter()
    .any(|at| record.pointer(at).is_some_and(|v| !v.is_null()));
    anyhow::ensure!(
        record["job"].as_str() != Some("cluster")
            || !mcmc
            || super::trace::vi_fit(record).is_some(),
        "{name} is the record of an MCMC fit of faba before 0.19, which faba no longer \
         reads: fit again with `faba cluster`"
    );
    Ok(())
}

/// The run record `name` in `dir`, parsed, or in one line why it cannot be.
pub(in crate::cluster) fn read(dir: &Path, name: &str) -> Result<serde_json::Value, String> {
    let unread = |e: &dyn std::fmt::Display| format!("{name} cannot be read: {e}");
    let text = std::fs::read_to_string(dir.join(name)).map_err(|e| unread(&e))?;
    serde_json::from_str(&text).map_err(|e| unread(&e))
}

/// Why there are no site blockmodel settings to draw.
#[derive(Debug, Clone, PartialEq)]
pub(in crate::cluster) enum Unsettled {
    /// A count fit's record: its model is not the site blockmodel's.
    CountFit,
    /// No record, or none that says them, in one line.
    Unrecorded(String),
}

impl std::fmt::Display for Unsettled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unsettled::CountFit => f.write_str(
                "a count fit: x_cg ~ Poisson(s_c r_g + sum_ab theta_ca lambda_ab psi_gb); \
                 no plate drawn yet",
            ),
            Unsettled::Unrecorded(why) => f.write_str(why),
        }
    }
}

/// The blockmodel's settings in `record`, the run record `name`, or why
/// there are none (a record without them, or a count fit's).
pub(in crate::cluster::browse) fn settings(
    (name, record): &(String, serde_json::Value),
) -> Result<Recorded, Unsettled> {
    let unread = |why| Err(Unsettled::Unrecorded(why));
    let Some(model) = record.pointer("/results/model") else {
        return unread(format!("{name} holds no blockmodel settings"));
    };
    // A count fit's record names its tables' modality, which is its kind.
    if record
        .pointer("/results/tables/modality")
        .and_then(|m| m.as_str())
        == Some(crate::cluster::command::COUNT_MODALITY)
    {
        return Err(Unsettled::CountFit);
    }
    serde_json::from_value(model.clone()).or_else(|e| unread(format!("{name} cannot be read: {e}")))
}
