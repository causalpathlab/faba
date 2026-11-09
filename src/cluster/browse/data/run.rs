//! The model the tables come from, as `faba cluster`'s blockmodel recorded
//! it beside them (`{m}_blockmodel.run.json`, see [`crate::run_record`]).

use crate::cluster::unified::Recorded;
use crate::cluster::{record_name, ModelKind};
use std::path::Path;

/// The run record of the run that wrote `modality`'s tables in `dir`, its
/// file name and contents, or in one line why there is none. Older runs
/// recorded the blockmodel in `{m}_cluster.run.json`, which a later crp
/// run replaces: read when there is no blockmodel record, saying so.
pub(in crate::cluster::browse) fn record(
    dir: &Path,
    modality: &str,
) -> Result<(String, serde_json::Value), String> {
    let own = record_name(modality, ModelKind::Blocks);
    let older = record_name(modality, ModelKind::Crp);
    let name = if dir.join(&own).is_file() {
        own
    } else if dir.join(&older).is_file() {
        log::warn!(
            "no {own} beside the tables: their model is read from {older}, as older runs \
             left it (a later `--model crp` run would have replaced it)"
        );
        older
    } else {
        return Err(format!("no run record ({own}) beside the tables"));
    };
    read(dir, &name).map(|record| (name, record))
}

/// The run record `name` in `dir`, parsed, or in one line why it cannot be.
pub(in crate::cluster) fn read(dir: &Path, name: &str) -> Result<serde_json::Value, String> {
    let unread = |e: &dyn std::fmt::Display| format!("{name} cannot be read: {e}");
    let text = std::fs::read_to_string(dir.join(name)).map_err(|e| unread(&e))?;
    serde_json::from_str(&text).map_err(|e| unread(&e))
}

/// Why a count fit's record draws no plate: its model is not the
/// blockmodel's.
pub(in crate::cluster::browse) const COUNT_FIT: &str = "a count fit: x_cg ~ Poisson(s_c r_g + \
     sum_ab theta_ca lambda_ab psi_gb); no plate drawn yet";

/// The blockmodel's settings in `record`, the run record `name`, or in one
/// line why there are none (an older record, a crp run's, or a count
/// fit's, see [`COUNT_FIT`]).
pub(in crate::cluster::browse) fn settings(
    (name, record): &(String, serde_json::Value),
) -> Result<Recorded, String> {
    let Some(model) = record.pointer("/results/model") else {
        return Err(format!("{name} holds no blockmodel settings"));
    };
    // A count fit's record names its tables' modality, which is its kind.
    if record
        .pointer("/results/tables/modality")
        .and_then(|m| m.as_str())
        == Some(crate::cluster::command::COUNT_MODALITY)
    {
        return Err(COUNT_FIT.to_string());
    }
    serde_json::from_value(model.clone()).map_err(|e| format!("{name} cannot be read: {e}"))
}
