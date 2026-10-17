//! The model the tables come from, as `faba cluster`'s blockmodel recorded
//! it beside them (`{m}_blockmodel.run.json`, see [`crate::run_record`]).

use crate::cluster::unified::Recorded;
use crate::cluster::{record_name, ModelKind};
use std::path::Path;

/// The settings of the run that wrote `modality`'s tables in `dir`, or in
/// one line why there are none (no record, or an older one). Older runs
/// recorded the blockmodel in `{m}_cluster.run.json`, which a later crp run
/// replaces: read when there is no blockmodel record, saying so.
pub(in crate::cluster::browse) fn model(dir: &Path, modality: &str) -> Result<Recorded, String> {
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
    let unread = |e: &dyn std::fmt::Display| format!("{name} cannot be read: {e}");
    let text = std::fs::read_to_string(dir.join(&name)).map_err(|e| unread(&e))?;
    let record: serde_json::Value = serde_json::from_str(&text).map_err(|e| unread(&e))?;
    let Some(model) = record.pointer("/results/model") else {
        return Err(format!("{name} holds no blockmodel settings"));
    };
    serde_json::from_value(model.clone()).map_err(|e| unread(&e))
}
