//! `docs`: the method write-ups, with this server's version.

use crate::mcp::args::Arguments;
use crate::mcp::run::command;
use anyhow::Result;
use legume_mcp::{Call, Content};
use serde_json::json;

pub(in crate::mcp) fn docs(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let mut argv = vec!["docs".to_string()];
    argv.extend(a.strings("query")?);
    let text = cx.run_reading(command(&argv)?)?;
    let version = json!({ "version": env!("CARGO_PKG_VERSION") });
    Ok(vec![Content::Text(text), Content::Json(version)])
}
