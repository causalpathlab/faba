//! `docs`: the method write-ups.

use crate::mcp::args::Arguments;
use crate::mcp::run::command;
use anyhow::Result;
use legume_mcp::{Call, Content};

pub(in crate::mcp) fn docs(a: &Arguments, cx: &Call) -> Result<Vec<Content>> {
    let mut argv = vec!["docs".to_string()];
    argv.extend(a.strings("query")?);
    let text = cx.run_reading(command(&argv)?)?;
    Ok(vec![Content::Text(text)])
}
