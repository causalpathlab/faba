//! The report's trace, under the Model section: each chain's log joint as
//! it ran (from `{m}_trace.parquet`), in its colour in the terminal, the
//! burn-in to mark, and the verdict and its advice (see
//! [`crate::cluster::progress::advice`]).

use super::super::data::Fit;
use crate::report::thousandth;
use crate::site_analysis::miami::palette::{color, Palette};
use serde::Serialize;

/// What the page reads of the trace: the chains and the verdict, or why
/// there is none.
#[derive(Serialize)]
#[serde(untagged)]
pub(super) enum Section {
    Drawn(Box<Drawn>),
    Omitted { why: String },
}

/// The chains, the sweeps run and the burn-in's, the split R-hat (none
/// when not known) and the verdict's words and advice.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Drawn {
    chains: Vec<Chain>,
    sweeps: usize,
    burn_in: usize,
    rhat: Option<f64>,
    converged: bool,
    said: String,
    advice: String,
}

/// A chain: its name, as the terminal numbers it, its colour there, and
/// its `[sweep, log joint]` points (some), the log joint to a thousandth.
#[derive(Serialize)]
struct Chain {
    name: String,
    color: String,
    points: Vec<(usize, f64)>,
}

/// The section for the fit `fit` (see [`super::super::Modules::fit`]).
pub(super) fn section(fit: &Fit) -> Section {
    let (t, verdict) = match (&fit.trace, fit.verdict()) {
        (Ok(t), Some(v)) => (t, v),
        (Err(why), _) => return Section::Omitted { why: why.clone() },
        (Ok(_), None) => {
            return Section::Omitted {
                why: "the run record gives no R-hat for the trace".to_string(),
            }
        }
    };
    // A chain with no points (none in the table) has nothing to draw.
    let chains = (t.chains.iter().enumerate())
        .filter(|(_, points)| !points.is_empty())
        .map(|(i, points)| {
            let (r, g, b) = color(&Palette::Category20, i);
            Chain {
                name: format!("chain {i}"),
                color: format!("#{r:02x}{g:02x}{b:02x}"),
                points: points.iter().map(|&(n, v)| (n, thousandth(v))).collect(),
            }
        })
        .collect();
    Section::Drawn(Box::new(Drawn {
        chains,
        sweeps: fit.span(),
        burn_in: t.burn_in,
        rhat: verdict.rhat.is_finite().then_some(verdict.rhat),
        converged: verdict.converged(),
        said: verdict.said(),
        advice: verdict.advice(),
    }))
}
