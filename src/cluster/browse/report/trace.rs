//! The report's trace, under the Model section: each restart's ELBO as it
//! ran (from `{m}_trace.parquet`), in its colour in the terminal, the end
//! of its k-means start to mark, and the verdict and its advice (see
//! [`crate::cluster::progress::advice`]).

use super::super::data::Fit;
use crate::report::thousandth;
use crate::site_analysis::miami::palette::{color, Palette};
use serde::Serialize;

/// What the page reads of the trace: the restarts and the verdict, or why
/// there is none.
#[derive(Serialize)]
#[serde(untagged)]
pub(super) enum Section {
    Drawn(Box<Drawn>),
    Omitted { why: String },
}

/// The restarts, the steps run and the k-means start's, and the verdict's
/// words and advice (none without a run record).
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Drawn {
    restarts: Vec<Restart>,
    steps: usize,
    warmup: usize,
    converged: bool,
    said: Option<String>,
    advice: Option<String>,
}

/// A restart: its name, as the terminal numbers it, its colour there, and
/// its `[step, ELBO]` points (some), the ELBO to a thousandth.
#[derive(Serialize)]
struct Restart {
    name: String,
    color: String,
    points: Vec<(usize, f64)>,
}

/// The section for the fit `fit` (see [`super::super::Modules::fit`]).
pub(super) fn section(fit: &Fit) -> Section {
    let t = match &fit.trace {
        Ok(t) => t,
        Err(why) => return Section::Omitted { why: why.clone() },
    };
    // A restart with no points (none in the table) has nothing to draw.
    let restarts = (t.restarts.iter().enumerate())
        .filter(|(_, points)| !points.is_empty())
        .map(|(i, points)| {
            let (r, g, b) = color(&Palette::Category20, i);
            Restart {
                name: format!("restart {i}"),
                color: format!("#{r:02x}{g:02x}{b:02x}"),
                points: points.iter().map(|&(n, v)| (n, thousandth(v))).collect(),
            }
        })
        .collect();
    let verdict = fit.verdict();
    Section::Drawn(Box::new(Drawn {
        restarts,
        steps: fit.span(),
        warmup: t.warmup,
        converged: verdict.is_none_or(|v| v.converged()),
        said: verdict.map(|v| v.said()),
        advice: verdict.map(|v| v.advice()),
    }))
}
