//! The end of a `faba cluster` fit on stderr: each chain's log joint as a
//! sparkline, all on one scale, then the verdict and its advice (see
//! [`super::advice`]).

use super::advice::Verdict;
use super::Progress;
use std::io::IsTerminal;

/// A sparkline's cells, low to high.
const BLOCKS: [char; 8] = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];

/// Cells a chain's sparkline takes at most.
const WIDTH: usize = 60;

/// The lines for `traces` and their `verdict`: what the sparklines span,
/// a sparkline per chain at most `width` cells wide, each cell the mean of
/// its points, every chain on one scale; then the verdict and its advice.
/// Each chain's first tenth of points is left out: there it climbs from
/// its start, which would flatten the rest.
pub fn lines(traces: &[Vec<(usize, f64)>], verdict: &Verdict, width: usize) -> Vec<String> {
    let settled: Vec<&[(usize, f64)]> = traces.iter().map(|t| &t[t.len() / 10..]).collect();
    let all = || settled.iter().flat_map(|t| t.iter());
    let (lo, hi) = all().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
        (lo.min(p.1), hi.max(p.1))
    });
    let from = all().map(|p| p.0).min().unwrap_or(0);
    let to = all().map(|p| p.0).max().unwrap_or(0);
    let mut out = vec![format!(
        " log joint per chain, sweeps {from} to {to}, {lo:.1} to {hi:.1}"
    )];
    out.extend(
        (settled.iter().enumerate())
            .map(|(i, t)| format!(" chain {i:<3}{}", sparkline(t, (lo, hi), width))),
    );
    out.push(format!(" {}", verdict.said()));
    out.push(format!(" {}", verdict.advice()));
    out
}

/// `points` in at most `width` cells, each the mean of its share of them,
/// its height on the scale `lo..=hi` (mid-height when that is one value).
fn sparkline(points: &[(usize, f64)], (lo, hi): (f64, f64), width: usize) -> String {
    let n = points.len();
    let cells = n.min(width);
    (0..cells)
        .map(|j| {
            let share = &points[j * n / cells..(j + 1) * n / cells];
            let mean = share.iter().map(|p| p.1).sum::<f64>() / share.len() as f64;
            let top = (BLOCKS.len() - 1) as f64;
            let at = if hi > lo {
                ((mean - lo) / (hi - lo) * top).round()
            } else {
                top / 2.0
            };
            BLOCKS[at as usize]
        })
        .collect()
}

/// Once a fit is summarised, its [`lines`] on stderr when it is a
/// terminal, else in the log (as warnings when it has not converged);
/// nothing when no chain was traced.
pub fn tell(progress: &Progress) {
    let traces = progress.traces();
    let Some(verdict) = progress.verdict() else {
        return;
    };
    if traces.iter().all(Vec::is_empty) {
        return;
    }
    let lines = lines(&traces, &verdict, WIDTH);
    if std::io::stderr().is_terminal() {
        lines.iter().for_each(|l| eprintln!("{l}"));
    } else {
        let level = match verdict.converged() {
            true => log::Level::Info,
            false => log::Level::Warn,
        };
        lines
            .iter()
            .for_each(|l| log::log!(level, "{}", l.trim_start()));
    }
}

#[cfg(test)]
mod tests;
