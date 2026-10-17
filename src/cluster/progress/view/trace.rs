//! The chains' log joint under their bars: a line per chain in its
//! label's colour, the burn-in's end marked; once done, whether a chain
//! lags the best.

use crate::site_analysis::miami::palette::{color, Palette};
use data_beans::interactive::ui::{DIM, HIGHLIGHT};
use ratatui::layout::Rect;
use ratatui::style::{Color, Style};
use ratatui::symbols::Marker;
use ratatui::text::{Line, Span};
use ratatui::widgets::{Axis, Chart, Dataset, GraphType};
use ratatui::Frame;

/// The fewest rows and columns a plot is drawn in: fewer and it goes.
const ROWS: u16 = 6;
const COLS: u16 = 30;

/// The burn-in's end down the plot: a plot's lines take a colour, not a
/// dim style.
pub(super) const MARK: Color = Color::DarkGray;

/// How far below the best chain a chain's mean log joint lags, in median
/// within-chain sds.
const FAR: f64 = 3.0;

/// Chain `i`'s colour, its label's and its trace's: the largest table's,
/// so that chains repeat a colour only past its 18.
pub(super) fn chain_style(i: usize) -> Style {
    let (r, g, b) = color(&Palette::Category20, i);
    Style::new().fg(Color::Rgb(r, g, b))
}

/// The median of `v` (taken in any order).
fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    (v[(n - 1) / 2] + v[n / 2]) / 2.0
}

/// The chains whose mean log joint over their traced points after the
/// `warmup` burn-in sweeps is more than [`FAR`] median within-chain sds
/// below the best chain's (the lower of two, or most stuck below one);
/// none to say with fewer than two chains, or a chain with fewer than two
/// points past its burn-in. A chain's level, not its mixing: that is the
/// R-hat's.
pub(super) fn lagging(traces: &[Vec<(usize, f64)>], warmup: usize) -> Option<Vec<usize>> {
    if traces.len() < 2 {
        return None;
    }
    let (means, sds): (Vec<f64>, Vec<f64>) = traces
        .iter()
        .map(|t| {
            let v: Vec<f64> = t.iter().filter(|p| p.0 > warmup).map(|p| p.1).collect();
            (v.len() >= 2).then(|| {
                let n = v.len() as f64;
                let mean = v.iter().sum::<f64>() / n;
                let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0);
                (mean, var.sqrt())
            })
        })
        .collect::<Option<Vec<_>>>()?
        .into_iter()
        .unzip();
    let best = means.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let bar = best - FAR * median(sds);
    Some((0..means.len()).filter(|&i| means[i] < bar).collect())
}

/// The line under the plot once every chain is done: whether a chain lags
/// after the burn-in, and which; none when there is nothing to say.
pub(super) fn verdict(traces: &[Vec<(usize, f64)>], warmup: usize) -> Option<Line<'static>> {
    let lag = lagging(traces, warmup)?;
    let names = lag.iter().map(usize::to_string).collect::<Vec<_>>();
    let (said, style) = match names.as_slice() {
        [] => ("no chain lags".to_string(), DIM),
        [k] => (format!("chain {k} lags"), HIGHLIGHT),
        _ => (format!("chains {} lag", names.join(", ")), HIGHLIGHT),
    };
    Some(Line::from(vec![
        Span::styled(" log joint after burn-in: ", DIM),
        Span::styled(said, style),
    ]))
}

/// The chains' `traces` in `area`, over sweeps `0..=total`: a line per
/// chain in its colour, a [`MARK`] line down at the end of the `warmup`
/// burn-in sweeps, the log joint's range that of each chain's points past
/// its first tenth (where it climbs from its start). Nothing when the area
/// is too small to read or nothing is traced yet.
pub(super) fn draw(
    frame: &mut Frame,
    area: Rect,
    traces: &[Vec<(usize, f64)>],
    (total, warmup): (usize, usize),
) {
    if area.height < ROWS || area.width < COLS {
        return;
    }
    let settled = traces.iter().flat_map(|t| &t[t.len() / 10..]).map(|p| p.1);
    let (lo, hi) = settled.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
        (lo.min(v), hi.max(v))
    });
    if lo > hi {
        return;
    }
    let pad = ((hi - lo) * 0.05).max(0.5);
    let (lo, hi) = (lo - pad, hi + pad);
    let points: Vec<Vec<(f64, f64)>> = traces
        .iter()
        .map(|t| t.iter().map(|&(n, v)| (n as f64, v)).collect())
        .collect();
    let mark = [(warmup as f64, lo), (warmup as f64, hi)];
    let line = |data, style| {
        Dataset::default()
            .marker(Marker::Braille)
            .graph_type(GraphType::Line)
            .style(style)
            .data(data)
    };
    let mut sets = Vec::new();
    if warmup > 0 && warmup < total {
        sets.push(line(&mark, Style::new().fg(MARK)));
    }
    sets.extend(
        points
            .iter()
            .enumerate()
            .map(|(i, p)| line(p, chain_style(i))),
    );
    let digits = usize::from(hi - lo < 10.0);
    let chart = Chart::new(sets)
        .x_axis(
            Axis::default()
                .title(Span::styled("sweep", DIM))
                .style(DIM)
                .bounds([0.0, total as f64])
                .labels(["0".to_string(), total.to_string()]),
        )
        .y_axis(
            Axis::default()
                .title(Span::styled("log joint", DIM))
                .style(DIM)
                .bounds([lo, hi])
                .labels([format!("{lo:.digits$}"), format!("{hi:.digits$}")]),
        );
    frame.render_widget(chart, area);
}
