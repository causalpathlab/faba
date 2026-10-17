//! The chains' log joint under their bars (and, saved, in `faba modules`'
//! trace pane): a line per chain in its label's colour, the burn-in's end
//! marked; once done, whether a chain lags the best.

use crate::cluster::progress::advice::{lagging, lags};
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

/// Chain `i`'s colour, its label's and its trace's: the largest table's,
/// so that chains repeat a colour only past its 18.
pub(in crate::cluster) fn chain_style(i: usize) -> Style {
    let (r, g, b) = color(&Palette::Category20, i);
    Style::new().fg(Color::Rgb(r, g, b))
}

/// The line under the plot once every chain is done: whether a chain lags
/// after the burn-in, and which; none when there is nothing to say.
pub(super) fn verdict(traces: &[Vec<(usize, f64)>], warmup: usize) -> Option<Line<'static>> {
    let lag = lagging(traces, warmup)?;
    let (said, style) = match lag.is_empty() {
        true => ("no chain lags".to_string(), DIM),
        false => (lags(&lag), HIGHLIGHT),
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
pub(in crate::cluster) fn draw(
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
