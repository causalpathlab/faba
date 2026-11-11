//! The sparklines' text: one line per restart on one scale, the climb left
//! out, the verdict and advice after.

use super::lines;
use crate::cluster::progress::advice::Verdict;

/// A restart's trace: its first 10 points a climb from far below, then
/// `level` at each of the next `n` points.
fn chain(n: usize, level: impl Fn(usize) -> f64) -> Vec<(usize, f64)> {
    (0..10)
        .map(|i| (i + 1, -1000.0))
        .chain((0..n).map(|i| (i + 11, level(i))))
        .collect()
}

/// The sparkline of `line`, its cells after the restart's name.
fn cells(line: &str) -> String {
    line.chars().filter(|c| ('▁'..='█').contains(c)).collect()
}

#[test]
fn each_restart_is_a_line_on_one_scale_after_its_climb() {
    let traces = [
        chain(90, |_| 10.0),
        chain(90, |_| 0.0),
        chain(90, |i| i as f64 / 89.0 * 10.0),
    ];
    let v = Verdict::vi(true, false, (90, 200), 3, (0.0, "read"));
    let out = lines(&traces, &v, 30);
    assert_eq!(out.len(), 1 + 3 + 2, "{out:#?}");
    assert!(
        out[0].contains("0.0 to 10.0"),
        "the climb left out: {}",
        out[0]
    );
    assert!(out[1].contains("restart 0") && out[3].contains("restart 2"));
    assert_eq!(
        cells(&out[1]),
        "█".repeat(30),
        "the highest restart at the top"
    );
    assert_eq!(cells(&out[2]), "▁".repeat(30), "the lowest at the bottom");
    let rising = cells(&out[3]);
    assert_eq!(rising.chars().count(), 30);
    assert!(rising.starts_with('▁') && rising.ends_with('█'), "{rising}");
    assert_eq!(out[4].trim(), v.said());
    assert_eq!(out[5].trim(), v.advice());
}

/// Fewer points than cells: a cell per point.
#[test]
fn a_short_trace_takes_a_cell_per_point() {
    let traces = [chain(5, |i| i as f64), chain(5, |i| 4.0 - i as f64)];
    let v = Verdict::vi(false, false, (5, 5), 2, (0.4, "read"));
    let out = lines(&traces, &v, 30);
    // A tenth of 15 points is one: the first climbing point is left out.
    assert_eq!(cells(&out[1]).chars().count(), 14);
    assert!(out[3].contains("ELBO still rising"), "{out:#?}");
    assert!(out[4].contains("--epochs 10"), "{out:#?}");
}

/// A flat trace (every chain at one level) draws mid-height, not a NaN.
#[test]
fn flat_traces_draw_at_mid_height() {
    let flat: Vec<(usize, f64)> = (1..=20).map(|n| (n, 3.0)).collect();
    let traces = [flat.clone(), flat];
    let v = Verdict::vi(true, false, (20, 40), 2, (0.0, "read"));
    let out = lines(&traces, &v, 10);
    assert_eq!(cells(&out[1]), "▄".repeat(10));
}
