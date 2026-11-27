//! The raw view's bands: equal, or by size within their cap.

use super::*;

fn widths(sizes: &[usize], by: RawBands, cap: (f64, f64)) -> Vec<f64> {
    bands(sizes, 100.0, by, cap)
        .iter()
        .map(|b| (b.1 * 1e6).round() / 1e6)
        .collect()
}

#[test]
fn equal_bands_share_the_length_evenly() {
    assert_eq!(
        widths(&[1, 9], RawBands::Equal, (0.5, 2.0)),
        vec![50.0, 50.0]
    );
}

#[test]
fn size_bands_follow_the_sizes_within_their_cap() {
    // Uncapped: as the sizes are.
    assert_eq!(
        widths(&[1, 3], RawBands::Size, (0.1, 10.0)),
        vec![25.0, 75.0]
    );
    // A small module kept at half an equal band, the rest to the others.
    let w = widths(&[1, 50, 49], RawBands::Size, (0.5, 2.0));
    assert!((w[0] - 100.0 / 3.0 * 0.5).abs() < 1e-6, "{w:?}");
    assert!((w.iter().sum::<f64>() - 100.0).abs() < 1e-4, "{w:?}");
    // None more than twice an equal band.
    let w = widths(&[100, 1, 1, 1], RawBands::Size, (0.5, 2.0));
    assert!(w[0] <= 50.0 + 1e-6, "{w:?}");
    assert!((w.iter().sum::<f64>() - 100.0).abs() < 1e-4, "{w:?}");
}

#[test]
fn crowded_names_are_spread_apart_within_the_length() {
    // Apart already: as they are.
    assert_eq!(
        spread(&[5.0, 20.0, 40.0], 10.0, 50.0),
        vec![5.0, 20.0, 40.0]
    );
    // Two near each other: the second pushed on.
    assert_eq!(spread(&[5.0, 8.0, 40.0], 10.0, 50.0), vec![5.0, 15.0, 40.0]);
    // Pushed past the end: back from it.
    assert_eq!(
        spread(&[30.0, 45.0, 48.0], 10.0, 50.0),
        vec![30.0, 40.0, 50.0]
    );
}

#[test]
fn names_too_many_for_the_length_are_spread_evenly_within_it() {
    assert_eq!(spread(&[1.0, 2.0, 3.0], 10.0, 15.0), vec![2.5, 7.5, 12.5]);
}
