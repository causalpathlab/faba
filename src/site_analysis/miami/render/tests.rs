use super::*;

#[test]
fn area_path_empty_when_flat() {
    assert!(area_path(&[0.0, 0.0], 0.0, 1.0, 10.0, 5.0, 1.0, "#abc").is_empty());
    let p = area_path(&[0.0, 1.0], 0.0, 1.0, 10.0, 5.0, 1.0, "#abc");
    assert!(p.contains("<path"));
}

#[test]
fn rgb_hex_formats() {
    assert_eq!(rgb_hex((255, 0, 16)), "#ff0010");
}

/// The panels' tracks take their own height, or share the height asked for
/// after the fixed parts; too little for a readable track is refused.
#[test]
fn the_tracks_share_the_height_asked_for() {
    let (dpi, fixed) = (100.0, 300.0);
    assert_eq!(track_height(None, 2.0, fixed, dpi).unwrap(), 85.0);
    // 7 inches is 700 dots: 400 for 2 panels' 4 tracks.
    assert_eq!(track_height(Some(7.0), 2.0, fixed, dpi).unwrap(), 100.0);
    let e = track_height(Some(4.0), 2.0, fixed, dpi)
        .unwrap_err()
        .to_string();
    assert!(e.contains("at least 4.3 inches"), "{e}");
}
