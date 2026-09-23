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
