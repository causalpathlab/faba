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

/// A cell type is one colour in every track: panels take their hue, not
/// their place.
#[test]
fn a_cell_type_keeps_its_colour_across_tracks() {
    let tmp = tempfile::tempdir().unwrap();
    let panel = |label: &str| PanelData {
        celltype: label.into(),
        hue: 0,
        epi_sites: vec![(10, 1.0)],
        depth_bins: vec![1.0; 4],
    };
    let opts = FigOpts {
        out_prefix: tmp.path().join("x").to_string_lossy().into(),
        width_in: 8.0,
        height_in: None,
        dpi: 100,
        palette: Palette::Auto,
        title: "GENE1".into(),
        top_label: "sites".into(),
    };
    let panels = [panel("a: CT1"), panel("b: CT1")];
    render_miami(&panels, &[], &BinEdges::new(0, 100, 4), &opts).unwrap();
    let svg = std::fs::read_to_string(tmp.path().join("x.miami.svg")).unwrap();
    let pal = palette::resolve(&Palette::Auto, 1);
    assert!(svg.contains(&rgb_hex(palette::color(&pal, 0))));
    assert!(!svg.contains(&rgb_hex(palette::color(&pal, 1))));
}

/// The left margin is 1.1 inches while the labels fit, else as wide as
/// the longest label needs.
#[test]
fn the_left_margin_fits_the_longest_label() {
    let (dpi, fs) = (100.0, 10.0);
    assert_eq!(
        left_margin(["all cells", "CT1"].into_iter(), dpi, fs),
        110.0
    );
    let long = "a long track label: CT1";
    let left = left_margin([long, "CT1"].into_iter(), dpi, fs);
    let end = label_x(dpi) + long.len() as f32 * fs * 0.55;
    assert!(
        end < left,
        "{long} ends at {end}, inside the plot from {left}"
    );
}
