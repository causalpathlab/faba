use super::*;

fn touch(dir: &Path, name: &str) {
    std::fs::write(dir.join(name), b"").unwrap();
}

#[test]
fn a_directory_gives_its_site_matrices_by_batch_and_its_site_table() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    for f in [
        "b2_m6a_site.zarr.zip",
        "b1_m6a_site.zarr.zip",
        "b1_m6a.zarr.zip",
        "b1_atoi_site.zarr.zip",
        "m6a_sites.parquet",
        "pipeline_summary.json",
    ] {
        touch(dir, f);
    }
    let m6a = output_sites(dir, None).unwrap();
    assert_eq!(&*m6a.modality, "m6a", "m6A first");
    let batches: Vec<&str> = m6a.matrices.iter().map(|(b, _)| &**b).collect();
    assert_eq!(batches, ["b1", "b2"], "site matrices only, by batch");
    assert!(m6a.site_table.unwrap().ends_with("m6a_sites.parquet"));
    let atoi = output_sites(dir, Some("atoi")).unwrap();
    assert_eq!(atoi.matrices.len(), 1);
    assert!(atoi.site_table.is_none());
    assert!(output_sites(dir, Some("apa")).is_err());

    // A run record names its folder; a matrix or another file does not.
    let record = dir.join("pipeline_summary.json");
    let named = crate::run_record::output_dir_of(&record.to_string_lossy()).unwrap();
    assert_eq!(named, dir);
    let matrix = dir.join("b1_m6a.zarr.zip");
    assert!(crate::run_record::output_dir_of(&matrix.to_string_lossy()).is_none());
}

#[test]
fn the_default_output_and_title_follow_the_site_tables() {
    let one: Vec<Box<str>> = vec!["out/m6a_sites.parquet".into()];
    assert_eq!(&*default_output(&one, ".pwm.tsv"), "m6a_sites.pwm.tsv");
    assert_eq!(title(&one), "m6a_sites.parquet");
    let two: Vec<Box<str>> = vec!["a/m6a_sites.parquet".into(), "b/atoi_sites.parquet".into()];
    assert_eq!(
        &*default_output(&two, ".metagene.tsv"),
        "m6a_sites_combined.metagene.tsv"
    );
    assert_eq!(title(&two), "m6a_sites.parquet + atoi_sites.parquet");
}

#[test]
fn site_matrices_come_from_the_first_modality_that_has_them() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    // m6A has only a site table; A-to-I has a matrix.
    for f in ["m6a_sites.parquet", "b1_atoi_site.zarr.zip"] {
        touch(dir, f);
    }
    assert_eq!(&*output_sites(dir, None).unwrap().modality, "m6a");
    let sites = output_site_matrices(dir, None).unwrap();
    assert_eq!(&*sites.modality, "atoi");
    assert_eq!(sites.matrices.len(), 1);
    // Asked for by name, a modality without matrices is an error.
    assert!(output_site_matrices(dir, Some("m6a")).is_err());
}
