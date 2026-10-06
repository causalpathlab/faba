use super::*;

fn dir_with(files: &[&str]) -> tempfile::TempDir {
    let tmp = tempfile::tempdir().unwrap();
    for f in files {
        std::fs::write(tmp.path().join(f), b"").unwrap();
    }
    std::fs::create_dir(tmp.path().join("sub")).unwrap();
    tmp
}

fn at(i: &mut Inputs, name: &str) {
    assert!(i.bams.list.select(name));
}

#[test]
fn the_browser_lists_folders_and_bams_only() {
    let tmp = dir_with(&[
        "sample_A.bam",
        "sample_A.bam.bai",
        "notes.txt",
        "sample_B.bam",
    ]);
    let i = Inputs::new(tmp.path().to_path_buf());
    let names: Vec<&str> = i.bams.list.shown().map(|e| e.name.as_str()).collect();
    assert_eq!(names, ["..", "sub", "sample_A.bam", "sample_B.bam"]);
    assert!(has_index(&tmp.path().join("sample_A.bam")));
    assert!(!has_index(&tmp.path().join("sample_B.bam")));
}

#[test]
fn space_cycles_fg_bg_off_and_argv_follows() {
    let tmp = dir_with(&["sample_A.bam", "sample_B.bam", "control_A.bam"]);
    let mut i = Inputs::new(tmp.path().to_path_buf());
    for n in ["sample_A.bam", "sample_B.bam", "control_A.bam"] {
        at(&mut i, n);
        i.toggle();
    }
    i.toggle(); // control_A.bam is highlighted: fg to bg
                // Round the cycle (bg, off, fg) keeps one entry, back as fg.
    at(&mut i, "sample_B.bam");
    for _ in 0..3 {
        i.toggle();
    }
    assert_eq!(i.picked.len(), 3);
    i.gff = Some(tmp.path().join("g.gff"));
    i.genome = Some(tmp.path().join("x.fa"));
    let argv = i.argv();
    let p = |n: &str| tmp.path().join(n).to_string_lossy().into_owned();
    assert_eq!(
        argv,
        [
            p("sample_A.bam"),
            p("sample_B.bam"),
            "--control-bam".into(),
            p("control_A.bam"),
            "-g".into(),
            p("g.gff"),
            "-f".into(),
            p("x.fa")
        ]
    );
    assert_eq!(i.picked.len(), 3);
}

#[test]
fn problems_name_what_blocks_the_run() {
    let tmp = dir_with(&["sample_A.bam"]);
    let mut i = Inputs::new(tmp.path().to_path_buf());
    let p = i.problems();
    assert!(p.iter().any(|s| s.contains("no fg BAM")));
    assert!(p.iter().any(|s| s.contains("GFF")) && p.iter().any(|s| s.contains("genome")));
    assert!(p.iter().any(|s| s.contains("no output folder")));
    at(&mut i, "sample_A.bam");
    i.toggle();
    assert_eq!(
        i.suggested_output(),
        tmp.path().join("faba_out").to_string_lossy()
    );
    i.output = tmp.path().to_string_lossy().into_owned(); // not empty
    assert!(i
        .output_problem()
        .unwrap()
        .contains("already contains files"));
}

#[test]
fn batch_names_are_the_pipelines_own() {
    let tmp = dir_with(&["sample_A.bam", "control_A.bam"]);
    let mut i = Inputs::new(tmp.path().to_path_buf());
    for n in ["sample_A.bam", "control_A.bam"] {
        at(&mut i, n);
        i.toggle();
    }
    i.toggle();
    let paths: Vec<Box<str>> = i
        .picked
        .iter()
        .map(|x| x.path.to_string_lossy().into_owned().into_boxed_str())
        .collect();
    let want = crate::common::uniq_batch_names(&paths).unwrap();
    let got = i.batch_names();
    assert_eq!(got.len(), 2);
    for ((path, name), (p, w)) in got.iter().zip(i.picked.iter().zip(want.iter())) {
        assert_eq!(path, &p.path);
        assert_eq!(name.as_str(), &**w);
    }
}

#[test]
fn the_browser_lists_what_known_snps_accepts() {
    let tmp = dir_with(&["a.vcf", "b.vcf.gz", "c.bcf", "d.parquet", "e.txt"]);
    let b = FileRow::KnownSnps.browser(tmp.path().to_path_buf());
    let names: Vec<&str> = b.list.shown().map(|e| e.name.as_str()).collect();
    assert_eq!(
        names,
        ["..", "sub", "a.vcf", "b.vcf.gz", "c.bcf", "d.parquet"]
    );
}

#[cfg(unix)]
#[test]
fn a_linked_bam_keeps_the_name_it_was_given() {
    let tmp = dir_with(&["sample_A.bam"]);
    let link = tmp.path().join("sub").join("linked_A.bam");
    std::os::unix::fs::symlink(tmp.path().join("sample_A.bam"), &link).unwrap();
    let spelled = tmp.path().join("sub/./../sub/linked_A.bam");
    assert_eq!(normalize(&spelled), link);
}
