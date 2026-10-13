use super::*;
use clap::FromArgMatches;

fn parse(args: &[&str]) -> ScanPwmArgs {
    let cmd = crate::faba_command()
        .find_subcommand("pwm")
        .unwrap()
        .clone();
    let m = cmd
        .try_get_matches_from(std::iter::once("pwm").chain(args.iter().copied()))
        .unwrap();
    ScanPwmArgs::from_arg_matches(&m).unwrap()
}

fn boxed(names: &[&str]) -> Vec<Box<str>> {
    names.iter().map(|&n| Box::from(n)).collect()
}

#[test]
fn bams_are_read_and_the_rest_are_sites() {
    let args = parse(&[
        "-s",
        "x.parquet",
        "a.bam,B.CRAM",
        "y.parquet",
        "--batch-process",
    ]);
    let (sites, bams) = args.inputs().unwrap().unwrap();
    assert_eq!(sites, boxed(&["x.parquet", "y.parquet"]), "-s first");
    assert_eq!(bams, boxed(&["a.bam", "B.CRAM"]));
    assert!(args.output.is_none(), "-o is optional");
}

#[test]
fn an_output_directory_or_its_run_record_gives_its_site_table() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    for f in ["m6a_sites.parquet", "pipeline_summary.json"] {
        std::fs::write(dir.join(f), b"").unwrap();
    }
    let table = dir.join("m6a_sites.parquet").to_string_lossy().into_owned();
    let record = dir.join("pipeline_summary.json");
    for given in [dir.to_path_buf(), record] {
        let args = parse(&[&given.to_string_lossy(), "s.bam"]);
        let (sites, bams) = args.inputs().unwrap().unwrap();
        assert_eq!(sites, boxed(&[&table]));
        assert_eq!(bams, boxed(&["s.bam"]));
        assert!(!args.batch_process, "the view by default");
    }
    let args = parse(&[&dir.to_string_lossy(), "--modality", "atoi"]);
    assert!(args.inputs().is_err(), "no A-to-I sites there");
}

#[test]
fn no_sites_in_batch_says_what_to_name() {
    let err = parse(&["a.bam", "--batch-process"]).inputs().unwrap_err();
    assert!(err.to_string().contains("faba output directory"), "{err}");
}

#[test]
fn sam_and_cram_are_reads_too() {
    assert!(is_bam("a.bam") && is_bam("b.CRAM") && is_bam("c.sam"));
    assert!(!is_bam("m6a_sites.parquet") && !is_bam("out/"));
}

#[test]
fn tables_called_against_different_genomes_need_one_named() {
    let tmp = tempfile::tempdir().unwrap();
    // Two output folders, each with its sites and a run record naming a genome.
    let table = |dir: &str, genome: &str| -> Box<str> {
        let d = tmp.path().join(dir);
        std::fs::create_dir_all(&d).unwrap();
        let fa = tmp.path().join(genome);
        std::fs::write(&fa, "").unwrap();
        let record = serde_json::json!({"inputs": {"genome": fa.to_str().unwrap()}});
        std::fs::write(d.join("pipeline_summary.json"), record.to_string()).unwrap();
        let sites = d.join("m6a_sites.parquet");
        std::fs::write(&sites, "").unwrap();
        sites.to_string_lossy().into()
    };
    let (a, b, c) = (
        table("a", "hg.fa"),
        table("b", "hg.fa"),
        table("c", "mm.fa"),
    );
    let same = genome_of(None, &[a.clone(), b]).unwrap().unwrap();
    assert!(same.ends_with("hg.fa"));
    let err = genome_of(None, &[a.clone(), c.clone()]).unwrap_err();
    assert!(err.to_string().contains("different genomes"), "{err}");
    assert_eq!(
        &*genome_of(Some("x.fa"), &[a, c]).unwrap().unwrap(),
        "x.fa",
        "--genome wins"
    );
}

#[test]
fn no_genome_recorded_is_none_to_ask_for() {
    let tmp = tempfile::tempdir().unwrap();
    let sites = tmp.path().join("m6a_sites.parquet");
    std::fs::write(&sites, "").unwrap();
    let sites: Box<str> = sites.to_string_lossy().into();
    assert!(genome_of(None, &[sites]).unwrap().is_none());
}
