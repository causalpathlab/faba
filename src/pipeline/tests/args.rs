use crate::pipeline::args::PipelineArgs;
use clap::Parser;

#[derive(Parser)]
struct Cli {
    #[command(flatten)]
    run: PipelineArgs,
}

fn parse(words: &[&str]) -> PipelineArgs {
    Cli::try_parse_from(std::iter::once("faba").chain(words.iter().copied()))
        .unwrap()
        .run
}

#[test]
fn nothing_is_required_until_batch_process() {
    let a = parse(&[]);
    assert!(a.bam_files.is_empty() && a.gff_file.is_none());
    assert!(!a.batch_process);
}

#[test]
fn batch_process_names_what_is_missing() {
    let a = parse(&["--batch-process", "a.bam", "-g", "g.gff"]);
    let e = a.check_batch().unwrap_err().to_string();
    assert!(
        e.contains("-f/--genome") && e.contains("-o/--output"),
        "{e}"
    );
    assert!(!e.contains("--gff"), "{e}");
    let a = parse(&[
        "--batch-process",
        "a.bam",
        "-g",
        "g.gff",
        "-f",
        "x.fa",
        "-o",
        "out",
    ]);
    a.check_batch().unwrap();
    assert_eq!((a.gff(), a.genome(), a.out()), ("g.gff", "x.fa", "out"));
}

#[test]
fn common_flags_carry_the_common_heading() {
    use clap::CommandFactory;
    let cmd = Cli::command();
    let heading = |long: &str| {
        cmd.get_arguments()
            .find(|a| a.get_long() == Some(long))
            .and_then(|a| a.get_help_heading())
            .map(str::to_string)
    };
    for long in [
        "max-threads",
        "backend",
        "no-zip",
        "cell-barcode-tag",
        "umi-tag",
        "min-mapping-quality",
    ] {
        assert_eq!(heading(long).as_deref(), Some("Common"), "{long}");
    }
}
