mod aggregate;
mod util;

use aggregate::run_aggregate;

use anyhow;
// use anyhow::{self};
use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(version, about, long_about=None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    commands: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// aggregate putative editing sites
    Aggregate(AggregateArgs),
    /// sifting through potential sites
    Sift,
}

#[derive(Args)]
struct SiftArgs {
    /// GFF file
    #[arg(short, long)]
    gff: Box<str>,
}

#[derive(Args)]
struct AggregateArgs {
    /// foreground BAM file
    #[arg(short, long)]
    fg_bam: Box<str>,

    /// background BAM file
    #[arg(short, long)]
    bg_bam: Box<str>,

    /// foreground BAI file (default: <FG_BAM>.bai)
    #[arg(long)]
    fg_bai: Option<Box<str>>,

    /// background BAI file (default: <BG_BAM>.bai)
    #[arg(long)]
    bg_bai: Option<Box<str>>,

    /// number of threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// block size (default: 10000)
    #[arg(long)]
    bsize: Option<usize>,

    /// output file header
    #[arg(short, long)]
    output: Option<Box<str>>,
}

/// main CLI for FAVA
///
fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match &cli.commands {
        Commands::Aggregate(args) => {
            //
            run_aggregate(args)?;
        }
        _ => {
            todo!("sifting");
        }
    }
    Ok(())
}
