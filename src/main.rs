mod aggregate;
mod depth;
mod sift;
mod util;

use anyhow;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(version, about, long_about=None)]
#[command(propagate_version = true)]
///
/// favabean
///
struct Cli {
    #[command(subcommand)]
    commands: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Sift through to collect statistics for putative editing sites
    Compare(sift::CaseControlArgs),
    /// Aggregate for global patterning
    Aggregate(aggregate::AggArgs),
    /// Depth
    Depth(depth::DepthArgs),
}

/// main CLI for FAVA
///
fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match &cli.commands {
        Commands::Compare(args) => {
            sift::compare::run(args)?;
        }
        Commands::Aggregate(args) => {
            //
        }

        Commands::Depth(args) => {
            //
        }
    }
    Ok(())
}
