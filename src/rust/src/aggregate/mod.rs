use anyhow;
use clap::Args;

#[derive(Args)]
pub struct AggArgs {
    /// GFF file
    #[arg(short, long)]
    gff: Box<str>,
    /// Output file header
    #[arg(short, long)]
    output: Box<str>,
}

pub fn run_agg(args: &AggArgs) -> anyhow::Result<()> {
    Ok(())
}
