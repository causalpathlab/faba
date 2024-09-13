pub mod run;

use clap::Args;

#[derive(Args)]
pub struct DepthArgs {
    /// GFF file
    #[arg(short, long)]
    gff: Box<str>,
    /// Output file header
    #[arg(short, long)]
    output: Box<str>,
}
