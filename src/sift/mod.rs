pub mod rules;
pub mod compare;

use clap::Args;

#[derive(Args)]
pub struct CaseControlArgs {
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

// a wrapper for input frequency stat argument
// pub struct FreqStat<'a> {
//     fg: &'a DnaFreq,
//     bg: &'a DnaFreq,
// }
