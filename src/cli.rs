#[cfg(feature = "cli")]
use clap::Parser;
use crate::{convert_file_to_ipsa, detect_format_from_path};
use std::path::PathBuf;

#[cfg(feature = "cli")]
#[derive(Parser)]
#[command(author, version, about = "Universal glycoproteomics PSM converter")]
pub struct Args {
    #[arg(short, long, help = "Input file path")]
    pub input: PathBuf,

    #[arg(short, long, help = "Output directory")]
    pub output: PathBuf,

    #[arg(long, help = "Show detection information")]
    pub detect_only: bool,
}

#[cfg(feature = "cli")]
pub fn run() -> anyhow::Result<()> {
    let args = Args::parse();

    if args.detect_only {
        let format = detect_format_from_path(&args.input)?;
        println!("Detected format: {}", format);
        return Ok(());
    }

    let stats = convert_file_to_ipsa(&args.input, &args.output)?;

    println!("✅ Conversion complete:");
    println!("   Format: {}", stats.format);
    println!("   Identifications: {}", stats.identifications_written);
    println!("   Modifications: {}", stats.modifications_written);
    println!("   Time: {}ms", stats.processing_time_ms);

    Ok(())
}

#[cfg(feature = "cli")]
pub fn main() -> anyhow::Result<()> {
    run()
}