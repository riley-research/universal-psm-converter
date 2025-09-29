#[cfg(feature = "cli")]
fn main() -> anyhow::Result<()> {
    universal_psm_converter::cli::main()
}

#[cfg(not(feature = "cli"))]
fn main() {
    println!("CLI feature not enabled. Use as library crate.");
}