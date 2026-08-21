use spring_rs::{decompress_aligned_paired, DecompressionOptions};
use std::fs::File;
use std::io::BufReader;

fn main() -> spring_rs::Result<()> {
    let archive = std::env::args()
        .nth(1)
        .ok_or_else(|| spring_rs::Error::InvalidInput("archive path is required".into()))?;
    let records = decompress_aligned_paired(
        BufReader::with_capacity(16 * 1024 * 1024, File::open(archive)?),
        &DecompressionOptions::default(),
    )?;
    println!(
        "decoded {} streams and {} records",
        records.streams.len(),
        records.streams.iter().map(Vec::len).sum::<usize>(),
    );
    Ok(())
}
