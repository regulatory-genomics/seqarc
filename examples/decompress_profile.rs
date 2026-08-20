use spring_rs::{decompress, DecompressionOptions};
use std::fs::File;
use std::io::BufReader;

fn main() -> spring_rs::Result<()> {
    let archive = std::env::args()
        .nth(1)
        .ok_or_else(|| spring_rs::Error::InvalidInput("archive path is required".into()))?;
    let started = std::time::Instant::now();
    let records = decompress(
        BufReader::with_capacity(16 * 1024 * 1024, File::open(archive)?),
        &DecompressionOptions::default(),
    )?;
    println!(
        "decoded {} streams and {} records in {:.3}s",
        records.streams.len(),
        records.streams.iter().map(Vec::len).sum::<usize>(),
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}
