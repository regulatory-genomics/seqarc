use spring_rs::{decompress, gzip_reader, DecompressionOptions, FastxReader, Format, Result};
use std::fs::File;
use std::io::BufReader;

fn last_record(path: &str) -> Result<spring_rs::ReadRecord> {
    FastxReader::new(
        BufReader::new(gzip_reader(BufReader::new(File::open(path)?))),
        Format::Fastq,
    )
    .last()
    .ok_or_else(|| spring_rs::Error::InvalidInput("empty FASTQ file".into()))?
}

fn main() -> Result<()> {
    let last = 6_989_780;
    let archive = std::env::var("SPRING_RS_ARCHIVE")
        .unwrap_or_else(|_| "data/20250816gan-2.spring-rs".into());
    let decoded = decompress(
        BufReader::new(File::open(archive)?),
        &DecompressionOptions {
            original_range: Some(last..=last),
            ..DecompressionOptions::default()
        },
    )?;
    assert_eq!(
        decoded.streams[0][0],
        last_record("data/20250816gan-2_R1.fastq.gz")?
    );
    assert_eq!(
        decoded.streams[1][0],
        last_record("data/20250816gan-2_R2.fastq.gz")?
    );
    println!("validated final paired record");
    Ok(())
}
