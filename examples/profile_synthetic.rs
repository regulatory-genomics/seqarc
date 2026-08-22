use seqarc::{compress_reordered_paired_fastq_files_profiled, gzip_writer, CompressionOptions};
use std::fs::File;
use std::io::Write;

fn main() -> seqarc::Result<()> {
    let pairs: usize = std::env::var("SEQARC_PROFILE_PAIRS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(10_000);
    let root = std::env::var("SEQARC_PROFILE_DIR")
        .unwrap_or_else(|_| "/tmp/opencode/seqarc-profile".into());
    std::fs::create_dir_all(&root)?;
    let input_1 = format!("{root}/reads.1.fastq.gz");
    let input_2 = format!("{root}/reads.2.fastq.gz");
    let output = format!("{root}/reads.seqarc");
    let work = format!("{root}/work");
    let mut first = gzip_writer(File::create(&input_1)?, 1);
    let mut second = gzip_writer(File::create(&input_2)?, 1);
    for index in 0..pairs {
        let sequence_1 = sequence(index, 0);
        let sequence_2 = sequence(index, 17);
        writeln!(first, "@read-{index}/1")?;
        first.write_all(&sequence_1)?;
        first.write_all(b"\n+\n")?;
        first.write_all(&vec![b'I'; sequence_1.len()])?;
        first.write_all(b"\n")?;
        writeln!(second, "@read-{index}/2")?;
        second.write_all(&sequence_2)?;
        second.write_all(b"\n+\n")?;
        second.write_all(&vec![b'J'; sequence_2.len()])?;
        second.write_all(b"\n")?;
    }
    first.finish()?;
    second.finish()?;
    let level = std::env::var("SEQARC_ZSTD_LEVEL")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(19);
    let threads = std::env::var("SEQARC_THREADS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(16);
    let options = CompressionOptions {
        zstd_level: level,
        threads,
        ..CompressionOptions::default()
    };
    let (metadata, profile) =
        compress_reordered_paired_fastq_files_profiled(input_1, input_2, output, &options, work)?;
    println!("{metadata:#?}\n{profile:#?}");
    Ok(())
}

fn sequence(index: usize, offset: usize) -> Vec<u8> {
    let mut sequence = (0..150)
        .map(|position| b"ACGT"[(position + offset + index / 100) % 4])
        .collect::<Vec<_>>();
    sequence[index % 150] = b"ACGT"[(index / 7 + offset) % 4];
    sequence
}
