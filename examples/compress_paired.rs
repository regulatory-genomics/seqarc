use spring_rs::{compress_paired_fastq_files, CompressionOptions, ReorderOptions};

fn main() -> spring_rs::Result<()> {
    let temp_dir = std::env::var("SPRING_RS_WORK_DIR").ok().map(Into::into);
    let reorder = ReorderOptions {
        kmer_size: env_usize("SPRING_RS_KMER", 21),
        window_size: env_usize("SPRING_RS_WINDOW", 31),
        max_bucket_records: env_usize("SPRING_RS_MAX_BUCKET", 4_096),
        temp_dir,
        ..ReorderOptions::default()
    };
    let options = CompressionOptions {
        reorder,
        ..CompressionOptions::default()
    };
    let input_1 = std::env::var("SPRING_RS_INPUT_1")
        .unwrap_or_else(|_| "data/20250816gan-2_R1.fastq.gz".into());
    let input_2 = std::env::var("SPRING_RS_INPUT_2")
        .unwrap_or_else(|_| "data/20250816gan-2_R2.fastq.gz".into());
    let output =
        std::env::var("SPRING_RS_OUTPUT").unwrap_or_else(|_| "data/20250816gan-2.spring-rs".into());
    let start = std::time::Instant::now();
    let metadata = compress_paired_fastq_files(input_1, input_2, output, &options)?;
    println!(
        "{metadata:?} using {} threads in {:.1}s",
        options.threads,
        start.elapsed().as_secs_f64()
    );
    Ok(())
}

fn env_usize(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}
