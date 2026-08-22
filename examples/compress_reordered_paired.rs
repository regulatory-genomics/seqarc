use seqarc::{compress_reordered_paired_fastq_files_profiled, CompressionOptions};

fn main() -> seqarc::Result<()> {
    let input_1 =
        std::env::var("SEQARC_INPUT_1").unwrap_or_else(|_| "data/20250816gan-2_R1.fastq.gz".into());
    let input_2 =
        std::env::var("SEQARC_INPUT_2").unwrap_or_else(|_| "data/20250816gan-2_R2.fastq.gz".into());
    let output = std::env::var("SEQARC_OUTPUT")
        .unwrap_or_else(|_| "data/20250816gan-2.consensus.seqarc".into());
    let work_dir = std::env::var("SEQARC_WORK_DIR")
        .unwrap_or_else(|_| "/tmp/opencode/seqarc-consensus-work".into());
    let mut options = CompressionOptions::default();
    options.reorder.max_bucket_records = std::env::var("SEQARC_MAX_BUCKET")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(options.reorder.max_bucket_records);
    options.reorder.kmer_size = std::env::var("SEQARC_KMER")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(options.reorder.kmer_size);
    options.reorder.window_size = std::env::var("SEQARC_WINDOW")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(options.reorder.window_size);
    if std::env::var_os("SEQARC_NO_QUALITY").is_some() {
        options.preserve_quality = false;
    }
    if std::env::var_os("SEQARC_NO_IDS").is_some() {
        options.preserve_ids = false;
    }
    let start = std::time::Instant::now();
    let (metadata, profile) = compress_reordered_paired_fastq_files_profiled(
        input_1, input_2, output, &options, work_dir,
    )?;
    println!(
        "{metadata:?} in {:.1}s (zstd level {}, sort threads {}, max bucket {})",
        start.elapsed().as_secs_f64(),
        options.zstd_level,
        options.reorder.sort_threads,
        options.reorder.max_bucket_records,
    );
    println!("profile: {profile:#?}");
    Ok(())
}
