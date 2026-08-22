use spring_rs::{
    compress_aligned_paired_fastq_files_with_index_profiled, AlignmentCompressionOptions,
};
use std::env;
use std::path::PathBuf;

fn main() -> spring_rs::Result<()> {
    let mut args = env::args_os().skip(1);
    let reference = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing reference FASTA".into()))?,
    );
    let index_prefix = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing index prefix".into()))?,
    );
    let input_1 = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing mate-1 FASTQ".into()))?,
    );
    let input_2 = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing mate-2 FASTQ".into()))?,
    );
    let output = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing output archive".into()))?,
    );
    if args.next().is_some() {
        return Err(spring_rs::Error::InvalidInput(
            "expected reference, index prefix, two FASTQ files, and output".into(),
        ));
    }
    let (_, profile) = compress_aligned_paired_fastq_files_with_index_profiled(
        reference,
        index_prefix,
        input_1,
        input_2,
        output,
        &AlignmentCompressionOptions::default(),
    )?;
    eprintln!("archive bytes: {}", profile.archive_bytes);
    eprintln!(
        "reference: {} bases, {} bytes compressed",
        profile.reference_bases, profile.reference_compressed_bytes
    );
    eprintln!(
        "chunk streams raw: ids {}, qualities {}, records {} bytes",
        profile.ids_raw_bytes, profile.qualities_raw_bytes, profile.records_raw_bytes
    );
    eprintln!(
        "chunk streams compressed: ids {}, qualities {}, records {} bytes",
        profile.ids_compressed_bytes,
        profile.qualities_compressed_bytes,
        profile.records_compressed_bytes
    );
    eprintln!(
        "mapping: {} mapped, {} raw, {} candidate hits, {} non-first selections, {} rejected",
        profile.mapped_reads,
        profile.raw_reads,
        profile.candidate_hits,
        profile.selected_non_first_hits,
        profile.rejected_alignments
    );
    Ok(())
}
