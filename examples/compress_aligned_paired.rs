use spring_rs::{compress_aligned_paired_fastq_files_with_index, AlignmentCompressionOptions};
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
    compress_aligned_paired_fastq_files_with_index(
        reference,
        index_prefix,
        input_1,
        input_2,
        output,
        &AlignmentCompressionOptions::default(),
    )?;
    Ok(())
}
