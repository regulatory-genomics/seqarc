//! A Rust-native, library-only sequencing read compressor.
//!
//! SeqArc archives are self-describing and support range decoding
//! without a command-line interface.

mod aligners;
mod alignment_archive;
mod fastx;
pub mod minibwa;
mod packed_sequence;
mod quality;
mod qvz;
mod rustar_index;

/// mimalloc reduces allocator overhead in the decode path, where chunks of
/// per-record allocations are freed while large section buffers churn.
#[global_allocator]
static GLOBAL_ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

use std::io::Read;
use std::ops::RangeInclusive;
use std::path::Path;
use std::time::Duration;

pub use aligners::{
    align_dna_batch, build_alignment_index, build_genome_index, build_genome_index_with_options,
    AlignerKind, GenomeAlignment, GenomeIndexBuildOptions,
};
pub use alignment_archive::{
    read_archive_info, read_archive_info_with_breakdown, verify_aligned_archive, AlignedRecordIter,
    ArchiveInfo, IndexStreamSubsectionSizes, MateSubsectionSizes, RecordSubsectionSizes,
    SectionBreakdown, SectionSizes, VerifyReport,
};
pub use fastx::{gzip_reader, gzip_writer, open_fastx_reader, FastxReader, FastxWriter, Format};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    InvalidFormat(String),
    InvalidInput(String),
    Checksum(String),
}

impl Error {
    /// Returns true when the error was caused by a failed checksum comparison.
    pub fn is_checksum(&self) -> bool {
        matches!(self, Self::Checksum(_))
    }
}

impl From<std::io::Error> for Error {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(error) => write!(f, "I/O error: {error}"),
            Self::InvalidFormat(error) => write!(f, "invalid archive: {error}"),
            Self::InvalidInput(error) => write!(f, "invalid input: {error}"),
            Self::Checksum(error) => write!(f, "checksum verification failed: {error}"),
        }
    }
}

impl std::error::Error for Error {}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadRecord {
    pub id: Vec<u8>,
    pub sequence: Vec<u8>,
    pub quality: Option<Vec<u8>>,
}

impl ReadRecord {
    pub fn new(
        id: impl Into<Vec<u8>>,
        sequence: impl Into<Vec<u8>>,
        quality: Option<Vec<u8>>,
    ) -> Self {
        Self {
            id: id.into(),
            sequence: sequence.into(),
            quality,
        }
    }

    pub(crate) fn validate(&self, preserve_quality: bool, allow_empty: bool) -> Result<()> {
        if self.sequence.is_empty() && !allow_empty {
            return Err(Error::InvalidInput("read sequence cannot be empty".into()));
        }
        if preserve_quality {
            let quality = self
                .quality
                .as_ref()
                .ok_or_else(|| Error::InvalidInput("quality is required".into()))?;
            if quality.len() != self.sequence.len() {
                return Err(Error::InvalidInput(
                    "quality length does not match sequence length".into(),
                ));
            }
            if !quality.iter().all(|&byte| (33..=126).contains(&byte)) {
                return Err(Error::InvalidInput(
                    "quality scores must be printable Phred+33 ASCII (bytes 33..=126)".into(),
                ));
            }
        }
        Ok(())
    }
}

pub(crate) fn apply_quality(record: &mut ReadRecord, mode: &QualityMode) {
    if let Some(quality) = record.quality.as_mut() {
        match mode {
            QualityMode::Lossless => {}
            QualityMode::IlluminaBinning => {
                for q in quality {
                    *q = match *q {
                        b'!'..=b'"' => b'!',
                        b'#'..=b'*' => b'\'',
                        b'+'..=b'4' => b'0',
                        b'5'..=b'9' => b'7',
                        b':'..=b'>' => b'<',
                        b'?'..=b'C' => b'B',
                        b'D'..=b'H' => b'F',
                        _ => b'I',
                    };
                }
            }
            QualityMode::Binary {
                threshold,
                high,
                low,
            } => {
                for q in quality {
                    *q = if q.saturating_sub(33) >= *threshold {
                        high.saturating_add(33)
                    } else {
                        low.saturating_add(33)
                    };
                }
            }
            // QVZ is a batch transform over a block of reads; it is applied
            // when the chunk quality stream is assembled, not per record.
            QualityMode::Qvz { .. } => {}
        }
    }
}

/// Truncate a read identifier to at most `max_chars` Unicode scalar values.
pub fn truncate_read_id(record: &mut ReadRecord, max_chars: usize) -> Result<()> {
    let id = std::str::from_utf8(&record.id)
        .map_err(|_| Error::InvalidInput("read identifier is not valid UTF-8".into()))?;
    let end = id
        .char_indices()
        .nth(max_chars)
        .map_or(id.len(), |(index, _)| index);
    record.id.truncate(end);
    Ok(())
}

/// Replace the read identifier with its `field`-th whitespace-delimited field.
///
/// Fields are indexed from zero, and runs of consecutive ASCII whitespace act
/// as one separator. Fails when the identifier has fewer fields than requested.
pub fn select_id_field(record: &mut ReadRecord, field: usize) -> Result<()> {
    let id = std::str::from_utf8(&record.id)
        .map_err(|_| Error::InvalidInput("read identifier is not valid UTF-8".into()))?;
    let selected = id.split_ascii_whitespace().nth(field).ok_or_else(|| {
        Error::InvalidInput(format!(
            "identifier field {field} is missing from read identifier {id:?}"
        ))
    })?;
    record.id = selected.as_bytes().to_vec();
    Ok(())
}

/// Return the identifier proper: the first whitespace-delimited token of the
/// read identifier, excluding any trailing description.
pub fn identifier_token(id: &[u8]) -> Result<&[u8]> {
    let id = std::str::from_utf8(id)
        .map_err(|_| Error::InvalidInput("read identifier is not valid UTF-8".into()))?;
    Ok(id
        .split_ascii_whitespace()
        .next()
        .unwrap_or_default()
        .as_bytes())
}

#[derive(Clone, Debug, PartialEq)]
pub enum QualityMode {
    Lossless,
    IlluminaBinning,
    Binary { threshold: u8, high: u8, low: u8 },
    /// QVZ rate-targeted lossy quantization (ported from Spring). `ratio` is
    /// the target entropy in bits per quality symbol.
    Qvz { ratio: f64 },
}

#[derive(Clone, Debug)]
pub struct AlignmentCompressionOptions {
    pub aligner: AlignerKind,
    pub chunk_records: usize,
    pub threads: usize,
    pub temp_dir: Option<std::path::PathBuf>,
    pub preserve_ids: bool,
    pub id_prefix_chars: Option<usize>,
    pub id_field: Option<usize>,
    pub preserve_quality: bool,
    pub quality_mode: QualityMode,
    pub zstd_level: i32,
    pub embedded_index_ranges: Vec<EmbeddedIndexRange>,
    pub index_read_files: Vec<std::path::PathBuf>,
    pub rustar_gtf: Option<std::path::PathBuf>,
    pub rustar_sjdb_overhang: u32,
    pub rustar_sa_sparse_d: u32,
    pub rustar_sa_index_nbases: u32,
    pub rustar_chr_bin_nbits: u32,
}

impl Default for AlignmentCompressionOptions {
    fn default() -> Self {
        Self {
            aligner: AlignerKind::Minibwa,
            chunk_records: 250_000,
            threads: 16,
            temp_dir: None,
            preserve_ids: true,
            id_prefix_chars: None,
            id_field: None,
            preserve_quality: true,
            quality_mode: QualityMode::Lossless,
            zstd_level: 19,
            embedded_index_ranges: Vec::new(),
            index_read_files: Vec::new(),
            rustar_gtf: None,
            rustar_sjdb_overhang: 100,
            rustar_sa_sparse_d: 1,
            rustar_sa_index_nbases: 14,
            rustar_chr_bin_nbits: 18,
        }
    }
}

#[derive(Clone, Debug)]
pub struct DecompressionOptions {
    pub original_range: Option<RangeInclusive<u64>>,
    pub output_order: OutputOrder,
    pub temp_dir: Option<std::path::PathBuf>,
    pub threads: usize,
    /// Skip chunks that fail their checksum instead of aborting the decode.
    pub best_effort: bool,
    /// Stop after emitting this many records; `None` decodes the whole archive.
    pub max_records: Option<u64>,
}

impl Default for DecompressionOptions {
    fn default() -> Self {
        Self {
            original_range: None,
            output_order: OutputOrder::Original,
            temp_dir: None,
            threads: 1,
            best_effort: false,
            max_records: None,
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OutputOrder {
    #[default]
    Original,
    Compression,
}

pub trait PairedRecordSink {
    fn write_pair(&mut self, original_id: u64, first: ReadRecord, second: ReadRecord)
        -> Result<()>;

    fn finish(&mut self) -> Result<()>;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignmentLayout {
    SingleEnd,
    PairedEnd,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmbeddedIndexMate {
    R1,
    R2,
}

impl EmbeddedIndexMate {
    pub(crate) fn stream(self) -> usize {
        match self {
            Self::R1 => 0,
            Self::R2 => 1,
        }
    }
}

/// A zero-based, half-open range extracted from a biological read before mapping.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EmbeddedIndexRange {
    pub mate: EmbeddedIndexMate,
    pub start: usize,
    pub end: usize,
}

impl EmbeddedIndexRange {
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

impl AlignmentLayout {
    pub fn streams(self) -> usize {
        match self {
            Self::SingleEnd => 1,
            Self::PairedEnd => 2,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SequenceRecordGroup {
    Single(ReadRecord),
    Paired(ReadRecord, ReadRecord),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignedRecordGroup {
    pub sequence: SequenceRecordGroup,
    pub index_reads: Vec<ReadRecord>,
}

pub trait AlignedRecordSink {
    fn write_group(&mut self, original_id: u64, group: AlignedRecordGroup) -> Result<()>;

    fn finish(&mut self) -> Result<()>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignedDecompressionMetadata {
    pub layout: AlignmentLayout,
    pub streams: usize,
    pub index_streams: usize,
    pub records: u64,
    pub output_order: OutputOrder,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DecompressionMetadata {
    pub streams: usize,
    pub pairs: u64,
    pub output_order: OutputOrder,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressionMetadata {
    pub streams: usize,
    pub records: u64,
    pub preserved_ids: bool,
    pub preserved_quality: bool,
    pub preserved_order: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlignedCompressionProfile {
    pub total: Duration,
    pub index: Duration,
    pub map_and_spool: Duration,
    pub reference: Duration,
    pub chunks: Duration,
    pub records: u64,
    pub chunk_count: u64,
    pub spool_bytes: u64,
    pub max_pending_chunks: u64,
    pub candidate_hits: u64,
    pub mapped_reads: u64,
    pub raw_reads: u64,
    pub selected_non_first_hits: u64,
    pub rejected_alignments: u64,
    pub reference_bases: u64,
    pub reference_raw_bytes: u64,
    pub reference_compressed_bytes: u64,
    pub ids_raw_bytes: u64,
    pub qualities_raw_bytes: u64,
    pub records_raw_bytes: u64,
    pub ids_compressed_bytes: u64,
    pub qualities_compressed_bytes: u64,
    pub records_compressed_bytes: u64,
    pub indexes_raw_bytes: u64,
    pub indexes_compressed_bytes: u64,
    pub chunks_compressed_bytes: u64,
    pub archive_bytes: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArchiveRecords {
    pub streams: Vec<Vec<ReadRecord>>,
}

/// Compress paired plain, gzip, or zstd FASTQ files against a genome reference.
///
/// Mapped reads store their reference offset, strand, CIGAR, and sequence noise.
/// Reads without a selected alignment are stored using the packed/raw fallback.
pub fn compress_aligned_paired_fastq_files(
    genome_fasta: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    compress_aligned_paired_fastq_files_profiled(genome_fasta, input_1, input_2, output, options)
        .map(|(metadata, _)| metadata)
}

/// Compress paired plain, gzip, or zstd FASTQ files and report archive details.
pub fn compress_aligned_paired_fastq_files_profiled(
    genome_fasta: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    alignment_archive::compress_paired_fastq_files(genome_fasta, input_1, input_2, output, options)
}

/// Compress paired plain, gzip, or zstd FASTQ files using an existing index.
pub fn compress_aligned_paired_fastq_files_with_index(
    index_file: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    compress_aligned_paired_fastq_files_with_index_profiled(
        index_file, input_1, input_2, output, options,
    )
    .map(|(metadata, _)| metadata)
}

/// Compress using an existing SeqArc alignment index and report archive size details.
pub fn compress_aligned_paired_fastq_files_with_index_profiled(
    index_file: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    alignment_archive::compress_paired_fastq_files_with_index(
        index_file, input_1, input_2, output, options,
    )
}

pub fn compress_aligned_single_fastq_file(
    genome_fasta: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    compress_aligned_single_fastq_file_profiled(genome_fasta, input, output, options)
        .map(|(metadata, _)| metadata)
}

pub fn compress_aligned_single_fastq_file_profiled(
    genome_fasta: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    alignment_archive::compress_single_fastq_file(genome_fasta, input, output, options)
}

pub fn compress_aligned_single_fastq_file_with_index(
    index_file: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    compress_aligned_single_fastq_file_with_index_profiled(index_file, input, output, options)
        .map(|(metadata, _)| metadata)
}

pub fn compress_aligned_single_fastq_file_with_index_profiled(
    index_file: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    alignment_archive::compress_single_fastq_file_with_index(index_file, input, output, options)
}

pub(crate) fn temporary_directory(root: Option<&Path>, prefix: &str) -> Result<tempfile::TempDir> {
    match root {
        Some(root) => {
            std::fs::create_dir_all(root)?;
            tempfile::Builder::new()
                .prefix(prefix)
                .tempdir_in(root)
                .map_err(Into::into)
        }
        None => tempfile::Builder::new()
            .prefix(prefix)
            .tempdir()
            .map_err(Into::into),
    }
}

/// Decompress an alignment archive into paired records.
pub fn decompress_aligned_paired<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    alignment_archive::decompress_paired(input, options)
}

pub fn decompress_aligned_single<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    alignment_archive::decompress_single(input, options)
}

/// Lazily decode aligned records in original logical-record order.
///
/// Up to `options.threads` chunks are decoded concurrently. At most twice that
/// many chunks are resident across active, completed, and currently yielded
/// work, applying backpressure when downstream consumption is slower. Exhaust
/// the iterator to validate all chunks and confirm that no bytes follow the
/// declared archive.
pub fn decompress_aligned_iter<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<AlignedRecordIter<R>> {
    alignment_archive::decompress_iter(input, options)
}

pub fn decompress_aligned_into<R: Read, S: AlignedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<AlignedDecompressionMetadata> {
    alignment_archive::decompress_into(input, options, sink)
}

pub fn decompress_aligned_paired_into<R: Read, S: PairedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<DecompressionMetadata> {
    alignment_archive::decompress_paired_into(input, options, sink)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::{BufReader, Cursor, Write};

    #[test]
    fn identifier_prefix_truncation_preserves_utf8_boundaries() {
        let mut record = ReadRecord::new("réad/1 suffix", b"ACGT".to_vec(), None);
        truncate_read_id(&mut record, 2).unwrap();
        assert_eq!(record.id, "ré".as_bytes());
        truncate_read_id(&mut record, 0).unwrap();
        assert!(record.id.is_empty());
    }

    #[test]
    fn identifier_field_selection_keeps_selected_whitespace_field() {
        let mut record = ReadRecord::new("READID 1:N:0:ATCG+CCAG", b"ACGT".to_vec(), None);
        select_id_field(&mut record, 1).unwrap();
        assert_eq!(record.id, b"1:N:0:ATCG+CCAG");
        let mut record = ReadRecord::new("READID 1:N:0:ATCG+CCAG", b"ACGT".to_vec(), None);
        select_id_field(&mut record, 0).unwrap();
        assert_eq!(record.id, b"READID");
    }

    #[test]
    fn identifier_field_selection_collapses_whitespace_runs() {
        let mut record = ReadRecord::new("a \t b\t\tc", b"ACGT".to_vec(), None);
        select_id_field(&mut record, 1).unwrap();
        assert_eq!(record.id, b"b");
    }

    #[test]
    fn identifier_field_selection_rejects_missing_fields() {
        let mut record = ReadRecord::new("only", b"ACGT".to_vec(), None);
        assert!(select_id_field(&mut record, 1).is_err());
        let mut invalid = ReadRecord::new(vec![0xff], b"ACGT".to_vec(), None);
        assert!(select_id_field(&mut invalid, 0).is_err());
    }

    #[test]
    fn identifier_token_returns_first_whitespace_field() {
        assert_eq!(
            identifier_token(b"read0 1:N:0:ATCG+CCAG").unwrap(),
            b"read0"
        );
        assert_eq!(identifier_token(b"read0\tcomment").unwrap(), b"read0");
        assert_eq!(identifier_token(b"read0").unwrap(), b"read0");
        assert_eq!(identifier_token(b"  spaced  id ").unwrap(), b"spaced");
        assert!(identifier_token(&[0xff]).is_err());
    }

    #[test]
    fn validate_rejects_quality_bytes_outside_printable_ascii() {
        let mut record = ReadRecord::new(
            "read",
            b"ACGTACGT".to_vec(),
            Some((33u8..=126).collect::<Vec<_>>()[..8].to_vec()),
        );
        assert!(record.validate(true, false).is_ok());

        // Control character below the Phred+33 floor.
        record.quality = Some(b"A\x01C\x01ACGT\x01".to_vec());
        assert!(record.validate(true, false).is_err());

        // DEL (127) and every byte above it are outside the printable range.
        record.quality = Some(b"A\x7FC\x7FACGT\x7F".to_vec());
        assert!(record.validate(true, false).is_err());
        record.quality = Some(b"AAAAACGT".to_vec());
        assert!(record.validate(true, false).is_ok());

        // FASTA-style records without quality remain valid.
        let record = ReadRecord::new("read", b"ACGT".to_vec(), None);
        assert!(record.validate(false, false).is_ok());
    }

    #[test]
    fn illumina_binning_matches_expected_table() {
        let mut record =
            ReadRecord::new("read", b"ACGTACGT".to_vec(), Some((33u8..=127).collect()));
        apply_quality(&mut record, &QualityMode::IlluminaBinning);
        let expected = (33u8..=127)
            .map(|quality| match quality {
                33..=34 => 33,
                35..=42 => 39,
                43..=52 => 48,
                53..=57 => 55,
                58..=62 => 60,
                63..=67 => 66,
                68..=72 => 70,
                _ => 73,
            })
            .collect::<Vec<_>>();
        assert_eq!(record.quality, Some(expected));
    }

    #[test]
    fn fastq_reader_and_writer_round_trip() {
        let input = b"@r1\nACGT\n+\nIIII\n";
        let records: Vec<_> = FastxReader::new(BufReader::new(&input[..]), Format::Fastq)
            .collect::<Result<Vec<_>>>()
            .unwrap();
        let mut output = Vec::new();
        let mut writer = FastxWriter::new(&mut output, Format::Fastq);
        for record in &records {
            writer.write_record(record).unwrap();
        }
        writer.finish().unwrap();
        assert_eq!(output, input);
    }

    #[test]
    fn fasta_reader_produces_records_without_quality() {
        let input = b">r1\nAC\nGT\n>r2\nNN\n";
        let records: Vec<_> = FastxReader::new(BufReader::new(&input[..]), Format::Fasta)
            .collect::<Result<Vec<_>>>()
            .unwrap();
        assert_eq!(records[0].sequence, b"ACGT");
        assert_eq!(records[1].quality, None);
    }

    #[test]
    fn aligned_paired_archive_round_trips_mapped_and_unmapped_reads() {
        let directory = tempfile::tempdir().unwrap();
        let genome = directory.path().join("genome.fa");
        let input_1 = directory.path().join("reads.1.fastq.gz");
        let input_2 = directory.path().join("reads.2.fastq.gz");
        let output = directory.path().join("reads.aln");
        let reference =
            b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
        std::fs::write(&genome, [b">chr1\n".as_slice(), reference, b"\n"].concat()).unwrap();
        let mut forward = reference[20..60].to_vec();
        forward[5] = match forward[5] {
            b'A' => b'C',
            _ => b'A',
        };
        let reverse = reference[20..60]
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => *base,
            })
            .collect::<Vec<_>>();
        let unknown = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let first_fastq = format!(
            "@pair0/1\n{}\n+\n{}\n@pair1/1\n{}\n+\n{}\n",
            String::from_utf8_lossy(&forward),
            "I".repeat(forward.len()),
            String::from_utf8_lossy(unknown),
            "#".repeat(unknown.len()),
        );
        let second_fastq = format!(
            "@pair0/2\n{}\n+\n{}\n@pair1/2\n{}\n+\n{}\n",
            String::from_utf8_lossy(&reverse),
            "J".repeat(reverse.len()),
            String::from_utf8_lossy(unknown),
            "$".repeat(unknown.len()),
        );
        for (path, contents) in [(&input_1, first_fastq), (&input_2, second_fastq)] {
            let encoder = gzip_writer(File::create(path).unwrap(), 1);
            let mut writer = encoder;
            writer.write_all(contents.as_bytes()).unwrap();
            writer.finish().unwrap();
        }

        let alignment_options = AlignmentCompressionOptions {
            chunk_records: 1,
            ..AlignmentCompressionOptions::default()
        };
        let (metadata, profile) = compress_aligned_paired_fastq_files_profiled(
            &genome,
            &input_1,
            &input_2,
            &output,
            &alignment_options,
        )
        .unwrap();
        assert_eq!(metadata.records, 2);
        assert_eq!(profile.records, 2);
        assert_eq!(profile.chunk_count, 2);
        assert_eq!(profile.mapped_reads + profile.raw_reads, 4);
        assert!(profile.spool_bytes > 0);
        assert!(profile.reference_bases > 0);
        assert!(profile.reference_compressed_bytes > 0);
        assert_eq!(
            profile.archive_bytes,
            46 + profile.reference_compressed_bytes
                + 8
                + 48 * profile.chunk_count
                + profile.chunks_compressed_bytes
        );
        let mut archive = Vec::new();
        File::open(&output)
            .unwrap()
            .read_to_end(&mut archive)
            .unwrap();
        let decoded = decompress_aligned_paired(
            Cursor::new(archive.clone()),
            &DecompressionOptions::default(),
        )
        .unwrap();

        assert_eq!(decoded.streams[0].len(), 2);
        assert_eq!(decoded.streams[1].len(), 2);
        assert_eq!(decoded.streams[0][0].sequence, forward);
        assert_eq!(decoded.streams[1][0].sequence, reverse);

        // The record limit applies to the paired collection decoder too.
        let limited = decompress_aligned_paired(
            Cursor::new(archive.clone()),
            &DecompressionOptions {
                max_records: Some(1),
                ..DecompressionOptions::default()
            },
        )
        .unwrap();
        assert_eq!(limited.streams[0].len(), 1);
        assert_eq!(limited.streams[1].len(), 1);

        struct VecPairSink {
            pairs: Vec<(u64, ReadRecord, ReadRecord)>,
        }

        impl PairedRecordSink for VecPairSink {
            fn write_pair(
                &mut self,
                original_id: u64,
                first: ReadRecord,
                second: ReadRecord,
            ) -> Result<()> {
                self.pairs.push((original_id, first, second));
                Ok(())
            }

            fn finish(&mut self) -> Result<()> {
                Ok(())
            }
        }

        let mut sink = VecPairSink { pairs: Vec::new() };
        decompress_aligned_paired_into(
            Cursor::new(&archive),
            &DecompressionOptions {
                max_records: Some(1),
                ..DecompressionOptions::default()
            },
            &mut sink,
        )
        .unwrap();
        assert_eq!(sink.pairs.len(), 1);
        assert_eq!(decoded.streams[0][1].sequence, unknown);
        assert_eq!(decoded.streams[1][1].sequence, unknown);
        assert_eq!(
            decoded.streams[0][0].quality.as_deref(),
            Some(vec![b'I'; forward.len()].as_slice())
        );
        assert_eq!(decoded.streams[0][0].id, b"pair0/1");

        let parallel_options = DecompressionOptions {
            original_range: Some(1..=1),
            threads: 4,
            ..DecompressionOptions::default()
        };
        let mut iter =
            decompress_aligned_iter(Cursor::new(archive.clone()), &parallel_options).unwrap();
        assert_eq!(iter.layout(), AlignmentLayout::PairedEnd);
        assert_eq!(iter.streams(), 2);
        assert_eq!(iter.total_chunks(), 2);
        assert_eq!(iter.decoded_chunks(), 0);
        let selected = iter.by_ref().collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(iter.decoded_chunks(), 2);
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].0, 1);
        assert!(selected[0].1.index_reads.is_empty());
        let SequenceRecordGroup::Paired(first, second) = &selected[0].1.sequence else {
            panic!("expected paired record")
        };
        assert_eq!(first.sequence, unknown);
        assert_eq!(second.sequence, unknown);

        let mut trailing = archive;
        trailing.push(0);
        let mut iter = decompress_aligned_iter(
            Cursor::new(trailing),
            &DecompressionOptions {
                threads: 2,
                ..DecompressionOptions::default()
            },
        )
        .unwrap();
        assert_eq!(iter.next().unwrap().unwrap().0, 0);
        assert_eq!(iter.next().unwrap().unwrap().0, 1);
        assert!(matches!(
            iter.next(),
            Some(Err(Error::InvalidFormat(message)))
                if message == "trailing aligned archive bytes"
        ));
        assert!(iter.next().is_none());
    }

    #[test]
    fn verify_and_best_effort_detect_and_recover_corrupt_chunks() {
        let directory = tempfile::tempdir().unwrap();
        let genome = directory.path().join("genome.fa");
        let input_1 = directory.path().join("reads.1.fastq.gz");
        let input_2 = directory.path().join("reads.2.fastq.gz");
        let output = directory.path().join("reads.aln");
        let reference =
            b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
        std::fs::write(&genome, [b">chr1\n".as_slice(), reference, b"\n"].concat()).unwrap();
        let mut forward = reference[20..60].to_vec();
        forward[5] = match forward[5] {
            b'A' => b'C',
            _ => b'A',
        };
        let reverse = reference[20..60]
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => *base,
            })
            .collect::<Vec<_>>();
        let unknown = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let first_fastq = format!(
            "@pair0/1\n{}\n+\n{}\n@pair1/1\n{}\n+\n{}\n",
            String::from_utf8_lossy(&forward),
            "I".repeat(forward.len()),
            String::from_utf8_lossy(unknown),
            "#".repeat(unknown.len()),
        );
        let second_fastq = format!(
            "@pair0/2\n{}\n+\n{}\n@pair1/2\n{}\n+\n{}\n",
            String::from_utf8_lossy(&reverse),
            "J".repeat(reverse.len()),
            String::from_utf8_lossy(unknown),
            "$".repeat(unknown.len()),
        );
        for (path, contents) in [(&input_1, first_fastq), (&input_2, second_fastq)] {
            let encoder = gzip_writer(File::create(path).unwrap(), 1);
            let mut writer = encoder;
            writer.write_all(contents.as_bytes()).unwrap();
            writer.finish().unwrap();
        }
        let alignment_options = AlignmentCompressionOptions {
            chunk_records: 1,
            ..AlignmentCompressionOptions::default()
        };
        compress_aligned_paired_fastq_files_profiled(
            &genome,
            &input_1,
            &input_2,
            &output,
            &alignment_options,
        )
        .unwrap();
        let mut archive = Vec::new();
        File::open(&output)
            .unwrap()
            .read_to_end(&mut archive)
            .unwrap();

        // A healthy archive verifies cleanly.
        let report = verify_aligned_archive(&mut Cursor::new(archive.clone())).unwrap();
        assert!(report.is_ok());
        assert_eq!(report.declared_chunks, 2);
        assert_eq!(report.scanned_chunks, 2);
        assert_eq!(report.records, 2);

        // Flip the final payload byte, which belongs to the last chunk.
        let mut damaged = archive.clone();
        let last = damaged.len() - 1;
        damaged[last] ^= 0xff;
        let report = verify_aligned_archive(&mut Cursor::new(damaged.clone())).unwrap();
        assert!(!report.is_ok());
        assert_eq!(report.bad_chunks, 1);
        assert_eq!(report.first_bad_chunk, Some(1));
        assert!(!report.truncated);
        assert_eq!(report.records, 1);

        // Strict decoding refuses the damaged archive.
        let strict = decompress_aligned_paired(
            Cursor::new(damaged.clone()),
            &DecompressionOptions::default(),
        );
        assert!(strict.is_err());

        // Best-effort decoding skips the damaged chunk and recovers the rest.
        let mut iter = decompress_aligned_iter(
            Cursor::new(damaged),
            &DecompressionOptions {
                best_effort: true,
                ..DecompressionOptions::default()
            },
        )
        .unwrap();
        let recovered = iter.by_ref().collect::<Result<Vec<_>>>().unwrap();
        assert_eq!(recovered.len(), 1);
        assert_eq!(iter.skipped_chunks(), 1);
        assert_eq!(iter.skipped_records(), 1);
        assert!(iter.next().is_none());
        let SequenceRecordGroup::Paired(first, second) = &recovered[0].1.sequence else {
            panic!("expected paired record")
        };
        assert_eq!(first.sequence, forward);
        assert_eq!(second.sequence, reverse);
    }

    #[test]
    fn gzip_reader_handles_concatenated_members() {
        use std::io::Write;

        let mut gzip = Vec::new();
        for bytes in [b"first".as_slice(), b"second".as_slice()] {
            let mut encoder = gzip_writer(Vec::new(), 6);
            encoder.write_all(bytes).unwrap();
            gzip.extend(encoder.finish().unwrap());
        }
        let mut decoded = String::new();
        gzip_reader(&gzip[..]).read_to_string(&mut decoded).unwrap();
        assert_eq!(decoded, "firstsecond");
    }
}
