//! SeqArc is a Rust archiver for sequencing reads.
//!
//! Instead of compressing raw FASTQ bytes, reads are aligned against a
//! reference genome and stored as a compact alignment-based representation,
//! compressed with Zstandard. Archives reconstruct the original reads exactly
//! while remaining far smaller than gzip-compressed FASTQ.
//!
//! Archives are self-describing: the global header and every chunk carry
//! CRC32C checksums, and records can be decoded by logical-record range. The
//! `seqarc` command-line binary in this crate exposes the full workflow
//! (index, compress, decompress, validate, verify, info); this library API
//! offers the same operations programmatically.

#![warn(missing_docs)]

/// RNA-seq aligner subset vendored from `scverse/rustar-aligner`.
///
/// See [`crate::aligners`] for the SeqArc-level API; the index container is
/// handled by [`crate::star_index`].
#[allow(missing_docs)]
pub(crate) mod star;

mod aligners;
mod alignment_archive;
mod fastx;
/// Read mapper ported from the C `minibwa` aligner.
///
/// Naming and layout follow the upstream C implementation; see the aligner
/// wrappers in [`crate::aligners`] for the SeqArc-level API.
#[allow(missing_docs)]
pub mod minibwa;
mod packed_sequence;
mod quality;
mod qvz;
mod star_index;

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

/// Library result type.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors returned by SeqArc archive operations.
#[derive(Debug)]
pub enum Error {
    /// Filesystem or stream failure.
    Io(std::io::Error),
    /// The archive bytes do not follow the expected format.
    InvalidFormat(String),
    /// The caller-supplied inputs are unusable.
    InvalidInput(String),
    /// A stored checksum did not match the data.
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

/// One decoded (or to-be-encoded) sequencing read.
///
/// Identifiers and sequences are raw bytes so non-UTF-8 archives round-trip.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadRecord {
    /// Read identifier without the leading `@`.
    pub id: Vec<u8>,
    /// Bases as ASCII (A, C, G, T, N, or other IUPAC codes).
    pub sequence: Vec<u8>,
    /// Phred+33 quality scores, present unless the archive drops quality.
    pub quality: Option<Vec<u8>>,
}

impl ReadRecord {
    /// Build a record from raw parts.
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

/// Quality-preservation strategy applied while compressing.
#[derive(Clone, Debug, PartialEq)]
pub enum QualityMode {
    /// Store Phred+33 scores unchanged.
    Lossless,
    /// Bin scores into the standard Illumina 8-bin scale.
    IlluminaBinning,
    /// Threshold scores into two values (`high`/`low`, Phred+33 encoded).
    Binary {
        /// Minimum raw Phred score (byte minus 33) that maps to `high`.
        threshold: u8,
        /// Stored byte for reads at or above `threshold`.
        high: u8,
        /// Stored byte for reads below `threshold`.
        low: u8,
    },
    /// QVZ rate-targeted lossy quantization (ported from Spring). `ratio` is
    /// the target entropy in bits per quality symbol.
    Qvz {
        /// Target entropy in bits per quality symbol.
        ratio: f64,
    },
}

/// Default worker count, derived from the machine's available parallelism.
pub fn default_thread_count() -> usize {
    std::thread::available_parallelism().map_or(1, |count| count.get())
}

/// Options for alignment-based FASTQ compression.
#[derive(Clone, Debug)]
pub struct AlignmentCompressionOptions {
    /// Genome aligner used to map reads.
    pub aligner: AlignerKind,
    /// Records grouped into one independently decodable archive chunk.
    pub chunk_records: usize,
    /// Worker threads for mapping and chunk encoding.
    pub threads: usize,
    /// Scratch directory for temporary files; defaults to the system temp dir.
    pub temp_dir: Option<std::path::PathBuf>,
    /// Store read identifiers instead of dropping them.
    pub preserve_ids: bool,
    /// Keep only the first N Unicode scalar values of each identifier.
    pub id_prefix_chars: Option<usize>,
    /// Keep only the N-th whitespace-delimited field of each identifier.
    pub id_field: Option<usize>,
    /// Store quality values instead of dropping them.
    pub preserve_quality: bool,
    /// Lossy quality transform applied before storage.
    pub quality_mode: QualityMode,
    /// Zstandard compression level for section payloads (up to 22).
    pub zstd_level: i32,
    /// Read regions extracted as embedded index reads before mapping.
    pub embedded_index_ranges: Vec<EmbeddedIndexRange>,
    /// Separate index-read FASTQ files (I1/I2) stored alongside the mates.
    pub index_read_files: Vec<std::path::PathBuf>,
    /// GTF annotation for star splice-junction index construction.
    pub star_gtf: Option<std::path::PathBuf>,
    /// star `sjdbOverhang` for splice-junction database construction.
    pub star_sjdb_overhang: u32,
    /// star `--sa-sparseD` suffix-array sparsity factor.
    pub star_sa_sparse_d: u32,
    /// star `--sa-index-nbases` suffix-array seed length.
    pub star_sa_index_nbases: u32,
    /// star `--chr-bin-nbits` chromosome address width as a base-2 exponent.
    pub star_chr_bin_nbits: u32,
}

impl Default for AlignmentCompressionOptions {
    fn default() -> Self {
        Self {
            aligner: AlignerKind::Minibwa,
            chunk_records: 250_000,
            threads: default_thread_count(),
            temp_dir: None,
            preserve_ids: true,
            id_prefix_chars: None,
            id_field: None,
            preserve_quality: true,
            quality_mode: QualityMode::Lossless,
            zstd_level: 19,
            embedded_index_ranges: Vec::new(),
            index_read_files: Vec::new(),
            star_gtf: None,
            star_sjdb_overhang: 100,
            star_sa_sparse_d: 1,
            star_sa_index_nbases: 14,
            star_chr_bin_nbits: 18,
        }
    }
}

/// Options for decoding an archive.
#[derive(Clone, Debug)]
pub struct DecompressionOptions {
    /// Inclusive zero-based logical-record range to decode; `None` means all.
    pub original_range: Option<RangeInclusive<u64>>,
    /// Order in which decoded records are delivered.
    pub output_order: OutputOrder,
    /// Scratch directory for temporary files; defaults to the system temp dir.
    pub temp_dir: Option<std::path::PathBuf>,
    /// Number of chunks decoded concurrently during streaming iteration.
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

/// Record delivery order for decompression.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum OutputOrder {
    /// Restore the original input order of the compressed records.
    #[default]
    Original,
    /// Deliver records in the (reordered) compression order, which is cheaper.
    Compression,
}

/// Receiver for decoded paired records, avoiding whole-archive buffering.
pub trait PairedRecordSink {
    /// Consume one pair; `original_id` is the zero-based input position.
    fn write_pair(&mut self, original_id: u64, first: ReadRecord, second: ReadRecord)
        -> Result<()>;

    /// Called once after the final pair.
    fn finish(&mut self) -> Result<()>;
}

/// Whether an archive stores single-end or paired-end reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlignmentLayout {
    /// One biological read per logical record.
    SingleEnd,
    /// Two mated reads per logical record.
    PairedEnd,
}

/// Which mate an embedded index range is extracted from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmbeddedIndexMate {
    /// Read 1 (forward mate).
    R1,
    /// Read 2 (reverse mate).
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
    /// Mate the region is cut from.
    pub mate: EmbeddedIndexMate,
    /// First included base offset within the read.
    pub start: usize,
    /// End offset (exclusive) within the read.
    pub end: usize,
}

impl EmbeddedIndexRange {
    /// Number of extracted bases.
    pub fn len(&self) -> usize {
        self.end.saturating_sub(self.start)
    }

    /// True when the range extracts no bases.
    pub fn is_empty(&self) -> bool {
        self.start >= self.end
    }
}

impl AlignmentLayout {
    /// Number of biological read streams the layout produces.
    pub fn streams(self) -> usize {
        match self {
            Self::SingleEnd => 1,
            Self::PairedEnd => 2,
        }
    }
}

/// The biological read payload of one logical record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SequenceRecordGroup {
    /// Single-end read.
    Single(ReadRecord),
    /// Paired-end reads (mate 1, mate 2).
    Paired(ReadRecord, ReadRecord),
}

/// One decoded logical record: biological reads plus any index reads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignedRecordGroup {
    /// The biological read or read pair.
    pub sequence: SequenceRecordGroup,
    /// Embedded or separate index reads (I1/I2, UMIs) for this record.
    pub index_reads: Vec<ReadRecord>,
}

/// Receiver for decoded logical records, avoiding whole-archive buffering.
pub trait AlignedRecordSink {
    /// Consume one record; `original_id` is the zero-based input position.
    fn write_group(&mut self, original_id: u64, group: AlignedRecordGroup) -> Result<()>;

    /// Called once after the final record.
    fn finish(&mut self) -> Result<()>;
}

/// Statistics returned by alignment-aware decompression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlignedDecompressionMetadata {
    /// Single- or paired-end layout of the archive.
    pub layout: AlignmentLayout,
    /// Number of biological read streams.
    pub streams: usize,
    /// Number of index-read streams.
    pub index_streams: usize,
    /// Logical records decoded.
    pub records: u64,
    /// Order the records were delivered in.
    pub output_order: OutputOrder,
}

/// Statistics returned by reference-free paired decompression.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DecompressionMetadata {
    /// Number of read streams (always 2 for paired archives).
    pub streams: usize,
    /// Pairs decoded.
    pub pairs: u64,
    /// Order the records were delivered in.
    pub output_order: OutputOrder,
}

/// Statistics returned by compression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressionMetadata {
    /// Number of biological read streams written to the archive.
    pub streams: usize,
    /// Logical records written to the archive.
    pub records: u64,
    /// Whether identifiers were stored.
    pub preserved_ids: bool,
    /// Whether quality values were stored.
    pub preserved_quality: bool,
    /// Whether original record order is restorable.
    pub preserved_order: bool,
}

/// Size and timing profile of an alignment-based compression run.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AlignedCompressionProfile {
    /// End-to-end wall time.
    pub total: Duration,
    /// Time spent loading or building the alignment index.
    pub index: Duration,
    /// Time spent mapping reads and spooling alignment data.
    pub map_and_spool: Duration,
    /// Time spent compacting and storing the embedded reference.
    pub reference: Duration,
    /// Time spent encoding and compressing chunks.
    pub chunks: Duration,
    /// Logical records compressed.
    pub records: u64,
    /// Chunks written to the archive.
    pub chunk_count: u64,
    /// Bytes written to temporary alignment spool files.
    pub spool_bytes: u64,
    /// Peak number of encoded chunks awaiting their in-order slot.
    pub max_pending_chunks: u64,
    /// Candidate alignments produced by the seed stage.
    pub candidate_hits: u64,
    /// Reads with a selected alignment.
    pub mapped_reads: u64,
    /// Reads stored through the raw fallback.
    pub raw_reads: u64,
    /// Alignments ranked below the selected hit.
    pub selected_non_first_hits: u64,
    /// Alignments discarded during selection.
    pub rejected_alignments: u64,
    /// Reference bases embedded in the archive.
    pub reference_bases: u64,
    /// Uncompressed size of the embedded reference.
    pub reference_raw_bytes: u64,
    /// Compressed size of the embedded reference.
    pub reference_compressed_bytes: u64,
    /// Uncompressed size of the identifier stream.
    pub ids_raw_bytes: u64,
    /// Uncompressed size of the quality stream.
    pub qualities_raw_bytes: u64,
    /// Uncompressed size of the record streams.
    pub records_raw_bytes: u64,
    /// Compressed size of the identifier stream.
    pub ids_compressed_bytes: u64,
    /// Compressed size of the quality stream.
    pub qualities_compressed_bytes: u64,
    /// Compressed size of the record streams.
    pub records_compressed_bytes: u64,
    /// Uncompressed size of the index-read streams.
    pub indexes_raw_bytes: u64,
    /// Compressed size of the index-read streams.
    pub indexes_compressed_bytes: u64,
    /// Compressed size of all chunk payloads.
    pub chunks_compressed_bytes: u64,
    /// Final archive size in bytes.
    pub archive_bytes: u64,
}

/// Decoded records grouped by stream (mate 1 first, then mate 2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArchiveRecords {
    /// One vector of records per read stream.
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

/// Compress single-end plain, gzip, or zstd FASTQ against a genome reference.
pub fn compress_aligned_single_fastq_file(
    genome_fasta: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    compress_aligned_single_fastq_file_profiled(genome_fasta, input, output, options)
        .map(|(metadata, _)| metadata)
}

/// Compress single-end FASTQ and report archive details.
pub fn compress_aligned_single_fastq_file_profiled(
    genome_fasta: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    alignment_archive::compress_single_fastq_file(genome_fasta, input, output, options)
}

/// Compress single-end FASTQ using an existing index.
pub fn compress_aligned_single_fastq_file_with_index(
    index_file: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    compress_aligned_single_fastq_file_with_index_profiled(index_file, input, output, options)
        .map(|(metadata, _)| metadata)
}

/// Compress single-end FASTQ with an existing SeqArc index and report details.
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

/// Decompress an alignment archive into single-end records.
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

/// Eagerly decode aligned records into an [`AlignedRecordSink`].
pub fn decompress_aligned_into<R: Read, S: AlignedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<AlignedDecompressionMetadata> {
    alignment_archive::decompress_into(input, options, sink)
}

/// Eagerly decode paired records into a [`PairedRecordSink`].
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
