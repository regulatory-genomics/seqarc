#[cfg(test)]
use crate::aligners::SelectedAlignment;
use crate::aligners::{
    self, compact_boundary, compact_reference, effective_cigar, reverse_complement, AlignmentIndex,
    AlignmentSelection, KeptRun, MappingProgressStats, ReferenceBases, ReferenceMask,
    MAP_TASK_READS,
};
#[cfg(test)]
use crate::minibwa::{self, Strand};
use crate::minibwa::{CigarKind, CigarOp};
#[cfg(test)]
use crate::packed_sequence::PackedSequence;
use crate::quality::{decode_stream_bytes, encode_stream_bytes};
use crate::{
    apply_quality, identifier_token, open_fastx_reader, select_id_field, truncate_read_id,
    AlignedCompressionProfile, AlignedDecompressionMetadata, AlignedRecordGroup, AlignedRecordSink,
    AlignmentLayout, ArchiveRecords, CompressionMetadata, DecompressionOptions, EmbeddedIndexMate,
    EmbeddedIndexRange, Error, FastxReader, Format, GenomeAlignment, PairedRecordSink, ReadRecord,
    Result, SequenceRecordGroup,
};
use bed_utils::extsort::{ExternalChunk, ExternalChunkBuilder};
use indicatif::{ProgressBar, ProgressStyle};
use rkyv::{Archive, Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Cursor, Read, Seek, SeekFrom, Write};
use std::ops::RangeInclusive;
use std::path::Path;
use std::sync::Arc;

const MAGIC: &[u8; 9] = b"SEQARC-AL";
const VERSION_MAJOR: u16 = 0;
const VERSION_MINOR: u16 = 28;
/// Header bytes covered by the preamble CRC: magic, version, flags, stream
/// count, position widths, block lengths, and chunk count.
const PREAMBLE_BYTES: usize = 42;
const MAP_TASKS_PER_WORKER: usize = 4;
const MAX_MAPPING_BATCH_RECORDS: usize = 131_072;
const PIPELINE_BATCHES: usize = 2;

#[derive(Clone, Debug)]
struct AlignedPair {
    first: ReadRecord,
    second: Option<ReadRecord>,
    first_hit: Option<GenomeAlignment>,
    second_hit: Option<GenomeAlignment>,
    first_selection: AlignmentSelection,
    second_selection: AlignmentSelection,
    indexes: Vec<SpoolIndexRead>,
    index_reads: Vec<ReadRecord>,
}

impl AlignedPair {
    fn layout(&self) -> AlignmentLayout {
        if self.second.is_some() {
            AlignmentLayout::PairedEnd
        } else {
            AlignmentLayout::SingleEnd
        }
    }

    fn record(&self, stream: usize) -> &ReadRecord {
        match stream {
            0 => &self.first,
            1 => self.second.as_ref().expect("paired stream is missing"),
            _ => unreachable!(),
        }
    }

    fn hit(&self, stream: usize) -> Option<&GenomeAlignment> {
        match stream {
            0 => self.first_hit.as_ref(),
            1 => self.second_hit.as_ref(),
            _ => unreachable!(),
        }
    }

    fn index_read(&self, stream: usize) -> Option<&ReadRecord> {
        self.index_reads.get(stream)
    }

    fn index_stream_count(&self) -> usize {
        self.index_reads.len()
    }
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolPair {
    first: SpoolRead,
    second: Option<SpoolRead>,
    indexes: Vec<SpoolIndexRead>,
    index_reads: Vec<SpoolIndexRead>,
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolIndexRead {
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
}

#[derive(Clone, Debug)]
struct IndexDescriptor {
    range: EmbeddedIndexRange,
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolRead {
    id: Vec<u8>,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
    alignment: Option<SpoolAlignment>,
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolAlignment {
    start: i64,
    query_start: i32,
    query_end: i32,
    is_reverse: bool,
    cigar: Vec<SpoolCigarOp>,
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolCigarOp {
    op: u8,
    len: u32,
}

struct EncodedChunk {
    pairs: u32,
    raw_len: u64,
    ids_len: u64,
    qualities_len: u64,
    records_len: u64,
    ids_compressed_len: u64,
    qualities_compressed_len: u64,
    records_compressed_len: u64,
    indexes_raw_len: u64,
    indexes_compressed_len: u64,
    compressed: Vec<u8>,
}

struct EncodedReference {
    block: Vec<u8>,
    raw_len: u64,
}

struct EncodedRecordStream {
    bytes: Vec<u8>,
    indexes_raw_len: u64,
    indexes_compressed_len: u64,
}

struct EncodedArchiveChunk {
    first_record: u64,
    record_count: usize,
    raw_len: usize,
    ids_len: usize,
    qualities_len: usize,
    records_len: usize,
    compressed: Vec<u8>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ArchiveOp {
    Match = 0,
    Literal = 1,
    RefAdvance = 2,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ArchiveRun {
    op: ArchiveOp,
    len: usize,
}

struct DecodedArchiveChunk {
    first_record: u64,
    next_record: usize,
    ids: Vec<Vec<Vec<u8>>>,
    qualities: Vec<Vec<Option<Vec<u8>>>>,
    sequences: Vec<Vec<Vec<u8>>>,
    index_reads: Vec<DecodedIndexStream>,
}

struct DecodedIndexStream {
    qualities: Vec<Option<Vec<u8>>>,
    sequences: Vec<Vec<u8>>,
}

struct DecodedRecordSections {
    biological: Vec<Vec<Vec<u8>>>,
    index_streams: Vec<DecodedIndexStream>,
}

struct ArchiveHeader {
    flags: u16,
    id_mode: IdMode,
    layout: AlignmentLayout,
    chunk_count: u64,
    reference: Vec<u8>,
    indexes: Vec<IndexDescriptor>,
    absolute_position_width: usize,
    delta_position_width: usize,
    index_streams: usize,
}

fn position_width(reference_len: usize) -> usize {
    let max = reference_len.saturating_sub(1);
    if max <= u8::MAX as usize {
        1
    } else if max <= u16::MAX as usize {
        2
    } else if max <= 0x00ff_ffff {
        3
    } else if max <= u32::MAX as usize {
        4
    } else {
        8
    }
}

fn delta_width(reference_len: usize) -> usize {
    let max_delta = reference_len.saturating_sub(1) as u64;
    let max_encoded = max_delta.saturating_mul(2);
    if max_encoded <= u8::MAX as u64 {
        1
    } else if max_encoded <= u16::MAX as u64 {
        2
    } else if max_encoded <= 0x00ff_ffff {
        3
    } else if max_encoded <= u32::MAX as u64 {
        4
    } else {
        8
    }
}

fn write_position(output: &mut Vec<u8>, value: u64, width: usize) -> Result<()> {
    let max = if width == 8 {
        u64::MAX
    } else {
        (1u64 << (width * 8)) - 1
    };
    if value > max {
        return Err(Error::InvalidInput(
            "position exceeds selected width".into(),
        ));
    }
    for byte in 0..width {
        output.push((value >> (byte * 8)) as u8);
    }
    Ok(())
}

fn read_position<T: AsRef<[u8]>>(input: &mut Cursor<T>, width: usize) -> Result<u64> {
    let mut value = 0u64;
    for byte in 0..width {
        value |= u64::from(read_u8(input)?) << (byte * 8);
    }
    Ok(value)
}

#[derive(Clone, Debug)]
struct V21NoiseEvent {
    delta: u16,
    symbol: u8,
    skip: u64,
}

struct V21AlignedFields {
    start: u64,
    reverse: bool,
    events: Vec<V21NoiseEvent>,
    emits: Vec<Vec<u8>>,
}

const V21_MATE_SECTION_COUNT: usize = 8;

/// Streaming decoder that yields aligned records lazily in original order.
///
/// Chunks decode concurrently on a bounded worker pool; exhaust the iterator
/// to validate every chunk and detect trailing bytes.
pub struct AlignedRecordIter<R: Read> {
    input: R,
    pool: Option<rayon::ThreadPool>,
    result_tx: Option<crossbeam_channel::Sender<(u64, Result<DecodedArchiveChunk>)>>,
    result_rx: Option<crossbeam_channel::Receiver<(u64, Result<DecodedArchiveChunk>)>>,
    reference: Arc<Vec<u8>>,
    id_mode: IdMode,
    layout: AlignmentLayout,
    preserve_quality: bool,
    indexes: Arc<Vec<IndexDescriptor>>,
    absolute_position_width: usize,
    delta_position_width: usize,
    index_streams: usize,
    original_range: Option<RangeInclusive<u64>>,
    best_effort: bool,
    total_chunks: u64,
    decoded_chunks: u64,
    chunks_remaining: u64,
    next_submit_chunk: u64,
    next_submit_record: u64,
    next_yield_chunk: u64,
    terminal_chunk: u64,
    outstanding_jobs: usize,
    target_jobs: usize,
    max_resident_chunks: usize,
    resident_bytes: u64,
    max_resident_bytes: u64,
    resident_estimates: BTreeMap<u64, u64>,
    pending: BTreeMap<u64, Result<DecodedArchiveChunk>>,
    skipped_chunks: BTreeMap<u64, u64>,
    skipped_chunk_total: u64,
    skipped_record_total: u64,
    current: Option<DecodedArchiveChunk>,
    input_finished: bool,
    done: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum IdMode {
    None = 0,
    Independent = 1,
    IdenticalMates = 2,
    SlashMate = 3,
}

impl IdMode {
    fn from_flags(flags: u16) -> Result<Self> {
        match (flags >> 3) & 3 {
            0 => Ok(Self::None),
            1 => Ok(Self::Independent),
            2 => Ok(Self::IdenticalMates),
            3 => Ok(Self::SlashMate),
            _ => unreachable!(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct IdModeDetector {
    mode: Option<IdMode>,
    preserve: bool,
    layout: AlignmentLayout,
}

impl IdModeDetector {
    fn new(preserve: bool, layout: AlignmentLayout) -> Self {
        Self {
            mode: None,
            preserve,
            layout,
        }
    }

    fn observe(&mut self, first: &[u8], second: Option<&[u8]>) {
        if !self.preserve {
            return;
        }
        if self.layout == AlignmentLayout::SingleEnd {
            self.mode = Some(IdMode::Independent);
            return;
        }
        let second = second.expect("paired identifier is missing");
        let observed = if first == second {
            IdMode::IdenticalMates
        } else if first.ends_with(b"/1")
            && second.ends_with(b"/2")
            && first[..first.len() - 2] == second[..second.len() - 2]
        {
            IdMode::SlashMate
        } else {
            IdMode::Independent
        };
        self.mode = Some(match (self.mode, observed) {
            (None, observed) => observed,
            (Some(IdMode::Independent), _) | (_, IdMode::Independent) => IdMode::Independent,
            (Some(IdMode::IdenticalMates), IdMode::IdenticalMates) => IdMode::IdenticalMates,
            (Some(IdMode::SlashMate), IdMode::SlashMate) => IdMode::SlashMate,
            _ => IdMode::Independent,
        });
    }

    fn finish(self) -> IdMode {
        if self.preserve {
            self.mode.unwrap_or(match self.layout {
                AlignmentLayout::SingleEnd => IdMode::Independent,
                AlignmentLayout::PairedEnd => IdMode::IdenticalMates,
            })
        } else {
            IdMode::None
        }
    }
}

fn index_build_options(
    options: &crate::AlignmentCompressionOptions,
    workspace: &Path,
) -> crate::GenomeIndexBuildOptions {
    crate::GenomeIndexBuildOptions {
        temp_dir: Some(workspace.to_path_buf()),
        aligner: options.aligner,
        threads: options.threads,
        star_gtf: options.star_gtf.clone(),
        star_sjdb_overhang: options.star_sjdb_overhang,
        star_sa_sparse_d: options.star_sa_sparse_d,
        star_sa_index_nbases: options.star_sa_index_nbases,
        star_chr_bin_nbits: options.star_chr_bin_nbits,
    }
}

pub(crate) fn compress_paired_fastq_files(
    genome_fasta: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let total_started = std::time::Instant::now();
    let workspace =
        crate::temporary_directory(options.temp_dir.as_deref(), "seqarc-aligned-index-")?;
    eprintln!(
        "[seqarc] stage=index_build start reference={}",
        genome_fasta.as_ref().display()
    );
    let index_started = std::time::Instant::now();
    let index_path = workspace.path().join("index.seqarc-index");
    crate::build_alignment_index(
        genome_fasta,
        &index_path,
        &index_build_options(options, workspace.path()),
    )?;
    let index = AlignmentIndex::load(&index_path)?;
    let index_elapsed = index_started.elapsed();
    eprintln!(
        "[seqarc] stage=index_build done elapsed={} reference_bases={}",
        display_duration(index_elapsed),
        index.reference_len()
    );
    let (metadata, mut profile) = compress_with_index(
        &index,
        input_1,
        Some(input_2),
        output,
        options,
        AlignmentLayout::PairedEnd,
    )?;
    profile.index = index_elapsed;
    profile.total = total_started.elapsed();
    eprintln!(
        "[seqarc] stage=total done elapsed={}",
        display_duration(profile.total)
    );
    Ok((metadata, profile))
}

pub(crate) fn compress_paired_fastq_files_with_index(
    index_file: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let total_started = std::time::Instant::now();
    let index_file = index_file.as_ref();
    eprintln!(
        "[seqarc] stage=index_load start path={}",
        index_file.display()
    );
    let index_started = std::time::Instant::now();
    let index = AlignmentIndex::load(index_file)?;
    let index_elapsed = index_started.elapsed();
    eprintln!(
        "[seqarc] stage=index_load done elapsed={} reference_bases={}",
        display_duration(index_elapsed),
        index.reference_len()
    );
    let (metadata, mut profile) = compress_with_index(
        &index,
        input_1,
        Some(input_2),
        output,
        options,
        AlignmentLayout::PairedEnd,
    )?;
    profile.index = index_elapsed;
    profile.total = total_started.elapsed();
    eprintln!(
        "[seqarc] stage=total done elapsed={}",
        display_duration(profile.total)
    );
    Ok((metadata, profile))
}

pub(crate) fn compress_single_fastq_file(
    genome_fasta: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let total_started = std::time::Instant::now();
    let workspace =
        crate::temporary_directory(options.temp_dir.as_deref(), "seqarc-aligned-index-")?;
    eprintln!(
        "[seqarc] stage=index_build start reference={}",
        genome_fasta.as_ref().display()
    );
    let index_started = std::time::Instant::now();
    let index_path = workspace.path().join("index.seqarc-index");
    crate::build_alignment_index(
        genome_fasta,
        &index_path,
        &index_build_options(options, workspace.path()),
    )?;
    let index = AlignmentIndex::load(&index_path)?;
    let index_elapsed = index_started.elapsed();
    eprintln!(
        "[seqarc] stage=index_build done elapsed={} reference_bases={}",
        display_duration(index_elapsed),
        index.reference_len()
    );
    let (metadata, mut profile) = compress_with_index(
        &index,
        input,
        Option::<&Path>::None,
        output,
        options,
        AlignmentLayout::SingleEnd,
    )?;
    profile.index = index_elapsed;
    profile.total = total_started.elapsed();
    eprintln!(
        "[seqarc] stage=total done elapsed={}",
        display_duration(profile.total)
    );
    Ok((metadata, profile))
}

pub(crate) fn compress_single_fastq_file_with_index(
    index_file: impl AsRef<Path>,
    input: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let total_started = std::time::Instant::now();
    let index_file = index_file.as_ref();
    eprintln!(
        "[seqarc] stage=index_load start path={}",
        index_file.display()
    );
    let index_started = std::time::Instant::now();
    let index = AlignmentIndex::load(index_file)?;
    let index_elapsed = index_started.elapsed();
    eprintln!(
        "[seqarc] stage=index_load done elapsed={} reference_bases={}",
        display_duration(index_elapsed),
        index.reference_len()
    );
    let (metadata, mut profile) = compress_with_index(
        &index,
        input,
        Option::<&Path>::None,
        output,
        options,
        AlignmentLayout::SingleEnd,
    )?;
    profile.index = index_elapsed;
    profile.total = total_started.elapsed();
    eprintln!(
        "[seqarc] stage=total done elapsed={}",
        display_duration(profile.total)
    );
    Ok((metadata, profile))
}

fn apply_identifier_policy(
    record: &mut ReadRecord,
    preserve_ids: bool,
    prefix_chars: Option<usize>,
    id_field: Option<usize>,
) -> Result<()> {
    if !preserve_ids {
        record.id.clear();
    } else {
        if let Some(id_field) = id_field {
            select_id_field(record, id_field)?;
        }
        if let Some(prefix_chars) = prefix_chars {
            truncate_read_id(record, prefix_chars)?;
        }
    }
    Ok(())
}

fn validate_index_ranges(ranges: &[EmbeddedIndexRange], layout: AlignmentLayout) -> Result<()> {
    if layout == AlignmentLayout::SingleEnd
        && ranges
            .iter()
            .any(|range| range.mate == EmbeddedIndexMate::R2)
    {
        return Err(Error::InvalidInput(
            "mate R2 index range requires paired-end input".into(),
        ));
    }
    for range in ranges {
        if range.is_empty() {
            return Err(Error::InvalidInput(
                "embedded index ranges must be nonempty".into(),
            ));
        }
    }
    for mate in [EmbeddedIndexMate::R1, EmbeddedIndexMate::R2] {
        let mut selected = ranges
            .iter()
            .filter(|range| range.mate == mate)
            .collect::<Vec<_>>();
        selected.sort_by_key(|range| range.start);
        if selected
            .windows(2)
            .any(|window| window[0].end > window[1].start)
        {
            return Err(Error::InvalidInput(
                "embedded index ranges overlap within a mate".into(),
            ));
        }
    }
    Ok(())
}

fn split_record_indexes(
    record: &mut ReadRecord,
    ranges: &[(usize, &EmbeddedIndexRange)],
    preserve_quality: bool,
) -> Result<Vec<(usize, SpoolIndexRead)>> {
    if ranges.is_empty() {
        return Ok(Vec::new());
    }
    let removed = ranges.iter().try_fold(0usize, |total, (_, range)| {
        if range.end > record.sequence.len() {
            return Err(Error::InvalidInput(
                "embedded index range exceeds read length".into(),
            ));
        }
        total
            .checked_add(range.len())
            .ok_or_else(|| Error::InvalidInput("embedded index length overflow".into()))
    })?;
    if removed > record.sequence.len() {
        return Err(Error::InvalidInput(
            "embedded index ranges exceed the biological read length".into(),
        ));
    }
    let quality = record.quality.as_ref().ok_or_else(|| {
        Error::InvalidInput("FASTQ record is missing index-range quality values".into())
    })?;
    if quality.len() != record.sequence.len() {
        return Err(Error::InvalidInput(
            "sequence and quality lengths differ before index extraction".into(),
        ));
    }

    let indexes = ranges
        .iter()
        .map(|(declaration, range)| {
            Ok((
                *declaration,
                SpoolIndexRead {
                    sequence: record.sequence[range.start..range.end].to_vec(),
                    quality: preserve_quality.then(|| quality[range.start..range.end].to_vec()),
                },
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut descending = ranges.to_vec();
    descending.sort_by_key(|(_, range)| std::cmp::Reverse(range.start));
    for (_, range) in descending {
        record.sequence.drain(range.start..range.end);
        record
            .quality
            .as_mut()
            .expect("quality was validated")
            .drain(range.start..range.end);
    }
    // The biological read may be empty when the embedded index ranges consume
    // the entire read; only the sequence/quality length agreement is required.
    record.validate(true, true)?;
    Ok(indexes)
}

fn split_pair_indexes(
    first: &mut ReadRecord,
    second: Option<&mut ReadRecord>,
    ranges: &[EmbeddedIndexRange],
    preserve_quality: bool,
) -> Result<Vec<SpoolIndexRead>> {
    let first_ranges = ranges
        .iter()
        .enumerate()
        .filter(|(_, range)| range.mate == EmbeddedIndexMate::R1)
        .collect::<Vec<_>>();
    let second_ranges = ranges
        .iter()
        .enumerate()
        .filter(|(_, range)| range.mate == EmbeddedIndexMate::R2)
        .collect::<Vec<_>>();
    let mut extracted = split_record_indexes(first, &first_ranges, preserve_quality)?;
    if !second_ranges.is_empty() {
        let second = second.ok_or_else(|| {
            Error::InvalidInput("mate R2 index range requires paired-end input".into())
        })?;
        extracted.extend(split_record_indexes(
            second,
            &second_ranges,
            preserve_quality,
        )?);
    }
    extracted.sort_by_key(|(declaration, _)| *declaration);
    if extracted.len() != ranges.len()
        || extracted
            .iter()
            .enumerate()
            .any(|(expected, (actual, _))| expected != *actual)
    {
        return Err(Error::InvalidInput(
            "embedded index extraction order mismatch".into(),
        ));
    }
    Ok(extracted.into_iter().map(|(_, index)| index).collect())
}

/// Reject unsupported alignment-compression option combinations up front.
fn validate_compress_options(
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<()> {
    validate_index_ranges(&options.embedded_index_ranges, layout)?;
    if options.index_read_files.len() > 2 {
        return Err(Error::InvalidInput(
            "at most two index read files are supported".into(),
        ));
    }
    if !options.preserve_ids && options.id_prefix_chars.is_some_and(|length| length != 0) {
        return Err(Error::InvalidInput(
            "identifier prefix cannot be retained when identifiers are disabled".into(),
        ));
    }
    if !options.preserve_ids && options.id_field.is_some() {
        return Err(Error::InvalidInput(
            "identifier field selection requires identifier preservation".into(),
        ));
    }
    Ok(())
}

/// Write the global preamble (header fields, embedded reference, index
/// metadata, and CRC) to `output`.
///
/// Returns the file opened read/write, positioned just past the preamble and
/// ready for chunk payloads, plus the offset of the chunk-count field the
/// caller must patch after writing chunks.
#[allow(clippy::too_many_arguments)]
fn write_archive_preamble(
    output: &Path,
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
    preserve_ids: bool,
    id_mode: IdMode,
    absolute_position_width: usize,
    delta_position_width: usize,
    reference_block: &[u8],
    index_metadata: &[u8],
) -> Result<(std::fs::File, u64)> {
    let mut archive = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(output)?;
    let header_size = 9 + 4 + 2 + 1 + 2 + 8 * 3 + 4;
    archive.write_all(MAGIC)?;
    write_u16(&mut archive, VERSION_MAJOR)?;
    write_u16(&mut archive, VERSION_MINOR)?;
    let mut flags = 0u16;
    if preserve_ids {
        flags |= 1;
    }
    if options.preserve_quality {
        flags |= 2;
    }
    if layout == AlignmentLayout::SingleEnd {
        flags |= 4;
    }
    flags |= (id_mode as u16) << 3;
    write_u16(&mut archive, flags)?;
    archive.write_all(&[options.index_read_files.len() as u8])?;
    archive.write_all(&[absolute_position_width as u8, delta_position_width as u8])?;
    write_u64(&mut archive, reference_block.len() as u64)?;
    write_u64(&mut archive, index_metadata.len() as u64)?;
    let chunk_count_pos = archive.stream_position()?;
    write_u64(&mut archive, 0)?;
    write_u32(&mut archive, 0)?;
    debug_assert_eq!(archive.stream_position()?, header_size as u64);
    archive.write_all(reference_block)?;
    archive.write_all(index_metadata)?;
    write_preamble_crc(&mut archive)?;
    Ok((archive, chunk_count_pos))
}

fn compress_with_index<I2: AsRef<Path>>(
    index: &AlignmentIndex,
    input_1: impl AsRef<Path>,
    input_2: Option<I2>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let pipeline_started = std::time::Instant::now();
    validate_compress_options(options, layout)?;
    let preserve_ids = options.preserve_ids && options.id_prefix_chars != Some(0);
    let mapping_threads = options.threads.max(1);
    let mapping_batch_records = mapping_threads
        .saturating_mul(MAP_TASKS_PER_WORKER)
        .saturating_mul(MAP_TASK_READS)
        .clamp(MAP_TASK_READS, MAX_MAPPING_BATCH_RECORDS);
    let mapping_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(mapping_threads)
        .build()
        .map_err(|error| Error::InvalidInput(format!("failed to create mapping pool: {error}")))?;
    let progress = mapping_progress();
    stage_log(&progress, format!(
        "[seqarc] stage=initialize layout={:?} streams={} index_streams={} threads={} mapping_batch_records={} map_task_reads={} chunk_records={} zstd_level={} preserve_ids={} id_prefix_chars={:?} id_field={:?} preserve_quality={}",
        layout,
        layout.streams(),
        options.index_read_files.len(),
        mapping_threads,
        mapping_batch_records,
        MAP_TASK_READS,
        options.chunk_records.max(1),
        options.zstd_level,
        preserve_ids,
        options.id_prefix_chars,
        options.id_field,
        options.preserve_quality
    ));
    let temp_dir =
        crate::temporary_directory(options.temp_dir.as_deref(), "seqarc-aligned-spool-")?;
    let spool_path = temp_dir.path().join("aligned-pairs.lz4");
    let mut spool_builder = ExternalChunkBuilder::new(File::create(&spool_path)?, 1)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let mut mask = ReferenceMask::new(index.reference_len());
    let mut id_detector = IdModeDetector::new(preserve_ids, layout);
    let mut profile = AlignedCompressionProfile::default();
    let input_1 = input_1.as_ref().to_path_buf();
    let input_2 = input_2.map(|path| path.as_ref().to_path_buf());
    let index_files = options.index_read_files.clone();
    let mapper_options = options.clone();
    let (parse_tx, parse_rx) =
        crossbeam_channel::bounded::<Result<Vec<AlignedPair>>>(PIPELINE_BATCHES);
    let (mapped_tx, mapped_rx) =
        crossbeam_channel::bounded::<Result<Vec<AlignedPair>>>(PIPELINE_BATCHES);
    stage_log(&progress, "[seqarc] stage=map_and_spool start");
    let map_started = std::time::Instant::now();
    let result = std::thread::scope(|scope| -> Result<()> {
        let producer_tx = parse_tx.clone();
        let producer_progress = progress.clone();
        let producer = scope.spawn(move || -> Result<()> {
            let mut readers: Vec<FastxReader<Box<dyn std::io::BufRead + Send>>> = Vec::new();
            readers.push(FastxReader::new(open_fastx_reader(input_1)?, Format::Fastq));
            if let Some(input_2) = input_2 {
                readers.push(FastxReader::new(open_fastx_reader(input_2)?, Format::Fastq));
            }
            for path in &index_files {
                readers.push(FastxReader::new(open_fastx_reader(path)?, Format::Fastq));
            }
            let biological_streams = layout.streams();
            let mut pairs = Vec::with_capacity(mapping_batch_records);
            let mut preview_identifier =
                mapper_options.id_prefix_chars.is_some() || mapper_options.id_field.is_some();
            let mut prepare_identifier = |record: &mut ReadRecord| -> Result<()> {
                let original = preview_identifier.then(|| record.id.clone());
                apply_identifier_policy(
                    record,
                    preserve_ids,
                    mapper_options.id_prefix_chars,
                    mapper_options.id_field,
                )?;
                if let Some(original) = original {
                    let old_header = format!("@{}", String::from_utf8_lossy(&original));
                    let new_header = format!("@{}", String::from_utf8_lossy(&record.id));
                    stage_log(
                        &producer_progress,
                        format!("[seqarc] id-transform preview: {old_header:?} -> {new_header:?}"),
                    );
                    preview_identifier = false;
                }
                Ok(())
            };
            loop {
                let next: Vec<Option<Result<ReadRecord>>> =
                    readers.iter_mut().map(|reader| reader.next()).collect();
                if next.iter().all(Option::is_none) {
                    break;
                }
                let mut records = Vec::with_capacity(next.len());
                for item in next {
                    let Some(record) = item else {
                        return Err(Error::InvalidInput(
                            "input streams have different record counts".into(),
                        ));
                    };
                    records.push(record?);
                }
                for record in records.iter().skip(biological_streams) {
                    if identifier_token(&record.id)? != identifier_token(&records[0].id)? {
                        return Err(Error::InvalidInput(
                            "index read identifier does not match mate 1 identifier".into(),
                        ));
                    }
                }
                for record in records.iter_mut() {
                    record.validate(true, false)?;
                    prepare_identifier(record)?;
                    if mapper_options.preserve_quality {
                        apply_quality(record, &mapper_options.quality_mode);
                    }
                }
                for record in records.iter().skip(biological_streams) {
                    if record.id != records[0].id {
                        return Err(Error::InvalidInput(
                            "transformed index read identifier does not match transformed mate 1 identifier".into(),
                        ));
                    }
                }
                let mut record_iter = records.into_iter();
                let mut first = record_iter.next().expect("mate 1 stream is required");
                let mut second = (biological_streams == 2)
                    .then(|| record_iter.next().expect("mate 2 stream is required"));
                let indexes = split_pair_indexes(
                    &mut first,
                    second.as_mut(),
                    &mapper_options.embedded_index_ranges,
                    mapper_options.preserve_quality,
                )?;
                let index_reads = record_iter.collect::<Vec<_>>();
                pairs.push(AlignedPair {
                    first,
                    second,
                    first_hit: None,
                    second_hit: None,
                    first_selection: AlignmentSelection::default(),
                    second_selection: AlignmentSelection::default(),
                    indexes,
                    index_reads,
                });
                if pairs.len() == mapping_batch_records {
                    producer_tx
                        .send(Ok(std::mem::take(&mut pairs)))
                        .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
                    pairs = Vec::with_capacity(mapping_batch_records);
                }
            }
            if !pairs.is_empty() {
                producer_tx
                    .send(Ok(pairs))
                    .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
            }
            Ok(())
        });
        drop(parse_tx);

        let mapper_tx = mapped_tx.clone();
        let mapper_progress = progress.clone();
        let mapper = scope.spawn(move || -> Result<()> {
            let mut progress_stats = MappingProgressStats::new(layout);
            while let Ok(batch) = parse_rx.recv() {
                let mapped = batch.and_then(|mut pairs| {
                    map_pairs(
                        index,
                        &mut pairs,
                        &mapping_pool,
                        &mapper_progress,
                        &mut progress_stats,
                    )?;
                    Ok(pairs)
                });
                if mapper_tx.send(mapped).is_err() {
                    break;
                }
            }
            Ok(())
        });
        drop(mapped_tx);

        while let Ok(batch) = mapped_rx.recv() {
            let pairs = batch?;
            update_mapping_profile(&mut profile, &pairs);
            spool_mapped_pairs(&pairs, &mut id_detector, &mut mask, &mut spool_builder)?;
        }
        producer
            .join()
            .map_err(|_| Error::InvalidInput("FASTQ producer panicked".into()))??;
        mapper
            .join()
            .map_err(|_| Error::InvalidInput("mapping worker panicked".into()))??;
        Ok(())
    });
    result?;

    let spool_chunk = spool_builder
        .finish()
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    drop(spool_chunk);
    let index_descriptors = options
        .embedded_index_ranges
        .iter()
        .cloned()
        .map(|range| IndexDescriptor { range })
        .collect::<Vec<_>>();
    profile.map_and_spool = map_started.elapsed();
    profile.spool_bytes = std::fs::metadata(&spool_path)?.len();
    let mapped_records = (profile.mapped_reads + profile.raw_reads) / layout.streams() as u64;
    stage_log(&progress, format!(
        "[seqarc] stage=map_and_spool done elapsed={} layout={:?} records={} mapped_reads={} raw_reads={} candidate_hits={} spool_bytes={}",
        display_duration(profile.map_and_spool),
        layout,
        mapped_records,
        profile.mapped_reads,
        profile.raw_reads,
        profile.candidate_hits,
        profile.spool_bytes
    ));
    progress.finish_and_clear();
    let progress = stage_progress("pruning reference");
    stage_log(&progress, "[seqarc] stage=reference start");
    let reference_started = std::time::Instant::now();
    let (compact_reference, kept_runs) = compact_reference(index, &mask)?;
    let id_mode = id_detector.finish();
    let absolute_position_width = position_width(compact_reference.len());
    let delta_position_width = delta_width(compact_reference.len());
    let reference_block = encode_reference(&compact_reference, options)?;
    profile.reference_bases = compact_reference.len() as u64;
    profile.reference_raw_bytes = reference_block.raw_len;
    profile.reference_compressed_bytes = reference_block.block.len() as u64;
    profile.reference = reference_started.elapsed();
    stage_log(&progress, format!(
        "[seqarc] stage=reference done elapsed={} original_bases={} retained_bases={} retained_runs={} raw_bytes={} compressed_bytes={}",
        display_duration(profile.reference),
        index.reference_len(),
        profile.reference_bases,
        kept_runs.len(),
        profile.reference_raw_bytes,
        profile.reference_compressed_bytes
    ));
    progress.finish_and_clear();
    let index_metadata = encode_index_metadata(&index_descriptors)?;
    // Read access is required to recompute and patch the preamble CRC.
    let (mut archive, chunk_count_pos) = write_archive_preamble(
        output.as_ref(),
        options,
        layout,
        preserve_ids,
        id_mode,
        absolute_position_width,
        delta_position_width,
        &reference_block.block,
        &index_metadata,
    )?;

    let mut written_chunks = 0u64;
    let mut replay_pair_count = 0u64;
    let chunk_records = options.chunk_records.max(1);
    let total_chunks = mapped_records.div_ceil(chunk_records as u64);
    let progress = chunk_progress(total_chunks);
    let worker_count = options.threads.max(1);
    stage_log(
        &progress,
        format!(
            "[seqarc] stage=chunks start workers={} chunk_records={} total_chunks={}",
            worker_count, chunk_records, total_chunks
        ),
    );
    let chunks_started = std::time::Instant::now();
    let (job_tx, job_rx) =
        crossbeam_channel::bounded::<Result<(usize, Vec<AlignedPair>)>>(worker_count);
    let (result_tx, result_rx) =
        crossbeam_channel::bounded::<Result<(usize, EncodedChunk)>>(worker_count);
    let reference_for_workers = compact_reference.as_slice();
    let index_descriptors_for_workers = index_descriptors.as_slice();
    std::thread::scope(|scope| -> Result<()> {
        let dispatcher = scope.spawn({
            let job_tx = job_tx.clone();
            let spool_path = spool_path.clone();
            let kept_runs = kept_runs.clone();
            move || -> Result<()> {
                let spool = ExternalChunk::<SpoolPair>::open(File::open(spool_path)?)
                    .map_err(|error| Error::InvalidFormat(error.to_string()))?;
                let mut pairs = Vec::with_capacity(chunk_records);
                let mut chunk_index = 0usize;
                for item in spool {
                    let item = item.map_err(|error| Error::InvalidFormat(error.to_string()))?;
                    pairs.push(item.into_aligned_pair(&kept_runs, index.reference_len())?);
                    if pairs.len() == chunk_records {
                        job_tx
                            .send(Ok((chunk_index, std::mem::take(&mut pairs))))
                            .map_err(|_| Error::InvalidInput("chunk workers exited".into()))?;
                        chunk_index += 1;
                        pairs = Vec::with_capacity(chunk_records);
                    }
                }
                if !pairs.is_empty() {
                    job_tx
                        .send(Ok((chunk_index, pairs)))
                        .map_err(|_| Error::InvalidInput("chunk workers exited".into()))?;
                }
                Ok(())
            }
        });
        drop(job_tx);
        let mut workers = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let job_rx = job_rx.clone();
            let result_tx = result_tx.clone();
            workers.push(scope.spawn(move || -> Result<()> {
                while let Ok(job) = job_rx.recv() {
                    let result = job.and_then(|(index, pairs)| {
                        encode_chunk(
                            reference_for_workers,
                            &pairs,
                            absolute_position_width,
                            delta_position_width,
                            options,
                            id_mode,
                            layout,
                            index_descriptors_for_workers,
                        )
                        .map(|chunk| (index, chunk))
                    });
                    if result_tx.send(result).is_err() {
                        break;
                    }
                }
                Ok(())
            }));
        }
        drop(job_rx);
        drop(result_tx);
        let mut pending = BTreeMap::new();
        let mut next_chunk = 0usize;
        while let Ok(result) = result_rx.recv() {
            let (index, chunk) = result?;
            progress.inc(1);
            pending.insert(index, chunk);
            profile.max_pending_chunks = profile.max_pending_chunks.max(pending.len() as u64);
            while let Some(chunk) = pending.remove(&next_chunk) {
                profile.ids_raw_bytes += chunk.ids_len;
                profile.qualities_raw_bytes += chunk.qualities_len;
                profile.records_raw_bytes += chunk.records_len;
                profile.ids_compressed_bytes += chunk.ids_compressed_len;
                profile.qualities_compressed_bytes += chunk.qualities_compressed_len;
                profile.records_compressed_bytes += chunk.records_compressed_len;
                profile.indexes_raw_bytes += chunk.indexes_raw_len;
                profile.indexes_compressed_bytes += chunk.indexes_compressed_len;
                profile.chunks_compressed_bytes += chunk.compressed.len() as u64;
                replay_pair_count += chunk.pairs as u64;
                write_encoded_chunk(&mut archive, chunk)?;
                written_chunks += 1;
                next_chunk += 1;
            }
        }
        dispatcher
            .join()
            .map_err(|_| Error::InvalidInput("chunk dispatcher panicked".into()))??;
        for worker in workers {
            worker
                .join()
                .map_err(|_| Error::InvalidInput("chunk worker panicked".into()))??;
        }
        if !pending.is_empty() {
            return Err(Error::InvalidFormat("missing chunk result".into()));
        }
        Ok(())
    })?;
    profile.chunks = chunks_started.elapsed();
    profile.chunk_count = written_chunks;
    stage_log(&progress, format!(
        "[seqarc] stage=chunks done elapsed={} chunks={} records={} raw_bytes={} compressed_bytes={} max_pending={}",
        display_duration(profile.chunks),
        profile.chunk_count,
        replay_pair_count,
        profile.ids_raw_bytes
            + profile.qualities_raw_bytes
            + profile.records_raw_bytes
            + profile.indexes_raw_bytes,
        profile.chunks_compressed_bytes,
        profile.max_pending_chunks
    ));

    let end = archive.stream_position()?;
    archive.seek(SeekFrom::Start(chunk_count_pos))?;
    write_u64(&mut archive, written_chunks)?;
    write_preamble_crc(&mut archive)?;
    archive.seek(SeekFrom::Start(end))?;
    archive.flush()?;
    profile.archive_bytes = end;
    profile.records = replay_pair_count;
    stage_log(
        &progress,
        format!(
            "[seqarc] stage=archive done elapsed={} layout={:?} records={} chunks={} archive_bytes={}",
            display_duration(pipeline_started.elapsed()),
            layout,
            profile.records,
            profile.chunk_count,
            profile.archive_bytes
        ),
    );
    progress.finish_and_clear();

    Ok((
        CompressionMetadata {
            streams: layout.streams(),
            records: replay_pair_count,
            preserved_ids: preserve_ids,
            preserved_quality: options.preserve_quality,
            preserved_order: true,
        },
        profile,
    ))
}

pub(crate) fn decompress_paired<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let records = decompress_iter(input, options)?;
    if records.layout() != AlignmentLayout::PairedEnd {
        return Err(Error::InvalidFormat(
            "single-end archive passed to paired decoder".into(),
        ));
    }
    let mut streams = vec![Vec::new(), Vec::new()];
    for _ in 0..records.index_streams() {
        streams.push(Vec::new());
    }
    let records: Box<dyn Iterator<Item = Result<(u64, AlignedRecordGroup)>>> = match options
        .max_records
        .and_then(|limit| usize::try_from(limit).ok())
    {
        Some(limit) => Box::new(records.take(limit)),
        None => Box::new(records),
    };
    for record in records {
        let (_, group) = record?;
        let AlignedRecordGroup {
            sequence,
            index_reads,
        } = group;
        match sequence {
            SequenceRecordGroup::Paired(first, second) => {
                streams[0].push(first);
                streams[1].push(second);
            }
            SequenceRecordGroup::Single(_) => unreachable!(),
        }
        for (stream, index_read) in streams.iter_mut().skip(2).zip(index_reads) {
            stream.push(index_read);
        }
    }
    Ok(ArchiveRecords { streams })
}

pub(crate) fn decompress_single<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let records = decompress_iter(input, options)?;
    if records.layout() != AlignmentLayout::SingleEnd {
        return Err(Error::InvalidFormat(
            "paired archive passed to single-end decoder".into(),
        ));
    }
    let mut stream = Vec::new();
    let mut index_streams = vec![Vec::new(); records.index_streams()];
    let records: Box<dyn Iterator<Item = Result<(u64, AlignedRecordGroup)>>> = match options
        .max_records
        .and_then(|limit| usize::try_from(limit).ok())
    {
        Some(limit) => Box::new(records.take(limit)),
        None => Box::new(records),
    };
    for record in records {
        let (_, group) = record?;
        let AlignedRecordGroup {
            sequence,
            index_reads,
        } = group;
        match sequence {
            SequenceRecordGroup::Single(record) => stream.push(record),
            SequenceRecordGroup::Paired(_, _) => unreachable!(),
        }
        for (stream, index_read) in index_streams.iter_mut().zip(index_reads) {
            stream.push(index_read);
        }
    }
    let mut streams = vec![stream];
    streams.extend(index_streams);
    Ok(ArchiveRecords { streams })
}

pub(crate) fn decompress_paired_into<R: Read, S: PairedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<crate::DecompressionMetadata> {
    let records = decompress_iter(input, options)?;
    if records.layout() != AlignmentLayout::PairedEnd {
        return Err(Error::InvalidFormat(
            "single-end archive passed to paired decoder".into(),
        ));
    }
    let mut pairs = 0;
    for record in records {
        if options.max_records.is_some_and(|limit| pairs >= limit) {
            break;
        }
        let (original_id, group) = record?;
        match group.sequence {
            SequenceRecordGroup::Paired(first, second) => {
                sink.write_pair(original_id, first, second)?;
                pairs += 1;
            }
            SequenceRecordGroup::Single(_) => unreachable!(),
        }
    }
    sink.finish()?;
    Ok(crate::DecompressionMetadata {
        streams: 2,
        pairs,
        output_order: crate::OutputOrder::Original,
    })
}

pub(crate) fn decompress_into<R: Read, S: AlignedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<AlignedDecompressionMetadata> {
    let decoded = decompress_iter(input, options)?;
    let layout = decoded.layout();
    let index_streams = decoded.index_streams();
    let mut records = 0;
    for record in decoded {
        if options.max_records.is_some_and(|limit| records >= limit) {
            break;
        }
        let (original_id, group) = record?;
        sink.write_group(original_id, group)?;
        records += 1;
    }
    sink.finish()?;
    Ok(AlignedDecompressionMetadata {
        layout,
        streams: layout.streams(),
        index_streams,
        records,
        output_order: crate::OutputOrder::Original,
    })
}

pub(crate) fn decompress_iter<R: Read>(
    mut input: R,
    options: &DecompressionOptions,
) -> Result<AlignedRecordIter<R>> {
    let header = read_header(&mut input)?;
    let decoder_threads = options.threads.max(1);
    let (pool, result_tx, result_rx) = if decoder_threads == 1 {
        (None, None, None)
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(decoder_threads)
            .build()
            .map_err(|error| {
                Error::InvalidInput(format!("failed to create decoder pool: {error}"))
            })?;
        let (result_tx, result_rx) = crossbeam_channel::unbounded();
        (Some(pool), Some(result_tx), Some(result_rx))
    };
    Ok(AlignedRecordIter {
        input,
        pool,
        result_tx,
        result_rx,
        reference: Arc::new(header.reference),
        id_mode: header.id_mode,
        layout: header.layout,
        preserve_quality: header.flags & 2 != 0,
        indexes: Arc::new(header.indexes),
        absolute_position_width: header.absolute_position_width,
        delta_position_width: header.delta_position_width,
        index_streams: header.index_streams,
        original_range: options.original_range.clone(),
        best_effort: options.best_effort,
        total_chunks: header.chunk_count,
        decoded_chunks: 0,
        chunks_remaining: header.chunk_count,
        next_submit_chunk: 0,
        next_submit_record: 0,
        next_yield_chunk: 0,
        terminal_chunk: header.chunk_count,
        outstanding_jobs: 0,
        target_jobs: decoder_threads,
        max_resident_chunks: if decoder_threads == 1 {
            1
        } else {
            decoder_threads.saturating_mul(2)
        },
        resident_bytes: 0,
        // Resident decoded chunks are bounded by memory rather than a chunk
        // count so that adding decoder threads does not inflate peak RSS.
        max_resident_bytes: 3 * 1024 * 1024 * 1024,
        resident_estimates: BTreeMap::new(),
        pending: BTreeMap::new(),
        skipped_chunks: BTreeMap::new(),
        skipped_chunk_total: 0,
        skipped_record_total: 0,
        current: None,
        input_finished: false,
        done: false,
    })
}

impl<R: Read> AlignedRecordIter<R> {
    /// Layout (single- or paired-end) declared by the archive.
    pub fn layout(&self) -> AlignmentLayout {
        self.layout
    }

    /// Number of biological read streams.
    pub fn streams(&self) -> usize {
        self.layout.streams()
    }

    /// Number of index-read streams.
    pub fn index_streams(&self) -> usize {
        self.index_streams
    }

    /// Total chunks declared by the archive header.
    pub fn total_chunks(&self) -> u64 {
        self.total_chunks
    }

    /// Number of chunks skipped in best-effort mode after failed checksums.
    pub fn skipped_chunks(&self) -> u64 {
        self.skipped_chunk_total
    }

    /// Estimated number of records lost in chunks skipped by best-effort mode.
    pub fn skipped_records(&self) -> u64 {
        self.skipped_record_total
    }

    /// Number of chunks fully decoded so far.
    pub fn decoded_chunks(&self) -> u64 {
        self.decoded_chunks
    }

    fn resident_chunks(&self) -> usize {
        self.outstanding_jobs + self.pending.len() + usize::from(self.current.is_some())
    }

    /// Records the memory estimate of a decoded chunk (in flight or pending).
    fn register_chunk_bytes(&mut self, chunk_index: u64, estimate: u64) {
        self.resident_bytes = self.resident_bytes.saturating_add(estimate);
        self.resident_estimates.insert(chunk_index, estimate);
    }

    /// Releases the memory estimate once a chunk has been fully consumed.
    fn release_chunk_bytes(&mut self, chunk_index: u64) {
        if let Some(estimate) = self.resident_estimates.remove(&chunk_index) {
            self.resident_bytes = self.resident_bytes.saturating_sub(estimate);
        }
    }

    /// Drops all residency accounting after pending chunks are discarded.
    fn drop_resident_state(&mut self) {
        self.resident_bytes = 0;
        self.resident_estimates.clear();
    }

    fn stop_input_with_error(&mut self, chunk_index: u64, error: Error) {
        self.pending.insert(chunk_index, Err(error));
        self.terminal_chunk = chunk_index + 1;
        self.chunks_remaining = 0;
        self.input_finished = true;
    }

    fn finish_input(&mut self) {
        if self.input_finished {
            return;
        }
        let mut trailing = [0; 1];
        match self.input.read(&mut trailing) {
            Ok(0) => {
                self.terminal_chunk = self.next_submit_chunk;
                self.input_finished = true;
            }
            Ok(_) => self.stop_input_with_error(
                self.next_submit_chunk,
                Error::InvalidFormat("trailing aligned archive bytes".into()),
            ),
            Err(error) => self.stop_input_with_error(self.next_submit_chunk, error.into()),
        }
    }

    fn submit_available(&mut self) {
        while !self.input_finished
            && self.outstanding_jobs < self.target_jobs
            && self.resident_chunks() < self.max_resident_chunks
            && self.resident_bytes < self.max_resident_bytes
        {
            if self.chunks_remaining == 0 {
                self.finish_input();
                break;
            }
            let chunk_index = self.next_submit_chunk;
            let first_record = self.next_submit_record;
            let header = match read_chunk_header(&mut self.input) {
                Ok(header) => header,
                Err(error) => {
                    self.stop_input_with_error(chunk_index, error);
                    break;
                }
            };
            let compressed =
                match read_chunk_payload(&mut self.input, &header, first_record, chunk_index) {
                    Ok(compressed) => compressed,
                    Err(error) => {
                        if self.best_effort && error.is_checksum() {
                            eprintln!("[seqarc] warning: skipping chunk {chunk_index}: {error}");
                            self.skipped_chunks.insert(chunk_index, first_record);
                            self.skipped_chunk_total += 1;
                            self.skipped_record_total += header.record_count as u64;
                            self.next_submit_chunk += 1;
                            self.next_submit_record = self
                                .next_submit_record
                                .saturating_add(header.record_count as u64);
                            self.chunks_remaining -= 1;
                            continue;
                        }
                        self.stop_input_with_error(chunk_index, error);
                        break;
                    }
                };
            let record_count = header.record_count;
            let raw_len = header.raw_len;
            let ids_len = header.ids_len;
            let qualities_len = header.qualities_len;
            let records_len = header.records_len;
            // Account for the chunk's decoded footprint (per-record Vec
            // materialization roughly quadruples the raw stream size).
            let estimate = u64::try_from(raw_len)
                .unwrap_or(u64::MAX)
                .saturating_mul(4)
                .saturating_add(u64::try_from(compressed.len()).unwrap_or(u64::MAX));
            let encoded = EncodedArchiveChunk {
                first_record,
                record_count,
                raw_len,
                ids_len,
                qualities_len,
                records_len,
                compressed,
            };
            let next_record = match first_record.checked_add(encoded.record_count as u64) {
                Some(next_record) => next_record,
                None => {
                    self.stop_input_with_error(
                        chunk_index,
                        Error::InvalidFormat("record count overflow".into()),
                    );
                    break;
                }
            };
            self.next_submit_chunk += 1;
            self.next_submit_record = next_record;
            self.chunks_remaining -= 1;
            self.register_chunk_bytes(chunk_index, estimate);
            if let (Some(pool), Some(sender)) = (&self.pool, &self.result_tx) {
                self.outstanding_jobs += 1;
                let sender = sender.clone();
                let reference = Arc::clone(&self.reference);
                let id_mode = self.id_mode;
                let layout = self.layout;
                let preserve_quality = self.preserve_quality;
                let indexes = Arc::clone(&self.indexes);
                let absolute_position_width = self.absolute_position_width;
                let delta_position_width = self.delta_position_width;
                let index_streams = self.index_streams;
                pool.spawn(move || {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        decode_archive_chunk(
                            encoded,
                            id_mode,
                            layout,
                            preserve_quality,
                            &reference,
                            &indexes,
                            absolute_position_width,
                            delta_position_width,
                            index_streams,
                        )
                    }))
                    .unwrap_or_else(|_| Err(Error::InvalidInput("chunk decoder panicked".into())));
                    let _ = sender.send((chunk_index, result));
                });
            } else {
                let result = decode_archive_chunk(
                    encoded,
                    self.id_mode,
                    self.layout,
                    self.preserve_quality,
                    &self.reference,
                    &self.indexes,
                    self.absolute_position_width,
                    self.delta_position_width,
                    self.index_streams,
                );
                if result.is_ok() {
                    self.decoded_chunks += 1;
                }
                self.pending.insert(chunk_index, result);
            }
            debug_assert!(self.resident_chunks() <= self.max_resident_chunks);
        }
    }

    fn harvest_ready_results(&mut self) -> Result<()> {
        loop {
            let result = match self.result_rx.as_ref() {
                Some(receiver) => receiver.try_recv(),
                None => return Ok(()),
            };
            match result {
                Ok((chunk_index, result)) => {
                    self.outstanding_jobs = self
                        .outstanding_jobs
                        .checked_sub(1)
                        .ok_or_else(|| Error::InvalidFormat("unexpected decoded chunk".into()))?;
                    if result.is_ok() {
                        self.decoded_chunks += 1;
                    }
                    self.pending.insert(chunk_index, result);
                    debug_assert!(self.resident_chunks() <= self.max_resident_chunks);
                }
                Err(crossbeam_channel::TryRecvError::Empty) => return Ok(()),
                Err(crossbeam_channel::TryRecvError::Disconnected) => {
                    return if self.outstanding_jobs == 0 {
                        Ok(())
                    } else {
                        Err(Error::InvalidFormat("chunk decoder exited".into()))
                    };
                }
            }
        }
    }

    fn maintain_parallelism(&mut self) -> Result<()> {
        self.harvest_ready_results()?;
        self.submit_available();
        Ok(())
    }

    fn load_next_chunk(&mut self) -> Result<bool> {
        loop {
            self.maintain_parallelism()?;
            if let Some(result) = self.pending.remove(&self.next_yield_chunk) {
                self.next_yield_chunk += 1;
                self.current = Some(result?);
                debug_assert!(self.resident_chunks() <= self.max_resident_chunks);
                return Ok(true);
            }
            if self.skipped_chunks.contains_key(&self.next_yield_chunk) {
                self.next_yield_chunk += 1;
                continue;
            }
            if self.input_finished && self.outstanding_jobs == 0 {
                if self.next_yield_chunk >= self.terminal_chunk {
                    return Ok(false);
                }
                return Err(Error::InvalidFormat("missing decoded chunk".into()));
            }
            let (chunk_index, result) = self
                .result_rx
                .as_ref()
                .ok_or_else(|| Error::InvalidFormat("missing chunk decoder".into()))?
                .recv()
                .map_err(|_| Error::InvalidFormat("chunk decoder exited".into()))?;
            self.outstanding_jobs = self
                .outstanding_jobs
                .checked_sub(1)
                .ok_or_else(|| Error::InvalidFormat("unexpected decoded chunk".into()))?;
            if result.is_ok() {
                self.decoded_chunks += 1;
            }
            self.pending.insert(chunk_index, result);
            debug_assert!(self.resident_chunks() <= self.max_resident_chunks);
        }
    }
}

impl DecodedArchiveChunk {
    fn next_group(
        &mut self,
        layout: AlignmentLayout,
        range: Option<&RangeInclusive<u64>>,
    ) -> Option<(u64, AlignedRecordGroup)> {
        while self.next_record < self.ids[0].len() {
            let index = self.next_record;
            self.next_record += 1;
            let original_id = self.first_record + index as u64;
            if range.is_some_and(|range| !range.contains(&original_id)) {
                continue;
            }
            let id = std::mem::take(&mut self.ids[0][index]);
            let index_reads = self
                .index_reads
                .iter_mut()
                .map(|stream| {
                    ReadRecord::new(
                        id.clone(),
                        std::mem::take(&mut stream.sequences[index]),
                        stream.qualities[index].take(),
                    )
                })
                .collect();
            let first = ReadRecord::new(
                id,
                std::mem::take(&mut self.sequences[0][index]),
                self.qualities[0][index].take(),
            );
            let sequence = if layout == AlignmentLayout::SingleEnd {
                SequenceRecordGroup::Single(first)
            } else {
                SequenceRecordGroup::Paired(
                    first,
                    ReadRecord::new(
                        std::mem::take(&mut self.ids[1][index]),
                        std::mem::take(&mut self.sequences[1][index]),
                        self.qualities[1][index].take(),
                    ),
                )
            };
            return Some((
                original_id,
                AlignedRecordGroup {
                    sequence,
                    index_reads,
                },
            ));
        }
        None
    }
}

impl<R: Read> Iterator for AlignedRecordIter<R> {
    type Item = Result<(u64, AlignedRecordGroup)>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if self.done {
                return None;
            }
            if let Err(error) = self.maintain_parallelism() {
                self.done = true;
                self.current = None;
                self.pending.clear();
                self.drop_resident_state();
                return Some(Err(error));
            }
            if let Some(current) = &mut self.current {
                if let Some(record) = current.next_group(self.layout, self.original_range.as_ref())
                {
                    return Some(Ok(record));
                }
                let exhausted_chunk = self.next_yield_chunk - 1;
                self.current = None;
                self.release_chunk_bytes(exhausted_chunk);
                debug_assert!(self.resident_chunks() <= self.max_resident_chunks);
            }
            match self.load_next_chunk() {
                Ok(true) => {}
                Ok(false) => {
                    self.done = true;
                    return None;
                }
                Err(error) => {
                    self.done = true;
                    self.current = None;
                    self.pending.clear();
                    self.drop_resident_state();
                    return Some(Err(error));
                }
            }
        }
    }
}

impl<R: Read> std::iter::FusedIterator for AlignedRecordIter<R> {}

struct ChunkHeader {
    bytes: [u8; 44],
    record_count: usize,
    compressed_len: usize,
    raw_len: usize,
    ids_len: usize,
    qualities_len: usize,
    records_len: usize,
    stored_crc: u32,
}

fn read_chunk_header<R: Read>(input: &mut R) -> Result<ChunkHeader> {
    let mut bytes = [0u8; 44];
    input.read_exact(&mut bytes)?;
    let mut fields = Cursor::new(&bytes[..]);
    let record_count = read_u32(&mut fields)? as usize;
    let compressed_len = usize_from(read_u64(&mut fields)?, "chunk compressed size")?;
    let raw_len = usize_from(read_u64(&mut fields)?, "chunk uncompressed size")?;
    let ids_len = usize_from(read_u64(&mut fields)?, "chunk ID block size")?;
    let qualities_len = usize_from(read_u64(&mut fields)?, "chunk quality section size")?;
    let records_len = usize_from(read_u64(&mut fields)?, "chunk record block size")?;
    let stored_crc = read_u32(input)?;
    Ok(ChunkHeader {
        bytes,
        record_count,
        compressed_len,
        raw_len,
        ids_len,
        qualities_len,
        records_len,
        stored_crc,
    })
}

fn read_chunk_payload<R: Read>(
    input: &mut R,
    header: &ChunkHeader,
    first_record: u64,
    chunk_index: u64,
) -> Result<Vec<u8>> {
    // Read the payload incrementally while hashing, so a corrupt
    // `compressed_len` cannot trigger a huge allocation up front; memory grows
    // only as bytes actually arrive. The initial capacity is capped so even
    // the reservation stays bounded.
    let mut compressed = Vec::with_capacity(header.compressed_len.min(4 * 1024 * 1024));
    let mut crc = crc32c::crc32c(&header.bytes);
    let mut buffer = [0u8; 64 * 1024];
    let mut remaining = header.compressed_len;
    while remaining > 0 {
        let want = remaining.min(buffer.len());
        let read = input.read(&mut buffer[..want])?;
        if read == 0 {
            return Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!(
                    "chunk {chunk_index} (record {first_record}): payload ended after {} of {} bytes",
                    compressed.len(),
                    header.compressed_len
                ),
            )));
        }
        crc = crc32c::crc32c_append(crc, &buffer[..read]);
        compressed.extend_from_slice(&buffer[..read]);
        remaining -= read;
    }
    if crc != header.stored_crc {
        return Err(Error::Checksum(format!(
            "chunk {chunk_index} (record {first_record}): stored CRC 0x{:08x} does not match \
             computed CRC 0x{:08x}",
            header.stored_crc, crc
        )));
    }
    Ok(compressed)
}

/// Result of a checksum-only archive scan.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct VerifyReport {
    /// Chunk count declared by the global header.
    pub declared_chunks: u64,
    /// Chunks whose 48-byte header and payload were successfully read.
    pub scanned_chunks: u64,
    /// Records contained in chunks whose checksum matched.
    pub records: u64,
    /// Chunks whose stored CRC32C did not match the computed CRC32C.
    pub bad_chunks: u64,
    /// Ordinal of the first chunk with a checksum mismatch, if any.
    pub first_bad_chunk: Option<u64>,
    /// True when the archive ended before all declared chunks were scanned.
    pub truncated: bool,
    /// True when bytes follow the declared final chunk.
    pub trailing_bytes: bool,
}

impl VerifyReport {
    /// True when no damage was detected.
    pub fn is_ok(&self) -> bool {
        self.bad_chunks == 0 && !self.truncated && !self.trailing_bytes
    }
}

/// Verify the global header, chunk framing, and per-chunk CRC32C checksums
/// without decompressing chunk payloads.
pub fn verify_aligned_archive<R: Read>(input: &mut R) -> Result<VerifyReport> {
    let header = read_header(input)?;
    let mut report = VerifyReport {
        declared_chunks: header.chunk_count,
        ..VerifyReport::default()
    };
    while report.scanned_chunks < header.chunk_count {
        let chunk_index = report.scanned_chunks;
        let chunk_header = match read_chunk_header(input) {
            Ok(chunk_header) => chunk_header,
            // A clean early EOF is truncation; any other error (I/O failure,
            // permissions) is reported to the caller instead of being
            // mislabeled as a truncated archive.
            Err(Error::Io(error)) if error.kind() == std::io::ErrorKind::UnexpectedEof => {
                report.truncated = true;
                return Ok(report);
            }
            Err(error) => return Err(error),
        };
        let mut crc = crc32c::crc32c(&chunk_header.bytes);
        let mut remaining = chunk_header.compressed_len;
        let mut buffer = [0u8; 64 * 1024];
        let mut truncated = false;
        while remaining > 0 {
            let want = remaining.min(buffer.len());
            match input.read(&mut buffer[..want]) {
                Ok(0) => {
                    truncated = true;
                    break;
                }
                Ok(read) => {
                    crc = crc32c::crc32c_append(crc, &buffer[..read]);
                    remaining -= read;
                }
                Err(error) => return Err(Error::Io(error)),
            }
        }
        if truncated {
            report.truncated = true;
            return Ok(report);
        }
        if crc == chunk_header.stored_crc {
            report.records += chunk_header.record_count as u64;
        } else {
            report.bad_chunks += 1;
            if report.first_bad_chunk.is_none() {
                report.first_bad_chunk = Some(chunk_index);
            }
        }
        report.scanned_chunks += 1;
    }
    let mut trailing = [0u8; 1];
    report.trailing_bytes = input.read(&mut trailing)? > 0;
    Ok(report)
}

/// Compressed and decompressed byte totals for one archive section.
///
/// `compressed` is what the section occupies on disk; `raw` is the total
/// length after decompression when the section stores compressed data.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SectionSizes {
    /// Bytes the section occupies on disk.
    pub compressed: u64,
    /// Total length after decompression.
    pub raw: u64,
}

impl std::ops::AddAssign for SectionSizes {
    fn add_assign(&mut self, rhs: Self) {
        self.compressed += rhs.compressed;
        self.raw += rhs.raw;
    }
}

/// The eight per-mate record subsections, in on-disk order.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MateSubsectionSizes {
    /// Read length column.
    pub lengths: SectionSizes,
    /// Reference position column.
    pub positions: SectionSizes,
    /// Mismatch/emit position column.
    pub noise_positions: SectionSizes,
    /// Mismatch symbol column.
    pub noise: SectionSizes,
    /// Reference skip (deletion) column.
    pub reference_skips: SectionSizes,
    /// Inserted-bases column.
    pub emits: SectionSizes,
    /// Strand-flags column.
    pub reverses: SectionSizes,
    /// Raw fallback sequences for unaligned reads.
    pub unaligned_sequences: SectionSizes,
}

/// The three record subsections per declared index-read stream, in on-disk
/// order.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct IndexStreamSubsectionSizes {
    /// Quality column.
    pub qualities: SectionSizes,
    /// Read length column.
    pub lengths: SectionSizes,
    /// Sequence column.
    pub sequences: SectionSizes,
}

/// Per-subsection record section sizes, in on-disk order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RecordSubsectionSizes {
    /// Per-record flags column.
    pub flags: SectionSizes,
    /// One entry per biological stream (one for single-end, two for paired).
    pub mates: Vec<MateSubsectionSizes>,
    /// One entry per embedded index range, in declaration order.
    pub embedded_index_ranges: Vec<SectionSizes>,
    /// One entry per declared index-read stream, in I1, I2 order.
    pub index_streams: Vec<IndexStreamSubsectionSizes>,
}

/// Aggregated on-disk storage sizes for every archive section.
///
/// Collected by [`read_archive_info_with_breakdown`] by walking chunk payload
/// directories; no payload is ever decompressed. The fixed preamble, chunk
/// headers, and reference/index-metadata blocks are included so the totals
/// reconcile with the archive file size.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SectionBreakdown {
    /// Fixed framing before the first chunk: the preamble fields plus the
    /// preamble CRC.
    pub preamble: u64,
    /// Chunk framing overhead: the 44 header fields plus the CRC per chunk.
    pub chunk_headers: u64,
    /// The 24-byte compressed section directory per chunk payload.
    pub payload_directories: u64,
    /// Identifier section sizes.
    pub ids: SectionSizes,
    /// Quality section sizes.
    pub qualities: SectionSizes,
    /// `compressed` is the verbatim record section bytes on disk; `raw` is
    /// the sum of the decompressed subsection lengths.
    pub records: SectionSizes,
    /// Compressed size of the embedded reference block.
    pub reference_block: u64,
    /// Size of the embedded index metadata block.
    pub index_metadata: u64,
    /// Per-subsection record sizes.
    pub subsections: RecordSubsectionSizes,
}

impl SectionBreakdown {
    /// Total on-disk bytes covered by the breakdown. For a complete,
    /// untruncated archive this equals the file size.
    pub fn accounted_bytes(&self) -> u64 {
        self.preamble
            + self.chunk_headers
            + self.payload_directories
            + self.ids.compressed
            + self.qualities.compressed
            + self.records.compressed
            + self.reference_block
            + self.index_metadata
    }
}

/// High-level archive summary readable without decoding chunk payloads.
///
/// The preamble (including its CRC) and every chunk header are read, but
/// chunk payloads are skipped without checksum verification; use
/// [`verify_aligned_archive`] for a full integrity check.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArchiveInfo {
    /// Major format version.
    pub version_major: u16,
    /// Minor format version.
    pub version_minor: u16,
    /// Single- or paired-end layout.
    pub layout: AlignmentLayout,
    /// Identifier storage mode: "none", "independent", "identical-mates", or
    /// "slash-mate".
    pub id_mode: &'static str,
    /// Whether identifiers are stored.
    pub preserve_ids: bool,
    /// Whether quality values are stored.
    pub preserve_quality: bool,
    /// Number of declared index-read streams.
    pub index_streams: usize,
    /// Embedded index read ranges in one-based inclusive coordinates.
    pub embedded_indexes: Vec<EmbeddedIndexRange>,
    /// Byte width of the absolute-position column.
    pub absolute_position_width: usize,
    /// Byte width of the delta-position column.
    pub delta_position_width: usize,
    /// Chunk count declared by the global header.
    pub declared_chunks: u64,
    /// Chunk headers successfully read; less than `declared_chunks` when the
    /// file is truncated.
    pub scanned_chunks: u64,
    /// Records declared by the scanned chunk headers.
    pub records: u64,
    /// Compressed size of the embedded reference block in bytes.
    pub reference_block_bytes: u64,
    /// Size of the embedded index metadata block in bytes.
    pub index_metadata_bytes: u64,
    /// True when the file ended before all declared chunk headers were read.
    pub truncated: bool,
    /// True when bytes follow the declared final chunk.
    pub trailing_bytes: bool,
    /// Per-section storage sizes, collected only by
    /// [`read_archive_info_with_breakdown`].
    pub breakdown: Option<SectionBreakdown>,
}

/// Read archive metadata from the header and chunk headers without decoding
/// any chunk payloads.
pub fn read_archive_info<R: Read + Seek>(input: &mut R) -> Result<ArchiveInfo> {
    read_archive_info_with_breakdown(input, false)
}

/// Like [`read_archive_info`], additionally walking each chunk payload's
/// section directory and the record section's subsection directory to
/// aggregate per-section storage sizes. Payloads are skipped with seeks and
/// never decompressed, so the extra cost is a handful of small reads per
/// chunk.
pub fn read_archive_info_with_breakdown<R: Read + Seek>(
    input: &mut R,
    collect_breakdown: bool,
) -> Result<ArchiveInfo> {
    let parsed = parse_preamble(input)?;
    let descriptors = decode_index_metadata(&parsed.index_metadata, parsed.layout)?;
    let embedded_indexes = descriptors.iter().map(|d| d.range.clone()).collect();
    let id_mode = match parsed.id_mode {
        IdMode::None => "none",
        IdMode::Independent => "independent",
        IdMode::IdenticalMates => "identical-mates",
        IdMode::SlashMate => "slash-mate",
    };
    let breakdown = collect_breakdown.then(|| SectionBreakdown {
        preamble: (PREAMBLE_BYTES + 4) as u64,
        reference_block: parsed.reference_block.len() as u64,
        index_metadata: parsed.index_metadata.len() as u64,
        subsections: RecordSubsectionSizes {
            mates: vec![MateSubsectionSizes::default(); parsed.layout.streams()],
            embedded_index_ranges: vec![SectionSizes::default(); descriptors.len()],
            index_streams: vec![IndexStreamSubsectionSizes::default(); parsed.index_streams],
            ..RecordSubsectionSizes::default()
        },
        ..SectionBreakdown::default()
    });
    let mut info = ArchiveInfo {
        version_major: parsed.version_major,
        version_minor: parsed.version_minor,
        layout: parsed.layout,
        id_mode,
        preserve_ids: parsed.flags & 1 != 0,
        preserve_quality: parsed.flags & 2 != 0,
        index_streams: parsed.index_streams,
        embedded_indexes,
        absolute_position_width: parsed.absolute_position_width,
        delta_position_width: parsed.delta_position_width,
        declared_chunks: parsed.chunk_count,
        scanned_chunks: 0,
        records: 0,
        reference_block_bytes: parsed.reference_block.len() as u64,
        index_metadata_bytes: parsed.index_metadata.len() as u64,
        truncated: false,
        trailing_bytes: false,
        breakdown,
    };
    // Walk the chunk headers, skipping each payload with a seek. Chunk
    // checksums are not verified here; `verify` covers that.
    while info.scanned_chunks < info.declared_chunks {
        let chunk_header = match read_chunk_header(input) {
            Ok(chunk_header) => chunk_header,
            Err(Error::Io(error)) if error.kind() == std::io::ErrorKind::UnexpectedEof => {
                info.truncated = true;
                break;
            }
            Err(error) => return Err(error),
        };
        info.records += chunk_header.record_count as u64;
        info.scanned_chunks += 1;
        // A corrupt length that cannot be used as a seek offset leaves the
        // framing unusable; stop and report what was scanned.
        let Ok(offset) = i64::try_from(chunk_header.compressed_len) else {
            info.truncated = true;
            break;
        };
        let payload = if let Some(breakdown) = info.breakdown.as_mut() {
            read_chunk_payload_sizes(
                input,
                &chunk_header,
                parsed.layout,
                descriptors.len(),
                parsed.index_streams,
                breakdown,
            )
        } else {
            input
                .seek(SeekFrom::Current(offset))
                .map(|_| ())
                .map_err(Error::Io)
        };
        match payload {
            Ok(()) => {}
            Err(Error::Io(error)) if error.kind() == std::io::ErrorKind::UnexpectedEof => {
                info.truncated = true;
                break;
            }
            Err(error) => return Err(error),
        }
    }
    let mut trailing = [0u8; 1];
    info.trailing_bytes = input.read(&mut trailing)? > 0;
    Ok(info)
}

/// Aggregate one chunk's on-disk section sizes into `breakdown`.
///
/// Reads the chunk's 24-byte payload directory and the record section's
/// subsection directory, skipping every compressed frame with a seek. The
/// reader must be positioned at the start of the chunk payload.
fn read_chunk_payload_sizes<R: Read + Seek>(
    input: &mut R,
    header: &ChunkHeader,
    layout: AlignmentLayout,
    descriptor_count: usize,
    index_stream_count: usize,
    breakdown: &mut SectionBreakdown,
) -> Result<()> {
    let ids_compressed = read_u64(input)?;
    let qualities_compressed = read_u64(input)?;
    let records_compressed = read_u64(input)?;
    let directory_total = 24u64
        .checked_add(ids_compressed)
        .and_then(|value| value.checked_add(qualities_compressed))
        .and_then(|value| value.checked_add(records_compressed))
        .ok_or_else(|| Error::InvalidFormat("compressed section size overflow".into()))?;
    if directory_total != header.compressed_len as u64 {
        return Err(Error::InvalidFormat(
            "chunk payload section sizes do not match the chunk header".into(),
        ));
    }
    // The record section is stored verbatim; the payload directory and the
    // chunk header must agree on its length.
    if records_compressed != header.records_len as u64 {
        return Err(Error::InvalidFormat(
            "record section length disagrees with the chunk header".into(),
        ));
    }
    breakdown.chunk_headers += 48;
    breakdown.payload_directories += 24;
    breakdown.ids.compressed += ids_compressed;
    breakdown.ids.raw += header.ids_len as u64;
    breakdown.qualities.compressed += qualities_compressed;
    breakdown.qualities.raw += header.qualities_len as u64;
    breakdown.records.compressed += records_compressed;
    skip_forward(input, ids_compressed)?;
    skip_forward(input, qualities_compressed)?;

    // The record section is a directory of independently compressed
    // subsection frames; read only each length pair and skip its frame.
    let subsection_count = read_u32(input)? as usize;
    let expected_subsections =
        1 + layout.streams() * V21_MATE_SECTION_COUNT + descriptor_count + index_stream_count * 3;
    if subsection_count != expected_subsections {
        return Err(Error::InvalidFormat(
            "record subsection count mismatch".into(),
        ));
    }
    let subsections = &mut breakdown.subsections;
    let embedded_start = 1 + layout.streams() * V21_MATE_SECTION_COUNT;
    let index_stream_start = embedded_start + descriptor_count;
    let mut subsection_compressed = 0u64;
    let mut subsection_raw = 0u64;
    for index in 0..subsection_count {
        let compressed_len = read_u64(input)?;
        let raw_len = read_u64(input)?;
        subsection_compressed += compressed_len;
        subsection_raw += raw_len;
        let sizes = SectionSizes {
            compressed: compressed_len,
            raw: raw_len,
        };
        if index == 0 {
            subsections.flags += sizes;
        } else if index < embedded_start {
            let entry = &mut subsections.mates[(index - 1) / V21_MATE_SECTION_COUNT];
            match (index - 1) % V21_MATE_SECTION_COUNT {
                0 => entry.lengths += sizes,
                1 => entry.positions += sizes,
                2 => entry.noise_positions += sizes,
                3 => entry.noise += sizes,
                4 => entry.reference_skips += sizes,
                5 => entry.emits += sizes,
                6 => entry.reverses += sizes,
                _ => entry.unaligned_sequences += sizes,
            }
        } else if index < index_stream_start {
            subsections.embedded_index_ranges[index - embedded_start] += sizes;
        } else {
            let entry = &mut subsections.index_streams[(index - index_stream_start) / 3];
            match (index - index_stream_start) % 3 {
                0 => entry.qualities += sizes,
                1 => entry.lengths += sizes,
                _ => entry.sequences += sizes,
            }
        }
        skip_forward(input, compressed_len)?;
    }
    // The directory and its frames must consume the record section exactly.
    let subsection_total = 4u64 + subsection_count as u64 * 16 + subsection_compressed;
    if subsection_total != records_compressed {
        return Err(Error::InvalidFormat(
            "record subsection directory does not consume the record section".into(),
        ));
    }
    breakdown.records.raw += subsection_raw;
    Ok(())
}

/// Seek forward over `count` payload bytes without reading them. Seeking past
/// EOF is permitted on files; the next chunk header read then reports the
/// truncation.
fn skip_forward<R: Read + Seek>(input: &mut R, count: u64) -> Result<()> {
    let offset = i64::try_from(count)
        .map_err(|_| Error::InvalidFormat("payload section size exceeds seek range".into()))?;
    input.seek(SeekFrom::Current(offset))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn decode_archive_chunk(
    encoded: EncodedArchiveChunk,
    id_mode: IdMode,
    layout: AlignmentLayout,
    preserve_quality: bool,
    reference: &[u8],
    index_descriptors: &[IndexDescriptor],
    absolute_position_width: usize,
    delta_position_width: usize,
    index_stream_count: usize,
) -> Result<DecodedArchiveChunk> {
    let (raw, records_section) = decode_chunk_sections(
        &encoded.compressed,
        encoded.ids_len,
        encoded.qualities_len,
        encoded.records_len,
    )?;
    let expected_raw_len = encoded
        .ids_len
        .checked_add(encoded.qualities_len)
        .and_then(|length| length.checked_add(encoded.records_len))
        .ok_or_else(|| Error::InvalidFormat("chunk block size overflow".into()))?;
    if expected_raw_len != encoded.raw_len {
        return Err(Error::InvalidFormat(
            "chunk block sizes do not match".into(),
        ));
    }
    let ids_end = encoded.ids_len;
    let qualities_end = ids_end + encoded.qualities_len;
    let qualities = decode_qualities(
        &raw[ids_end..qualities_end],
        encoded.record_count,
        preserve_quality,
        layout,
    )?;
    let (sequences, index_reads) = {
        let sections = decode_v21_records(
            records_section,
            encoded.record_count,
            reference,
            layout,
            index_descriptors,
            preserve_quality,
            index_stream_count,
            absolute_position_width,
            delta_position_width,
        )?;
        (sections.biological, sections.index_streams)
    };
    for mate in 0..layout.streams() {
        for record in 0..encoded.record_count {
            if qualities[mate][record]
                .as_ref()
                .is_some_and(|quality| quality.len() != sequences[mate][record].len())
            {
                return Err(Error::InvalidFormat(
                    "reconstructed sequence and quality lengths differ".into(),
                ));
            }
        }
    }
    for stream in &index_reads {
        for record in 0..encoded.record_count {
            if stream.qualities[record]
                .as_ref()
                .is_some_and(|quality| quality.len() != stream.sequences[record].len())
            {
                return Err(Error::InvalidFormat(
                    "reconstructed sequence and quality lengths differ".into(),
                ));
            }
        }
    }
    Ok(DecodedArchiveChunk {
        first_record: encoded.first_record,
        next_record: 0,
        ids: decode_ids(&raw[..ids_end], encoded.record_count, id_mode, layout)?,
        qualities,
        sequences,
        index_reads,
    })
}

thread_local! {
    /// Reused zstd decompression context per worker thread, avoiding context
    /// reallocation for every decoded section.
    static DECOMPRESSOR: std::cell::RefCell<zstd::zstd_safe::DCtx<'static>> =
        std::cell::RefCell::new(zstd::zstd_safe::DCtx::create());
}

/// Decompresses a single zstd frame into `destination` without intermediate
/// allocations. Frames carrying an embedded content size (all frames written
/// by the current encoder) are decoded one-shot; legacy streaming frames
/// without one fall back to the streaming decoder.
fn decompress_frame_into(source: &[u8], destination: &mut [u8]) -> Result<()> {
    DECOMPRESSOR.with(|context| {
        let mut context = context.borrow_mut();
        let content_size = zstd::zstd_safe::get_frame_content_size(source)
            .map_err(|_| Error::InvalidFormat("invalid zstd frame header".into()))?;
        if content_size == Some(destination.len() as u64) {
            context
                .decompress(destination, source)
                .map_err(|error| Error::InvalidFormat(format!("zstd decompression failed: {error}")))?;
            Ok(())
        } else {
            let decoded = zstd::stream::decode_all(Cursor::new(source))?;
            if decoded.len() != destination.len() {
                return Err(Error::InvalidFormat("zstd frame length mismatch".into()));
            }
            destination.copy_from_slice(&decoded);
            Ok(())
        }
    })
}

/// Decompresses the chunk ID and quality sections in place and returns them
/// alongside a borrowed view of the (already framed) record subsections. The
/// returned buffer holds ids followed by qualities; the record section is not
/// copied.
fn decode_chunk_sections(
    compressed: &[u8],
    ids_len: usize,
    qualities_len: usize,
    records_len: usize,
) -> Result<(Vec<u8>, &[u8])> {
    let mut header = Cursor::new(compressed);
    let ids_compressed_len = usize_from(read_u64(&mut header)?, "compressed ID section")?;
    let qualities_compressed_len =
        usize_from(read_u64(&mut header)?, "compressed quality section")?;
    let records_compressed_len = usize_from(read_u64(&mut header)?, "compressed record section")?;
    let expected = 24usize
        .checked_add(ids_compressed_len)
        .and_then(|value| value.checked_add(qualities_compressed_len))
        .and_then(|value| value.checked_add(records_compressed_len))
        .ok_or_else(|| Error::InvalidFormat("compressed section size overflow".into()))?;
    if expected != compressed.len() {
        return Err(Error::InvalidFormat(
            "compressed section sizes do not match".into(),
        ));
    }
    let ids_source = &compressed[24..24 + ids_compressed_len];
    let qualities_source = &compressed[24 + ids_compressed_len..expected - records_len];
    let records_source = &compressed[expected - records_len..];
    if records_source.len() != records_len {
        return Err(Error::InvalidFormat(
            "chunk record section length mismatch".into(),
        ));
    }
    let mut raw = vec![0u8; ids_len + qualities_len];
    decompress_frame_into(ids_source, &mut raw[..ids_len])?;
    decompress_frame_into(qualities_source, &mut raw[ids_len..])?;
    Ok((raw, records_source))
}

fn map_pairs(
    index: &AlignmentIndex,
    pairs: &mut [AlignedPair],
    mapping_pool: &rayon::ThreadPool,
    progress: &ProgressBar,
    progress_stats: &mut MappingProgressStats,
) -> Result<()> {
    let first_records = pairs.iter().map(|pair| pair.record(0)).collect::<Vec<_>>();
    let first_hits = aligners::map_records(
        index,
        &first_records,
        mapping_pool,
        &alignment_candidate_cost,
    )?;
    progress_stats.observe(0, &first_hits, progress);
    for (pair, first) in pairs.iter_mut().zip(first_hits) {
        pair.first_hit = first.hit;
        pair.first_selection = first.selection;
    }
    if pairs.first().is_some_and(|pair| pair.second.is_some()) {
        let second_records = pairs.iter().map(|pair| pair.record(1)).collect::<Vec<_>>();
        let second_hits = aligners::map_records(
            index,
            &second_records,
            mapping_pool,
            &alignment_candidate_cost,
        )?;
        progress_stats.observe(1, &second_hits, progress);
        for (pair, second) in pairs.iter_mut().zip(second_hits) {
            pair.second_hit = second.hit;
            pair.second_selection = second.selection;
        }
    }
    Ok(())
}

fn update_mapping_profile(profile: &mut AlignedCompressionProfile, pairs: &[AlignedPair]) {
    for pair in pairs {
        for stream in 0..pair.layout().streams() {
            let (hit, selection) = if stream == 0 {
                (&pair.first_hit, pair.first_selection)
            } else {
                (&pair.second_hit, pair.second_selection)
            };
            profile.candidate_hits += u64::from(selection.candidate_hits);
            profile.selected_non_first_hits += u64::from(selection.selected_non_first);
            profile.rejected_alignments += u64::from(selection.rejected_alignment);
            if hit.is_some() {
                profile.mapped_reads += 1;
            } else {
                profile.raw_reads += 1;
            }
        }
    }
}

fn spool_mapped_pairs(
    pairs: &[AlignedPair],
    id_detector: &mut IdModeDetector,
    mask: &mut ReferenceMask,
    spool_builder: &mut ExternalChunkBuilder<SpoolPair>,
) -> Result<()> {
    for pair in pairs {
        id_detector.observe(
            &pair.first.id,
            pair.second.as_ref().map(|record| record.id.as_slice()),
        );
        if let Some(hit) = &pair.first_hit {
            mask.mark_hit(hit)?;
        }
        if let Some(hit) = &pair.second_hit {
            mask.mark_hit(hit)?;
        }
        let spool_pair = SpoolPair::from_records(pair)?;
        spool_builder
            .add(spool_pair)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
    }
    Ok(())
}

impl SpoolPair {
    fn from_records(pair: &AlignedPair) -> Result<Self> {
        Ok(Self {
            first: SpoolRead::from_record(&pair.first, pair.first_hit.clone())?,
            second: pair
                .second
                .as_ref()
                .map(|record| SpoolRead::from_record(record, pair.second_hit.clone()))
                .transpose()?,
            indexes: pair.indexes.clone(),
            index_reads: pair
                .index_reads
                .iter()
                .map(|record| SpoolIndexRead {
                    sequence: record.sequence.clone(),
                    quality: record.quality.clone(),
                })
                .collect::<Vec<_>>(),
        })
    }

    fn into_aligned_pair(self, runs: &[KeptRun], reference_len: usize) -> Result<AlignedPair> {
        Ok(AlignedPair {
            first: self.first.record()?,
            second: self.second.as_ref().map(SpoolRead::record).transpose()?,
            first_hit: self
                .first
                .alignment
                .map(|hit| hit.into_genome(runs, reference_len))
                .transpose()?,
            second_hit: self
                .second
                .and_then(|read| read.alignment)
                .map(|hit| hit.into_genome(runs, reference_len))
                .transpose()?,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: self.indexes,
            index_reads: self
                .index_reads
                .iter()
                .map(|index| {
                    Ok(ReadRecord::new(
                        Vec::new(),
                        index.sequence.clone(),
                        index.quality.clone(),
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
        })
    }
}

impl SpoolRead {
    fn from_record(record: &ReadRecord, hit: Option<GenomeAlignment>) -> Result<Self> {
        Ok(Self {
            id: record.id.clone(),
            sequence: record.sequence.clone(),
            quality: record.quality.clone(),
            alignment: hit.map(SpoolAlignment::from_genome).transpose()?,
        })
    }

    fn record(&self) -> Result<ReadRecord> {
        Ok(ReadRecord::new(
            self.id.clone(),
            self.sequence.clone(),
            self.quality.clone(),
        ))
    }
}

impl SpoolAlignment {
    fn from_genome(hit: GenomeAlignment) -> Result<Self> {
        Ok(Self {
            start: hit.start,
            query_start: hit.query_start,
            query_end: hit.query_end,
            is_reverse: hit.is_reverse,
            cigar: hit
                .cigar
                .into_iter()
                .map(|op| {
                    Ok(SpoolCigarOp {
                        op: cigar_code(op.kind)?,
                        len: op.len,
                    })
                })
                .collect::<Result<Vec<_>>>()?,
        })
    }

    fn into_genome(self, runs: &[KeptRun], reference_len: usize) -> Result<GenomeAlignment> {
        let original_start = nonnegative_usize(self.start, "alignment start")?;
        let start = compact_boundary(runs, original_start, reference_len)?;
        let mut original_position = original_start;
        let mut cigar: Vec<CigarOp> = Vec::with_capacity(self.cigar.len());
        for op in self.cigar {
            let kind = cigar_kind(op.op)?;
            let len = match kind {
                CigarKind::Match | CigarKind::Eq | CigarKind::Diff => {
                    let end = original_position
                        .checked_add(op.len as usize)
                        .ok_or_else(|| Error::InvalidInput("alignment position overflow".into()))?;
                    let compact_start = compact_boundary(runs, original_position, reference_len)?;
                    let compact_end = compact_boundary(runs, end, reference_len)?;
                    let compact_len = compact_end.checked_sub(compact_start).ok_or_else(|| {
                        Error::InvalidFormat("compact alignment position underflow".into())
                    })?;
                    if compact_len != op.len as usize {
                        return Err(Error::InvalidFormat(
                            "aligned reference bases were pruned".into(),
                        ));
                    }
                    original_position = end;
                    op.len
                }
                CigarKind::Del | CigarKind::RefSkip => {
                    let end = original_position
                        .checked_add(op.len as usize)
                        .ok_or_else(|| Error::InvalidInput("alignment position overflow".into()))?;
                    let compact_start = compact_boundary(runs, original_position, reference_len)?;
                    let compact_end = compact_boundary(runs, end, reference_len)?;
                    original_position = end;
                    u32::try_from(compact_end.checked_sub(compact_start).ok_or_else(|| {
                        Error::InvalidFormat("compact alignment position underflow".into())
                    })?)
                    .map_err(|_| Error::InvalidInput("compact CIGAR length exceeds u32".into()))?
                }
                CigarKind::Ins | CigarKind::SoftClip | CigarKind::HardClip | CigarKind::Pad => {
                    op.len
                }
            };
            if len == 0 && matches!(kind, CigarKind::Del | CigarKind::RefSkip) {
                continue;
            }
            if let Some(previous) = cigar.last_mut().filter(|previous| previous.kind == kind) {
                previous.len = previous
                    .len
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("CIGAR length overflow".into()))?;
            } else {
                cigar.push(CigarOp { kind, len });
            }
        }
        Ok(GenomeAlignment {
            start: start as i64,
            query_start: self.query_start,
            query_end: self.query_end,
            is_reverse: self.is_reverse,
            cigar,
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn encode_chunk(
    reference: &[u8],
    pairs: &[AlignedPair],
    absolute_position_width: usize,
    delta_position_width: usize,
    options: &crate::AlignmentCompressionOptions,
    id_mode: IdMode,
    layout: AlignmentLayout,
    index_descriptors: &[IndexDescriptor],
) -> Result<EncodedChunk> {
    let mut ids = Vec::new();
    write_id_stream(&mut ids, pairs, id_mode, layout)?;
    let mut qualities = Vec::new();
    write_quality_stream(&mut qualities, pairs, index_descriptors, options, layout)?;
    let records_stream = write_v21_record_stream(
        reference,
        pairs,
        absolute_position_width,
        delta_position_width,
        layout,
        index_descriptors,
        options.preserve_quality,
        options.index_read_files.len(),
        options.zstd_level,
    )?;
    let ids_compressed = zstd::bulk::compress(&ids, options.zstd_level)?;
    let qualities_compressed = zstd::bulk::compress(&qualities, options.zstd_level)?;
    // `records` already contains independently compressed subsection frames.
    let records_len = records_stream.bytes.len() as u64;
    let records_compressed = records_stream.bytes;
    let indexes_raw_len = records_stream.indexes_raw_len;
    let indexes_compressed_len = records_stream.indexes_compressed_len;
    let mut compressed = Vec::with_capacity(
        24 + ids_compressed.len() + qualities_compressed.len() + records_compressed.len(),
    );
    write_u64(&mut compressed, ids_compressed.len() as u64)?;
    write_u64(&mut compressed, qualities_compressed.len() as u64)?;
    write_u64(&mut compressed, records_compressed.len() as u64)?;
    compressed.extend(&ids_compressed);
    compressed.extend(&qualities_compressed);
    compressed.extend(&records_compressed);
    Ok(EncodedChunk {
        pairs: u32::try_from(pairs.len())
            .map_err(|_| Error::InvalidInput("chunk contains too many records".into()))?,
        raw_len: (ids.len() + qualities.len()) as u64 + records_len,
        ids_len: ids.len() as u64,
        qualities_len: qualities.len() as u64,
        records_len,
        ids_compressed_len: ids_compressed.len() as u64,
        qualities_compressed_len: qualities_compressed.len() as u64,
        records_compressed_len: records_compressed.len() as u64,
        indexes_raw_len,
        indexes_compressed_len,
        compressed,
    })
}

fn write_encoded_chunk(archive: &mut File, chunk: EncodedChunk) -> Result<()> {
    let mut header = Vec::with_capacity(44);
    write_u32(&mut header, chunk.pairs)?;
    write_u64(&mut header, chunk.compressed.len() as u64)?;
    write_u64(&mut header, chunk.raw_len)?;
    write_u64(&mut header, chunk.ids_len)?;
    write_u64(&mut header, chunk.qualities_len)?;
    write_u64(&mut header, chunk.records_len)?;
    debug_assert_eq!(header.len(), 44);
    let crc = chunk_crc32c(&header, &chunk.compressed);
    archive.write_all(&header)?;
    write_u32(archive, crc)?;
    archive.write_all(&chunk.compressed)?;
    Ok(())
}

/// CRC32C over the 44-byte chunk header fields followed by the compressed payload.
fn chunk_crc32c(header: &[u8], payload: &[u8]) -> u32 {
    let crc = crc32c::crc32c(header);
    crc32c::crc32c_append(crc, payload)
}

fn write_id_stream(
    output: &mut Vec<u8>,
    pairs: &[AlignedPair],
    mode: IdMode,
    layout: AlignmentLayout,
) -> Result<()> {
    for mate in 0..layout.streams() {
        let mut block = Vec::new();
        if mode != IdMode::None && !(mode != IdMode::Independent && mate == 1) {
            let ids = pairs.iter().map(|pair| {
                let id = &pair.record(mate).id;
                if mode == IdMode::SlashMate {
                    &id[..id.len().saturating_sub(2)]
                } else {
                    id.as_slice()
                }
            });
            write_prefix_ids(&mut block, ids)?;
        }
        write_u64(output, block.len() as u64)?;
        output.extend(block);
    }
    Ok(())
}

fn write_prefix_ids<'a, I>(output: &mut Vec<u8>, ids: I) -> Result<()>
where
    I: IntoIterator<Item = &'a [u8]>,
{
    let mut previous: Vec<u8> = Vec::new();
    for id in ids {
        let common = previous
            .iter()
            .zip(id)
            .take_while(|(left, right)| left == right)
            .count();
        write_varint(output, common as u64)?;
        write_varint(output, (id.len() - common) as u64)?;
        output.extend(&id[common..]);
        previous.clear();
        previous.extend(id);
    }
    Ok(())
}

fn write_quality_stream(
    output: &mut Vec<u8>,
    pairs: &[AlignedPair],
    descriptors: &[IndexDescriptor],
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<()> {
    for mate in 0..layout.streams() {
        let stream = if options.preserve_quality {
            let mut rows = pairs
                .iter()
                .map(|pair| reconstruct_quality_row(pair, mate, descriptors))
                .collect::<Result<Vec<_>>>()?;
            if let crate::QualityMode::Qvz { ratio } = options.quality_mode {
                // Quantize the whole chunk as one training block, mirroring
                // Spring's per-block QVZ application. The decoder never sees
                // the codebooks; the quantized values are stored losslessly.
                let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
                crate::qvz::quantize(&mut slices, ratio)?;
            }
            let values = rows.iter().map(Vec::as_slice).collect::<Vec<_>>();
            encode_stream_bytes(&values)?
        } else {
            Vec::new()
        };
        write_u64(output, stream.len() as u64)?;
        output.extend(stream);
    }
    Ok(())
}

fn reconstruct_quality_row(
    pair: &AlignedPair,
    mate: usize,
    descriptors: &[IndexDescriptor],
) -> Result<Vec<u8>> {
    if pair.indexes.len() != descriptors.len() {
        return Err(Error::InvalidInput(
            "embedded index count does not match metadata".into(),
        ));
    }
    let biological = pair
        .record(mate)
        .quality
        .as_deref()
        .ok_or_else(|| Error::InvalidInput("biological read quality is missing".into()))?;
    let mut selected = descriptors
        .iter()
        .enumerate()
        .filter(|(_, descriptor)| descriptor.range.mate.stream() == mate)
        .collect::<Vec<_>>();
    selected.sort_by_key(|(_, descriptor)| descriptor.range.start);
    let removed = selected.iter().try_fold(0usize, |total, (_, descriptor)| {
        total
            .checked_add(descriptor.range.len())
            .ok_or_else(|| Error::InvalidInput("embedded index length overflow".into()))
    })?;
    let mut output = Vec::with_capacity(biological.len().saturating_add(removed));
    let mut biological_start = 0usize;
    let mut prior_removed = 0usize;
    for (descriptor_index, descriptor) in selected {
        let boundary = descriptor
            .range
            .start
            .checked_sub(prior_removed)
            .ok_or_else(|| Error::InvalidInput("embedded index boundary underflow".into()))?;
        if boundary < biological_start || boundary > biological.len() {
            return Err(Error::InvalidInput(
                "embedded index quality insertion is out of bounds".into(),
            ));
        }
        output.extend_from_slice(&biological[biological_start..boundary]);
        let quality = pair.indexes[descriptor_index]
            .quality
            .as_deref()
            .ok_or_else(|| Error::InvalidInput("embedded index quality is missing".into()))?;
        if quality.len() != descriptor.range.len() {
            return Err(Error::InvalidInput(
                "embedded index quality length does not match its range".into(),
            ));
        }
        output.extend_from_slice(quality);
        biological_start = boundary;
        prior_removed = prior_removed
            .checked_add(descriptor.range.len())
            .ok_or_else(|| Error::InvalidInput("embedded index length overflow".into()))?;
    }
    output.extend_from_slice(&biological[biological_start..]);
    Ok(output)
}

#[allow(clippy::too_many_arguments)]
fn write_v21_record_stream<R: ReferenceBases + ?Sized>(
    reference: &R,
    pairs: &[AlignedPair],
    absolute_position_width: usize,
    delta_position_width: usize,
    layout: AlignmentLayout,
    descriptors: &[IndexDescriptor],
    preserve_quality: bool,
    index_stream_count: usize,
    zstd_level: i32,
) -> Result<EncodedRecordStream> {
    let mut flags = Vec::with_capacity(pairs.len());
    for pair in pairs {
        flags.push(match layout {
            AlignmentLayout::SingleEnd => {
                if pair.first_hit.is_some() {
                    0
                } else {
                    3
                }
            }
            AlignmentLayout::PairedEnd => {
                match (pair.first_hit.is_some(), pair.second_hit.is_some()) {
                    (true, true) => 0,
                    (true, false) => 1,
                    (false, true) => 2,
                    (false, false) => 3,
                }
            }
        });
    }
    let mut sections =
        Vec::with_capacity(1 + layout.streams() * V21_MATE_SECTION_COUNT + descriptors.len());
    sections.push(flags.clone());
    for mate in 0..layout.streams() {
        let mut lengths = Vec::with_capacity(pairs.len() * 2);
        let mut positions = Vec::new();
        let mut noise_positions = Vec::new();
        let mut noise = Vec::new();
        let mut skips = Vec::new();
        let mut emits = Vec::new();
        let mut reverses = Vec::new();
        let mut unaligned = Vec::new();
        let mut oriented_scratch = Vec::new();
        for (pair, &flag) in pairs.iter().zip(&flags) {
            let hit = pair.hit(mate);
            let aligned = matches!((flag, mate), (0, _) | (1, 0) | (2, 1));
            if aligned != hit.is_some() {
                return Err(Error::InvalidInput(
                    "alignment flag does not match mate state".into(),
                ));
            }
            let record = pair.record(mate);
            write_u16(
                &mut lengths,
                u16::try_from(record.sequence.len())
                    .map_err(|_| Error::InvalidInput("read length exceeds u16 capacity".into()))?,
            )?;
            if let Some(hit) = hit {
                let fields = collect_v21_aligned_fields(
                    reference,
                    record,
                    hit,
                    &mut oriented_scratch,
                )?;
                if mate == 1 && flag == 0 {
                    let first = pair.first_hit.as_ref().ok_or_else(|| {
                        Error::InvalidInput("paired R2 alignment is missing R1 alignment".into())
                    })?;
                    let first_start = first.start;
                    let second_start = hit.start;
                    let delta = second_start.checked_sub(first_start).ok_or_else(|| {
                        Error::InvalidInput("R2 relative alignment position overflow".into())
                    })?;
                    let encoded = zigzag_encode(delta);
                    write_position(&mut positions, encoded, delta_position_width)?;
                    write_u8(&mut reverses, u8::from(fields.reverse == first.is_reverse))?;
                } else {
                    write_position(&mut positions, fields.start, absolute_position_width)?;
                    write_u8(&mut reverses, u8::from(fields.reverse))?;
                }
                for event in fields.events {
                    write_u16(&mut noise_positions, event.delta)?;
                    noise.push(b'0' + event.symbol);
                    if event.symbol == 4 {
                        write_u64(&mut skips, event.skip)?;
                    }
                }
                noise.push(b'\n');
                for emit in fields.emits {
                    if emit.contains(&b'\n') {
                        return Err(Error::InvalidInput("emit sequence contains newline".into()));
                    }
                    emits.extend(emit);
                    emits.push(b'\n');
                }
            } else {
                unaligned.extend_from_slice(&record.sequence);
            }
        }
        sections.extend([
            lengths,
            positions,
            noise_positions,
            noise,
            skips,
            emits,
            reverses,
            unaligned,
        ]);
    }
    for (descriptor_index, descriptor) in descriptors.iter().enumerate() {
        let mut section = Vec::new();
        for pair in pairs {
            let index = pair
                .indexes
                .get(descriptor_index)
                .ok_or_else(|| Error::InvalidInput("embedded index fragment is missing".into()))?;
            if index.sequence.len() != descriptor.range.len() {
                return Err(Error::InvalidInput("index read length mismatch".into()));
            }
            section.extend_from_slice(&index.sequence);
        }
        sections.push(section);
    }
    if pairs
        .iter()
        .any(|pair| pair.index_stream_count() != index_stream_count)
    {
        return Err(Error::InvalidInput(
            "index read count changed during compression".into(),
        ));
    }
    for stream in 0..index_stream_count {
        let mut quality_stream = Vec::new();
        if preserve_quality {
            let values = pairs
                .iter()
                .map(|pair| {
                    let index = pair
                        .index_read(stream)
                        .ok_or_else(|| Error::InvalidInput("index read is missing".into()))?;
                    index
                        .quality
                        .as_deref()
                        .ok_or_else(|| Error::InvalidInput("index read quality is missing".into()))
                })
                .collect::<Result<Vec<_>>>()?;
            quality_stream = encode_stream_bytes(&values)?;
        }
        sections.push(quality_stream);
        let mut lengths = Vec::with_capacity(pairs.len() * 2);
        let mut sequences = Vec::new();
        for pair in pairs {
            let index = pair
                .index_read(stream)
                .ok_or_else(|| Error::InvalidInput("index read is missing".into()))?;
            write_u16(
                &mut lengths,
                u16::try_from(index.sequence.len()).map_err(|_| {
                    Error::InvalidInput("index read length exceeds u16 capacity".into())
                })?,
            )?;
            sequences.extend_from_slice(&index.sequence);
        }
        sections.push(lengths);
        sections.push(sequences);
    }

    let index_start = 1 + layout.streams() * V21_MATE_SECTION_COUNT + descriptors.len();
    let mut output = Vec::new();
    write_u32(
        &mut output,
        u32::try_from(sections.len())
            .map_err(|_| Error::InvalidInput("too many record sections".into()))?,
    )?;
    let mut indexes_raw_len = 0u64;
    let mut indexes_compressed_len = 0u64;
    for (section_index, section) in sections.iter().enumerate() {
        let compressed = zstd::bulk::compress(section, zstd_level)?;
        if section_index >= index_start {
            indexes_raw_len += section.len() as u64;
            indexes_compressed_len += compressed.len() as u64;
        }
        write_u64(&mut output, compressed.len() as u64)?;
        write_u64(&mut output, section.len() as u64)?;
        output.extend(compressed);
    }
    Ok(EncodedRecordStream {
        bytes: output,
        indexes_raw_len,
        indexes_compressed_len,
    })
}

fn collect_v21_aligned_fields<R: ReferenceBases + ?Sized>(
    reference: &R,
    record: &ReadRecord,
    hit: &GenomeAlignment,
    oriented: &mut Vec<u8>,
) -> Result<V21AlignedFields> {
    let start = nonnegative_usize(hit.start, "alignment start")?;
    // The caller reuses this buffer across records; per-read allocation is
    // avoided and the copy doubles as the reverse-complement workspace.
    oriented.clear();
    oriented.extend_from_slice(&record.sequence);
    if hit.is_reverse {
        reverse_complement(oriented);
    }
    let cigar = effective_cigar(hit, record.sequence.len())?;
    let mut events = Vec::new();
    let mut emit_payloads = Vec::new();
    let mut ref_pos = start;
    let mut query_pos = 0usize;
    let mut event_cursor = 0usize;
    for op in cigar {
        let length = op.len as usize;
        match op.kind {
            CigarKind::Match | CigarKind::Eq | CigarKind::Diff => {
                let ref_end = ref_pos
                    .checked_add(length)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                let query_end = query_pos
                    .checked_add(length)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if ref_end > reference.len() || query_end > oriented.len() {
                    return Err(Error::InvalidFormat(
                        "aligned operation exceeds bounds".into(),
                    ));
                }
                for offset in 0..length {
                    let actual = oriented[query_pos + offset];
                    let expected = reference.base(ref_pos + offset);
                    if actual != expected {
                        events.push(V21NoiseEvent {
                            delta: event_delta(query_pos + offset, event_cursor)?,
                            symbol: substitution_symbol(expected, actual)?,
                            skip: 0,
                        });
                        event_cursor = query_pos + offset + 1;
                    }
                }
                ref_pos = ref_end;
                query_pos = query_end;
            }
            CigarKind::Ins | CigarKind::SoftClip => {
                let end = query_pos
                    .checked_add(length)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if end > oriented.len() {
                    return Err(Error::InvalidFormat(
                        "literal operation exceeds query".into(),
                    ));
                }
                events.push(V21NoiseEvent {
                    delta: event_delta(query_pos, event_cursor)?,
                    symbol: 5,
                    skip: 0,
                });
                emit_payloads.push(oriented[query_pos..end].to_vec());
                query_pos = end;
                event_cursor = end;
            }
            CigarKind::Del | CigarKind::RefSkip => {
                ref_pos = ref_pos
                    .checked_add(length)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "reference skip exceeds reference".into(),
                    ));
                }
                events.push(V21NoiseEvent {
                    delta: event_delta(query_pos, event_cursor)?,
                    symbol: 4,
                    skip: length as u64,
                });
                event_cursor = query_pos;
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }
    if query_pos != oriented.len() {
        return Err(Error::InvalidFormat(
            "alignment does not consume read".into(),
        ));
    }
    Ok(V21AlignedFields {
        start: start as u64,
        reverse: hit.is_reverse,
        events,
        emits: emit_payloads,
    })
}

#[cfg(test)]
fn encode_v21_aligned_read<R: ReferenceBases + ?Sized>(
    reference: &R,
    record: &ReadRecord,
    hit: &GenomeAlignment,
) -> Result<Vec<u8>> {
    let fields = collect_v21_aligned_fields(reference, record, hit, &mut Vec::new())?;
    let mut body = Vec::new();
    write_u64(&mut body, fields.start)?;
    write_u8(&mut body, u8::from(fields.reverse))?;
    write_u32(&mut body, fields.events.len() as u32)?;
    for event in fields.events {
        write_u16(&mut body, event.delta)?;
        write_u8(&mut body, event.symbol)?;
        if event.symbol == 4 {
            write_u64(&mut body, event.skip)?;
        }
    }
    let emit_payload = fields.emits.join(&b'\n');
    write_u32(&mut body, emit_payload.len() as u32)?;
    body.extend(emit_payload);
    Ok(body)
}

fn event_delta(position: usize, cursor: usize) -> Result<u16> {
    let delta = position
        .checked_sub(cursor)
        .ok_or_else(|| Error::InvalidInput("noise events are out of order".into()))?;
    u16::try_from(delta)
        .map_err(|_| Error::InvalidInput("noise position delta exceeds u16 capacity".into()))
}

fn substitution_symbol(reference: u8, actual: u8) -> Result<u8> {
    let alternatives = *b"ACGTN";
    let reference_index = alternatives
        .iter()
        .position(|&base| base == reference)
        .ok_or_else(|| Error::InvalidInput("noncanonical reference base".into()))?;
    let actual_index = alternatives
        .iter()
        .position(|&base| base == actual)
        .ok_or_else(|| Error::InvalidInput("noncanonical aligned read base".into()))?;
    if reference_index == actual_index {
        return Err(Error::InvalidInput(
            "substitution symbol matches reference".into(),
        ));
    }
    Ok((0..alternatives.len())
        .filter(|&index| index != reference_index)
        .position(|index| index == actual_index)
        .unwrap() as u8)
}

fn substitution_base(reference: u8, symbol: u8) -> Result<u8> {
    let alternatives = *b"ACGTN";
    if !alternatives.contains(&reference) {
        return Err(Error::InvalidFormat("noncanonical reference base".into()));
    }
    alternatives
        .iter()
        .copied()
        .filter(|&base| base != reference)
        .nth(symbol as usize)
        .ok_or_else(|| Error::InvalidFormat("invalid substitution symbol".into()))
}

fn encode_record_columns<R: ReferenceBases + ?Sized>(
    programs: &mut Vec<u8>,
    literals: &mut Vec<u8>,
    sequence: &[u8],
    mapped: Option<(&GenomeAlignment, &R)>,
    runs: &[ArchiveRun],
) -> Result<()> {
    let mut oriented = sequence.to_vec();
    let (mut ref_pos, reverse) = if let Some((hit, _)) = mapped {
        (
            nonnegative_usize(hit.start, "alignment start")?,
            hit.is_reverse,
        )
    } else {
        (0, false)
    };
    if reverse {
        reverse_complement(&mut oriented);
    }

    write_varint(programs, runs.len() as u64)?;
    for chunk in runs.chunks(4) {
        let mut packed = 0;
        for (index, run) in chunk.iter().enumerate() {
            packed |= (run.op as u8) << (6 - index * 2);
        }
        programs.push(packed);
    }
    let mut query_pos = 0usize;
    for run in runs {
        if run.len == 0 {
            return Err(Error::InvalidInput("invalid archive operation run".into()));
        }
        write_varint(programs, run.len as u64)?;
        match run.op {
            ArchiveOp::Match => {
                let (_, reference) = mapped.ok_or_else(|| {
                    Error::InvalidInput("unmapped body contains a match operation".into())
                })?;
                let ref_end = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                let query_end = query_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if ref_end > reference.len() || query_end > oriented.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds bounds".into(),
                    ));
                }
                ref_pos = ref_end;
                query_pos = query_end;
            }
            ArchiveOp::Literal => {
                query_pos = query_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if query_pos > oriented.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds query".into(),
                    ));
                }
            }
            ArchiveOp::RefAdvance => {
                let (_, reference) = mapped.ok_or_else(|| {
                    Error::InvalidInput("unmapped body contains a reference advance".into())
                })?;
                ref_pos = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds reference".into(),
                    ));
                }
            }
        }
    }
    if query_pos != oriented.len() {
        return Err(Error::InvalidFormat(
            "archive operations do not consume reduced query".into(),
        ));
    }

    let (mut ref_pos, mut query_pos) = if let Some((hit, _)) = mapped {
        (nonnegative_usize(hit.start, "alignment start")?, 0)
    } else {
        (0, 0)
    };
    for run in runs {
        match run.op {
            ArchiveOp::Match => {
                let (_, reference) = mapped.unwrap();
                let end = ref_pos + run.len;
                let mismatches = (0..run.len)
                    .filter_map(|offset| {
                        (reference.base(ref_pos + offset) != oriented[query_pos + offset])
                            .then_some((offset, oriented[query_pos + offset]))
                    })
                    .collect::<Vec<_>>();
                write_varint(programs, mismatches.len() as u64)?;
                let mut previous = 0usize;
                for (position, base) in mismatches {
                    write_varint(programs, (position - previous) as u64)?;
                    programs.push(base);
                    previous = position;
                }
                ref_pos = end;
                query_pos += run.len;
            }
            ArchiveOp::Literal => query_pos += run.len,
            ArchiveOp::RefAdvance => ref_pos += run.len,
        }
    }

    let mut query_pos = 0usize;
    for run in runs {
        match run.op {
            ArchiveOp::Match => {
                query_pos += run.len;
            }
            ArchiveOp::Literal => {
                let end = query_pos + run.len;
                literals.extend_from_slice(&oriented[query_pos..end]);
                query_pos = end;
            }
            ArchiveOp::RefAdvance => {}
        }
    }
    Ok(())
}

fn alignment_candidate_core_cost<R: ReferenceBases + ?Sized>(
    sequence: &[u8],
    mapped: Option<(&GenomeAlignment, &R)>,
    runs: &[ArchiveRun],
) -> Result<usize> {
    let mut programs = Vec::new();
    let mut literals = Vec::new();
    encode_record_columns(&mut programs, &mut literals, sequence, mapped, runs)?;
    Ok(programs.len() + literals.len())
}

fn alignment_candidate_cost<R: ReferenceBases + ?Sized>(
    record: &ReadRecord,
    alignment: Option<&GenomeAlignment>,
    reference: &R,
) -> Result<usize> {
    let runs = if let Some(alignment) = alignment {
        archive_normalize(effective_cigar(alignment, record.sequence.len())?)?
    } else {
        vec![ArchiveRun {
            op: ArchiveOp::Literal,
            len: record.sequence.len(),
        }]
    };
    alignment_candidate_core_cost(
        &record.sequence,
        alignment.map(|alignment| (alignment, reference)),
        &runs,
    )
}

#[cfg(test)]
fn select_alignment<R: ReferenceBases + ?Sized>(
    record: &ReadRecord,
    hits: Vec<minibwa::Hit>,
    reference: &R,
) -> Result<SelectedAlignment> {
    aligners::select_minibwa_alignment(record, hits, reference, &alignment_candidate_cost)
}

#[allow(clippy::too_many_arguments)]
fn decode_v21_records(
    input: &[u8],
    count: usize,
    reference: &[u8],
    layout: AlignmentLayout,
    descriptors: &[IndexDescriptor],
    preserve_quality: bool,
    index_stream_count: usize,
    absolute_position_width: usize,
    delta_position_width: usize,
) -> Result<DecodedRecordSections> {
    let mut cursor = Cursor::new(input);
    let section_count = read_u32(&mut cursor)? as usize;
    let expected_sections = layout
        .streams()
        .checked_mul(V21_MATE_SECTION_COUNT)
        .and_then(|value| value.checked_add(1))
        .and_then(|value| value.checked_add(descriptors.len()))
        .and_then(|value| value.checked_add(index_stream_count.saturating_mul(3)))
        .ok_or_else(|| Error::InvalidFormat("record section count overflow".into()))?;
    if section_count != expected_sections {
        return Err(Error::InvalidFormat("record section count mismatch".into()));
    }
    let mut sections = Vec::with_capacity(section_count);
    for _ in 0..section_count {
        let compressed_len = usize_from(read_u64(&mut cursor)?, "record subsection")?;
        let raw_len = usize_from(read_u64(&mut cursor)?, "record subsection")?;
        let compressed_start = cursor.position() as usize;
        let compressed_end = compressed_start
            .checked_add(compressed_len)
            .ok_or_else(|| Error::InvalidFormat("record subsection overflow".into()))?;
        if compressed_end > input.len() {
            return Err(Error::InvalidFormat("record subsection out of bounds".into()));
        }
        cursor.set_position(compressed_end as u64);
        let mut raw = vec![0u8; raw_len];
        decompress_frame_into(&input[compressed_start..compressed_end], &mut raw)?;
        sections.push(raw);
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat(
            "trailing record section directory bytes".into(),
        ));
    }

    let flags = &sections[0];
    if flags.len() != count
        || flags.iter().any(|&flag| {
            flag > 3 || (layout == AlignmentLayout::SingleEnd && flag != 0 && flag != 3)
        })
    {
        return Err(Error::InvalidFormat(
            "invalid alignment flag section".into(),
        ));
    }
    let mut all = Vec::with_capacity(layout.streams());
    let mut first_alignment = None;
    for mate in 0..layout.streams() {
        let start = 1 + mate * V21_MATE_SECTION_COUNT;
        let (sequences, alignments) = decode_v21_mate_sections(
            &sections[start..start + V21_MATE_SECTION_COUNT],
            flags,
            mate,
            reference,
            first_alignment.as_ref(),
            absolute_position_width,
            delta_position_width,
        )?;
        if mate == 0 {
            first_alignment = Some(alignments);
        }
        all.push(sequences);
    }

    let mut index_sections = Vec::with_capacity(descriptors.len());
    for (descriptor_index, descriptor) in descriptors.iter().enumerate() {
        let section = &sections[1 + layout.streams() * V21_MATE_SECTION_COUNT + descriptor_index];
        let expected_len = descriptor
            .range
            .len()
            .checked_mul(count)
            .ok_or_else(|| Error::InvalidFormat("index subsection length overflow".into()))?;
        if section.len() != expected_len {
            return Err(Error::InvalidFormat(
                "index subsection length mismatch".into(),
            ));
        }
        let indexes = section
            .chunks_exact(descriptor.range.len())
            .map(<[u8]>::to_vec)
            .collect::<Vec<_>>();
        index_sections.push(indexes);
    }
    for (mate, sequences) in all.iter_mut().enumerate() {
        let mut selected = descriptors
            .iter()
            .enumerate()
            .filter(|(_, descriptor)| descriptor.range.mate.stream() == mate)
            .collect::<Vec<_>>();
        selected.sort_by_key(|(_, descriptor)| descriptor.range.start);
        for record in 0..count {
            let rebuilt = {
                let biological = &sequences[record];
                let mut rebuilt = Vec::with_capacity(biological.len());
                let mut biological_start = 0usize;
                let mut removed = 0usize;
                for (descriptor_index, descriptor) in &selected {
                    let boundary =
                        descriptor.range.start.checked_sub(removed).ok_or_else(|| {
                            Error::InvalidFormat("index insertion underflow".into())
                        })?;
                    if boundary < biological_start || boundary > biological.len() {
                        return Err(Error::InvalidFormat("index insertion exceeds read".into()));
                    }
                    let index = &index_sections[*descriptor_index][record];
                    if index.len() != descriptor.range.len() {
                        return Err(Error::InvalidFormat("index read length mismatch".into()));
                    }
                    rebuilt.extend_from_slice(&biological[biological_start..boundary]);
                    rebuilt.extend_from_slice(index);
                    biological_start = boundary;
                    removed = removed
                        .checked_add(descriptor.range.len())
                        .ok_or_else(|| Error::InvalidFormat("index length overflow".into()))?;
                }
                rebuilt.extend_from_slice(&biological[biological_start..]);
                rebuilt
            };
            sequences[record] = rebuilt;
        }
    }
    let mut index_streams = Vec::with_capacity(index_stream_count);
    for stream in 0..index_stream_count {
        let start = 1 + layout.streams() * V21_MATE_SECTION_COUNT + descriptors.len() + stream * 3;
        let quality_block = &sections[start];
        let lengths_block = &sections[start + 1];
        let sequence_block = &sections[start + 2];
        let qualities = if preserve_quality {
            decode_stream_bytes(quality_block, count)?
                .into_iter()
                .map(Some)
                .collect()
        } else {
            if !quality_block.is_empty() {
                return Err(Error::InvalidFormat(
                    "unexpected index quality bytes".into(),
                ));
            }
            vec![None; count]
        };
        if lengths_block.len()
            != count
                .checked_mul(2)
                .ok_or_else(|| Error::InvalidFormat("index read length section overflow".into()))?
        {
            return Err(Error::InvalidFormat(
                "index read length section mismatch".into(),
            ));
        }
        let mut lengths_cursor = Cursor::new(lengths_block.clone());
        let mut offset = 0usize;
        let mut sequences = Vec::with_capacity(count);
        for _ in 0..count {
            let length = read_u16(&mut lengths_cursor)? as usize;
            if length == 0 {
                return Err(Error::InvalidFormat(
                    "index read sequence cannot be empty".into(),
                ));
            }
            let end = offset
                .checked_add(length)
                .ok_or_else(|| Error::InvalidFormat("index read length overflow".into()))?;
            if end > sequence_block.len() {
                return Err(Error::InvalidFormat(
                    "index sequence section mismatch".into(),
                ));
            }
            sequences.push(sequence_block[offset..end].to_vec());
            offset = end;
        }
        if offset != sequence_block.len() {
            return Err(Error::InvalidFormat("trailing index read bytes".into()));
        }
        index_streams.push(DecodedIndexStream {
            qualities,
            sequences,
        });
    }
    Ok(DecodedRecordSections {
        biological: all,
        index_streams,
    })
}

type MateAlignments = Vec<Option<(usize, bool)>>;

fn decode_v21_mate_sections(
    sections: &[Vec<u8>],
    flags: &[u8],
    mate: usize,
    reference: &[u8],
    first_alignment: Option<&MateAlignments>,
    absolute_position_width: usize,
    delta_position_width: usize,
) -> Result<(Vec<Vec<u8>>, MateAlignments)> {
    if sections.len() != V21_MATE_SECTION_COUNT || sections[0].len() != flags.len() * 2 {
        return Err(Error::InvalidFormat("read length section mismatch".into()));
    }
    let aligned_count = flags
        .iter()
        .filter(|&&flag| matches!((flag, mate), (0, _) | (1, 0) | (2, 1)))
        .count();
    let delta_count = if mate == 1 {
        flags.iter().filter(|&&flag| flag == 0).count()
    } else {
        0
    };
    let expected_position_bytes = (aligned_count - delta_count) * absolute_position_width
        + delta_count * delta_position_width;
    if sections[1].len() != expected_position_bytes || sections[6].len() != aligned_count {
        return Err(Error::InvalidFormat(
            "aligned metadata section mismatch".into(),
        ));
    }
    let noise_payload = validated_newline_rows(&sections[3], aligned_count, true, "noise")?;
    let mut emit_count = 0usize;
    let mut event_count = 0usize;
    let mut skip_count = 0usize;
    // Single pass over the noise payload: three separate filter passes here
    // showed up as measurable decode overhead.
    for row in noise_payload.split(|&byte| byte == b'\n') {
        event_count += row.len();
        for &symbol in row {
            match symbol {
                b'4' => skip_count += 1,
                b'5' => emit_count += 1,
                _ => {}
            }
        }
    }
    if sections[2].len() != event_count * 2 {
        return Err(Error::InvalidFormat(
            "noise position section mismatch".into(),
        ));
    }
    if sections[4].len() != skip_count * 8 {
        return Err(Error::InvalidFormat(
            "reference skip section mismatch".into(),
        ));
    }
    let emit_payload = validated_newline_rows(&sections[5], emit_count, false, "emit")?;
    let mut noise_rows = noise_payload.split(|&byte| byte == b'\n');
    let mut emit_rows = emit_payload.split(|&byte| byte == b'\n');
    let mut lengths = Cursor::new(&sections[0][..]);
    let mut positions = Cursor::new(&sections[1][..]);
    let mut noise_positions = Cursor::new(&sections[2][..]);
    let mut skips = Cursor::new(&sections[4][..]);
    let mut reverses = Cursor::new(&sections[6][..]);
    let mut unaligned = Cursor::new(&sections[7][..]);
    let mut noise_row_index = 0usize;
    let mut sequences = Vec::with_capacity(flags.len());
    let mut alignments = vec![None; flags.len()];
    // Reused across records; the aligned-record loop allocates otherwise.
    let mut events: Vec<V21NoiseEvent> = Vec::new();
    let mut emits: Vec<&[u8]> = Vec::new();
    for (record_index, &flag) in flags.iter().enumerate() {
        let length = read_u16(&mut lengths)? as usize;
        let aligned = matches!((flag, mate), (0, _) | (1, 0) | (2, 1));
        if aligned {
            let encoded_position = read_position(
                &mut positions,
                if mate == 1 && flag == 0 {
                    delta_position_width
                } else {
                    absolute_position_width
                },
            )?;
            let start = if mate == 1 && flag == 0 {
                let first = first_alignment
                    .and_then(|alignments| alignments.get(record_index))
                    .and_then(|alignment| *alignment)
                    .ok_or_else(|| {
                        Error::InvalidFormat("paired R2 alignment is missing R1 alignment".into())
                    })?;
                let delta = zigzag_decode(encoded_position)?;
                let first_start = i64::try_from(first.0).map_err(|_| {
                    Error::InvalidFormat("R1 alignment position exceeds signed range".into())
                })?;
                let second_start = first_start.checked_add(delta).ok_or_else(|| {
                    Error::InvalidFormat("R2 relative alignment position overflow".into())
                })?;
                nonnegative_usize(second_start, "R2 alignment position")?
            } else {
                usize_from(encoded_position, "alignment position")?
            };
            let encoded_reverse = read_u8(&mut reverses)?;
            if encoded_reverse > 1 {
                return Err(Error::InvalidFormat(
                    "invalid reverse-complement flag".into(),
                ));
            }
            let reverse = if mate == 1 && flag == 0 {
                let first_reverse = first_alignment
                    .and_then(|alignments| alignments.get(record_index))
                    .and_then(|alignment| *alignment)
                    .ok_or_else(|| {
                        Error::InvalidFormat(
                            "paired R2 orientation is missing R1 orientation".into(),
                        )
                    })?
                    .1;
                if encoded_reverse == 1 {
                    first_reverse
                } else {
                    !first_reverse
                }
            } else {
                encoded_reverse != 0
            };
            let noise_row = noise_rows
                .next()
                .ok_or_else(|| Error::InvalidFormat("noise row count mismatch".into()))?;
            events.clear();
            emits.clear();
            for &encoded_symbol in noise_row {
                if !(b'0'..=b'5').contains(&encoded_symbol) {
                    return Err(Error::InvalidFormat("unknown noise symbol".into()));
                }
                let symbol = encoded_symbol - b'0';
                let delta = read_u16(&mut noise_positions)?;
                let skip = if symbol == 4 {
                    read_u64(&mut skips)?
                } else {
                    0
                };
                if symbol == 5 {
                    let emit_row = emit_rows
                        .next()
                        .ok_or_else(|| Error::InvalidFormat("emit row count mismatch".into()))?;
                    emits.push(emit_row);
                }
                events.push(V21NoiseEvent {
                    delta,
                    symbol,
                    skip,
                });
            }
            noise_row_index += 1;
            sequences.push(reconstruct_v21_aligned_read(
                reference, length, start, reverse, &events, &emits,
            )?);
            alignments[record_index] = Some((start, reverse));
        } else {
            let mut sequence = vec![0; length];
            unaligned.read_exact(&mut sequence)?;
            sequences.push(sequence);
        }
    }
    if lengths.position() != lengths.get_ref().len() as u64
        || positions.position() != positions.get_ref().len() as u64
        || noise_positions.position() != noise_positions.get_ref().len() as u64
        || skips.position() != skips.get_ref().len() as u64
        || reverses.position() != reverses.get_ref().len() as u64
        || unaligned.position() != unaligned.get_ref().len() as u64
        || noise_row_index != aligned_count
    {
        return Err(Error::InvalidFormat(
            "trailing read subsection bytes".into(),
        ));
    }
    Ok((sequences, alignments))
}

/// Validates that `input` holds exactly `expected` newline-separated rows and
/// returns the payload between the outer newlines. Callers iterate rows
/// lazily with `split` instead of materializing a row vector per section.
fn validated_newline_rows<'a>(
    input: &'a [u8],
    expected: usize,
    allow_empty: bool,
    name: &str,
) -> Result<&'a [u8]> {
    if expected == 0 {
        return if input.is_empty() {
            Ok(&[])
        } else {
            Err(Error::InvalidFormat(format!("unexpected {name} rows")))
        };
    }
    if input.last() != Some(&b'\n') {
        return Err(Error::InvalidFormat(format!("unterminated {name} rows")));
    }
    let payload = &input[..input.len() - 1];
    let mut rows = 0usize;
    let mut all_nonempty = true;
    for row in payload.split(|&byte| byte == b'\n') {
        rows += 1;
        all_nonempty &= !row.is_empty();
    }
    if rows != expected || (!allow_empty && !all_nonempty) {
        return Err(Error::InvalidFormat(format!("{name} row count mismatch")));
    }
    Ok(payload)
}

#[cfg(test)]
fn decode_v21_aligned_read(reference: &[u8], length: usize, input: &[u8]) -> Result<Vec<u8>> {
    let mut cursor = Cursor::new(input.to_vec());
    let start = usize_from(read_u64(&mut cursor)?, "alignment position")?;
    let reverse = read_u8(&mut cursor)?;
    if reverse > 1 {
        return Err(Error::InvalidFormat(
            "invalid reverse-complement flag".into(),
        ));
    }
    let event_count = usize_from(read_u32(&mut cursor)? as u64, "noise event count")?;
    let mut events = Vec::with_capacity(event_count);
    let mut emit_count = 0usize;
    for _ in 0..event_count {
        let delta = read_u16(&mut cursor)?;
        let symbol = read_u8(&mut cursor)?;
        let skip = if symbol == 4 {
            read_u64(&mut cursor)?
        } else {
            0
        };
        if symbol == 5 {
            emit_count = emit_count
                .checked_add(1)
                .ok_or_else(|| Error::InvalidFormat("emit count overflow".into()))?;
        } else if symbol > 4 {
            return Err(Error::InvalidFormat("unknown noise symbol".into()));
        }
        events.push(V21NoiseEvent {
            delta,
            symbol,
            skip,
        });
    }
    let emit_len = usize_from(read_u32(&mut cursor)? as u64, "emit payload")?;
    let mut emit_payload = vec![0; emit_len];
    cursor.read_exact(&mut emit_payload)?;
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat(
            "trailing aligned read body bytes".into(),
        ));
    }
    let emits = if emit_count == 0 {
        if !emit_payload.is_empty() {
            return Err(Error::InvalidFormat("unexpected emit payload".into()));
        }
        Vec::new()
    } else {
        let values = emit_payload
            .split(|&byte| byte == b'\n')
            .collect::<Vec<_>>();
        if values.len() != emit_count || values.iter().any(|value| value.is_empty()) {
            return Err(Error::InvalidFormat(
                "emit count does not match noise symbols".into(),
            ));
        }
        values
    };
    reconstruct_v21_aligned_read(
        reference,
        length,
        start,
        reverse != 0,
        &events,
        &emits,
    )
}

fn reconstruct_v21_aligned_read(
    reference: &[u8],
    length: usize,
    start: usize,
    reverse: bool,
    events: &[V21NoiseEvent],
    emits: &[&[u8]],
) -> Result<Vec<u8>> {
    // Fast path: without noise events the read is a contiguous reference
    // slice (the general path below writes every byte exactly once, so the
    // zero-initialized buffer it replaced was pure overhead).
    if events.is_empty() {
        if !emits.is_empty() {
            return Err(Error::InvalidFormat(
                "aligned read length mismatch".into(),
            ));
        }
        let end = start
            .checked_add(length)
            .ok_or_else(|| Error::InvalidFormat("final reference span overflow".into()))?;
        if end > reference.len() {
            return Err(Error::InvalidFormat(
                "final reference span exceeds reference".into(),
            ));
        }
        let mut sequence = reference[start..end].to_vec();
        if reverse {
            reverse_complement(&mut sequence);
        }
        return Ok(sequence);
    }
    // The construction is strictly left-to-right: reference gaps, event
    // effects, and the tail all extend the sequence monotonically, so the
    // sequence is built by extending a capacity-reserved buffer instead of
    // writing into a zero-initialized one.
    let mut sequence = Vec::with_capacity(length);
    let mut ref_pos = start;
    let mut emit_index = 0usize;
    for event in events {
        let reference_end = ref_pos
            .checked_add(event.delta as usize)
            .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
        let target = sequence
            .len()
            .checked_add(event.delta as usize)
            .ok_or_else(|| Error::InvalidFormat("noise position overflow".into()))?;
        if target > length || reference_end > reference.len() {
            return Err(Error::InvalidFormat(format!(
                "aligned read event exceeds bounds: query={} delta={} target={target} length={length} reference={ref_pos} reference_target={reference_end} reference_length={} symbol={}",
                sequence.len(),
                event.delta,
                reference.len(),
                event.symbol
            )));
        }
        sequence.extend_from_slice(&reference[ref_pos..reference_end]);
        ref_pos = reference_end;
        match event.symbol {
            0..=3 => {
                if sequence.len() >= length || ref_pos >= reference.len() {
                    return Err(Error::InvalidFormat("substitution exceeds bounds".into()));
                }
                sequence.push(substitution_base(reference[ref_pos], event.symbol)?);
                ref_pos += 1;
            }
            4 => {
                ref_pos = ref_pos
                    .checked_add(event.skip as usize)
                    .ok_or_else(|| Error::InvalidFormat("reference skip overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "reference skip exceeds reference".into(),
                    ));
                }
            }
            5 => {
                let emit = emits[emit_index];
                emit_index += 1;
                let end = sequence
                    .len()
                    .checked_add(emit.len())
                    .ok_or_else(|| Error::InvalidFormat("emit sequence overflow".into()))?;
                if end > length {
                    return Err(Error::InvalidFormat(
                        "emit sequence exceeds read length".into(),
                    ));
                }
                sequence.extend_from_slice(emit);
            }
            _ => unreachable!(),
        }
    }
    if sequence.len() < length {
        let end = ref_pos
            .checked_add(length - sequence.len())
            .ok_or_else(|| Error::InvalidFormat("final reference span overflow".into()))?;
        if end > reference.len() {
            return Err(Error::InvalidFormat(
                "final reference span exceeds reference".into(),
            ));
        }
        sequence.extend_from_slice(&reference[ref_pos..end]);
    }
    if sequence.len() != length || emit_index != emits.len() {
        return Err(Error::InvalidFormat("aligned read length mismatch".into()));
    }
    if reverse {
        reverse_complement(&mut sequence);
    }
    Ok(sequence)
}

fn decode_ids(
    input: &[u8],
    count: usize,
    mode: IdMode,
    layout: AlignmentLayout,
) -> Result<Vec<Vec<Vec<u8>>>> {
    let mut cursor = Cursor::new(input);
    let mut blocks = vec![Vec::new(); layout.streams()];
    for block in &mut blocks {
        let length = usize_from(read_u64(&mut cursor)?, "ID block")?;
        *block = vec![0; length];
        cursor.read_exact(block)?;
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat("trailing ID block bytes".into()));
    }
    let first = if mode == IdMode::None {
        if blocks.iter().any(|block| !block.is_empty()) {
            return Err(Error::InvalidFormat("unexpected ID bytes".into()));
        }
        vec![Vec::new(); count]
    } else {
        decode_prefix_ids(&blocks[0], count)?
    };
    if layout == AlignmentLayout::SingleEnd {
        if !matches!(mode, IdMode::None | IdMode::Independent) {
            return Err(Error::InvalidFormat(
                "mate-specific ID mode in single-end archive".into(),
            ));
        }
        return Ok(vec![first]);
    }
    if mode == IdMode::SlashMate && !blocks[1].is_empty() {
        return Err(Error::InvalidFormat("unexpected mate-2 ID bytes".into()));
    }
    let (first, second) = match mode {
        IdMode::None => (first, vec![Vec::new(); count]),
        IdMode::Independent => (first, decode_prefix_ids(&blocks[1], count)?),
        IdMode::IdenticalMates => {
            if !blocks[1].is_empty() {
                return Err(Error::InvalidFormat("unexpected mate-2 ID bytes".into()));
            }
            let second = first.clone();
            (first, second)
        }
        // Each identifier is allocated exactly once per mate: mate 1 gets the
        // "/1" suffix appended in place, mate 2 reuses a clone of mate 1 with
        // the suffix replaced (no reallocation, same length).
        IdMode::SlashMate => {
            let first = first
                .into_iter()
                .map(|mut id| {
                    id.extend_from_slice(b"/1");
                    id
                })
                .collect::<Vec<_>>();
            let second = first
                .iter()
                .map(|id| {
                    let mut value = id.clone();
                    let base = value.len() - 2;
                    value.truncate(base);
                    value.extend_from_slice(b"/2");
                    value
                })
                .collect();
            (first, second)
        }
    };
    Ok(vec![first, second])
}

fn decode_prefix_ids(input: &[u8], count: usize) -> Result<Vec<Vec<u8>>> {
    let mut cursor = Cursor::new(input);
    let mut result: Vec<Vec<u8>> = Vec::with_capacity(count);
    for _ in 0..count {
        let common = usize_from(read_varint(&mut cursor)?, "ID prefix")?;
        let suffix_len = usize_from(read_varint(&mut cursor)?, "ID suffix")?;
        if common > result.last().map_or(0, |value| value.len()) {
            return Err(Error::InvalidFormat("ID prefix out of bounds".into()));
        }
        let empty: Vec<u8> = Vec::new();
        let previous = result.last().unwrap_or(&empty);
        let mut value = previous[..common].to_vec();
        let start = cursor.position() as usize;
        let end = start
            .checked_add(suffix_len)
            .ok_or_else(|| Error::InvalidFormat("ID suffix overflow".into()))?;
        if end > cursor.get_ref().len() {
            return Err(Error::InvalidFormat("ID suffix out of bounds".into()));
        }
        value.extend_from_slice(&cursor.get_ref()[start..end]);
        cursor.set_position(end as u64);
        result.push(value);
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat("trailing ID bytes".into()));
    }
    Ok(result)
}

fn decode_qualities(
    input: &[u8],
    count: usize,
    preserve: bool,
    layout: AlignmentLayout,
) -> Result<Vec<Vec<Option<Vec<u8>>>>> {
    let mut cursor = Cursor::new(input);
    let mut result = Vec::with_capacity(layout.streams());
    for _ in 0..layout.streams() {
        let length = usize_from(read_u64(&mut cursor)?, "quality stream")?;
        let mut stream = vec![0; length];
        cursor.read_exact(&mut stream)?;
        if preserve {
            result.push(
                decode_stream_bytes(&stream, count)?
                    .into_iter()
                    .map(Some)
                    .collect(),
            );
        } else if length != 0 {
            return Err(Error::InvalidFormat("unexpected quality bytes".into()));
        } else {
            result.push(vec![None; count]);
        }
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat(
            "trailing quality section bytes".into(),
        ));
    }
    Ok(result)
}

fn archive_normalize(cigar: impl IntoIterator<Item = CigarOp>) -> Result<Vec<ArchiveRun>> {
    let mut runs: Vec<ArchiveRun> = Vec::new();
    for cigar_op in cigar {
        if cigar_op.len == 0 || matches!(cigar_op.kind, CigarKind::HardClip | CigarKind::Pad) {
            continue;
        }
        let op = match cigar_op.kind {
            CigarKind::Match | CigarKind::Eq | CigarKind::Diff => ArchiveOp::Match,
            CigarKind::Ins | CigarKind::SoftClip => ArchiveOp::Literal,
            CigarKind::Del | CigarKind::RefSkip => ArchiveOp::RefAdvance,
            CigarKind::HardClip | CigarKind::Pad => unreachable!(),
        };
        if let Some(previous) = runs.last_mut().filter(|previous| previous.op == op) {
            previous.len = previous
                .len
                .checked_add(cigar_op.len as usize)
                .ok_or_else(|| Error::InvalidInput("archive operation length overflow".into()))?;
        } else {
            runs.push(ArchiveRun {
                op,
                len: cigar_op.len as usize,
            });
        }
    }
    Ok(runs)
}

fn encode_index_metadata(descriptors: &[IndexDescriptor]) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    write_u64(&mut output, descriptors.len() as u64)?;
    for descriptor in descriptors {
        write_u8(
            &mut output,
            match descriptor.range.mate {
                EmbeddedIndexMate::R1 => 0,
                EmbeddedIndexMate::R2 => 1,
            },
        )?;
        write_u64(&mut output, descriptor.range.start as u64)?;
        write_u64(&mut output, descriptor.range.len() as u64)?;
    }
    Ok(output)
}

fn decode_index_metadata(input: &[u8], layout: AlignmentLayout) -> Result<Vec<IndexDescriptor>> {
    if input.len() < 8 {
        return Err(Error::InvalidFormat(
            "embedded index metadata is too short".into(),
        ));
    }
    let mut input = Cursor::new(input.to_vec());
    let count = usize_from(read_u64(&mut input)?, "embedded index range count")?;
    let minimum_descriptors = count
        .checked_mul(17)
        .ok_or_else(|| Error::InvalidFormat("embedded index metadata size overflow".into()))?;
    if minimum_descriptors > input.get_ref().len() - input.position() as usize {
        return Err(Error::InvalidFormat(
            "embedded index range count exceeds metadata".into(),
        ));
    }
    let mut descriptors = Vec::with_capacity(count);
    for _ in 0..count {
        let mate = match read_u8(&mut input)? {
            0 => EmbeddedIndexMate::R1,
            1 => EmbeddedIndexMate::R2,
            _ => return Err(Error::InvalidFormat("unknown embedded index mate".into())),
        };
        let start = usize_from(read_u64(&mut input)?, "embedded index start")?;
        let length = usize_from(read_u64(&mut input)?, "embedded index length")?;
        let end = start
            .checked_add(length)
            .ok_or_else(|| Error::InvalidFormat("embedded index range overflow".into()))?;
        descriptors.push(IndexDescriptor {
            range: EmbeddedIndexRange { mate, start, end },
        });
    }
    if input.position() != input.get_ref().len() as u64 {
        return Err(Error::InvalidFormat(
            "trailing embedded index metadata bytes".into(),
        ));
    }
    let ranges = descriptors
        .iter()
        .map(|descriptor| descriptor.range.clone())
        .collect::<Vec<_>>();
    validate_index_ranges(&ranges, layout).map_err(|error| match error {
        Error::InvalidInput(message) | Error::InvalidFormat(message) => {
            Error::InvalidFormat(message)
        }
        error => error,
    })?;
    Ok(descriptors)
}

fn read_header<R: Read>(input: &mut R) -> Result<ArchiveHeader> {
    let parsed = parse_preamble(input)?;
    Ok(ArchiveHeader {
        flags: parsed.flags,
        id_mode: parsed.id_mode,
        layout: parsed.layout,
        chunk_count: parsed.chunk_count,
        reference: decode_reference(&parsed.reference_block)?,
        indexes: decode_index_metadata(&parsed.index_metadata, parsed.layout)?,
        absolute_position_width: parsed.absolute_position_width,
        delta_position_width: parsed.delta_position_width,
        index_streams: parsed.index_streams,
    })
}

/// The archive preamble, parsed and CRC-verified, with the embedded blocks
/// retained as raw bytes so callers can avoid decoding the reference.
struct ParsedPreamble {
    version_major: u16,
    version_minor: u16,
    flags: u16,
    id_mode: IdMode,
    layout: AlignmentLayout,
    absolute_position_width: usize,
    delta_position_width: usize,
    index_streams: usize,
    chunk_count: u64,
    reference_block: Vec<u8>,
    index_metadata: Vec<u8>,
}

fn parse_preamble<R: Read>(input: &mut R) -> Result<ParsedPreamble> {
    let mut magic = [0; 9];
    input.read_exact(&mut magic)?;
    if magic != *MAGIC {
        return Err(Error::InvalidFormat("bad alignment archive magic".into()));
    }
    let mut version = [0; 4];
    input.read_exact(&mut version)?;
    let mut version_fields = Cursor::new(&version[..]);
    let version_major = read_u16(&mut version_fields)?;
    let version_minor = read_u16(&mut version_fields)?;
    if version_major != VERSION_MAJOR || version_minor != VERSION_MINOR {
        return Err(Error::InvalidFormat(format!(
            "unsupported alignment archive version {version_major}.{version_minor}; expected \
             {VERSION_MAJOR}.{VERSION_MINOR}"
        )));
    }
    let mut fields = [0u8; PREAMBLE_BYTES - 9 - 4];
    input.read_exact(&mut fields)?;
    let mut parsed = Cursor::new(&fields[..]);
    let flags = read_u16(&mut parsed)?;
    if flags & !0x1f != 0 {
        return Err(Error::InvalidFormat(
            "unknown alignment archive flags".into(),
        ));
    }
    let index_stream_count = usize::from(read_u8(&mut parsed)?);
    if index_stream_count > 2 {
        return Err(Error::InvalidFormat(
            "index stream count must be at most two".into(),
        ));
    }
    let absolute_position_width = usize::from(read_u8(&mut parsed)?);
    let delta_position_width = usize::from(read_u8(&mut parsed)?);
    if !matches!(absolute_position_width, 1 | 2 | 3 | 4 | 8)
        || !matches!(delta_position_width, 1 | 2 | 3 | 4 | 8)
    {
        return Err(Error::InvalidFormat("invalid position width".into()));
    }
    let id_mode = IdMode::from_flags(flags)?;
    let layout = if flags & 4 != 0 {
        AlignmentLayout::SingleEnd
    } else {
        AlignmentLayout::PairedEnd
    };
    if (flags & 1 != 0) != (id_mode != IdMode::None) {
        return Err(Error::InvalidFormat(
            "identifier flag and mode disagree".into(),
        ));
    }
    if layout == AlignmentLayout::SingleEnd
        && !matches!(id_mode, IdMode::None | IdMode::Independent)
    {
        return Err(Error::InvalidFormat(
            "mate-specific ID mode in single-end archive".into(),
        ));
    }
    let reference_len = usize_from(read_u64(&mut parsed)?, "reference block")?;
    let index_metadata_len = usize_from(read_u64(&mut parsed)?, "embedded index metadata block")?;
    let chunks = read_u64(&mut parsed)?;
    let stored_crc = read_u32(input)?;
    let reference_block = read_incremental(input, reference_len, "reference block")?;
    let index_metadata =
        read_incremental(input, index_metadata_len, "embedded index metadata block")?;
    let mut preamble = [0u8; PREAMBLE_BYTES];
    preamble[..magic.len()].copy_from_slice(&magic);
    preamble[magic.len()..magic.len() + version.len()].copy_from_slice(&version);
    preamble[magic.len() + version.len()..].copy_from_slice(&fields);
    let computed_crc = preamble_crc32c(&preamble, &reference_block, &index_metadata);
    if computed_crc != stored_crc {
        return Err(Error::Checksum(format!(
            "archive preamble: stored CRC 0x{stored_crc:08x} does not match computed CRC \
             0x{computed_crc:08x}"
        )));
    }
    Ok(ParsedPreamble {
        version_major,
        version_minor,
        flags,
        id_mode,
        layout,
        absolute_position_width,
        delta_position_width,
        index_streams: index_stream_count,
        chunk_count: chunks,
        reference_block,
        index_metadata,
    })
}

/// CRC32C over the archive preamble: the fixed header bytes followed by the
/// embedded reference block and the index metadata block.
fn preamble_crc32c(preamble: &[u8], reference_block: &[u8], index_metadata: &[u8]) -> u32 {
    let crc = crc32c::crc32c(preamble);
    let crc = crc32c::crc32c_append(crc, reference_block);
    crc32c::crc32c_append(crc, index_metadata)
}

/// Recompute the preamble CRC from the bytes on disk and patch it into the
/// header. Used after the preamble is first written and again after the final
/// chunk count is patched, since the chunk count is covered by the CRC.
fn write_preamble_crc(archive: &mut File) -> Result<()> {
    archive.seek(SeekFrom::Start(0))?;
    let mut preamble = [0u8; PREAMBLE_BYTES];
    archive.read_exact(&mut preamble)?;
    let reference_len = u64::from_le_bytes(preamble[18..26].try_into().unwrap());
    let metadata_len = u64::from_le_bytes(preamble[26..34].try_into().unwrap());
    let mut crc = crc32c::crc32c(&preamble);
    // The CRC field itself sits between the covered header bytes and the
    // reference block, so resume hashing after it.
    archive.seek(SeekFrom::Start(PREAMBLE_BYTES as u64 + 4))?;
    let mut remaining = reference_len + metadata_len;
    let mut buffer = vec![0u8; 64 * 1024];
    while remaining > 0 {
        let want = remaining.min(buffer.len() as u64) as usize;
        archive.read_exact(&mut buffer[..want])?;
        crc = crc32c::crc32c_append(crc, &buffer[..want]);
        remaining -= want as u64;
    }
    archive.seek(SeekFrom::Start(PREAMBLE_BYTES as u64))?;
    write_u32(archive, crc)?;
    // Restore the position to the first byte after the preamble (header
    // fields, CRC field, reference block, and metadata) so callers that are
    // mid-write continue appending in the right place.
    archive.seek(SeekFrom::Start(
        PREAMBLE_BYTES as u64 + 4 + reference_len + metadata_len,
    ))?;
    Ok(())
}

fn encode_reference(
    reference: &[u8],
    options: &crate::AlignmentCompressionOptions,
) -> Result<EncodedReference> {
    let mut raw = Vec::with_capacity(8 + reference.len().div_ceil(4));
    write_u64(&mut raw, reference.len() as u64)?;
    for (index, &base) in reference.iter().enumerate() {
        let bits = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                return Err(Error::InvalidInput(
                    "reference contains a noncanonical base".into(),
                ));
            }
        };
        if index % 4 == 0 {
            raw.push(0);
        }
        *raw.last_mut().unwrap() |= bits << (6 - (index % 4) * 2);
    }
    Ok(EncodedReference {
        block: zstd::bulk::compress(&raw, options.zstd_level)?,
        raw_len: raw.len() as u64,
    })
}

fn decode_reference(block: &[u8]) -> Result<Vec<u8>> {
    let frame_len = zstd::zstd_safe::find_frame_compressed_size(block)
        .map_err(|_| Error::InvalidFormat("invalid reference zstd frame".into()))?;
    if frame_len != block.len() {
        return Err(Error::InvalidFormat(
            "trailing reference compressed bytes".into(),
        ));
    }

    let invalid_payload = || Error::InvalidFormat("invalid reference payload".into());
    // Decompress the packed payload exactly once, sizing the buffer from the
    // frame's embedded content size when present; legacy frames without one
    // fall back to the streaming decoder.
    let payload = match zstd::zstd_safe::get_frame_content_size(block)
        .map_err(|_| invalid_payload())?
    {
        Some(size) => {
            let mut payload =
                vec![0u8; usize::try_from(size).map_err(|_| invalid_payload())?];
            decompress_frame_into(block, &mut payload)?;
            payload
        }
        None => zstd::stream::decode_all(Cursor::new(block)).map_err(|_| invalid_payload())?,
    };
    let count_bytes: [u8; 8] = payload
        .get(..8)
        .and_then(|bytes| bytes.try_into().ok())
        .ok_or_else(invalid_payload)?;
    let base_count = usize_from(u64::from_le_bytes(count_bytes), "reference base count")?;
    let packed_len = base_count.div_ceil(4);
    if payload.len() != 8 + packed_len {
        return Err(Error::InvalidFormat(
            "reference packed payload length mismatch".into(),
        ));
    }
    let last_packed = payload[payload.len() - 1];
    let remainder = base_count % 4;
    if remainder != 0 && last_packed & ((1 << ((4 - remainder) * 2)) - 1) != 0 {
        return Err(Error::InvalidFormat(
            "nonzero reference padding bits".into(),
        ));
    }

    let mut reference = Vec::with_capacity(base_count);
    for &packed in &payload[8..] {
        for shift in [6, 4, 2, 0] {
            if reference.len() == base_count {
                break;
            }
            reference.push(b"ACGT"[((packed >> shift) & 3) as usize]);
        }
    }
    Ok(reference)
}

#[cfg(test)]
fn encode_packed(sequence: &[u8]) -> Result<Vec<u8>> {
    let packed = PackedSequence::try_from_bytes(sequence)?;
    let mut output = Vec::new();
    write_u32(&mut output, packed.len())?;
    if packed.is_packed() {
        write_u8(&mut output, 1)?;
        write_u64(&mut output, packed.encoded_size() as u64)?;
        output.extend(packed.encoded_bytes());
    } else {
        let sparse = encode_sparse_body(sequence)?;
        let (mode, body) = if sparse.len() < sequence.len() {
            (2, sparse.as_slice())
        } else {
            (0, sequence)
        };
        write_u8(&mut output, mode)?;
        write_u64(&mut output, body.len() as u64)?;
        output.extend(body);
    }
    Ok(output)
}
#[cfg(test)]
fn read_packed(input: &mut Cursor<Vec<u8>>) -> Result<PackedSequence> {
    let length = read_u32(input)?;
    let mode = read_u8(input)?;
    let body_len = usize_from(read_u64(input)?, "packed sequence body")?;
    let mut body = vec![0; body_len];
    input.read_exact(&mut body)?;
    match mode {
        0 => PackedSequence::from_encoded(length, false, body),
        1 => PackedSequence::from_encoded(length, true, body),
        2 => PackedSequence::try_from_bytes(&decode_sparse_body(length as usize, &body)?),
        _ => Err(Error::InvalidFormat(
            "unknown packed sequence encoding".into(),
        )),
    }
}

#[cfg(test)]
fn encode_sparse_body(sequence: &[u8]) -> Result<Vec<u8>> {
    let mut packed = vec![0; sequence.len().div_ceil(4)];
    let mut exceptions = Vec::new();
    for (index, &base) in sequence.iter().enumerate() {
        let bits = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                exceptions.push((index, base));
                0
            }
        };
        packed[index / 4] |= bits << (6 - (index % 4) * 2);
    }
    write_varint(&mut packed, exceptions.len() as u64)?;
    let mut previous = 0usize;
    for (position, base) in exceptions {
        write_varint(&mut packed, (position - previous) as u64)?;
        packed.push(base);
        previous = position;
    }
    Ok(packed)
}

#[cfg(test)]
fn decode_sparse_body(length: usize, body: &[u8]) -> Result<Vec<u8>> {
    let packed_len = length.div_ceil(4);
    if body.len() < packed_len {
        return Err(Error::InvalidFormat("truncated sparse sequence".into()));
    }
    let mut sequence = (0..length)
        .map(|index| b"ACGT"[((body[index / 4] >> (6 - (index % 4) * 2)) & 3) as usize])
        .collect::<Vec<_>>();
    let mut cursor = Cursor::new(body[packed_len..].to_vec());
    let count = usize_from(read_varint(&mut cursor)?, "sparse sequence exceptions")?;
    let mut previous = 0usize;
    for index in 0..count {
        let position = previous
            .checked_add(usize_from(
                read_varint(&mut cursor)?,
                "sparse sequence exception offset",
            )?)
            .ok_or_else(|| Error::InvalidFormat("sparse sequence offset overflow".into()))?;
        if position >= length || (index != 0 && position <= previous) {
            return Err(Error::InvalidFormat(
                "sparse sequence exception out of bounds".into(),
            ));
        }
        sequence[position] = read_u8(&mut cursor)?;
        previous = position;
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat(
            "trailing sparse sequence bytes".into(),
        ));
    }
    Ok(sequence)
}
fn cigar_code(op: CigarKind) -> Result<u8> {
    Ok(match op {
        CigarKind::Match => 0,
        CigarKind::Ins => 1,
        CigarKind::Del => 2,
        CigarKind::RefSkip => 3,
        CigarKind::SoftClip => 4,
        CigarKind::HardClip => 5,
        CigarKind::Pad => 6,
        CigarKind::Eq => 7,
        CigarKind::Diff => 8,
    })
}
fn cigar_kind(code: u8) -> Result<CigarKind> {
    Ok(match code {
        0 => CigarKind::Match,
        1 => CigarKind::Ins,
        2 => CigarKind::Del,
        3 => CigarKind::RefSkip,
        4 => CigarKind::SoftClip,
        5 => CigarKind::HardClip,
        6 => CigarKind::Pad,
        7 => CigarKind::Eq,
        8 => CigarKind::Diff,
        _ => return Err(Error::InvalidFormat("unknown CIGAR operation".into())),
    })
}

fn write_u8(output: &mut impl Write, value: u8) -> Result<()> {
    output.write_all(&[value])?;
    Ok(())
}
fn write_u16(output: &mut impl Write, value: u16) -> Result<()> {
    output.write_all(&value.to_le_bytes())?;
    Ok(())
}
fn write_u32(output: &mut impl Write, value: u32) -> Result<()> {
    output.write_all(&value.to_le_bytes())?;
    Ok(())
}
fn write_u64(output: &mut impl Write, value: u64) -> Result<()> {
    output.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_varint(output: &mut Vec<u8>, mut value: u64) -> Result<()> {
    while value >= 0x80 {
        output.push((value as u8) | 0x80);
        value >>= 7;
    }
    output.push(value as u8);
    Ok(())
}

fn read_varint<T: AsRef<[u8]>>(input: &mut Cursor<T>) -> Result<u64> {
    let mut value = 0u64;
    for shift in (0..=63).step_by(7) {
        let byte = read_u8(input)?;
        let bits = (byte & 0x7f) as u64;
        value |= bits << shift;
        if byte & 0x80 == 0 {
            if shift == 63 && bits > 1 {
                return Err(Error::InvalidFormat("varint overflow".into()));
            }
            return Ok(value);
        }
    }
    Err(Error::InvalidFormat("unterminated varint".into()))
}

fn zigzag_encode(value: i64) -> u64 {
    ((value << 1) ^ (value >> 63)) as u64
}

fn zigzag_decode(value: u64) -> Result<i64> {
    let decoded = (value >> 1) as i64;
    if value & 1 == 0 {
        Ok(decoded)
    } else {
        decoded
            .checked_neg()
            .and_then(|value| value.checked_sub(1))
            .ok_or_else(|| Error::InvalidFormat("ZigZag value overflow".into()))
    }
}
fn read_u8(input: &mut impl Read) -> Result<u8> {
    let mut bytes = [0; 1];
    input.read_exact(&mut bytes)?;
    Ok(bytes[0])
}
fn read_u16(input: &mut impl Read) -> Result<u16> {
    let mut bytes = [0; 2];
    input.read_exact(&mut bytes)?;
    Ok(u16::from_le_bytes(bytes))
}
fn read_u32(input: &mut impl Read) -> Result<u32> {
    let mut bytes = [0; 4];
    input.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

/// Read exactly `len` bytes, growing the buffer incrementally so a corrupt or
/// attacker-controlled length cannot trigger a huge allocation up front.
fn read_incremental(input: &mut impl Read, len: usize, label: &str) -> Result<Vec<u8>> {
    let mut data = Vec::with_capacity(len.min(4 * 1024 * 1024));
    let mut buffer = [0u8; 64 * 1024];
    let mut remaining = len;
    while remaining > 0 {
        let want = remaining.min(buffer.len());
        let read = input.read(&mut buffer[..want])?;
        if read == 0 {
            return Err(Error::Io(std::io::Error::new(
                std::io::ErrorKind::UnexpectedEof,
                format!("{label} ended after {} of {} bytes", data.len(), len),
            )));
        }
        data.extend_from_slice(&buffer[..read]);
        remaining -= read;
    }
    Ok(data)
}

fn read_u64(input: &mut impl Read) -> Result<u64> {
    let mut bytes = [0; 8];
    input.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}
fn usize_from(value: u64, what: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| Error::InvalidFormat(format!("{what} is too large")))
}
fn nonnegative_usize(value: i64, what: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| Error::InvalidFormat(format!("invalid {what}")))
}

fn display_duration(duration: std::time::Duration) -> String {
    let seconds = duration.as_secs_f64();
    if seconds >= 60.0 {
        format!("{:.2}m", seconds / 60.0)
    } else {
        format!("{seconds:.2}s")
    }
}

fn mapping_progress() -> ProgressBar {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} mapping reads [{elapsed_precise}] {human_pos} reads @ {per_sec} | {msg}",
        )
        .expect("valid mapping progress template"),
    );
    progress.enable_steady_tick(std::time::Duration::from_millis(250));
    progress
}

fn stage_progress(message: &'static str) -> ProgressBar {
    let progress = ProgressBar::new_spinner();
    progress.set_style(
        ProgressStyle::with_template("{spinner:.green} {msg} [{elapsed_precise}]")
            .expect("valid stage progress template"),
    );
    progress.set_message(message);
    progress.enable_steady_tick(std::time::Duration::from_millis(250));
    progress
}

fn chunk_progress(total: u64) -> ProgressBar {
    let progress = ProgressBar::new(total);
    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} compressing chunks [{elapsed_precise}] {pos}/{len} [{wide_bar:.cyan/blue}] {percent}% ETA {eta}",
        )
        .expect("valid chunk progress template"),
    );
    progress.enable_steady_tick(std::time::Duration::from_millis(250));
    progress
}

fn stage_log(progress: &ProgressBar, message: impl AsRef<str>) {
    if progress.is_hidden() {
        eprintln!("{}", message.as_ref());
    } else {
        progress.println(message.as_ref());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compress_reference_payload(payload: &[u8]) -> Vec<u8> {
        zstd::stream::encode_all(Cursor::new(payload), 1).unwrap()
    }

    fn progress_selection(candidate_hits: u32) -> SelectedAlignment {
        SelectedAlignment {
            hit: None,
            selection: AlignmentSelection {
                candidate_hits,
                rejected_alignment: candidate_hits != 0,
                ..AlignmentSelection::default()
            },
        }
    }

    fn single_quality_pair(quality: &[u8]) -> AlignedPair {
        AlignedPair {
            first: ReadRecord::new("read", vec![b'A'; quality.len()], Some(quality.to_vec())),
            second: None,
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: Vec::new(),
            index_reads: Vec::new(),
        }
    }

    #[test]
    fn mapping_progress_accumulates_candidate_hits_per_mate() {
        let progress = ProgressBar::hidden();
        let mut stats = MappingProgressStats::new(AlignmentLayout::PairedEnd);
        assert_eq!(stats.message(), "R1 aligned 0.00% | R2 aligned 0.00%");

        stats.observe(
            0,
            &[progress_selection(1), progress_selection(0)],
            &progress,
        );
        stats.observe(0, &[progress_selection(2)], &progress);
        stats.observe(
            1,
            &[
                progress_selection(0),
                progress_selection(1),
                progress_selection(1),
                progress_selection(1),
            ],
            &progress,
        );

        assert_eq!(stats.processed, [3, 4]);
        assert_eq!(stats.aligned, [2, 3]);
        assert_eq!(stats.message(), "R1 aligned 66.67% | R2 aligned 75.00%");
        assert_eq!(progress.position(), 7);
    }

    #[test]
    fn single_end_mapping_progress_omits_r2() {
        let progress = ProgressBar::hidden();
        let mut stats = MappingProgressStats::new(AlignmentLayout::SingleEnd);
        stats.observe(0, &[progress_selection(1)], &progress);
        assert_eq!(stats.message(), "R1 aligned 100.00%");
        assert_eq!(stats.aligned[0], 1);
    }

    #[test]
    fn reference_block_round_trips_canonical_lengths_zero_through_nine() {
        let reference = b"ACGTACGTA";
        let expected_packed = [
            &[][..],
            &[0x00],
            &[0x10],
            &[0x18],
            &[0x1b],
            &[0x1b, 0x00],
            &[0x1b, 0x10],
            &[0x1b, 0x18],
            &[0x1b, 0x1b],
            &[0x1b, 0x1b, 0x00],
        ];
        for length in 0..=reference.len() {
            let encoded = encode_reference(
                &reference[..length],
                &crate::AlignmentCompressionOptions::default(),
            )
            .unwrap();
            let payload = zstd::stream::decode_all(Cursor::new(&encoded.block)).unwrap();
            assert_eq!(&payload[..8], &(length as u64).to_le_bytes());
            assert_eq!(&payload[8..], expected_packed[length]);
            assert_eq!(encoded.raw_len, (8 + length.div_ceil(4)) as u64);
            assert_eq!(
                decode_reference(&encoded.block).unwrap(),
                &reference[..length]
            );
        }
    }

    #[test]
    fn reference_encoder_rejects_noncanonical_bases() {
        for reference in [b"N".as_slice(), b"a", b"ACGU", b"ACGT\0"] {
            assert!(matches!(
                encode_reference(reference, &crate::AlignmentCompressionOptions::default()),
                Err(Error::InvalidInput(_))
            ));
        }
    }

    #[test]
    fn compact_reference_rejects_noncanonical_bases() {
        let reference = b"AN";
        let mut mask = ReferenceMask::new(reference.len());
        mask.mark_range(0, reference.len()).unwrap();
        assert!(matches!(
            compact_reference(reference.as_slice(), &mask),
            Err(Error::InvalidInput(_))
        ));
    }

    #[test]
    fn reference_decoder_rejects_malformed_lengths_and_padding() {
        for payload in [
            vec![0; 7],
            u64::MAX.to_le_bytes().to_vec(),
            [5u64.to_le_bytes().as_slice(), &[0x1b]].concat(),
            [4u64.to_le_bytes().as_slice(), &[0x1b, 0x00]].concat(),
            [1u64.to_le_bytes().as_slice(), &[0x01]].concat(),
        ] {
            assert!(decode_reference(&compress_reference_payload(&payload)).is_err());
        }
    }

    #[test]
    fn reference_decoder_rejects_malformed_or_trailing_compressed_data() {
        let encoded = encode_reference(b"ACGT", &crate::AlignmentCompressionOptions::default())
            .unwrap()
            .block;

        assert!(decode_reference(&encoded[..encoded.len() - 1]).is_err());

        let mut trailing_byte = encoded.clone();
        trailing_byte.push(0);
        assert!(decode_reference(&trailing_byte).is_err());

        let mut concatenated = encoded.clone();
        concatenated.extend(compress_reference_payload(&0u64.to_le_bytes()));
        assert!(decode_reference(&concatenated).is_err());
    }

    #[test]
    fn embedded_indexes_split_and_reconstruct_full_quality_rows() {
        let ranges = vec![
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 6,
                end: 8,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R2,
                start: 1,
                end: 4,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 1,
                end: 3,
            },
        ];
        validate_index_ranges(&ranges, AlignmentLayout::PairedEnd).unwrap();
        let first_sequence = b"ACGTACGTAA".to_vec();
        let second_sequence = b"TTGCAACCGG".to_vec();
        let first_quality = b"ABCDEFGHIJ".to_vec();
        let second_quality = b"JKLMNOPQRS".to_vec();
        let mut first =
            ReadRecord::new("first", first_sequence.clone(), Some(first_quality.clone()));
        let mut second = ReadRecord::new(
            "second",
            second_sequence.clone(),
            Some(second_quality.clone()),
        );
        let extracted = split_pair_indexes(&mut first, Some(&mut second), &ranges, true).unwrap();

        assert_eq!(extracted[0].sequence, b"GT");
        assert_eq!(extracted[1].sequence, b"TGC");
        assert_eq!(extracted[2].sequence, b"CG");
        assert_eq!(first.sequence, b"ATACAA");
        assert_eq!(second.sequence, b"TAACCGG");

        let descriptors = ranges
            .iter()
            .cloned()
            .map(|range| IndexDescriptor { range })
            .collect::<Vec<_>>();
        let pair = AlignedPair {
            first,
            second: Some(second),
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: extracted,
            index_reads: Vec::new(),
        };
        assert_eq!(
            reconstruct_quality_row(&pair, 0, &descriptors).unwrap(),
            first_quality
        );
        assert_eq!(
            reconstruct_quality_row(&pair, 1, &descriptors).unwrap(),
            second_quality
        );
    }

    #[test]
    fn embedded_indexes_may_consume_the_entire_read() {
        let ranges = vec![EmbeddedIndexRange {
            mate: EmbeddedIndexMate::R1,
            start: 0,
            end: 10,
        }];
        validate_index_ranges(&ranges, AlignmentLayout::SingleEnd).unwrap();
        let original_sequence = b"ACGTACGTAA".to_vec();
        let original_quality = b"ABCDEFGHIJ".to_vec();
        let mut record = ReadRecord::new(
            "read",
            original_sequence.clone(),
            Some(original_quality.clone()),
        );
        let extracted = split_pair_indexes(&mut record, None, &ranges, true).unwrap();

        assert_eq!(extracted.len(), 1);
        assert_eq!(extracted[0].sequence, original_sequence);
        assert_eq!(
            extracted[0].quality.as_deref(),
            Some(original_quality.as_slice())
        );
        assert!(record.sequence.is_empty());
        assert!(record.quality.as_ref().is_some_and(Vec::is_empty));

        let descriptors = ranges
            .iter()
            .cloned()
            .map(|range| IndexDescriptor { range })
            .collect::<Vec<_>>();
        let pair = AlignedPair {
            first: record,
            second: None,
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: extracted,
            index_reads: Vec::new(),
        };
        assert_eq!(
            reconstruct_quality_row(&pair, 0, &descriptors).unwrap(),
            original_quality
        );
    }

    #[test]
    fn aligned_quality_stream_round_trips_one_fixed_width_raw_run() {
        let pairs = [single_quality_pair(b"AB"), single_quality_pair(b"CD")];
        let options = crate::AlignmentCompressionOptions {
            preserve_quality: true,
            ..crate::AlignmentCompressionOptions::default()
        };
        let mut encoded = Vec::new();
        write_quality_stream(
            &mut encoded,
            &pairs,
            &[],
            &options,
            AlignmentLayout::SingleEnd,
        )
        .unwrap();

        let expected_stream = [
            2u32.to_le_bytes().as_slice(),
            2u32.to_le_bytes().as_slice(),
            b"ABCD",
        ]
        .concat();
        assert_eq!(&encoded[..8], &(expected_stream.len() as u64).to_le_bytes());
        assert_eq!(&encoded[8..], expected_stream);
        assert_eq!(
            decode_qualities(&encoded, pairs.len(), true, AlignmentLayout::SingleEnd).unwrap(),
            vec![vec![Some(b"AB".to_vec()), Some(b"CD".to_vec())]]
        );
    }

    #[test]
    fn aligned_quality_stream_round_trips_variable_width_raw_runs() {
        let pairs = [
            single_quality_pair(b"A"),
            single_quality_pair(b"BC"),
            single_quality_pair(b"DE"),
            single_quality_pair(b"FGH"),
        ];
        let options = crate::AlignmentCompressionOptions {
            preserve_quality: true,
            ..crate::AlignmentCompressionOptions::default()
        };
        let mut encoded = Vec::new();
        write_quality_stream(
            &mut encoded,
            &pairs,
            &[],
            &options,
            AlignmentLayout::SingleEnd,
        )
        .unwrap();

        assert_eq!(
            decode_qualities(&encoded, pairs.len(), true, AlignmentLayout::SingleEnd).unwrap(),
            vec![vec![
                Some(b"A".to_vec()),
                Some(b"BC".to_vec()),
                Some(b"DE".to_vec()),
                Some(b"FGH".to_vec()),
            ]]
        );
    }

    #[test]
    fn aligned_quality_section_keeps_zero_length_per_disabled_stream() {
        let pair = AlignedPair {
            first: ReadRecord::new("first", b"A".to_vec(), None),
            second: Some(ReadRecord::new("second", b"T".to_vec(), None)),
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: Vec::new(),
            index_reads: Vec::new(),
        };
        let options = crate::AlignmentCompressionOptions {
            preserve_quality: false,
            ..crate::AlignmentCompressionOptions::default()
        };
        let mut encoded = Vec::new();
        write_quality_stream(
            &mut encoded,
            &[pair],
            &[],
            &options,
            AlignmentLayout::PairedEnd,
        )
        .unwrap();

        assert_eq!(encoded, [0u8; 16]);
        assert_eq!(
            decode_qualities(&encoded, 1, false, AlignmentLayout::PairedEnd).unwrap(),
            vec![vec![None], vec![None]]
        );
    }

    #[test]
    fn aligned_quality_decoder_rejects_truncated_trailing_and_count_mismatch() {
        let stream = encode_stream_bytes(&[b"AB".as_slice()]).unwrap();
        let section = |stream: &[u8]| {
            let mut encoded = Vec::new();
            write_u64(&mut encoded, stream.len() as u64).unwrap();
            encoded.extend_from_slice(stream);
            encoded
        };

        let encoded = section(&stream);
        assert!(decode_qualities(
            &encoded[..encoded.len() - 1],
            1,
            true,
            AlignmentLayout::SingleEnd
        )
        .is_err());
        assert!(decode_qualities(&encoded, 2, true, AlignmentLayout::SingleEnd).is_err());

        let mut excessive_count = stream.clone();
        excessive_count[4..8].copy_from_slice(&2u32.to_le_bytes());
        assert!(decode_qualities(
            &section(&excessive_count),
            1,
            true,
            AlignmentLayout::SingleEnd
        )
        .is_err());

        let mut trailing_stream = stream;
        trailing_stream.push(0);
        assert!(decode_qualities(
            &section(&trailing_stream),
            1,
            true,
            AlignmentLayout::SingleEnd
        )
        .is_err());

        let mut trailing_section = encoded;
        trailing_section.push(0);
        assert!(decode_qualities(&trailing_section, 1, true, AlignmentLayout::SingleEnd).is_err());
    }

    #[test]
    fn chunk_round_trips_embedded_indexes() {
        for preserve_quality in [false, true] {
            let ranges = vec![
                EmbeddedIndexRange {
                    mate: EmbeddedIndexMate::R1,
                    start: 5,
                    end: 7,
                },
                EmbeddedIndexRange {
                    mate: EmbeddedIndexMate::R1,
                    start: 1,
                    end: 3,
                },
            ];
            let original_sequence = b"ACGTACGTAA".to_vec();
            let original_quality = b"ABCDEFGHIJ".to_vec();
            let mut record = ReadRecord::new(
                "read",
                original_sequence.clone(),
                Some(original_quality.clone()),
            );
            let indexes = split_pair_indexes(&mut record, None, &ranges, preserve_quality).unwrap();
            let descriptors = ranges
                .into_iter()
                .map(|range| IndexDescriptor { range })
                .collect::<Vec<_>>();
            let pair = AlignedPair {
                first: record,
                second: None,
                first_hit: None,
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes,
                index_reads: Vec::new(),
            };
            let options = crate::AlignmentCompressionOptions {
                preserve_quality,
                zstd_level: 1,
                ..crate::AlignmentCompressionOptions::default()
            };
            let encoded = encode_chunk(
                &[],
                &[pair],
                1,
                1,
                &options,
                IdMode::Independent,
                AlignmentLayout::SingleEnd,
                &descriptors,
            )
            .unwrap();
            let archive_chunk = EncodedArchiveChunk {
                first_record: 0,
                record_count: encoded.pairs as usize,
                raw_len: encoded.raw_len as usize,
                ids_len: encoded.ids_len as usize,
                qualities_len: encoded.qualities_len as usize,
                records_len: encoded.records_len as usize,
                compressed: encoded.compressed,
            };
            let mut decoded = decode_archive_chunk(
                archive_chunk,
                IdMode::Independent,
                AlignmentLayout::SingleEnd,
                preserve_quality,
                &[],
                &descriptors,
                1,
                1,
                0,
            )
            .unwrap();
            let (_, group) = decoded
                .next_group(AlignmentLayout::SingleEnd, None)
                .unwrap();
            assert!(group.index_reads.is_empty());
            let SequenceRecordGroup::Single(record) = group.sequence else {
                unreachable!();
            };
            assert_eq!(record.sequence, original_sequence);
            assert_eq!(
                record.quality,
                preserve_quality.then_some(original_quality.clone())
            );
        }
    }

    #[test]
    fn breakdown_aggregates_subsections_across_chunks() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        let index_path = directory.path().join("genome.seqarc-index");
        let reference: Vec<u8> = b"ACGT".iter().copied().cycle().take(300).collect();
        std::fs::write(
            &fasta,
            [b">chr1\n".as_slice(), reference.as_slice(), b"\n"].concat(),
        )
        .unwrap();
        crate::build_alignment_index(
            &fasta,
            &index_path,
            &crate::aligners::GenomeIndexBuildOptions::default(),
        )
        .unwrap();

        let first_path = directory.path().join("reads.1.fastq");
        let second_path = directory.path().join("reads.2.fastq");
        let archive_path = directory.path().join("reads.seqarc");
        let mut first = Vec::new();
        let mut second = Vec::new();
        for index in 0..6u8 {
            let start = usize::from(index) * 8;
            // Pairs 0-2 align to the reference; pairs 3-5 (poly-G, absent from
            // the reference) do not.
            let (first_sequence, second_sequence) = if index < 3 {
                (
                    reference[start..start + 12].to_vec(),
                    reference[start + 100..start + 112].to_vec(),
                )
            } else {
                (vec![b'A'; 12], vec![b'G'; 12])
            };
            for (buffer, mate, sequence) in [
                (&mut first, b"/1".as_slice(), &first_sequence),
                (&mut second, b"/2".as_slice(), &second_sequence),
            ] {
                buffer.extend_from_slice(b"@r");
                buffer.push(b'0' + index);
                buffer.extend_from_slice(mate);
                buffer.extend_from_slice(b"\n");
                buffer.extend_from_slice(sequence);
                buffer.extend_from_slice(b"\n+\n");
                buffer.extend(std::iter::repeat_n(b'I', sequence.len()));
                buffer.extend_from_slice(b"\n");
            }
        }
        std::fs::write(&first_path, &first).unwrap();
        std::fs::write(&second_path, &second).unwrap();

        let options = crate::AlignmentCompressionOptions {
            threads: 1,
            chunk_records: 3,
            // The whole R1 read is embedded index, mirroring barcode layouts.
            embedded_index_ranges: vec![EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 0,
                end: 12,
            }],
            ..crate::AlignmentCompressionOptions::default()
        };
        compress_paired_fastq_files_with_index(
            &index_path,
            &first_path,
            &second_path,
            &archive_path,
            &options,
        )
        .unwrap();

        let mut archive = File::open(&archive_path).unwrap();
        let info = read_archive_info_with_breakdown(&mut archive, true).unwrap();
        assert_eq!(info.records, 6);
        assert_eq!(info.declared_chunks, 2);
        let breakdown = info.breakdown.as_ref().unwrap();
        let subsections = &breakdown.subsections;

        // Six records across two chunks: one flags byte and one u16 length
        // per record per mate, and 12 index bases per record.
        assert_eq!(subsections.flags.raw, 6);
        assert_eq!(subsections.mates[0].lengths.raw, 12);
        assert_eq!(subsections.mates[1].lengths.raw, 12);
        // R1 is fully consumed by the index, so it stores no sequences.
        assert_eq!(subsections.mates[0].unaligned_sequences.raw, 0);
        assert_eq!(subsections.mates[0].positions.raw, 0);
        assert_eq!(subsections.embedded_index_ranges.len(), 1);
        assert_eq!(subsections.embedded_index_ranges[0].raw, 72);

        // Every subsection's sizes must sum to the record-section totals,
        // which include the subsection directory framing.
        let mut compressed_sum = subsections.flags.compressed;
        let mut raw_sum = subsections.flags.raw;
        for mate in &subsections.mates {
            for sizes in [
                mate.lengths,
                mate.positions,
                mate.noise_positions,
                mate.noise,
                mate.reference_skips,
                mate.emits,
                mate.reverses,
                mate.unaligned_sequences,
            ] {
                compressed_sum += sizes.compressed;
                raw_sum += sizes.raw;
            }
        }
        for sizes in &subsections.embedded_index_ranges {
            compressed_sum += sizes.compressed;
            raw_sum += sizes.raw;
        }
        for stream in &subsections.index_streams {
            for sizes in [stream.qualities, stream.lengths, stream.sequences] {
                compressed_sum += sizes.compressed;
                raw_sum += sizes.raw;
            }
        }
        let subsection_count = 1 + 2 * V21_MATE_SECTION_COUNT + 1;
        let framing = info.declared_chunks * (4 + 16 * subsection_count as u64);
        assert_eq!(compressed_sum + framing, breakdown.records.compressed);
        assert_eq!(raw_sum, breakdown.records.raw);
    }

    #[test]
    fn chunk_round_trips_fully_consumed_reads() {
        for preserve_quality in [false, true] {
            let ranges = vec![EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 0,
                end: 10,
            }];
            let original_sequence = b"ACGTACGTAA".to_vec();
            let original_quality = b"ABCDEFGHIJ".to_vec();
            let mut record = ReadRecord::new(
                "read",
                original_sequence.clone(),
                Some(original_quality.clone()),
            );
            let indexes = split_pair_indexes(&mut record, None, &ranges, preserve_quality).unwrap();
            assert!(record.sequence.is_empty());
            let descriptors = ranges
                .into_iter()
                .map(|range| IndexDescriptor { range })
                .collect::<Vec<_>>();
            let pair = AlignedPair {
                first: record,
                second: None,
                first_hit: None,
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes,
                index_reads: Vec::new(),
            };
            let options = crate::AlignmentCompressionOptions {
                preserve_quality,
                zstd_level: 1,
                ..crate::AlignmentCompressionOptions::default()
            };
            let encoded = encode_chunk(
                &[],
                &[pair],
                1,
                1,
                &options,
                IdMode::Independent,
                AlignmentLayout::SingleEnd,
                &descriptors,
            )
            .unwrap();
            let archive_chunk = EncodedArchiveChunk {
                first_record: 0,
                record_count: encoded.pairs as usize,
                raw_len: encoded.raw_len as usize,
                ids_len: encoded.ids_len as usize,
                qualities_len: encoded.qualities_len as usize,
                records_len: encoded.records_len as usize,
                compressed: encoded.compressed,
            };
            let mut decoded = decode_archive_chunk(
                archive_chunk,
                IdMode::Independent,
                AlignmentLayout::SingleEnd,
                preserve_quality,
                &[],
                &descriptors,
                1,
                1,
                0,
            )
            .unwrap();
            let (_, group) = decoded
                .next_group(AlignmentLayout::SingleEnd, None)
                .unwrap();
            assert!(group.index_reads.is_empty());
            let SequenceRecordGroup::Single(record) = group.sequence else {
                unreachable!();
            };
            assert_eq!(record.sequence, original_sequence);
            assert_eq!(
                record.quality,
                preserve_quality.then_some(original_quality.clone())
            );
        }
    }

    #[test]
    fn index_read_files_round_trip_through_chunks() {
        for (preserve_ids, preserve_quality) in [(true, true), (true, false), (false, true)] {
            let pairs = vec![
                AlignedPair {
                    first: ReadRecord::new("read0", b"ACGT".to_vec(), Some(b"ABCD".to_vec())),
                    second: None,
                    first_hit: None,
                    second_hit: None,
                    first_selection: AlignmentSelection::default(),
                    second_selection: AlignmentSelection::default(),
                    indexes: Vec::new(),
                    index_reads: vec![
                        ReadRecord::new("idx0", b"AACCGG".to_vec(), Some(b"abcdef".to_vec())),
                        ReadRecord::new("idx1", b"TTGGCC".to_vec(), Some(b"ghijkl".to_vec())),
                    ],
                },
                AlignedPair {
                    first: ReadRecord::new("read1", b"TTTT".to_vec(), Some(b"EFGH".to_vec())),
                    second: None,
                    first_hit: None,
                    second_hit: None,
                    first_selection: AlignmentSelection::default(),
                    second_selection: AlignmentSelection::default(),
                    indexes: Vec::new(),
                    index_reads: vec![
                        ReadRecord::new("idx2", b"ACGTAC".to_vec(), Some(b"mnopqr".to_vec())),
                        ReadRecord::new("idx3", b"GGTTAC".to_vec(), Some(b"stuvwx".to_vec())),
                    ],
                },
            ];
            let mode = if preserve_ids {
                IdMode::Independent
            } else {
                IdMode::None
            };
            let options = crate::AlignmentCompressionOptions {
                preserve_ids,
                id_prefix_chars: if preserve_ids { None } else { Some(0) },
                preserve_quality,
                zstd_level: 1,
                index_read_files: vec![
                    std::path::PathBuf::from("i1.fastq"),
                    std::path::PathBuf::from("i2.fastq"),
                ],
                ..crate::AlignmentCompressionOptions::default()
            };
            let encoded = encode_chunk(
                &[],
                &pairs,
                1,
                1,
                &options,
                mode,
                AlignmentLayout::SingleEnd,
                &[],
            )
            .unwrap();
            assert!(encoded.indexes_raw_len > 0);
            let archive_chunk = EncodedArchiveChunk {
                first_record: 0,
                record_count: encoded.pairs as usize,
                raw_len: encoded.raw_len as usize,
                ids_len: encoded.ids_len as usize,
                qualities_len: encoded.qualities_len as usize,
                records_len: encoded.records_len as usize,
                compressed: encoded.compressed,
            };
            let mut decoded = decode_archive_chunk(
                archive_chunk,
                mode,
                AlignmentLayout::SingleEnd,
                preserve_quality,
                &[],
                &[],
                1,
                1,
                2,
            )
            .unwrap();
            let (_, group) = decoded
                .next_group(AlignmentLayout::SingleEnd, None)
                .unwrap();
            assert_eq!(group.index_reads.len(), 2);
            let SequenceRecordGroup::Single(record) = group.sequence else {
                unreachable!();
            };
            assert_eq!(record.sequence, b"ACGT");
            assert_eq!(record.quality, preserve_quality.then_some(b"ABCD".to_vec()));
            assert_eq!(
                record.id,
                if preserve_ids {
                    &b"read0"[..]
                } else {
                    &b""[..]
                }
            );
            assert_eq!(group.index_reads[0].sequence, b"AACCGG");
            assert_eq!(group.index_reads[1].sequence, b"TTGGCC");
            assert_eq!(
                group.index_reads[0].quality,
                preserve_quality.then_some(b"abcdef".to_vec())
            );
            assert_eq!(
                group.index_reads[0].id,
                if preserve_ids {
                    &b"read0"[..]
                } else {
                    &b""[..]
                }
            );
            let (_, group) = decoded
                .next_group(AlignmentLayout::SingleEnd, None)
                .unwrap();
            let SequenceRecordGroup::Single(record) = group.sequence else {
                unreachable!();
            };
            assert_eq!(record.sequence, b"TTTT");
            assert_eq!(group.index_reads[0].sequence, b"ACGTAC");
            assert_eq!(group.index_reads[1].sequence, b"GGTTAC");
            assert_eq!(
                group.index_reads[1].id,
                if preserve_ids {
                    &b"read1"[..]
                } else {
                    &b""[..]
                }
            );
            assert!(decoded
                .next_group(AlignmentLayout::SingleEnd, None)
                .is_none());
        }
    }

    #[test]
    fn reads_v21_header_layout() {
        let options = crate::AlignmentCompressionOptions::default();
        let reference = encode_reference(b"ACGT", &options).unwrap().block;
        let range = EmbeddedIndexRange {
            mate: EmbeddedIndexMate::R1,
            start: 1,
            end: 3,
        };
        let descriptors = vec![IndexDescriptor {
            range: range.clone(),
        }];
        let metadata = encode_index_metadata(&descriptors).unwrap();
        let mut header = Vec::new();
        header.extend(MAGIC);
        write_u16(&mut header, VERSION_MAJOR).unwrap();
        write_u16(&mut header, VERSION_MINOR).unwrap();
        write_u16(&mut header, 4).unwrap();
        header.push(0);
        header.extend([1, 1]);
        write_u64(&mut header, reference.len() as u64).unwrap();
        write_u64(&mut header, metadata.len() as u64).unwrap();
        write_u64(&mut header, 9).unwrap();
        assert_eq!(header.len(), PREAMBLE_BYTES);
        let crc = preamble_crc32c(&header, &reference, &metadata);
        write_u32(&mut header, crc).unwrap();
        header.extend(reference);
        header.extend(metadata);
        let decoded = read_header(&mut Cursor::new(header)).unwrap();
        assert_eq!(decoded.chunk_count, 9);
        assert_eq!(decoded.index_streams, 0);
        assert_eq!(decoded.indexes.len(), 1);
        assert_eq!(decoded.indexes[0].range, range);
    }

    #[test]
    fn preamble_crc_detects_corruption() {
        let options = crate::AlignmentCompressionOptions::default();
        let reference = encode_reference(b"ACGT", &options).unwrap().block;
        let metadata = encode_index_metadata(&[]).unwrap();
        let build = |flip: Option<usize>| -> Vec<u8> {
            let mut header = Vec::new();
            header.extend(MAGIC);
            write_u16(&mut header, VERSION_MAJOR).unwrap();
            write_u16(&mut header, VERSION_MINOR).unwrap();
            write_u16(&mut header, 4).unwrap();
            header.push(0);
            header.extend([1, 1]);
            write_u64(&mut header, reference.len() as u64).unwrap();
            write_u64(&mut header, metadata.len() as u64).unwrap();
            write_u64(&mut header, 9).unwrap();
            let crc = preamble_crc32c(&header, &reference, &metadata);
            write_u32(&mut header, crc).unwrap();
            header.extend(&reference);
            header.extend(&metadata);
            if let Some(index) = flip {
                header[index] ^= 0xff;
            }
            header
        };
        read_header(&mut Cursor::new(build(None))).unwrap();
        // Flips in magic, version, or flags are caught by earlier validation,
        // so only positions guarded solely by the CRC are exercised here.
        for flip in [
            34,                                   // chunk count low byte
            35,                                   // chunk count second byte
            42,                                   // stored CRC field
            PREAMBLE_BYTES + 3,                   // reference block payload
            PREAMBLE_BYTES + reference.len() + 2, // index metadata payload
        ] {
            let error = match read_header(&mut Cursor::new(build(Some(flip)))) {
                Err(error) => error,
                Ok(_) => panic!("flip at {flip} was not detected"),
            };
            assert!(
                error.is_checksum(),
                "flip at {flip} did not fail the preamble CRC: {error}"
            );
        }
    }

    #[test]
    fn archive_operations_normalize_and_pack_in_archive_order() {
        let normalized = archive_normalize([
            CigarOp {
                kind: CigarKind::Eq,
                len: 2,
            },
            CigarOp {
                kind: CigarKind::Diff,
                len: 3,
            },
            CigarOp {
                kind: CigarKind::Ins,
                len: 1,
            },
            CigarOp {
                kind: CigarKind::SoftClip,
                len: 2,
            },
            CigarOp {
                kind: CigarKind::Del,
                len: 4,
            },
            CigarOp {
                kind: CigarKind::RefSkip,
                len: 5,
            },
            CigarOp {
                kind: CigarKind::HardClip,
                len: 9,
            },
            CigarOp {
                kind: CigarKind::Pad,
                len: 9,
            },
            CigarOp {
                kind: CigarKind::Match,
                len: 0,
            },
        ])
        .unwrap();
        assert_eq!(
            normalized,
            vec![
                ArchiveRun {
                    op: ArchiveOp::Match,
                    len: 5,
                },
                ArchiveRun {
                    op: ArchiveOp::Literal,
                    len: 3,
                },
                ArchiveRun {
                    op: ArchiveOp::RefAdvance,
                    len: 9,
                },
            ]
        );
        assert_eq!(
            [ArchiveOp::Match, ArchiveOp::Literal, ArchiveOp::RefAdvance].map(|op| op as u8),
            [0, 1, 2]
        );
    }

    #[test]
    fn chunk_header_is_48_bytes_with_crc() {
        let payload = vec![7u8; 13];
        let mut encoded = Vec::new();
        write_u32(&mut encoded, 2).unwrap();
        for value in [13u64, 0, 0, 0, 0] {
            write_u64(&mut encoded, value).unwrap();
        }
        let header = encoded.clone();
        let crc = chunk_crc32c(&header, &payload);
        write_u32(&mut encoded, crc).unwrap();
        encoded.extend(&payload);
        assert_eq!(encoded.len(), 48 + 13);
        let chunk = read_chunk_header(&mut Cursor::new(&encoded)).unwrap();
        assert_eq!(chunk.record_count, 2);
        assert_eq!(chunk.compressed_len, 13);
        let decoded = read_chunk_payload(&mut Cursor::new(&encoded[48..]), &chunk, 7, 3).unwrap();
        assert_eq!(decoded, payload);
        // A flipped payload byte must fail the checksum.
        let mut damaged = encoded.clone();
        let last = damaged.len() - 1;
        damaged[last] ^= 1;
        let error = read_chunk_payload(&mut Cursor::new(&damaged[48..]), &chunk, 7, 3).unwrap_err();
        assert!(error.is_checksum());
        // A flipped header byte must fail the checksum too.
        let mut damaged = encoded.clone();
        damaged[0] ^= 1;
        let error = read_chunk_header(&mut Cursor::new(&damaged))
            .and_then(|header| read_chunk_payload(&mut Cursor::new(&damaged[48..]), &header, 7, 3))
            .unwrap_err();
        assert!(error.is_checksum());
    }

    #[test]
    fn chunk_checksum_error_names_the_chunk() {
        let payload = vec![1u8; 4];
        let mut encoded = Vec::new();
        write_u32(&mut encoded, 1).unwrap();
        for value in [4u64, 0, 0, 0, 0] {
            write_u64(&mut encoded, value).unwrap();
        }
        let crc = chunk_crc32c(&encoded, &payload);
        write_u32(&mut encoded, crc).unwrap();
        encoded.extend(&payload);
        let mut damaged = encoded;
        let last = damaged.len() - 1;
        damaged[last] ^= 0xff;
        let header = read_chunk_header(&mut Cursor::new(&damaged)).unwrap();
        let error =
            read_chunk_payload(&mut Cursor::new(&damaged[48..]), &header, 9, 4).unwrap_err();
        assert!(
            error.to_string().contains("chunk 4 (record 9)"),
            "unexpected message: {error}"
        );
    }

    #[test]
    fn corrupt_chunk_length_fails_without_huge_allocation() {
        let payload = vec![7u8; 13];
        let mut encoded = Vec::new();
        write_u32(&mut encoded, 2).unwrap();
        for value in [13u64, 0, 0, 0, 0] {
            write_u64(&mut encoded, value).unwrap();
        }
        let header = encoded.clone();
        let crc = chunk_crc32c(&header, &payload);
        write_u32(&mut encoded, crc).unwrap();
        encoded.extend(&payload);
        // Corrupt compressed_len to usize::MAX; the reader must hit the end of
        // the input and return an error instead of attempting the allocation.
        let mut damaged = encoded;
        damaged[4..12].copy_from_slice(&u64::MAX.to_le_bytes());
        let corrupt = read_chunk_header(&mut Cursor::new(&damaged)).unwrap();
        let error = read_chunk_payload(&mut Cursor::new(&payload), &corrupt, 0, 0).unwrap_err();
        assert!(matches!(error, Error::Io(_)), "unexpected error: {error}");
    }

    #[test]
    fn corrupt_reference_length_fails_without_huge_allocation() {
        let options = crate::AlignmentCompressionOptions::default();
        let reference = encode_reference(b"ACGT", &options).unwrap().block;
        let metadata = encode_index_metadata(&[]).unwrap();
        let mut header = Vec::new();
        header.extend(MAGIC);
        write_u16(&mut header, VERSION_MAJOR).unwrap();
        write_u16(&mut header, VERSION_MINOR).unwrap();
        write_u16(&mut header, 4).unwrap();
        header.push(0);
        header.extend([1, 1]);
        // Declare a reference block far larger than the input provides.
        write_u64(&mut header, u64::MAX).unwrap();
        write_u64(&mut header, metadata.len() as u64).unwrap();
        write_u64(&mut header, 0).unwrap();
        let crc = preamble_crc32c(&header, &reference, &metadata);
        write_u32(&mut header, crc).unwrap();
        header.extend(&reference);
        header.extend(&metadata);
        // The read must fail at the end of the input before the CRC check;
        // with the old pre-CRC allocation this attempted a usize::MAX buffer.
        let error = match read_header(&mut Cursor::new(&header)) {
            Err(error) => error,
            Ok(_) => panic!("expected corrupt reference length to fail"),
        };
        assert!(matches!(error, Error::Io(_)), "unexpected error: {error}");
    }

    #[test]
    fn verify_reports_io_errors_instead_of_truncation() {
        let options = crate::AlignmentCompressionOptions::default();
        let reference = encode_reference(b"ACGT", &options).unwrap().block;
        let metadata = encode_index_metadata(&[]).unwrap();
        let mut header = Vec::new();
        header.extend(MAGIC);
        write_u16(&mut header, VERSION_MAJOR).unwrap();
        write_u16(&mut header, VERSION_MINOR).unwrap();
        write_u16(&mut header, 4).unwrap();
        header.push(0);
        header.extend([1, 1]);
        write_u64(&mut header, reference.len() as u64).unwrap();
        write_u64(&mut header, metadata.len() as u64).unwrap();
        write_u64(&mut header, 1).unwrap();
        let crc = preamble_crc32c(&header, &reference, &metadata);
        write_u32(&mut header, crc).unwrap();
        header.extend(&reference);
        header.extend(&metadata);
        // Serves the archive bytes, then fails with an I/O error instead of
        // reporting a clean end of file.
        struct FailAfterEof<'a> {
            data: Cursor<&'a [u8]>,
        }
        impl Read for FailAfterEof<'_> {
            fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
                let read = self.data.read(buf)?;
                if read == 0 {
                    return Err(std::io::Error::from(std::io::ErrorKind::PermissionDenied));
                }
                Ok(read)
            }
        }
        let error = verify_aligned_archive(&mut FailAfterEof {
            data: Cursor::new(&header),
        })
        .unwrap_err();
        assert!(matches!(error, Error::Io(_)), "unexpected error: {error}");
        // A clean early end of file during the chunk scan is still reported
        // as truncation; append a partial chunk header and cut there.
        let mut truncated_archive = header.clone();
        truncated_archive.extend(std::iter::repeat_n(0u8, 10));
        let report = verify_aligned_archive(&mut Cursor::new(&truncated_archive)).unwrap();
        assert!(report.truncated);
    }

    #[test]
    fn compact_boundaries_include_only_retained_bases() {
        let runs = vec![
            KeptRun {
                original_start: 5,
                compact_start: 0,
                length: 3,
            },
            KeptRun {
                original_start: 10,
                compact_start: 3,
                length: 2,
            },
        ];
        for (original, compact) in [
            (0, 0),
            (5, 0),
            (6, 1),
            (8, 3),
            (9, 3),
            (10, 3),
            (12, 5),
            (20, 5),
        ] {
            assert_eq!(compact_boundary(&runs, original, 20).unwrap(), compact);
        }
        assert!(compact_boundary(&runs, 21, 20).is_err());
    }

    #[test]
    fn effective_cigar_is_canonical() {
        let hit = GenomeAlignment {
            start: 0,
            query_start: 2,
            query_end: 8,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 1,
                },
                CigarOp {
                    kind: CigarKind::Eq,
                    len: 2,
                },
                CigarOp {
                    kind: CigarKind::Diff,
                    len: 1,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 2,
                },
                CigarOp {
                    kind: CigarKind::HardClip,
                    len: 4,
                },
                CigarOp {
                    kind: CigarKind::Pad,
                    len: 1,
                },
                CigarOp {
                    kind: CigarKind::Del,
                    len: 0,
                },
            ],
            is_reverse: false,
        };
        assert_eq!(
            effective_cigar(&hit, 10).unwrap(),
            vec![
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 3,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 5,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: 2,
                },
            ]
        );
    }

    #[test]
    fn uncovered_deletions_and_reference_skips_are_pruned() {
        let reference = [b"AAA".as_slice(), &[b'C'; 10], b"GGG".as_slice()].concat();
        for kind in [CigarKind::Del, CigarKind::RefSkip] {
            for reverse in [false, true] {
                let hit = GenomeAlignment {
                    start: 0,
                    query_start: 0,
                    query_end: 6,
                    cigar: vec![
                        CigarOp {
                            kind: CigarKind::Match,
                            len: 3,
                        },
                        CigarOp { kind, len: 10 },
                        CigarOp {
                            kind: CigarKind::Match,
                            len: 3,
                        },
                    ],
                    is_reverse: reverse,
                };
                let mut mask = ReferenceMask::new(reference.len());
                mask.mark_hit(&hit).unwrap();
                let (compact, runs) = compact_reference(reference.as_slice(), &mask).unwrap();
                assert_eq!(compact, b"AAAGGG");

                let transformed = SpoolAlignment::from_genome(hit)
                    .unwrap()
                    .into_genome(&runs, reference.len())
                    .unwrap();
                assert_eq!(
                    transformed.cigar,
                    vec![CigarOp {
                        kind: CigarKind::Match,
                        len: 6,
                    }]
                );
            }
        }
    }

    #[test]
    fn reference_skip_retains_bases_used_by_an_overlapping_alignment() {
        let reference = [b"AAA".as_slice(), &[b'C'; 10], b"GGG".as_slice()].concat();
        let skipped = GenomeAlignment {
            start: 0,
            query_start: 0,
            query_end: 6,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 3,
                },
                CigarOp {
                    kind: CigarKind::RefSkip,
                    len: 10,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 3,
                },
            ],
            is_reverse: false,
        };
        let overlapping = GenomeAlignment {
            start: 6,
            query_start: 0,
            query_end: 4,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 4,
            }],
            is_reverse: false,
        };
        let mut mask = ReferenceMask::new(reference.len());
        mask.mark_hit(&skipped).unwrap();
        mask.mark_hit(&overlapping).unwrap();
        let (compact, runs) = compact_reference(reference.as_slice(), &mask).unwrap();
        assert_eq!(compact, b"AAACCCCGGG");

        let skipped = SpoolAlignment::from_genome(skipped)
            .unwrap()
            .into_genome(&runs, reference.len())
            .unwrap();
        assert_eq!(
            skipped.cigar,
            vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 3,
                },
                CigarOp {
                    kind: CigarKind::RefSkip,
                    len: 4,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 3,
                },
            ]
        );
        let overlapping = SpoolAlignment::from_genome(overlapping)
            .unwrap()
            .into_genome(&runs, reference.len())
            .unwrap();
        assert_eq!(overlapping.start, 3);
    }

    #[test]
    fn leading_and_trailing_skips_do_not_retain_reference() {
        let reference = [b"CCCCC".as_slice(), b"AAA", b"GGGGG"].concat();
        let hit = GenomeAlignment {
            start: 0,
            query_start: 0,
            query_end: 3,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Del,
                    len: 5,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 3,
                },
                CigarOp {
                    kind: CigarKind::RefSkip,
                    len: 5,
                },
            ],
            is_reverse: false,
        };
        let mut mask = ReferenceMask::new(reference.len());
        mask.mark_hit(&hit).unwrap();
        let (compact, runs) = compact_reference(reference.as_slice(), &mask).unwrap();
        assert_eq!(compact, b"AAA");
        let transformed = SpoolAlignment::from_genome(hit)
            .unwrap()
            .into_genome(&runs, reference.len())
            .unwrap();
        assert_eq!(transformed.start, 0);
        assert_eq!(
            transformed.cigar,
            vec![CigarOp {
                kind: CigarKind::Match,
                len: 3,
            }]
        );
    }

    #[test]
    fn id_modes_round_trip_losslessly() {
        for (mode, first_id, second_id, expected_first, expected_second) in [
            (
                IdMode::IdenticalMates,
                b"read".to_vec(),
                b"read".to_vec(),
                b"read".to_vec(),
                b"read".to_vec(),
            ),
            (
                IdMode::SlashMate,
                b"read/1".to_vec(),
                b"read/2".to_vec(),
                b"read/1".to_vec(),
                b"read/2".to_vec(),
            ),
            (
                IdMode::Independent,
                b"first".to_vec(),
                b"second".to_vec(),
                b"first".to_vec(),
                b"second".to_vec(),
            ),
        ] {
            let pair = AlignedPair {
                first: ReadRecord::new(first_id, b"A".to_vec(), None),
                second: Some(ReadRecord::new(second_id, b"A".to_vec(), None)),
                first_hit: None,
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
                index_reads: Vec::new(),
            };
            let mut encoded = Vec::new();
            write_id_stream(&mut encoded, &[pair], mode, AlignmentLayout::PairedEnd).unwrap();
            let decoded = decode_ids(&encoded, 1, mode, AlignmentLayout::PairedEnd).unwrap();
            assert_eq!(decoded[0][0], expected_first);
            assert_eq!(decoded[1][0], expected_second);
        }
    }

    #[test]
    fn alignment_cost_counts_raw_literal_bytes() {
        assert_eq!(
            alignment_candidate_core_cost(
                b"ACGTNN",
                None::<(&GenomeAlignment, &[u8])>,
                &[ArchiveRun {
                    op: ArchiveOp::Literal,
                    len: 6,
                }],
            )
            .unwrap(),
            9
        );
    }

    #[test]
    fn packed_sequence_uses_sparse_exceptions_when_smaller() {
        let sparse = b"ACGTACGTACGTACGTACGTNCGTACGTACGTACGTACGT";
        let encoded = encode_packed(sparse).unwrap();
        assert_eq!(encoded[4], 2);
        assert_eq!(
            read_packed(&mut Cursor::new(encoded)).unwrap().decode(),
            sparse
        );

        let raw = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
        let encoded = encode_packed(raw).unwrap();
        assert_eq!(encoded[4], 0);
        assert_eq!(
            read_packed(&mut Cursor::new(encoded)).unwrap().decode(),
            raw
        );
    }

    #[test]
    fn packed_and_raw_sequence_modes_round_trip() {
        let mut packed = Vec::new();
        write_u32(&mut packed, 4).unwrap();
        write_u8(&mut packed, 1).unwrap();
        write_u64(&mut packed, 1).unwrap();
        packed.push(0x1b);
        assert_eq!(
            read_packed(&mut Cursor::new(packed)).unwrap().decode(),
            b"ACGT"
        );

        let mut raw = Vec::new();
        write_u32(&mut raw, 4).unwrap();
        write_u8(&mut raw, 0).unwrap();
        write_u64(&mut raw, 4).unwrap();
        raw.extend(b"ACNT");
        assert_eq!(
            read_packed(&mut Cursor::new(raw)).unwrap().decode(),
            b"ACNT"
        );
    }

    #[test]
    fn rejects_non_v0_28_archives() {
        for (major, minor) in [
            (0, 0),
            (0, 7),
            (0, 26),
            (0, 27),
            (0, 29),
            (1, 0),
            (1, 28),
            (2, 3),
            (u16::MAX, u16::MAX),
        ] {
            let mut header = Vec::new();
            header.extend(MAGIC);
            write_u16(&mut header, major).unwrap();
            write_u16(&mut header, minor).unwrap();
            let message = match read_header(&mut Cursor::new(header)) {
                Err(Error::InvalidFormat(message)) => message,
                Err(Error::Checksum(message)) => {
                    panic!("unexpected checksum error: {message}")
                }
                Err(Error::InvalidInput(message)) => {
                    panic!("unexpected input error: {message}")
                }
                Err(Error::Io(error)) => panic!("unexpected I/O error: {error}"),
                Ok(_) => panic!("expected version {major}.{minor} to be rejected"),
            };
            assert_eq!(
                message,
                format!(
                    "unsupported alignment archive version {major}.{minor}; expected {}.{}",
                    VERSION_MAJOR, VERSION_MINOR
                )
            );
        }
    }

    #[test]
    fn rejects_foreign_alignment_magic() {
        let mut header = vec![0; MAGIC.len()];
        write_u16(&mut header, VERSION_MAJOR).unwrap();
        write_u16(&mut header, VERSION_MINOR).unwrap();
        assert!(matches!(
            read_header(&mut Cursor::new(header)),
            Err(Error::InvalidFormat(message)) if message == "bad alignment archive magic"
        ));
    }

    #[test]
    fn v21_aligned_read_round_trips_noise_skips_emits_and_reverse() {
        let reference = b"AACCGGTTAACCGGTT";
        for reverse in [false, true] {
            let mut oriented = b"AATCGAACT".to_vec();
            let hit = GenomeAlignment {
                start: 0,
                query_start: 0,
                query_end: oriented.len() as i32,
                is_reverse: reverse,
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 3,
                    },
                    CigarOp {
                        kind: CigarKind::Del,
                        len: 2,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 2,
                    },
                    CigarOp {
                        kind: CigarKind::Ins,
                        len: 2,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 2,
                    },
                ],
            };
            if reverse {
                reverse_complement(&mut oriented);
            }
            let record = ReadRecord::new("read", oriented.clone(), None);
            let body = encode_v21_aligned_read(reference.as_slice(), &record, &hit).unwrap();
            assert_eq!(
                decode_v21_aligned_read(reference, record.sequence.len(), &body).unwrap(),
                record.sequence
            );
        }
    }

    #[test]
    fn v21_noise_deltas_keep_trailing_emit_after_match() {
        let reference = vec![b'C'; 57];
        let sequence = [reference.as_slice(), &[b'A'; 43]].concat();
        let record = ReadRecord::new("read", sequence.clone(), None);
        let hit = GenomeAlignment {
            start: 0,
            query_start: 0,
            query_end: 57,
            is_reverse: false,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 57,
            }],
        };
        let encoded = encode_v21_aligned_read(reference.as_slice(), &record, &hit).unwrap();
        assert_eq!(
            decode_v21_aligned_read(reference.as_slice(), sequence.len(), &encoded).unwrap(),
            sequence
        );
    }

    #[test]
    fn v21_noise_deltas_do_not_repeat_prefix_after_reference_skip() {
        let reference = vec![b'C'; 12];
        let sequence = b"CCCCCACC".to_vec();
        let record = ReadRecord::new("read", sequence.clone(), None);
        let hit = GenomeAlignment {
            start: 0,
            query_start: 0,
            query_end: sequence.len() as i32,
            is_reverse: false,
            cigar: vec![
                CigarOp {
                    kind: CigarKind::Match,
                    len: 3,
                },
                CigarOp {
                    kind: CigarKind::Del,
                    len: 1,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 2,
                },
                CigarOp {
                    kind: CigarKind::Ins,
                    len: 1,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: 2,
                },
            ],
        };
        let encoded = encode_v21_aligned_read(reference.as_slice(), &record, &hit).unwrap();
        assert_eq!(
            decode_v21_aligned_read(reference.as_slice(), sequence.len(), &encoded).unwrap(),
            sequence
        );
    }

    #[test]
    fn v21_record_sections_are_independent_zstd_frames_and_cover_partial_pairs() {
        let reference = b"AACCGGTT";
        let hit = |start| GenomeAlignment {
            start,
            query_start: 0,
            query_end: 4,
            is_reverse: false,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 4,
            }],
        };
        let pairs = vec![
            AlignedPair {
                first: ReadRecord::new("a/1", b"AACC".to_vec(), None),
                second: Some(ReadRecord::new("a/2", b"NNNN".to_vec(), None)),
                first_hit: Some(hit(0)),
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
                index_reads: Vec::new(),
            },
            AlignedPair {
                first: ReadRecord::new("b/1", b"NNNN".to_vec(), None),
                second: Some(ReadRecord::new("b/2", b"GGTT".to_vec(), None)),
                first_hit: None,
                second_hit: Some(hit(4)),
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
                index_reads: Vec::new(),
            },
        ];
        let encoded = write_v21_record_stream(
            reference.as_slice(),
            &pairs,
            1,
            1,
            AlignmentLayout::PairedEnd,
            &[],
            true,
            0,
            1,
        )
        .unwrap()
        .bytes;
        let mut cursor = Cursor::new(encoded.clone());
        let section_count = 1 + 2 * V21_MATE_SECTION_COUNT;
        assert_eq!(read_u32(&mut cursor).unwrap(), section_count as u32);
        for _ in 0..section_count {
            let compressed = read_u64(&mut cursor).unwrap() as usize;
            let raw = read_u64(&mut cursor).unwrap() as usize;
            let mut frame = vec![0; compressed];
            cursor.read_exact(&mut frame).unwrap();
            assert_eq!(
                zstd::stream::decode_all(Cursor::new(frame)).unwrap().len(),
                raw
            );
        }
        assert_eq!(cursor.position(), encoded.len() as u64);
        let (decoded, index_streams) = {
            let sections = decode_v21_records(
                &encoded,
                pairs.len(),
                reference,
                AlignmentLayout::PairedEnd,
                &[],
                true,
                0,
                1,
                1,
            )
            .unwrap();
            (sections.biological, sections.index_streams)
        };
        assert!(index_streams.is_empty());
        assert_eq!(decoded[0][0], b"AACC");
        assert_eq!(decoded[1][0], b"NNNN");
        assert_eq!(decoded[0][1], b"NNNN");
        assert_eq!(decoded[1][1], b"GGTT");
    }

    #[test]
    fn v21_paired_r2_position_and_reverse_flag_are_relative_to_r1() {
        let reference = b"AACCGGTT";
        let hit = |start, is_reverse| GenomeAlignment {
            start,
            query_start: 0,
            query_end: 4,
            is_reverse,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 4,
            }],
        };
        let pairs = vec![
            AlignedPair {
                first: ReadRecord::new("a/1", b"CCAA".to_vec(), None),
                second: Some(ReadRecord::new("a/2", b"AACC".to_vec(), None)),
                first_hit: Some(hit(4, true)),
                second_hit: Some(hit(0, false)),
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
                index_reads: Vec::new(),
            },
            AlignedPair {
                first: ReadRecord::new("b/1", b"CCAA".to_vec(), None),
                second: Some(ReadRecord::new("b/2", b"GGTT".to_vec(), None)),
                first_hit: Some(hit(4, true)),
                second_hit: Some(hit(0, true)),
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
                index_reads: Vec::new(),
            },
        ];
        let encoded = write_v21_record_stream(
            reference.as_slice(),
            &pairs,
            1,
            1,
            AlignmentLayout::PairedEnd,
            &[],
            true,
            0,
            1,
        )
        .unwrap()
        .bytes;
        let mut cursor = Cursor::new(encoded.clone());
        let section_count = read_u32(&mut cursor).unwrap() as usize;
        let mut sections = Vec::with_capacity(section_count);
        for _ in 0..section_count {
            let compressed_len = read_u64(&mut cursor).unwrap() as usize;
            let raw_len = read_u64(&mut cursor).unwrap() as usize;
            let mut compressed = vec![0; compressed_len];
            cursor.read_exact(&mut compressed).unwrap();
            let raw = zstd::stream::decode_all(Cursor::new(compressed)).unwrap();
            assert_eq!(raw.len(), raw_len);
            sections.push(raw);
        }
        assert_eq!(sections[10], [7, 7]);
        assert_eq!(sections[15], [0, 1]);

        let (decoded, index_streams) = {
            let sections = decode_v21_records(
                &encoded,
                pairs.len(),
                reference,
                AlignmentLayout::PairedEnd,
                &[],
                true,
                0,
                1,
                1,
            )
            .unwrap();
            (sections.biological, sections.index_streams)
        };
        assert!(index_streams.is_empty());
        assert_eq!(decoded[0][0], b"CCAA");
        assert_eq!(decoded[1][0], b"AACC");
        assert_eq!(decoded[0][1], b"CCAA");
        assert_eq!(decoded[1][1], b"GGTT");
    }

    #[test]
    fn position_width_depends_on_pruned_reference_length() {
        for (length, width) in [
            (0, 1),
            (1, 1),
            (256, 1),
            (257, 2),
            (65_536, 2),
            (65_537, 3),
            (16_777_216, 3),
            (16_777_217, 4),
            (u32::MAX as usize + 1, 4),
            (u32::MAX as usize + 2, 8),
        ] {
            assert_eq!(position_width(length), width);
        }
        for (length, width) in [
            (1, 1),
            (128, 1),
            (129, 2),
            (32_768, 2),
            (32_769, 3),
            (8_388_608, 3),
            (8_388_609, 4),
        ] {
            assert_eq!(delta_width(length), width);
        }
    }

    #[test]
    fn positions_round_trip_at_each_supported_width() {
        for width in [1, 2, 3, 4, 8] {
            let max = if width == 8 {
                u64::MAX
            } else {
                (1u64 << (width * 8)) - 1
            };
            let values = [0, max / 2, max];
            let mut encoded = Vec::new();
            for value in values {
                write_position(&mut encoded, value, width).unwrap();
            }
            let mut cursor = Cursor::new(encoded);
            for value in values {
                assert_eq!(read_position(&mut cursor, width).unwrap(), value);
            }
        }
    }

    #[test]
    fn alignment_selection_uses_smallest_candidate_payload() {
        let reference = [vec![b'A'; 40], vec![b'C'; 40]].concat();
        let record = ReadRecord::new("read", vec![b'C'; 40], None);
        let hit = |start| minibwa::Hit {
            tid: 0,
            contig: Some("genome".to_owned()),
            ref_start: start,
            ref_end: start + 40,
            query_start: 0,
            query_end: 40,
            strand: Strand::Forward,
            mapq: 60,
            score: if start == 40 { 40 } else { 1 },
            n_sub: 0,
            blen: 40,
            mlen: if start == 40 { 40 } else { 0 },
            proper_pair: false,
            is_primary: start == 0,
            is_secondary: false,
            is_supplementary: false,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 40,
            }],
        };
        let reference = reference.as_slice();
        let selected = select_alignment(&record, vec![hit(0), hit(40)], reference).unwrap();
        assert_eq!(selected.hit.unwrap().start, 40);
        assert!(selected.selection.selected_non_first);
        assert!(!selected.selection.rejected_alignment);
    }

    #[test]
    fn alignment_selection_rejects_payload_larger_than_raw_sequence() {
        let reference = vec![b'A'; 40];
        let record = ReadRecord::new("read", vec![b'N'; 40], None);
        let hit = minibwa::Hit {
            tid: 0,
            contig: Some("genome".to_owned()),
            ref_start: 0,
            ref_end: 40,
            query_start: 0,
            query_end: 40,
            strand: Strand::Forward,
            mapq: 1,
            score: 1,
            n_sub: 0,
            blen: 40,
            mlen: 0,
            proper_pair: false,
            is_primary: true,
            is_secondary: false,
            is_supplementary: false,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 40,
            }],
        };
        let reference = reference.as_slice();
        let selected = select_alignment(&record, vec![hit], reference).unwrap();
        assert!(selected.hit.is_none());
        assert!(selected.selection.rejected_alignment);
    }
}
