use crate::index_read::{self, IndexReadEncoding, IndexReadError, Whitelist, WhitelistBuilder};
use crate::minibwa;
use crate::minibwa::{Aligner, CigarKind, CigarOp, Index, Opts, Strand, ThreadBuf};
#[cfg(test)]
use crate::packed_sequence::PackedSequence;
use crate::quality::{decode_stream_bytes, encode_stream_bytes};
use crate::{
    apply_quality, open_fastx_reader, truncate_read_id, AlignedCompressionProfile,
    AlignedDecompressionMetadata, AlignedRecordGroup, AlignedRecordSink, AlignmentLayout,
    ArchiveRecords, CompressionMetadata, DecompressionOptions, EmbeddedIndexMate,
    EmbeddedIndexRange, Error, FastxReader, Format, GenomeAlignment, PairedRecordSink, ReadRecord,
    Result,
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
const VERSION: u32 = 20;
const MAP_TASK_READS: usize = 256;
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
}

#[derive(Clone, Copy, Debug, Default)]
struct AlignmentSelection {
    candidate_hits: u32,
    selected_non_first: bool,
    rejected_alignment: bool,
}

struct SelectedAlignment {
    hit: Option<GenomeAlignment>,
    selection: AlignmentSelection,
}

struct MappingProgressStats {
    layout: AlignmentLayout,
    processed: [u64; 2],
    aligned: [u64; 2],
}

impl MappingProgressStats {
    fn new(layout: AlignmentLayout) -> Self {
        Self {
            layout,
            processed: [0; 2],
            aligned: [0; 2],
        }
    }

    fn observe(&mut self, mate: usize, selected: &[SelectedAlignment], progress: &ProgressBar) {
        self.processed[mate] = self.processed[mate].saturating_add(selected.len() as u64);
        self.aligned[mate] = self.aligned[mate].saturating_add(
            selected
                .iter()
                .filter(|selected| selected.selection.candidate_hits != 0)
                .count() as u64,
        );
        progress.set_message(self.message());
        progress.inc(selected.len() as u64);
    }

    fn message(&self) -> String {
        let percent = |mate: usize| {
            if self.processed[mate] == 0 {
                0.0
            } else {
                self.aligned[mate] as f64 * 100.0 / self.processed[mate] as f64
            }
        };
        if self.layout == AlignmentLayout::PairedEnd {
            format!(
                "R1 aligned {:.2}% | R2 aligned {:.2}%",
                percent(0),
                percent(1)
            )
        } else {
            format!("R1 aligned {:.2}%", percent(0))
        }
    }
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolPair {
    first: SpoolRead,
    second: Option<SpoolRead>,
    indexes: Vec<SpoolIndexRead>,
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

#[derive(Clone, Debug)]
struct IndexEncodingDescriptor {
    range: EmbeddedIndexRange,
    local_whitelist: Whitelist,
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

struct ReferenceMask {
    bits: Vec<u64>,
    length: usize,
}

#[derive(Clone, Copy, Debug)]
struct KeptRun {
    original_start: usize,
    compact_start: usize,
    length: usize,
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
    compressed: Vec<u8>,
}

struct EncodedReference {
    block: Vec<u8>,
    raw_len: u64,
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
    Whitelist = 3,
}

impl ArchiveOp {
    fn from_code(code: u8) -> Result<Self> {
        match code {
            0 => Ok(Self::Match),
            1 => Ok(Self::Literal),
            2 => Ok(Self::RefAdvance),
            3 => Ok(Self::Whitelist),
            _ => Err(Error::InvalidFormat("unknown archive operation".into())),
        }
    }

    fn consumes_query(self) -> bool {
        matches!(self, Self::Match | Self::Literal)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ArchiveRun {
    op: ArchiveOp,
    len: usize,
    descriptors: Vec<usize>,
}

struct DecodedArchiveChunk {
    first_record: u64,
    next_record: usize,
    ids: Vec<Vec<Vec<u8>>>,
    qualities: Vec<Vec<Option<Vec<u8>>>>,
    sequences: Vec<Vec<Vec<u8>>>,
}

struct ArchiveHeader {
    flags: u16,
    id_mode: IdMode,
    layout: AlignmentLayout,
    chunk_count: u64,
    reference: Vec<u8>,
    indexes: Vec<IndexDescriptor>,
    global_whitelist: Whitelist,
}

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
    global_whitelist: Arc<Whitelist>,
    original_range: Option<RangeInclusive<u64>>,
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
    pending: BTreeMap<u64, Result<DecodedArchiveChunk>>,
    current: Option<DecodedArchiveChunk>,
    input_finished: bool,
    done: bool,
}

trait ReferenceBases: Sync {
    fn len(&self) -> usize;
    fn base(&self, position: usize) -> u8;
}

impl ReferenceBases for Index {
    fn len(&self) -> usize {
        self.reference_len()
    }

    #[inline(always)]
    fn base(&self, position: usize) -> u8 {
        self.reference_base_at(position)
    }
}

impl ReferenceBases for [u8] {
    fn len(&self) -> usize {
        <[u8]>::len(self)
    }

    #[inline(always)]
    fn base(&self, position: usize) -> u8 {
        self[position]
    }
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
    let index = crate::build_genome_index_with_options(
        genome_fasta,
        workspace.path().join("index.seqarc-index"),
        &crate::GenomeIndexBuildOptions {
            temp_dir: Some(workspace.path().to_path_buf()),
        },
    )?;
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
    let index =
        Index::load(index_file, false).map_err(|error| Error::InvalidInput(error.to_string()))?;
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
    let index = crate::build_genome_index_with_options(
        genome_fasta,
        workspace.path().join("index.seqarc-index"),
        &crate::GenomeIndexBuildOptions {
            temp_dir: Some(workspace.path().to_path_buf()),
        },
    )?;
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
    let index =
        Index::load(index_file, false).map_err(|error| Error::InvalidInput(error.to_string()))?;
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
) -> Result<()> {
    if !preserve_ids {
        record.id.clear();
    } else if let Some(prefix_chars) = prefix_chars {
        truncate_read_id(record, prefix_chars)?;
    }
    Ok(())
}

fn index_encode_error(error: IndexReadError) -> Error {
    Error::InvalidInput(error.to_string())
}

fn index_decode_error(error: IndexReadError) -> Error {
    Error::InvalidFormat(error.to_string())
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
    if removed >= record.sequence.len() {
        return Err(Error::InvalidInput(
            "embedded index ranges leave an empty biological read".into(),
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
    record.validate(true)?;
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

fn compress_with_index<I2: AsRef<Path>>(
    index: &Index,
    input_1: impl AsRef<Path>,
    input_2: Option<I2>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let pipeline_started = std::time::Instant::now();
    validate_index_ranges(&options.embedded_index_ranges, layout)?;
    if !options.preserve_ids && options.id_prefix_chars.is_some_and(|length| length != 0) {
        return Err(Error::InvalidInput(
            "identifier prefix cannot be retained when identifiers are disabled".into(),
        ));
    }
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
        "[seqarc] stage=initialize layout={:?} streams={} threads={} mapping_batch_records={} map_task_reads={} chunk_records={} zstd_level={} preserve_ids={} id_prefix_chars={:?} preserve_quality={}",
        layout,
        layout.streams(),
        mapping_threads,
        mapping_batch_records,
        MAP_TASK_READS,
        options.chunk_records.max(1),
        options.zstd_level,
        preserve_ids,
        options.id_prefix_chars,
        options.preserve_quality
    ));
    let temp_dir =
        crate::temporary_directory(options.temp_dir.as_deref(), "seqarc-aligned-spool-")?;
    let spool_path = temp_dir.path().join("aligned-pairs.lz4");
    let mut spool_builder = ExternalChunkBuilder::new(File::create(&spool_path)?, 1)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let mut mask = ReferenceMask::new(index.reference_len());
    let mut id_detector = IdModeDetector::new(preserve_ids, layout);
    let mut whitelist_builders = (0..options.embedded_index_ranges.len())
        .map(|_| WhitelistBuilder::new())
        .collect::<Vec<_>>();
    let mut profile = AlignedCompressionProfile::default();
    let input_1 = input_1.as_ref().to_path_buf();
    let input_2 = input_2.map(|path| path.as_ref().to_path_buf());
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
            let mut first = FastxReader::new(open_fastx_reader(input_1)?, Format::Fastq);
            let mut pairs = Vec::with_capacity(mapping_batch_records);
            let mut preview_identifier = mapper_options.id_prefix_chars.is_some();
            let mut prepare_identifier = |record: &mut ReadRecord| -> Result<()> {
                let original = preview_identifier.then(|| record.id.clone());
                apply_identifier_policy(record, preserve_ids, mapper_options.id_prefix_chars)?;
                if let Some(original) = original {
                    let old_header = format!("@{}", String::from_utf8_lossy(&original));
                    let new_header = format!("@{}", String::from_utf8_lossy(&record.id));
                    stage_log(
                        &producer_progress,
                        format!("[seqarc] id-prefix preview: {old_header:?} -> {new_header:?}"),
                    );
                    preview_identifier = false;
                }
                Ok(())
            };
            if let Some(input_2) = input_2 {
                let mut second = FastxReader::new(open_fastx_reader(input_2)?, Format::Fastq);
                loop {
                    match (first.next(), second.next()) {
                        (None, None) => break,
                        (Some(Ok(mut first)), Some(Ok(mut second))) => {
                            first.validate(true)?;
                            second.validate(true)?;
                            prepare_identifier(&mut first)?;
                            prepare_identifier(&mut second)?;
                            if mapper_options.preserve_quality {
                                apply_quality(&mut first, &mapper_options.quality_mode);
                                apply_quality(&mut second, &mapper_options.quality_mode);
                            }
                            let indexes = split_pair_indexes(
                                &mut first,
                                Some(&mut second),
                                &mapper_options.embedded_index_ranges,
                                mapper_options.preserve_quality,
                            )?;
                            pairs.push(AlignedPair {
                                first,
                                second: Some(second),
                                first_hit: None,
                                second_hit: None,
                                first_selection: AlignmentSelection::default(),
                                second_selection: AlignmentSelection::default(),
                                indexes,
                            });
                        }
                        (None, Some(_)) | (Some(_), None) => {
                            return Err(Error::InvalidInput(
                                "paired streams have different record counts".into(),
                            ));
                        }
                        (Some(Err(error)), _) | (_, Some(Err(error))) => return Err(error),
                    }
                    if pairs.len() == mapping_batch_records {
                        producer_tx
                            .send(Ok(std::mem::take(&mut pairs)))
                            .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
                        pairs = Vec::with_capacity(mapping_batch_records);
                    }
                }
            } else {
                for first in first {
                    let mut first = first?;
                    first.validate(true)?;
                    prepare_identifier(&mut first)?;
                    if mapper_options.preserve_quality {
                        apply_quality(&mut first, &mapper_options.quality_mode);
                    }
                    let indexes = split_pair_indexes(
                        &mut first,
                        None,
                        &mapper_options.embedded_index_ranges,
                        mapper_options.preserve_quality,
                    )?;
                    pairs.push(AlignedPair {
                        first,
                        second: None,
                        first_hit: None,
                        second_hit: None,
                        first_selection: AlignmentSelection::default(),
                        second_selection: AlignmentSelection::default(),
                        indexes,
                    });
                    if pairs.len() == mapping_batch_records {
                        producer_tx
                            .send(Ok(std::mem::take(&mut pairs)))
                            .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
                        pairs = Vec::with_capacity(mapping_batch_records);
                    }
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
            for pair in &pairs {
                if pair.indexes.len() != whitelist_builders.len() {
                    return Err(Error::InvalidInput(
                        "embedded index count changed during ingestion".into(),
                    ));
                }
                for (builder, index) in whitelist_builders.iter_mut().zip(&pair.indexes) {
                    builder.add(&index.sequence).map_err(index_encode_error)?;
                }
            }
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
    let local_whitelists = whitelist_builders
        .into_iter()
        .map(|builder| builder.finish().map_err(index_encode_error))
        .collect::<Result<Vec<_>>>()?;
    let global_whitelist =
        index_read::merge_whitelists(&local_whitelists).map_err(index_encode_error)?;
    let index_descriptors = options
        .embedded_index_ranges
        .iter()
        .cloned()
        .zip(local_whitelists)
        .map(|(range, local_whitelist)| IndexEncodingDescriptor {
            range,
            local_whitelist,
        })
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
    let index_metadata = encode_index_metadata(&index_descriptors, &global_whitelist)?;
    let mut archive = File::create(output)?;
    let header_size = 9 + 4 + 2 + 8 * 3;
    archive.write_all(MAGIC)?;
    write_u32(&mut archive, VERSION)?;
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
    write_u64(&mut archive, reference_block.block.len() as u64)?;
    write_u64(&mut archive, index_metadata.len() as u64)?;
    let chunk_count_pos = archive.stream_position()?;
    write_u64(&mut archive, 0)?;
    debug_assert_eq!(archive.stream_position()?, header_size as u64);
    archive.write_all(&reference_block.block)?;
    archive.write_all(&index_metadata)?;

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
    let global_whitelist_for_workers = &global_whitelist;
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
                            options,
                            id_mode,
                            layout,
                            index_descriptors_for_workers,
                            global_whitelist_for_workers,
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
    archive.seek(SeekFrom::Start(end))?;
    archive.flush()?;
    profile.archive_bytes = end;
    profile.records = replay_pair_count;
    profile.indexes_raw_bytes = 0;
    profile.indexes_compressed_bytes = 0;
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
    for record in records {
        let (_, group) = record?;
        match group {
            AlignedRecordGroup::Paired(first, second) => {
                streams[0].push(first);
                streams[1].push(second);
            }
            AlignedRecordGroup::Single(_) => unreachable!(),
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
    for record in records {
        let (_, group) = record?;
        match group {
            AlignedRecordGroup::Single(record) => stream.push(record),
            AlignedRecordGroup::Paired(_, _) => unreachable!(),
        }
    }
    Ok(ArchiveRecords {
        streams: vec![stream],
    })
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
        let (original_id, group) = record?;
        match group {
            AlignedRecordGroup::Paired(first, second) => {
                sink.write_pair(original_id, first, second)?;
                pairs += 1;
            }
            AlignedRecordGroup::Single(_) => unreachable!(),
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
    let mut records = 0;
    for record in decoded {
        let (original_id, group) = record?;
        sink.write_group(original_id, group)?;
        records += 1;
    }
    sink.finish()?;
    Ok(AlignedDecompressionMetadata {
        layout,
        streams: layout.streams(),
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
        global_whitelist: Arc::new(header.global_whitelist),
        original_range: options.original_range.clone(),
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
        pending: BTreeMap::new(),
        current: None,
        input_finished: false,
        done: false,
    })
}

impl<R: Read> AlignedRecordIter<R> {
    pub fn layout(&self) -> AlignmentLayout {
        self.layout
    }

    pub fn streams(&self) -> usize {
        self.layout.streams()
    }

    pub fn total_chunks(&self) -> u64 {
        self.total_chunks
    }

    pub fn decoded_chunks(&self) -> u64 {
        self.decoded_chunks
    }

    fn resident_chunks(&self) -> usize {
        self.outstanding_jobs + self.pending.len() + usize::from(self.current.is_some())
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
        {
            if self.chunks_remaining == 0 {
                self.finish_input();
                break;
            }
            let chunk_index = self.next_submit_chunk;
            let first_record = self.next_submit_record;
            let encoded = match read_encoded_archive_chunk(&mut self.input, first_record) {
                Ok(encoded) => encoded,
                Err(error) => {
                    self.stop_input_with_error(chunk_index, error);
                    break;
                }
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
            if let (Some(pool), Some(sender)) = (&self.pool, &self.result_tx) {
                self.outstanding_jobs += 1;
                let sender = sender.clone();
                let reference = Arc::clone(&self.reference);
                let id_mode = self.id_mode;
                let layout = self.layout;
                let preserve_quality = self.preserve_quality;
                let indexes = Arc::clone(&self.indexes);
                let global_whitelist = Arc::clone(&self.global_whitelist);
                pool.spawn(move || {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        decode_archive_chunk(
                            encoded,
                            id_mode,
                            layout,
                            preserve_quality,
                            &reference,
                            &indexes,
                            &global_whitelist,
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
                    &self.global_whitelist,
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
            let first = ReadRecord::new(
                std::mem::take(&mut self.ids[0][index]),
                std::mem::take(&mut self.sequences[0][index]),
                self.qualities[0][index].take(),
            );
            let group = if layout == AlignmentLayout::SingleEnd {
                AlignedRecordGroup::Single(first)
            } else {
                AlignedRecordGroup::Paired(
                    first,
                    ReadRecord::new(
                        std::mem::take(&mut self.ids[1][index]),
                        std::mem::take(&mut self.sequences[1][index]),
                        self.qualities[1][index].take(),
                    ),
                )
            };
            return Some((original_id, group));
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
                return Some(Err(error));
            }
            if let Some(current) = &mut self.current {
                if let Some(record) = current.next_group(self.layout, self.original_range.as_ref())
                {
                    return Some(Ok(record));
                }
                self.current = None;
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
                    return Some(Err(error));
                }
            }
        }
    }
}

impl<R: Read> std::iter::FusedIterator for AlignedRecordIter<R> {}

fn read_encoded_archive_chunk<R: Read>(
    input: &mut R,
    first_record: u64,
) -> Result<EncodedArchiveChunk> {
    let record_count = read_u32(input)? as usize;
    let compressed_len = usize_from(read_u64(input)?, "chunk compressed size")?;
    let raw_len = usize_from(read_u64(input)?, "chunk uncompressed size")?;
    let ids_len = usize_from(read_u64(input)?, "chunk ID block size")?;
    let qualities_len = usize_from(read_u64(input)?, "chunk quality section size")?;
    let records_len = usize_from(read_u64(input)?, "chunk record block size")?;
    let mut compressed = vec![0; compressed_len];
    input.read_exact(&mut compressed)?;
    Ok(EncodedArchiveChunk {
        first_record,
        record_count,
        raw_len,
        ids_len,
        qualities_len,
        records_len,
        compressed,
    })
}

fn decode_archive_chunk(
    encoded: EncodedArchiveChunk,
    id_mode: IdMode,
    layout: AlignmentLayout,
    preserve_quality: bool,
    reference: &[u8],
    index_descriptors: &[IndexDescriptor],
    global_whitelist: &Whitelist,
) -> Result<DecodedArchiveChunk> {
    let raw = decode_chunk_sections(
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
    if raw.len() != encoded.raw_len || expected_raw_len != raw.len() {
        return Err(Error::InvalidFormat(
            "chunk block sizes do not match".into(),
        ));
    }
    let ids_end = encoded.ids_len;
    let qualities_end = ids_end + encoded.qualities_len;
    let records_end = qualities_end + encoded.records_len;
    let qualities = decode_qualities(
        &raw[ids_end..qualities_end],
        encoded.record_count,
        preserve_quality,
        layout,
    )?;
    let sequences = decode_records(
        &raw[qualities_end..records_end],
        encoded.record_count,
        reference,
        layout,
        index_descriptors,
        global_whitelist,
    )?;
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
    Ok(DecodedArchiveChunk {
        first_record: encoded.first_record,
        next_record: 0,
        ids: decode_ids(&raw[..ids_end], encoded.record_count, id_mode, layout)?,
        qualities,
        sequences,
    })
}

fn decode_chunk_sections(
    compressed: &[u8],
    ids_len: usize,
    qualities_len: usize,
    records_len: usize,
) -> Result<Vec<u8>> {
    let mut input = Cursor::new(compressed.to_vec());
    let ids_compressed_len = usize_from(read_u64(&mut input)?, "compressed ID section")?;
    let qualities_compressed_len = usize_from(read_u64(&mut input)?, "compressed quality section")?;
    let records_compressed_len = usize_from(read_u64(&mut input)?, "compressed record section")?;
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
    let mut decode_section = |compressed_len: usize, raw_len: usize| -> Result<Vec<u8>> {
        let mut section = vec![0; compressed_len];
        input.read_exact(&mut section)?;
        let raw = zstd::stream::decode_all(Cursor::new(section))?;
        if raw.len() != raw_len {
            return Err(Error::InvalidFormat("chunk section length mismatch".into()));
        }
        Ok(raw)
    };
    let ids = decode_section(ids_compressed_len, ids_len)?;
    let qualities = decode_section(qualities_compressed_len, qualities_len)?;
    let records = decode_section(records_compressed_len, records_len)?;
    let mut raw = Vec::with_capacity(ids_len + qualities_len + records_len);
    raw.extend(ids);
    raw.extend(qualities);
    raw.extend(records);
    Ok(raw)
}

fn map_pairs(
    index: &Index,
    pairs: &mut [AlignedPair],
    mapping_pool: &rayon::ThreadPool,
    progress: &ProgressBar,
    progress_stats: &mut MappingProgressStats,
) -> Result<()> {
    let first_hits = map_mate(index, pairs, 0, mapping_pool)?;
    progress_stats.observe(0, &first_hits, progress);
    for (pair, first) in pairs.iter_mut().zip(first_hits) {
        pair.first_hit = first.hit;
        pair.first_selection = first.selection;
    }
    if pairs.first().is_some_and(|pair| pair.second.is_some()) {
        let second_hits = map_mate(index, pairs, 1, mapping_pool)?;
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

fn map_mate(
    index: &Index,
    pairs: &[AlignedPair],
    mate: usize,
    mapping_pool: &rayon::ThreadPool,
) -> Result<Vec<SelectedAlignment>> {
    let opts = Opts::new().set_paired(false);
    let aligner = Aligner::new(index, &opts);
    mapping_pool.install(|| {
        use rayon::prelude::*;
        pairs
            .par_chunks(MAP_TASK_READS)
            .map_init(ThreadBuf::new, |buffer, batch| {
                let records = batch.iter().map(|pair| pair.record(mate));
                align_batch(&aligner, buffer, index, records)
            })
            .collect::<Result<Vec<_>>>()
            .map(|batches| batches.into_iter().flatten().collect())
    })
}

impl ReferenceMask {
    fn new(length: usize) -> Self {
        Self {
            bits: vec![0; length.div_ceil(64)],
            length,
        }
    }

    fn mark_hit(&mut self, hit: &GenomeAlignment) -> Result<()> {
        let mut position = nonnegative_usize(hit.start, "alignment start")?;
        for op in &hit.cigar {
            let length = op.len as usize;
            match op.kind {
                CigarKind::Match | CigarKind::Eq | CigarKind::Diff => {
                    let end = position.checked_add(length).ok_or_else(|| {
                        Error::InvalidInput("reference mask position overflow".into())
                    })?;
                    self.mark_range(position, end)?;
                    position = end;
                }
                CigarKind::Del | CigarKind::RefSkip => {
                    position = position.checked_add(length).ok_or_else(|| {
                        Error::InvalidInput("reference mask position overflow".into())
                    })?;
                    if position > self.length {
                        return Err(Error::InvalidFormat(
                            "reference mask range out of bounds".into(),
                        ));
                    }
                }
                CigarKind::Ins | CigarKind::SoftClip | CigarKind::HardClip | CigarKind::Pad => {}
            }
        }
        Ok(())
    }

    fn mark_range(&mut self, start: usize, end: usize) -> Result<()> {
        if start > end || end > self.length {
            return Err(Error::InvalidFormat(
                "reference mask range out of bounds".into(),
            ));
        }
        for position in start..end {
            self.bits[position / 64] |= 1u64 << (position % 64);
        }
        Ok(())
    }

    fn contains(&self, position: usize) -> bool {
        position < self.length && (self.bits[position / 64] & (1u64 << (position % 64))) != 0
    }
}

fn compact_reference<R: ReferenceBases + ?Sized>(
    reference: &R,
    mask: &ReferenceMask,
) -> Result<(Vec<u8>, Vec<KeptRun>)> {
    let mut compact = Vec::new();
    let mut runs: Vec<KeptRun> = Vec::new();
    for position in 0..reference.len() {
        if !mask.contains(position) {
            continue;
        }
        let base = reference.base(position);
        if !matches!(base, b'A' | b'C' | b'G' | b'T') {
            return Err(Error::InvalidInput(
                "reference contains a noncanonical base".into(),
            ));
        }
        let compact_position = compact.len();
        compact.push(base);
        if let Some(run) = runs.last_mut() {
            if run.original_start + run.length == position
                && run.compact_start + run.length == compact_position
            {
                run.length += 1;
                continue;
            }
        }
        runs.push(KeptRun {
            original_start: position,
            compact_start: compact_position,
            length: 1,
        });
    }
    Ok((compact, runs))
}

fn compact_boundary(runs: &[KeptRun], original: usize, reference_len: usize) -> Result<usize> {
    if original > reference_len {
        return Err(Error::InvalidFormat(
            "alignment position exceeds reference".into(),
        ));
    }
    let index = runs.partition_point(|run| run.original_start <= original);
    Ok(index
        .checked_sub(1)
        .and_then(|index| runs.get(index))
        .map_or(0, |run| {
            run.compact_start + (original - run.original_start).min(run.length)
        }))
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

fn encode_chunk(
    reference: &[u8],
    pairs: &[AlignedPair],
    options: &crate::AlignmentCompressionOptions,
    id_mode: IdMode,
    layout: AlignmentLayout,
    index_descriptors: &[IndexEncodingDescriptor],
    global_whitelist: &Whitelist,
) -> Result<EncodedChunk> {
    let mut ids = Vec::new();
    write_id_stream(&mut ids, pairs, id_mode, layout)?;
    let mut qualities = Vec::new();
    write_quality_stream(&mut qualities, pairs, index_descriptors, options, layout)?;
    let records = write_record_stream(
        reference,
        pairs,
        layout,
        index_descriptors,
        global_whitelist,
    )?;
    let ids_compressed = zstd::stream::encode_all(Cursor::new(&ids), options.zstd_level)?;
    let qualities_compressed =
        zstd::stream::encode_all(Cursor::new(&qualities), options.zstd_level)?;
    let records_compressed = zstd::stream::encode_all(Cursor::new(&records), options.zstd_level)?;
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
        raw_len: (ids.len() + qualities.len() + records.len()) as u64,
        ids_len: ids.len() as u64,
        qualities_len: qualities.len() as u64,
        records_len: records.len() as u64,
        ids_compressed_len: ids_compressed.len() as u64,
        qualities_compressed_len: qualities_compressed.len() as u64,
        records_compressed_len: records_compressed.len() as u64,
        compressed,
    })
}

fn write_encoded_chunk(archive: &mut File, chunk: EncodedChunk) -> Result<()> {
    write_u32(archive, chunk.pairs)?;
    write_u64(archive, chunk.compressed.len() as u64)?;
    write_u64(archive, chunk.raw_len)?;
    write_u64(archive, chunk.ids_len)?;
    write_u64(archive, chunk.qualities_len)?;
    write_u64(archive, chunk.records_len)?;
    archive.write_all(&chunk.compressed)?;
    Ok(())
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
    descriptors: &[IndexEncodingDescriptor],
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<()> {
    for mate in 0..layout.streams() {
        let stream = if options.preserve_quality {
            let rows = pairs
                .iter()
                .map(|pair| reconstruct_quality_row(pair, mate, descriptors))
                .collect::<Result<Vec<_>>>()?;
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
    descriptors: &[IndexEncodingDescriptor],
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

fn write_record_stream<R: ReferenceBases + ?Sized>(
    reference: &R,
    pairs: &[AlignedPair],
    layout: AlignmentLayout,
    descriptors: &[IndexEncodingDescriptor],
    global_whitelist: &Whitelist,
) -> Result<Vec<u8>> {
    if global_whitelist.len() > u16::MAX as usize + 1 {
        return Err(Error::InvalidInput(
            "global whitelist exceeds u16 ordinal capacity".into(),
        ));
    }
    let mut output = Vec::new();
    let mut starts = [0i64; 2];
    for mate in 0..layout.streams() {
        let bitmap_len = pairs.len().div_ceil(8);
        let mut mapped_bitmap = vec![0; bitmap_len];
        let mut reverse_bitmap = vec![0; bitmap_len];
        let mut positions = Vec::new();
        let mut programs = Vec::new();
        let mut raw_bits = BitWriter::default();
        let mut mismatch_bits = BitWriter::default();
        let mut ordinals = Vec::new();
        let mut mismatch_details = Vec::new();
        let mut literals = Vec::new();
        for (index, pair) in pairs.iter().enumerate() {
            let record = pair.record(mate);
            let hit = pair.hit(mate);
            if let Some(hit) = hit {
                mapped_bitmap[index / 8] |= 1 << (index % 8);
                let encoded_reverse = if mate == 1 {
                    pair.first_hit
                        .as_ref()
                        .map_or(hit.is_reverse, |first| hit.is_reverse ^ first.is_reverse)
                } else {
                    hit.is_reverse
                };
                if encoded_reverse {
                    reverse_bitmap[index / 8] |= 1 << (index % 8);
                }
                let base = if mate == 1 {
                    pair.first_hit
                        .as_ref()
                        .map_or(starts[mate], |first| first.start)
                } else {
                    starts[mate]
                };
                let delta = hit
                    .start
                    .checked_sub(base)
                    .ok_or_else(|| Error::InvalidInput("alignment offset overflow".into()))?;
                write_varint(&mut positions, zigzag_encode(delta))?;
                starts[mate] = hit.start;
                let core = archive_normalize(effective_cigar(hit, record.sequence.len())?)?;
                let runs = insert_whitelist_runs(
                    core,
                    descriptor_boundaries(
                        descriptors,
                        mate,
                        record.sequence.len(),
                        hit.is_reverse,
                    )?,
                    record.sequence.len(),
                )?;
                encode_record_columns(
                    &mut programs,
                    &mut raw_bits,
                    &mut mismatch_bits,
                    &mut ordinals,
                    &mut mismatch_details,
                    &mut literals,
                    &record.sequence,
                    Some((hit, reference)),
                    &runs,
                    &pair.indexes,
                    descriptors,
                    global_whitelist,
                )?;
            } else {
                let runs = insert_whitelist_runs(
                    vec![ArchiveRun {
                        op: ArchiveOp::Literal,
                        len: record.sequence.len(),
                        descriptors: Vec::new(),
                    }],
                    descriptor_boundaries(descriptors, mate, record.sequence.len(), false)?,
                    record.sequence.len(),
                )?;
                encode_record_columns(
                    &mut programs,
                    &mut raw_bits,
                    &mut mismatch_bits,
                    &mut ordinals,
                    &mut mismatch_details,
                    &mut literals,
                    &record.sequence,
                    None::<(&GenomeAlignment, &R)>,
                    &runs,
                    &pair.indexes,
                    descriptors,
                    global_whitelist,
                )?;
            }
        }
        output.extend(mapped_bitmap);
        output.extend(reverse_bitmap);
        write_u64(&mut output, positions.len() as u64)?;
        output.extend(positions);
        write_u64(&mut output, programs.len() as u64)?;
        write_u64(&mut output, raw_bits.len() as u64)?;
        write_u64(&mut output, raw_bits.bytes.len() as u64)?;
        write_u64(&mut output, mismatch_bits.len() as u64)?;
        write_u64(&mut output, mismatch_bits.bytes.len() as u64)?;
        write_u64(&mut output, ordinals.len() as u64)?;
        write_u64(&mut output, mismatch_details.len() as u64)?;
        write_u64(&mut output, literals.len() as u64)?;
        output.extend(programs);
        output.extend(raw_bits.finish());
        output.extend(mismatch_bits.finish());
        output.extend(ordinals);
        output.extend(mismatch_details);
        output.extend(literals);
    }
    Ok(output)
}

fn encode_record_columns<R: ReferenceBases + ?Sized>(
    programs: &mut Vec<u8>,
    raw_bits: &mut BitWriter,
    mismatch_bits: &mut BitWriter,
    ordinals: &mut Vec<u8>,
    mismatch_details: &mut Vec<u8>,
    literals: &mut Vec<u8>,
    sequence: &[u8],
    mapped: Option<(&GenomeAlignment, &R)>,
    runs: &[ArchiveRun],
    indexes: &[SpoolIndexRead],
    descriptors: &[IndexEncodingDescriptor],
    global_whitelist: &Whitelist,
) -> Result<()> {
    let mut classifications = Vec::new();
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
        if run.op == ArchiveOp::Whitelist {
            if run.len != 0 || run.descriptors.is_empty() {
                return Err(Error::InvalidInput(
                    "whitelist operation has invalid descriptors".into(),
                ));
            }
            continue;
        }
        if run.len == 0 || !run.descriptors.is_empty() {
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
            ArchiveOp::Whitelist => unreachable!(),
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
            ArchiveOp::Whitelist => {}
        }
    }

    for run in runs.iter().filter(|run| run.op == ArchiveOp::Whitelist) {
        for &descriptor_index in &run.descriptors {
            let descriptor = descriptors.get(descriptor_index).ok_or_else(|| {
                Error::InvalidInput("whitelist descriptor is out of bounds".into())
            })?;
            let index = indexes
                .get(descriptor_index)
                .ok_or_else(|| Error::InvalidInput("embedded index fragment is missing".into()))?;
            if index.sequence.len() != descriptor.range.len() {
                return Err(Error::InvalidInput(
                    "embedded index sequence length does not match its range".into(),
                ));
            }
            let classification = index_read::classify_with_global(
                &index.sequence,
                &descriptor.local_whitelist,
                global_whitelist,
            )
            .map_err(index_encode_error)?;
            let classification = if descriptor.range.len() > 64
                && matches!(classification, IndexReadEncoding::OneMismatch { .. })
            {
                IndexReadEncoding::CanonicalRaw
            } else {
                classification
            };
            classifications.push((descriptor_index, classification));
        }
    }
    for (_, classification) in &classifications {
        raw_bits.write(
            usize::from(matches!(
                classification,
                IndexReadEncoding::CanonicalRaw | IndexReadEncoding::ExceptionalRaw
            )),
            1,
        )?;
    }
    for &(descriptor_index, ref classification) in &classifications {
        let descriptor = &descriptors[descriptor_index];
        match *classification {
            IndexReadEncoding::Exact { global_ordinal } => {
                let ordinal = u16::try_from(global_ordinal)
                    .map_err(|_| Error::InvalidInput("whitelist ordinal exceeds u16".into()))?;
                mismatch_bits.write(0, 1)?;
                ordinals.extend_from_slice(&ordinal.to_le_bytes());
            }
            IndexReadEncoding::OneMismatch {
                global_ordinal,
                position,
                original,
            } => {
                let selected = global_whitelist.sequence(global_ordinal).ok_or_else(|| {
                    Error::InvalidInput("whitelist ordinal is out of bounds".into())
                })?;
                if descriptor.range.len() > 64
                    || selected.len() != descriptor.range.len()
                    || position >= selected.len()
                {
                    return Err(Error::InvalidInput(
                        "whitelist mismatch does not match descriptor".into(),
                    ));
                }
                let original = canonical_bits(original)?;
                if original == canonical_bits(selected[position])? {
                    return Err(Error::InvalidInput(
                        "whitelist replacement does not differ".into(),
                    ));
                }
                let ordinal = u16::try_from(global_ordinal)
                    .map_err(|_| Error::InvalidInput("whitelist ordinal exceeds u16".into()))?;
                mismatch_bits.write(1, 1)?;
                ordinals.extend_from_slice(&ordinal.to_le_bytes());
                mismatch_details.push(((position as u8) << 2) | original);
            }
            IndexReadEncoding::CanonicalRaw | IndexReadEncoding::ExceptionalRaw => {}
        }
    }

    let mut classification_index = 0;
    let mut query_pos = 0;
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
            ArchiveOp::Whitelist => {
                for &descriptor_index in &run.descriptors {
                    let (_, classification) =
                        classifications.get(classification_index).ok_or_else(|| {
                            Error::InvalidInput("missing whitelist classification".into())
                        })?;
                    classification_index += 1;
                    if matches!(
                        classification,
                        IndexReadEncoding::CanonicalRaw | IndexReadEncoding::ExceptionalRaw
                    ) {
                        let mut fragment = indexes[descriptor_index].sequence.clone();
                        if reverse {
                            reverse_complement(&mut fragment);
                        }
                        literals.extend(fragment);
                    }
                }
            }
        }
    }
    if classification_index != classifications.len() {
        return Err(Error::InvalidInput(
            "literal source consumption mismatch".into(),
        ));
    }
    Ok(())
}

fn alignment_candidate_core_cost<R: ReferenceBases + ?Sized>(
    sequence: &[u8],
    mapped: Option<(&GenomeAlignment, &R)>,
    runs: &[ArchiveRun],
) -> Result<usize> {
    let mut programs = Vec::new();
    let mut raw_bits = BitWriter::default();
    let mut mismatch_bits = BitWriter::default();
    let mut ordinals = Vec::new();
    let mut mismatch_details = Vec::new();
    let mut literals = Vec::new();
    let empty_whitelist = index_read::merge_whitelists(&[]).map_err(index_encode_error)?;
    encode_record_columns(
        &mut programs,
        &mut raw_bits,
        &mut mismatch_bits,
        &mut ordinals,
        &mut mismatch_details,
        &mut literals,
        sequence,
        mapped,
        runs,
        &[],
        &[],
        &empty_whitelist,
    )?;
    Ok(programs.len() + literals.len())
}

fn align_batch<'a, I, R>(
    aligner: &Aligner<'_>,
    buffer: &mut ThreadBuf,
    reference: &R,
    records: I,
) -> Result<Vec<SelectedAlignment>>
where
    I: IntoIterator<Item = &'a ReadRecord>,
    R: ReferenceBases + ?Sized,
{
    let records = records.into_iter().collect::<Vec<_>>();
    let names = (0..records.len())
        .map(|index| format!("read-{index}"))
        .collect::<Vec<_>>();
    let sequences = records
        .iter()
        .map(|record| record.sequence.as_slice())
        .collect::<Vec<_>>();
    aligner
        .map_many(buffer, &names, &sequences)
        .map_err(|error| Error::InvalidInput(error.to_string()))
        .and_then(|hits| {
            records
                .into_iter()
                .zip(hits)
                .map(|(record, hits)| select_alignment(record, hits, reference))
                .collect()
        })
}

fn select_alignment<R: ReferenceBases + ?Sized>(
    record: &ReadRecord,
    hits: Vec<minibwa::Hit>,
    reference: &R,
) -> Result<SelectedAlignment> {
    let raw_cost = alignment_candidate_core_cost(
        &record.sequence,
        None::<(&GenomeAlignment, &R)>,
        &[ArchiveRun {
            op: ArchiveOp::Literal,
            len: record.sequence.len(),
            descriptors: Vec::new(),
        }],
    )?;
    let candidate_hits = u32::try_from(hits.len()).unwrap_or(u32::MAX);
    let mut best: Option<(usize, i32, usize, GenomeAlignment)> = None;
    for (index, hit) in hits.into_iter().enumerate() {
        let alignment = GenomeAlignment {
            start: hit.ref_start,
            query_start: hit.query_start,
            query_end: hit.query_end,
            cigar: hit.cigar,
            is_reverse: hit.strand == Strand::Reverse,
        };
        let encoded_result = effective_cigar(&alignment, record.sequence.len())
            .and_then(archive_normalize)
            .and_then(|runs| {
                alignment_candidate_core_cost(
                    &record.sequence,
                    Some((&alignment, reference)),
                    &runs,
                )
            });
        if encoded_result.is_err() {
            continue;
        }
        let candidate = (encoded_result.unwrap(), hit.score, index, alignment);
        if best.as_ref().is_none_or(|current| {
            candidate.0 < current.0
                || (candidate.0 == current.0 && candidate.1 > current.1)
                || (candidate.0 == current.0 && candidate.1 == current.1 && candidate.2 < current.2)
        }) {
            best = Some(candidate);
        }
    }
    let use_alignment = best
        .as_ref()
        .is_some_and(|candidate| candidate.0 < raw_cost);
    let selected_non_first =
        use_alignment && best.as_ref().is_some_and(|candidate| candidate.2 != 0);
    let rejected_alignment = candidate_hits != 0 && !use_alignment;
    Ok(SelectedAlignment {
        hit: if use_alignment {
            best.map(|candidate| candidate.3)
        } else {
            None
        },
        selection: AlignmentSelection {
            candidate_hits,
            selected_non_first,
            rejected_alignment,
        },
    })
}

#[cfg(test)]
fn encode_record_body<R: ReferenceBases + ?Sized>(
    output: &mut Vec<u8>,
    record: &ReadRecord,
    hit: &GenomeAlignment,
    reference: &R,
) -> Result<()> {
    let runs = archive_normalize(effective_cigar(hit, record.sequence.len())?)?;
    let empty_whitelist = index_read::merge_whitelists(&[]).map_err(index_encode_error)?;
    encode_operation_body(
        output,
        &record.sequence,
        Some((hit, reference)),
        &runs,
        &[],
        &[],
        &empty_whitelist,
    )
}

enum WhitelistFragment {
    Raw(usize),
    Bases(Vec<u8>),
}

struct RecordProgram {
    runs: Vec<ArchiveRun>,
    matched: Vec<Option<Vec<u8>>>,
    reduced_len: usize,
    reverse: bool,
    whitelist_positions: Vec<usize>,
    group_sizes: Vec<usize>,
    descriptor_groups: Vec<Vec<usize>>,
    fragments: Vec<(usize, WhitelistFragment)>,
}

fn decode_record_program(
    input: &mut Cursor<Vec<u8>>,
    mapped: Option<(&[u8], usize, bool)>,
) -> Result<RecordProgram> {
    let count = usize_from(read_varint(input)?, "archive operation count")?;
    let packed_len = count.div_ceil(4);
    let mut packed = vec![0; packed_len];
    input.read_exact(&mut packed)?;
    if count % 4 != 0 {
        let padding_bits = (4 - count % 4) * 2;
        if packed.last().copied().unwrap_or(0) & ((1 << padding_bits) - 1) != 0 {
            return Err(Error::InvalidFormat(
                "nonzero archive operation padding bits".into(),
            ));
        }
    }
    let mut runs = Vec::with_capacity(count);
    for index in 0..count {
        let op = ArchiveOp::from_code((packed[index / 4] >> (6 - index % 4 * 2)) & 3)?;
        let len = if op == ArchiveOp::Whitelist {
            0
        } else {
            let len = usize_from(read_varint(input)?, "archive operation length")?;
            if len == 0 {
                return Err(Error::InvalidFormat("zero-length archive operation".into()));
            }
            len
        };
        runs.push(ArchiveRun {
            op,
            len,
            descriptors: Vec::new(),
        });
    }
    if mapped.is_none()
        && runs
            .iter()
            .any(|run| matches!(run.op, ArchiveOp::Match | ArchiveOp::RefAdvance))
    {
        return Err(Error::InvalidFormat(
            "unmapped body contains mapped operations".into(),
        ));
    }
    let (mut ref_pos, reverse) = mapped.map_or((0, false), |(_, start, reverse)| (start, reverse));
    let mut reduced_len = 0usize;
    let mut matched = (0..runs.len()).map(|_| None).collect::<Vec<_>>();
    let mut whitelist_positions = Vec::new();
    let mut group_sizes = Vec::new();
    for (run_index, run) in runs.iter().enumerate() {
        match run.op {
            ArchiveOp::Match => {
                let reference = mapped.unwrap().0;
                let end = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if end > reference.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds reference".into(),
                    ));
                }
                let mut bases = reference[ref_pos..end].to_vec();
                let mismatch_count = usize_from(read_varint(input)?, "mismatch count")?;
                let mut previous: usize = 0;
                for mismatch in 0..mismatch_count {
                    let position = previous
                        .checked_add(usize_from(read_varint(input)?, "mismatch position")?)
                        .ok_or_else(|| Error::InvalidFormat("mismatch position overflow".into()))?;
                    if position >= run.len || (mismatch != 0 && position <= previous) {
                        return Err(Error::InvalidFormat(
                            "mismatch position out of bounds".into(),
                        ));
                    }
                    bases[position] = read_u8(input)?;
                    previous = position;
                }
                ref_pos = end;
                reduced_len = reduced_len
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reduced query length overflow".into()))?;
                matched[run_index] = Some(bases);
            }
            ArchiveOp::Literal => {
                reduced_len = reduced_len
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reduced query length overflow".into()))?;
            }
            ArchiveOp::RefAdvance => {
                ref_pos = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if ref_pos > mapped.unwrap().0.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds reference".into(),
                    ));
                }
            }
            ArchiveOp::Whitelist => {
                whitelist_positions.push(reduced_len);
                group_sizes.push(0);
            }
        }
    }
    if let Some((reference, _, _)) = mapped {
        let _ = reference;
    }
    Ok(RecordProgram {
        runs,
        matched,
        reduced_len,
        reverse,
        whitelist_positions,
        group_sizes,
        descriptor_groups: Vec::new(),
        fragments: Vec::new(),
    })
}

fn expected_whitelist_groups(
    reduced_len: usize,
    descriptors: &[IndexDescriptor],
    mate: usize,
    reverse: bool,
) -> Result<Vec<(usize, Vec<usize>)>> {
    let mut selected = descriptors
        .iter()
        .enumerate()
        .filter(|(_, descriptor)| descriptor.range.mate.stream() == mate)
        .collect::<Vec<_>>();
    selected.sort_by_key(|(_, descriptor)| descriptor.range.start);
    let mut prior_removed = 0usize;
    let mut groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for (index, descriptor) in selected {
        let boundary = descriptor
            .range
            .start
            .checked_sub(prior_removed)
            .ok_or_else(|| Error::InvalidFormat("embedded index boundary underflow".into()))?;
        if boundary > reduced_len {
            return Err(Error::InvalidFormat(
                "embedded index boundary exceeds reduced read".into(),
            ));
        }
        if let Some((_, group)) = groups
            .last_mut()
            .filter(|(position, _)| *position == boundary)
        {
            group.push(index);
        } else {
            groups.push((boundary, vec![index]));
        }
        prior_removed = prior_removed
            .checked_add(descriptor.range.len())
            .ok_or_else(|| Error::InvalidFormat("embedded index length overflow".into()))?;
    }
    if reverse {
        groups.reverse();
        for (position, group) in &mut groups {
            *position = reduced_len - *position;
            group.reverse();
        }
    }
    Ok(groups)
}

fn decode_record_literals(program: &RecordProgram, input: &mut Cursor<Vec<u8>>) -> Result<Vec<u8>> {
    enum Piece {
        Bases(Vec<u8>),
        Literal(usize),
    }
    let mut pieces = Vec::new();
    let mut group_index = 0;
    for (index, run) in program.runs.iter().enumerate() {
        match run.op {
            ArchiveOp::Match => {
                pieces.push(Piece::Bases(program.matched[index].clone().ok_or_else(
                    || Error::InvalidFormat("missing matched sequence".into()),
                )?))
            }
            ArchiveOp::Literal => pieces.push(Piece::Literal(run.len)),
            ArchiveOp::RefAdvance => {}
            ArchiveOp::Whitelist => {
                for descriptor_index in &program.descriptor_groups[group_index] {
                    let (_, fragment) = program
                        .fragments
                        .iter()
                        .find(|(index, _)| index == descriptor_index)
                        .ok_or_else(|| Error::InvalidFormat("missing whitelist fragment".into()))?;
                    match fragment {
                        WhitelistFragment::Raw(length) => pieces.push(Piece::Literal(*length)),
                        WhitelistFragment::Bases(sequence) => {
                            pieces.push(Piece::Bases(sequence.clone()))
                        }
                    }
                }
                group_index += 1;
            }
        }
    }
    if program.fragments.len()
        != program
            .descriptor_groups
            .iter()
            .map(Vec::len)
            .sum::<usize>()
    {
        return Err(Error::InvalidFormat("unused whitelist fragment".into()));
    }
    let mut sequence = Vec::new();
    let mut index = 0;
    while index < pieces.len() {
        match &pieces[index] {
            Piece::Bases(bases) => {
                sequence.extend_from_slice(bases);
                index += 1;
            }
            Piece::Literal(_) => {
                let start = index;
                let mut length = 0usize;
                while index < pieces.len() {
                    let Piece::Literal(value) = pieces[index] else {
                        break;
                    };
                    length = length.checked_add(value).ok_or_else(|| {
                        Error::InvalidFormat("literal run length overflow".into())
                    })?;
                    index += 1;
                }
                let literal = take_literal_bytes(input, length)?;
                let mut offset = 0usize;
                for piece in &pieces[start..index] {
                    let Piece::Literal(source_length) = *piece else {
                        unreachable!();
                    };
                    let end = offset.checked_add(source_length).ok_or_else(|| {
                        Error::InvalidFormat("literal source length overflow".into())
                    })?;
                    sequence.extend_from_slice(&literal[offset..end]);
                    offset = end;
                }
                debug_assert_eq!(offset, literal.len());
            }
        }
    }
    if program.reverse {
        reverse_complement(&mut sequence);
    }
    Ok(sequence)
}

fn decode_records(
    input: &[u8],
    count: usize,
    reference: &[u8],
    layout: AlignmentLayout,
    descriptors: &[IndexDescriptor],
    global_whitelist: &Whitelist,
) -> Result<Vec<Vec<Vec<u8>>>> {
    let mut cursor = Cursor::new(input.to_vec());
    let mut starts = [0i64; 2];
    let mut first_starts = vec![None; count];
    let mut first_reverses = vec![None; count];
    let mut all = (0..layout.streams())
        .map(|_| Vec::with_capacity(count))
        .collect::<Vec<_>>();
    for mate in 0..layout.streams() {
        let bitmap_len = count.div_ceil(8);
        let mut mapped_bitmap = vec![0; bitmap_len];
        let mut reverse_bitmap = vec![0; bitmap_len];
        cursor.read_exact(&mut mapped_bitmap)?;
        cursor.read_exact(&mut reverse_bitmap)?;
        if count % 8 != 0 {
            let padding = (8 - count % 8) as u8;
            let mask = !((1u8 << (8 - padding)) - 1);
            if mapped_bitmap.last().copied().unwrap_or(0) & mask != 0
                || reverse_bitmap.last().copied().unwrap_or(0) & mask != 0
            {
                return Err(Error::InvalidFormat(
                    "nonzero alignment bitmap padding".into(),
                ));
            }
        }
        let positions_len = usize_from(read_u64(&mut cursor)?, "alignment positions")?;
        let mut positions = vec![0; positions_len];
        cursor.read_exact(&mut positions)?;
        let mut positions = Cursor::new(positions);
        let mut mapped_values = Vec::with_capacity(count);
        for ordinal in 0..count {
            let mapped = mapped_bitmap[ordinal / 8] & (1 << (ordinal % 8)) != 0;
            let encoded_reverse = reverse_bitmap[ordinal / 8] & (1 << (ordinal % 8)) != 0;
            if mapped {
                let delta = zigzag_decode(read_varint(&mut positions)?)?;
                let base = if mate == 1 {
                    first_starts[ordinal].unwrap_or(starts[mate])
                } else {
                    starts[mate]
                };
                let start = base
                    .checked_add(delta)
                    .ok_or_else(|| Error::InvalidFormat("alignment offset overflow".into()))?;
                starts[mate] = start;
                let reverse = if mate == 1 {
                    first_reverses[ordinal].map_or(encoded_reverse, |first| encoded_reverse ^ first)
                } else {
                    encoded_reverse
                };
                if mate == 0 {
                    first_starts[ordinal] = Some(start);
                    first_reverses[ordinal] = Some(reverse);
                }
                mapped_values.push(Some((start, reverse)));
            } else {
                mapped_values.push(None);
            }
        }
        if positions.position() != positions.get_ref().len() as u64 {
            return Err(Error::InvalidFormat(
                "trailing alignment position bytes".into(),
            ));
        }

        let programs_len = usize_from(read_u64(&mut cursor)?, "record programs")?;
        let raw_bit_count = usize_from(read_u64(&mut cursor)?, "raw bitmap bit count")?;
        let raw_byte_len = usize_from(read_u64(&mut cursor)?, "raw bitmap byte length")?;
        let mismatch_bit_count = usize_from(read_u64(&mut cursor)?, "mismatch bitmap bit count")?;
        let mismatch_byte_len = usize_from(read_u64(&mut cursor)?, "mismatch bitmap byte length")?;
        let ordinal_byte_len = usize_from(read_u64(&mut cursor)?, "ordinal column byte length")?;
        let mismatch_detail_byte_len =
            usize_from(read_u64(&mut cursor)?, "mismatch detail byte length")?;
        let literals_len = usize_from(read_u64(&mut cursor)?, "record literals")?;
        if raw_byte_len != raw_bit_count.div_ceil(8)
            || mismatch_byte_len != mismatch_bit_count.div_ceil(8)
        {
            return Err(Error::InvalidFormat("record bitmap length mismatch".into()));
        }
        let payload_len = programs_len
            .checked_add(raw_byte_len)
            .and_then(|length| length.checked_add(mismatch_byte_len))
            .and_then(|length| length.checked_add(ordinal_byte_len))
            .and_then(|length| length.checked_add(mismatch_detail_byte_len))
            .and_then(|length| length.checked_add(literals_len))
            .ok_or_else(|| Error::InvalidFormat("record section length overflow".into()))?;
        let payload_end = (cursor.position() as usize)
            .checked_add(payload_len)
            .ok_or_else(|| Error::InvalidFormat("record section length overflow".into()))?;
        if payload_end > cursor.get_ref().len() {
            return Err(Error::InvalidFormat("truncated record section".into()));
        }
        let take = |cursor: &mut Cursor<Vec<u8>>, length: usize| -> Result<Vec<u8>> {
            let start = cursor.position() as usize;
            let end = start
                .checked_add(length)
                .ok_or_else(|| Error::InvalidFormat("record section length overflow".into()))?;
            if end > cursor.get_ref().len() {
                return Err(Error::InvalidFormat("truncated record section".into()));
            }
            let value = cursor.get_ref()[start..end].to_vec();
            cursor.set_position(end as u64);
            Ok(value)
        };
        let programs = take(&mut cursor, programs_len)?;
        let raw_bytes = take(&mut cursor, raw_byte_len)?;
        let mismatch_bytes = take(&mut cursor, mismatch_byte_len)?;
        let ordinal_bytes = take(&mut cursor, ordinal_byte_len)?;
        let mismatch_detail_bytes = take(&mut cursor, mismatch_detail_byte_len)?;
        let literal_bytes = take(&mut cursor, literals_len)?;

        let mut program_cursor = Cursor::new(programs);
        let mut decoded = Vec::with_capacity(count);
        for mapped in &mapped_values {
            let mapping = match mapped {
                Some((start, reverse)) => Some((
                    reference,
                    nonnegative_usize(*start, "alignment start")?,
                    *reverse,
                )),
                None => None,
            };
            decoded.push(decode_record_program(&mut program_cursor, mapping)?);
        }
        if program_cursor.position() != program_cursor.get_ref().len() as u64 {
            return Err(Error::InvalidFormat("trailing record program bytes".into()));
        }
        let expected_descriptors = descriptors
            .iter()
            .enumerate()
            .filter(|(_, descriptor)| descriptor.range.mate.stream() == mate)
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        let expected_occurrences = count
            .checked_mul(expected_descriptors.len())
            .ok_or_else(|| Error::InvalidFormat("whitelist occurrence count overflow".into()))?;
        if raw_bit_count != expected_occurrences {
            return Err(Error::InvalidFormat("raw bitmap count mismatch".into()));
        }
        let mut raw = BitReader::new(&raw_bytes, raw_bit_count)?;
        let mut mismatch = BitReader::new(&mismatch_bytes, mismatch_bit_count)?;
        let raw_count = raw_bytes
            .iter()
            .map(|byte| byte.count_ones() as usize)
            .sum::<usize>();
        let nonraw_count = expected_occurrences
            .checked_sub(raw_count)
            .ok_or_else(|| Error::InvalidFormat("raw bitmap count exceeds descriptors".into()))?;
        if mismatch_bit_count != nonraw_count {
            return Err(Error::InvalidFormat(
                "mismatch bitmap count mismatch".into(),
            ));
        }
        let expected_ordinal_bytes = nonraw_count
            .checked_mul(2)
            .ok_or_else(|| Error::InvalidFormat("ordinal column length overflow".into()))?;
        if ordinal_byte_len != expected_ordinal_bytes {
            return Err(Error::InvalidFormat(
                "ordinal column length mismatch".into(),
            ));
        }
        let mismatch_count = mismatch_bytes
            .iter()
            .map(|byte| byte.count_ones() as usize)
            .sum::<usize>();
        if mismatch_detail_byte_len != mismatch_count {
            return Err(Error::InvalidFormat(
                "mismatch detail column length mismatch".into(),
            ));
        }
        let mut ordinal_offset = 0usize;
        let mut detail_offset = 0usize;
        let mut literal_cursor = Cursor::new(literal_bytes);
        for program in &mut decoded {
            let expected_groups =
                expected_whitelist_groups(program.reduced_len, descriptors, mate, program.reverse)?;
            if program.whitelist_positions
                != expected_groups
                    .iter()
                    .map(|group| group.0)
                    .collect::<Vec<_>>()
                || program.group_sizes.len() != expected_groups.len()
            {
                return Err(Error::InvalidFormat(
                    "whitelist operation group mismatch".into(),
                ));
            }
            program.descriptor_groups = expected_groups.into_iter().map(|group| group.1).collect();
            for group in &program.descriptor_groups {
                for &descriptor_index in group {
                    let descriptor = descriptors.get(descriptor_index).ok_or_else(|| {
                        Error::InvalidFormat("whitelist descriptor is out of bounds".into())
                    })?;
                    let is_raw = raw.read(1)? != 0;
                    let fragment = if is_raw {
                        WhitelistFragment::Raw(descriptor.range.len())
                    } else {
                        let mode = mismatch.read(1)?;
                        let ordinal_end = ordinal_offset.checked_add(2).ok_or_else(|| {
                            Error::InvalidFormat("ordinal column offset overflow".into())
                        })?;
                        let ordinal = u16::from_le_bytes(
                            ordinal_bytes[ordinal_offset..ordinal_end]
                                .try_into()
                                .expect("ordinal length was validated"),
                        ) as usize;
                        ordinal_offset = ordinal_end;
                        let mut sequence = global_whitelist
                            .sequence(ordinal)
                            .ok_or_else(|| {
                                Error::InvalidFormat("whitelist ordinal is out of bounds".into())
                            })?
                            .to_vec();
                        if sequence.len() != descriptor.range.len() {
                            return Err(Error::InvalidFormat(
                                "whitelist sequence length does not match descriptor".into(),
                            ));
                        }
                        if mode == 1 {
                            let detail = mismatch_detail_bytes[detail_offset];
                            detail_offset += 1;
                            let position = (detail >> 2) as usize;
                            if position >= sequence.len() {
                                return Err(Error::InvalidFormat(
                                    "whitelist mismatch position is out of bounds".into(),
                                ));
                            }
                            let replacement = bits_to_base(detail & 3);
                            if replacement == sequence[position] {
                                return Err(Error::InvalidFormat(
                                    "whitelist replacement does not differ".into(),
                                ));
                            }
                            sequence[position] = replacement;
                        } else if mode != 0 {
                            return Err(Error::InvalidFormat(
                                "invalid whitelist mismatch mode".into(),
                            ));
                        }
                        if program.reverse {
                            reverse_complement(&mut sequence);
                        }
                        WhitelistFragment::Bases(sequence)
                    };
                    program.fragments.push((descriptor_index, fragment));
                }
            }
            let sequence = decode_record_literals(program, &mut literal_cursor)?;
            all[mate].push(sequence);
        }
        if raw.remaining() != 0
            || mismatch.remaining() != 0
            || ordinal_offset != ordinal_bytes.len()
            || detail_offset != mismatch_detail_bytes.len()
        {
            return Err(Error::InvalidFormat("trailing record column data".into()));
        }
        if literal_cursor.position() != literal_cursor.get_ref().len() as u64 {
            return Err(Error::InvalidFormat("trailing record literal bytes".into()));
        }
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat("trailing record bytes".into()));
    }
    Ok(all)
}

fn decode_ids(
    input: &[u8],
    count: usize,
    mode: IdMode,
    layout: AlignmentLayout,
) -> Result<Vec<Vec<Vec<u8>>>> {
    let mut cursor = Cursor::new(input.to_vec());
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
    let second = match mode {
        IdMode::None => vec![Vec::new(); count],
        IdMode::Independent => decode_prefix_ids(&blocks[1], count)?,
        IdMode::IdenticalMates => {
            if !blocks[1].is_empty() {
                return Err(Error::InvalidFormat("unexpected mate-2 ID bytes".into()));
            }
            first.clone()
        }
        IdMode::SlashMate => {
            if !blocks[1].is_empty() {
                return Err(Error::InvalidFormat("unexpected mate-2 ID bytes".into()));
            }
            first
                .iter()
                .map(|id| {
                    let mut value = id.clone();
                    value.extend_from_slice(b"/2");
                    value
                })
                .collect()
        }
    };
    let first = if mode == IdMode::SlashMate {
        first
            .iter()
            .map(|id| {
                let mut value = id.clone();
                value.extend_from_slice(b"/1");
                value
            })
            .collect()
    } else {
        first
    };
    Ok(vec![first, second])
}

fn decode_prefix_ids(input: &[u8], count: usize) -> Result<Vec<Vec<u8>>> {
    let mut cursor = Cursor::new(input.to_vec());
    let mut previous = Vec::new();
    let mut result = Vec::with_capacity(count);
    for _ in 0..count {
        let common = usize_from(read_varint(&mut cursor)?, "ID prefix")?;
        let suffix_len = usize_from(read_varint(&mut cursor)?, "ID suffix")?;
        if common > previous.len() {
            return Err(Error::InvalidFormat("ID prefix out of bounds".into()));
        }
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
        previous = value.clone();
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
    let mut cursor = Cursor::new(input.to_vec());
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
                descriptors: Vec::new(),
            });
        }
    }
    Ok(runs)
}

fn descriptor_boundaries(
    descriptors: &[IndexEncodingDescriptor],
    mate: usize,
    reduced_len: usize,
    reverse: bool,
) -> Result<Vec<(Vec<usize>, usize)>> {
    let mut selected = descriptors
        .iter()
        .enumerate()
        .filter(|(_, descriptor)| descriptor.range.mate.stream() == mate)
        .collect::<Vec<_>>();
    selected.sort_by_key(|(_, descriptor)| descriptor.range.start);
    let mut prior_removed = 0usize;
    let mut boundaries: Vec<(Vec<usize>, usize)> = Vec::with_capacity(selected.len());
    for (descriptor_index, descriptor) in selected {
        let forward = descriptor
            .range
            .start
            .checked_sub(prior_removed)
            .ok_or_else(|| Error::InvalidInput("embedded index boundary underflow".into()))?;
        if forward > reduced_len {
            return Err(Error::InvalidInput(
                "embedded index boundary exceeds reduced read".into(),
            ));
        }
        if let Some((group, _)) = boundaries
            .last_mut()
            .filter(|(_, boundary)| *boundary == forward)
        {
            group.push(descriptor_index);
        } else {
            boundaries.push((vec![descriptor_index], forward));
        }
        prior_removed = prior_removed
            .checked_add(descriptor.range.len())
            .ok_or_else(|| Error::InvalidInput("embedded index length overflow".into()))?;
    }
    if reverse {
        boundaries.reverse();
        for (group, boundary) in &mut boundaries {
            *boundary = reduced_len - *boundary;
            group.reverse();
        }
    }
    Ok(boundaries)
}

fn insert_whitelist_runs(
    core: Vec<ArchiveRun>,
    boundaries: Vec<(Vec<usize>, usize)>,
    reduced_len: usize,
) -> Result<Vec<ArchiveRun>> {
    let mut output = Vec::with_capacity(core.len() + boundaries.len());
    let mut boundary_index = 0usize;
    let mut query_pos = 0usize;
    for run in core {
        if !run.op.consumes_query() {
            output.push(run);
            continue;
        }
        let mut remaining = run.len;
        while boundary_index < boundaries.len()
            && boundaries[boundary_index].1 <= query_pos + remaining
        {
            let (descriptors, boundary) = &boundaries[boundary_index];
            if *boundary < query_pos {
                return Err(Error::InvalidInput(
                    "embedded index boundaries are not ordered".into(),
                ));
            }
            let prefix = *boundary - query_pos;
            if prefix != 0 {
                output.push(ArchiveRun {
                    op: run.op,
                    len: prefix,
                    descriptors: Vec::new(),
                });
                remaining -= prefix;
                query_pos = *boundary;
            }
            output.push(ArchiveRun {
                op: ArchiveOp::Whitelist,
                len: 0,
                descriptors: descriptors.clone(),
            });
            boundary_index += 1;
        }
        if remaining != 0 {
            output.push(ArchiveRun {
                op: run.op,
                len: remaining,
                descriptors: Vec::new(),
            });
            query_pos += remaining;
        }
    }
    while boundary_index < boundaries.len() && boundaries[boundary_index].1 == query_pos {
        output.push(ArchiveRun {
            op: ArchiveOp::Whitelist,
            len: 0,
            descriptors: boundaries[boundary_index].0.clone(),
        });
        boundary_index += 1;
    }
    if query_pos != reduced_len || boundary_index != boundaries.len() {
        return Err(Error::InvalidInput(
            "archive operations do not cover embedded index boundaries".into(),
        ));
    }
    Ok(output)
}

#[cfg(test)]
fn encode_operation_body<R: ReferenceBases + ?Sized>(
    output: &mut Vec<u8>,
    sequence: &[u8],
    mapped: Option<(&GenomeAlignment, &R)>,
    runs: &[ArchiveRun],
    indexes: &[SpoolIndexRead],
    descriptors: &[IndexEncodingDescriptor],
    global_whitelist: &Whitelist,
) -> Result<()> {
    if indexes.len() != descriptors.len() {
        return Err(Error::InvalidInput(
            "embedded index count does not match metadata".into(),
        ));
    }
    write_varint(output, runs.len() as u64)?;
    for chunk in runs.chunks(4) {
        let mut packed = 0u8;
        for (index, run) in chunk.iter().enumerate() {
            packed |= (run.op as u8) << (6 - index * 2);
        }
        output.push(packed);
    }
    for run in runs {
        if run.op == ArchiveOp::Whitelist {
            if run.len != 0 || run.descriptors.is_empty() {
                return Err(Error::InvalidInput(
                    "whitelist operation has invalid descriptors".into(),
                ));
            }
        } else {
            if run.len == 0 || !run.descriptors.is_empty() {
                return Err(Error::InvalidInput("invalid archive operation run".into()));
            }
            write_varint(output, run.len as u64)?;
        }
    }

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
    let mut query_pos = 0usize;
    for run in runs {
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
                let mut mismatches = Vec::new();
                for offset in 0..run.len {
                    if reference.base(ref_pos + offset) != oriented[query_pos + offset] {
                        mismatches.push((offset, oriented[query_pos + offset]));
                    }
                }
                write_varint(output, mismatches.len() as u64)?;
                let mut previous = 0usize;
                for (position, base) in mismatches {
                    write_varint(output, (position - previous) as u64)?;
                    output.push(base);
                    previous = position;
                }
                ref_pos = ref_end;
                query_pos = query_end;
            }
            ArchiveOp::Literal => {
                let end = query_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if end > oriented.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds query".into(),
                    ));
                }
                query_pos = end;
            }
            ArchiveOp::RefAdvance => {
                if mapped.is_none() {
                    return Err(Error::InvalidInput(
                        "unmapped body contains a reference advance".into(),
                    ));
                }
                ref_pos = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                if ref_pos > mapped.unwrap().1.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds reference".into(),
                    ));
                }
            }
            ArchiveOp::Whitelist => {}
        }
    }
    if query_pos != oriented.len() {
        return Err(Error::InvalidFormat(
            "archive operations do not consume reduced query".into(),
        ));
    }

    let ordinal_bits = bit_width(global_whitelist.len());
    let mut classifications = Vec::new();
    for run in runs.iter().filter(|run| run.op == ArchiveOp::Whitelist) {
        for &descriptor_index in &run.descriptors {
            let descriptor = descriptors.get(descriptor_index).ok_or_else(|| {
                Error::InvalidInput("whitelist descriptor is out of bounds".into())
            })?;
            let index = indexes
                .get(descriptor_index)
                .ok_or_else(|| Error::InvalidInput("embedded index fragment is missing".into()))?;
            if index.sequence.len() != descriptor.range.len() {
                return Err(Error::InvalidInput(
                    "embedded index sequence length does not match its range".into(),
                ));
            }
            let classification = index_read::classify_with_global(
                &index.sequence,
                &descriptor.local_whitelist,
                global_whitelist,
            )
            .map_err(index_encode_error)?;
            classifications.push((descriptor_index, classification));
        }
    }

    let mut bits = BitWriter::default();
    for (_, classification) in &classifications {
        bits.write(
            usize::from(matches!(
                classification,
                IndexReadEncoding::CanonicalRaw | IndexReadEncoding::ExceptionalRaw
            )),
            1,
        )?;
    }
    for &(descriptor_index, ref classification) in &classifications {
        let descriptor = &descriptors[descriptor_index];
        match *classification {
            IndexReadEncoding::Exact { global_ordinal } => {
                if global_whitelist
                    .sequence(global_ordinal)
                    .is_none_or(|sequence| sequence.len() != descriptor.range.len())
                {
                    return Err(Error::InvalidInput(
                        "whitelist sequence length does not match descriptor".into(),
                    ));
                }
                bits.write(0, 1)?;
                bits.write(global_ordinal, ordinal_bits)?;
            }
            IndexReadEncoding::OneMismatch {
                global_ordinal,
                position,
                original,
            } => {
                let selected = global_whitelist.sequence(global_ordinal).ok_or_else(|| {
                    Error::InvalidInput("whitelist ordinal is out of bounds".into())
                })?;
                if selected.len() != descriptor.range.len() || position >= selected.len() {
                    return Err(Error::InvalidInput(
                        "whitelist mismatch does not match descriptor".into(),
                    ));
                }
                if canonical_bits(original)? == canonical_bits(selected[position])? {
                    return Err(Error::InvalidInput(
                        "whitelist replacement does not differ".into(),
                    ));
                }
                bits.write(1, 1)?;
                bits.write(global_ordinal, ordinal_bits)?;
                bits.write(position, bit_width(descriptor.range.len()))?;
                bits.write(canonical_bits(original)? as usize, 2)?;
            }
            IndexReadEncoding::CanonicalRaw | IndexReadEncoding::ExceptionalRaw => {}
        }
    }
    write_varint(output, bits.len() as u64)?;
    output.extend(bits.finish());

    let mut classification_index = 0usize;
    let mut query_pos = 0usize;
    for run in runs {
        match run.op {
            ArchiveOp::Match => {
                query_pos = query_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
            }
            ArchiveOp::Literal => {
                let end = query_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if end > oriented.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds query".into(),
                    ));
                }
                output.extend_from_slice(&oriented[query_pos..end]);
                query_pos = end;
            }
            ArchiveOp::RefAdvance => {}
            ArchiveOp::Whitelist => {
                for &descriptor_index in &run.descriptors {
                    let (classified_descriptor, classification) =
                        classifications.get(classification_index).ok_or_else(|| {
                            Error::InvalidInput("missing whitelist classification".into())
                        })?;
                    if *classified_descriptor != descriptor_index {
                        return Err(Error::InvalidInput(
                            "whitelist classification order mismatch".into(),
                        ));
                    }
                    classification_index += 1;
                    if matches!(
                        classification,
                        IndexReadEncoding::CanonicalRaw | IndexReadEncoding::ExceptionalRaw
                    ) {
                        let mut fragment = indexes[descriptor_index].sequence.clone();
                        if reverse {
                            reverse_complement(&mut fragment);
                        }
                        output.extend(fragment);
                    }
                }
            }
        }
    }
    if query_pos != oriented.len() || classification_index != classifications.len() {
        return Err(Error::InvalidInput(
            "literal source consumption mismatch".into(),
        ));
    }
    Ok(())
}

#[cfg(test)]
fn decode_operation_body(
    input: &mut Cursor<Vec<u8>>,
    mapped: Option<(&[u8], usize)>,
    reverse: bool,
    descriptors: &[IndexDescriptor],
    mate: usize,
    global_whitelist: &Whitelist,
) -> Result<Vec<u8>> {
    let count = usize_from(read_varint(input)?, "archive operation count")?;
    let packed_len = count.div_ceil(4);
    let mut packed = vec![0; packed_len];
    input.read_exact(&mut packed)?;
    if count % 4 != 0 {
        let padding_bits = (4 - count % 4) * 2;
        if packed.last().copied().unwrap_or(0) & ((1 << padding_bits) - 1) != 0 {
            return Err(Error::InvalidFormat(
                "nonzero archive operation padding bits".into(),
            ));
        }
    }
    let mut runs = Vec::with_capacity(count);
    for index in 0..count {
        let code = (packed[index / 4] >> (6 - (index % 4) * 2)) & 3;
        let op = ArchiveOp::from_code(code)?;
        let len = if op == ArchiveOp::Whitelist {
            0
        } else {
            let len = usize_from(read_varint(input)?, "archive operation length")?;
            if len == 0 {
                return Err(Error::InvalidFormat("zero-length archive operation".into()));
            }
            len
        };
        runs.push(ArchiveRun {
            op,
            len,
            descriptors: Vec::new(),
        });
    }
    if mapped.is_none()
        && runs
            .iter()
            .any(|run| matches!(run.op, ArchiveOp::Match | ArchiveOp::RefAdvance))
    {
        return Err(Error::InvalidFormat(
            "unmapped body contains mapped operations".into(),
        ));
    }

    let mut matched = (0..runs.len()).map(|_| None).collect::<Vec<_>>();
    let mut ref_pos = mapped.map_or(0, |(_, start)| start);
    let mut reduced_len = 0usize;
    let mut whitelist_positions = Vec::new();
    for (run_index, run) in runs.iter().enumerate() {
        match run.op {
            ArchiveOp::Match => {
                let reference = mapped.unwrap().0;
                let end = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if end > reference.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds reference".into(),
                    ));
                }
                let mut bases = reference[ref_pos..end].to_vec();
                let mismatch_count = usize_from(read_varint(input)?, "mismatch count")?;
                let mut previous = 0usize;
                for mismatch in 0..mismatch_count {
                    let position = previous
                        .checked_add(usize_from(read_varint(input)?, "mismatch position")?)
                        .ok_or_else(|| Error::InvalidFormat("mismatch position overflow".into()))?;
                    if position >= run.len || (mismatch != 0 && position <= previous) {
                        return Err(Error::InvalidFormat(
                            "mismatch position out of bounds".into(),
                        ));
                    }
                    bases[position] = read_u8(input)?;
                    previous = position;
                }
                ref_pos = end;
                reduced_len = reduced_len
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reduced query length overflow".into()))?;
                matched[run_index] = Some(bases);
            }
            ArchiveOp::Literal => {
                reduced_len = reduced_len
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reduced query length overflow".into()))?;
            }
            ArchiveOp::RefAdvance => {
                ref_pos = ref_pos
                    .checked_add(run.len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if ref_pos > mapped.unwrap().0.len() {
                    return Err(Error::InvalidFormat(
                        "archive operation exceeds reference".into(),
                    ));
                }
            }
            ArchiveOp::Whitelist => {
                whitelist_positions.push(reduced_len);
            }
        }
    }

    let mut selected = descriptors
        .iter()
        .enumerate()
        .filter(|(_, descriptor)| descriptor.range.mate.stream() == mate)
        .collect::<Vec<_>>();
    selected.sort_by_key(|(_, descriptor)| descriptor.range.start);
    let mut prior_removed = 0usize;
    let mut expected_groups: Vec<(usize, Vec<usize>)> = Vec::new();
    for (descriptor_index, descriptor) in selected {
        let forward = descriptor
            .range
            .start
            .checked_sub(prior_removed)
            .ok_or_else(|| Error::InvalidFormat("embedded index boundary underflow".into()))?;
        if forward > reduced_len {
            return Err(Error::InvalidFormat(
                "embedded index boundary exceeds reduced read".into(),
            ));
        }
        if let Some((_, group)) = expected_groups
            .last_mut()
            .filter(|(boundary, _)| *boundary == forward)
        {
            group.push(descriptor_index);
        } else {
            expected_groups.push((forward, vec![descriptor_index]));
        }
        prior_removed = prior_removed
            .checked_add(descriptor.range.len())
            .ok_or_else(|| Error::InvalidFormat("embedded index length overflow".into()))?;
    }
    if reverse {
        expected_groups.reverse();
        for (boundary, group) in &mut expected_groups {
            *boundary = reduced_len - *boundary;
            group.reverse();
        }
    }
    let expected_positions = expected_groups
        .iter()
        .map(|(boundary, _)| *boundary)
        .collect::<Vec<_>>();
    if whitelist_positions != expected_positions {
        return Err(Error::InvalidFormat(
            "whitelist operation count or order mismatch".into(),
        ));
    }
    for (run, (_, group)) in runs
        .iter_mut()
        .filter(|run| run.op == ArchiveOp::Whitelist)
        .zip(&expected_groups)
    {
        run.descriptors = group.clone();
    }

    let bit_count = usize_from(read_varint(input)?, "whitelist bit count")?;
    let mut bit_bytes = vec![0; bit_count.div_ceil(8)];
    input.read_exact(&mut bit_bytes)?;
    let mut bits = BitReader::new(&bit_bytes, bit_count)?;
    let ordinal_bits = bit_width(global_whitelist.len());
    let descriptor_order = expected_groups
        .iter()
        .flat_map(|(_, group)| group.iter().copied())
        .collect::<Vec<_>>();
    let raw_bitmap = (0..descriptor_order.len())
        .map(|_| bits.read(1).map(|value| value != 0))
        .collect::<Result<Vec<_>>>()?;

    enum WhitelistFragment {
        Raw(usize),
        Bases(Vec<u8>),
    }
    let mut fragments = (0..descriptors.len()).map(|_| None).collect::<Vec<_>>();
    for (&descriptor_index, raw) in descriptor_order.iter().zip(raw_bitmap) {
        let descriptor = descriptors
            .get(descriptor_index)
            .ok_or_else(|| Error::InvalidFormat("whitelist descriptor is out of bounds".into()))?;
        let length = descriptor.range.len();
        let fragment = if raw {
            WhitelistFragment::Raw(length)
        } else {
            let mode = bits.read(1)?;
            let ordinal = bits.read(ordinal_bits)?;
            let mut sequence = global_whitelist
                .sequence(ordinal)
                .ok_or_else(|| Error::InvalidFormat("whitelist ordinal is out of bounds".into()))?
                .to_vec();
            if sequence.len() != length {
                return Err(Error::InvalidFormat(
                    "whitelist sequence length does not match descriptor".into(),
                ));
            }
            if mode == 1 {
                let position = bits.read(bit_width(length))?;
                if position >= length {
                    return Err(Error::InvalidFormat(
                        "whitelist mismatch position is out of bounds".into(),
                    ));
                }
                let original = bits_to_base(bits.read(2)? as u8);
                if original == sequence[position] {
                    return Err(Error::InvalidFormat(
                        "whitelist replacement does not differ".into(),
                    ));
                }
                sequence[position] = original;
            }
            if reverse {
                reverse_complement(&mut sequence);
            }
            WhitelistFragment::Bases(sequence)
        };
        if fragments[descriptor_index].replace(fragment).is_some() {
            return Err(Error::InvalidFormat(
                "duplicate whitelist descriptor".into(),
            ));
        }
    }
    if bits.remaining() != 0 {
        return Err(Error::InvalidFormat("trailing whitelist bits".into()));
    }

    enum Piece {
        Bases(Vec<u8>),
        Literal(usize),
    }
    let mut pieces = Vec::new();
    for (run_index, run) in runs.iter().enumerate() {
        match run.op {
            ArchiveOp::Match => {
                pieces.push(Piece::Bases(matched[run_index].take().ok_or_else(
                    || Error::InvalidFormat("missing matched sequence".into()),
                )?))
            }
            ArchiveOp::Literal => pieces.push(Piece::Literal(run.len)),
            ArchiveOp::RefAdvance => {}
            ArchiveOp::Whitelist => {
                for &descriptor_index in &run.descriptors {
                    match fragments[descriptor_index]
                        .take()
                        .ok_or_else(|| Error::InvalidFormat("missing whitelist fragment".into()))?
                    {
                        WhitelistFragment::Raw(length) => pieces.push(Piece::Literal(length)),
                        WhitelistFragment::Bases(sequence) => pieces.push(Piece::Bases(sequence)),
                    }
                }
            }
        }
    }
    if fragments.into_iter().any(|fragment| fragment.is_some()) {
        return Err(Error::InvalidFormat("unused whitelist fragment".into()));
    }

    let mut sequence = Vec::new();
    let mut piece_index = 0usize;
    while piece_index < pieces.len() {
        match &pieces[piece_index] {
            Piece::Bases(bases) => {
                sequence.extend_from_slice(bases);
                piece_index += 1;
            }
            Piece::Literal(_) => {
                let start = piece_index;
                let mut length = 0usize;
                while piece_index < pieces.len() {
                    let Piece::Literal(source_length) = pieces[piece_index] else {
                        break;
                    };
                    length = length.checked_add(source_length).ok_or_else(|| {
                        Error::InvalidFormat("literal run length overflow".into())
                    })?;
                    piece_index += 1;
                }
                let literal = take_literal_bytes(input, length)?;
                let mut offset = 0usize;
                for piece in &pieces[start..piece_index] {
                    let Piece::Literal(source_length) = *piece else {
                        unreachable!();
                    };
                    let end = offset.checked_add(source_length).ok_or_else(|| {
                        Error::InvalidFormat("literal source length overflow".into())
                    })?;
                    if end > literal.len() {
                        return Err(Error::InvalidFormat("truncated literal run".into()));
                    }
                    sequence.extend_from_slice(&literal[offset..end]);
                    offset = end;
                }
                if offset != literal.len() {
                    return Err(Error::InvalidFormat(
                        "literal source splitting mismatch".into(),
                    ));
                }
            }
        }
    }
    if reverse {
        reverse_complement(&mut sequence);
    }
    Ok(sequence)
}

#[cfg(test)]
fn bit_width(values: usize) -> usize {
    if values <= 1 {
        0
    } else {
        (usize::BITS - (values - 1).leading_zeros()) as usize
    }
}

fn canonical_bits(base: u8) -> Result<u8> {
    match base {
        b'A' => Ok(0),
        b'C' => Ok(1),
        b'G' => Ok(2),
        b'T' => Ok(3),
        _ => Err(Error::InvalidInput("noncanonical whitelist base".into())),
    }
}

fn bits_to_base(bits: u8) -> u8 {
    b"ACGT"[bits as usize]
}

#[derive(Default)]
struct BitWriter {
    bytes: Vec<u8>,
    bits: usize,
}

impl BitWriter {
    fn write(&mut self, value: usize, width: usize) -> Result<()> {
        if width < usize::BITS as usize && value >= (1usize << width) {
            return Err(Error::InvalidInput(
                "bit field value is out of range".into(),
            ));
        }
        for shift in (0..width).rev() {
            if self.bits % 8 == 0 {
                self.bytes.push(0);
            }
            let byte = self.bytes.len() - 1;
            self.bytes[byte] |= (((value >> shift) & 1) as u8) << (7 - self.bits % 8);
            self.bits += 1;
        }
        Ok(())
    }

    fn len(&self) -> usize {
        self.bits
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }
}

struct BitReader<'a> {
    bytes: &'a [u8],
    bits: usize,
    position: usize,
}

impl<'a> BitReader<'a> {
    fn new(bytes: &'a [u8], bits: usize) -> Result<Self> {
        if bytes.len() != bits.div_ceil(8) {
            return Err(Error::InvalidFormat(
                "whitelist bitstream length mismatch".into(),
            ));
        }
        if bits % 8 != 0 && bytes.last().copied().unwrap_or(0) & ((1 << (8 - bits % 8)) - 1) != 0 {
            return Err(Error::InvalidFormat(
                "nonzero whitelist bitstream padding".into(),
            ));
        }
        Ok(Self {
            bytes,
            bits,
            position: 0,
        })
    }

    fn read(&mut self, width: usize) -> Result<usize> {
        if width > self.remaining() {
            return Err(Error::InvalidFormat("truncated whitelist bitstream".into()));
        }
        let mut value = 0usize;
        for _ in 0..width {
            value = (value << 1)
                | ((self.bytes[self.position / 8] >> (7 - self.position % 8)) & 1) as usize;
            self.position += 1;
        }
        Ok(value)
    }

    fn remaining(&self) -> usize {
        self.bits - self.position
    }
}

#[cfg(test)]
fn decode_mapped_body(
    input: &mut Cursor<Vec<u8>>,
    reference: &[u8],
    start: i64,
    is_reverse: bool,
) -> Result<Vec<u8>> {
    let whitelist = index_read::merge_whitelists(&[]).map_err(index_decode_error)?;
    decode_operation_body(
        input,
        Some((reference, nonnegative_usize(start, "alignment start")?)),
        is_reverse,
        &[],
        0,
        &whitelist,
    )
}

fn effective_cigar(hit: &GenomeAlignment, read_len: usize) -> Result<Vec<CigarOp>> {
    let query_start = usize::try_from(hit.query_start)
        .map_err(|_| Error::InvalidFormat("negative query start".into()))?;
    let query_end = usize::try_from(hit.query_end)
        .map_err(|_| Error::InvalidFormat("negative query end".into()))?;
    if query_start > query_end || query_end > read_len {
        return Err(Error::InvalidFormat(
            "query alignment range out of bounds".into(),
        ));
    }
    let (prefix, suffix) = if hit.is_reverse {
        (read_len - query_end, query_start)
    } else {
        (query_start, read_len - query_end)
    };
    let mut cigar = Vec::with_capacity(hit.cigar.len() + 2);
    if prefix != 0 {
        cigar.push(CigarOp {
            kind: CigarKind::SoftClip,
            len: prefix as u32,
        });
    }
    cigar.extend_from_slice(&hit.cigar);
    if suffix != 0 {
        cigar.push(CigarOp {
            kind: CigarKind::SoftClip,
            len: suffix as u32,
        });
    }
    normalize_cigar(cigar)
}

fn normalize_cigar(cigar: impl IntoIterator<Item = CigarOp>) -> Result<Vec<CigarOp>> {
    let mut normalized: Vec<CigarOp> = Vec::new();
    for op in cigar {
        let kind = match op.kind {
            CigarKind::Eq | CigarKind::Diff => CigarKind::Match,
            CigarKind::HardClip | CigarKind::Pad => continue,
            kind => kind,
        };
        if op.len == 0 {
            continue;
        }
        if let Some(previous) = normalized
            .last_mut()
            .filter(|previous| previous.kind == kind)
        {
            previous.len = previous
                .len
                .checked_add(op.len)
                .ok_or_else(|| Error::InvalidInput("CIGAR length overflow".into()))?;
        } else {
            normalized.push(CigarOp { kind, len: op.len });
        }
    }
    Ok(normalized)
}

fn encode_index_metadata(
    descriptors: &[IndexEncodingDescriptor],
    global_whitelist: &Whitelist,
) -> Result<Vec<u8>> {
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
    let whitelist = index_read::encode_whitelist(global_whitelist).map_err(index_encode_error)?;
    write_u64(&mut output, whitelist.len() as u64)?;
    output.extend(whitelist);
    Ok(output)
}

fn decode_index_metadata(
    input: &[u8],
    layout: AlignmentLayout,
) -> Result<(Vec<IndexDescriptor>, Whitelist)> {
    if input.len() < 16 {
        return Err(Error::InvalidFormat(
            "embedded index metadata is too short".into(),
        ));
    }
    let mut input = Cursor::new(input.to_vec());
    let count = usize_from(read_u64(&mut input)?, "embedded index range count")?;
    let minimum_descriptors = count
        .checked_mul(17)
        .and_then(|length| length.checked_add(8))
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
    let whitelist_len = usize_from(read_u64(&mut input)?, "embedded index whitelist")?;
    let whitelist_start = input.position() as usize;
    let whitelist_end = whitelist_start
        .checked_add(whitelist_len)
        .ok_or_else(|| Error::InvalidFormat("embedded index whitelist overflow".into()))?;
    if whitelist_end > input.get_ref().len() {
        return Err(Error::InvalidFormat(
            "embedded index whitelist exceeds metadata".into(),
        ));
    }
    let global_whitelist =
        index_read::decode_whitelist(&input.get_ref()[whitelist_start..whitelist_end])
            .map_err(index_decode_error)?;
    input.set_position(whitelist_end as u64);
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
    Ok((descriptors, global_whitelist))
}

fn read_header<R: Read>(input: &mut R) -> Result<ArchiveHeader> {
    let mut magic = [0; 9];
    input.read_exact(&mut magic)?;
    if magic != *MAGIC {
        return Err(Error::InvalidFormat("bad alignment archive magic".into()));
    }
    let version = read_u32(input)?;
    if version != VERSION {
        return Err(Error::InvalidFormat(
            "unsupported alignment archive version".into(),
        ));
    }
    let flags = read_u16(input)?;
    if flags & !0x1f != 0 {
        return Err(Error::InvalidFormat(
            "unknown alignment archive flags".into(),
        ));
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
    let reference_len = usize_from(read_u64(input)?, "reference block")?;
    let index_metadata_len = usize_from(read_u64(input)?, "embedded index metadata block")?;
    let chunks = read_u64(input)?;
    let mut reference_block = vec![0; reference_len];
    input.read_exact(&mut reference_block)?;
    let mut index_metadata = vec![0; index_metadata_len];
    input.read_exact(&mut index_metadata)?;
    let (indexes, global_whitelist) = decode_index_metadata(&index_metadata, layout)?;
    Ok(ArchiveHeader {
        flags,
        id_mode,
        layout,
        chunk_count: chunks,
        reference: decode_reference(&reference_block)?,
        indexes,
        global_whitelist,
    })
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
        block: zstd::stream::encode_all(Cursor::new(&raw), options.zstd_level)?,
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
    let mut decoder = zstd::stream::read::Decoder::new(Cursor::new(block))?.single_frame();
    let mut count_bytes = [0; 8];
    decoder
        .read_exact(&mut count_bytes)
        .map_err(|_| invalid_payload())?;
    let base_count = usize_from(u64::from_le_bytes(count_bytes), "reference base count")?;
    let packed_len = base_count.div_ceil(4);
    let mut decoded_packed_len = 0usize;
    let mut last_packed = 0u8;
    let mut buffer = [0; 8192];
    loop {
        let read = decoder.read(&mut buffer).map_err(|_| invalid_payload())?;
        if read == 0 {
            break;
        }
        decoded_packed_len = decoded_packed_len
            .checked_add(read)
            .ok_or_else(|| Error::InvalidFormat("reference payload length overflow".into()))?;
        last_packed = buffer[read - 1];
    }
    if decoded_packed_len != packed_len {
        return Err(Error::InvalidFormat(
            "reference packed payload length mismatch".into(),
        ));
    }
    let remainder = base_count % 4;
    if remainder != 0 && last_packed & ((1 << ((4 - remainder) * 2)) - 1) != 0 {
        return Err(Error::InvalidFormat(
            "nonzero reference padding bits".into(),
        ));
    }

    let mut reference = Vec::new();
    reference
        .try_reserve_exact(base_count)
        .map_err(|_| Error::InvalidFormat("reference base count is too large".into()))?;
    let mut decoder = zstd::stream::read::Decoder::new(Cursor::new(block))?.single_frame();
    decoder
        .read_exact(&mut count_bytes)
        .map_err(|_| invalid_payload())?;
    for _ in 0..packed_len {
        let mut packed = [0];
        decoder
            .read_exact(&mut packed)
            .map_err(|_| invalid_payload())?;
        for shift in [6, 4, 2, 0] {
            if reference.len() == base_count {
                break;
            }
            reference.push(b"ACGT"[((packed[0] >> shift) & 3) as usize]);
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

fn take_literal_bytes(input: &mut Cursor<Vec<u8>>, length: usize) -> Result<Vec<u8>> {
    let start = usize::try_from(input.position())
        .map_err(|_| Error::InvalidFormat("record literal offset is too large".into()))?;
    let end = start
        .checked_add(length)
        .ok_or_else(|| Error::InvalidFormat("record literal length overflow".into()))?;
    if end > input.get_ref().len() {
        return Err(Error::InvalidFormat(
            "truncated record literal bytes".into(),
        ));
    }
    let literal = input.get_ref()[start..end].to_vec();
    input.set_position(end as u64);
    Ok(literal)
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
fn reverse_complement(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => *base,
        };
    }
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

fn read_varint(input: &mut Cursor<Vec<u8>>) -> Result<u64> {
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

    struct RecordDirectory {
        start: usize,
        values: [usize; 8],
        payload_start: usize,
    }

    fn record_directory(encoded: &[u8], count: usize) -> RecordDirectory {
        let mut cursor = Cursor::new(encoded.to_vec());
        cursor.set_position((count.div_ceil(8) * 2) as u64);
        let positions_len = read_u64(&mut cursor).unwrap() as usize;
        cursor.set_position(cursor.position() + positions_len as u64);
        let start = cursor.position() as usize;
        let mut values = [0; 8];
        for value in &mut values {
            *value = read_u64(&mut cursor).unwrap() as usize;
        }
        RecordDirectory {
            start,
            values,
            payload_start: cursor.position() as usize,
        }
    }

    fn record_literals(encoded: &[u8], count: usize) -> &[u8] {
        let directory = record_directory(encoded, count);
        let literal_start = directory.payload_start
            + directory.values[0]
            + directory.values[2]
            + directory.values[4]
            + directory.values[5]
            + directory.values[6];
        &encoded[literal_start..literal_start + directory.values[7]]
    }

    fn single_index_record_stream(
        index_sequence: &[u8],
        whitelist_sequences: &[Vec<u8>],
    ) -> (Vec<u8>, Vec<IndexDescriptor>, Whitelist) {
        let range = EmbeddedIndexRange {
            mate: EmbeddedIndexMate::R1,
            start: 0,
            end: index_sequence.len(),
        };
        let mut builder = WhitelistBuilder::new();
        for sequence in whitelist_sequences {
            builder.add(sequence).unwrap();
        }
        let local = builder.finish().unwrap();
        let global = index_read::merge_whitelists(&[local.clone()]).unwrap();
        let encoding = vec![IndexEncodingDescriptor {
            range: range.clone(),
            local_whitelist: local,
        }];
        let descriptors = vec![IndexDescriptor { range }];
        let pair = AlignedPair {
            first: ReadRecord::new("read", b"G".to_vec(), None),
            second: None,
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: vec![SpoolIndexRead {
                sequence: index_sequence.to_vec(),
                quality: None,
            }],
        };
        let encoded = write_record_stream(
            b"".as_slice(),
            &[pair],
            AlignmentLayout::SingleEnd,
            &encoding,
            &global,
        )
        .unwrap();
        (encoded, descriptors, global)
    }

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

    fn encoding_descriptors(ranges: &[EmbeddedIndexRange]) -> Vec<IndexEncodingDescriptor> {
        ranges
            .iter()
            .cloned()
            .map(|range| {
                let mut builder = WhitelistBuilder::new();
                builder.add(&vec![b'A'; range.len()]).unwrap();
                IndexEncodingDescriptor {
                    range,
                    local_whitelist: builder.finish().unwrap(),
                }
            })
            .collect()
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
            .map(|range| {
                let mut builder = WhitelistBuilder::new();
                builder.add(b"AC").unwrap();
                IndexEncodingDescriptor {
                    range,
                    local_whitelist: builder.finish().unwrap(),
                }
            })
            .collect::<Vec<_>>();
        let pair = AlignedPair {
            first,
            second: Some(second),
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: extracted,
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
    fn v20_chunk_round_trips_embedded_indexes() {
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
            let local_whitelists = ranges
                .iter()
                .enumerate()
                .map(|(index, _)| {
                    let mut builder = WhitelistBuilder::new();
                    builder.add(&indexes[index].sequence).unwrap();
                    builder.finish().unwrap()
                })
                .collect::<Vec<_>>();
            let global_whitelist = index_read::merge_whitelists(&local_whitelists).unwrap();
            let encoding_descriptors = ranges
                .iter()
                .cloned()
                .zip(local_whitelists)
                .map(|(range, local_whitelist)| IndexEncodingDescriptor {
                    range,
                    local_whitelist,
                })
                .collect::<Vec<_>>();
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
            };
            let options = crate::AlignmentCompressionOptions {
                preserve_quality,
                zstd_level: 1,
                ..crate::AlignmentCompressionOptions::default()
            };
            let encoded = encode_chunk(
                &[],
                &[pair],
                &options,
                IdMode::Independent,
                AlignmentLayout::SingleEnd,
                &encoding_descriptors,
                &global_whitelist,
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
                &global_whitelist,
            )
            .unwrap();
            let (_, group) = decoded
                .next_group(AlignmentLayout::SingleEnd, None)
                .unwrap();
            let AlignedRecordGroup::Single(record) = group else {
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
    fn reads_v20_header_layout() {
        let options = crate::AlignmentCompressionOptions::default();
        let reference = encode_reference(b"ACGT", &options).unwrap().block;
        let range = EmbeddedIndexRange {
            mate: EmbeddedIndexMate::R1,
            start: 1,
            end: 3,
        };
        let mut builder = WhitelistBuilder::new();
        builder.add(b"CG").unwrap();
        let local_whitelist = builder.finish().unwrap();
        let global_whitelist = index_read::merge_whitelists(&[local_whitelist.clone()]).unwrap();
        let descriptors = vec![IndexEncodingDescriptor {
            range: range.clone(),
            local_whitelist,
        }];
        let metadata = encode_index_metadata(&descriptors, &global_whitelist).unwrap();
        let mut header = Vec::new();
        header.extend(MAGIC);
        write_u32(&mut header, VERSION).unwrap();
        write_u16(&mut header, 4).unwrap();
        write_u64(&mut header, reference.len() as u64).unwrap();
        write_u64(&mut header, metadata.len() as u64).unwrap();
        write_u64(&mut header, 9).unwrap();
        header.extend(reference);
        header.extend(metadata);
        let decoded = read_header(&mut Cursor::new(header)).unwrap();
        assert_eq!(decoded.chunk_count, 9);
        assert_eq!(decoded.indexes.len(), 1);
        assert_eq!(decoded.indexes[0].range, range);
        assert_eq!(
            index_read::encode_whitelist(&decoded.global_whitelist).unwrap(),
            index_read::encode_whitelist(&global_whitelist).unwrap()
        );
    }

    #[test]
    fn v20_metadata_serializes_shared_whitelist_entries_once() {
        let ranges = [
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 0,
                end: 2,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 2,
                end: 4,
            },
        ];
        let local_whitelists = ranges
            .iter()
            .map(|_| {
                let mut builder = WhitelistBuilder::new();
                builder.add(b"CG").unwrap();
                builder.finish().unwrap()
            })
            .collect::<Vec<_>>();
        let global_whitelist = index_read::merge_whitelists(&local_whitelists).unwrap();
        let descriptors = ranges
            .into_iter()
            .zip(local_whitelists)
            .map(|(range, local_whitelist)| IndexEncodingDescriptor {
                range,
                local_whitelist,
            })
            .collect::<Vec<_>>();

        let metadata = encode_index_metadata(&descriptors, &global_whitelist).unwrap();
        let whitelist_payload = index_read::encode_whitelist(&global_whitelist).unwrap();
        assert_eq!(metadata.len(), 8 + 2 * 17 + 8 + whitelist_payload.len());
        assert_eq!(
            metadata.windows(2).filter(|bytes| *bytes == b"CG").count(),
            1
        );

        let (decoded, decoded_whitelist) =
            decode_index_metadata(&metadata, AlignmentLayout::SingleEnd).unwrap();
        assert_eq!(decoded.len(), 2);
        assert_eq!(
            index_read::encode_whitelist(&decoded_whitelist).unwrap(),
            whitelist_payload
        );
    }

    #[test]
    fn v20_empty_index_metadata_includes_an_empty_global_whitelist() {
        let global_whitelist = index_read::merge_whitelists(&[]).unwrap();
        let metadata = encode_index_metadata(&[], &global_whitelist).unwrap();
        assert_eq!(metadata.len(), 17);

        let (descriptors, decoded_whitelist) =
            decode_index_metadata(&metadata, AlignmentLayout::SingleEnd).unwrap();
        assert!(descriptors.is_empty());
        assert_eq!(
            index_read::encode_whitelist(&decoded_whitelist).unwrap(),
            index_read::encode_whitelist(&global_whitelist).unwrap()
        );
    }

    #[test]
    fn local_and_global_whitelist_ordinal_differences_round_trip() {
        let mut first_builder = WhitelistBuilder::new();
        for sequence in [b"TTTT".as_slice(), b"TTTT", b"AAAA"] {
            first_builder.add(sequence).unwrap();
        }
        let first_local = first_builder.finish().unwrap();
        let mut second_builder = WhitelistBuilder::new();
        second_builder.add(b"AAAA").unwrap();
        let second_local = second_builder.finish().unwrap();
        let global_whitelist =
            index_read::merge_whitelists(&[first_local.clone(), second_local.clone()]).unwrap();
        let encoding_descriptors = vec![
            IndexEncodingDescriptor {
                range: EmbeddedIndexRange {
                    mate: EmbeddedIndexMate::R1,
                    start: 0,
                    end: 4,
                },
                local_whitelist: first_local,
            },
            IndexEncodingDescriptor {
                range: EmbeddedIndexRange {
                    mate: EmbeddedIndexMate::R1,
                    start: 4,
                    end: 8,
                },
                local_whitelist: second_local,
            },
        ];
        let descriptors = encoding_descriptors
            .iter()
            .map(|descriptor| IndexDescriptor {
                range: descriptor.range.clone(),
            })
            .collect::<Vec<_>>();
        let pair = AlignedPair {
            first: ReadRecord::new("read", b"G".to_vec(), None),
            second: None,
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: vec![
                SpoolIndexRead {
                    sequence: b"TTTT".to_vec(),
                    quality: None,
                },
                SpoolIndexRead {
                    sequence: b"AAAA".to_vec(),
                    quality: None,
                },
            ],
        };

        let encoded = write_record_stream(
            b"".as_slice(),
            &[pair],
            AlignmentLayout::SingleEnd,
            &encoding_descriptors,
            &global_whitelist,
        )
        .unwrap();
        let decoded = decode_records(
            &encoded,
            1,
            &[],
            AlignmentLayout::SingleEnd,
            &descriptors,
            &global_whitelist,
        )
        .unwrap();
        assert_eq!(decoded[0][0], b"TTTTAAAAG");
    }

    #[test]
    fn v20_record_directory_and_u16_ordinals_are_byte_aligned() {
        let sequences = (0..257)
            .map(|mut value| {
                let mut sequence = vec![b'A'; 5];
                for position in (0..5).rev() {
                    sequence[position] = b"ACGT"[value % 4];
                    value /= 4;
                }
                sequence
            })
            .collect::<Vec<_>>();
        let (encoded, descriptors, global) =
            single_index_record_stream(&sequences[256], &sequences);
        let directory = record_directory(&encoded, 1);
        assert_eq!(directory.payload_start - directory.start, 64);
        assert_eq!(directory.values[1..7], [1, 1, 1, 1, 2, 0]);

        let ordinal_start = directory.payload_start
            + directory.values[0]
            + directory.values[2]
            + directory.values[4];
        assert_eq!(&encoded[ordinal_start..ordinal_start + 2], &[0x00, 0x01]);
        assert_eq!(
            decode_records(
                &encoded,
                1,
                &[],
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
            .unwrap()[0][0],
            [sequences[256].as_slice(), b"G"].concat()
        );
    }

    #[test]
    fn v20_decoder_rejects_odd_truncated_and_extra_ordinal_columns() {
        let sequence = b"AAAA".to_vec();
        let (encoded, descriptors, global) =
            single_index_record_stream(&sequence, std::slice::from_ref(&sequence));
        let directory = record_directory(&encoded, 1);
        let ordinal_len_offset = directory.start + 5 * 8;
        let ordinal_start = directory.payload_start
            + directory.values[0]
            + directory.values[2]
            + directory.values[4];
        let decode = |bytes: Vec<u8>| {
            decode_records(
                &bytes,
                1,
                &[],
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
        };

        let mut odd = encoded.clone();
        odd[ordinal_len_offset..ordinal_len_offset + 8].copy_from_slice(&1u64.to_le_bytes());
        assert!(decode(odd).is_err());

        let mut truncated = encoded.clone();
        truncated.remove(ordinal_start + 1);
        assert!(decode(truncated).is_err());

        let mut extra = encoded;
        extra[ordinal_len_offset..ordinal_len_offset + 8].copy_from_slice(&4u64.to_le_bytes());
        assert!(decode(extra).is_err());
    }

    #[test]
    fn v20_mismatch_details_use_one_validated_byte() {
        let whitelist_sequence = b"AAAAA".to_vec();
        let (encoded, descriptors, global) =
            single_index_record_stream(b"AATAA", std::slice::from_ref(&whitelist_sequence));
        let directory = record_directory(&encoded, 1);
        assert_eq!(directory.values[3..7], [1, 1, 2, 1]);
        let detail_start = directory.payload_start
            + directory.values[0]
            + directory.values[2]
            + directory.values[4]
            + directory.values[5];
        assert_eq!(encoded[detail_start], (2 << 2) | 3);
        assert_eq!(
            decode_records(
                &encoded,
                1,
                &[],
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
            .unwrap()[0][0],
            b"AATAAG"
        );

        let decode = |bytes: Vec<u8>| {
            decode_records(
                &bytes,
                1,
                &[],
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
        };
        let mut invalid_position = encoded.clone();
        invalid_position[detail_start] = (63 << 2) | 3;
        assert!(decode(invalid_position).is_err());
        let mut unchanged_base = encoded;
        unchanged_base[detail_start] = 2 << 2;
        assert!(decode(unchanged_base).is_err());
    }

    #[test]
    fn v20_long_one_mismatch_falls_back_to_raw_before_columns() {
        let whitelist_sequence = vec![b'A'; 65];
        let mut index_sequence = whitelist_sequence.clone();
        index_sequence[64] = b'C';
        let (encoded, descriptors, global) =
            single_index_record_stream(&index_sequence, std::slice::from_ref(&whitelist_sequence));
        let directory = record_directory(&encoded, 1);
        assert_eq!(directory.values[1..7], [1, 1, 0, 0, 0, 0]);
        let raw_start = directory.payload_start + directory.values[0];
        assert_eq!(encoded[raw_start], 0x80);
        assert_eq!(
            decode_records(
                &encoded,
                1,
                &[],
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
            .unwrap()[0][0],
            [index_sequence.as_slice(), b"G"].concat()
        );
    }

    #[test]
    fn v20_record_bitmap_padding_must_be_zero() {
        let sequence = b"AAAA".to_vec();
        let (encoded, descriptors, global) =
            single_index_record_stream(&sequence, std::slice::from_ref(&sequence));
        let directory = record_directory(&encoded, 1);
        let raw_start = directory.payload_start + directory.values[0];
        let mismatch_start = raw_start + directory.values[2];
        let decode = |bytes: Vec<u8>| {
            decode_records(
                &bytes,
                1,
                &[],
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
        };

        let mut raw_padding = encoded.clone();
        raw_padding[raw_start] |= 1;
        assert!(decode(raw_padding).is_err());
        let mut mismatch_padding = encoded;
        mismatch_padding[mismatch_start] |= 1;
        assert!(decode(mismatch_padding).is_err());
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
                    descriptors: Vec::new(),
                },
                ArchiveRun {
                    op: ArchiveOp::Literal,
                    len: 3,
                    descriptors: Vec::new(),
                },
                ArchiveRun {
                    op: ArchiveOp::RefAdvance,
                    len: 9,
                    descriptors: Vec::new(),
                },
            ]
        );
        assert_eq!(
            [
                ArchiveOp::Match,
                ArchiveOp::Literal,
                ArchiveOp::RefAdvance,
                ArchiveOp::Whitelist,
            ]
            .map(|op| op as u8),
            [0, 1, 2, 3]
        );
        for code in 0..=3 {
            assert_eq!(ArchiveOp::from_code(code).unwrap() as u8, code);
        }
    }

    #[test]
    fn adjacent_whitelist_ranges_share_one_operation() {
        let ranges = (0..5)
            .map(|start| EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start,
                end: start + 1,
            })
            .collect::<Vec<_>>();
        let descriptors = encoding_descriptors(&ranges);
        let runs = insert_whitelist_runs(
            vec![ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            }],
            descriptor_boundaries(&descriptors, 0, 1, false).unwrap(),
            1,
        )
        .unwrap();

        assert_eq!(
            runs,
            vec![
                ArchiveRun {
                    op: ArchiveOp::Whitelist,
                    len: 0,
                    descriptors: vec![0, 1, 2, 3, 4],
                },
                ArchiveRun {
                    op: ArchiveOp::Literal,
                    len: 1,
                    descriptors: Vec::new(),
                },
            ]
        );
    }

    #[test]
    fn nonadjacent_whitelist_ranges_remain_separate_operations() {
        let ranges = vec![
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 0,
                end: 1,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 2,
                end: 3,
            },
        ];
        let descriptors = encoding_descriptors(&ranges);
        let runs = insert_whitelist_runs(
            vec![ArchiveRun {
                op: ArchiveOp::Literal,
                len: 2,
                descriptors: Vec::new(),
            }],
            descriptor_boundaries(&descriptors, 0, 2, false).unwrap(),
            2,
        )
        .unwrap();

        assert_eq!(
            runs.iter()
                .filter(|run| run.op == ArchiveOp::Whitelist)
                .map(|run| run.descriptors.as_slice())
                .collect::<Vec<_>>(),
            vec![&[0][..], &[1][..]]
        );
    }

    #[test]
    fn reverse_adjacent_whitelist_group_reverses_descriptor_order() {
        let ranges = (0..5)
            .map(|start| EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start,
                end: start + 1,
            })
            .collect::<Vec<_>>();
        let descriptors = encoding_descriptors(&ranges);

        assert_eq!(
            descriptor_boundaries(&descriptors, 0, 1, true).unwrap(),
            vec![(vec![4, 3, 2, 1, 0], 1)]
        );
    }

    #[test]
    fn archive_opcode_padding_must_be_zero() {
        let mut body = Vec::new();
        write_varint(&mut body, 1).unwrap();
        body.push(1);
        write_varint(&mut body, 1).unwrap();
        write_varint(&mut body, 0).unwrap();
        write_varint(&mut body, 0).unwrap();
        write_varint(&mut body, 0).unwrap();
        let whitelist = index_read::merge_whitelists(&[]).unwrap();
        assert!(matches!(
            decode_operation_body(
                &mut Cursor::new(body),
                Some((b"A", 0)),
                false,
                &[],
                0,
                &whitelist,
            ),
            Err(Error::InvalidFormat(message))
                if message == "nonzero archive operation padding bits"
        ));
    }

    #[test]
    fn whitelist_modes_and_adjacent_unmapped_ranges_round_trip() {
        let ranges = [
            (0, 4, b"ACGT".as_slice(), b"ACGT".as_slice()),
            (4, 8, b"ACGA".as_slice(), b"ACGT".as_slice()),
            (8, 12, b"TTAA".as_slice(), b"CCCC".as_slice()),
            (12, 16, b"ANNT".as_slice(), b"GGGG".as_slice()),
            (16, 20, b"AGGA".as_slice(), b"ACGT".as_slice()),
        ];
        let mut encoding = Vec::new();
        let mut indexes = Vec::new();
        let mut locals = Vec::new();
        let mut original = Vec::new();
        for &(start, end, sequence, whitelist_sequence) in &ranges {
            let mut builder = WhitelistBuilder::new();
            builder.add(whitelist_sequence).unwrap();
            let local = builder.finish().unwrap();
            locals.push(local.clone());
            encoding.push(IndexEncodingDescriptor {
                range: EmbeddedIndexRange {
                    mate: EmbeddedIndexMate::R1,
                    start,
                    end,
                },
                local_whitelist: local,
            });
            indexes.push(SpoolIndexRead {
                sequence: sequence.to_vec(),
                quality: None,
            });
            original.extend_from_slice(sequence);
        }
        original.push(b'G');
        let global = index_read::merge_whitelists(&locals).unwrap();
        let descriptors = encoding
            .iter()
            .map(|descriptor| IndexDescriptor {
                range: descriptor.range.clone(),
            })
            .collect::<Vec<_>>();
        let pair = AlignedPair {
            first: ReadRecord::new("read", b"G".to_vec(), None),
            second: None,
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes,
        };
        let encoded = write_record_stream(
            b"".as_slice(),
            &[pair],
            AlignmentLayout::SingleEnd,
            &encoding,
            &global,
        )
        .unwrap();
        assert_eq!(record_literals(&encoded, 1), b"TTAAANNTAGGAG");
        let decoded = decode_records(
            &encoded,
            1,
            &[],
            AlignmentLayout::SingleEnd,
            &descriptors,
            &global,
        )
        .unwrap();
        assert_eq!(decoded[0][0], original);
    }

    #[test]
    fn whitelist_placement_round_trips_forward_and_reverse_mapped_reads() {
        let ranges = vec![
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 1,
                end: 3,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 5,
                end: 7,
            },
        ];
        let original = b"ACGTACTGAA";
        for reverse in [false, true] {
            let mut record = ReadRecord::new("read", original.to_vec(), Some(vec![b'I'; 10]));
            let indexes = split_pair_indexes(&mut record, None, &ranges, false).unwrap();
            let mut locals = Vec::new();
            for index in &indexes {
                let mut builder = WhitelistBuilder::new();
                builder.add(&index.sequence).unwrap();
                locals.push(builder.finish().unwrap());
            }
            let global = index_read::merge_whitelists(&locals).unwrap();
            let encoding = ranges
                .iter()
                .cloned()
                .zip(locals)
                .map(|(range, local_whitelist)| IndexEncodingDescriptor {
                    range,
                    local_whitelist,
                })
                .collect::<Vec<_>>();
            let descriptors = ranges
                .iter()
                .cloned()
                .map(|range| IndexDescriptor { range })
                .collect::<Vec<_>>();
            let mut reference = record.sequence.clone();
            if reverse {
                reverse_complement(&mut reference);
            }
            let hit = GenomeAlignment {
                start: 0,
                query_start: 0,
                query_end: record.sequence.len() as i32,
                is_reverse: reverse,
                cigar: vec![CigarOp {
                    kind: CigarKind::Match,
                    len: record.sequence.len() as u32,
                }],
            };
            let pair = AlignedPair {
                first: record,
                second: None,
                first_hit: Some(hit),
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes,
            };
            let encoded = write_record_stream(
                reference.as_slice(),
                &[pair],
                AlignmentLayout::SingleEnd,
                &encoding,
                &global,
            )
            .unwrap();
            let decoded = decode_records(
                &encoded,
                1,
                &reference,
                AlignmentLayout::SingleEnd,
                &descriptors,
                &global,
            )
            .unwrap();
            assert_eq!(decoded[0][0], original);
        }
    }

    #[test]
    fn whitelist_decoder_rejects_bad_ordinal_position_padding_and_count() {
        let range = EmbeddedIndexRange {
            mate: EmbeddedIndexMate::R1,
            start: 0,
            end: 3,
        };
        let mut builder = WhitelistBuilder::new();
        builder.add(b"ACG").unwrap();
        builder.add(b"TGC").unwrap();
        builder.add(b"GAT").unwrap();
        let whitelist = builder.finish().unwrap();
        let descriptors = vec![IndexDescriptor { range }];
        let decode = |mode_bits: &[u8], bit_count: usize, op_count: usize| {
            let mut body = Vec::new();
            write_varint(&mut body, op_count as u64).unwrap();
            if op_count != 0 {
                body.push(0xc0);
            }
            write_varint(&mut body, bit_count as u64).unwrap();
            body.extend_from_slice(mode_bits);
            decode_operation_body(
                &mut Cursor::new(body),
                None,
                false,
                &descriptors,
                0,
                &whitelist,
            )
        };
        assert!(decode(&[0b0011_0000], 4, 1).is_err());
        assert!(decode(&[0b0100_1110], 8, 1).is_err());
        assert!(decode(&[0b0000_0001], 4, 1).is_err());
        assert!(decode(&[], 0, 1).is_err());
        assert!(decode(&[0], 1, 1).is_err());
        assert!(decode(&[], 0, 0).is_err());
    }

    #[test]
    fn literal_and_raw_whitelist_sources_share_one_payload_across_ref_advance() {
        let ranges = vec![
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 1,
                end: 3,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 3,
                end: 5,
            },
        ];
        let mut locals = Vec::new();
        for _ in &ranges {
            let mut builder = WhitelistBuilder::new();
            builder.add(b"AA").unwrap();
            locals.push(builder.finish().unwrap());
        }
        let global = index_read::merge_whitelists(&locals).unwrap();
        let encoding = ranges
            .iter()
            .cloned()
            .zip(locals)
            .map(|(range, local_whitelist)| IndexEncodingDescriptor {
                range,
                local_whitelist,
            })
            .collect::<Vec<_>>();
        let descriptors = ranges
            .into_iter()
            .map(|range| IndexDescriptor { range })
            .collect::<Vec<_>>();
        let indexes = vec![
            SpoolIndexRead {
                sequence: b"CC".to_vec(),
                quality: None,
            },
            SpoolIndexRead {
                sequence: b"NN".to_vec(),
                quality: None,
            },
        ];
        let runs = vec![
            ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            },
            ArchiveRun {
                op: ArchiveOp::Whitelist,
                len: 0,
                descriptors: vec![0, 1],
            },
            ArchiveRun {
                op: ArchiveOp::RefAdvance,
                len: 1,
                descriptors: Vec::new(),
            },
            ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            },
        ];
        let hit = GenomeAlignment {
            start: 0,
            query_start: 0,
            query_end: 2,
            cigar: Vec::new(),
            is_reverse: false,
        };
        let mut body = Vec::new();
        encode_operation_body(
            &mut body,
            b"GT",
            Some((&hit, b"A".as_slice())),
            &runs,
            &indexes,
            &encoding,
            &global,
        )
        .unwrap();

        let mut parsed = Cursor::new(body.clone());
        let operation_count = read_varint(&mut parsed).unwrap() as usize;
        assert_eq!(operation_count, 4);
        let mut packed = vec![0; operation_count.div_ceil(4)];
        parsed.read_exact(&mut packed).unwrap();
        for index in 0..operation_count {
            let code = (packed[index / 4] >> (6 - index % 4 * 2)) & 3;
            if ArchiveOp::from_code(code).unwrap() != ArchiveOp::Whitelist {
                read_varint(&mut parsed).unwrap();
            }
        }
        let bit_count = read_varint(&mut parsed).unwrap() as usize;
        assert_eq!(bit_count, 2);
        let mut bit_bytes = vec![0; bit_count.div_ceil(8)];
        parsed.read_exact(&mut bit_bytes).unwrap();
        assert_eq!(bit_bytes, [0b1100_0000]);
        let literal_start = parsed.position() as usize;
        assert_eq!(&parsed.get_ref()[literal_start..], b"GCCNNT");
        parsed.set_position(parsed.get_ref().len() as u64);
        assert_eq!(parsed.position(), parsed.get_ref().len() as u64);

        assert_eq!(
            decode_operation_body(
                &mut Cursor::new(body),
                Some((b"A", 0)),
                false,
                &descriptors,
                0,
                &global,
            )
            .unwrap(),
            b"GCCNNT"
        );
    }

    #[test]
    fn exact_and_mismatch_whitelist_fragments_split_mixed_literal_runs() {
        let ranges = (0..5)
            .map(|index| EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 1 + index * 2,
                end: 3 + index * 2,
            })
            .collect::<Vec<_>>();
        let mut locals = Vec::new();
        for _ in &ranges {
            let mut builder = WhitelistBuilder::new();
            builder.add(b"AA").unwrap();
            locals.push(builder.finish().unwrap());
        }
        let global = index_read::merge_whitelists(&locals).unwrap();
        let encoding = ranges
            .iter()
            .cloned()
            .zip(locals)
            .map(|(range, local_whitelist)| IndexEncodingDescriptor {
                range,
                local_whitelist,
            })
            .collect::<Vec<_>>();
        let descriptors = ranges
            .into_iter()
            .map(|range| IndexDescriptor { range })
            .collect::<Vec<_>>();
        let indexes = [b"CC".as_slice(), b"AA", b"NN", b"AC", b"GG"]
            .into_iter()
            .map(|sequence| SpoolIndexRead {
                sequence: sequence.to_vec(),
                quality: None,
            })
            .collect::<Vec<_>>();
        let runs = vec![
            ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            },
            ArchiveRun {
                op: ArchiveOp::Whitelist,
                len: 0,
                descriptors: vec![0, 1, 2, 3, 4],
            },
            ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            },
        ];
        let mut body = Vec::new();
        encode_operation_body(
            &mut body,
            b"GT",
            None::<(&GenomeAlignment, &[u8])>,
            &runs,
            &indexes,
            &encoding,
            &global,
        )
        .unwrap();

        let mut parsed = Cursor::new(body.clone());
        assert_eq!(read_varint(&mut parsed).unwrap(), 3);
        let mut packed = [0];
        parsed.read_exact(&mut packed).unwrap();
        read_varint(&mut parsed).unwrap();
        read_varint(&mut parsed).unwrap();
        let bit_count = read_varint(&mut parsed).unwrap() as usize;
        assert_eq!(bit_count, 10);
        let mut bit_bytes = vec![0; bit_count.div_ceil(8)];
        parsed.read_exact(&mut bit_bytes).unwrap();
        let literal_start = parsed.position() as usize;
        assert_eq!(&parsed.get_ref()[literal_start..], b"GCCNNGGT");
        parsed.set_position(parsed.get_ref().len() as u64);
        assert_eq!(parsed.position(), parsed.get_ref().len() as u64);

        assert_eq!(
            decode_operation_body(
                &mut Cursor::new(body),
                None,
                false,
                &descriptors,
                0,
                &global,
            )
            .unwrap(),
            b"GCCAANNACGGT"
        );
    }

    #[test]
    fn reverse_raw_whitelist_fragments_use_oriented_literal_bytes() {
        let ranges = vec![
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 1,
                end: 3,
            },
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R1,
                start: 3,
                end: 5,
            },
        ];
        let mut locals = Vec::new();
        for _ in &ranges {
            let mut builder = WhitelistBuilder::new();
            builder.add(b"AA").unwrap();
            locals.push(builder.finish().unwrap());
        }
        let global = index_read::merge_whitelists(&locals).unwrap();
        let encoding = ranges
            .iter()
            .cloned()
            .zip(locals)
            .map(|(range, local_whitelist)| IndexEncodingDescriptor {
                range,
                local_whitelist,
            })
            .collect::<Vec<_>>();
        let descriptors = ranges
            .into_iter()
            .map(|range| IndexDescriptor { range })
            .collect::<Vec<_>>();
        let indexes = vec![
            SpoolIndexRead {
                sequence: b"CC".to_vec(),
                quality: None,
            },
            SpoolIndexRead {
                sequence: b"NN".to_vec(),
                quality: None,
            },
        ];
        let runs = vec![
            ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            },
            ArchiveRun {
                op: ArchiveOp::Whitelist,
                len: 0,
                descriptors: vec![1, 0],
            },
            ArchiveRun {
                op: ArchiveOp::Literal,
                len: 1,
                descriptors: Vec::new(),
            },
        ];
        let hit = GenomeAlignment {
            start: 0,
            query_start: 0,
            query_end: 2,
            cigar: Vec::new(),
            is_reverse: true,
        };
        let mut body = Vec::new();
        encode_operation_body(
            &mut body,
            b"GT",
            Some((&hit, b"".as_slice())),
            &runs,
            &indexes,
            &encoding,
            &global,
        )
        .unwrap();
        let mut parsed = Cursor::new(body.clone());
        let operation_count = read_varint(&mut parsed).unwrap() as usize;
        parsed.set_position(parsed.position() + operation_count.div_ceil(4) as u64);
        for index in 0..operation_count {
            if !matches!(index, 1) {
                read_varint(&mut parsed).unwrap();
            }
        }
        let bit_count = read_varint(&mut parsed).unwrap() as usize;
        parsed.set_position(parsed.position() + bit_count.div_ceil(8) as u64);
        assert_eq!(&parsed.get_ref()[parsed.position() as usize..], b"ANNGGC");
        assert_eq!(
            decode_operation_body(
                &mut Cursor::new(body),
                Some((b"", 0)),
                true,
                &descriptors,
                0,
                &global,
            )
            .unwrap(),
            b"GCCNNT"
        );
    }

    #[test]
    fn chunk_header_is_44_bytes() {
        let mut encoded = Vec::new();
        write_u32(&mut encoded, 2).unwrap();
        for value in [0u64, 0, 0, 0, 0] {
            write_u64(&mut encoded, value).unwrap();
        }
        assert_eq!(encoded.len(), 44);
        let chunk = read_encoded_archive_chunk(&mut Cursor::new(encoded), 7).unwrap();
        assert_eq!(chunk.first_record, 7);
        assert_eq!(chunk.record_count, 2);
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

                let mut sequence = compact.clone();
                if reverse {
                    reverse_complement(&mut sequence);
                }
                let record = ReadRecord::new("read", sequence.clone(), None);
                let mut encoded = Vec::new();
                encode_record_body(&mut encoded, &record, &transformed, compact.as_slice())
                    .unwrap();
                let decoded = decode_mapped_body(
                    &mut Cursor::new(encoded),
                    &compact,
                    transformed.start,
                    reverse,
                )
                .unwrap();
                assert_eq!(decoded, sequence);
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

        for (sequence, hit) in [(b"AAAGGG".as_slice(), skipped), (b"CCCC", overlapping)] {
            let record = ReadRecord::new("read", sequence.to_vec(), None);
            let mut encoded = Vec::new();
            encode_record_body(&mut encoded, &record, &hit, compact.as_slice()).unwrap();
            assert_eq!(
                decode_mapped_body(&mut Cursor::new(encoded), &compact, hit.start, false).unwrap(),
                sequence
            );
        }
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
            };
            let mut encoded = Vec::new();
            write_id_stream(&mut encoded, &[pair], mode, AlignmentLayout::PairedEnd).unwrap();
            let decoded = decode_ids(&encoded, 1, mode, AlignmentLayout::PairedEnd).unwrap();
            assert_eq!(decoded[0][0], expected_first);
            assert_eq!(decoded[1][0], expected_second);
        }
    }

    #[test]
    fn v20_literal_column_is_unframed_raw_source_bytes() {
        let whitelist = index_read::merge_whitelists(&[]).unwrap();
        for sequence in [b"ACGTAC".as_slice(), b"ACNTac".as_slice()] {
            let pair = AlignedPair {
                first: ReadRecord::new("read", sequence.to_vec(), None),
                second: None,
                first_hit: None,
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
            };
            let encoded = write_record_stream(
                b"".as_slice(),
                &[pair],
                AlignmentLayout::SingleEnd,
                &[],
                &whitelist,
            )
            .unwrap();
            assert_eq!(record_literals(&encoded, 1), sequence);
            assert_eq!(
                decode_records(
                    &encoded,
                    1,
                    &[],
                    AlignmentLayout::SingleEnd,
                    &[],
                    &whitelist,
                )
                .unwrap()[0][0],
                sequence
            );
        }
    }

    #[test]
    fn v20_mapped_literal_and_raw_whitelist_bytes_use_alignment_orientation() {
        let range = EmbeddedIndexRange {
            mate: EmbeddedIndexMate::R1,
            start: 1,
            end: 3,
        };
        let mut builder = WhitelistBuilder::new();
        builder.add(b"AA").unwrap();
        let local = builder.finish().unwrap();
        let global = index_read::merge_whitelists(std::slice::from_ref(&local)).unwrap();
        let encoding = vec![IndexEncodingDescriptor {
            range: range.clone(),
            local_whitelist: local,
        }];
        let descriptors = vec![IndexDescriptor { range }];

        for (reverse, reference, expected_literals) in [
            (false, b"A".as_slice(), b"GCCC".as_slice()),
            (true, b"T", b"GGGC"),
        ] {
            let hit = GenomeAlignment {
                start: 0,
                query_start: 0,
                query_end: 3,
                is_reverse: reverse,
                cigar: vec![
                    CigarOp {
                        kind: CigarKind::Ins,
                        len: 1,
                    },
                    CigarOp {
                        kind: CigarKind::Match,
                        len: 1,
                    },
                    CigarOp {
                        kind: CigarKind::Ins,
                        len: 1,
                    },
                ],
            };
            let pair = AlignedPair {
                first: ReadRecord::new("read", b"GAC".to_vec(), None),
                second: None,
                first_hit: Some(hit),
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: vec![SpoolIndexRead {
                    sequence: b"CC".to_vec(),
                    quality: None,
                }],
            };
            let encoded = write_record_stream(
                reference,
                &[pair],
                AlignmentLayout::SingleEnd,
                &encoding,
                &global,
            )
            .unwrap();
            assert_eq!(record_literals(&encoded, 1), expected_literals);
            assert_eq!(
                decode_records(
                    &encoded,
                    1,
                    reference,
                    AlignmentLayout::SingleEnd,
                    &descriptors,
                    &global,
                )
                .unwrap()[0][0],
                b"GCCAC"
            );
        }
    }

    #[test]
    fn v20_decoder_rejects_truncated_and_extra_literal_columns() {
        let whitelist = index_read::merge_whitelists(&[]).unwrap();
        let pair = AlignedPair {
            first: ReadRecord::new("read", b"ACGN".to_vec(), None),
            second: None,
            first_hit: None,
            second_hit: None,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
            indexes: Vec::new(),
        };
        let encoded = write_record_stream(
            b"".as_slice(),
            &[pair],
            AlignmentLayout::SingleEnd,
            &[],
            &whitelist,
        )
        .unwrap();
        let directory = record_directory(&encoded, 1);
        let literal_len_offset = directory.start + 7 * 8;
        let literal_start = encoded.len() - directory.values[7];
        let decode = |bytes: &[u8]| {
            decode_records(bytes, 1, &[], AlignmentLayout::SingleEnd, &[], &whitelist)
        };

        let mut truncated = encoded.clone();
        truncated[literal_len_offset..literal_len_offset + 8].copy_from_slice(&3u64.to_le_bytes());
        truncated.remove(literal_start + 3);
        assert!(matches!(
            decode(&truncated),
            Err(Error::InvalidFormat(message)) if message == "truncated record literal bytes"
        ));

        let mut extra = encoded;
        extra[literal_len_offset..literal_len_offset + 8].copy_from_slice(&5u64.to_le_bytes());
        extra.push(b'X');
        assert!(matches!(
            decode(&extra),
            Err(Error::InvalidFormat(message)) if message == "trailing record literal bytes"
        ));
    }

    #[test]
    fn v20_alignment_cost_counts_raw_literal_bytes() {
        assert_eq!(
            alignment_candidate_core_cost(
                b"ACGTNN",
                None::<(&GenomeAlignment, &[u8])>,
                &[ArchiveRun {
                    op: ArchiveOp::Literal,
                    len: 6,
                    descriptors: Vec::new(),
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
    fn rejects_non_v20_archives() {
        for version in [
            0,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            21,
            u32::MAX,
        ] {
            let mut header = Vec::new();
            header.extend(MAGIC);
            write_u32(&mut header, version).unwrap();
            assert!(matches!(
                read_header(&mut Cursor::new(header)),
                Err(Error::InvalidFormat(message))
                    if message == "unsupported alignment archive version"
            ));
        }
    }

    #[test]
    fn rejects_foreign_alignment_magic() {
        let mut header = vec![0; MAGIC.len()];
        write_u32(&mut header, VERSION).unwrap();
        assert!(matches!(
            read_header(&mut Cursor::new(header)),
            Err(Error::InvalidFormat(message)) if message == "bad alignment archive magic"
        ));
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

    #[test]
    fn record_stream_round_trips_mate_relative_coordinates() {
        let reference =
            b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
        let first_sequence = reference[10..30].to_vec();
        let mut second_sequence = reference[40..60].to_vec();
        reverse_complement(&mut second_sequence);
        let alignment = |start, reverse| GenomeAlignment {
            start,
            query_start: 0,
            query_end: 20,
            cigar: vec![CigarOp {
                kind: CigarKind::Match,
                len: 20,
            }],
            is_reverse: reverse,
        };
        let pairs = vec![
            AlignedPair {
                first: ReadRecord::new("pair0/1", first_sequence.clone(), None),
                second: Some(ReadRecord::new("pair0/2", second_sequence.clone(), None)),
                first_hit: Some(alignment(10, false)),
                second_hit: Some(alignment(40, true)),
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
            },
            AlignedPair {
                first: ReadRecord::new("pair1/1", b"ACGTNCGT".to_vec(), None),
                second: Some(ReadRecord::new("pair1/2", b"NNNNNNNN".to_vec(), None)),
                first_hit: None,
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
                indexes: Vec::new(),
            },
        ];
        let reference_slice = reference.as_slice();
        let whitelist = index_read::merge_whitelists(&[]).unwrap();
        let encoded = write_record_stream(
            reference_slice,
            &pairs,
            AlignmentLayout::PairedEnd,
            &[],
            &whitelist,
        )
        .unwrap();
        let decoded = decode_records(
            &encoded,
            pairs.len(),
            reference,
            AlignmentLayout::PairedEnd,
            &[],
            &whitelist,
        )
        .unwrap();
        assert_eq!(decoded[0][0], first_sequence);
        assert_eq!(decoded[1][0], second_sequence);
        assert_eq!(decoded[0][1], b"ACGTNCGT");
        assert_eq!(decoded[1][1], b"NNNNNNNN");
    }
}
