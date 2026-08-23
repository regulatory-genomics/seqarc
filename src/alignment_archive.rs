use crate::minibwa;
use crate::minibwa::{Aligner, CigarKind, CigarOp, Index, Opts, Strand, ThreadBuf};
use crate::packed_sequence::PackedSequence;
use crate::quality_blocks::{decode_stream_bytes, encode_stream_bytes};
use crate::{
    apply_quality, gzip_reader, AlignedCompressionProfile, AlignedDecompressionMetadata,
    AlignedRecordGroup, AlignedRecordSink, AlignmentLayout, ArchiveRecords, CompressionMetadata,
    DecompressionOptions, Error, FastxReader, Format, GenomeAlignment, PairedRecordSink,
    ReadRecord, Result,
};
use bed_utils::extsort::{ExternalChunk, ExternalChunkBuilder};
use indicatif::{ProgressBar, ProgressStyle};
use rkyv::{Archive, Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom, Write};
use std::path::Path;

const MAGIC: &[u8; 9] = b"SEQARC-AL";
const VERSION: u32 = 11;
const BATCH_SIZE: usize = 256;
const PARSE_BATCH_PAIRS: usize = 4096;

#[derive(Clone, Debug)]
struct AlignedPair {
    first: ReadRecord,
    second: Option<ReadRecord>,
    first_hit: Option<GenomeAlignment>,
    second_hit: Option<GenomeAlignment>,
    first_selection: AlignmentSelection,
    second_selection: AlignmentSelection,
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

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolPair {
    first: SpoolRead,
    second: Option<SpoolRead>,
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

fn compress_with_index<I2: AsRef<Path>>(
    index: &Index,
    input_1: impl AsRef<Path>,
    input_2: Option<I2>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<(CompressionMetadata, AlignedCompressionProfile)> {
    let pipeline_started = std::time::Instant::now();
    let mapping_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads.max(1))
        .build()
        .map_err(|error| Error::InvalidInput(format!("failed to create mapping pool: {error}")))?;
    let progress = mapping_progress();
    stage_log(&progress, format!(
        "[seqarc] stage=initialize layout={:?} streams={} threads={} chunk_records={} zstd_level={} preserve_ids={} preserve_quality={}",
        layout,
        layout.streams(),
        options.threads.max(1),
        options.chunk_records.max(1),
        options.zstd_level,
        options.preserve_ids,
        options.preserve_quality
    ));
    let temp_dir =
        crate::temporary_directory(options.temp_dir.as_deref(), "seqarc-aligned-spool-")?;
    let spool_path = temp_dir.path().join("aligned-pairs.lz4");
    let mut spool_builder = ExternalChunkBuilder::new(File::create(&spool_path)?, 1)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let mut mask = ReferenceMask::new(index.reference_len());
    let mut id_detector = IdModeDetector::new(options.preserve_ids, layout);
    let mut profile = AlignedCompressionProfile::default();
    let input_1 = input_1.as_ref().to_path_buf();
    let input_2 = input_2.map(|path| path.as_ref().to_path_buf());
    let mapper_options = options.clone();
    let (parse_tx, parse_rx) = crossbeam_channel::bounded::<Result<Vec<AlignedPair>>>(4);
    let (mapped_tx, mapped_rx) = crossbeam_channel::bounded::<Result<Vec<AlignedPair>>>(4);
    stage_log(&progress, "[seqarc] stage=map_and_spool start");
    let map_started = std::time::Instant::now();
    let result = std::thread::scope(|scope| -> Result<()> {
        let producer_tx = parse_tx.clone();
        let producer = scope.spawn(move || -> Result<()> {
            let mut first = FastxReader::new(
                BufReader::new(gzip_reader(File::open(input_1)?)),
                Format::Fastq,
            );
            let mut pairs = Vec::with_capacity(PARSE_BATCH_PAIRS);
            if let Some(input_2) = input_2 {
                let mut second = FastxReader::new(
                    BufReader::new(gzip_reader(File::open(input_2)?)),
                    Format::Fastq,
                );
                loop {
                    match (first.next(), second.next()) {
                        (None, None) => break,
                        (Some(Ok(mut first)), Some(Ok(mut second))) => {
                            first.validate(true)?;
                            second.validate(true)?;
                            if mapper_options.preserve_quality {
                                apply_quality(&mut first, &mapper_options.quality_mode);
                                apply_quality(&mut second, &mapper_options.quality_mode);
                            }
                            pairs.push(AlignedPair {
                                first,
                                second: Some(second),
                                first_hit: None,
                                second_hit: None,
                                first_selection: AlignmentSelection::default(),
                                second_selection: AlignmentSelection::default(),
                            });
                        }
                        (None, Some(_)) | (Some(_), None) => {
                            return Err(Error::InvalidInput(
                                "paired streams have different record counts".into(),
                            ));
                        }
                        (Some(Err(error)), _) | (_, Some(Err(error))) => return Err(error),
                    }
                    if pairs.len() == PARSE_BATCH_PAIRS {
                        producer_tx
                            .send(Ok(std::mem::take(&mut pairs)))
                            .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
                        pairs = Vec::with_capacity(PARSE_BATCH_PAIRS);
                    }
                }
            } else {
                for first in first {
                    let mut first = first?;
                    first.validate(true)?;
                    if mapper_options.preserve_quality {
                        apply_quality(&mut first, &mapper_options.quality_mode);
                    }
                    pairs.push(AlignedPair {
                        first,
                        second: None,
                        first_hit: None,
                        second_hit: None,
                        first_selection: AlignmentSelection::default(),
                        second_selection: AlignmentSelection::default(),
                    });
                    if pairs.len() == PARSE_BATCH_PAIRS {
                        producer_tx
                            .send(Ok(std::mem::take(&mut pairs)))
                            .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
                        pairs = Vec::with_capacity(PARSE_BATCH_PAIRS);
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
            while let Ok(batch) = parse_rx.recv() {
                let mapped = batch.and_then(|mut pairs| {
                    map_pairs(index, &mut pairs, &mapping_pool, &mapper_progress)?;
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
    let mut archive = File::create(output)?;
    let header_size = 9 + 4 + 2 + 8 * 2;
    archive.write_all(MAGIC)?;
    write_u32(&mut archive, VERSION)?;
    let mut flags = 0u16;
    if options.preserve_ids {
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
    let chunk_count_pos = archive.stream_position()?;
    write_u64(&mut archive, 0)?;
    debug_assert_eq!(archive.stream_position()?, header_size as u64);
    archive.write_all(&reference_block.block)?;

    let mut written_chunks = 0u64;
    let mut replay_pair_count = 0u64;
    progress.set_message("compressing chunks".to_owned());
    let worker_count = options.threads.max(1);
    stage_log(
        &progress,
        format!(
            "[seqarc] stage=chunks start workers={} chunk_records={}",
            worker_count,
            options.chunk_records.max(1)
        ),
    );
    let chunks_started = std::time::Instant::now();
    let (job_tx, job_rx) =
        crossbeam_channel::bounded::<Result<(usize, Vec<AlignedPair>)>>(worker_count);
    let (result_tx, result_rx) =
        crossbeam_channel::bounded::<Result<(usize, EncodedChunk)>>(worker_count);
    let reference_for_workers = compact_reference.as_slice();
    std::thread::scope(|scope| -> Result<()> {
        let dispatcher = scope.spawn({
            let job_tx = job_tx.clone();
            let spool_path = spool_path.clone();
            let kept_runs = kept_runs.clone();
            move || -> Result<()> {
                let spool = ExternalChunk::<SpoolPair>::open(File::open(spool_path)?)
                    .map_err(|error| Error::InvalidFormat(error.to_string()))?;
                let mut pairs = Vec::with_capacity(options.chunk_records.max(1));
                let mut chunk_index = 0usize;
                for item in spool {
                    let item = item.map_err(|error| Error::InvalidFormat(error.to_string()))?;
                    pairs.push(item.into_aligned_pair(&kept_runs)?);
                    if pairs.len() == options.chunk_records.max(1) {
                        job_tx
                            .send(Ok((chunk_index, std::mem::take(&mut pairs))))
                            .map_err(|_| Error::InvalidInput("chunk workers exited".into()))?;
                        chunk_index += 1;
                        pairs = Vec::with_capacity(options.chunk_records.max(1));
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
                        encode_chunk(reference_for_workers, &pairs, options, id_mode, layout)
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
        profile.ids_raw_bytes + profile.qualities_raw_bytes + profile.records_raw_bytes,
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
            preserved_ids: options.preserve_ids,
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
    let mut streams = vec![Vec::new(), Vec::new()];
    decode_aligned(input, options, |_, group| match group {
        AlignedRecordGroup::Paired(first, second) => {
            streams[0].push(first);
            streams[1].push(second);
            Ok(())
        }
        AlignedRecordGroup::Single(_) => Err(Error::InvalidFormat(
            "single-end archive passed to paired decoder".into(),
        )),
    })?;
    Ok(ArchiveRecords { streams })
}

pub(crate) fn decompress_single<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let mut stream = Vec::new();
    decode_aligned(input, options, |_, group| match group {
        AlignedRecordGroup::Single(record) => {
            stream.push(record);
            Ok(())
        }
        AlignedRecordGroup::Paired(_, _) => Err(Error::InvalidFormat(
            "paired archive passed to single-end decoder".into(),
        )),
    })?;
    Ok(ArchiveRecords {
        streams: vec![stream],
    })
}

pub(crate) fn decompress_paired_into<R: Read, S: PairedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<crate::DecompressionMetadata> {
    let mut pairs = 0;
    decode_aligned(input, options, |original_id, group| match group {
        AlignedRecordGroup::Paired(first, second) => {
            sink.write_pair(original_id, first, second)?;
            pairs += 1;
            Ok(())
        }
        AlignedRecordGroup::Single(_) => Err(Error::InvalidFormat(
            "single-end archive passed to paired decoder".into(),
        )),
    })?;
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
    let mut records = 0;
    let layout = decode_aligned(input, options, |original_id, group| {
        sink.write_group(original_id, group)?;
        records += 1;
        Ok(())
    })?;
    sink.finish()?;
    Ok(AlignedDecompressionMetadata {
        layout,
        streams: layout.streams(),
        records,
        output_order: crate::OutputOrder::Original,
    })
}

fn decode_aligned<R: Read, F: FnMut(u64, AlignedRecordGroup) -> Result<()>>(
    mut input: R,
    options: &DecompressionOptions,
    mut emit: F,
) -> Result<AlignmentLayout> {
    let (flags, id_mode, layout, chunk_count, reference) = read_header(&mut input)?;
    let preserve_quality = flags & 2 != 0;
    let mut first_pair = 0u64;
    for _ in 0..chunk_count {
        let pair_count = read_u32(&mut input)? as usize;
        let compressed_len = usize_from(read_u64(&mut input)?, "chunk compressed size")?;
        let raw_len = usize_from(read_u64(&mut input)?, "chunk uncompressed size")?;
        let ids_len = usize_from(read_u64(&mut input)?, "chunk ID block size")?;
        let qualities_len = usize_from(read_u64(&mut input)?, "chunk quality block size")?;
        let records_len = usize_from(read_u64(&mut input)?, "chunk record block size")?;
        let mut compressed = vec![0; compressed_len];
        input.read_exact(&mut compressed)?;
        let raw = decode_chunk_sections(&compressed, ids_len, qualities_len, records_len)?;
        if raw.len() != raw_len || ids_len + qualities_len + records_len != raw.len() {
            return Err(Error::InvalidFormat(
                "chunk block sizes do not match".into(),
            ));
        }
        let ids_end = ids_len;
        let qualities_end = ids_end + qualities_len;
        let ids = decode_ids(&raw[..ids_end], pair_count, id_mode, layout)?;
        let qualities = decode_qualities(
            &raw[ids_end..qualities_end],
            pair_count,
            preserve_quality,
            layout,
        )?;
        let records = &raw[qualities_end..];
        let decoded = decode_records(records, pair_count, first_pair, &reference, layout)?;
        for index in 0..pair_count {
            let original_id = decoded[index].0;
            if options
                .original_range
                .as_ref()
                .is_some_and(|range| !range.contains(&original_id))
            {
                continue;
            }
            let first = ReadRecord::new(
                ids[0][index].clone(),
                decoded[index].1[0].clone(),
                qualities[0][index].clone(),
            );
            let group = if layout == AlignmentLayout::SingleEnd {
                AlignedRecordGroup::Single(first)
            } else {
                AlignedRecordGroup::Paired(
                    first,
                    ReadRecord::new(
                        ids[1][index].clone(),
                        decoded[index].1[1].clone(),
                        qualities[1][index].clone(),
                    ),
                )
            };
            emit(original_id, group)?;
        }
        first_pair = first_pair
            .checked_add(pair_count as u64)
            .ok_or_else(|| Error::InvalidFormat("pair count overflow".into()))?;
    }
    let mut trailing = [0; 1];
    if input.read(&mut trailing)? != 0 {
        return Err(Error::InvalidFormat(
            "trailing aligned archive bytes".into(),
        ));
    }
    Ok(layout)
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
) -> Result<()> {
    let first_hits = map_mate(index, pairs, 0, mapping_pool, progress)?;
    for (pair, first) in pairs.iter_mut().zip(first_hits) {
        pair.first_hit = first.hit;
        pair.first_selection = first.selection;
    }
    if pairs.first().is_some_and(|pair| pair.second.is_some()) {
        let second_hits = map_mate(index, pairs, 1, mapping_pool, progress)?;
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
    progress: &ProgressBar,
) -> Result<Vec<SelectedAlignment>> {
    let opts = Opts::new().set_paired(false);
    let aligner = Aligner::new(index, &opts);
    mapping_pool.install(|| {
        use rayon::prelude::*;
        pairs
            .par_chunks(BATCH_SIZE)
            .map_init(ThreadBuf::new, |buffer, batch| {
                let records = batch.iter().map(|pair| pair.record(mate));
                let result = align_batch(&aligner, buffer, index, records);
                progress.inc(batch.len() as u64);
                result
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
                CigarKind::Match
                | CigarKind::Eq
                | CigarKind::Diff
                | CigarKind::Del
                | CigarKind::RefSkip => {
                    let end = position.checked_add(length).ok_or_else(|| {
                        Error::InvalidInput("reference mask position overflow".into())
                    })?;
                    self.mark_range(position, end)?;
                    position = end;
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

fn compact_position(runs: &[KeptRun], original: usize) -> Result<usize> {
    let index = runs.partition_point(|run| run.original_start <= original);
    let run = index
        .checked_sub(1)
        .and_then(|index| runs.get(index))
        .filter(|run| original < run.original_start + run.length)
        .ok_or_else(|| Error::InvalidFormat("alignment starts in pruned reference".into()))?;
    Ok(run.compact_start + (original - run.original_start))
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
        })
    }

    fn into_aligned_pair(self, runs: &[KeptRun]) -> Result<AlignedPair> {
        Ok(AlignedPair {
            first: self.first.record()?,
            second: self.second.as_ref().map(SpoolRead::record).transpose()?,
            first_hit: self
                .first
                .alignment
                .map(|hit| hit.into_genome(runs))
                .transpose()?,
            second_hit: self
                .second
                .and_then(|read| read.alignment)
                .map(|hit| hit.into_genome(runs))
                .transpose()?,
            first_selection: AlignmentSelection::default(),
            second_selection: AlignmentSelection::default(),
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

    fn into_genome(self, runs: &[KeptRun]) -> Result<GenomeAlignment> {
        let start = compact_position(runs, nonnegative_usize(self.start, "alignment start")?)?;
        let cigar = self
            .cigar
            .into_iter()
            .map(|op| {
                Ok(CigarOp {
                    kind: cigar_kind(op.op)?,
                    len: op.len,
                })
            })
            .collect::<Result<Vec<_>>>()?;
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
) -> Result<EncodedChunk> {
    let mut ids = Vec::new();
    write_id_stream(&mut ids, pairs, id_mode, layout)?;
    let mut qualities = Vec::new();
    write_quality_stream(&mut qualities, pairs, options, layout)?;
    let records = write_record_stream(reference, pairs, layout)?;
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
    options: &crate::AlignmentCompressionOptions,
    layout: AlignmentLayout,
) -> Result<()> {
    for mate in 0..layout.streams() {
        let values = pairs
            .iter()
            .map(|pair| pair.record(mate).quality.as_deref().unwrap_or(&[]))
            .collect::<Vec<_>>();
        let block = if options.preserve_quality {
            encode_stream_bytes(&values)?
        } else {
            Vec::new()
        };
        write_u64(output, block.len() as u64)?;
        output.extend(block);
    }
    Ok(())
}

fn write_record_stream<R: ReferenceBases + ?Sized>(
    reference: &R,
    pairs: &[AlignedPair],
    layout: AlignmentLayout,
) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    let mut starts = [0i64; 2];
    for mate in 0..layout.streams() {
        let bitmap_len = pairs.len().div_ceil(8);
        let mut mapped_bitmap = vec![0; bitmap_len];
        let mut reverse_bitmap = vec![0; bitmap_len];
        let mut positions = Vec::new();
        let mut bodies = Vec::new();
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
                encode_record_body(&mut bodies, record, hit, reference)?;
            } else {
                write_packed(&mut bodies, &record.sequence)?;
            }
        }
        output.extend(mapped_bitmap);
        output.extend(reverse_bitmap);
        write_u64(&mut output, positions.len() as u64)?;
        output.extend(positions);
        output.extend(bodies);
    }
    Ok(output)
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
    let raw_cost = encode_packed(&record.sequence)?.len();
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
        let mut encoded = Vec::new();
        if encode_record(&mut encoded, record, Some(&alignment), reference, &mut 0).is_err() {
            continue;
        }
        let candidate = (encoded.len(), hit.score, index, alignment);
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

fn encode_record<R: ReferenceBases + ?Sized>(
    output: &mut Vec<u8>,
    record: &ReadRecord,
    hit: Option<&GenomeAlignment>,
    reference: &R,
    previous_start: &mut i64,
) -> Result<()> {
    if let Some(hit) = hit {
        let delta = hit
            .start
            .checked_sub(*previous_start)
            .ok_or_else(|| Error::InvalidInput("alignment offset overflow".into()))?;
        *previous_start = hit.start;
        write_varint(output, zigzag_encode(delta))?;
        encode_record_body(output, record, hit, reference)?;
    } else {
        write_packed(output, &record.sequence)?;
    }
    Ok(())
}

fn encode_record_body<R: ReferenceBases + ?Sized>(
    output: &mut Vec<u8>,
    record: &ReadRecord,
    hit: &GenomeAlignment,
    reference: &R,
) -> Result<()> {
    let cigar = effective_cigar(hit, record.sequence.len())?;
    write_varint(output, cigar.len() as u64)?;
    for chunk in cigar.chunks(2) {
        let first = cigar_code(chunk[0].kind)? << 4;
        let second = if chunk.len() == 2 {
            cigar_code(chunk[1].kind)?
        } else {
            0
        };
        write_u8(output, first | second)?;
    }
    for op in &cigar {
        write_varint(output, op.len as u64)?;
    }
    encode_noise(output, &record.sequence, hit, &cigar, reference)
}

fn decode_records(
    input: &[u8],
    count: usize,
    first_pair: u64,
    reference: &[u8],
    layout: AlignmentLayout,
) -> Result<Vec<(u64, Vec<Vec<u8>>)>> {
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
        let positions_len = usize_from(read_u64(&mut cursor)?, "alignment positions")?;
        let mut positions = vec![0; positions_len];
        cursor.read_exact(&mut positions)?;
        let mut positions = Cursor::new(positions);
        for ordinal in 0..count {
            let mapped = mapped_bitmap[ordinal / 8] & (1 << (ordinal % 8)) != 0;
            let encoded_reverse = reverse_bitmap[ordinal / 8] & (1 << (ordinal % 8)) != 0;
            let sequence = if mapped {
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
                decode_mapped_body(&mut cursor, reference, start, reverse)?
            } else {
                read_packed(&mut cursor)?.decode()
            };
            all[mate].push(sequence);
        }
        if positions.position() != positions.get_ref().len() as u64 {
            return Err(Error::InvalidFormat(
                "trailing alignment position bytes".into(),
            ));
        }
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat("trailing record bytes".into()));
    }
    Ok((0..count)
        .map(|index| {
            (
                first_pair + index as u64,
                all.iter().map(|stream| stream[index].clone()).collect(),
            )
        })
        .collect())
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
        let length = usize_from(read_u64(&mut cursor)?, "quality block")?;
        let mut block = vec![0; length];
        cursor.read_exact(&mut block)?;
        if preserve {
            result.push(
                decode_stream_bytes(&block, count)?
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
        return Err(Error::InvalidFormat("trailing quality block bytes".into()));
    }
    Ok(result)
}

fn decode_mapped_body(
    input: &mut Cursor<Vec<u8>>,
    reference: &[u8],
    start: i64,
    is_reverse: bool,
) -> Result<Vec<u8>> {
    let count = usize_from(read_varint(input)?, "CIGAR count")?;
    let mut op_codes = Vec::with_capacity(count);
    for index in 0..count.div_ceil(2) {
        let packed = read_u8(input)?;
        op_codes.push(packed >> 4);
        if index * 2 + 1 < count {
            op_codes.push(packed & 0x0f);
        }
    }
    let mut cigar = Vec::with_capacity(count);
    for op in op_codes {
        cigar.push(CigarOp {
            kind: cigar_kind(op)?,
            len: u32::try_from(read_varint(input)?)
                .map_err(|_| Error::InvalidFormat("CIGAR length too large".into()))?,
        });
    }
    let mut sequence = reconstruct_with_noise(
        input,
        reference,
        nonnegative_usize(start, "alignment start")?,
        &cigar,
    )?;
    if is_reverse {
        reverse_complement(&mut sequence);
    }
    Ok(sequence)
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
    Ok(cigar)
}

fn encode_noise<R: ReferenceBases + ?Sized>(
    output: &mut Vec<u8>,
    sequence: &[u8],
    hit: &GenomeAlignment,
    cigar: &[CigarOp],
    reference: &R,
) -> Result<()> {
    let mut oriented = sequence.to_vec();
    if hit.is_reverse {
        reverse_complement(&mut oriented);
    }
    let mut ref_pos = nonnegative_usize(hit.start, "alignment start")?;
    let mut query_pos = 0usize;
    for op in cigar {
        let len = op.len as usize;
        match op.kind {
            CigarKind::Match | CigarKind::Eq | CigarKind::Diff => {
                let ref_end = ref_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                let query_end = query_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if ref_end > reference.len() || query_end > oriented.len() {
                    return Err(Error::InvalidFormat("CIGAR exceeds bounds".into()));
                }
                let mut mismatches = Vec::new();
                for offset in 0..len {
                    let reference_base = reference.base(ref_pos + offset);
                    if reference_base != oriented[query_pos + offset] {
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
            CigarKind::Ins => {
                let end = query_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if end > oriented.len() {
                    return Err(Error::InvalidFormat("CIGAR exceeds query bounds".into()));
                }
                output.extend(&oriented[query_pos..end]);
                query_pos = end;
            }
            CigarKind::SoftClip => {
                let end = query_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if end > oriented.len() {
                    return Err(Error::InvalidFormat("CIGAR exceeds query bounds".into()));
                }
                write_softclip(output, &oriented[query_pos..end])?;
                query_pos = end;
            }
            CigarKind::Del | CigarKind::RefSkip => {
                ref_pos = ref_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "CIGAR exceeds reference bounds".into(),
                    ));
                }
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }
    if query_pos != oriented.len() {
        return Err(Error::InvalidFormat("CIGAR does not consume query".into()));
    }
    Ok(())
}

fn reconstruct_with_noise(
    input: &mut Cursor<Vec<u8>>,
    reference: &[u8],
    start: usize,
    cigar: &[CigarOp],
) -> Result<Vec<u8>> {
    let mut sequence = Vec::new();
    let mut ref_pos = start;
    for op in cigar {
        let len = op.len as usize;
        match op.kind {
            CigarKind::Match | CigarKind::Eq | CigarKind::Diff => {
                let end = ref_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if end > reference.len() {
                    return Err(Error::InvalidFormat(
                        "CIGAR exceeds reference bounds".into(),
                    ));
                }
                let sequence_start = sequence.len();
                sequence.extend_from_slice(&reference[ref_pos..end]);
                let mismatch_count = usize_from(read_varint(input)?, "mismatch count")?;
                let mut previous = 0usize;
                for _ in 0..mismatch_count {
                    let position = previous
                        .checked_add(usize_from(read_varint(input)?, "mismatch position")?)
                        .ok_or_else(|| Error::InvalidFormat("mismatch position overflow".into()))?;
                    if position >= len {
                        return Err(Error::InvalidFormat(
                            "mismatch position out of bounds".into(),
                        ));
                    }
                    sequence[sequence_start + position] = read_u8(input)?;
                    previous = position;
                }
                ref_pos = end;
            }
            CigarKind::Ins => {
                let end = sequence.len().checked_add(len).ok_or_else(|| {
                    Error::InvalidFormat("inserted sequence length overflow".into())
                })?;
                let mut bases = vec![0; len];
                input.read_exact(&mut bases)?;
                sequence.extend_from_slice(&bases);
                if sequence.len() != end {
                    return Err(Error::InvalidFormat(
                        "inserted sequence length mismatch".into(),
                    ));
                }
            }
            CigarKind::SoftClip => {
                let bases = read_softclip(input, len)?;
                sequence.extend_from_slice(&bases);
            }
            CigarKind::Del | CigarKind::RefSkip => {
                ref_pos = ref_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "CIGAR exceeds reference bounds".into(),
                    ));
                }
            }
            CigarKind::HardClip | CigarKind::Pad => {}
        }
    }
    Ok(sequence)
}

fn read_header<R: Read>(input: &mut R) -> Result<(u16, IdMode, AlignmentLayout, u64, Vec<u8>)> {
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
    let chunks = read_u64(input)?;
    let mut reference_block = vec![0; reference_len];
    input.read_exact(&mut reference_block)?;
    Ok((
        flags,
        id_mode,
        layout,
        chunks,
        decode_reference(&reference_block)?,
    ))
}

fn encode_reference(
    reference: &[u8],
    options: &crate::AlignmentCompressionOptions,
) -> Result<EncodedReference> {
    let mut packed = vec![0; reference.len().div_ceil(4)];
    let mut runs: Vec<(usize, Vec<u8>)> = Vec::new();
    for (index, &base) in reference.iter().enumerate() {
        let bits = match base {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        };
        if let Some(bits) = bits {
            packed[index / 4] |= bits << (6 - (index % 4) * 2);
        } else if let Some((start, bytes)) = runs.last_mut() {
            if *start + bytes.len() == index {
                bytes.push(base);
                continue;
            }
            runs.push((index, vec![base]));
        } else {
            runs.push((index, vec![base]));
        }
    }

    let mut raw = Vec::new();
    write_u64(&mut raw, reference.len() as u64)?;
    write_u64(&mut raw, packed.len() as u64)?;
    write_u64(&mut raw, runs.len() as u64)?;
    raw.extend(packed);
    let mut previous_end = 0usize;
    for (start, bytes) in runs {
        let length = bytes.len();
        write_u64(
            &mut raw,
            start
                .checked_sub(previous_end)
                .ok_or_else(|| Error::InvalidInput("reference run order overflow".into()))?
                as u64,
        )?;
        write_u64(&mut raw, length as u64)?;
        raw.extend(bytes);
        previous_end = start + length;
    }

    let compressed = zstd::stream::encode_all(Cursor::new(&raw), options.zstd_level)?;
    let mut block = Vec::new();
    write_u64(&mut block, raw.len() as u64)?;
    block.extend(compressed);
    Ok(EncodedReference {
        block,
        raw_len: raw.len() as u64,
    })
}

fn decode_reference(block: &[u8]) -> Result<Vec<u8>> {
    let mut input = Cursor::new(block.to_vec());
    let raw_len = usize_from(read_u64(&mut input)?, "reference payload")?;
    let mut compressed = Vec::new();
    input.read_to_end(&mut compressed)?;
    let raw = zstd::stream::decode_all(Cursor::new(compressed))?;
    if raw.len() != raw_len {
        return Err(Error::InvalidFormat(
            "reference payload length mismatch".into(),
        ));
    }
    let mut input = Cursor::new(raw);
    let base_count = usize_from(read_u64(&mut input)?, "reference base count")?;
    let packed_len = usize_from(read_u64(&mut input)?, "reference packed length")?;
    let run_count = usize_from(read_u64(&mut input)?, "reference ambiguity runs")?;
    if packed_len != base_count.div_ceil(4) {
        return Err(Error::InvalidFormat(
            "invalid reference packed length".into(),
        ));
    }
    let mut packed = vec![0; packed_len];
    input.read_exact(&mut packed)?;
    let mut reference = (0..base_count)
        .map(|index| b"ACGT"[((packed[index / 4] >> (6 - (index % 4) * 2)) & 3) as usize])
        .collect::<Vec<_>>();
    let mut previous_end = 0usize;
    for _ in 0..run_count {
        let start = previous_end
            .checked_add(usize_from(read_u64(&mut input)?, "reference run offset")?)
            .ok_or_else(|| Error::InvalidFormat("reference run offset overflow".into()))?;
        let length = usize_from(read_u64(&mut input)?, "reference run length")?;
        let end = start
            .checked_add(length)
            .ok_or_else(|| Error::InvalidFormat("reference run length overflow".into()))?;
        if end > reference.len() {
            return Err(Error::InvalidFormat(
                "reference ambiguity run out of bounds".into(),
            ));
        }
        input.read_exact(&mut reference[start..end])?;
        previous_end = end;
    }
    if input.position() != input.get_ref().len() as u64 {
        return Err(Error::InvalidFormat(
            "trailing reference payload bytes".into(),
        ));
    }
    Ok(reference)
}

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
fn write_packed(output: &mut Vec<u8>, sequence: &[u8]) -> Result<()> {
    output.extend(encode_packed(sequence)?);
    Ok(())
}
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

fn write_softclip(output: &mut Vec<u8>, sequence: &[u8]) -> Result<()> {
    let packed = PackedSequence::try_from_bytes(sequence)?;
    if packed.is_packed() {
        write_u8(output, 0)?;
        output.extend(packed.encoded_bytes());
    } else {
        let sparse = encode_sparse_body(sequence)?;
        if sparse.len() < sequence.len() {
            write_u8(output, 2)?;
            output.extend(sparse);
        } else {
            write_u8(output, 1)?;
            output.extend(sequence);
        }
    }
    Ok(())
}

fn read_softclip(input: &mut Cursor<Vec<u8>>, length: usize) -> Result<Vec<u8>> {
    match read_u8(input)? {
        0 => {
            let body_len = length.div_ceil(4);
            let mut body = vec![0; body_len];
            input.read_exact(&mut body)?;
            PackedSequence::from_encoded(length as u32, true, body).map(|value| value.decode())
        }
        1 => {
            let mut sequence = vec![0; length];
            input.read_exact(&mut sequence)?;
            Ok(sequence)
        }
        2 => {
            let packed_len = length.div_ceil(4);
            let mut body = vec![0; packed_len];
            input.read_exact(&mut body)?;
            let count = usize_from(read_varint(input)?, "sparse sequence exceptions")?;
            write_varint(&mut body, count as u64)?;
            for _ in 0..count {
                write_varint(&mut body, read_varint(input)?)?;
                body.push(read_u8(input)?);
            }
            decode_sparse_body(length, &body)
        }
        _ => Err(Error::InvalidFormat("unknown soft-clip encoding".into())),
    }
}

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
            "{spinner:.green} mapping reads [{elapsed_precise}] {human_pos} reads @ {per_sec}",
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

    #[test]
    fn reference_block_round_trips_ambiguous_bases() {
        let reference = b"ACGTNNNNacgtRYACGT";
        let encoded =
            encode_reference(reference, &crate::AlignmentCompressionOptions::default()).unwrap();
        assert_eq!(decode_reference(&encoded.block).unwrap(), reference);
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
            };
            let mut encoded = Vec::new();
            write_id_stream(&mut encoded, &[pair], mode, AlignmentLayout::PairedEnd).unwrap();
            let decoded = decode_ids(&encoded, 1, mode, AlignmentLayout::PairedEnd).unwrap();
            assert_eq!(decoded[0][0], expected_first);
            assert_eq!(decoded[1][0], expected_second);
        }
    }

    #[test]
    fn softclip_payload_uses_two_bit_encoding_with_raw_fallback() {
        for sequence in [b"ACGTAC".as_slice(), b"ACNTac".as_slice()] {
            let mut encoded = Vec::new();
            write_softclip(&mut encoded, sequence).unwrap();
            let decoded = read_softclip(&mut Cursor::new(encoded), sequence.len()).unwrap();
            assert_eq!(decoded, sequence);
        }
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
    fn rejects_pre_v11_archives() {
        for version in [7, 8, 9, 10] {
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
    fn v11_record_stream_round_trips_mate_relative_coordinates() {
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
            },
            AlignedPair {
                first: ReadRecord::new("pair1/1", b"ACGTNCGT".to_vec(), None),
                second: Some(ReadRecord::new("pair1/2", b"NNNNNNNN".to_vec(), None)),
                first_hit: None,
                second_hit: None,
                first_selection: AlignmentSelection::default(),
                second_selection: AlignmentSelection::default(),
            },
        ];
        let reference_slice = reference.as_slice();
        let encoded =
            write_record_stream(reference_slice, &pairs, AlignmentLayout::PairedEnd).unwrap();
        let decoded = decode_records(
            &encoded,
            pairs.len(),
            0,
            reference,
            AlignmentLayout::PairedEnd,
        )
        .unwrap();
        assert_eq!(decoded[0].1[0], first_sequence);
        assert_eq!(decoded[0].1[1], second_sequence);
        assert_eq!(decoded[1].1[0], b"ACGTNCGT");
        assert_eq!(decoded[1].1[1], b"NNNNNNNN");
    }
}
