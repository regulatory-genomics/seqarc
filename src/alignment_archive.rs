use crate::packed_sequence::PackedSequence;
use crate::quality_blocks::{decode_stream_bytes, encode_stream_bytes};
use crate::{
    apply_quality, gzip_reader, ArchiveRecords, CompressionMetadata, DecompressionOptions, Error,
    FastxReader, Format, GenomeAlignment, PairedRecordSink, ReadRecord, Result,
};
use bed_utils::extsort::{ExternalChunk, ExternalChunkBuilder};
use indicatif::ProgressBar;
use minibwa::{CigarOp, CigarOpKind, Index, Options, Query};
use rkyv::{Archive, Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom, Write};
use std::path::Path;

const MAGIC: &[u8; 9] = b"SPRNG-AL2";
const VERSION: u32 = 6;
const CHUNK_HEADER_SIZE: u64 = 44;
const BATCH_SIZE: usize = 256;
const PARSE_BATCH_PAIRS: usize = 4096;

#[derive(Clone, Debug)]
struct AlignedPair {
    first: ReadRecord,
    second: ReadRecord,
    first_hit: Option<GenomeAlignment>,
    second_hit: Option<GenomeAlignment>,
}

#[derive(Archive, Serialize, Deserialize, Clone, Debug)]
struct SpoolPair {
    first: SpoolRead,
    second: SpoolRead,
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

#[derive(Clone, Copy, Debug)]
struct ChunkIndexEntry {
    offset: u64,
    first_pair: u64,
    pairs: u32,
    compressed_size: u64,
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
}

impl IdModeDetector {
    fn new(preserve: bool) -> Self {
        Self {
            mode: None,
            preserve,
        }
    }

    fn observe(&mut self, first: &[u8], second: &[u8]) {
        if !self.preserve {
            return;
        }
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
            self.mode.unwrap_or(IdMode::IdenticalMates)
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
) -> Result<CompressionMetadata> {
    let reference = read_reference(genome_fasta)?;
    let index = build_index(&reference)?;
    compress_with_index(&index, &reference, input_1, input_2, output, options)
}

pub(crate) fn compress_paired_fastq_files_with_index(
    genome_fasta: impl AsRef<Path>,
    index_prefix: impl AsRef<Path>,
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    let reference = read_reference(genome_fasta)?;
    let index =
        Index::load(index_prefix, false).map_err(|error| Error::InvalidInput(error.to_string()))?;
    compress_with_index(&index, &reference, input_1, input_2, output, options)
}

fn compress_with_index(
    index: &Index,
    reference: &[u8],
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &crate::AlignmentCompressionOptions,
) -> Result<CompressionMetadata> {
    let mapping_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.threads.max(1))
        .build()
        .map_err(|error| Error::InvalidInput(format!("failed to create mapping pool: {error}")))?;
    let progress = ProgressBar::new_spinner();
    progress.enable_steady_tick(std::time::Duration::from_millis(250));
    progress.set_message("aligned 0 reads");
    let temp_dir = tempfile::tempdir()?;
    let spool_path = temp_dir.path().join("aligned-pairs.lz4");
    let mut spool_builder = ExternalChunkBuilder::new(File::create(&spool_path)?, 1)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let mut mask = ReferenceMask::new(reference.len());
    let mut id_detector = IdModeDetector::new(options.preserve_ids);
    let input_1 = input_1.as_ref().to_path_buf();
    let input_2 = input_2.as_ref().to_path_buf();
    let mapper_options = options.clone();
    let (parse_tx, parse_rx) = crossbeam_channel::bounded::<Result<Vec<AlignedPair>>>(4);
    let (mapped_tx, mapped_rx) = crossbeam_channel::bounded::<Result<Vec<AlignedPair>>>(4);
    let result = std::thread::scope(|scope| -> Result<()> {
        let producer_tx = parse_tx.clone();
        let producer = scope.spawn(move || -> Result<()> {
            let mut first = FastxReader::new(
                BufReader::new(gzip_reader(File::open(input_1)?)),
                Format::Fastq,
            );
            let mut second = FastxReader::new(
                BufReader::new(gzip_reader(File::open(input_2)?)),
                Format::Fastq,
            );
            let mut pairs = Vec::with_capacity(PARSE_BATCH_PAIRS);
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
                            second,
                            first_hit: None,
                            second_hit: None,
                        });
                        if pairs.len() == PARSE_BATCH_PAIRS {
                            producer_tx
                                .send(Ok(std::mem::take(&mut pairs)))
                                .map_err(|_| Error::InvalidInput("mapping worker exited".into()))?;
                            pairs = Vec::with_capacity(PARSE_BATCH_PAIRS);
                        }
                    }
                    (None, Some(_)) | (Some(_), None) => {
                        return Err(Error::InvalidInput(
                            "paired streams have different record counts".into(),
                        ));
                    }
                    (Some(Err(error)), _) | (_, Some(Err(error))) => return Err(error),
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
                    map_pairs(&index, &mut pairs, &mapping_pool, &mapper_progress)?;
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
    progress.set_message("pruning reference".to_owned());
    let (compact_reference, kept_runs) = compact_reference(reference, &mask);
    let id_mode = id_detector.finish();
    let reference_block = encode_reference(&compact_reference, options)?;
    let mut archive = File::create(output)?;
    let header_size = 9 + 4 + 2 + 8 * 4;
    archive.write_all(MAGIC)?;
    write_u32(&mut archive, VERSION)?;
    let mut flags = 0u16;
    if options.preserve_ids {
        flags |= 1;
    }
    if options.preserve_quality {
        flags |= 2;
    }
    if options.include_chunk_index {
        flags |= 4;
    }
    flags |= (id_mode as u16) << 3;
    write_u16(&mut archive, flags)?;
    write_u64(&mut archive, reference_block.len() as u64)?;
    let chunk_count_pos = archive.stream_position()?;
    write_u64(&mut archive, 0)?;
    let index_offset_pos = archive.stream_position()?;
    write_u64(&mut archive, 0)?;
    let index_length_pos = archive.stream_position()?;
    write_u64(&mut archive, 0)?;
    debug_assert_eq!(archive.stream_position()?, header_size as u64);
    archive.write_all(&reference_block)?;

    let mut chunks = Vec::new();
    let mut replay_pairs = Vec::with_capacity(options.chunk_pairs.max(1));
    let mut replay_pair_count = 0u64;
    let spool = ExternalChunk::<SpoolPair>::open(File::open(&spool_path)?)
        .map_err(|error| Error::InvalidFormat(error.to_string()))?;
    for item in spool {
        let item = item.map_err(|error| Error::InvalidFormat(error.to_string()))?;
        replay_pairs.push(item.into_aligned_pair(&kept_runs)?);
        if replay_pairs.len() == options.chunk_pairs.max(1) {
            chunks.push(write_chunk(
                &mut archive,
                &compact_reference,
                replay_pair_count,
                &replay_pairs,
                options,
                id_mode,
            )?);
            replay_pair_count += replay_pairs.len() as u64;
            replay_pairs.clear();
        }
    }
    if !replay_pairs.is_empty() {
        chunks.push(write_chunk(
            &mut archive,
            &compact_reference,
            replay_pair_count,
            &replay_pairs,
            options,
            id_mode,
        )?);
        replay_pair_count += replay_pairs.len() as u64;
    }

    let (index_offset, index_length) = if options.include_chunk_index {
        let offset = archive.stream_position()?;
        archive.write_all(b"AL2IDX\0\0")?;
        write_u64(&mut archive, chunks.len() as u64)?;
        for entry in &chunks {
            write_u64(&mut archive, entry.offset)?;
            write_u64(&mut archive, entry.first_pair)?;
            write_u32(&mut archive, entry.pairs)?;
            write_u64(&mut archive, entry.compressed_size)?;
        }
        (offset, archive.stream_position()? - offset)
    } else {
        (0, 0)
    };
    let end = archive.stream_position()?;
    archive.seek(SeekFrom::Start(chunk_count_pos))?;
    write_u64(&mut archive, chunks.len() as u64)?;
    archive.seek(SeekFrom::Start(index_offset_pos))?;
    write_u64(&mut archive, index_offset)?;
    archive.seek(SeekFrom::Start(index_length_pos))?;
    write_u64(&mut archive, index_length)?;
    archive.seek(SeekFrom::Start(end))?;
    archive.flush()?;
    progress.finish_with_message(format!("aligned {} reads", progress.position()));

    Ok(CompressionMetadata {
        streams: 2,
        records: replay_pair_count,
        preserved_ids: options.preserve_ids,
        preserved_quality: options.preserve_quality,
        preserved_order: true,
    })
}

pub(crate) fn decompress_paired<R: Read>(
    input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let mut streams = vec![Vec::new(), Vec::new()];
    decode_paired(input, options, |_, first, second| {
        streams[0].push(first);
        streams[1].push(second);
        Ok(())
    })?;
    Ok(ArchiveRecords { streams })
}

pub(crate) fn decompress_paired_into<R: Read, S: PairedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<crate::DecompressionMetadata> {
    let mut pairs = 0;
    decode_paired(input, options, |original_id, first, second| {
        sink.write_pair(original_id, first, second)?;
        pairs += 1;
        Ok(())
    })?;
    sink.finish()?;
    Ok(crate::DecompressionMetadata {
        streams: 2,
        pairs,
        output_order: crate::OutputOrder::Original,
    })
}

fn decode_paired<R: Read, F: FnMut(u64, ReadRecord, ReadRecord) -> Result<()>>(
    mut input: R,
    options: &DecompressionOptions,
    mut emit: F,
) -> Result<()> {
    let (flags, id_mode, chunk_count, index_offset, reference, mut position) =
        read_header(&mut input)?;
    let preserve_quality = flags & 2 != 0;
    let has_index = flags & 4 != 0;
    let mut first_pair = 0u64;
    for _ in 0..chunk_count {
        if has_index && position == index_offset {
            return Err(Error::InvalidFormat(
                "chunk count reaches index early".into(),
            ));
        }
        let pair_count = read_u32(&mut input)? as usize;
        let compressed_len = usize_from(read_u64(&mut input)?, "chunk compressed size")?;
        let raw_len = usize_from(read_u64(&mut input)?, "chunk uncompressed size")?;
        let ids_len = usize_from(read_u64(&mut input)?, "chunk ID block size")?;
        let qualities_len = usize_from(read_u64(&mut input)?, "chunk quality block size")?;
        let records_len = usize_from(read_u64(&mut input)?, "chunk record block size")?;
        let mut compressed = vec![0; compressed_len];
        input.read_exact(&mut compressed)?;
        position = position
            .checked_add(CHUNK_HEADER_SIZE + compressed_len as u64)
            .ok_or_else(|| Error::InvalidFormat("archive position overflow".into()))?;
        let raw = zstd::stream::decode_all(Cursor::new(compressed))?;
        if raw.len() != raw_len || ids_len + qualities_len + records_len != raw.len() {
            return Err(Error::InvalidFormat(
                "chunk block sizes do not match".into(),
            ));
        }
        let ids_end = ids_len;
        let qualities_end = ids_end + qualities_len;
        let ids = decode_ids(&raw[..ids_end], pair_count, id_mode)?;
        let qualities =
            decode_qualities(&raw[ids_end..qualities_end], pair_count, preserve_quality)?;
        let records = &raw[qualities_end..];
        let decoded = decode_records(records, pair_count, first_pair, &reference)?;
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
                decoded[index].1.clone(),
                qualities[0][index].clone(),
            );
            let second = ReadRecord::new(
                ids[1][index].clone(),
                decoded[index].2.clone(),
                qualities[1][index].clone(),
            );
            emit(original_id, first, second)?;
        }
        first_pair = first_pair
            .checked_add(pair_count as u64)
            .ok_or_else(|| Error::InvalidFormat("pair count overflow".into()))?;
    }
    if has_index && position != index_offset {
        return Err(Error::InvalidFormat("chunk index offset mismatch".into()));
    }
    Ok(())
}

fn map_pairs(
    index: &Index,
    pairs: &mut [AlignedPair],
    mapping_pool: &rayon::ThreadPool,
    progress: &ProgressBar,
) -> Result<()> {
    let first_hits = map_mate(index, pairs, 0, mapping_pool, progress)?;
    let second_hits = map_mate(index, pairs, 1, mapping_pool, progress)?;
    for ((pair, first_hit), second_hit) in pairs.iter_mut().zip(first_hits).zip(second_hits) {
        pair.first_hit = first_hit;
        pair.second_hit = second_hit;
    }
    Ok(())
}

fn spool_mapped_pairs(
    pairs: &[AlignedPair],
    id_detector: &mut IdModeDetector,
    mask: &mut ReferenceMask,
    spool_builder: &mut ExternalChunkBuilder<SpoolPair>,
) -> Result<()> {
    for pair in pairs {
        id_detector.observe(&pair.first.id, &pair.second.id);
        if let Some(hit) = &pair.first_hit {
            mask.mark_hit(hit)?;
        }
        if let Some(hit) = &pair.second_hit {
            mask.mark_hit(hit)?;
        }
        let spool_pair =
            SpoolPair::from_records(pair, pair.first_hit.clone(), pair.second_hit.clone())?;
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
) -> Result<Vec<Option<GenomeAlignment>>> {
    mapping_pool.install(|| {
        use rayon::prelude::*;
        pairs
            .par_chunks(BATCH_SIZE)
            .map(|batch| {
                let sequences = batch.iter().map(|pair| {
                    if mate == 0 {
                        &pair.first.sequence
                    } else {
                        &pair.second.sequence
                    }
                });
                let result = align_batch(index, sequences);
                progress.inc(batch.len() as u64);
                progress.set_message(format!("aligned {} reads", progress.position()));
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
            match op.op {
                CigarOpKind::Match
                | CigarOpKind::Equal
                | CigarOpKind::Diff
                | CigarOpKind::Del
                | CigarOpKind::Skip => {
                    let end = position.checked_add(length).ok_or_else(|| {
                        Error::InvalidInput("reference mask position overflow".into())
                    })?;
                    self.mark_range(position, end)?;
                    position = end;
                }
                CigarOpKind::Ins
                | CigarOpKind::SoftClip
                | CigarOpKind::HardClip
                | CigarOpKind::Pad => {}
                CigarOpKind::Unknown(_) => {
                    return Err(Error::InvalidInput("unsupported CIGAR operation".into()))
                }
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

fn compact_reference(reference: &[u8], mask: &ReferenceMask) -> (Vec<u8>, Vec<KeptRun>) {
    let mut compact = Vec::new();
    let mut runs: Vec<KeptRun> = Vec::new();
    for (position, &base) in reference.iter().enumerate() {
        if !mask.contains(position) {
            continue;
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
    (compact, runs)
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
    fn from_records(
        pair: &AlignedPair,
        first_hit: Option<GenomeAlignment>,
        second_hit: Option<GenomeAlignment>,
    ) -> Result<Self> {
        Ok(Self {
            first: SpoolRead::from_record(&pair.first, first_hit)?,
            second: SpoolRead::from_record(&pair.second, second_hit)?,
        })
    }

    fn into_aligned_pair(self, runs: &[KeptRun]) -> Result<AlignedPair> {
        Ok(AlignedPair {
            first: self.first.record()?,
            second: self.second.record()?,
            first_hit: self
                .first
                .alignment
                .map(|hit| hit.into_genome(runs))
                .transpose()?,
            second_hit: self
                .second
                .alignment
                .map(|hit| hit.into_genome(runs))
                .transpose()?,
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
                        op: cigar_code(op.op)?,
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
                    op: cigar_kind(op.op)?,
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

fn write_chunk(
    archive: &mut File,
    reference: &[u8],
    first_pair: u64,
    pairs: &[AlignedPair],
    options: &crate::AlignmentCompressionOptions,
    id_mode: IdMode,
) -> Result<ChunkIndexEntry> {
    let mut ids = Vec::new();
    write_id_stream(&mut ids, pairs, id_mode)?;
    let mut qualities = Vec::new();
    write_quality_stream(&mut qualities, pairs, options)?;
    let records = write_record_stream(reference, pairs)?;
    let mut raw = Vec::with_capacity(ids.len() + qualities.len() + records.len());
    raw.extend(&ids);
    raw.extend(&qualities);
    raw.extend(&records);
    let zstd_threads = u32::try_from(options.zstd_threads.max(1))
        .map_err(|_| Error::InvalidInput("too many zstd threads".into()))?;
    let mut encoder = zstd::stream::Encoder::new(Vec::new(), options.zstd_level)?;
    encoder.multithread(zstd_threads)?;
    encoder.write_all(&raw)?;
    let compressed = encoder.finish()?;
    let offset = archive.stream_position()?;
    write_u32(archive, pairs.len() as u32)?;
    write_u64(archive, compressed.len() as u64)?;
    write_u64(archive, raw.len() as u64)?;
    write_u64(archive, ids.len() as u64)?;
    write_u64(archive, qualities.len() as u64)?;
    write_u64(archive, records.len() as u64)?;
    archive.write_all(&compressed)?;
    Ok(ChunkIndexEntry {
        offset,
        first_pair,
        pairs: pairs.len() as u32,
        compressed_size: compressed.len() as u64,
    })
}

fn write_id_stream(output: &mut Vec<u8>, pairs: &[AlignedPair], mode: IdMode) -> Result<()> {
    for mate in 0..2 {
        let mut block = Vec::new();
        if mode != IdMode::None && !(mode != IdMode::Independent && mate == 1) {
            let ids = pairs.iter().map(|pair| {
                let id = if mate == 0 {
                    &pair.first.id
                } else {
                    &pair.second.id
                };
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
) -> Result<()> {
    for mate in 0..2 {
        let values = pairs
            .iter()
            .map(|pair| {
                if mate == 0 {
                    pair.first.quality.as_deref().unwrap_or(&[])
                } else {
                    pair.second.quality.as_deref().unwrap_or(&[])
                }
            })
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

fn write_record_stream(reference: &[u8], pairs: &[AlignedPair]) -> Result<Vec<u8>> {
    let mut output = Vec::new();
    let mut starts = [0i64; 2];
    for mate in 0..2 {
        let bitmap_len = pairs.len().div_ceil(8);
        let mut mapped_bitmap = vec![0; bitmap_len];
        let mut reverse_bitmap = vec![0; bitmap_len];
        for (index, pair) in pairs.iter().enumerate() {
            let hit = if mate == 0 {
                pair.first_hit.as_ref()
            } else {
                pair.second_hit.as_ref()
            };
            if let Some(hit) = hit {
                mapped_bitmap[index / 8] |= 1 << (index % 8);
                if hit.is_reverse {
                    reverse_bitmap[index / 8] |= 1 << (index % 8);
                }
            }
        }
        output.extend(mapped_bitmap);
        output.extend(reverse_bitmap);
        for pair in pairs {
            let record = if mate == 0 { &pair.first } else { &pair.second };
            let hit = if mate == 0 {
                pair.first_hit.as_ref()
            } else {
                pair.second_hit.as_ref()
            };
            encode_record(&mut output, record, hit, reference, &mut starts[mate])?;
        }
    }
    Ok(output)
}

fn align_batch<'a, I>(index: &Index, sequences: I) -> Result<Vec<Option<GenomeAlignment>>>
where
    I: IntoIterator<Item = &'a Vec<u8>>,
{
    let sequences = sequences.into_iter().collect::<Vec<_>>();
    let names = (0..sequences.len())
        .map(|index| format!("read-{index}"))
        .collect::<Vec<_>>();
    let queries = names
        .iter()
        .zip(&sequences)
        .map(|(name, sequence)| Query {
            name,
            sequence,
            methylation_type: 0,
        })
        .collect::<Vec<_>>();
    let mut options = Options::default();
    options.set_paired_end(false);
    minibwa::map_batch(index, &options, &queries, None)
        .map(|hits| {
            hits.into_iter()
                .map(|hits| {
                    hits.into_iter().next().map(|hit| GenomeAlignment {
                        start: hit.target_start,
                        query_start: hit.query_start,
                        query_end: hit.query_end,
                        cigar: hit.cigar,
                        is_reverse: hit.is_reverse,
                    })
                })
                .collect()
        })
        .map_err(|error| Error::InvalidInput(error.to_string()))
}

fn encode_record(
    output: &mut Vec<u8>,
    record: &ReadRecord,
    hit: Option<&GenomeAlignment>,
    reference: &[u8],
    previous_start: &mut i64,
) -> Result<()> {
    if let Some(hit) = hit {
        let delta = hit
            .start
            .checked_sub(*previous_start)
            .ok_or_else(|| Error::InvalidInput("alignment offset overflow".into()))?;
        *previous_start = hit.start;
        let cigar = effective_cigar(hit, record.sequence.len())?;
        write_varint(output, zigzag_encode(delta))?;
        write_varint(output, cigar.len() as u64)?;
        for chunk in cigar.chunks(2) {
            let first = cigar_code(chunk[0].op)? << 4;
            let second = if chunk.len() == 2 {
                cigar_code(chunk[1].op)?
            } else {
                0
            };
            write_u8(output, first | second)?;
        }
        for op in &cigar {
            write_varint(output, op.len as u64)?;
        }
        encode_noise(output, &record.sequence, hit, &cigar, reference)?;
    } else {
        write_packed(output, &record.sequence)?;
    }
    Ok(())
}

fn decode_records(
    input: &[u8],
    count: usize,
    first_pair: u64,
    reference: &[u8],
) -> Result<Vec<(u64, Vec<u8>, Vec<u8>)>> {
    let mut cursor = Cursor::new(input.to_vec());
    let mut out = Vec::with_capacity(count);
    let mut starts = [0i64; 2];
    let mut all = Vec::new();
    for mate in 0..2 {
        let mut stream = Vec::with_capacity(count);
        let bitmap_len = count.div_ceil(8);
        let mut mapped_bitmap = vec![0; bitmap_len];
        let mut reverse_bitmap = vec![0; bitmap_len];
        cursor.read_exact(&mut mapped_bitmap)?;
        cursor.read_exact(&mut reverse_bitmap)?;
        for ordinal in 0..count {
            let mapped = mapped_bitmap[ordinal / 8] & (1 << (ordinal % 8)) != 0;
            let reverse = reverse_bitmap[ordinal / 8] & (1 << (ordinal % 8)) != 0;
            let sequence = if mapped {
                let delta = zigzag_decode(read_varint(&mut cursor)?)?;
                starts[mate] = starts[mate]
                    .checked_add(delta)
                    .ok_or_else(|| Error::InvalidFormat("alignment offset overflow".into()))?;
                decode_mapped_body(&mut cursor, reference, starts[mate], reverse)?
            } else {
                read_packed(&mut cursor)?.decode()
            };
            stream.push(sequence);
        }
        all.push(stream);
    }
    if cursor.position() != cursor.get_ref().len() as u64 {
        return Err(Error::InvalidFormat("trailing record bytes".into()));
    }
    for index in 0..count {
        out.push((
            first_pair + index as u64,
            all[0][index].clone(),
            all[1][index].clone(),
        ));
    }
    Ok(out)
}

fn decode_ids(input: &[u8], count: usize, mode: IdMode) -> Result<[Vec<Vec<u8>>; 2]> {
    let mut cursor = Cursor::new(input.to_vec());
    let mut blocks = [Vec::new(), Vec::new()];
    for mate in 0..2 {
        let length = usize_from(read_u64(&mut cursor)?, "ID block")?;
        blocks[mate] = vec![0; length];
        cursor.read_exact(&mut blocks[mate])?;
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
    Ok([first, second])
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
) -> Result<[Vec<Option<Vec<u8>>>; 2]> {
    let mut cursor = Cursor::new(input.to_vec());
    let mut result = [Vec::with_capacity(count), Vec::with_capacity(count)];
    for mate in 0..2 {
        let length = usize_from(read_u64(&mut cursor)?, "quality block")?;
        let mut block = vec![0; length];
        cursor.read_exact(&mut block)?;
        if preserve {
            result[mate] = decode_stream_bytes(&block, count)?
                .into_iter()
                .map(Some)
                .collect();
        } else if length != 0 {
            return Err(Error::InvalidFormat("unexpected quality bytes".into()));
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
            op: cigar_kind(op)?,
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
            op: CigarOpKind::SoftClip,
            len: prefix as u32,
        });
    }
    cigar.extend_from_slice(&hit.cigar);
    if suffix != 0 {
        cigar.push(CigarOp {
            op: CigarOpKind::SoftClip,
            len: suffix as u32,
        });
    }
    Ok(cigar)
}

fn encode_noise(
    output: &mut Vec<u8>,
    sequence: &[u8],
    hit: &GenomeAlignment,
    cigar: &[CigarOp],
    reference: &[u8],
) -> Result<()> {
    let mut oriented = sequence.to_vec();
    if hit.is_reverse {
        reverse_complement(&mut oriented);
    }
    let mut ref_pos = nonnegative_usize(hit.start, "alignment start")?;
    let mut query_pos = 0usize;
    for op in cigar {
        let len = op.len as usize;
        match op.op {
            CigarOpKind::Match | CigarOpKind::Equal | CigarOpKind::Diff => {
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
                    if reference[ref_pos + offset] != oriented[query_pos + offset] {
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
            CigarOpKind::Ins | CigarOpKind::SoftClip => {
                let end = query_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("query position overflow".into()))?;
                if end > oriented.len() {
                    return Err(Error::InvalidFormat("CIGAR exceeds query bounds".into()));
                }
                output.extend(&oriented[query_pos..end]);
                query_pos = end;
            }
            CigarOpKind::Del | CigarOpKind::Skip => {
                ref_pos = ref_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidInput("reference position overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "CIGAR exceeds reference bounds".into(),
                    ));
                }
            }
            CigarOpKind::HardClip | CigarOpKind::Pad => {}
            CigarOpKind::Unknown(_) => {
                return Err(Error::InvalidInput("unsupported CIGAR operation".into()))
            }
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
        match op.op {
            CigarOpKind::Match | CigarOpKind::Equal | CigarOpKind::Diff => {
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
            CigarOpKind::Ins | CigarOpKind::SoftClip => {
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
            CigarOpKind::Del | CigarOpKind::Skip => {
                ref_pos = ref_pos
                    .checked_add(len)
                    .ok_or_else(|| Error::InvalidFormat("reference position overflow".into()))?;
                if ref_pos > reference.len() {
                    return Err(Error::InvalidFormat(
                        "CIGAR exceeds reference bounds".into(),
                    ));
                }
            }
            CigarOpKind::HardClip | CigarOpKind::Pad => {}
            CigarOpKind::Unknown(_) => {
                return Err(Error::InvalidFormat("unsupported CIGAR operation".into()))
            }
        }
    }
    Ok(sequence)
}

fn read_header<R: Read>(input: &mut R) -> Result<(u16, IdMode, u64, u64, Vec<u8>, u64)> {
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
    let reference_len = usize_from(read_u64(input)?, "reference block")?;
    let chunks = read_u64(input)?;
    let index_offset = read_u64(input)?;
    let _index_len = read_u64(input)?;
    let mut reference_block = vec![0; reference_len];
    input.read_exact(&mut reference_block)?;
    let position = 9 + 4 + 2 + 8 * 4 + reference_len as u64;
    Ok((
        flags,
        id_mode,
        chunks,
        index_offset,
        decode_reference(&reference_block)?,
        position,
    ))
}

fn read_reference(path: impl AsRef<Path>) -> Result<Vec<u8>> {
    let path = path.as_ref();
    let input = File::open(path)?;
    let reader: Box<dyn Read> = if path.extension().is_some_and(|extension| extension == "gz") {
        Box::new(gzip_reader(input))
    } else {
        Box::new(input)
    };
    let mut reference = Vec::new();
    for record in FastxReader::new(BufReader::new(reader), Format::Fasta) {
        reference.extend(record?.sequence);
    }
    if reference.is_empty() {
        return Err(Error::InvalidInput("genome FASTA is empty".into()));
    }
    Ok(reference)
}

fn build_index(reference: &[u8]) -> Result<Index> {
    let mut fasta = tempfile::NamedTempFile::new()?;
    fasta.write_all(b">genome\n")?;
    fasta.write_all(reference)?;
    fasta.write_all(b"\n")?;
    fasta.as_file_mut().sync_all()?;
    let prefix = tempfile::tempdir()?;
    Index::build(fasta.path(), prefix.path().join("index"), 4, false)
        .map_err(|error| Error::InvalidInput(error.to_string()))
}

fn encode_reference(
    reference: &[u8],
    options: &crate::AlignmentCompressionOptions,
) -> Result<Vec<u8>> {
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

    let zstd_threads = u32::try_from(options.zstd_threads.max(1))
        .map_err(|_| Error::InvalidInput("too many zstd threads".into()))?;
    let mut encoder = zstd::stream::Encoder::new(Vec::new(), options.zstd_level)?;
    encoder.multithread(zstd_threads)?;
    encoder.write_all(&raw)?;
    let compressed = encoder.finish()?;
    let mut block = Vec::new();
    write_u64(&mut block, raw.len() as u64)?;
    block.extend(compressed);
    Ok(block)
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
    write_u8(&mut output, u8::from(packed.is_packed()))?;
    write_u64(&mut output, packed.encoded_size() as u64)?;
    output.extend(packed.encoded_bytes());
    Ok(output)
}
fn write_packed(output: &mut Vec<u8>, sequence: &[u8]) -> Result<()> {
    output.extend(encode_packed(sequence)?);
    Ok(())
}
fn read_packed(input: &mut Cursor<Vec<u8>>) -> Result<PackedSequence> {
    let length = read_u32(input)?;
    let packed = read_u8(input)? != 0;
    let body_len = usize_from(read_u64(input)?, "packed sequence body")?;
    let mut body = vec![0; body_len];
    input.read_exact(&mut body)?;
    PackedSequence::from_encoded(length, packed, body)
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
fn cigar_code(op: CigarOpKind) -> Result<u8> {
    Ok(match op {
        CigarOpKind::Match => 0,
        CigarOpKind::Ins => 1,
        CigarOpKind::Del => 2,
        CigarOpKind::Skip => 3,
        CigarOpKind::SoftClip => 4,
        CigarOpKind::HardClip => 5,
        CigarOpKind::Pad => 6,
        CigarOpKind::Equal => 7,
        CigarOpKind::Diff => 8,
        CigarOpKind::Unknown(_) => {
            return Err(Error::InvalidInput("unsupported CIGAR operation".into()))
        }
    })
}
fn cigar_kind(code: u8) -> Result<CigarOpKind> {
    Ok(match code {
        0 => CigarOpKind::Match,
        1 => CigarOpKind::Ins,
        2 => CigarOpKind::Del,
        3 => CigarOpKind::Skip,
        4 => CigarOpKind::SoftClip,
        5 => CigarOpKind::HardClip,
        6 => CigarOpKind::Pad,
        7 => CigarOpKind::Equal,
        8 => CigarOpKind::Diff,
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_block_round_trips_ambiguous_bases() {
        let reference = b"ACGTNNNNacgtRYACGT";
        let encoded =
            encode_reference(reference, &crate::AlignmentCompressionOptions::default()).unwrap();
        assert_eq!(decode_reference(&encoded).unwrap(), reference);
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
                second: ReadRecord::new(second_id, b"A".to_vec(), None),
                first_hit: None,
                second_hit: None,
            };
            let mut encoded = Vec::new();
            write_id_stream(&mut encoded, &[pair], mode).unwrap();
            let decoded = decode_ids(&encoded, 1, mode).unwrap();
            assert_eq!(decoded[0][0], expected_first);
            assert_eq!(decoded[1][0], expected_second);
        }
    }
}
