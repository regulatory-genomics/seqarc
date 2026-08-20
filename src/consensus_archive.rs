use crate::external_bytes::{copy_to as copy_external_bytes, Lz4ByteWriter};
use crate::packed_sequence::PackedSequence;
#[cfg(test)]
use crate::record_store::PairStore;
#[cfg(test)]
use crate::record_store::StoreStats;
use crate::reorder::{ConsensusCluster, ConsensusRead};
use crate::sequence_spool::{SequenceSpoolReader, SequenceSpoolWriter, SortableRead};
use crate::{
    ArchiveRecords, CompressionMetadata, CompressionOptions, CompressionProfile,
    DecompressionOptions, Error, ReadRecord, Result,
};
use bed_utils::extsort::{ExternalChunk, ExternalChunkBuilder};
use crossbeam_channel::bounded;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

const MAGIC: &[u8; 9] = b"SPRNG-RS4";
const VERSION: u32 = 5;
const LONG_DISTANCE_WINDOW_LOG: u32 = 27;
const ORDER_ENTRY_SIZE: usize = 9;
const MAX_OVERLAP_MISMATCHES: usize = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairedIdPattern {
    Independent = 0,
    Identical = 1,
    SlashMate = 2,
}

impl PairedIdPattern {
    fn detect(first: &[u8], second: &[u8]) -> Self {
        if first == second {
            Self::Identical
        } else if first.len() == second.len()
            && first.ends_with(b"/1")
            && second.ends_with(b"/2")
            && first[..first.len() - 2] == second[..second.len() - 2]
        {
            Self::SlashMate
        } else {
            Self::Independent
        }
    }

    fn derive_second(self, first: &[u8]) -> Result<Vec<u8>> {
        match self {
            Self::Independent => Err(Error::InvalidFormat(
                "cannot derive an independent mate ID".into(),
            )),
            Self::Identical => Ok(first.to_vec()),
            Self::SlashMate if first.ends_with(b"/1") => {
                let mut second = first.to_vec();
                *second.last_mut().unwrap() = b'2';
                Ok(second)
            }
            Self::SlashMate => Err(Error::InvalidFormat(
                "invalid mate-1 ID for /1-/2 pattern".into(),
            )),
        }
    }
}

impl TryFrom<u8> for PairedIdPattern {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self> {
        match value {
            0 => Ok(Self::Independent),
            1 => Ok(Self::Identical),
            2 => Ok(Self::SlashMate),
            _ => Err(Error::InvalidFormat("unknown paired ID pattern".into())),
        }
    }
}

pub fn read_reordered<R: Read>(
    mut input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let mut magic = [0; 9];
    input.read_exact(&mut magic)?;
    if magic != *MAGIC {
        return Err(Error::InvalidFormat("bad RS4 magic".into()));
    }
    let version = read_u32(&mut input)?;
    if version != VERSION {
        return Err(Error::InvalidFormat("unsupported RS4 version".into()));
    }
    let streams = read_u16(&mut input)? as usize;
    let flags = read_u16(&mut input)?;
    let raw_len = usize_from(read_u64(&mut input)?, "payload")?;
    let encoded_len = usize_from(read_u64(&mut input)?, "payload")?;
    let expected = read_u32(&mut input)?;
    let mut encoded = input.take(encoded_len as u64);
    let mut decoder = zstd::stream::read::Decoder::new(&mut encoded)?;
    read_payload_collect(
        &mut decoder,
        raw_len as u64,
        expected,
        streams,
        flags,
        options,
    )
}

fn read_payload_collect<R: Read>(
    payload: &mut R,
    raw_len: u64,
    expected: u32,
    streams: usize,
    flags: u16,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let mut decoded = (0..streams)
        .map(|_| Vec::<(u64, ReadRecord)>::new())
        .collect::<Vec<_>>();
    let (records, id_pattern) = decode_payload_records(
        payload,
        raw_len,
        expected,
        streams,
        flags,
        |stream, original_id, record| {
            decoded[stream].push((original_id, record));
            Ok(())
        },
    )?;
    if flags & 1 != 0 && id_pattern != PairedIdPattern::Independent {
        let first_ids = decoded[0]
            .iter()
            .map(|(id, record)| (*id, record.id.clone()))
            .collect::<HashMap<_, _>>();
        for (original_id, record) in &mut decoded[1] {
            let first_id = first_ids
                .get(original_id)
                .ok_or_else(|| Error::InvalidFormat("paired ID original ID mismatch".into()))?;
            record.id = id_pattern.derive_second(first_id)?;
        }
    }
    if options.output_order == crate::OutputOrder::Original {
        for stream in &mut decoded {
            stream.sort_unstable_by_key(|(original_id, _)| *original_id);
        }
    }
    let streams: Vec<Vec<ReadRecord>> = decoded
        .into_iter()
        .map(|stream| stream.into_iter().map(|(_, record)| record).collect())
        .collect();
    if streams
        .iter()
        .any(|stream: &Vec<ReadRecord>| stream.len() as u64 != records)
    {
        return Err(Error::InvalidFormat("RS4 record count mismatch".into()));
    }
    finish_decoding(streams, options)
}

pub(crate) fn read_paired_into<R: Read, S: crate::PairedRecordSink>(
    mut input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<crate::DecompressionMetadata> {
    let mut magic = [0; 9];
    input.read_exact(&mut magic)?;
    if magic != *MAGIC {
        return Err(Error::InvalidFormat("bad RS4 magic".into()));
    }
    if read_u32(&mut input)? != VERSION {
        return Err(Error::InvalidFormat("unsupported RS4 version".into()));
    }
    let streams = read_u16(&mut input)? as usize;
    let flags = read_u16(&mut input)?;
    let raw_len = read_u64(&mut input)?;
    let encoded_len = read_u64(&mut input)?;
    let expected = read_u32(&mut input)?;
    if streams != 2 {
        return Err(Error::InvalidFormat(
            "paired decompression requires two streams".into(),
        ));
    }
    let mut encoded = input.take(encoded_len);
    let mut decoder = zstd::stream::read::Decoder::new(&mut encoded)?;
    let temp_root = options.temp_dir.clone().unwrap_or_else(std::env::temp_dir);
    std::fs::create_dir_all(&temp_root)?;
    let workspace = unique_decode_workspace(&temp_root)?;
    let result = (|| -> Result<crate::DecompressionMetadata> {
        let mut run_writers = [
            RunWriter::new(workspace.join("mate-0"))?,
            RunWriter::new(workspace.join("mate-1"))?,
        ];
        let (pairs, id_pattern) = decode_payload_records(
            &mut decoder,
            raw_len,
            expected,
            streams,
            flags,
            |stream, original_id, record| run_writers[stream].push(original_id, record),
        )?;
        let [first_writer, second_writer] = run_writers;
        let [first_runs, second_runs] = [first_writer.finish()?, second_writer.finish()?];
        merge_paired_runs(
            first_runs,
            second_runs,
            &workspace,
            pairs,
            flags,
            id_pattern,
            options,
            sink,
        )?;
        sink.finish()?;
        Ok(crate::DecompressionMetadata {
            streams: 2,
            pairs,
            output_order: options.output_order,
        })
    })();
    let _ = std::fs::remove_dir_all(&workspace);
    result
}

fn decode_payload_records<R: Read, F: FnMut(usize, u64, ReadRecord) -> Result<()>>(
    payload: &mut R,
    raw_len: u64,
    expected: u32,
    streams: usize,
    flags: u16,
    mut visit: F,
) -> Result<(u64, PairedIdPattern)> {
    let mut payload = ChecksummedReader::new(payload, raw_len);
    let records = payload.u64()?;
    let payload_streams = payload.u16()? as usize;
    if streams == 0 || streams != payload_streams {
        return Err(Error::InvalidFormat("RS4 stream count mismatch".into()));
    }
    let id_pattern = PairedIdPattern::try_from(payload.u8()?)?;
    let has_ids = flags & 1 != 0;
    let has_quality = flags & 2 != 0;
    if id_pattern != PairedIdPattern::Independent && (!has_ids || streams != 2) {
        return Err(Error::InvalidFormat(
            "paired ID pattern requires two ID-bearing streams".into(),
        ));
    }
    for stream_index in 0..streams {
        let multi_groups = payload.u64()?;
        let multi_records = payload.u64()?;
        let singleton_buckets = payload.u64()?;
        let singleton_records = payload.u64()?;
        if multi_records + singleton_records != records {
            return Err(Error::InvalidFormat(
                "RS4 section record count mismatch".into(),
            ));
        }
        let ids_enabled =
            has_ids && (stream_index == 0 || id_pattern == PairedIdPattern::Independent);
        let mut decoded_multi = 0u64;
        for _ in 0..multi_groups {
            decoded_multi += decode_payload_group(
                &mut payload,
                stream_index,
                multi_records,
                0,
                ids_enabled,
                has_quality,
                &mut visit,
            )?;
        }
        let mut decoded_singletons = 0u64;
        for _ in 0..singleton_buckets {
            decoded_singletons += decode_payload_group(
                &mut payload,
                stream_index,
                singleton_records,
                2,
                ids_enabled,
                has_quality,
                &mut visit,
            )?;
        }
        if decoded_multi != multi_records || decoded_singletons != singleton_records {
            return Err(Error::InvalidFormat(
                "RS4 section group count mismatch".into(),
            ));
        }
    }
    payload.finish(expected)?;
    Ok((records, id_pattern))
}

fn decode_payload_group<R: Read, F: FnMut(usize, u64, ReadRecord) -> Result<()>>(
    payload: &mut ChecksummedReader<'_, R>,
    stream_index: usize,
    section_records: u64,
    expected_kind: u8,
    ids_enabled: bool,
    has_quality: bool,
    visit: &mut F,
) -> Result<u64> {
    let group_records = payload.varint()?;
    if group_records == 0 || group_records > section_records {
        return Err(Error::InvalidFormat(
            "invalid RS4 group record count".into(),
        ));
    }
    let kind = payload.u8()?;
    let singleton = expected_kind == 2;
    if (singleton && !(kind == 1 || kind == 2 || kind == 3)) || (!singleton && kind != 0) {
        return Err(Error::InvalidFormat(format!(
            "invalid RS4 section group kind {kind} (expected {expected_kind})"
        )));
    }
    let reference = if kind == 0 {
        read_reference_stream(payload)?
    } else {
        Vec::new()
    };
    let fixed_length = if kind == 2 || kind == 3 {
        Some(
            u32::try_from(payload.varint()?)
                .map_err(|_| Error::InvalidFormat("sequence too large".into()))?,
        )
    } else {
        None
    };
    let mut records_in_group = Vec::with_capacity(usize_from(group_records, "group")?);
    let mut previous_id = 0u64;
    for _ in 0..group_records {
        let original_id = zigzag_decode(payload.varint()?, previous_id)?;
        previous_id = original_id;
        let sequence = if kind == 1 {
            read_packed_sequence_stream(payload)?
        } else if kind == 2 || kind == 3 {
            let length = fixed_length.unwrap();
            let encoded_length = if kind == 2 {
                (length as usize).div_ceil(4)
            } else {
                length as usize
            };
            let mut bytes = vec![0; encoded_length];
            payload.read_exact(&mut bytes)?;
            PackedSequence::from_encoded(length, kind == 2, bytes)?.decode()
        } else {
            let read = read_consensus_read_stream(payload, original_id)?;
            crate::reorder::decode_consensus_read(
                &ConsensusCluster {
                    reference: reference.clone(),
                    reads: Vec::new(),
                },
                &read,
            )?
        };
        records_in_group.push((original_id, sequence));
    }
    let ids = if ids_enabled {
        read_group_metadata(payload, group_records, false)?
    } else {
        vec![Vec::new(); usize_from(group_records, "group")?]
    };
    let qualities = if has_quality {
        read_group_metadata(payload, group_records, true)?
    } else {
        vec![Vec::new(); usize_from(group_records, "group")?]
    };
    for (index, (original_id, sequence)) in records_in_group.into_iter().enumerate() {
        visit(
            stream_index,
            original_id,
            ReadRecord {
                id: ids[index].clone(),
                sequence,
                quality: if has_quality {
                    Some(qualities[index].clone())
                } else {
                    None
                },
            },
        )?;
    }
    Ok(group_records)
}

fn unique_decode_workspace(root: &Path) -> Result<std::path::PathBuf> {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    for _ in 0..1000 {
        let path = root.join(format!(
            "spring-rs-decode-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        match std::fs::create_dir(&path) {
            Ok(()) => return Ok(path),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(Error::Io(error)),
        }
    }
    Err(Error::InvalidInput(
        "unable to allocate decode workspace".into(),
    ))
}

fn finish_decoding(
    mut streams: Vec<Vec<ReadRecord>>,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let range = options.original_range.as_ref();
    if let Some(range) = range {
        for stream in &mut streams {
            *stream = stream
                .iter()
                .enumerate()
                .filter(|(index, _)| range.contains(&(*index as u64)))
                .map(|(_, record)| record.clone())
                .collect();
        }
    }
    Ok(ArchiveRecords { streams })
}

struct RunWriter {
    prefix: std::path::PathBuf,
    next_run: usize,
    records: Vec<(u64, ReadRecord)>,
    paths: Vec<std::path::PathBuf>,
    count: u64,
}

struct RunSet {
    paths: Vec<std::path::PathBuf>,
    records: u64,
}

const RUN_MERGE_FAN_IN: usize = 64;

impl RunWriter {
    const RUN_RECORDS: usize = 65_536;

    fn new(prefix: std::path::PathBuf) -> Result<Self> {
        std::fs::create_dir_all(&prefix)?;
        Ok(Self {
            prefix,
            next_run: 0,
            records: Vec::with_capacity(Self::RUN_RECORDS),
            paths: Vec::new(),
            count: 0,
        })
    }

    fn push(&mut self, original_id: u64, record: ReadRecord) -> Result<()> {
        self.records.push((original_id, record));
        self.count = self
            .count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many decoded records".into()))?;
        if self.records.len() == Self::RUN_RECORDS {
            self.flush_run_with_paths()?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<RunSet> {
        self.flush_run_with_paths()?;
        Ok(RunSet {
            paths: self.paths,
            records: self.count,
        })
    }

    fn flush_run_with_paths(&mut self) -> Result<()> {
        if self.records.is_empty() {
            return Ok(());
        }
        self.records.sort_unstable_by_key(|(id, _)| *id);
        let path = self.prefix.join(format!("run-{:08}.bin", self.next_run));
        self.next_run += 1;
        let mut file = std::io::BufWriter::new(File::create(path.clone())?);
        file.write_all(&(self.records.len() as u64).to_le_bytes())?;
        for (original_id, record) in self.records.drain(..) {
            write_run_record(&mut file, original_id, &record)?;
        }
        file.flush()?;
        self.paths.push(path);
        Ok(())
    }
}

fn write_run_record<W: Write>(writer: &mut W, original_id: u64, record: &ReadRecord) -> Result<()> {
    writer.write_all(&original_id.to_le_bytes())?;
    for field in [&record.id, &record.sequence] {
        writer.write_all(&(field.len() as u64).to_le_bytes())?;
        writer.write_all(field)?;
    }
    match &record.quality {
        Some(quality) => {
            writer.write_all(&[1])?;
            writer.write_all(&(quality.len() as u64).to_le_bytes())?;
            writer.write_all(quality)?;
        }
        None => writer.write_all(&[0])?,
    }
    Ok(())
}

struct RunReader {
    file: BufReader<File>,
    remaining: u64,
}

impl RunReader {
    fn open(path: &Path) -> Result<Self> {
        let mut file = BufReader::new(File::open(path)?);
        let mut count = [0u8; 8];
        file.read_exact(&mut count)?;
        Ok(Self {
            file,
            remaining: u64::from_le_bytes(count),
        })
    }

    fn next(&mut self) -> Result<Option<(u64, ReadRecord)>> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let mut original_id = [0u8; 8];
        self.file.read_exact(&mut original_id)?;
        let id = read_run_field(&mut self.file)?;
        let sequence = read_run_field(&mut self.file)?;
        let mut present = [0u8; 1];
        self.file.read_exact(&mut present)?;
        let quality = match present[0] {
            0 => None,
            1 => Some(read_run_field(&mut self.file)?),
            _ => return Err(Error::InvalidFormat("invalid RS4 run quality flag".into())),
        };
        self.remaining -= 1;
        Ok(Some((
            u64::from_le_bytes(original_id),
            ReadRecord {
                id,
                sequence,
                quality,
            },
        )))
    }
}

fn read_run_field<R: Read>(reader: &mut R) -> Result<Vec<u8>> {
    let mut length = [0u8; 8];
    reader.read_exact(&mut length)?;
    let length = usize_from(u64::from_le_bytes(length), "run field")?;
    let mut field = vec![0; length];
    reader.read_exact(&mut field)?;
    Ok(field)
}

struct HeapEntry {
    original_id: u64,
    run: usize,
    record: ReadRecord,
}

impl PartialEq for HeapEntry {
    fn eq(&self, other: &Self) -> bool {
        (self.original_id, self.run) == (other.original_id, other.run)
    }
}

impl Eq for HeapEntry {}

impl PartialOrd for HeapEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for HeapEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        (other.original_id, other.run).cmp(&(self.original_id, self.run))
    }
}

struct RunMerger {
    readers: Vec<RunReader>,
    heap: std::collections::BinaryHeap<HeapEntry>,
    count: u64,
}

impl RunMerger {
    fn open(paths: &[std::path::PathBuf]) -> Result<Self> {
        let mut readers = Vec::with_capacity(paths.len());
        let mut heap = std::collections::BinaryHeap::new();
        for path in paths {
            readers.push(RunReader::open(path)?);
        }
        for run in 0..readers.len() {
            if let Some((original_id, record)) = readers[run].next()? {
                heap.push(HeapEntry {
                    original_id,
                    run,
                    record,
                });
            }
        }
        Ok(Self {
            readers,
            heap,
            count: 0,
        })
    }

    fn next(&mut self) -> Result<Option<(u64, ReadRecord)>> {
        let Some(entry) = self.heap.pop() else {
            return Ok(None);
        };
        if let Some((original_id, record)) = self.readers[entry.run].next()? {
            self.heap.push(HeapEntry {
                original_id,
                run: entry.run,
                record,
            });
        }
        self.count += 1;
        Ok(Some((entry.original_id, entry.record)))
    }
}

fn merge_run_batch(paths: &[std::path::PathBuf], output: &Path) -> Result<u64> {
    let mut merger = RunMerger::open(paths)?;
    let mut output = std::io::BufWriter::new(File::create(output)?);
    output.write_all(&0u64.to_le_bytes())?;
    let mut count = 0u64;
    while let Some((original_id, record)) = merger.next()? {
        write_run_record(&mut output, original_id, &record)?;
        count = count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many merged run records".into()))?;
    }
    output.flush()?;
    let mut output = output
        .into_inner()
        .map_err(|error| Error::Io(error.into_error()))?;
    output.seek(SeekFrom::Start(0))?;
    output.write_all(&count.to_le_bytes())?;
    output.flush()?;
    Ok(count)
}

fn compact_run_set(set: &mut RunSet, workspace: &Path) -> Result<()> {
    if set.paths.len() <= RUN_MERGE_FAN_IN {
        return Ok(());
    }
    let merge_root = workspace.join("merged-runs");
    std::fs::create_dir_all(&merge_root)?;
    let mut current = std::mem::take(&mut set.paths);
    let mut pass = 0usize;
    while current.len() > RUN_MERGE_FAN_IN {
        let mut next = Vec::with_capacity(current.len().div_ceil(RUN_MERGE_FAN_IN));
        for (batch, input_paths) in current.chunks(RUN_MERGE_FAN_IN).enumerate() {
            let output = merge_root.join(format!("pass-{pass:04}-{batch:08}.bin"));
            merge_run_batch(input_paths, &output)?;
            for path in input_paths {
                std::fs::remove_file(path)?;
            }
            next.push(output);
        }
        current = next;
        pass += 1;
    }
    set.paths = current;
    Ok(())
}

fn merge_paired_runs<S: crate::PairedRecordSink>(
    mut first_set: RunSet,
    mut second_set: RunSet,
    workspace: &Path,
    pairs: u64,
    flags: u16,
    id_pattern: PairedIdPattern,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<()> {
    if first_set.records != pairs || second_set.records != pairs {
        return Err(Error::InvalidFormat(
            "paired stream record count mismatch".into(),
        ));
    }
    compact_run_set(&mut first_set, workspace)?;
    compact_run_set(&mut second_set, workspace)?;
    let mut first = RunMerger::open(&first_set.paths)?;
    let mut second = RunMerger::open(&second_set.paths)?;
    for expected_id in 0..pairs {
        let (first_id, first_record) = first
            .next()?
            .ok_or_else(|| Error::InvalidFormat("truncated first mate runs".into()))?;
        let (second_id, mut second_record) = second
            .next()?
            .ok_or_else(|| Error::InvalidFormat("truncated second mate runs".into()))?;
        if first_id != expected_id || second_id != expected_id {
            return Err(Error::InvalidFormat(
                "paired run original ID mismatch".into(),
            ));
        }
        if flags & 1 != 0 && id_pattern != PairedIdPattern::Independent {
            second_record.id = id_pattern.derive_second(&first_record.id)?;
        }
        if options
            .original_range
            .as_ref()
            .is_none_or(|range| range.contains(&expected_id))
        {
            sink.write_pair(expected_id, first_record, second_record)?;
        }
    }
    if first.next()?.is_some() || second.next()?.is_some() {
        return Err(Error::InvalidFormat("trailing paired run records".into()));
    }
    Ok(())
}

#[cfg(test)]
pub fn write_reordered_store_profiled<W: Write>(
    store: &mut PairStore,
    mut output: W,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<(CompressionMetadata, CompressionProfile)> {
    let temp_dir = temp_dir.as_ref();
    let mut prepared = PreparedRecordsBuilder::create(temp_dir)?;
    for id in 0..store.count() {
        let (first, second) = store.read_pair(id)?;
        prepared.push(&first, &second, options)?;
    }
    let mut result = write_prepared_profiled(prepared.finish()?, &mut output, options, temp_dir)?;
    result.1.store_stats = StoreStats {
        index_reads: store.count().saturating_mul(5),
        blob_reads: store.count().saturating_mul(8),
        bytes_read: store.stats().bytes_read,
    };
    return Ok(result);
    #[allow(unreachable_code)]
    {
        let temp_dir = temp_dir.as_ref();
        std::fs::create_dir_all(temp_dir)?;
        let assignment_start = std::time::Instant::now();
        let assignments = [
            crate::reorder::write_store_assignments(store, &options.reorder, temp_dir, 0)?,
            crate::reorder::write_store_assignments(store, &options.reorder, temp_dir, 1)?,
        ];
        let assignment_time = assignment_start.elapsed();
        let assignment_bytes = assignments.iter().try_fold(0u64, |total, path| {
            Ok::<_, Error>(total + std::fs::metadata(path)?.len())
        })?;

        let payload_start = std::time::Instant::now();
        let (record_columns, id_pattern) = build_store_record_columns(store, options, temp_dir)?;
        let payload_path = temp_dir.join("spring-rs3.payload");
        let payload_file = std::fs::OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(&payload_path)?;
        let mut payload = std::io::BufWriter::with_capacity(16 * 1024 * 1024, payload_file);
        put_u64_writer(&mut payload, store.count())?;
        put_u16_writer(&mut payload, 2)?;
        payload.write_all(&[1])?;
        payload.write_all(&[id_pattern as u8])?;
        for (mate, (assignment_path, columns)) in assignments.iter().zip(record_columns).enumerate()
        {
            let stats = write_store_clusters(
                &mut payload,
                store,
                assignment_path,
                mate,
                &options.reorder,
                options.threads,
                temp_dir,
            )?;
            store.add_stats(stats);
            let ids_stored =
                options.preserve_ids && (mate == 0 || id_pattern == PairedIdPattern::Independent);
            write_store_record_columns(
                &mut payload,
                columns,
                store.count(),
                ids_stored,
                options.preserve_quality,
                65_536,
            )?;
        }
        payload.flush()?;
        let mut payload_file = payload
            .into_inner()
            .map_err(|error| Error::Io(error.into_error()))?;
        let raw_len = payload_file.seek(SeekFrom::End(0))?;
        payload_file.seek(SeekFrom::Start(0))?;
        let payload_time = payload_start.elapsed();
        let checksum_start = std::time::Instant::now();
        let checksum = checksum_file(std::fs::File::open(&payload_path)?)?;
        let checksum_time = checksum_start.elapsed();
        let encoded_path = temp_dir.join("spring-rs3.payload.zst");
        let encoded_file = std::fs::OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(&encoded_path)?;
        let zstd_start = std::time::Instant::now();
        let mut encoder = zstd::stream::Encoder::new(encoded_file, options.zstd_level)?;
        encoder.multithread(worker_count(options.threads)?)?;
        if raw_len >= 1 << LONG_DISTANCE_WINDOW_LOG {
            encoder.set_parameter(zstd::zstd_safe::CParameter::WindowLog(
                LONG_DISTANCE_WINDOW_LOG,
            ))?;
            encoder.set_parameter(zstd::zstd_safe::CParameter::EnableLongDistanceMatching(
                true,
            ))?;
        }
        encoder.set_pledged_src_size(Some(raw_len))?;
        let mut transfer_buffer = vec![0; 16 * 1024 * 1024];
        let mut payload_reader = std::fs::File::open(&payload_path)?;
        loop {
            let bytes_read = payload_reader.read(&mut transfer_buffer)?;
            if bytes_read == 0 {
                break;
            }
            encoder.write_all(&transfer_buffer[..bytes_read])?;
        }
        let mut encoded_file = encoder.finish()?;
        encoded_file.flush()?;
        let encoded_len = encoded_file.seek(SeekFrom::End(0))?;
        let zstd_time = zstd_start.elapsed();
        encoded_file.seek(SeekFrom::Start(0))?;
        output.write_all(MAGIC)?;
        put_u32_writer(&mut output, VERSION)?;
        put_u16_writer(&mut output, 2)?;
        put_u16_writer(&mut output, flags(options))?;
        put_u64_writer(&mut output, raw_len)?;
        put_u64_writer(&mut output, encoded_len)?;
        put_u32_writer(&mut output, checksum)?;
        std::io::copy(&mut encoded_file, &mut output)?;
        let metadata = CompressionMetadata {
            streams: 2,
            records: store.count(),
            preserved_ids: options.preserve_ids,
            preserved_quality: options.preserve_quality,
            preserved_order: true,
        };
        let profile = CompressionProfile {
            assignments: assignment_time,
            consensus_and_payload: payload_time,
            checksum: checksum_time,
            zstd: zstd_time,
            store_stats: store.stats(),
            assignment_bytes,
            payload_bytes: raw_len,
            archive_bytes: 37 + encoded_len,
            ..CompressionProfile::default()
        };
        Ok((metadata, profile))
    }
}

#[cfg(test)]
fn write_reordered_store_profiled_v2<W: Write>(
    store: &mut PairStore,
    mut output: W,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<(CompressionMetadata, CompressionProfile)> {
    let temp_dir = temp_dir.as_ref();
    let mut metadata_builder = PreparedRecordsBuilder::create_metadata(temp_dir)?;
    let mut reads = [Vec::new(), Vec::new()];
    for id in 0..store.count() {
        let (first, second) = store.read_pair(id)?;
        let prepared = prepare_pair(first, second, id, options)?;
        let sorted = metadata_builder.push_metadata(prepared, options)?;
        reads[0].push(sorted[0].clone());
        reads[1].push(sorted[1].clone());
    }
    let metadata = metadata_builder.finish_metadata()?;
    let sections = [0usize, 1].map(|mate| {
        let (tx, rx) = crossbeam_channel::bounded(reads[mate].len().max(1));
        for read in std::mem::take(&mut reads[mate]) {
            tx.send(read)
                .map_err(|_| Error::InvalidInput("sorter exited".into()))?;
        }
        drop(tx);
        sort_and_write_mate(rx, store.count(), mate, options, temp_dir)
    });
    let sections = match sections {
        [Ok(a), Ok(b)] => [a, b],
        [Err(e), _] | [_, Err(e)] => return Err(e),
    };
    let mut archive = std::io::Cursor::new(Vec::new());
    let (metadata, mut profile) =
        write_streamed_profiled(metadata, sections, &mut archive, options, temp_dir)?;
    output.write_all(archive.get_ref())?;
    profile.store_stats = StoreStats {
        index_reads: store.count().saturating_mul(5),
        blob_reads: store.count().saturating_mul(8),
        bytes_read: store.stats().bytes_read,
    };
    Ok((metadata, profile))
}

#[allow(dead_code)]
pub(crate) fn write_prepared_profiled<W: Write>(
    prepared: PreparedRecords,
    mut output: W,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<(CompressionMetadata, CompressionProfile)> {
    return write_prepared_profiled_v2(prepared, output, options, temp_dir);
    #[allow(unreachable_code)]
    {
        let temp_dir = temp_dir.as_ref();
        let PreparedRecords {
            columns,
            spool_paths,
            max_sortable_bytes,
            id_pattern,
            count,
        } = prepared;
        let assignment_bytes = spool_paths.iter().try_fold(0u64, |total, path| {
            Ok::<_, Error>(total + std::fs::metadata(path)?.len())
        })?;
        let payload_path = temp_dir.join("spring-rs3.payload");
        let payload_file = std::fs::OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(&payload_path)?;
        let mut payload = std::io::BufWriter::with_capacity(16 * 1024 * 1024, payload_file);
        put_u64_writer(&mut payload, count)?;
        put_u16_writer(&mut payload, 2)?;
        payload.write_all(&[1])?;
        payload.write_all(&[id_pattern as u8])?;

        let mut assignment_time = std::time::Duration::default();
        let mut consensus_metrics = ConsensusMetrics::default();
        let payload_start = std::time::Instant::now();
        for (mate, ((spool_path, columns), max_sortable_bytes)) in spool_paths
            .iter()
            .zip(columns)
            .zip(max_sortable_bytes)
            .enumerate()
        {
            const SORT_MEMORY_BYTES: usize = 512 * 1024 * 1024;
            let memory_limited_records = SORT_MEMORY_BYTES / max_sortable_bytes.max(1);
            let sort_chunk_size = options
                .reorder
                .sort_chunk_size
                .min(memory_limited_records.max(1))
                .max(1);
            let sorter = bed_utils::extsort::ExternalSorterBuilder::new()
                .with_chunk_size(sort_chunk_size)
                .with_tmp_dir(temp_dir)
                .with_compression(3)
                .num_threads(options.reorder.sort_threads.max(1))
                .build()
                .map_err(Error::Io)?;
            let mut reader = SequenceSpoolReader::open(spool_path)?;
            if reader.record_count() != count {
                return Err(Error::InvalidFormat(
                    "sequence spool record count mismatch".into(),
                ));
            }
            let source_error = std::rc::Rc::new(std::cell::RefCell::new(None));
            let input_error = source_error.clone();
            let input = std::iter::from_fn(move || match reader.next() {
                Some(Ok(read)) => Some(read),
                Some(Err(error)) => {
                    *input_error.borrow_mut() = Some(error);
                    None
                }
                None => None,
            });
            let sort_start = std::time::Instant::now();
            let sorted = sorter
                .sort_by(input, SortableRead::compare)
                .map_err(|error| Error::InvalidInput(error.to_string()))?;
            assignment_time += sort_start.elapsed();
            if let Some(error) = source_error.borrow_mut().take() {
                return Err(error);
            }
            let order_path = temp_dir.join(format!("prepared.{mate}.order"));
            let mut order = std::fs::File::create(&order_path)?;
            order.set_len(count.saturating_mul(ORDER_ENTRY_SIZE as u64))?;
            let metrics = write_sorted_clusters(
                &mut payload,
                sorted,
                mate,
                &options.reorder,
                options.threads,
                temp_dir,
                count,
                &mut order,
            )?;
            order.flush()?;
            append_compact_order(&mut payload, &order_path, count, temp_dir)?;
            std::fs::remove_file(order_path)?;
            consensus_metrics.merge(metrics);
            let ids_stored =
                options.preserve_ids && (mate == 0 || id_pattern == PairedIdPattern::Independent);
            write_prepared_record_columns(
                &mut payload,
                columns,
                count,
                ids_stored,
                options.preserve_quality,
            )?;
            std::fs::remove_file(spool_path)?;
        }
        payload.flush()?;
        let mut payload_file = payload
            .into_inner()
            .map_err(|error| Error::Io(error.into_error()))?;
        let raw_len = payload_file.seek(SeekFrom::End(0))?;
        payload_file.seek(SeekFrom::Start(0))?;
        let payload_time = payload_start.elapsed();
        let checksum_start = std::time::Instant::now();
        let checksum = checksum_file(std::fs::File::open(&payload_path)?)?;
        let checksum_time = checksum_start.elapsed();
        let encoded_path = temp_dir.join("spring-rs3.payload.zst");
        let encoded_file = std::fs::OpenOptions::new()
            .create(true)
            .read(true)
            .write(true)
            .truncate(true)
            .open(&encoded_path)?;
        let zstd_start = std::time::Instant::now();
        let mut encoder = zstd::stream::Encoder::new(encoded_file, options.zstd_level)?;
        encoder.multithread(worker_count(options.threads)?)?;
        if raw_len >= 1 << LONG_DISTANCE_WINDOW_LOG {
            encoder.set_parameter(zstd::zstd_safe::CParameter::WindowLog(
                LONG_DISTANCE_WINDOW_LOG,
            ))?;
            encoder.set_parameter(zstd::zstd_safe::CParameter::EnableLongDistanceMatching(
                true,
            ))?;
        }
        encoder.set_pledged_src_size(Some(raw_len))?;
        let mut transfer_buffer = vec![0; 16 * 1024 * 1024];
        let mut payload_reader = std::fs::File::open(&payload_path)?;
        loop {
            let bytes_read = payload_reader.read(&mut transfer_buffer)?;
            if bytes_read == 0 {
                break;
            }
            encoder.write_all(&transfer_buffer[..bytes_read])?;
        }
        let mut encoded_file = encoder.finish()?;
        encoded_file.flush()?;
        let encoded_len = encoded_file.seek(SeekFrom::End(0))?;
        let zstd_time = zstd_start.elapsed();
        encoded_file.seek(SeekFrom::Start(0))?;
        output.write_all(MAGIC)?;
        put_u32_writer(&mut output, VERSION)?;
        put_u16_writer(&mut output, 2)?;
        put_u16_writer(&mut output, flags(options))?;
        put_u64_writer(&mut output, raw_len)?;
        put_u64_writer(&mut output, encoded_len)?;
        put_u32_writer(&mut output, checksum)?;
        std::io::copy(&mut encoded_file, &mut output)?;
        let metadata = CompressionMetadata {
            streams: 2,
            records: count,
            preserved_ids: options.preserve_ids,
            preserved_quality: options.preserve_quality,
            preserved_order: true,
        };
        let profile = CompressionProfile {
            assignments: assignment_time,
            consensus_and_payload: payload_time.saturating_sub(assignment_time),
            checksum: checksum_time,
            zstd: zstd_time,
            assignment_bytes,
            payload_bytes: raw_len,
            archive_bytes: 37 + encoded_len,
            consensus_producer: consensus_metrics.producer,
            consensus_producer_wait: consensus_metrics.producer_wait,
            consensus_workers: consensus_metrics.workers,
            consensus_writer: consensus_metrics.writer,
            consensus_writer_active: consensus_metrics.writer_active,
            consensus_writer_wait: consensus_metrics.writer_wait,
            consensus_max_pending_buckets: consensus_metrics.max_pending_buckets,
            consensus_multi_clusters: consensus_metrics.multi_clusters,
            consensus_singletons: consensus_metrics.singletons,
            consensus_buckets: consensus_metrics.buckets,
            consensus_overlap_comparisons: consensus_metrics.overlap_comparisons,
            consensus_max_bucket_records: consensus_metrics.max_bucket_records,
            consensus_max_groups: consensus_metrics.max_groups,
            ..CompressionProfile::default()
        };
        Ok((metadata, profile))
    }
}

#[allow(dead_code)]
fn write_prepared_profiled_v2<W: Write>(
    prepared: PreparedRecords,
    mut output: W,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<(CompressionMetadata, CompressionProfile)> {
    let temp_dir = temp_dir.as_ref();
    let PreparedRecords {
        columns,
        spool_paths,
        id_pattern,
        count,
        ..
    } = prepared;
    let (ids0, qualities0) = load_prepared_metadata(&columns[0])?;
    let (ids1, qualities1) = load_prepared_metadata(&columns[1])?;
    let metadata = PreparedMetadata {
        columns,
        id_pattern,
        count,
        ids: [ids0, ids1],
        qualities: [qualities0, qualities1],
    };
    let sections = [0usize, 1].map(|mate| -> Result<PreparedMateSection> {
        let mut reader = SequenceSpoolReader::open(&spool_paths[mate])?;
        let mut reads = Vec::new();
        while let Some(read) = reader.next() {
            reads.push(read?);
        }
        let (tx, rx) = crossbeam_channel::bounded(reads.len().max(1));
        for read in reads {
            tx.send(read)
                .map_err(|_| Error::InvalidInput("sorter exited".into()))?;
        }
        drop(tx);
        sort_and_write_mate(rx, count, mate, options, temp_dir)
    });
    let sections = match sections {
        [Ok(first), Ok(second)] => [first, second],
        [Err(error), _] | [_, Err(error)] => return Err(error),
    };
    let mut archive = std::io::Cursor::new(Vec::new());
    let result = write_streamed_profiled(metadata, sections, &mut archive, options, temp_dir);
    if result.is_ok() {
        output.write_all(archive.get_ref())?;
    }
    for path in spool_paths {
        let _ = std::fs::remove_file(path);
    }
    result
}

fn worker_count(workers: usize) -> Result<u32> {
    u32::try_from(workers).map_err(|_| Error::InvalidInput("thread count exceeds u32".into()))
}

pub(crate) struct PreparedMateSection {
    path: std::path::PathBuf,
    order_path: std::path::PathBuf,
    assignment_time: std::time::Duration,
    metrics: ConsensusMetrics,
}

pub(crate) fn sort_and_write_mate(
    input: crossbeam_channel::Receiver<SortableRead>,
    expected_count: u64,
    mate: usize,
    options: &CompressionOptions,
    temp_dir: &Path,
) -> Result<PreparedMateSection> {
    let sorter = bed_utils::extsort::ExternalSorterBuilder::new()
        .with_chunk_size(options.reorder.sort_chunk_size.max(1))
        .with_tmp_dir(temp_dir)
        .with_compression(3)
        .num_threads(options.reorder.sort_threads.div_ceil(2).max(1))
        .build()
        .map_err(Error::Io)?;
    let sort_started = std::time::Instant::now();
    let sorted = sorter
        .sort_by(input, SortableRead::compare)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let assignment_time = sort_started.elapsed();
    let path = temp_dir.join(format!("mate.{mate}.section"));
    let order_path = temp_dir.join(format!("mate.{mate}.order"));
    let mut order = std::fs::File::create(&order_path)?;
    let file = std::fs::File::create(&path)?;
    let mut output = std::io::BufWriter::with_capacity(16 * 1024 * 1024, file);
    let metrics = write_sorted_clusters(
        &mut output,
        sorted,
        mate,
        &options.reorder,
        options.threads.div_ceil(2).max(1),
        temp_dir,
        expected_count,
        &mut order,
    )?;
    order.set_len(
        metrics
            .source_records
            .saturating_mul(ORDER_ENTRY_SIZE as u64),
    )?;
    output.flush()?;
    order.flush()?;
    append_compact_order(&mut output, &order_path, metrics.source_records, temp_dir)?;
    output.flush()?;
    Ok(PreparedMateSection {
        path,
        order_path,
        assignment_time,
        metrics,
    })
}

pub(crate) fn write_streamed_profiled<W: Write + Seek>(
    prepared: PreparedMetadata,
    sections: [PreparedMateSection; 2],
    mut output: W,
    options: &CompressionOptions,
    _temp_dir: &Path,
) -> Result<(CompressionMetadata, CompressionProfile)> {
    let PreparedMetadata {
        id_pattern,
        count,
        ids,
        qualities,
        columns: _,
    } = prepared;
    output.write_all(&[0; 37])?;
    let encoded_start = output.stream_position()?;
    let zstd_started = std::time::Instant::now();
    let mut encoder = zstd::stream::Encoder::new(&mut output, options.zstd_level)?;
    encoder.multithread(worker_count(options.threads)?)?;
    encoder.set_parameter(zstd::zstd_safe::CParameter::WindowLog(
        LONG_DISTANCE_WINDOW_LOG,
    ))?;
    encoder.set_parameter(zstd::zstd_safe::CParameter::EnableLongDistanceMatching(
        true,
    ))?;
    let mut payload = PayloadChecksumWriter::new(encoder);
    put_u64_writer(&mut payload, count)?;
    put_u16_writer(&mut payload, 2)?;
    payload.write_all(&[id_pattern as u8])?;
    let payload_started = std::time::Instant::now();
    let mut assignment_time = std::time::Duration::default();
    let mut consensus_metrics = ConsensusMetrics::default();
    for (mate, section) in sections.into_iter().enumerate() {
        if section.metrics.source_records != count {
            return Err(Error::InvalidFormat(
                "sorted mate stream record count mismatch".into(),
            ));
        }
        assignment_time += section.assignment_time;
        consensus_metrics.merge(section.metrics);
        let mut section_file = std::io::BufReader::new(std::fs::File::open(&section.path)?);
        write_singleton_tail_section(
            &mut section_file,
            &mut payload,
            count,
            mate,
            &ids[mate],
            &qualities[mate],
            options,
            id_pattern,
        )?;
        std::fs::remove_file(section.path)?;
        std::fs::remove_file(section.order_path)?;
    }
    payload.flush()?;
    let (encoder, raw_len, checksum) = payload.finish();
    let output_ref = encoder.finish()?;
    output_ref.flush()?;
    let payload_time = payload_started.elapsed();
    let zstd_time = zstd_started.elapsed();
    let encoded_end = output.stream_position()?;
    let encoded_len = encoded_end
        .checked_sub(encoded_start)
        .ok_or_else(|| Error::InvalidInput("invalid encoded payload length".into()))?;
    output.seek(SeekFrom::Start(0))?;
    output.write_all(MAGIC)?;
    put_u32_writer(&mut output, VERSION)?;
    put_u16_writer(&mut output, 2)?;
    put_u16_writer(&mut output, flags(options))?;
    put_u64_writer(&mut output, raw_len)?;
    put_u64_writer(&mut output, encoded_len)?;
    put_u32_writer(&mut output, checksum)?;
    output.seek(SeekFrom::Start(encoded_end))?;
    Ok((
        CompressionMetadata {
            streams: 2,
            records: count,
            preserved_ids: options.preserve_ids,
            preserved_quality: options.preserve_quality,
            preserved_order: true,
        },
        CompressionProfile {
            assignments: assignment_time,
            consensus_and_payload: payload_time,
            checksum: std::time::Duration::default(),
            zstd: zstd_time,
            payload_bytes: raw_len,
            archive_bytes: 37 + encoded_len,
            consensus_producer: consensus_metrics.producer,
            consensus_producer_wait: consensus_metrics.producer_wait,
            consensus_workers: consensus_metrics.workers,
            consensus_writer: consensus_metrics.writer,
            consensus_writer_active: consensus_metrics.writer_active,
            consensus_writer_wait: consensus_metrics.writer_wait,
            consensus_max_pending_buckets: consensus_metrics.max_pending_buckets,
            consensus_multi_clusters: consensus_metrics.multi_clusters,
            consensus_singletons: consensus_metrics.singletons,
            consensus_buckets: consensus_metrics.buckets,
            consensus_overlap_comparisons: consensus_metrics.overlap_comparisons,
            consensus_max_bucket_records: consensus_metrics.max_bucket_records,
            consensus_max_groups: consensus_metrics.max_groups,
            ..CompressionProfile::default()
        },
    ))
}

struct SingletonTailCounts {
    multi_groups: u64,
    multi_records: u64,
    singleton_buckets: u64,
    singleton_records: u64,
    keys: Vec<(u32, bool)>,
}

fn write_singleton_tail_section<R: Read + Seek, W: Write>(
    input: &mut R,
    output: &mut W,
    count: u64,
    mate: usize,
    ids: &[Vec<u8>],
    qualities: &[Vec<u8>],
    options: &CompressionOptions,
    id_pattern: PairedIdPattern,
) -> Result<()> {
    let counts = scan_singleton_tail(input, count)?;
    put_u64_writer(output, counts.multi_groups)?;
    put_u64_writer(output, counts.multi_records)?;
    put_u64_writer(output, counts.singleton_buckets)?;
    put_u64_writer(output, counts.singleton_records)?;

    input.seek(SeekFrom::Start(0))?;
    let mut records = 0u64;
    while records < count {
        let group_records = read_varint(input)?;
        let kind = read_one(input)?;
        if kind != 0 {
            skip_source_singleton(input, group_records)?;
            records += group_records;
            continue;
        }
        let reference = read_reference_stream(input)?;
        put_varint(output, group_records)?;
        output.write_all(&[0])?;
        put_reference_v8(output, &reference)?;
        let mut group_ids = Vec::with_capacity(usize_from(group_records, "group")?);
        let mut previous_id = 0u64;
        for _ in 0..group_records {
            let delta = read_varint(input)?;
            let original_id = zigzag_decode(delta, previous_id)?;
            previous_id = original_id;
            group_ids.push(original_id);
            put_varint(output, delta)?;
            output.write_all(&[read_one(input)?])?;
            put_varint(output, read_varint(input)?)?;
            put_varint(output, read_varint(input)?)?;
            let mismatch_count = read_varint(input)?;
            put_varint(output, mismatch_count)?;
            for _ in 0..mismatch_count {
                put_varint(output, read_varint(input)?)?;
                output.write_all(&[read_one(input)?])?;
            }
        }
        write_group_metadata(
            output, &group_ids, ids, qualities, options, mate, id_pattern,
        )?;
        records += group_records;
    }

    for &(length, packed) in &counts.keys {
        input.seek(SeekFrom::Start(0))?;
        let mut pending = None;
        let mut records = 0u64;
        while records < count {
            let group_records = read_varint(input)?;
            let kind = read_one(input)?;
            if kind == 0 {
                skip_source_multi(input, group_records)?;
            } else {
                if group_records != 1 {
                    return Err(Error::InvalidFormat(
                        "invalid singleton group record count".into(),
                    ));
                }
                let original_id = zigzag_decode(read_varint(input)?, 0)?;
                let sequence = read_packed_sequence(input)?;
                if sequence.len() == length && sequence.is_packed() == packed {
                    let fits = pending
                        .as_ref()
                        .is_some_and(|batch: &PendingSingletonBatch| {
                            batch.ids.len() < MAX_BATCH_RECORDS
                                && batch.body_len.saturating_add(sequence.encoded_size())
                                    <= MAX_BATCH_BYTES
                        });
                    if !fits {
                        flush_singleton_tail_batch(
                            output,
                            pending.take(),
                            ids,
                            qualities,
                            options,
                            mate,
                            id_pattern,
                        )?;
                        pending = Some(PendingSingletonBatch {
                            ids: Vec::with_capacity(MAX_BATCH_RECORDS),
                            bodies: Vec::with_capacity(MAX_BATCH_RECORDS),
                            length,
                            packed,
                            body_len: 0,
                        });
                    }
                    let batch = pending.as_mut().unwrap();
                    batch.ids.push(original_id);
                    batch.body_len += sequence.encoded_size();
                    batch.bodies.push(sequence.encoded_bytes().to_vec());
                }
            }
            records += group_records;
        }
        flush_singleton_tail_batch(output, pending, ids, qualities, options, mate, id_pattern)?;
    }
    Ok(())
}

const MAX_BATCH_RECORDS: usize = 4096;
const MAX_BATCH_BYTES: usize = 4 * 1024 * 1024;

fn scan_singleton_tail<R: Read + Seek>(input: &mut R, count: u64) -> Result<SingletonTailCounts> {
    input.seek(SeekFrom::Start(0))?;
    let mut records = 0u64;
    let mut multi_groups = 0;
    let mut multi_records = 0;
    let mut singleton_buckets = 0;
    let mut singleton_records = 0;
    let mut singleton_batches = BTreeMap::<(u32, bool), (usize, usize)>::new();
    while records < count {
        let group_records = read_varint(input)?;
        let kind = read_one(input)?;
        if group_records == 0 || records + group_records > count {
            return Err(Error::InvalidFormat(
                "invalid RS4 source group count".into(),
            ));
        }
        if kind == 0 {
            skip_source_multi(input, group_records)?;
            multi_groups += 1;
            multi_records += group_records;
        } else if kind == 1 && group_records == 1 {
            let _ = read_varint(input)?;
            let sequence = read_packed_sequence(input)?;
            let key = (sequence.len(), sequence.is_packed());
            let batch = singleton_batches.entry(key).or_default();
            if batch.0 != 0
                && (batch.0 >= MAX_BATCH_RECORDS
                    || batch.1.saturating_add(sequence.encoded_size()) > MAX_BATCH_BYTES)
            {
                singleton_buckets += 1;
                *batch = (0, 0);
            }
            batch.0 += 1;
            batch.1 = batch.1.saturating_add(sequence.encoded_size());
            singleton_records += 1;
        } else {
            return Err(Error::InvalidFormat("invalid RS4 source group kind".into()));
        }
        records += group_records;
    }
    singleton_buckets += singleton_batches
        .values()
        .filter(|(count, _)| *count != 0)
        .count() as u64;
    Ok(SingletonTailCounts {
        multi_groups,
        multi_records,
        singleton_buckets,
        singleton_records,
        keys: singleton_batches.keys().copied().collect(),
    })
}

fn skip_source_singleton<R: Read>(input: &mut R, count: u64) -> Result<()> {
    if count != 1 {
        return Err(Error::InvalidFormat(
            "invalid singleton group record count".into(),
        ));
    }
    let _ = read_varint(input)?;
    let _ = read_packed_sequence(input)?;
    Ok(())
}

fn skip_source_multi<R: Read>(input: &mut R, count: u64) -> Result<()> {
    let _ = read_reference_stream(input)?;
    for _ in 0..count {
        let _ = read_varint(input)?;
        let _ = read_one(input)?;
        let _ = read_varint(input)?;
        let _ = read_varint(input)?;
        let mismatch_count = read_varint(input)?;
        for _ in 0..mismatch_count {
            let _ = read_varint(input)?;
            let _ = read_one(input)?;
        }
    }
    Ok(())
}

#[cfg(test)]
fn write_store_clusters<W: Write>(
    payload: &mut W,
    store: &PairStore,
    assignments: &Path,
    mate: usize,
    reorder: &crate::ReorderOptions,
    threads: usize,
    temp_dir: &Path,
) -> Result<StoreStats> {
    let limit = reorder.max_bucket_records.max(1);
    let body_path = temp_dir.join(format!("clusters.{mate}.multi.body"));
    let mut multi_body =
        std::io::BufWriter::with_capacity(8 * 1024 * 1024, std::fs::File::create(&body_path)?);
    let mut singleton_sequences = StoredColumn::create(temp_dir, mate, "singleton-sequences")?;
    let order_path = temp_dir.join(format!("clusters.{mate}.order"));
    let mut order = std::fs::File::create(&order_path)?;
    order.set_len(store.count().saturating_mul(ORDER_ENTRY_SIZE as u64))?;
    let workers = threads.max(1);
    let queue_capacity = workers.div_ceil(2).max(1);
    let in_flight = workers + queue_capacity;
    let (job_tx, job_rx) = bounded::<(usize, Vec<AssignedRead>)>(queue_capacity);
    let (done_tx, done_rx) = bounded::<Result<ConsensusBucket>>(queue_capacity);
    let (permit_tx, permit_rx) = bounded::<()>(in_flight);
    for _ in 0..in_flight {
        permit_tx
            .send(())
            .map_err(|_| Error::InvalidInput("failed to initialize consensus pipeline".into()))?;
    }
    let worker_stores = (0..workers)
        .map(|_| PairStore::open(store.root()))
        .collect::<Result<Vec<_>>>()?;

    let (cluster_count, _pending, stats) =
        std::thread::scope(|scope| -> Result<((u32, u32), u64, StoreStats)> {
            let producer_errors = done_tx.clone();
            scope.spawn(move || {
                let result = produce_consensus_buckets(assignments, limit, permit_rx, job_tx);
                if let Err(error) = result {
                    let _ = producer_errors.send(Err(error));
                }
            });

            for mut worker_store in worker_stores {
                let job_rx = job_rx.clone();
                let done_tx = done_tx.clone();
                scope.spawn(move || {
                    while let Ok((index, ids)) = job_rx.recv() {
                        worker_store.reset_stats();
                        let result =
                            cluster_bucket(&mut worker_store, mate, reorder.kmer_size, &ids)
                                .and_then(|clusters| {
                                    let (
                                        multi_bytes,
                                        multi_cluster_count,
                                        singletons,
                                        order_entries,
                                    ) = partition_clusters_v8(&clusters)?;
                                    Ok(ConsensusBucket {
                                        index,
                                        multi_cluster_count,
                                        multi_bytes,
                                        singletons,
                                        order_entries,
                                        stats: worker_store.stats(),
                                        overlap_comparisons: 0,
                                        source_records: ids.len() as u64,
                                        groups: clusters.len() as u64,
                                        consensus_groups: clusters,
                                    })
                                });
                        if done_tx.send(result).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(job_rx);
            drop(done_tx);

            let result = write_consensus_buckets(
                done_rx,
                permit_tx,
                &mut multi_body,
                &mut singleton_sequences,
                &mut order,
            )?;
            Ok((result.counts, result.max_pending_buckets, result.stats))
        })?;
    multi_body.flush()?;
    drop(multi_body);
    put_u32_writer(payload, cluster_count.0)?;
    let mut multi_body = std::fs::File::open(body_path)?;
    std::io::copy(&mut multi_body, payload)?;
    put_u32_writer(payload, cluster_count.1)?;
    singleton_sequences.write_packed_dna_to(payload, cluster_count.1 as usize)?;
    order.flush()?;
    append_compact_order(payload, &order_path, store.count(), temp_dir)?;
    std::fs::remove_file(order_path)?;
    Ok(stats)
}

struct ConsensusBucket {
    index: usize,
    multi_cluster_count: u32,
    multi_bytes: Vec<u8>,
    singletons: Vec<(u64, PackedSequence)>,
    order_entries: Vec<(u64, OrderEntry)>,
    #[cfg(test)]
    stats: StoreStats,
    overlap_comparisons: u64,
    source_records: u64,
    groups: u64,
    consensus_groups: Vec<ConsensusCluster>,
}

#[derive(Default)]
struct ConsensusMetrics {
    producer: std::time::Duration,
    producer_wait: std::time::Duration,
    workers: std::time::Duration,
    writer: std::time::Duration,
    writer_active: std::time::Duration,
    writer_wait: std::time::Duration,
    max_pending_buckets: u64,
    multi_clusters: u64,
    singletons: u64,
    buckets: u64,
    overlap_comparisons: u64,
    max_bucket_records: u64,
    max_groups: u64,
    source_records: u64,
}

impl ConsensusMetrics {
    fn merge(&mut self, other: Self) {
        self.producer += other.producer;
        self.producer_wait += other.producer_wait;
        self.workers += other.workers;
        self.writer += other.writer;
        self.writer_active += other.writer_active;
        self.writer_wait += other.writer_wait;
        self.max_pending_buckets = self.max_pending_buckets.max(other.max_pending_buckets);
        self.multi_clusters += other.multi_clusters;
        self.singletons += other.singletons;
        self.buckets += other.buckets;
        self.overlap_comparisons += other.overlap_comparisons;
        self.max_bucket_records = self.max_bucket_records.max(other.max_bucket_records);
        self.max_groups = self.max_groups.max(other.max_groups);
        self.source_records += other.source_records;
    }
}

#[cfg(test)]
#[derive(Clone, Copy)]
struct AssignedRead {
    id: u64,
    position: u32,
    reverse: bool,
    anchored: bool,
}

struct AssignedSequence {
    id: u64,
    position: u32,
    reverse: bool,
    anchored: bool,
    sequence: Vec<u8>,
}

#[cfg(test)]
trait SingletonSequenceSink {
    fn push_sequence(&mut self, sequence: PackedSequence) -> Result<()>;
}

#[cfg(test)]
#[cfg(test)]
impl SingletonSequenceSink for StoredColumn {
    fn push_sequence(&mut self, sequence: PackedSequence) -> Result<()> {
        self.push(&sequence.decode())
    }
}

#[derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize, Clone, Debug)]
struct SingletonRecord {
    id: u64,
    sequence: PackedSequence,
}

struct SingletonRecordWriter {
    path: std::path::PathBuf,
    builder: ExternalChunkBuilder<SingletonRecord>,
    count: usize,
    common_length: Option<u32>,
    variable_length: bool,
    data_len: u64,
    clean_records: usize,
}

impl SingletonRecordWriter {
    fn create(temp_dir: &Path, mate: usize) -> Result<Self> {
        let path = temp_dir.join(format!("records.{mate}.singletons.lz4"));
        let builder = ExternalChunkBuilder::new(std::fs::File::create(&path)?, 1)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        Ok(Self {
            path,
            builder,
            count: 0,
            common_length: None,
            variable_length: false,
            data_len: 0,
            clean_records: 0,
        })
    }

    fn push_singleton(&mut self, id: u64, sequence: PackedSequence) -> Result<()> {
        let length = sequence.len();
        self.variable_length |= self.common_length.is_some_and(|common| common != length);
        self.common_length.get_or_insert(length);
        self.data_len = self
            .data_len
            .checked_add(length as u64)
            .ok_or_else(|| Error::InvalidInput("singleton sequence column too large".into()))?;
        self.clean_records += sequence.is_packed() as usize;
        self.count = self
            .count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many singleton sequences".into()))?;
        self.builder
            .add(SingletonRecord { id, sequence })
            .map_err(|error| Error::InvalidInput(error.to_string()))
    }

    fn finish_to<W: Write>(self, payload: &mut W, expected_count: usize) -> Result<()> {
        let Self {
            path,
            builder,
            count,
            common_length,
            variable_length,
            data_len,
            clean_records,
        } = self;
        if count != expected_count {
            return Err(Error::InvalidInput(
                "singleton sequence count mismatch".into(),
            ));
        }
        let chunk = builder
            .finish()
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        drop(chunk);
        put_u32_writer(payload, count as u32)?;
        if variable_length {
            payload.write_all(&[2])?;
            put_u64_writer(payload, data_len)?;
            for_each_singleton(&path, count, |record| {
                put_u32_writer(payload, record.sequence.len())
            })?;
            for_each_singleton(&path, count, |record| {
                payload.write_all(&record.sequence.decode())?;
                Ok(())
            })?;
        } else {
            let width = common_length.unwrap_or(0) as usize;
            let packed_width = width.div_ceil(4);
            let clean = clean_records == count;
            if clean {
                payload.write_all(&[4])?;
                put_u32_writer(payload, width as u32)?;
                put_u64_writer(payload, (packed_width as u64) * count as u64)?;
            } else {
                let ambiguous = count - clean_records;
                let encoded_len =
                    count
                        .checked_add(clean_records.checked_mul(packed_width).ok_or_else(|| {
                            Error::InvalidInput("mixed DNA column too large".into())
                        })?)
                        .and_then(|length| length.checked_add(ambiguous.checked_mul(width)?))
                        .ok_or_else(|| Error::InvalidInput("mixed DNA column too large".into()))?;
                payload.write_all(&[5])?;
                put_u32_writer(payload, width as u32)?;
                put_u64_writer(payload, encoded_len as u64)?;
            }
            for_each_singleton(&path, count, |record| {
                if !clean {
                    payload.write_all(&[(!record.sequence.is_packed()) as u8])?;
                }
                if record.sequence.is_packed() {
                    payload.write_all(record.sequence.encoded_bytes())?;
                } else {
                    payload.write_all(&record.sequence.decode())?;
                }
                Ok(())
            })?;
        }
        std::fs::remove_file(path)?;
        Ok(())
    }
}

fn for_each_singleton<F>(path: &Path, expected_count: usize, mut visit: F) -> Result<()>
where
    F: FnMut(SingletonRecord) -> Result<()>,
{
    let chunk = ExternalChunk::<SingletonRecord>::open(std::fs::File::open(path)?)
        .map_err(|error| Error::InvalidFormat(error.to_string()))?;
    let mut count = 0usize;
    for record in chunk {
        let record = record.map_err(|error| Error::InvalidFormat(error.to_string()))?;
        visit(record)?;
        count += 1;
    }
    if count != expected_count {
        return Err(Error::InvalidFormat(
            "singleton chunk count mismatch".into(),
        ));
    }
    Ok(())
}

fn write_sorted_clusters<W, I, E>(
    payload: &mut W,
    sorted: I,
    mate: usize,
    reorder: &crate::ReorderOptions,
    threads: usize,
    temp_dir: &Path,
    expected_count: u64,
    order_output: &mut std::fs::File,
) -> Result<ConsensusMetrics>
where
    W: Write + Seek,
    I: Iterator<Item = std::result::Result<SortableRead, E>> + Send,
    E: std::fmt::Display + Send,
{
    let limit = reorder.max_bucket_records.max(1);
    let workers = threads.max(1);
    let queue_capacity = workers.div_ceil(2).max(1);
    let in_flight = workers + queue_capacity;
    let (job_tx, job_rx) = bounded::<(usize, Vec<AssignedSequence>)>(queue_capacity);
    let (done_tx, done_rx) = bounded::<Result<ConsensusBucket>>(queue_capacity);
    let (permit_tx, permit_rx) = bounded::<()>(in_flight);
    for _ in 0..in_flight {
        permit_tx
            .send(())
            .map_err(|_| Error::InvalidInput("failed to initialize consensus pipeline".into()))?;
    }

    let worker_nanos = Arc::new(AtomicU64::new(0));
    let producer_nanos = Arc::new(AtomicU64::new(0));
    let producer_wait_nanos = Arc::new(AtomicU64::new(0));
    let (cluster_count, metrics) =
        std::thread::scope(|scope| -> Result<((u32, u32), ConsensusMetrics)> {
            let producer_errors = done_tx.clone();
            let producer_metric = producer_nanos.clone();
            let producer_wait_metric = producer_wait_nanos.clone();
            scope.spawn(move || {
                let started = std::time::Instant::now();
                let result = produce_sorted_consensus_buckets(
                    sorted,
                    limit,
                    permit_rx,
                    job_tx,
                    expected_count,
                    &producer_wait_metric,
                );
                producer_metric.fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                if let Err(error) = result {
                    let _ = producer_errors.send(Err(error));
                }
            });

            for _ in 0..workers {
                let job_rx = job_rx.clone();
                let done_tx = done_tx.clone();
                let worker_nanos = worker_nanos.clone();
                scope.spawn(move || {
                    while let Ok((index, reads)) = job_rx.recv() {
                        let started = std::time::Instant::now();
                        let result = process_consensus_bucket(index, reorder.kmer_size, reads);
                        worker_nanos
                            .fetch_add(started.elapsed().as_nanos() as u64, Ordering::Relaxed);
                        if done_tx.send(result).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(job_rx);
            drop(done_tx);

            let writer_started = std::time::Instant::now();
            let result = write_stream_consensus_buckets(done_rx, permit_tx, payload)?;
            let metrics = ConsensusMetrics {
                producer: std::time::Duration::from_nanos(producer_nanos.load(Ordering::Relaxed)),
                producer_wait: std::time::Duration::from_nanos(
                    producer_wait_nanos.load(Ordering::Relaxed),
                ),
                workers: std::time::Duration::from_nanos(worker_nanos.load(Ordering::Relaxed)),
                writer: writer_started.elapsed(),
                writer_active: result.writer_active,
                writer_wait: result.writer_wait,
                multi_clusters: result.counts.0 as u64,
                singletons: result.counts.1 as u64,
                max_pending_buckets: result.max_pending_buckets,
                buckets: result.buckets,
                overlap_comparisons: result.overlap_comparisons,
                max_bucket_records: result.max_bucket_records,
                max_groups: result.max_groups,
                source_records: result.source_records,
            };
            Ok((result.counts, metrics))
        })?;
    let _ = (cluster_count, order_output, mate, temp_dir);
    Ok(metrics)
}

fn process_consensus_bucket(
    index: usize,
    kmer_size: usize,
    mut reads: Vec<AssignedSequence>,
) -> Result<ConsensusBucket> {
    let source_records = reads.len() as u64;
    if reads.len() == 1 {
        let read = reads.pop().unwrap();
        return Ok(ConsensusBucket {
            index,
            multi_cluster_count: 0,
            multi_bytes: Vec::new(),
            singletons: vec![(read.id, PackedSequence::try_from_bytes(&read.sequence)?)],
            order_entries: vec![(
                read.id,
                OrderEntry {
                    cluster: 0,
                    read: 0,
                    singleton: true,
                },
            )],
            #[cfg(test)]
            stats: StoreStats::default(),
            overlap_comparisons: 0,
            source_records,
            groups: 1,
            consensus_groups: vec![ConsensusCluster {
                reference: read.sequence.clone(),
                reads: vec![ConsensusRead {
                    original_id: read.id,
                    offset: 0,
                    reverse: false,
                    length: u32::try_from(read.sequence.len())
                        .map_err(|_| Error::InvalidInput("read length too large".into()))?,
                    mismatches: Vec::new(),
                }],
            }],
        });
    }
    let (clusters, overlap_comparisons) = cluster_sequences_profiled(kmer_size, reads)?;
    let groups = clusters.len() as u64;
    let (multi_bytes, multi_cluster_count, singletons, order_entries) =
        partition_clusters_v8(&clusters)?;
    Ok(ConsensusBucket {
        index,
        multi_cluster_count,
        multi_bytes,
        singletons,
        order_entries,
        #[cfg(test)]
        stats: StoreStats::default(),
        overlap_comparisons,
        source_records,
        groups,
        consensus_groups: clusters,
    })
}

fn produce_sorted_consensus_buckets<I, E>(
    sorted: I,
    limit: usize,
    permit_rx: crossbeam_channel::Receiver<()>,
    job_tx: crossbeam_channel::Sender<(usize, Vec<AssignedSequence>)>,
    expected_count: u64,
    wait_nanos: &AtomicU64,
) -> Result<()>
where
    I: Iterator<Item = std::result::Result<SortableRead, E>>,
    E: std::fmt::Display,
{
    let mut bucket = Vec::with_capacity(limit);
    let mut index = 0usize;
    let mut record_count = 0u64;
    let mut id_sum = 0u128;
    let mut id_square_sum = 0u128;
    let mut id_xor = 0u64;
    for read in sorted {
        let read = read.map_err(|error| {
            Error::InvalidInput(format!("cannot read sorted sequence spool: {error}"))
        })?;
        if expected_count != u64::MAX && read.original_id >= expected_count {
            return Err(Error::InvalidFormat(
                "sorted sequence spool contains an out-of-range ID".into(),
            ));
        }
        record_count += 1;
        id_sum += read.original_id as u128;
        id_square_sum += (read.original_id as u128) * (read.original_id as u128);
        id_xor ^= read.original_id;
        if bucket.len() == limit {
            let wait_started = std::time::Instant::now();
            permit_rx.recv().map_err(|_| {
                Error::InvalidInput("consensus workers exited before all buckets were sent".into())
            })?;
            wait_nanos.fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            let wait_started = std::time::Instant::now();
            job_tx
                .send((index, std::mem::take(&mut bucket)))
                .map_err(|_| {
                    Error::InvalidInput(
                        "consensus workers exited before all buckets were sent".into(),
                    )
                })?;
            wait_nanos.fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
            bucket = Vec::with_capacity(limit);
            index += 1;
        }
        bucket.push(AssignedSequence {
            id: read.original_id,
            position: read.position,
            reverse: read.reverse,
            anchored: read.hash != u64::MAX,
            sequence: read.sequence.decode(),
        });
    }
    if !bucket.is_empty() {
        let wait_started = std::time::Instant::now();
        permit_rx.recv().map_err(|_| {
            Error::InvalidInput("consensus workers exited before all buckets were sent".into())
        })?;
        wait_nanos.fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
        let wait_started = std::time::Instant::now();
        job_tx.send((index, bucket)).map_err(|_| {
            Error::InvalidInput("consensus workers exited before all buckets were sent".into())
        })?;
        wait_nanos.fetch_add(wait_started.elapsed().as_nanos() as u64, Ordering::Relaxed);
    }
    if expected_count != u64::MAX
        && !valid_id_permutation(expected_count, record_count, id_sum, id_square_sum, id_xor)
    {
        return Err(Error::InvalidFormat(
            "sorted sequence spool is not a complete ID permutation".into(),
        ));
    }
    Ok(())
}

fn valid_id_permutation(
    expected_count: u64,
    record_count: u64,
    id_sum: u128,
    id_square_sum: u128,
    id_xor: u64,
) -> bool {
    let n = expected_count as u128;
    let expected_sum = n * n.saturating_sub(1) / 2;
    let expected_square_sum = if n == 0 {
        0
    } else {
        n * (n - 1) * (2 * n - 1) / 6
    };
    record_count == expected_count
        && id_sum == expected_sum
        && id_square_sum == expected_square_sum
        && id_xor == xor_zero_to(expected_count.saturating_sub(1))
}

fn xor_zero_to(value: u64) -> u64 {
    match value % 4 {
        0 => value,
        1 => 1,
        2 => value + 1,
        _ => 0,
    }
}

#[cfg(test)]
fn produce_consensus_buckets(
    assignments: &Path,
    limit: usize,
    permit_rx: crossbeam_channel::Receiver<()>,
    job_tx: crossbeam_channel::Sender<(usize, Vec<AssignedRead>)>,
) -> Result<()> {
    let mut bucket = Vec::with_capacity(limit);
    let mut index = 0usize;
    for assignment in crate::reorder::AssignmentReader::open(assignments)? {
        let assignment = assignment?;
        if bucket.len() == limit {
            permit_rx.recv().map_err(|_| {
                Error::InvalidInput("consensus workers exited before all buckets were sent".into())
            })?;
            job_tx
                .send((index, std::mem::take(&mut bucket)))
                .map_err(|_| {
                    Error::InvalidInput(
                        "consensus workers exited before all buckets were sent".into(),
                    )
                })?;
            bucket = Vec::with_capacity(limit);
            index += 1;
        }
        bucket.push(AssignedRead {
            id: assignment.read_id,
            position: assignment.position,
            reverse: assignment.reverse,
            anchored: assignment.hash != u64::MAX,
        });
    }
    if !bucket.is_empty() {
        permit_rx.recv().map_err(|_| {
            Error::InvalidInput("consensus workers exited before all buckets were sent".into())
        })?;
        job_tx.send((index, bucket)).map_err(|_| {
            Error::InvalidInput("consensus workers exited before all buckets were sent".into())
        })?;
    }
    Ok(())
}

fn write_stream_consensus_buckets<W: Write>(
    done_rx: crossbeam_channel::Receiver<Result<ConsensusBucket>>,
    permit_tx: crossbeam_channel::Sender<()>,
    multi_output: &mut W,
) -> Result<ConsensusWriteResult> {
    let mut pending = BTreeMap::new();
    let mut next_to_write = 0usize;
    let mut multi_cluster_count = 0u32;
    let mut singleton_count = 0u32;
    let mut max_pending_buckets = 0u64;
    let mut writer_active = std::time::Duration::default();
    let mut writer_wait = std::time::Duration::default();
    let mut buckets = 0u64;
    let mut overlap_comparisons = 0u64;
    let mut max_bucket_records = 0u64;
    let mut max_groups = 0u64;
    let mut source_records = 0u64;
    while let Ok(result) = done_rx.recv() {
        let wait_started = std::time::Instant::now();
        let bucket = result?;
        writer_wait += wait_started.elapsed();
        let active_started = std::time::Instant::now();
        buckets += 1;
        overlap_comparisons += bucket.overlap_comparisons;
        max_bucket_records = max_bucket_records.max(bucket.source_records);
        max_groups = max_groups.max(bucket.groups);
        source_records += bucket.source_records;
        pending.insert(bucket.index, bucket);
        max_pending_buckets = max_pending_buckets.max(pending.len() as u64);
        while let Some(bucket) = pending.remove(&next_to_write) {
            for group in &bucket.consensus_groups {
                write_group_descriptor(multi_output, group)?;
                if group.reads.len() == 1 {
                    singleton_count += 1;
                } else {
                    multi_cluster_count += 1;
                }
            }
            let _ = permit_tx.send(());
            next_to_write += 1;
        }
        writer_active += active_started.elapsed();
    }
    if !pending.is_empty() {
        return Err(Error::InvalidInput(
            "consensus workers exited before all buckets were written".into(),
        ));
    }
    Ok(ConsensusWriteResult {
        counts: (multi_cluster_count, singleton_count),
        max_pending_buckets,
        #[cfg(test)]
        stats: StoreStats::default(),
        writer_active,
        writer_wait,
        buckets,
        overlap_comparisons,
        max_bucket_records,
        max_groups,
        source_records,
    })
}

#[cfg(test)]
fn write_consensus_buckets<W: Write, S: SingletonSequenceSink>(
    done_rx: crossbeam_channel::Receiver<Result<ConsensusBucket>>,
    permit_tx: crossbeam_channel::Sender<()>,
    multi_output: &mut W,
    singleton_sequences: &mut S,
    order_output: &mut std::fs::File,
) -> Result<ConsensusWriteResult> {
    let mut pending = BTreeMap::new();
    let mut next_to_write = 0usize;
    let mut multi_cluster_count = 0u32;
    let mut singleton_count = 0u32;
    let mut max_pending_buckets = 0u64;
    #[cfg(test)]
    let mut stats = StoreStats::default();
    let mut buckets = 0u64;
    let mut overlap_comparisons = 0u64;
    let mut max_bucket_records = 0u64;
    let mut max_groups = 0u64;
    let mut source_records = 0u64;
    let mut writer_active = std::time::Duration::default();
    let mut writer_wait = std::time::Duration::default();
    loop {
        let wait_started = std::time::Instant::now();
        let result = match done_rx.recv() {
            Ok(result) => result,
            Err(_) => break,
        };
        writer_wait += wait_started.elapsed();
        let active_started = std::time::Instant::now();
        let bucket = result?;
        buckets += 1;
        overlap_comparisons += bucket.overlap_comparisons;
        max_bucket_records = max_bucket_records.max(bucket.source_records);
        max_groups = max_groups.max(bucket.groups);
        source_records += bucket.source_records;
        pending.insert(bucket.index, bucket);
        max_pending_buckets = max_pending_buckets.max(pending.len() as u64);
        while let Some(bucket) = pending.remove(&next_to_write) {
            for (original_id, entry) in &bucket.order_entries {
                let offset = original_id
                    .checked_mul(ORDER_ENTRY_SIZE as u64)
                    .ok_or_else(|| Error::InvalidInput("order index too large".into()))?;
                order_output.seek(SeekFrom::Start(offset))?;
                order_output.write_all(&[entry.singleton as u8])?;
                let cluster = if entry.singleton {
                    singleton_count
                        .checked_add(entry.cluster)
                        .ok_or_else(|| Error::InvalidInput("singleton index overflow".into()))?
                } else {
                    multi_cluster_count
                        .checked_add(entry.cluster)
                        .ok_or_else(|| Error::InvalidInput("cluster index overflow".into()))?
                };
                order_output.write_all(&cluster.to_le_bytes())?;
                order_output.write_all(&entry.read.to_le_bytes())?;
            }
            multi_output.write_all(&bucket.multi_bytes)?;
            multi_cluster_count = multi_cluster_count
                .checked_add(bucket.multi_cluster_count)
                .ok_or_else(|| Error::InvalidInput("too many consensus clusters".into()))?;
            for (_, sequence) in bucket.singletons {
                singleton_sequences.push_sequence(sequence)?;
                singleton_count = singleton_count
                    .checked_add(1)
                    .ok_or_else(|| Error::InvalidInput("too many consensus clusters".into()))?;
            }
            #[cfg(test)]
            {
                stats.index_reads += bucket.stats.index_reads;
                stats.blob_reads += bucket.stats.blob_reads;
                stats.bytes_read += bucket.stats.bytes_read;
            }
            let _ = permit_tx.send(());
            next_to_write += 1;
        }
        writer_active += active_started.elapsed();
    }
    if !pending.is_empty() {
        return Err(Error::InvalidInput(
            "consensus workers exited before all buckets were written".into(),
        ));
    }
    Ok(ConsensusWriteResult {
        counts: (multi_cluster_count, singleton_count),
        max_pending_buckets,
        #[cfg(test)]
        stats,
        writer_active,
        writer_wait,
        buckets,
        overlap_comparisons,
        max_bucket_records,
        max_groups,
        source_records,
    })
}

struct ConsensusWriteResult {
    counts: (u32, u32),
    max_pending_buckets: u64,
    #[cfg(test)]
    stats: StoreStats,
    writer_active: std::time::Duration,
    writer_wait: std::time::Duration,
    buckets: u64,
    overlap_comparisons: u64,
    max_bucket_records: u64,
    max_groups: u64,
    source_records: u64,
}

#[derive(Clone)]
struct ClusterInput {
    id: u64,
    sequence: Vec<u8>,
    reverse: bool,
    offset: i32,
}

struct ClusterGroup {
    anchor: i32,
    minimum_offset: i32,
    reference: Vec<u8>,
    counts: Vec<[u32; 5]>,
    members: Vec<ClusterInput>,
}

#[cfg(test)]
fn cluster_bucket(
    store: &mut PairStore,
    mate: usize,
    kmer_size: usize,
    assignments: &[AssignedRead],
) -> Result<Vec<ConsensusCluster>> {
    let mut reads = Vec::with_capacity(assignments.len());
    for assignment in assignments {
        reads.push(AssignedSequence {
            id: assignment.id,
            position: assignment.position,
            reverse: assignment.reverse,
            anchored: assignment.anchored,
            sequence: store.read_sequence(assignment.id, mate)?,
        });
    }
    cluster_sequences_profiled(kmer_size, reads).map(|(clusters, _)| clusters)
}

fn cluster_sequences_profiled(
    kmer_size: usize,
    assignments: Vec<AssignedSequence>,
) -> Result<(Vec<ConsensusCluster>, u64)> {
    let mut groups: Vec<ClusterGroup> = Vec::new();
    let mut groups_by_anchor: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    let mut overlap_comparisons = 0u64;
    let max_read_length = assignments
        .iter()
        .map(|assignment| assignment.sequence.len())
        .max()
        .unwrap_or(0) as i32;
    for assignment in assignments {
        let mut sequence = assignment.sequence;
        if assignment.reverse {
            reverse_complement_in_place(&mut sequence);
        }
        let anchor = if assignment.anchored {
            let position = assignment.position as usize;
            if position
                .checked_add(kmer_size)
                .is_none_or(|end| end > sequence.len())
            {
                return Err(Error::InvalidInput(
                    "minimizer position exceeds read length".into(),
                ));
            }
            if assignment.reverse {
                (sequence.len() - kmer_size - position) as i32
            } else {
                position as i32
            }
        } else {
            0
        };
        let lower_anchor = anchor.saturating_sub(max_read_length);
        let upper_anchor = anchor.saturating_add(max_read_length);
        let mut best: Option<(usize, i32, usize, usize)> = None;
        let candidate_indices = groups_by_anchor
            .range(lower_anchor..=upper_anchor)
            .flat_map(|(_, indices)| indices.iter().copied())
            .collect::<Vec<_>>();
        for index in candidate_indices {
            let group = &groups[index];
            overlap_comparisons += 1;
            let offset = group.anchor - anchor;
            let consensus_offset = offset - group.minimum_offset;
            let (overlap, distance) =
                overlap_distance(&group.reference, &sequence, consensus_offset);
            let minimum_overlap = group.reference.len().min(sequence.len()).div_ceil(2);
            if overlap >= minimum_overlap && distance <= MAX_OVERLAP_MISMATCHES {
                let replace = best.is_none_or(|(_, _, best_overlap, best_distance)| {
                    distance < best_distance
                        || (distance == best_distance && overlap > best_overlap)
                });
                if replace {
                    best = Some((index, offset, overlap, distance));
                }
            }
        }
        match best {
            Some((index, offset, _, _)) => {
                let group = &mut groups[index];
                add_cluster_member(
                    group,
                    ClusterInput {
                        id: assignment.id,
                        sequence,
                        reverse: assignment.reverse,
                        offset,
                    },
                )?;
            }
            None => {
                let counts = sequence
                    .iter()
                    .map(|&base| {
                        let mut count = [0u32; 5];
                        count[base_index(base)] = 1;
                        count
                    })
                    .collect();
                let index = groups.len();
                groups.push(ClusterGroup {
                    anchor,
                    minimum_offset: 0,
                    reference: sequence.clone(),
                    counts,
                    members: vec![ClusterInput {
                        id: assignment.id,
                        sequence,
                        reverse: assignment.reverse,
                        offset: 0,
                    }],
                });
                groups_by_anchor.entry(anchor).or_default().push(index);
            }
        }
    }
    Ok((
        groups
            .into_iter()
            .map(finish_cluster)
            .collect::<Result<Vec<_>>>()?,
        overlap_comparisons,
    ))
}

fn add_cluster_member(group: &mut ClusterGroup, input: ClusterInput) -> Result<()> {
    let old_minimum = group.minimum_offset;
    let old_end = old_minimum as i64 + group.counts.len() as i64;
    let new_minimum = old_minimum.min(input.offset);
    let new_end = old_end.max(input.offset as i64 + input.sequence.len() as i64);
    let new_len = usize::try_from(new_end - new_minimum as i64)
        .map_err(|_| Error::InvalidInput("consensus reference too large".into()))?;
    if new_minimum < old_minimum {
        let prefix = usize::try_from(old_minimum - new_minimum)
            .map_err(|_| Error::InvalidInput("consensus offset too large".into()))?;
        let mut counts = vec![[0u32; 5]; prefix];
        counts.append(&mut group.counts);
        group.counts = counts;
        let mut reference = vec![0; prefix];
        reference.append(&mut group.reference);
        group.reference = reference;
    }
    group.counts.resize(new_len, [0u32; 5]);
    group.reference.resize(new_len, b'A');
    let start = usize::try_from(input.offset - new_minimum)
        .map_err(|_| Error::InvalidInput("negative consensus offset".into()))?;
    for (position, &base) in input.sequence.iter().enumerate() {
        let index = start + position;
        group.counts[index][base_index(base)] += 1;
        group.reference[index] = consensus_base(&group.counts[index]);
    }
    group.minimum_offset = new_minimum;
    group.members.push(input);
    Ok(())
}

fn overlap_distance(reference: &[u8], sequence: &[u8], offset: i32) -> (usize, usize) {
    let reference_start = offset.max(0) as usize;
    let sequence_start = (-offset).max(0) as usize;
    if reference_start >= reference.len() || sequence_start >= sequence.len() {
        return (0, 0);
    }
    let overlap = (reference.len() - reference_start).min(sequence.len() - sequence_start);
    let mut distance = 0usize;
    for (left, right) in reference[reference_start..reference_start + overlap]
        .iter()
        .zip(&sequence[sequence_start..sequence_start + overlap])
    {
        distance += (left != right) as usize;
        if distance > MAX_OVERLAP_MISMATCHES {
            break;
        }
    }
    (overlap, distance)
}

fn consensus_base(count: &[u32; 5]) -> u8 {
    let mut best = 0;
    for index in 1..5 {
        if count[index] > count[best] {
            best = index;
        }
    }
    [b'A', b'C', b'G', b'T', b'N'][best]
}

fn finish_cluster(mut group: ClusterGroup) -> Result<ConsensusCluster> {
    let mut reads = Vec::with_capacity(group.members.len());
    for input in &mut group.members {
        let offset = usize::try_from(input.offset - group.minimum_offset)
            .map_err(|_| Error::InvalidInput("negative consensus offset".into()))?;
        let mismatches = group.reference[offset..offset + input.sequence.len()]
            .iter()
            .zip(&input.sequence)
            .enumerate()
            .filter_map(|(position, (&expected, actual))| {
                (expected != *actual).then_some((position as u32, *actual))
            })
            .collect();
        reads.push(ConsensusRead {
            original_id: input.id,
            offset: u32::try_from(offset)
                .map_err(|_| Error::InvalidInput("consensus offset too large".into()))?,
            reverse: input.reverse,
            length: u32::try_from(input.sequence.len())
                .map_err(|_| Error::InvalidInput("read length too large".into()))?,
            mismatches,
        });
    }
    Ok(ConsensusCluster {
        reference: group.reference,
        reads,
    })
}
fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}
fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

fn reverse_complement_in_place(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        };
    }
}

type SingletonSequence = (u64, PackedSequence);

#[derive(Clone, Copy)]
struct OrderEntry {
    cluster: u32,
    read: u32,
    singleton: bool,
}

fn partition_clusters_v8(
    clusters: &[ConsensusCluster],
) -> Result<(Vec<u8>, u32, Vec<SingletonSequence>, Vec<(u64, OrderEntry)>)> {
    let mut multi_bytes = Vec::new();
    let mut multi_cluster_count = 0u32;
    let mut singletons = Vec::new();
    let mut order = Vec::new();
    let mut multi_index = 0u32;
    for cluster in clusters {
        if cluster.reads.len() == 1 {
            let read = &cluster.reads[0];
            if read.offset != 0
                || read.length as usize != cluster.reference.len()
                || !read.mismatches.is_empty()
            {
                return Err(Error::InvalidInput(
                    "singleton consensus invariant violated".into(),
                ));
            }
            let sequence = if read.reverse {
                reverse_complement(&cluster.reference)
            } else {
                cluster.reference.clone()
            };
            let singleton_index = singletons.len() as u32;
            singletons.push((read.original_id, PackedSequence::try_from_bytes(&sequence)?));
            order.push((
                read.original_id,
                OrderEntry {
                    cluster: singleton_index,
                    read: 0,
                    singleton: true,
                },
            ));
        } else {
            put_cluster_bodies_v8(&mut multi_bytes, std::slice::from_ref(cluster))?;
            for (read_index, read) in cluster.reads.iter().enumerate() {
                order.push((
                    read.original_id,
                    OrderEntry {
                        cluster: multi_index,
                        read: read_index as u32,
                        singleton: false,
                    },
                ));
            }
            multi_cluster_count = multi_cluster_count
                .checked_add(1)
                .ok_or_else(|| Error::InvalidInput("too many consensus clusters".into()))?;
            multi_index += 1;
        }
    }
    Ok((multi_bytes, multi_cluster_count, singletons, order))
}

fn put_cluster_bodies_v8<W: Write>(writer: &mut W, clusters: &[ConsensusCluster]) -> Result<()> {
    for cluster in clusters {
        if cluster.reads.len() < 2 {
            return Err(Error::InvalidInput("RS3 multi cluster is too small".into()));
        }
        put_reference_v8(writer, &cluster.reference)?;
        put_u32_writer(writer, cluster.reads.len() as u32)?;
        for read in &cluster.reads {
            writer.write_all(&[read.reverse as u8])?;
            put_varint(writer, read.offset as u64)?;
            put_varint(writer, read.length as u64)?;
            put_varint(writer, read.mismatches.len() as u64)?;
            let mut previous_position = 0u32;
            for &(position, base) in &read.mismatches {
                if position < previous_position {
                    return Err(Error::InvalidInput(
                        "consensus mismatches are not ordered".into(),
                    ));
                }
                put_varint(writer, (position - previous_position) as u64)?;
                writer.write_all(&[base])?;
                previous_position = position;
            }
        }
    }
    Ok(())
}

fn write_group_descriptor<W: Write>(writer: &mut W, group: &ConsensusCluster) -> Result<()> {
    put_varint(writer, group.reads.len() as u64)?;
    let singleton = group.reads.len() == 1;
    writer.write_all(&[singleton as u8])?;
    if singleton {
        let sequence = crate::reorder::decode_consensus_read(group, &group.reads[0])?;
        put_varint(writer, zigzag_delta(group.reads[0].original_id, 0))?;
        write_packed_sequence(writer, &PackedSequence::try_from_bytes(&sequence)?)?;
        return Ok(());
    }
    put_reference_v8(writer, &group.reference)?;
    let mut previous_id = 0;
    for read in &group.reads {
        put_varint(writer, zigzag_delta(read.original_id, previous_id))?;
        previous_id = read.original_id;
        writer.write_all(&[read.reverse as u8])?;
        put_varint(writer, read.offset as u64)?;
        put_varint(writer, read.length as u64)?;
        put_varint(writer, read.mismatches.len() as u64)?;
        let mut previous_position = 0u32;
        for &(position, base) in &read.mismatches {
            put_varint(writer, (position - previous_position) as u64)?;
            writer.write_all(&[base])?;
            previous_position = position;
        }
    }
    Ok(())
}

fn copy_group_descriptors<R: Read, W: Write>(
    input: &mut R,
    output: &mut W,
    count: u64,
    mate: usize,
    ids: &[Vec<u8>],
    qualities: &[Vec<u8>],
    options: &CompressionOptions,
    id_pattern: PairedIdPattern,
) -> Result<()> {
    const MAX_BATCH_RECORDS: usize = 4096;
    const MAX_BATCH_BYTES: usize = 4 * 1024 * 1024;
    let mut pending = None;
    let mut records = 0u64;
    while records < count {
        let group_count = read_varint(input)?;
        let kind = read_one(input)?;
        let reference = if kind == 0 {
            read_reference_stream(input)?
        } else if kind == 1 {
            Vec::new()
        } else {
            return Err(Error::InvalidFormat("invalid RS4 group kind".into()));
        };
        if kind == 1 && group_count != 1 {
            return Err(Error::InvalidFormat(
                "invalid singleton group record count".into(),
            ));
        }
        if kind == 0 {
            flush_singleton_batch(
                output,
                pending.take(),
                ids,
                qualities,
                options,
                mate,
                id_pattern,
            )?;
            put_varint(output, group_count)?;
            output.write_all(&[kind])?;
            put_reference_v8(output, &reference)?;
        }
        let mut group_ids = Vec::with_capacity(usize_from(group_count, "group")?);
        let mut previous_id = 0u64;
        for _ in 0..group_count {
            let original_id = zigzag_decode(read_varint(input)?, previous_id)?;
            previous_id = original_id;
            group_ids.push(original_id);
            if kind == 1 {
                let sequence = read_packed_sequence(input)?;
                let packed = sequence.is_packed();
                if !packed {
                    flush_singleton_batch(
                        output,
                        pending.take(),
                        ids,
                        qualities,
                        options,
                        mate,
                        id_pattern,
                    )?;
                    put_varint(output, 1)?;
                    output.write_all(&[1])?;
                    put_varint(output, zigzag_delta(original_id, 0))?;
                    write_packed_sequence(output, &sequence)?;
                    write_group_metadata(
                        output,
                        &[original_id],
                        ids,
                        qualities,
                        options,
                        mate,
                        id_pattern,
                    )?;
                    records += 1;
                    continue;
                }
                let fits = pending
                    .as_ref()
                    .is_some_and(|batch: &PendingSingletonBatch| {
                        batch.packed == packed
                            && batch.length == sequence.len()
                            && batch.ids.len() < MAX_BATCH_RECORDS
                            && batch.body_len.saturating_add(sequence.encoded_size())
                                <= MAX_BATCH_BYTES
                    });
                if !fits {
                    flush_singleton_batch(
                        output,
                        pending.take(),
                        ids,
                        qualities,
                        options,
                        mate,
                        id_pattern,
                    )?;
                    pending = Some(PendingSingletonBatch {
                        ids: Vec::with_capacity(MAX_BATCH_RECORDS),
                        bodies: Vec::with_capacity(MAX_BATCH_RECORDS),
                        length: sequence.len(),
                        packed,
                        body_len: 0,
                    });
                }
                let batch = pending.as_mut().expect("singleton batch initialized");
                batch.body_len += sequence.encoded_size();
                batch.ids.push(original_id);
                batch.bodies.push(sequence.encoded_bytes().to_vec());
            } else {
                put_varint(
                    output,
                    zigzag_delta(
                        original_id,
                        if group_ids.len() == 1 {
                            0
                        } else {
                            group_ids[group_ids.len() - 2]
                        },
                    ),
                )?;
                let reverse = read_one(input)?;
                output.write_all(&[reverse])?;
                let offset = read_varint(input)?;
                let length = read_varint(input)?;
                let mismatch_count = read_varint(input)?;
                put_varint(output, offset)?;
                put_varint(output, length)?;
                put_varint(output, mismatch_count)?;
                for _ in 0..mismatch_count {
                    let delta = read_varint(input)?;
                    let base = read_one(input)?;
                    put_varint(output, delta)?;
                    output.write_all(&[base])?;
                }
            }
        }
        if kind == 0 {
            write_group_metadata(
                output, &group_ids, ids, qualities, options, mate, id_pattern,
            )?;
            records += group_count;
        } else {
            records += 1;
        }
    }
    flush_singleton_batch(output, pending, ids, qualities, options, mate, id_pattern)?;
    Ok(())
}

struct PendingSingletonBatch {
    ids: Vec<u64>,
    bodies: Vec<Vec<u8>>,
    length: u32,
    packed: bool,
    body_len: usize,
}

fn flush_singleton_batch<W: Write>(
    output: &mut W,
    batch: Option<PendingSingletonBatch>,
    ids: &[Vec<u8>],
    qualities: &[Vec<u8>],
    options: &CompressionOptions,
    mate: usize,
    id_pattern: PairedIdPattern,
) -> Result<()> {
    let Some(batch) = batch else { return Ok(()) };
    if batch.ids.len() < 2 {
        put_varint(output, 1)?;
        output.write_all(&[1])?;
        put_varint(output, zigzag_delta(batch.ids[0], 0))?;
        write_packed_sequence(
            output,
            &PackedSequence::from_encoded(batch.length, batch.packed, batch.bodies[0].clone())?,
        )?;
        return write_group_metadata(
            output, &batch.ids, ids, qualities, options, mate, id_pattern,
        );
    }
    put_varint(output, batch.ids.len() as u64)?;
    output.write_all(&[2])?;
    put_varint(output, batch.length as u64)?;
    let mut previous_id = 0;
    for (&id, body) in batch.ids.iter().zip(&batch.bodies) {
        put_varint(output, zigzag_delta(id, previous_id))?;
        previous_id = id;
        output.write_all(body)?;
    }
    write_group_metadata(
        output, &batch.ids, ids, qualities, options, mate, id_pattern,
    )
}

fn flush_singleton_tail_batch<W: Write>(
    output: &mut W,
    batch: Option<PendingSingletonBatch>,
    ids: &[Vec<u8>],
    qualities: &[Vec<u8>],
    options: &CompressionOptions,
    mate: usize,
    id_pattern: PairedIdPattern,
) -> Result<()> {
    let Some(batch) = batch else { return Ok(()) };
    let mut entries = batch.ids.into_iter().zip(batch.bodies).collect::<Vec<_>>();
    entries.sort_unstable_by_key(|(id, _)| *id);
    let (ids_sorted, bodies_sorted): (Vec<_>, Vec<_>) = entries.into_iter().unzip();
    flush_singleton_batch(
        output,
        Some(PendingSingletonBatch {
            ids: ids_sorted,
            bodies: bodies_sorted,
            length: batch.length,
            packed: batch.packed,
            body_len: batch.body_len,
        }),
        ids,
        qualities,
        options,
        mate,
        id_pattern,
    )
}

fn put_reference_v8<W: Write>(writer: &mut W, reference: &[u8]) -> Result<()> {
    let clean = reference
        .iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'));
    if !clean {
        writer.write_all(&[0])?;
        put_bytes_writer(writer, reference, true)?;
        return Ok(());
    }
    let packed_width = reference.len().div_ceil(4);
    let mut packed = vec![0u8; packed_width];
    for (index, &base) in reference.iter().enumerate() {
        let value = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => unreachable!(),
        };
        packed[index / 4] |= value << (6 - 2 * (index % 4));
    }
    writer.write_all(&[1])?;
    put_u32_writer(writer, reference.len() as u32)?;
    put_u64_writer(writer, packed.len() as u64)?;
    writer.write_all(&packed)?;
    Ok(())
}

fn put_varint<W: Write>(writer: &mut W, mut value: u64) -> Result<()> {
    let mut bytes = [0u8; 10];
    let mut length = 0;
    loop {
        bytes[length] = (value as u8 & 0x7f) | if value >= 0x80 { 0x80 } else { 0 };
        length += 1;
        value >>= 7;
        if value == 0 {
            writer.write_all(&bytes[..length])?;
            return Ok(());
        }
    }
}

fn varint_len(value: u64) -> usize {
    if value == 0 {
        1
    } else {
        ((u64::BITS - value.leading_zeros()) as usize).div_ceil(7)
    }
}

fn zigzag_delta(value: u64, previous: u64) -> u64 {
    let delta = value as i128 - previous as i128;
    ((delta << 1) ^ (delta >> 127)) as u64
}

fn zigzag_decode(value: u64, previous: u64) -> Result<u64> {
    let delta = ((value >> 1) as i128) ^ -((value & 1) as i128);
    let decoded = previous as i128 + delta;
    u64::try_from(decoded)
        .map_err(|_| Error::InvalidFormat("RS4 original ID delta underflow".into()))
}

fn write_packed_sequence<W: Write>(writer: &mut W, sequence: &PackedSequence) -> Result<()> {
    writer.write_all(&[sequence.is_packed() as u8])?;
    put_varint(writer, sequence.len() as u64)?;
    put_varint(writer, sequence.encoded_size() as u64)?;
    writer.write_all(sequence.encoded_bytes())?;
    Ok(())
}

fn read_packed_sequence_stream<R: Read>(reader: &mut ChecksummedReader<'_, R>) -> Result<Vec<u8>> {
    let packed = reader.u8()? != 0;
    let length = u32::try_from(reader.varint()?)
        .map_err(|_| Error::InvalidFormat("sequence too large".into()))?;
    let encoded_length = usize_from(reader.varint()?, "sequence")?;
    let mut bytes = vec![0; encoded_length];
    reader.read_exact(&mut bytes)?;
    Ok(PackedSequence::from_encoded(length, packed, bytes)?.decode())
}

fn read_packed_sequence<R: Read>(reader: &mut R) -> Result<PackedSequence> {
    let mut mode = [0];
    reader.read_exact(&mut mode)?;
    let length = u32::try_from(read_varint(reader)?)
        .map_err(|_| Error::InvalidFormat("sequence too large".into()))?;
    let encoded_length = usize_from(read_varint(reader)?, "sequence")?;
    let mut bytes = vec![0; encoded_length];
    reader.read_exact(&mut bytes)?;
    PackedSequence::from_encoded(length, mode[0] != 0, bytes)
}

fn write_group_metadata<W: Write>(
    output: &mut W,
    group_ids: &[u64],
    ids: &[Vec<u8>],
    qualities: &[Vec<u8>],
    options: &CompressionOptions,
    mate: usize,
    id_pattern: PairedIdPattern,
) -> Result<()> {
    let store_id =
        options.preserve_ids && (mate == 0 || id_pattern == PairedIdPattern::Independent);
    if store_id {
        write_metadata_values(
            output,
            group_ids.iter().map(|&id| {
                ids.get(id as usize)
                    .map(Vec::as_slice)
                    .ok_or_else(|| Error::InvalidFormat("missing ID metadata".into()))
            }),
        )?;
    }
    if options.preserve_quality {
        write_quality_metadata(
            output,
            group_ids.iter().map(|&id| {
                qualities
                    .get(id as usize)
                    .map(Vec::as_slice)
                    .ok_or_else(|| Error::InvalidFormat("missing quality metadata".into()))
            }),
        )?;
    }
    Ok(())
}

fn write_metadata_values<'a, W: Write, I: Iterator<Item = Result<&'a [u8]>>>(
    output: &mut W,
    values: I,
) -> Result<()> {
    let values: Vec<_> = values.collect::<Result<Vec<_>>>()?;
    let fixed = values
        .first()
        .map(|value| values.iter().all(|other| other.len() == value.len()))
        .unwrap_or(true);
    let fixed_size = 1
        + varint_len(values.first().map_or(0, |value| value.len() as u64))
        + values.iter().map(|value| value.len()).sum::<usize>();
    let variable_size = 1
        + values
            .iter()
            .map(|value| varint_len(value.len() as u64))
            .sum::<usize>()
        + values.iter().map(|value| value.len()).sum::<usize>();

    let prefix_length = values
        .first()
        .map(|first| {
            values.iter().skip(1).fold(first.len(), |length, value| {
                length.min(
                    value
                        .iter()
                        .zip(first.iter())
                        .take(length)
                        .take_while(|(left, right)| left == right)
                        .count(),
                )
            })
        })
        .unwrap_or(0);
    let prefix_size = 1
        + varint_len(prefix_length as u64)
        + prefix_length
        + values
            .iter()
            .map(|value| {
                let suffix_length = value.len() - prefix_length;
                varint_len(suffix_length as u64) + suffix_length
            })
            .sum::<usize>();

    let fallback_size = if fixed { fixed_size } else { variable_size };
    if prefix_size < fallback_size {
        output.write_all(&[2])?;
        put_varint(output, prefix_length as u64)?;
        if let Some(first) = values.first() {
            output.write_all(&first[..prefix_length])?;
        }
        for value in &values {
            put_varint(output, (value.len() - prefix_length) as u64)?;
        }
        for value in values {
            output.write_all(&value[prefix_length..])?;
        }
    } else if fixed && fixed_size <= variable_size {
        output.write_all(&[1])?;
        put_varint(output, values.first().map_or(0, |value| value.len() as u64))?;
        for value in values {
            output.write_all(value)?;
        }
    } else {
        output.write_all(&[0])?;
        for value in &values {
            put_varint(output, value.len() as u64)?;
        }
        for value in values {
            output.write_all(value)?;
        }
    }
    Ok(())
}

const QUALITY_METADATA_PACKED: u8 = 1;
const QUALITY_METADATA_BLOCK_RECORDS: usize = 4096;

fn write_quality_metadata<'a, W: Write, I: Iterator<Item = Result<&'a [u8]>>>(
    output: &mut W,
    values: I,
) -> Result<()> {
    let values = values.collect::<Result<Vec<_>>>()?;
    let encoded =
        crate::quality_blocks::encode_blocks_bytes(&values, QUALITY_METADATA_BLOCK_RECORDS)?;
    output.write_all(&[QUALITY_METADATA_PACKED])?;
    put_varint(output, encoded.len() as u64)?;
    output.write_all(&encoded)?;
    Ok(())
}

fn read_group_metadata<R: Read>(
    reader: &mut ChecksummedReader<'_, R>,
    count: u64,
    packed: bool,
) -> Result<Vec<Vec<u8>>> {
    let count = usize_from(count, "group")?;
    if packed {
        if reader.u8()? != QUALITY_METADATA_PACKED {
            return Err(Error::InvalidFormat(
                "unsupported RS4 quality metadata encoding".into(),
            ));
        }
        let encoded_length = usize_from(reader.varint()?, "quality metadata")?;
        let mut encoded = vec![0; encoded_length];
        reader.read_exact(&mut encoded)?;
        let values = crate::quality_blocks::QualityBlocks::decode(&encoded)?;
        if values.len() != count {
            return Err(Error::InvalidFormat(
                "RS4 quality block record count mismatch".into(),
            ));
        }
        return Ok(values);
    }
    let mode = reader.u8()?;
    let mut lengths = match mode {
        0 => (0..count)
            .map(|_| {
                reader
                    .varint()
                    .and_then(|value| usize_from(value, "metadata"))
            })
            .collect::<Result<Vec<_>>>()?,
        1 => vec![usize_from(reader.varint()?, "metadata")?; count],
        2 => {
            let prefix_length = usize_from(reader.varint()?, "metadata")?;
            let mut prefix = vec![0; prefix_length];
            reader.read_exact(&mut prefix)?;
            let suffix_lengths = (0..count)
                .map(|_| {
                    reader
                        .varint()
                        .and_then(|value| usize_from(value, "metadata"))
                })
                .collect::<Result<Vec<_>>>()?;
            let mut values = Vec::with_capacity(count);
            for suffix_length in suffix_lengths {
                let mut value = prefix.clone();
                let start = value.len();
                let end = start
                    .checked_add(suffix_length)
                    .ok_or_else(|| Error::InvalidFormat("RS4 metadata value too large".into()))?;
                value.resize(end, 0);
                reader.read_exact(&mut value[start..])?;
                values.push(value);
            }
            return Ok(values);
        }
        _ => {
            return Err(Error::InvalidFormat(
                "unsupported RS4 metadata encoding".into(),
            ))
        }
    };
    let mut values = Vec::with_capacity(count);
    for length in lengths.drain(..) {
        let mut value = vec![0; length];
        reader.read_exact(&mut value)?;
        values.push(value);
    }
    Ok(values)
}

fn append_compact_order<W: Write>(
    payload: &mut W,
    order_path: &Path,
    records: u64,
    temp_dir: &Path,
) -> Result<()> {
    let compact_path = temp_dir.join(format!(
        "{}.compact",
        order_path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| Error::InvalidInput("invalid order path".into()))?
    ));
    let mut input = std::fs::File::open(order_path)?;
    let mut compact = std::fs::File::create(&compact_path)?;
    let mut slot = [0u8; ORDER_ENTRY_SIZE];
    for _ in 0..records {
        input.read_exact(&mut slot)?;
        let singleton = slot[0] != 0;
        let cluster = u32::from_le_bytes(slot[1..5].try_into().unwrap());
        let read = u32::from_le_bytes(slot[5..9].try_into().unwrap());
        put_varint(&mut compact, ((cluster as u64) << 1) | singleton as u64)?;
        if !singleton {
            put_varint(&mut compact, read as u64)?;
        }
    }
    compact.flush()?;
    let mut compact = std::fs::File::open(&compact_path)?;
    put_u64_writer(payload, compact.metadata()?.len())?;
    std::io::copy(&mut compact, payload)?;
    std::fs::remove_file(compact_path)?;
    Ok(())
}

#[cfg(test)]
struct RecordColumns {
    ids: StoredColumn,
    qualities: StoredColumn,
}

struct PreparedColumnBuilder {
    ids: PreparedIdWriter,
    qualities: PreparedQualityWriter,
}

struct PreparedColumn {
    ids: PreparedIdColumn,
    qualities: PreparedQuality,
}

struct PreparedIdWriter {
    path: std::path::PathBuf,
    builder: ExternalChunkBuilder<Vec<u8>>,
    count: u64,
    common_length: Option<u32>,
    variable_length: bool,
    data_len: u64,
}

struct PreparedIdColumn {
    path: std::path::PathBuf,
    count: u64,
    common_length: Option<u32>,
    variable_length: bool,
    data_len: u64,
}

impl PreparedIdWriter {
    fn create(temp_dir: &Path, mate: usize) -> Result<Self> {
        let path = temp_dir.join(format!("records.{mate}.ids.lz4"));
        let builder = ExternalChunkBuilder::new(std::fs::File::create(&path)?, 1)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        Ok(Self {
            path,
            builder,
            count: 0,
            common_length: None,
            variable_length: false,
            data_len: 0,
        })
    }

    fn push(&mut self, id: &[u8]) -> Result<()> {
        let length =
            u32::try_from(id.len()).map_err(|_| Error::InvalidInput("ID exceeds 4 GiB".into()))?;
        self.variable_length |= self.common_length.is_some_and(|common| common != length);
        self.common_length.get_or_insert(length);
        self.data_len = self
            .data_len
            .checked_add(length as u64)
            .ok_or_else(|| Error::InvalidInput("ID column too large".into()))?;
        self.count += 1;
        self.builder
            .add(id.to_vec())
            .map_err(|error| Error::InvalidInput(error.to_string()))
    }

    fn finish(self) -> Result<PreparedIdColumn> {
        let chunk = self
            .builder
            .finish()
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        drop(chunk);
        Ok(PreparedIdColumn {
            path: self.path,
            count: self.count,
            common_length: self.common_length,
            variable_length: self.variable_length,
            data_len: self.data_len,
        })
    }
}

impl PreparedIdColumn {
    fn write_to<W: Write>(self, payload: &mut W, enabled: bool, count: usize) -> Result<()> {
        if !enabled {
            std::fs::remove_file(self.path)?;
            payload.write_all(&[0])?;
            return Ok(());
        }
        if self.count != count as u64 {
            return Err(Error::InvalidInput("ID column count mismatch".into()));
        }
        let path = self.path;
        if self.variable_length {
            payload.write_all(&[2])?;
            put_u64_writer(payload, self.data_len)?;
            for_each_id(&path, self.count, |id| {
                put_u32_writer(payload, id.len() as u32)
            })?;
            for_each_id(&path, self.count, |id| {
                payload.write_all(&id)?;
                Ok(())
            })?;
        } else {
            let width = self.common_length.unwrap_or(0);
            payload.write_all(&[1])?;
            put_u32_writer(payload, width)?;
            put_u64_writer(payload, self.data_len)?;
            for_each_id(&path, self.count, |id| {
                payload.write_all(&id)?;
                Ok(())
            })?;
        }
        std::fs::remove_file(path)?;
        Ok(())
    }
}

fn for_each_id<F>(path: &Path, expected_count: u64, mut visit: F) -> Result<()>
where
    F: FnMut(Vec<u8>) -> Result<()>,
{
    let chunk = ExternalChunk::<Vec<u8>>::open(std::fs::File::open(path)?)
        .map_err(|error| Error::InvalidFormat(error.to_string()))?;
    let mut count = 0u64;
    for id in chunk {
        let id = id.map_err(|error| Error::InvalidFormat(error.to_string()))?;
        visit(id)?;
        count += 1;
    }
    if count != expected_count {
        return Err(Error::InvalidFormat("ID chunk count mismatch".into()));
    }
    Ok(())
}

fn load_prepared_metadata(column: &PreparedColumn) -> Result<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
    let mut ids = Vec::new();
    if column.ids.count != 0 {
        for_each_id(&column.ids.path, column.ids.count, |id| {
            ids.push(id);
            Ok(())
        })?;
    }
    let qualities = if column.qualities.block_count == 0 {
        Vec::new()
    } else {
        let mut bytes = Vec::new();
        copy_external_bytes(&column.qualities.path, &mut bytes)?;
        let mut encoded = Vec::with_capacity(8 + bytes.len());
        encoded.extend_from_slice(&(PreparedQualityWriter::BLOCK_RECORDS as u32).to_le_bytes());
        encoded.extend_from_slice(&column.qualities.block_count.to_le_bytes());
        encoded.extend_from_slice(&bytes);
        let values = crate::quality_blocks::QualityBlocks::decode(&encoded)?;
        values
    };
    Ok((ids, qualities))
}

struct PreparedQualityWriter {
    output: Lz4ByteWriter,
    values: Vec<Vec<u8>>,
    width: Option<usize>,
    block_count: u32,
    encoded_len: u64,
}

struct PreparedQuality {
    path: std::path::PathBuf,
    block_count: u32,
    encoded_len: u64,
}

impl PreparedColumnBuilder {
    fn create(temp_dir: &Path, mate: usize) -> Result<Self> {
        Ok(Self {
            ids: PreparedIdWriter::create(temp_dir, mate)?,
            qualities: PreparedQualityWriter::create(temp_dir, mate)?,
        })
    }

    fn finish(self) -> Result<PreparedColumn> {
        Ok(PreparedColumn {
            ids: self.ids.finish()?,
            qualities: self.qualities.finish()?,
        })
    }
}

impl PreparedQualityWriter {
    const BLOCK_RECORDS: usize = 65_536;

    fn create(temp_dir: &Path, mate: usize) -> Result<Self> {
        Ok(Self {
            output: Lz4ByteWriter::create(
                temp_dir.join(format!("records.{mate}.quality-blocks.lz4")),
                1,
            )?,
            values: Vec::with_capacity(Self::BLOCK_RECORDS),
            width: None,
            block_count: 0,
            encoded_len: 0,
        })
    }

    fn push(&mut self, quality: &[u8]) -> Result<()> {
        if !self.values.is_empty()
            && (self.values.len() == Self::BLOCK_RECORDS
                || self.width.is_some_and(|width| width != quality.len()))
        {
            self.flush_block()?;
        }
        self.values.push(quality.to_vec());
        self.width = Some(quality.len());
        Ok(())
    }

    fn flush_block(&mut self) -> Result<()> {
        if self.values.is_empty() {
            return Ok(());
        }
        let bytes = crate::quality_blocks::encode_block_bytes(&self.values)?;
        self.encoded_len = self
            .encoded_len
            .checked_add(bytes.len() as u64)
            .ok_or_else(|| Error::InvalidInput("quality blocks too large".into()))?;
        self.output.write_all(&bytes)?;
        self.block_count = self
            .block_count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many quality blocks".into()))?;
        self.values.clear();
        self.width = None;
        Ok(())
    }

    fn finish(mut self) -> Result<PreparedQuality> {
        self.flush_block()?;
        Ok(PreparedQuality {
            path: self.output.finish()?,
            block_count: self.block_count,
            encoded_len: self.encoded_len,
        })
    }
}

#[cfg(test)]
impl RecordColumns {
    fn create(temp_dir: &Path, mate: usize) -> Result<Self> {
        Ok(Self {
            ids: StoredColumn::create(temp_dir, mate, "ids")?,
            qualities: StoredColumn::create(temp_dir, mate, "qualities")?,
        })
    }
}

#[allow(dead_code)]
pub(crate) struct PreparedRecordsBuilder {
    columns: [PreparedColumnBuilder; 2],
    spools: Option<[SequenceSpoolWriter; 2]>,
    id_pattern: Option<PairedIdPattern>,
    max_sortable_bytes: [usize; 2],
    count: u64,
}

pub(crate) struct PreparedPair {
    ids: [Vec<u8>; 2],
    qualities: [Vec<u8>; 2],
    reads: [SortableRead; 2],
}

pub(crate) fn prepare_pair(
    first: ReadRecord,
    second: ReadRecord,
    original_id: u64,
    options: &CompressionOptions,
) -> Result<PreparedPair> {
    let records = [first, second];
    let mut reads = Vec::with_capacity(2);
    for record in &records {
        let assignment = crate::reorder::primary_minimizer(
            &record.sequence,
            options.reorder.kmer_size,
            options.reorder.window_size,
        );
        let (hash, position, reverse) = assignment.unwrap_or((u64::MAX, 0, false));
        reads.push(SortableRead {
            hash,
            original_id,
            position,
            reverse,
            sequence: PackedSequence::try_from_bytes(&record.sequence)?,
        });
    }
    let [first, second] = records;
    let [first_read, second_read] = reads.try_into().unwrap();
    Ok(PreparedPair {
        ids: [first.id, second.id],
        qualities: [
            first.quality.unwrap_or_default(),
            second.quality.unwrap_or_default(),
        ],
        reads: [first_read, second_read],
    })
}

#[allow(dead_code)]
pub(crate) struct PreparedRecords {
    columns: [PreparedColumn; 2],
    spool_paths: [std::path::PathBuf; 2],
    max_sortable_bytes: [usize; 2],
    id_pattern: PairedIdPattern,
    count: u64,
}

pub(crate) struct PreparedMetadata {
    columns: [PreparedColumn; 2],
    id_pattern: PairedIdPattern,
    count: u64,
    ids: [Vec<Vec<u8>>; 2],
    qualities: [Vec<Vec<u8>>; 2],
}

impl PreparedRecordsBuilder {
    #[allow(dead_code)]
    pub(crate) fn create(temp_dir: &Path) -> Result<Self> {
        std::fs::create_dir_all(temp_dir)?;
        Ok(Self {
            columns: [
                PreparedColumnBuilder::create(temp_dir, 0)?,
                PreparedColumnBuilder::create(temp_dir, 1)?,
            ],
            spools: Some([
                SequenceSpoolWriter::create(temp_dir.join("sequences.0.spool"))?,
                SequenceSpoolWriter::create(temp_dir.join("sequences.1.spool"))?,
            ]),
            id_pattern: None,
            max_sortable_bytes: [0; 2],
            count: 0,
        })
    }

    #[cfg(test)]
    pub(crate) fn push(
        &mut self,
        first: &ReadRecord,
        second: &ReadRecord,
        options: &CompressionOptions,
    ) -> Result<()> {
        if self.count == u32::MAX as u64 {
            return Err(Error::InvalidInput(
                "RS3 supports at most u32::MAX paired records".into(),
            ));
        }
        if first.sequence.len() > i32::MAX as usize || second.sequence.len() > i32::MAX as usize {
            return Err(Error::InvalidInput(
                "short-read consensus supports sequences up to i32::MAX bases".into(),
            ));
        }
        let prepared = prepare_pair(first.clone(), second.clone(), self.count, options)?;
        self.push_prepared(prepared, options)
    }

    #[allow(dead_code)]
    pub(crate) fn push_prepared(
        &mut self,
        prepared: PreparedPair,
        options: &CompressionOptions,
    ) -> Result<()> {
        let PreparedPair {
            ids,
            qualities,
            reads,
        } = prepared;
        if options.preserve_ids {
            let current = PairedIdPattern::detect(&ids[0], &ids[1]);
            self.id_pattern = Some(match self.id_pattern {
                None => current,
                Some(previous) if previous == current => previous,
                Some(_) => PairedIdPattern::Independent,
            });
            self.columns[0].ids.push(&ids[0])?;
            self.columns[1].ids.push(&ids[1])?;
        }
        if options.preserve_quality {
            self.columns[0].qualities.push(&qualities[0])?;
            self.columns[1].qualities.push(&qualities[1])?;
        }
        if let Some(spools) = &mut self.spools {
            for (mate, (spool, read)) in spools.iter_mut().zip(reads).enumerate() {
                self.max_sortable_bytes[mate] = self.max_sortable_bytes[mate]
                    .max(std::mem::size_of::<SortableRead>() + read.sequence.encoded_size());
                spool.push(&read)?;
            }
        }
        self.count = self
            .count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many paired records".into()))?;
        Ok(())
    }

    #[allow(dead_code)]
    pub(crate) fn finish(self) -> Result<PreparedRecords> {
        let PreparedRecordsBuilder {
            columns,
            spools,
            id_pattern,
            max_sortable_bytes,
            count,
        } = self;
        let [first, second] = spools
            .ok_or_else(|| Error::InvalidInput("metadata builder has no sequence spools".into()))?;
        let [first_columns, second_columns] = columns;
        Ok(PreparedRecords {
            columns: [first_columns.finish()?, second_columns.finish()?],
            spool_paths: [first.finish()?, second.finish()?],
            max_sortable_bytes,
            id_pattern: id_pattern.unwrap_or(PairedIdPattern::Independent),
            count,
        })
    }

    pub(crate) fn create_metadata(temp_dir: &Path) -> Result<Self> {
        std::fs::create_dir_all(temp_dir)?;
        Ok(Self {
            columns: [
                PreparedColumnBuilder::create(temp_dir, 0)?,
                PreparedColumnBuilder::create(temp_dir, 1)?,
            ],
            spools: None,
            id_pattern: None,
            max_sortable_bytes: [0; 2],
            count: 0,
        })
    }

    pub(crate) fn push_metadata(
        &mut self,
        prepared: PreparedPair,
        options: &CompressionOptions,
    ) -> Result<[SortableRead; 2]> {
        let PreparedPair {
            ids,
            qualities,
            reads,
        } = prepared;
        if options.preserve_ids {
            let current = PairedIdPattern::detect(&ids[0], &ids[1]);
            self.id_pattern = Some(match self.id_pattern {
                None => current,
                Some(previous) if previous == current => previous,
                Some(_) => PairedIdPattern::Independent,
            });
            self.columns[0].ids.push(&ids[0])?;
            self.columns[1].ids.push(&ids[1])?;
        }
        if options.preserve_quality {
            self.columns[0].qualities.push(&qualities[0])?;
            self.columns[1].qualities.push(&qualities[1])?;
        }
        self.count = self
            .count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many paired records".into()))?;
        Ok(reads)
    }

    pub(crate) fn finish_metadata(self) -> Result<PreparedMetadata> {
        let [first, second] = self.columns;
        let columns = [first.finish()?, second.finish()?];
        let (first_ids, first_qualities) = load_prepared_metadata(&columns[0])?;
        let (second_ids, second_qualities) = load_prepared_metadata(&columns[1])?;
        Ok(PreparedMetadata {
            columns,
            id_pattern: self.id_pattern.unwrap_or(PairedIdPattern::Independent),
            count: self.count,
            ids: [first_ids, second_ids],
            qualities: [first_qualities, second_qualities],
        })
    }
}

#[cfg(test)]
fn build_store_record_columns(
    store: &mut PairStore,
    options: &CompressionOptions,
    temp_dir: &Path,
) -> Result<([RecordColumns; 2], PairedIdPattern)> {
    let mut columns = [
        RecordColumns::create(temp_dir, 0)?,
        RecordColumns::create(temp_dir, 1)?,
    ];
    let mut id_pattern = None;
    for pair_id in 0..store.count() {
        let (first, second) =
            store.read_pair_metadata(pair_id, options.preserve_ids, options.preserve_quality)?;
        if options.preserve_ids {
            let current = PairedIdPattern::detect(&first.id, &second.id);
            id_pattern = Some(match id_pattern {
                None => current,
                Some(previous) if previous == current => previous,
                Some(_) => PairedIdPattern::Independent,
            });
            columns[0].ids.push(&first.id)?;
            columns[1].ids.push(&second.id)?;
        }
        if options.preserve_quality {
            columns[0].qualities.push(&first.quality)?;
            columns[1].qualities.push(&second.quality)?;
        }
    }
    Ok((columns, id_pattern.unwrap_or(PairedIdPattern::Independent)))
}

#[cfg(test)]
fn write_store_record_columns<W: Write + Seek>(
    payload: &mut W,
    columns: RecordColumns,
    count: u64,
    ids: bool,
    quality: bool,
    quality_blocks_records: usize,
) -> Result<()> {
    let count =
        u32::try_from(count).map_err(|_| Error::InvalidInput("too many records for RS3".into()))?;
    put_u32_writer(payload, count)?;
    let RecordColumns {
        ids: id_column,
        qualities,
    } = columns;
    id_column.write_to(payload, ids, count as usize)?;
    let quality_length_position = payload.stream_position()?;
    put_u64_writer(payload, 0)?;
    let quality_start = payload.stream_position()?;
    let mut blocks =
        crate::quality_blocks::QualityBlockWriter::new(&mut *payload, quality_blocks_records)?;
    if quality {
        qualities.write_quality_blocks(&mut blocks, count as usize, quality_blocks_records)?;
    }
    blocks.finish()?;
    let quality_end = payload.stream_position()?;
    let quality_len = quality_end
        .checked_sub(quality_start)
        .ok_or_else(|| Error::InvalidInput("invalid quality section length".into()))?;
    payload.seek(SeekFrom::Start(quality_length_position))?;
    put_u64_writer(payload, quality_len)?;
    payload.seek(SeekFrom::Start(quality_end))?;
    Ok(())
}

fn write_prepared_record_columns<W: Write>(
    payload: &mut W,
    columns: PreparedColumn,
    count: u64,
    ids: bool,
    quality: bool,
) -> Result<()> {
    let count =
        u32::try_from(count).map_err(|_| Error::InvalidInput("too many records for RS3".into()))?;
    put_u32_writer(payload, count)?;
    columns.ids.write_to(payload, ids, count as usize)?;
    let quality_len = 8u64
        .checked_add(if quality {
            columns.qualities.encoded_len
        } else {
            0
        })
        .ok_or_else(|| Error::InvalidInput("quality section too large".into()))?;
    put_u64_writer(payload, quality_len)?;
    payload.write_all(&(PreparedQualityWriter::BLOCK_RECORDS as u32).to_le_bytes())?;
    payload.write_all(
        &if quality {
            columns.qualities.block_count
        } else {
            0
        }
        .to_le_bytes(),
    )?;
    if quality {
        copy_external_bytes(&columns.qualities.path, payload)?;
    }
    std::fs::remove_file(columns.qualities.path)?;
    Ok(())
}

#[cfg(test)]
struct StoredColumn {
    data: std::io::BufWriter<std::fs::File>,
    lengths: std::io::BufWriter<std::fs::File>,
    common_length: Option<u32>,
    variable_length: bool,
    data_len: u64,
    clean_dna: bool,
    clean_dna_records: u64,
}

#[cfg(test)]
impl StoredColumn {
    fn create(temp_dir: &Path, mate: usize, name: &str) -> Result<Self> {
        Ok(Self {
            data: std::io::BufWriter::with_capacity(
                8 * 1024 * 1024,
                create_temporary_column(&temp_dir.join(format!("records.{mate}.{name}.data")))?,
            ),
            lengths: std::io::BufWriter::with_capacity(
                8 * 1024 * 1024,
                create_temporary_column(&temp_dir.join(format!("records.{mate}.{name}.lengths")))?,
            ),
            common_length: None,
            variable_length: false,
            data_len: 0,
            clean_dna: true,
            clean_dna_records: 0,
        })
    }

    fn push(&mut self, bytes: &[u8]) -> Result<()> {
        let length = u32::try_from(bytes.len())
            .map_err(|_| Error::InvalidInput("record field exceeds 4 GiB".into()))?;
        self.lengths.write_all(&length.to_le_bytes())?;
        self.data.write_all(bytes)?;
        let clean_dna = bytes
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'));
        self.clean_dna &= clean_dna;
        self.clean_dna_records += clean_dna as u64;
        self.data_len = self
            .data_len
            .checked_add(length as u64)
            .ok_or_else(|| Error::InvalidInput("record column too large".into()))?;
        match self.common_length {
            None => self.common_length = Some(length),
            Some(common) if common != length => self.variable_length = true,
            Some(_) => {}
        }
        Ok(())
    }

    fn write_to<W: Write>(mut self, payload: &mut W, enabled: bool, count: usize) -> Result<()> {
        if !enabled {
            payload.write_all(&[0])?;
            return Ok(());
        }
        self.data.flush()?;
        self.data.get_mut().seek(SeekFrom::Start(0))?;
        if self.variable_length {
            payload.write_all(&[2])?;
            put_u64_writer(payload, self.data_len)?;
            self.lengths.flush()?;
            self.lengths.get_mut().seek(SeekFrom::Start(0))?;
            std::io::copy(self.lengths.get_mut(), payload)?;
        } else {
            payload.write_all(&[1])?;
            let width = self.common_length.unwrap_or(0);
            put_u32_writer(payload, width)?;
            let expected = (width as u64)
                .checked_mul(count as u64)
                .ok_or_else(|| Error::InvalidInput("record column too large".into()))?;
            if expected != self.data_len {
                return Err(Error::InvalidInput("record column length mismatch".into()));
            }
            put_u64_writer(payload, self.data_len)?;
        }
        std::io::copy(self.data.get_mut(), payload)?;
        Ok(())
    }

    #[cfg(test)]
    fn write_quality_blocks<W: Write + Seek>(
        mut self,
        blocks: &mut crate::quality_blocks::QualityBlockWriter<W>,
        count: usize,
        block_records: usize,
    ) -> Result<()> {
        let block_records = block_records.max(1);
        self.lengths.flush()?;
        self.data.flush()?;
        self.lengths.get_mut().seek(SeekFrom::Start(0))?;
        self.data.get_mut().seek(SeekFrom::Start(0))?;
        if self.offsets_len() != count {
            return Err(Error::InvalidInput("record column length mismatch".into()));
        }
        let mut values = Vec::with_capacity(block_records.min(count));
        let mut width = None;
        for _ in 0..count {
            let mut length = [0u8; 4];
            self.lengths.get_mut().read_exact(&mut length)?;
            let length = u32::from_le_bytes(length) as usize;
            if !values.is_empty()
                && (values.len() == block_records || width.is_some_and(|width| width != length))
            {
                blocks.push(&values)?;
                values.clear();
            }
            let mut value = vec![0; length];
            self.data.get_mut().read_exact(&mut value)?;
            values.push(value);
            width = Some(length);
        }
        if !values.is_empty() {
            blocks.push(&values)?;
        }
        Ok(())
    }

    #[cfg(test)]
    fn offsets_len(&self) -> usize {
        let len = self
            .lengths
            .get_ref()
            .metadata()
            .map(|m| m.len() as usize)
            .unwrap_or(0);
        len / std::mem::size_of::<u32>()
    }

    #[cfg(test)]
    fn write_packed_dna_to<W: Write>(mut self, payload: &mut W, count: usize) -> Result<()> {
        if self.variable_length {
            return self.write_to(payload, true, count);
        }
        let width = self.common_length.unwrap_or(0) as usize;
        let expected = width
            .checked_mul(count)
            .ok_or_else(|| Error::InvalidInput("DNA column too large".into()))?;
        if expected as u64 != self.data_len {
            return Err(Error::InvalidInput("DNA column length mismatch".into()));
        }
        let packed_width = width.div_ceil(4);
        let packed_len = packed_width
            .checked_mul(count)
            .ok_or_else(|| Error::InvalidInput("packed DNA column too large".into()))?;
        if self.clean_dna {
            payload.write_all(&[4])?;
            put_u32_writer(payload, width as u32)?;
            put_u64_writer(payload, packed_len as u64)?;
        } else {
            let clean = usize::try_from(self.clean_dna_records)
                .map_err(|_| Error::InvalidInput("too many DNA records".into()))?;
            let ambiguous = count
                .checked_sub(clean)
                .ok_or_else(|| Error::InvalidInput("invalid DNA record count".into()))?;
            let encoded_len = count
                .checked_add(
                    clean
                        .checked_mul(packed_width)
                        .ok_or_else(|| Error::InvalidInput("mixed DNA column too large".into()))?,
                )
                .and_then(|length| length.checked_add(ambiguous.checked_mul(width)?))
                .ok_or_else(|| Error::InvalidInput("mixed DNA column too large".into()))?;
            payload.write_all(&[5])?;
            put_u32_writer(payload, width as u32)?;
            put_u64_writer(payload, encoded_len as u64)?;
        }
        self.data.flush()?;
        self.data.get_mut().seek(SeekFrom::Start(0))?;
        let mut sequence = vec![0; width];
        let mut packed = vec![0; packed_width];
        for _ in 0..count {
            self.data.get_mut().read_exact(&mut sequence)?;
            let clean = sequence
                .iter()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'));
            if !self.clean_dna {
                payload.write_all(&[(!clean) as u8])?;
                if !clean {
                    payload.write_all(&sequence)?;
                    continue;
                }
            }
            packed.fill(0);
            for (index, &base) in sequence.iter().enumerate() {
                let value = match base {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => unreachable!(),
                };
                packed[index / 4] |= value << (6 - 2 * (index % 4));
            }
            payload.write_all(&packed)?;
        }
        Ok(())
    }
}

#[cfg(test)]
fn create_temporary_column(path: &Path) -> Result<std::fs::File> {
    Ok(std::fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(path)?)
}

fn checksum_file<R: Read>(mut reader: R) -> Result<u32> {
    let mut hash = 2166136261u32;
    let mut buffer = [0u8; 1024 * 1024];
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        for &byte in &buffer[..read] {
            hash = (hash ^ byte as u32).wrapping_mul(16777619);
        }
    }
    Ok(hash)
}

struct PayloadChecksumWriter<W> {
    inner: W,
    length: u64,
    checksum: u32,
}

impl<W> PayloadChecksumWriter<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            length: 0,
            checksum: 2166136261,
        }
    }

    fn finish(self) -> (W, u64, u32) {
        (self.inner, self.length, self.checksum)
    }
}

impl<W: Write> Write for PayloadChecksumWriter<W> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(bytes)?;
        self.length = self
            .length
            .checked_add(written as u64)
            .ok_or_else(|| std::io::Error::other("uncompressed payload length overflow"))?;
        for &byte in &bytes[..written] {
            self.checksum = (self.checksum ^ byte as u32).wrapping_mul(16777619);
        }
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

fn flags(options: &CompressionOptions) -> u16 {
    (options.preserve_ids as u16) | ((options.preserve_quality as u16) << 1)
}
fn put_u32_writer<W: Write>(writer: &mut W, value: u32) -> Result<()> {
    writer.write_all(&value.to_le_bytes()).map_err(Into::into)
}
fn put_u64_writer<W: Write>(writer: &mut W, value: u64) -> Result<()> {
    writer.write_all(&value.to_le_bytes()).map_err(Into::into)
}

fn put_bytes_writer<W: Write>(writer: &mut W, bytes: &[u8], enabled: bool) -> Result<()> {
    let length = if enabled { bytes.len() } else { 0 };
    put_u32_writer(
        writer,
        u32::try_from(length).map_err(|_| Error::InvalidInput("field too large".into()))?,
    )?;
    if enabled {
        writer.write_all(bytes)?;
    }
    Ok(())
}
fn put_u16_writer<W: Write>(writer: &mut W, value: u16) -> Result<()> {
    writer.write_all(&value.to_le_bytes()).map_err(Into::into)
}
fn usize_from(value: u64, field: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| Error::InvalidFormat(format!("{field} too large")))
}
fn read_u16<R: Read>(reader: &mut R) -> Result<u16> {
    let mut bytes = [0; 2];
    reader.read_exact(&mut bytes)?;
    Ok(u16::from_le_bytes(bytes))
}
fn read_u32<R: Read>(reader: &mut R) -> Result<u32> {
    let mut bytes = [0; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}
fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut bytes = [0; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

struct ChecksummedReader<'a, R> {
    inner: &'a mut R,
    remaining: u64,
    checksum: u32,
}

impl<'a, R: Read> ChecksummedReader<'a, R> {
    fn new(inner: &'a mut R, remaining: u64) -> Self {
        Self {
            inner,
            remaining,
            checksum: 2166136261,
        }
    }
    fn read_exact(&mut self, bytes: &mut [u8]) -> Result<()> {
        if bytes.len() as u64 > self.remaining {
            return Err(Error::InvalidFormat("RS4 payload length mismatch".into()));
        }
        self.inner.read_exact(bytes)?;
        self.remaining -= bytes.len() as u64;
        for &byte in bytes.iter() {
            self.checksum = (self.checksum ^ byte as u32).wrapping_mul(16777619);
        }
        Ok(())
    }
    fn u8(&mut self) -> Result<u8> {
        let mut b = [0];
        self.read_exact(&mut b)?;
        Ok(b[0])
    }
    fn u16(&mut self) -> Result<u16> {
        let mut b = [0; 2];
        self.read_exact(&mut b)?;
        Ok(u16::from_le_bytes(b))
    }
    fn u32(&mut self) -> Result<u32> {
        let mut b = [0; 4];
        self.read_exact(&mut b)?;
        Ok(u32::from_le_bytes(b))
    }
    fn u64(&mut self) -> Result<u64> {
        let mut b = [0; 8];
        self.read_exact(&mut b)?;
        Ok(u64::from_le_bytes(b))
    }
    fn varint(&mut self) -> Result<u64> {
        let mut value = 0u64;
        for shift in (0..=63).step_by(7) {
            let byte = self.u8()?;
            if shift == 63 && byte > 1 {
                return Err(Error::InvalidFormat("RS4 varint overflow".into()));
            }
            value |= ((byte & 0x7f) as u64) << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
        }
        Err(Error::InvalidFormat("unterminated RS4 varint".into()))
    }
    fn bytes_u64(&mut self) -> Result<Vec<u8>> {
        let len = usize_from(self.u64()?, "field")?;
        let mut value = vec![0; len];
        self.read_exact(&mut value)?;
        Ok(value)
    }
    fn finish(mut self, expected: u32) -> Result<()> {
        let mut buffer = [0u8; 64 * 1024];
        while self.remaining != 0 {
            let length = (self.remaining as usize).min(buffer.len());
            self.read_exact(&mut buffer[..length])?;
        }
        if self.checksum != expected {
            return Err(Error::InvalidFormat("RS4 checksum failure".into()));
        }
        Ok(())
    }
}

impl<R: Read> Read for ChecksummedReader<'_, R> {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        if self.remaining == 0 {
            return Ok(0);
        }
        let length = bytes.len().min(self.remaining as usize);
        let read = self.inner.read(&mut bytes[..length])?;
        self.remaining -= read as u64;
        for &byte in &bytes[..read] {
            self.checksum = (self.checksum ^ byte as u32).wrapping_mul(16777619);
        }
        Ok(read)
    }
}

fn read_one<R: Read>(reader: &mut R) -> Result<u8> {
    let mut b = [0];
    reader.read_exact(&mut b)?;
    Ok(b[0])
}
fn read_bytes_u64<R: Read>(reader: &mut R) -> Result<Vec<u8>> {
    let len = usize_from(read_u64(reader)?, "field")?;
    let mut value = vec![0; len];
    reader.read_exact(&mut value)?;
    Ok(value)
}
fn read_varint<R: Read>(reader: &mut R) -> Result<u64> {
    let mut value = 0u64;
    for shift in (0..=63).step_by(7) {
        let byte = read_one(reader)?;
        if shift == 63 && byte > 1 {
            return Err(Error::InvalidFormat("RS4 varint overflow".into()));
        }
        value |= ((byte & 0x7f) as u64) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    Err(Error::InvalidFormat("unterminated RS4 varint".into()))
}
fn read_reference_stream<R: Read>(reader: &mut R) -> Result<Vec<u8>> {
    match read_one(reader)? {
        0 => {
            let len = usize_from(read_u32(reader)? as u64, "reference")?;
            let mut value = vec![0; len];
            reader.read_exact(&mut value)?;
            Ok(value)
        }
        1 => {
            let width = usize_from(read_u32(reader)? as u64, "reference")?;
            let data_len = usize_from(read_u64(reader)?, "reference")?;
            let mut packed = vec![0; data_len];
            reader.read_exact(&mut packed)?;
            if data_len != width.div_ceil(4) {
                return Err(Error::InvalidFormat("invalid packed RS4 reference".into()));
            }
            Ok((0..width)
                .map(|i| {
                    [b'A', b'C', b'G', b'T'][((packed[i / 4] >> (6 - 2 * (i % 4))) & 3) as usize]
                })
                .collect())
        }
        _ => Err(Error::InvalidFormat("invalid RS4 reference mode".into())),
    }
}
fn read_consensus_read_stream<R: Read>(
    reader: &mut ChecksummedReader<'_, R>,
    original_id: u64,
) -> Result<ConsensusRead> {
    let reverse = reader.u8()? != 0;
    let offset = u32::try_from(read_varint(reader)?)
        .map_err(|_| Error::InvalidFormat("read offset too large".into()))?;
    let length = u32::try_from(read_varint(reader)?)
        .map_err(|_| Error::InvalidFormat("read length too large".into()))?;
    let mismatch_count = usize_from(read_varint(reader)?, "mismatch count")?;
    let mut mismatches = Vec::with_capacity(mismatch_count);
    let mut position: u32 = 0;
    for _ in 0..mismatch_count {
        position = position
            .checked_add(
                u32::try_from(read_varint(reader)?)
                    .map_err(|_| Error::InvalidFormat("mismatch position too large".into()))?,
            )
            .ok_or_else(|| Error::InvalidFormat("mismatch position overflow".into()))?;
        mismatches.push((position, reader.u8()?));
    }
    Ok(ConsensusRead {
        original_id,
        offset,
        reverse,
        length,
        mismatches,
    })
}

struct Cursor<'a> {
    data: &'a [u8],
    position: usize,
}
impl<'a> Cursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, position: 0 }
    }
    fn take(&mut self, length: usize) -> Result<&'a [u8]> {
        if self.data.len().saturating_sub(self.position) < length {
            return Err(Error::InvalidFormat("truncated RS3 payload".into()));
        }
        let value = &self.data[self.position..self.position + length];
        self.position += length;
        Ok(value)
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
    fn bytes(&mut self) -> Result<Vec<u8>> {
        let length = usize_from(self.u32()? as u64, "field")?;
        Ok(self.take(length)?.to_vec())
    }
    fn multi_clusters_v8(&mut self) -> Result<Vec<ConsensusCluster>> {
        let count = self.u32()? as usize;
        let mut clusters = Vec::with_capacity(count);
        for _ in 0..count {
            let reference = self.reference_v8()?;
            let reads_count = self.u32()? as usize;
            if reads_count < 2 {
                return Err(Error::InvalidFormat("invalid RS3 multi cluster".into()));
            }
            let mut reads = Vec::with_capacity(reads_count);
            for _ in 0..reads_count {
                let reverse = self.u8()? != 0;
                let offset = u32::try_from(self.varint()?)
                    .map_err(|_| Error::InvalidFormat("RS3 offset too large".into()))?;
                let length = u32::try_from(self.varint()?)
                    .map_err(|_| Error::InvalidFormat("RS3 read length too large".into()))?;
                let mismatch_count = usize_from(self.varint()?, "mismatch count")?;
                let mut mismatches = Vec::with_capacity(mismatch_count);
                let mut position = 0u32;
                for _ in 0..mismatch_count {
                    let delta = u32::try_from(self.varint()?).map_err(|_| {
                        Error::InvalidFormat("RS3 mismatch position too large".into())
                    })?;
                    position = position.checked_add(delta).ok_or_else(|| {
                        Error::InvalidFormat("RS3 mismatch position overflow".into())
                    })?;
                    mismatches.push((position, self.u8()?));
                }
                let end = (offset as usize)
                    .checked_add(length as usize)
                    .ok_or_else(|| Error::InvalidFormat("RS3 read range overflow".into()))?;
                if end > reference.len() {
                    return Err(Error::InvalidFormat("RS3 read exceeds reference".into()));
                }
                reads.push(ConsensusRead {
                    original_id: 0,
                    offset,
                    reverse,
                    length,
                    mismatches,
                });
            }
            clusters.push(ConsensusCluster { reference, reads });
        }
        Ok(clusters)
    }

    fn clusters_v8(&mut self) -> Result<Vec<ConsensusCluster>> {
        let mut clusters = self.multi_clusters_v8()?;
        let singleton_count = self.u32()? as usize;
        let sequences = self.column(singleton_count, true, "singleton sequence")?;
        for reference in sequences {
            let length = u32::try_from(reference.len())
                .map_err(|_| Error::InvalidFormat("RS3 reference too large".into()))?;
            clusters.push(ConsensusCluster {
                reference,
                reads: vec![ConsensusRead {
                    original_id: 0,
                    offset: 0,
                    reverse: false,
                    length,
                    mismatches: Vec::new(),
                }],
            });
        }
        Ok(clusters)
    }

    fn order_v8(&mut self, records: usize) -> Result<Vec<OrderEntry>> {
        let byte_len = usize_from(self.u64()?, "RS3 order index")?;
        let bytes = self.take(byte_len)?;
        let mut entries = Vec::with_capacity(records);
        let mut order = Cursor::new(bytes);
        for _ in 0..records {
            let code = order.varint()?;
            let singleton = code & 1 != 0;
            let cluster = u32::try_from(code >> 1)
                .map_err(|_| Error::InvalidFormat("RS3 order cluster too large".into()))?;
            let read = if singleton {
                0
            } else {
                u32::try_from(order.varint()?)
                    .map_err(|_| Error::InvalidFormat("RS3 order read too large".into()))?
            };
            entries.push(OrderEntry {
                singleton,
                cluster,
                read,
            });
        }
        if order.position != bytes.len() {
            return Err(Error::InvalidFormat("trailing RS3 order index data".into()));
        }
        Ok(entries)
    }

    fn reference_v8(&mut self) -> Result<Vec<u8>> {
        match self.u8()? {
            0 => self.bytes(),
            1 => {
                let width = self.u32()? as usize;
                let packed_width = width.div_ceil(4);
                let data_len = usize_from(self.u64()?, "RS3 reference")?;
                if data_len != packed_width {
                    return Err(Error::InvalidFormat("invalid RS3 packed reference".into()));
                }
                let packed = self.take(data_len)?;
                let bases = [b'A', b'C', b'G', b'T'];
                let mut reference = Vec::with_capacity(width);
                for position in 0..width {
                    let value = (packed[position / 4] >> (6 - 2 * (position % 4))) & 3;
                    reference.push(bases[value as usize]);
                }
                Ok(reference)
            }
            _ => Err(Error::InvalidFormat("invalid RS3 reference mode".into())),
        }
    }

    fn varint(&mut self) -> Result<u64> {
        let mut value = 0u64;
        for shift in (0..=63).step_by(7) {
            let byte = self.u8()?;
            if shift == 63 && byte > 1 {
                return Err(Error::InvalidFormat("RS3 varint overflow".into()));
            }
            value |= ((byte & 0x7f) as u64) << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
        }
        Err(Error::InvalidFormat("unterminated RS3 varint".into()))
    }
    fn record_metadata_v8(&mut self, ids: bool, quality: bool) -> Result<Vec<ReadRecord>> {
        let count = self.u32()? as usize;
        let id_values = self.column(count, ids, "ID")?;
        let quality_len = usize_from(self.u64()?, "quality blocks")?;
        let quality_bytes = self.take(quality_len)?;
        let quality_values = crate::quality_blocks::QualityBlocks::decode(quality_bytes)?;
        if quality_values.len() != if quality { count } else { 0 } {
            return Err(Error::InvalidFormat(
                "RS3 quality block record count mismatch".into(),
            ));
        }
        let quality_values = if quality {
            quality_values
        } else {
            (0..count).map(|_| Vec::new()).collect()
        };
        Ok(id_values
            .into_iter()
            .zip(quality_values)
            .map(|(id, quality_value)| ReadRecord {
                id,
                sequence: Vec::new(),
                quality: quality.then_some(quality_value),
            })
            .collect())
    }

    fn columnar_records(&mut self, ids: bool, quality: bool) -> Result<Vec<ReadRecord>> {
        let count = self.u32()? as usize;
        let id_values = self.column(count, ids, "ID")?;
        let quality_values = self.column(count, quality, "quality")?;
        Ok(id_values
            .into_iter()
            .zip(quality_values)
            .map(|(id, quality_value)| ReadRecord {
                id,
                sequence: Vec::new(),
                quality: quality.then_some(quality_value),
            })
            .collect())
    }

    fn column(&mut self, count: usize, expected: bool, field: &str) -> Result<Vec<Vec<u8>>> {
        let mode = self.u8()?;
        if mode == 0 {
            if expected {
                return Err(Error::InvalidFormat(format!("missing RS3 {field} column")));
            }
            return Ok((0..count).map(|_| Vec::new()).collect());
        }
        if !expected {
            return Err(Error::InvalidFormat(format!(
                "unexpected RS3 {field} column"
            )));
        }
        match mode {
            1 => {
                let width = self.u32()? as usize;
                let data_len = usize_from(self.u64()?, field)?;
                let expected_len = width
                    .checked_mul(count)
                    .ok_or_else(|| Error::InvalidFormat(format!("RS3 {field} column too large")))?;
                if data_len != expected_len {
                    return Err(Error::InvalidFormat(format!(
                        "invalid RS3 {field} column length"
                    )));
                }
                let data = self.take(data_len)?;
                if width == 0 {
                    return Ok((0..count).map(|_| Vec::new()).collect());
                }
                Ok(data.chunks_exact(width).map(<[u8]>::to_vec).collect())
            }
            2 => {
                let data_len = usize_from(self.u64()?, field)?;
                let lengths = (0..count)
                    .map(|_| self.u32().map(|length| length as usize))
                    .collect::<Result<Vec<_>>>()?;
                let expected_len = lengths.iter().try_fold(0usize, |total, &length| {
                    total.checked_add(length).ok_or_else(|| {
                        Error::InvalidFormat(format!("RS3 {field} column too large"))
                    })
                })?;
                if data_len != expected_len {
                    return Err(Error::InvalidFormat(format!(
                        "invalid RS3 {field} column length"
                    )));
                }
                let data = self.take(data_len)?;
                let mut position = 0usize;
                Ok(lengths
                    .into_iter()
                    .map(|length| {
                        let value = data[position..position + length].to_vec();
                        position += length;
                        value
                    })
                    .collect())
            }
            3 => {
                let width = self.u32()? as usize;
                let block_records = self.u32()? as usize;
                if block_records == 0 {
                    return Err(Error::InvalidFormat(format!(
                        "invalid RS3 {field} transpose block"
                    )));
                }
                let data_len = usize_from(self.u64()?, field)?;
                let expected_len = width
                    .checked_mul(count)
                    .ok_or_else(|| Error::InvalidFormat(format!("RS3 {field} column too large")))?;
                if data_len != expected_len {
                    return Err(Error::InvalidFormat(format!(
                        "invalid RS3 {field} column length"
                    )));
                }
                let data = self.take(data_len)?;
                let mut values = (0..count).map(|_| vec![0; width]).collect::<Vec<_>>();
                let mut source = 0usize;
                for start in (0..count).step_by(block_records) {
                    let rows = (count - start).min(block_records);
                    for column in 0..width {
                        for row in 0..rows {
                            values[start + row][column] = data[source + column * rows + row];
                        }
                    }
                    source += rows * width;
                }
                Ok(values)
            }
            4 => {
                let width = self.u32()? as usize;
                let packed_width = width.div_ceil(4);
                let data_len = usize_from(self.u64()?, field)?;
                let expected_len = packed_width
                    .checked_mul(count)
                    .ok_or_else(|| Error::InvalidFormat(format!("RS3 {field} column too large")))?;
                if data_len != expected_len {
                    return Err(Error::InvalidFormat(format!(
                        "invalid RS3 {field} column length"
                    )));
                }
                let data = self.take(data_len)?;
                let bases = [b'A', b'C', b'G', b'T'];
                let mut values = Vec::with_capacity(count);
                for packed in data.chunks_exact(packed_width.max(1)).take(count) {
                    let mut sequence = Vec::with_capacity(width);
                    for position in 0..width {
                        let value = (packed[position / 4] >> (6 - 2 * (position % 4))) & 3;
                        sequence.push(bases[value as usize]);
                    }
                    values.push(sequence);
                }
                if width == 0 {
                    values.resize_with(count, Vec::new);
                }
                Ok(values)
            }
            5 => {
                let width = self.u32()? as usize;
                let packed_width = width.div_ceil(4);
                let data_len = usize_from(self.u64()?, field)?;
                let data = self.take(data_len)?;
                let mut source = Cursor::new(data);
                let bases = [b'A', b'C', b'G', b'T'];
                let mut values = Vec::with_capacity(count);
                for _ in 0..count {
                    match source.u8()? {
                        0 => {
                            let packed = source.take(packed_width)?;
                            let mut sequence = Vec::with_capacity(width);
                            for position in 0..width {
                                let value = (packed[position / 4] >> (6 - 2 * (position % 4))) & 3;
                                sequence.push(bases[value as usize]);
                            }
                            values.push(sequence);
                        }
                        1 => values.push(source.take(width)?.to_vec()),
                        _ => {
                            return Err(Error::InvalidFormat(format!(
                                "invalid RS3 {field} DNA tag"
                            )))
                        }
                    }
                }
                if source.position != data.len() {
                    return Err(Error::InvalidFormat(format!(
                        "trailing RS3 {field} DNA data"
                    )));
                }
                Ok(values)
            }
            _ => Err(Error::InvalidFormat(format!(
                "invalid RS3 {field} column mode"
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlap_distance_stops_after_rejection_threshold() {
        assert_eq!(overlap_distance(b"AAAAAAAA", b"AAAACCCC", 0), (8, 4));
        assert_eq!(overlap_distance(b"AAAAAAAA", b"CCCCCCCC", 0), (8, 5));
    }

    #[test]
    fn validates_complete_sorted_id_permutations() {
        assert!(valid_id_permutation(4, 4, 6, 14, 0));
        assert!(!valid_id_permutation(4, 3, 3, 5, 3));
        assert!(!valid_id_permutation(4, 4, 5, 9, 1));
    }

    #[test]
    fn prepared_pipeline_matches_store_backed_archive_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let old_dir = directory.path().join("old");
        let new_dir = directory.path().join("new");
        std::fs::create_dir_all(&old_dir).unwrap();
        std::fs::create_dir_all(&new_dir).unwrap();
        let records = [
            (
                ReadRecord::new("r0/1", b"ACGTACG".to_vec(), Some(b"DDD9DDD".to_vec())),
                ReadRecord::new("r0/2", b"TGCATGC".to_vec(), Some(b"DDDD9DD".to_vec())),
            ),
            (
                ReadRecord::new("r1/1", b"ACGTACA".to_vec(), Some(b"DD,9DDD".to_vec())),
                ReadRecord::new("r1/2", b"TGCATGA".to_vec(), Some(b"DD!D9DD".to_vec())),
            ),
            (
                ReadRecord::new("r2/1", b"ACNTACG".to_vec(), Some(b"DDD9DDD".to_vec())),
                ReadRecord::new("r2/2", b"TGCNTGC".to_vec(), Some(b"DDDD9DD".to_vec())),
            ),
        ];
        let options = CompressionOptions {
            threads: 1,
            zstd_level: 3,
            reorder: crate::ReorderOptions {
                kmer_size: 3,
                window_size: 3,
                max_bucket_records: 2,
                sort_chunk_size: 1,
                sort_threads: 1,
                temp_dir: None,
            },
            ..CompressionOptions::default()
        };

        let mut store = PairStore::create(old_dir.join("store")).unwrap();
        let mut prepared = PreparedRecordsBuilder::create(&new_dir).unwrap();
        for (first, second) in &records {
            store.append(first, second).unwrap();
            prepared.push(first, second, &options).unwrap();
        }
        let mut old_archive = Vec::new();
        write_reordered_store_profiled(&mut store, &mut old_archive, &options, &old_dir).unwrap();
        let mut new_archive = Vec::new();
        write_prepared_profiled(
            prepared.finish().unwrap(),
            &mut new_archive,
            &options,
            &new_dir,
        )
        .unwrap();

        let old_decoded = crate::decompress(
            std::io::Cursor::new(old_archive),
            &crate::DecompressionOptions::default(),
        )
        .unwrap();
        let new_decoded = crate::decompress(
            std::io::Cursor::new(new_archive),
            &crate::DecompressionOptions::default(),
        )
        .unwrap();
        assert_eq!(new_decoded, old_decoded);
    }

    #[test]
    fn singleton_descriptor_batches_clean_sequences_and_falls_back() {
        let sequences = [
            PackedSequence::from_bytes(b"ACGT"),
            PackedSequence::from_bytes(b"TGCA"),
            PackedSequence::from_bytes(b"ACG"),
            PackedSequence::from_bytes(b"ACGN"),
        ];
        let mut input = Vec::new();
        for (id, sequence) in sequences.iter().enumerate() {
            put_varint(&mut input, 1).unwrap();
            input.write_all(&[1]).unwrap();
            put_varint(&mut input, zigzag_delta(id as u64, 0)).unwrap();
            write_packed_sequence(&mut input, sequence).unwrap();
        }
        let options = CompressionOptions {
            preserve_ids: false,
            preserve_quality: false,
            ..CompressionOptions::default()
        };
        let mut output = Vec::new();
        copy_group_descriptors(
            &mut std::io::Cursor::new(input),
            &mut output,
            sequences.len() as u64,
            0,
            &[],
            &[],
            &options,
            PairedIdPattern::Independent,
        )
        .unwrap();

        let mut output = std::io::Cursor::new(output);
        assert_eq!(read_varint(&mut output).unwrap(), 2);
        assert_eq!(read_one(&mut output).unwrap(), 2);
        assert_eq!(read_varint(&mut output).unwrap(), 4);
        assert_eq!(
            zigzag_decode(read_varint(&mut output).unwrap(), 0).unwrap(),
            0
        );
        assert_eq!(read_one(&mut output).unwrap(), 0x1b);
        assert_eq!(
            zigzag_decode(read_varint(&mut output).unwrap(), 0).unwrap(),
            1
        );
        assert_eq!(read_one(&mut output).unwrap(), 0xe4);
        assert_eq!(read_varint(&mut output).unwrap(), 1);
        assert_eq!(read_one(&mut output).unwrap(), 1);
        assert_eq!(read_varint(&mut output).unwrap(), 4);
        assert_eq!(read_packed_sequence(&mut output).unwrap(), sequences[2]);
        assert_eq!(read_varint(&mut output).unwrap(), 1);
        assert_eq!(read_one(&mut output).unwrap(), 1);
        assert_eq!(read_varint(&mut output).unwrap(), 6);
        assert_eq!(read_packed_sequence(&mut output).unwrap(), sequences[3]);
    }

    #[test]
    fn similar_store_reads_share_consensus_and_split_at_limit() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        for sequence in [
            b"ACGTACGT".as_slice(),
            b"ACGTACGA",
            b"ACGTACGG",
            b"TTTTTTTT",
        ] {
            let record = ReadRecord::new("id", sequence.to_vec(), None);
            store.append(&record, &record).unwrap();
        }
        let assignments = [
            AssignedRead {
                id: 0,
                position: 0,
                reverse: false,
                anchored: false,
            },
            AssignedRead {
                id: 1,
                position: 0,
                reverse: false,
                anchored: false,
            },
            AssignedRead {
                id: 2,
                position: 0,
                reverse: false,
                anchored: false,
            },
        ];
        let clusters = cluster_bucket(&mut store, 0, 3, &assignments).unwrap();
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].reads.len(), 3);
        for read in &clusters[0].reads {
            let decoded = crate::reorder::decode_consensus_read(&clusters[0], read).unwrap();
            assert_eq!(
                decoded,
                store.read_first(read.original_id).unwrap().sequence
            );
        }
    }

    #[test]
    fn overlap_clustering_builds_positional_consensus() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        let template = (0..200).map(|index| b"ACGT"[index % 4]).collect::<Vec<_>>();
        let first = ReadRecord::new("first", template[..150].to_vec(), None);
        let second = ReadRecord::new("second", template[50..].to_vec(), None);
        store.append(&first, &first).unwrap();
        store.append(&second, &second).unwrap();
        let assignments = [
            AssignedRead {
                id: 0,
                position: 100,
                reverse: false,
                anchored: true,
            },
            AssignedRead {
                id: 1,
                position: 50,
                reverse: false,
                anchored: true,
            },
        ];
        let clusters = cluster_bucket(&mut store, 0, 3, &assignments).unwrap();
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].reference.len(), 200);
        assert_eq!(clusters[0].reads[0].offset, 0);
        assert_eq!(clusters[0].reads[1].offset, 50);
        for read in &clusters[0].reads {
            assert_eq!(
                crate::reorder::decode_consensus_read(&clusters[0], read).unwrap(),
                if read.original_id == 0 {
                    first.sequence.clone()
                } else {
                    second.sequence.clone()
                }
            );
        }
    }

    #[test]
    fn parallel_store_consensus_matches_serial_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        for index in 0..48 {
            let mut first = b"ACGTACGTACGTACGT".to_vec();
            let length = first.len();
            first[index % length] = b"ACGT"[(index / 3) % 4];
            let mut second = reverse_complement(&first);
            second[(index * 5) % length] = b"ACGT"[(index / 7) % 4];
            store
                .append(
                    &ReadRecord::new(format!("r{index}/1"), first, None),
                    &ReadRecord::new(format!("r{index}/2"), second, None),
                )
                .unwrap();
        }
        let options = crate::ReorderOptions {
            kmer_size: 3,
            window_size: 3,
            sort_chunk_size: 5,
            sort_threads: 2,
            max_bucket_records: 3,
            ..crate::ReorderOptions::default()
        };
        for mate in 0..2 {
            let assignments = crate::reorder::write_store_assignments(
                &mut store,
                &options,
                directory.path(),
                mate,
            )
            .unwrap();
            let mut serial = Vec::new();
            let serial_stats = write_store_clusters(
                &mut serial,
                &store,
                &assignments,
                mate,
                &options,
                1,
                directory.path(),
            )
            .unwrap();
            let mut parallel = Vec::new();
            let parallel_stats = write_store_clusters(
                &mut parallel,
                &store,
                &assignments,
                mate,
                &options,
                8,
                directory.path(),
            )
            .unwrap();
            assert_eq!(serial, parallel);
            assert_eq!(serial_stats, parallel_stats);
        }
    }

    #[test]
    fn record_columns_support_fixed_and_variable_widths() {
        let directory = tempfile::tempdir().unwrap();
        for (name, values) in [
            ("fixed", vec![b"AAAA".as_slice(), b"CCCC", b"GGGG"]),
            ("variable", vec![b"A".as_slice(), b"CCCC", b""]),
        ] {
            let mut column = StoredColumn::create(directory.path(), 0, name).unwrap();
            for value in &values {
                column.push(value).unwrap();
            }
            let mut encoded = Vec::new();
            column.write_to(&mut encoded, true, values.len()).unwrap();
            let mut cursor = Cursor::new(&encoded);
            assert_eq!(cursor.column(values.len(), true, name).unwrap(), values);
        }
    }

    fn metadata_round_trip(values: &[&[u8]], expected_mode: u8) -> Vec<Vec<u8>> {
        let mut encoded = Vec::new();
        write_metadata_values(&mut encoded, values.iter().copied().map(Ok)).unwrap();
        assert_eq!(encoded[0], expected_mode);
        let length = encoded.len() as u64;
        let mut input = std::io::Cursor::new(encoded);
        let mut reader = ChecksummedReader::new(&mut input, length);
        read_group_metadata(&mut reader, values.len() as u64, false).unwrap()
    }

    #[test]
    fn group_metadata_compresses_common_prefixes() {
        let values = [b"sample-000/1".as_slice(), b"sample-001/1", b"sample-002/1"];
        assert_eq!(
            metadata_round_trip(&values, 2),
            values
                .iter()
                .map(|value| value.to_vec())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn group_metadata_falls_back_when_prefix_is_not_useful() {
        let values = [b"a".as_slice(), b"bbbb", b"ccc"];
        assert_eq!(
            metadata_round_trip(&values, 0),
            values
                .iter()
                .map(|value| value.to_vec())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn group_metadata_preserves_arbitrary_bytes() {
        let values = [b"\x80\x00prefix-a".as_slice(), b"\x80\x00prefix-b"];
        assert_eq!(
            metadata_round_trip(&values, 2),
            values
                .iter()
                .map(|value| value.to_vec())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn paired_id_patterns_are_losslessly_derivable() {
        assert_eq!(
            PairedIdPattern::detect(b"read/1", b"read/2"),
            PairedIdPattern::SlashMate
        );
        assert_eq!(
            PairedIdPattern::SlashMate.derive_second(b"read/1").unwrap(),
            b"read/2"
        );
        assert_eq!(
            PairedIdPattern::detect(b"same", b"same"),
            PairedIdPattern::Identical
        );
        assert_eq!(
            PairedIdPattern::detect(b"left", b"right"),
            PairedIdPattern::Independent
        );
    }

    #[test]
    fn packed_dna_columns_round_trip_and_ambiguous_bases_fall_back() {
        let directory = tempfile::tempdir().unwrap();
        for (name, values, expected_mode) in [
            ("packed", vec![b"ACGTAC".as_slice(), b"TTGGCA"], 4u8),
            ("ambiguous", vec![b"ACGTAN".as_slice(), b"TTGGCA"], 5u8),
        ] {
            let mut column = StoredColumn::create(directory.path(), 0, name).unwrap();
            for value in &values {
                column.push(value).unwrap();
            }
            let mut encoded = Vec::new();
            column
                .write_packed_dna_to(&mut encoded, values.len())
                .unwrap();
            assert_eq!(encoded[0], expected_mode);
            let mut cursor = Cursor::new(&encoded);
            assert_eq!(cursor.column(values.len(), true, name).unwrap(), values);
        }
    }

    #[test]
    fn singleton_clusters_are_separated() {
        let clusters = (0..100)
            .map(|original_id| ConsensusCluster {
                reference: vec![b'A'; 150],
                reads: vec![ConsensusRead {
                    original_id,
                    offset: 0,
                    reverse: false,
                    length: 150,
                    mismatches: Vec::new(),
                }],
            })
            .collect::<Vec<_>>();
        let (multi, count, singletons, order) = partition_clusters_v8(&clusters).unwrap();
        assert!(multi.is_empty());
        assert_eq!(count, 0);
        assert_eq!(singletons.len(), clusters.len());
        assert_eq!(order.len(), clusters.len());
    }

    #[test]
    fn run_compaction_uses_bounded_fan_in() {
        let directory = tempfile::tempdir().unwrap();
        let runs_dir = directory.path().join("runs");
        std::fs::create_dir_all(&runs_dir).unwrap();
        let mut paths = Vec::new();
        for original_id in (0..(RUN_MERGE_FAN_IN as u64 + 1)).rev() {
            let path = runs_dir.join(format!("input-{original_id}.bin"));
            let mut file = std::io::BufWriter::new(File::create(&path).unwrap());
            file.write_all(&1u64.to_le_bytes()).unwrap();
            write_run_record(
                &mut file,
                original_id,
                &ReadRecord::new("id", b"ACGT".to_vec(), None),
            )
            .unwrap();
            file.flush().unwrap();
            paths.push(path);
        }
        let mut set = RunSet {
            paths,
            records: RUN_MERGE_FAN_IN as u64 + 1,
        };
        compact_run_set(&mut set, directory.path()).unwrap();
        assert!(set.paths.len() <= RUN_MERGE_FAN_IN);
        let mut merger = RunMerger::open(&set.paths).unwrap();
        for expected_id in 0..(RUN_MERGE_FAN_IN as u64 + 1) {
            assert_eq!(merger.next().unwrap().unwrap().0, expected_id);
        }
        assert!(merger.next().unwrap().is_none());
    }
}
