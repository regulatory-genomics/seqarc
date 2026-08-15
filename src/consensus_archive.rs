use crate::codec;
use crate::record_store::{PairStore, StoreStats};
use crate::reorder::{ConsensusCluster, ConsensusRead};
use crate::{
    ArchiveRecords, CompressionMetadata, CompressionOptions, CompressionProfile,
    DecompressionOptions, Error, ReadRecord, Result,
};
use crossbeam_channel::bounded;
use std::collections::BTreeMap;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

const MAGIC: &[u8; 9] = b"SPRNG-RS3";
const VERSION: u32 = 7;
const LONG_DISTANCE_WINDOW_LOG: u32 = 27;

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
        return Err(Error::InvalidFormat("bad RS3 magic".into()));
    }
    let version = read_u32(&mut input)?;
    if version != VERSION {
        return Err(Error::InvalidFormat("unsupported RS3 version".into()));
    }
    let streams = read_u16(&mut input)? as usize;
    let flags = read_u16(&mut input)?;
    let raw_len = usize_from(read_u64(&mut input)?, "payload")?;
    let encoded_len = usize_from(read_u64(&mut input)?, "payload")?;
    let expected = read_u32(&mut input)?;
    let mut encoded = vec![0; encoded_len];
    input.read_exact(&mut encoded)?;
    let mut decoder = codec::Decoder::new()?;
    let payload = decoder.decode(&encoded, raw_len)?;
    if checksum(&payload) != expected {
        return Err(Error::InvalidFormat("RS3 checksum failure".into()));
    }
    read_payload_v7(&payload, streams, flags, options)
}

fn read_payload_v7(
    payload: &[u8],
    streams: usize,
    flags: u16,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let mut cursor = Cursor::new(payload);
    let records = usize_from(cursor.u64()?, "record count")?;
    let payload_streams = cursor.u16()? as usize;
    if streams == 0 || streams != payload_streams {
        return Err(Error::InvalidFormat("RS3 stream count mismatch".into()));
    }
    let _quality_mode = cursor.u8()?;
    let id_pattern = PairedIdPattern::try_from(cursor.u8()?)?;
    let has_ids = flags & 1 != 0;
    let has_quality = flags & 2 != 0;
    if id_pattern != PairedIdPattern::Independent && (!has_ids || streams != 2) {
        return Err(Error::InvalidFormat(
            "paired ID pattern requires two ID-bearing streams".into(),
        ));
    }
    let mut decoded_streams = Vec::with_capacity(streams);
    for stream_index in 0..streams {
        let clusters = cursor.clusters_v7()?;
        let ids_stored =
            has_ids && (stream_index == 0 || id_pattern == PairedIdPattern::Independent);
        let mut metadata = cursor.columnar_records(ids_stored, has_quality)?;
        if metadata.len() != records {
            return Err(Error::InvalidFormat("RS3 record count mismatch".into()));
        }
        let mut seen = vec![false; records];
        for cluster in &clusters {
            for read in &cluster.reads {
                let original_id = usize_from(read.original_id, "original read ID")?;
                if original_id >= records || seen[original_id] {
                    return Err(Error::InvalidFormat("invalid RS3 original read ID".into()));
                }
                metadata[original_id].sequence =
                    crate::reorder::decode_consensus_read(cluster, read)?;
                seen[original_id] = true;
            }
        }
        if seen.iter().any(|seen| !seen) {
            return Err(Error::InvalidFormat("incomplete RS3 residuals".into()));
        }
        decoded_streams.push(metadata);
    }
    if has_ids && id_pattern != PairedIdPattern::Independent {
        let (first, rest) = decoded_streams.split_at_mut(1);
        for (first, second) in first[0].iter().zip(&mut rest[0]) {
            second.id = id_pattern.derive_second(&first.id)?;
        }
    }
    if cursor.position != payload.len() {
        return Err(Error::InvalidFormat("trailing RS3 payload data".into()));
    }
    finish_decoding(decoded_streams, options)
}

fn finish_decoding(
    mut streams: Vec<Vec<ReadRecord>>,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    if let Some(range) = options.range.as_ref() {
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

pub fn write_reordered_store_profiled<W: Write>(
    store: &mut PairStore,
    mut output: W,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<(CompressionMetadata, CompressionProfile)> {
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
    let mut payload = std::fs::OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(&payload_path)?;
    put_u64_writer(&mut payload, store.count())?;
    put_u16_writer(&mut payload, 2)?;
    payload.write_all(&[0])?;
    payload.write_all(&[id_pattern as u8])?;
    for (mate, (assignment_path, columns)) in assignments.iter().zip(record_columns).enumerate() {
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
        )?;
    }
    payload.flush()?;
    let raw_len = payload.seek(SeekFrom::End(0))?;
    payload.seek(SeekFrom::Start(0))?;
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
    loop {
        let bytes_read = payload.read(&mut transfer_buffer)?;
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

fn worker_count(workers: usize) -> Result<u32> {
    u32::try_from(workers).map_err(|_| Error::InvalidInput("thread count exceeds u32".into()))
}

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
    let mut multi_body = std::fs::File::create(&body_path)?;
    let singleton_ids_path = temp_dir.join(format!("clusters.{mate}.singleton.ids"));
    let mut singleton_ids = create_temporary_column(&singleton_ids_path)?;
    let mut singleton_sequences = StoredColumn::create(temp_dir, mate, "singleton-sequences")?;
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

    let (cluster_count, stats) = std::thread::scope(|scope| -> Result<((u32, u32), StoreStats)> {
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
                    let result = cluster_bucket(&mut worker_store, mate, reorder.kmer_size, &ids)
                        .and_then(|clusters| {
                            let (multi_bytes, multi_cluster_count, singletons) =
                                partition_clusters_v7(&clusters)?;
                            Ok(ConsensusBucket {
                                index,
                                multi_cluster_count,
                                multi_bytes,
                                singletons,
                                stats: worker_store.stats(),
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

        write_consensus_buckets(
            done_rx,
            permit_tx,
            &mut multi_body,
            &mut singleton_ids,
            &mut singleton_sequences,
        )
    })?;
    multi_body.flush()?;
    drop(multi_body);
    put_u32_writer(payload, cluster_count.0)?;
    let mut multi_body = std::fs::File::open(body_path)?;
    std::io::copy(&mut multi_body, payload)?;
    put_u32_writer(payload, cluster_count.1)?;
    singleton_ids.flush()?;
    let singleton_id_bytes = singleton_ids.seek(SeekFrom::End(0))?;
    singleton_ids.seek(SeekFrom::Start(0))?;
    put_u64_writer(payload, singleton_id_bytes)?;
    std::io::copy(&mut singleton_ids, payload)?;
    singleton_sequences.write_packed_dna_to(payload, cluster_count.1 as usize)?;
    Ok(stats)
}

struct ConsensusBucket {
    index: usize,
    multi_cluster_count: u32,
    multi_bytes: Vec<u8>,
    singletons: Vec<(u64, Vec<u8>)>,
    stats: StoreStats,
}

#[derive(Clone, Copy)]
struct AssignedRead {
    id: u64,
    position: u32,
    reverse: bool,
    anchored: bool,
}

fn produce_consensus_buckets(
    assignments: &Path,
    limit: usize,
    permit_rx: crossbeam_channel::Receiver<()>,
    job_tx: crossbeam_channel::Sender<(usize, Vec<AssignedRead>)>,
) -> Result<()> {
    let mut bucket = Vec::with_capacity(limit);
    let mut current_hash = None;
    let mut index = 0usize;
    for assignment in crate::reorder::AssignmentReader::open(assignments)? {
        let assignment = assignment?;
        if current_hash.is_some_and(|hash| hash != assignment.hash) || bucket.len() == limit {
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
        current_hash = Some(assignment.hash);
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

fn write_consensus_buckets<W: Write>(
    done_rx: crossbeam_channel::Receiver<Result<ConsensusBucket>>,
    permit_tx: crossbeam_channel::Sender<()>,
    multi_output: &mut W,
    singleton_ids: &mut std::fs::File,
    singleton_sequences: &mut StoredColumn,
) -> Result<((u32, u32), StoreStats)> {
    let mut pending = BTreeMap::new();
    let mut next_to_write = 0usize;
    let mut multi_cluster_count = 0u32;
    let mut singleton_count = 0u32;
    let mut previous_singleton_id = None;
    let mut stats = StoreStats::default();
    while let Ok(result) = done_rx.recv() {
        let bucket = result?;
        pending.insert(bucket.index, bucket);
        while let Some(bucket) = pending.remove(&next_to_write) {
            multi_output.write_all(&bucket.multi_bytes)?;
            multi_cluster_count = multi_cluster_count
                .checked_add(bucket.multi_cluster_count)
                .ok_or_else(|| Error::InvalidInput("too many consensus clusters".into()))?;
            for (original_id, sequence) in bucket.singletons {
                let encoded_id = previous_singleton_id
                    .map_or(original_id, |previous| zigzag_delta(original_id, previous));
                put_varint(singleton_ids, encoded_id)?;
                singleton_sequences.push(&sequence)?;
                previous_singleton_id = Some(original_id);
                singleton_count = singleton_count
                    .checked_add(1)
                    .ok_or_else(|| Error::InvalidInput("too many consensus clusters".into()))?;
            }
            stats.index_reads += bucket.stats.index_reads;
            stats.blob_reads += bucket.stats.blob_reads;
            stats.bytes_read += bucket.stats.bytes_read;
            let _ = permit_tx.send(());
            next_to_write += 1;
        }
    }
    if !pending.is_empty() {
        return Err(Error::InvalidInput(
            "consensus workers exited before all buckets were written".into(),
        ));
    }
    Ok(((multi_cluster_count, singleton_count), stats))
}

#[derive(Clone)]
struct ClusterInput {
    id: u64,
    sequence: Vec<u8>,
    reverse: bool,
    offset: i32,
}

fn cluster_bucket(
    store: &mut PairStore,
    mate: usize,
    kmer_size: usize,
    assignments: &[AssignedRead],
) -> Result<Vec<ConsensusCluster>> {
    let mut groups: Vec<(i32, Vec<ClusterInput>)> = Vec::new();
    for assignment in assignments {
        let mut sequence = store.read_sequence(assignment.id, mate)?;
        if assignment.reverse {
            sequence = reverse_complement(&sequence);
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
        let mut best = None;
        for (index, (group_anchor, group)) in groups.iter().enumerate() {
            let offset = *group_anchor - anchor;
            let (overlap, distance) = overlap_distance(&group[0].sequence, &sequence, offset);
            let minimum_overlap = group[0].sequence.len().min(sequence.len()).div_ceil(2);
            if overlap >= minimum_overlap && distance <= 4 {
                best = Some((index, offset));
                break;
            }
        }
        match best {
            Some((index, offset)) => groups[index].1.push(ClusterInput {
                id: assignment.id,
                sequence,
                reverse: assignment.reverse,
                offset,
            }),
            None => groups.push((
                anchor,
                vec![ClusterInput {
                    id: assignment.id,
                    sequence,
                    reverse: assignment.reverse,
                    offset: 0,
                }],
            )),
        }
    }
    groups
        .into_iter()
        .map(|(_, group)| finish_cluster(group))
        .collect()
}

fn overlap_distance(reference: &[u8], sequence: &[u8], offset: i32) -> (usize, usize) {
    let reference_start = offset.max(0) as usize;
    let sequence_start = (-offset).max(0) as usize;
    if reference_start >= reference.len() || sequence_start >= sequence.len() {
        return (0, 0);
    }
    let overlap = (reference.len() - reference_start).min(sequence.len() - sequence_start);
    let distance = reference[reference_start..reference_start + overlap]
        .iter()
        .zip(&sequence[sequence_start..sequence_start + overlap])
        .filter(|(left, right)| left != right)
        .count();
    (overlap, distance)
}

fn finish_cluster(mut group: Vec<ClusterInput>) -> Result<ConsensusCluster> {
    let minimum_offset = group.iter().map(|input| input.offset).min().unwrap_or(0);
    let maximum_end = group
        .iter()
        .map(|input| input.offset as i64 + input.sequence.len() as i64)
        .max()
        .unwrap_or(0);
    let reference_len = usize::try_from(maximum_end - minimum_offset as i64)
        .map_err(|_| Error::InvalidInput("consensus reference too large".into()))?;
    let mut counts = vec![[0u32; 5]; reference_len];
    for input in &group {
        let start = usize::try_from(input.offset - minimum_offset)
            .map_err(|_| Error::InvalidInput("negative consensus offset".into()))?;
        for (position, &base) in input.sequence.iter().enumerate() {
            counts[start + position][base_index(base)] += 1;
        }
    }
    let reference: Vec<u8> = counts
        .into_iter()
        .map(|count| {
            let mut best = 0;
            for index in 1..5 {
                if count[index] > count[best] {
                    best = index;
                }
            }
            [b'A', b'C', b'G', b'T', b'N'][best]
        })
        .collect();
    let mut reads = Vec::with_capacity(group.len());
    for input in &mut group {
        let offset = usize::try_from(input.offset - minimum_offset)
            .map_err(|_| Error::InvalidInput("negative consensus offset".into()))?;
        let mismatches = reference[offset..offset + input.sequence.len()]
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
    Ok(ConsensusCluster { reference, reads })
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

type SingletonSequence = (u64, Vec<u8>);

fn partition_clusters_v7(
    clusters: &[ConsensusCluster],
) -> Result<(Vec<u8>, u32, Vec<SingletonSequence>)> {
    let mut multi_bytes = Vec::new();
    let mut multi_cluster_count = 0u32;
    let mut singletons = Vec::new();
    for cluster in clusters {
        if cluster.reads.len() == 1 {
            let read = &cluster.reads[0];
            if read.offset != 0
                || read.length as usize != cluster.reference.len()
                || !read.mismatches.is_empty()
            {
                return Err(Error::InvalidInput(
                    "V5 singleton consensus invariant violated".into(),
                ));
            }
            let sequence = if read.reverse {
                reverse_complement(&cluster.reference)
            } else {
                cluster.reference.clone()
            };
            singletons.push((read.original_id, sequence));
        } else {
            put_cluster_bodies_v7(&mut multi_bytes, std::slice::from_ref(cluster))?;
            multi_cluster_count = multi_cluster_count
                .checked_add(1)
                .ok_or_else(|| Error::InvalidInput("too many consensus clusters".into()))?;
        }
    }
    Ok((multi_bytes, multi_cluster_count, singletons))
}

fn put_cluster_bodies_v7<W: Write>(writer: &mut W, clusters: &[ConsensusCluster]) -> Result<()> {
    for cluster in clusters {
        if cluster.reads.len() < 2 {
            return Err(Error::InvalidInput("V7 multi cluster is too small".into()));
        }
        put_reference_v7(writer, &cluster.reference)?;
        put_u32_writer(writer, cluster.reads.len() as u32)?;
        let mut previous_id = 0u64;
        for (index, read) in cluster.reads.iter().enumerate() {
            if index == 0 {
                put_varint(writer, read.original_id)?;
            } else {
                put_varint(writer, zigzag_delta(read.original_id, previous_id))?;
            }
            previous_id = read.original_id;
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

fn put_reference_v7<W: Write>(writer: &mut W, reference: &[u8]) -> Result<()> {
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

fn zigzag_delta(value: u64, previous: u64) -> u64 {
    if value >= previous {
        (value - previous) << 1
    } else {
        ((previous - value) << 1) - 1
    }
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

struct RecordColumns {
    ids: StoredColumn,
    qualities: StoredColumn,
}

impl RecordColumns {
    fn create(temp_dir: &Path, mate: usize) -> Result<Self> {
        Ok(Self {
            ids: StoredColumn::create(temp_dir, mate, "ids")?,
            qualities: StoredColumn::create(temp_dir, mate, "qualities")?,
        })
    }
}

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

fn write_store_record_columns<W: Write>(
    payload: &mut W,
    columns: RecordColumns,
    count: u64,
    ids: bool,
    quality: bool,
) -> Result<()> {
    let count =
        u32::try_from(count).map_err(|_| Error::InvalidInput("too many records for RS3".into()))?;
    put_u32_writer(payload, count)?;
    columns.ids.write_to(payload, ids, count as usize)?;
    columns
        .qualities
        .write_to(payload, quality, count as usize)?;
    Ok(())
}

struct StoredColumn {
    data: std::fs::File,
    lengths: std::fs::File,
    common_length: Option<u32>,
    variable_length: bool,
    data_len: u64,
    clean_dna: bool,
    clean_dna_records: u64,
}

impl StoredColumn {
    fn create(temp_dir: &Path, mate: usize, name: &str) -> Result<Self> {
        Ok(Self {
            data: create_temporary_column(&temp_dir.join(format!("records.{mate}.{name}.data")))?,
            lengths: create_temporary_column(
                &temp_dir.join(format!("records.{mate}.{name}.lengths")),
            )?,
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
        self.data.seek(SeekFrom::Start(0))?;
        if self.variable_length {
            payload.write_all(&[2])?;
            put_u64_writer(payload, self.data_len)?;
            self.lengths.flush()?;
            self.lengths.seek(SeekFrom::Start(0))?;
            std::io::copy(&mut self.lengths, payload)?;
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
        std::io::copy(&mut self.data, payload)?;
        Ok(())
    }

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
        self.data.seek(SeekFrom::Start(0))?;
        let mut sequence = vec![0; width];
        let mut packed = vec![0; packed_width];
        for _ in 0..count {
            self.data.read_exact(&mut sequence)?;
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
fn checksum(bytes: &[u8]) -> u32 {
    bytes.iter().fold(2166136261, |hash, byte| {
        (hash ^ *byte as u32).wrapping_mul(16777619)
    })
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
    fn multi_clusters_v7(&mut self) -> Result<Vec<ConsensusCluster>> {
        let count = self.u32()? as usize;
        let mut clusters = Vec::with_capacity(count);
        for _ in 0..count {
            let reference = self.reference_v7()?;
            let reads_count = self.u32()? as usize;
            if reads_count < 2 {
                return Err(Error::InvalidFormat("invalid V7 multi cluster".into()));
            }
            let mut reads = Vec::with_capacity(reads_count);
            let mut previous_id = 0u64;
            for index in 0..reads_count {
                let encoded_id = self.varint()?;
                let original_id = if index == 0 {
                    encoded_id
                } else {
                    decode_zigzag_delta(encoded_id, previous_id)?
                };
                previous_id = original_id;
                let reverse = self.u8()? != 0;
                let offset = u32::try_from(self.varint()?)
                    .map_err(|_| Error::InvalidFormat("V7 offset too large".into()))?;
                let length = u32::try_from(self.varint()?)
                    .map_err(|_| Error::InvalidFormat("V7 read length too large".into()))?;
                let mismatch_count = usize_from(self.varint()?, "mismatch count")?;
                let mut mismatches = Vec::with_capacity(mismatch_count);
                let mut position = 0u32;
                for _ in 0..mismatch_count {
                    let delta = u32::try_from(self.varint()?).map_err(|_| {
                        Error::InvalidFormat("V7 mismatch position too large".into())
                    })?;
                    position = position.checked_add(delta).ok_or_else(|| {
                        Error::InvalidFormat("V7 mismatch position overflow".into())
                    })?;
                    mismatches.push((position, self.u8()?));
                }
                let end = (offset as usize)
                    .checked_add(length as usize)
                    .ok_or_else(|| Error::InvalidFormat("V7 read range overflow".into()))?;
                if end > reference.len() {
                    return Err(Error::InvalidFormat("V7 read exceeds reference".into()));
                }
                reads.push(ConsensusRead {
                    original_id,
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

    fn clusters_v7(&mut self) -> Result<Vec<ConsensusCluster>> {
        let mut clusters = self.multi_clusters_v7()?;
        let singleton_count = self.u32()? as usize;
        let id_bytes_len = usize_from(self.u64()?, "singleton ID column")?;
        let id_bytes = self.take(id_bytes_len)?;
        let mut id_cursor = Cursor::new(id_bytes);
        let mut singleton_ids = Vec::with_capacity(singleton_count);
        let mut previous_id = 0u64;
        for index in 0..singleton_count {
            let encoded = id_cursor.varint()?;
            let original_id = if index == 0 {
                encoded
            } else {
                decode_zigzag_delta(encoded, previous_id)?
            };
            singleton_ids.push(original_id);
            previous_id = original_id;
        }
        if id_cursor.position != id_bytes.len() {
            return Err(Error::InvalidFormat("trailing V7 singleton ID data".into()));
        }
        let sequences = self.column(singleton_count, true, "singleton sequence")?;
        for (original_id, reference) in singleton_ids.into_iter().zip(sequences) {
            let length = u32::try_from(reference.len())
                .map_err(|_| Error::InvalidFormat("V7 reference too large".into()))?;
            clusters.push(ConsensusCluster {
                reference,
                reads: vec![ConsensusRead {
                    original_id,
                    offset: 0,
                    reverse: false,
                    length,
                    mismatches: Vec::new(),
                }],
            });
        }
        Ok(clusters)
    }

    fn reference_v7(&mut self) -> Result<Vec<u8>> {
        match self.u8()? {
            0 => self.bytes(),
            1 => {
                let width = self.u32()? as usize;
                let packed_width = width.div_ceil(4);
                let data_len = usize_from(self.u64()?, "V7 reference")?;
                if data_len != packed_width {
                    return Err(Error::InvalidFormat("invalid V7 packed reference".into()));
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
            _ => Err(Error::InvalidFormat("invalid V7 reference mode".into())),
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

fn decode_zigzag_delta(encoded: u64, previous: u64) -> Result<u64> {
    let magnitude = (encoded >> 1) + (encoded & 1);
    if encoded & 1 == 0 {
        previous
            .checked_add(magnitude)
            .ok_or_else(|| Error::InvalidFormat("RS3 read ID overflow".into()))
    } else {
        previous
            .checked_sub(magnitude)
            .ok_or_else(|| Error::InvalidFormat("RS3 read ID underflow".into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let (multi, count, singletons) = partition_clusters_v7(&clusters).unwrap();
        assert!(multi.is_empty());
        assert_eq!(count, 0);
        assert_eq!(singletons.len(), clusters.len());
    }
}
