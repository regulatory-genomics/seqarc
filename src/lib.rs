//! A Rust-native, library-only sequencing read compressor.
//!
//! The archive format is intentionally independent from the original C++
//! SPRING format. RS4 archives are self-describing and support range decoding
//! without a command-line interface.

mod consensus_archive;
mod external_bytes;
mod fastx;
mod packed_sequence;
mod quality_blocks;
mod record_store;
pub mod reorder;
#[allow(dead_code)]
mod sequence_spool;

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::ops::RangeInclusive;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Duration;

pub use fastx::{gzip_reader, gzip_writer, FastxReader, FastxWriter, Format};
pub use record_store::{PairStore, StoreStats};
pub use reorder::{
    canonical_minimizers, decode_consensus_read, ConsensusCluster, ConsensusRead, MinimizerEntry,
    OrderedRead,
};

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    InvalidFormat(String),
    InvalidInput(String),
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

    pub(crate) fn validate(&self, preserve_quality: bool) -> Result<()> {
        if self.sequence.is_empty() {
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
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum QualityMode {
    Lossless,
    IlluminaBinning,
    Binary { threshold: u8, high: u8, low: u8 },
}

#[derive(Clone, Debug)]
pub struct CompressionOptions {
    pub preserve_ids: bool,
    pub preserve_quality: bool,
    pub preserve_order: bool,
    pub quality_mode: QualityMode,
    pub threads: usize,
    pub zstd_threads: usize,
    pub zstd_level: i32,
    pub reorder: ReorderOptions,
}

#[derive(Clone, Debug)]
pub struct ReorderOptions {
    pub kmer_size: usize,
    pub window_size: usize,
    pub max_bucket_records: usize,
    pub sort_chunk_size: usize,
    pub sort_threads: usize,
    pub temp_dir: Option<std::path::PathBuf>,
}

impl Default for ReorderOptions {
    fn default() -> Self {
        Self {
            kmer_size: 21,
            window_size: 31,
            max_bucket_records: 4_096,
            sort_chunk_size: 10_000_000,
            sort_threads: 16,
            temp_dir: None,
        }
    }
}

impl Default for CompressionOptions {
    fn default() -> Self {
        Self {
            preserve_ids: true,
            preserve_quality: true,
            preserve_order: true,
            quality_mode: QualityMode::Lossless,
            threads: 16,
            zstd_threads: 16,
            zstd_level: 19,
            reorder: ReorderOptions::default(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct DecompressionOptions {
    pub original_range: Option<RangeInclusive<u64>>,
    pub output_order: OutputOrder,
    pub temp_dir: Option<std::path::PathBuf>,
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

#[derive(Clone, Debug, Default)]
pub struct CompressionProfile {
    pub ingest: Duration,
    pub assignments: Duration,
    pub consensus_and_payload: Duration,
    pub checksum: Duration,
    pub zstd: Duration,
    pub store_stats: StoreStats,
    pub store_bytes: u64,
    pub assignment_bytes: u64,
    pub payload_bytes: u64,
    pub archive_bytes: u64,
    pub consensus_producer: Duration,
    pub consensus_producer_wait: Duration,
    pub consensus_workers: Duration,
    pub consensus_writer: Duration,
    pub consensus_writer_active: Duration,
    pub consensus_writer_wait: Duration,
    pub consensus_max_pending_buckets: u64,
    pub consensus_multi_clusters: u64,
    pub consensus_singletons: u64,
    pub consensus_buckets: u64,
    pub consensus_overlap_comparisons: u64,
    pub consensus_max_bucket_records: u64,
    pub consensus_max_groups: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArchiveRecords {
    pub streams: Vec<Vec<ReadRecord>>,
}

/// Compress paired gzipped FASTQ files into a native archive.
///
/// Uses the store-backed RS4 consensus path.
pub fn compress_paired_fastq_files(
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &CompressionOptions,
) -> Result<CompressionMetadata> {
    let temp_dir =
        options.reorder.temp_dir.clone().unwrap_or_else(|| {
            std::env::temp_dir().join(format!("spring-rs-{}", std::process::id()))
        });
    compress_reordered_paired_fastq_files(input_1, input_2, output, options, temp_dir)
}

/// Reads paired FASTQ streams and writes an RS4 global-reorder archive.
///
/// Minimizer entries spill through `bed-utils`, records are stored on disk,
/// and consensus buckets are processed by a bounded worker pool.
pub fn compress_reordered_paired_fastq_files(
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<CompressionMetadata> {
    compress_reordered_paired_fastq_files_profiled(input_1, input_2, output, options, temp_dir)
        .map(|(metadata, _)| metadata)
}

pub fn compress_reordered_paired_fastq_files_profiled(
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &CompressionOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<(CompressionMetadata, CompressionProfile)> {
    let temp_root = temp_dir.as_ref();
    reorder::validate_options(&options.reorder)?;
    let output = output.as_ref();
    let partial_output =
        output.with_extension(match output.extension().and_then(|value| value.to_str()) {
            Some(extension) => format!("{extension}.partial"),
            None => "partial".into(),
        });
    let mut guard = CompressionGuard::new(create_workspace(temp_root)?, partial_output);
    let temp_dir = guard.workspace().to_path_buf();
    let ingest_start = std::time::Instant::now();
    let mut prepared = consensus_archive::PreparedRecordsBuilder::create_metadata(&temp_dir)?;
    let input_1 = input_1.as_ref().to_path_buf();
    let input_2 = input_2.as_ref().to_path_buf();
    let workers = options.threads.max(1);
    let queue_capacity = workers.saturating_mul(2).max(1);
    let (job_tx, job_rx) =
        crossbeam_channel::bounded::<(usize, Vec<(u64, ReadRecord, ReadRecord)>)>(queue_capacity);
    let (result_tx, result_rx) = crossbeam_channel::bounded::<
        Result<(usize, Vec<consensus_archive::PreparedPair>)>,
    >(queue_capacity);
    let (sort_0_tx, sort_0_rx) = crossbeam_channel::bounded(queue_capacity);
    let (sort_1_tx, sort_1_rx) = crossbeam_channel::bounded(queue_capacity);
    let mut pending = BTreeMap::new();
    let mut next_batch = 0usize;
    let sections = std::thread::scope(
        |scope| -> Result<[consensus_archive::PreparedMateSection; 2]> {
            let sort_options = options.clone();
            let sort_temp_dir = temp_dir.clone();
            let sort_0_handle = scope.spawn(move || {
                consensus_archive::sort_and_write_mate(
                    sort_0_rx,
                    u64::MAX,
                    0,
                    &sort_options,
                    &sort_temp_dir,
                )
            });
            let sort_options = options.clone();
            let sort_temp_dir = temp_dir.clone();
            let sort_1_handle = scope.spawn(move || {
                consensus_archive::sort_and_write_mate(
                    sort_1_rx,
                    u64::MAX,
                    1,
                    &sort_options,
                    &sort_temp_dir,
                )
            });
            let producer_results = result_tx.clone();
            let producer_options = options.clone();
            scope.spawn(move || {
                let result = (|| -> Result<()> {
                    let first = FastxReader::new(
                        BufReader::with_capacity(
                            8 * 1024 * 1024,
                            gzip_reader(BufReader::with_capacity(
                                1024 * 1024,
                                File::open(input_1)?,
                            )),
                        ),
                        Format::Fastq,
                    );
                    let second = FastxReader::new(
                        BufReader::with_capacity(
                            8 * 1024 * 1024,
                            gzip_reader(BufReader::with_capacity(
                                1024 * 1024,
                                File::open(input_2)?,
                            )),
                        ),
                        Format::Fastq,
                    );
                    let mut first = first;
                    let mut second = second;
                    let mut batch = Vec::with_capacity(2048);
                    let mut batch_index = 0usize;
                    let mut original_id = 0u64;
                    loop {
                        match (first.next(), second.next()) {
                            (None, None) => break,
                            (Some(Ok(mut first)), Some(Ok(mut second))) => {
                                first.validate(true)?;
                                second.validate(true)?;
                                if producer_options.preserve_quality {
                                    apply_quality(&mut first, &producer_options.quality_mode);
                                    apply_quality(&mut second, &producer_options.quality_mode);
                                }
                                batch.push((original_id, first, second));
                                original_id += 1;
                                if batch.len() == 2048 {
                                    job_tx
                                        .send((batch_index, std::mem::take(&mut batch)))
                                        .map_err(|_| {
                                            Error::InvalidInput("preparation workers exited".into())
                                        })?;
                                    batch_index += 1;
                                    batch = Vec::with_capacity(2048);
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
                    if !batch.is_empty() {
                        job_tx.send((batch_index, batch)).map_err(|_| {
                            Error::InvalidInput("preparation workers exited".into())
                        })?;
                    }
                    Ok(())
                })();
                if let Err(error) = result {
                    let _ = producer_results.send(Err(error));
                }
            });
            for _ in 0..workers {
                let job_rx = job_rx.clone();
                let result_tx = result_tx.clone();
                let worker_options = options.clone();
                scope.spawn(move || {
                    while let Ok((batch_index, batch)) = job_rx.recv() {
                        let result = batch
                            .into_iter()
                            .map(|(original_id, first, second)| {
                                consensus_archive::prepare_pair(
                                    first,
                                    second,
                                    original_id,
                                    &worker_options,
                                )
                            })
                            .collect::<Result<Vec<_>>>()
                            .map(|pairs| (batch_index, pairs));
                        if result_tx.send(result).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(job_rx);
            drop(result_tx);
            while let Ok(result) = result_rx.recv() {
                let (batch_index, pairs) = result?;
                pending.insert(batch_index, pairs);
                while let Some(pairs) = pending.remove(&next_batch) {
                    for pair in pairs {
                        let [read_0, read_1] = prepared.push_metadata(pair, options)?;
                        sort_0_tx
                            .send(read_0)
                            .map_err(|_| Error::InvalidInput("mate-1 sorter exited".into()))?;
                        sort_1_tx
                            .send(read_1)
                            .map_err(|_| Error::InvalidInput("mate-2 sorter exited".into()))?;
                    }
                    next_batch += 1;
                }
            }
            drop(sort_0_tx);
            drop(sort_1_tx);
            let section_0 = sort_0_handle
                .join()
                .map_err(|_| Error::InvalidInput("mate-1 sorter panicked".into()))??;
            let section_1 = sort_1_handle
                .join()
                .map_err(|_| Error::InvalidInput("mate-2 sorter panicked".into()))??;
            Ok([section_0, section_1])
        },
    )?;
    let prepared = prepared.finish_metadata()?;
    let ingest = ingest_start.elapsed();
    let archive = File::create(guard.partial_output())?;
    let (metadata, mut profile) = consensus_archive::write_streamed_profiled(
        prepared, sections, archive, options, &temp_dir,
    )?;
    profile.store_bytes = directory_size(&temp_dir)?;
    profile.ingest = ingest;
    guard.commit(output)?;
    Ok((metadata, profile))
}

struct CompressionGuard {
    workspace: Option<std::path::PathBuf>,
    partial_output: Option<std::path::PathBuf>,
}

impl CompressionGuard {
    fn new(workspace: std::path::PathBuf, partial_output: std::path::PathBuf) -> Self {
        Self {
            workspace: Some(workspace),
            partial_output: Some(partial_output),
        }
    }

    fn workspace(&self) -> &Path {
        self.workspace.as_deref().unwrap()
    }

    fn partial_output(&self) -> &Path {
        self.partial_output.as_deref().unwrap()
    }

    fn commit(&mut self, output: &Path) -> Result<()> {
        if let Some(workspace) = self.workspace.take() {
            if let Err(error) = std::fs::remove_dir_all(&workspace) {
                self.workspace = Some(workspace);
                return Err(Error::Io(error));
            }
        }
        let partial_output = self.partial_output.take().unwrap();
        if let Err(error) = std::fs::rename(&partial_output, output) {
            self.partial_output = Some(partial_output);
            return Err(Error::Io(error));
        }
        Ok(())
    }
}

impl Drop for CompressionGuard {
    fn drop(&mut self) {
        if let Some(workspace) = self.workspace.take() {
            let _ = std::fs::remove_dir_all(workspace);
        }
        if let Some(partial_output) = self.partial_output.take() {
            let _ = std::fs::remove_file(partial_output);
        }
    }
}

fn create_workspace(root: &Path) -> Result<std::path::PathBuf> {
    static NEXT_WORKSPACE: AtomicU64 = AtomicU64::new(0);

    std::fs::create_dir_all(root)?;
    for _ in 0..1000 {
        let sequence = NEXT_WORKSPACE.fetch_add(1, Ordering::Relaxed);
        let workspace = root.join(format!("run-{}-{sequence}", std::process::id()));
        match std::fs::create_dir(&workspace) {
            Ok(()) => return Ok(workspace),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(Error::Io(error)),
        }
    }
    Err(Error::InvalidInput(
        "unable to allocate a unique compression workspace".into(),
    ))
}

fn directory_size(path: &Path) -> Result<u64> {
    let mut total = 0;
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        if entry.file_type()?.is_file() {
            total += entry.metadata()?.len();
        }
    }
    Ok(total)
}

/// Decompress an archive, optionally selecting an inclusive zero-based record range.
pub fn decompress<R: Read>(input: R, options: &DecompressionOptions) -> Result<ArchiveRecords> {
    consensus_archive::read_reordered(input, options)
}

/// Decode paired records into a sink without requiring callers to collect the
/// output into `ArchiveRecords`.
///
/// Records are delivered in original order.
pub fn decompress_paired_into<R: Read, S: PairedRecordSink>(
    input: R,
    options: &DecompressionOptions,
    sink: &mut S,
) -> Result<DecompressionMetadata> {
    consensus_archive::read_paired_into(input, options, sink)
}

/// Copy decoded records to FASTQ or FASTA writers.
pub fn write_fastx<W: Write>(records: &[ReadRecord], output: W, format: Format) -> Result<()> {
    let mut writer = FastxWriter::new(output, format);
    for record in records {
        writer.write_record(record)?;
    }
    writer.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;
    use std::io::Cursor;

    fn record(n: u8) -> ReadRecord {
        ReadRecord::new(
            format!("read{n}"),
            b"ACGTN".to_vec(),
            Some(b"IIIII".to_vec()),
        )
    }

    #[test]
    fn defaults_use_sixteen_workers_and_zstd_threads() {
        let options = CompressionOptions::default();
        assert_eq!(options.threads, 16);
        assert_eq!(options.zstd_threads, 16);
        assert_eq!(options.zstd_level, 19);
    }

    #[test]
    fn illumina_binning_matches_spring_table() {
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
        write_fastx(&records, &mut output, Format::Fastq).unwrap();
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
    fn store_backed_reordered_archive_round_trips() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        let first = [
            ReadRecord::new("a", b"ACGTACGT".to_vec(), Some(b"IIIIIIII".to_vec())),
            ReadRecord::new("b", b"ACGTACGA".to_vec(), Some(b"HHHHHHHH".to_vec())),
        ];
        let second = [
            ReadRecord::new("a/2", b"TGCATGCA".to_vec(), Some(b"JJJJJJJJ".to_vec())),
            ReadRecord::new("b/2", b"TGCATGCT".to_vec(), Some(b"GGGGGGGG".to_vec())),
        ];
        for (left, right) in first.iter().zip(&second) {
            store.append(left, right).unwrap();
        }
        let mut archive = Vec::new();
        let (_, profile) = consensus_archive::write_reordered_store_profiled(
            &mut store,
            &mut archive,
            &CompressionOptions::default(),
            directory.path(),
        )
        .unwrap();
        assert_eq!(profile.store_stats.index_reads, 10);
        assert_eq!(profile.store_stats.blob_reads, 16);
        let decoded = decompress(
            Cursor::new(archive),
            &DecompressionOptions {
                original_range: Some(1..=1),
                ..DecompressionOptions::default()
            },
        )
        .unwrap();
        assert_eq!(decoded.streams[0], vec![first[1].clone()]);
        assert_eq!(decoded.streams[1], vec![second[1].clone()]);
    }

    #[test]
    fn rs4_round_trips_singletons_ids_qualities_and_original_order() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        let first = [
            ReadRecord::new("z/1", b"ACGT".to_vec(), Some(b"IIII".to_vec())),
            ReadRecord::new("a/1", b"NNNN".to_vec(), Some(b"####".to_vec())),
            ReadRecord::new("m/1", b"TGCA".to_vec(), Some(b"JJJJ".to_vec())),
        ];
        let second = [
            ReadRecord::new("z/2", b"GATC".to_vec(), Some(b"HHHH".to_vec())),
            ReadRecord::new("a/2", b"CCCC".to_vec(), Some(b"$$$$".to_vec())),
            ReadRecord::new("m/2", b"ATAT".to_vec(), Some(b"GGGG".to_vec())),
        ];
        for (left, right) in first.iter().zip(&second) {
            store.append(left, right).unwrap();
        }
        let mut archive = Vec::new();
        consensus_archive::write_reordered_store_profiled(
            &mut store,
            &mut archive,
            &CompressionOptions::default(),
            directory.path(),
        )
        .unwrap();
        let decoded = decompress(Cursor::new(archive), &DecompressionOptions::default()).unwrap();
        assert_eq!(decoded.streams, vec![first.to_vec(), second.to_vec()]);
    }

    #[test]
    fn rs4_round_trips_prefixed_arbitrary_byte_ids() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        let first = [
            ReadRecord::new(b"\x80sample-000/1".to_vec(), b"ACGT".to_vec(), None),
            ReadRecord::new(b"\x80sample-001/1".to_vec(), b"NNNN".to_vec(), None),
            ReadRecord::new(b"\x80sample-002/1".to_vec(), b"TGCA".to_vec(), None),
        ];
        let second = [
            ReadRecord::new(b"\x80sample-000/2".to_vec(), b"TGCA".to_vec(), None),
            ReadRecord::new(b"\x80sample-001/2".to_vec(), b"CCCC".to_vec(), None),
            ReadRecord::new(b"\x80sample-002/2".to_vec(), b"GATC".to_vec(), None),
        ];
        for (left, right) in first.iter().zip(&second) {
            store.append(left, right).unwrap();
        }
        let mut archive = Vec::new();
        consensus_archive::write_reordered_store_profiled(
            &mut store,
            &mut archive,
            &CompressionOptions {
                preserve_quality: false,
                ..CompressionOptions::default()
            },
            directory.path(),
        )
        .unwrap();
        let decoded = decompress(Cursor::new(archive), &DecompressionOptions::default()).unwrap();
        assert_eq!(decoded.streams, vec![first.to_vec(), second.to_vec()]);
    }

    #[test]
    fn store_backed_rs3_round_trips_with_multiple_zstd_workers() {
        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        for index in 0..32 {
            let first = record(index);
            let second = record(index.wrapping_add(32));
            store.append(&first, &second).unwrap();
        }
        let mut archive = Vec::new();
        consensus_archive::write_reordered_store_profiled(
            &mut store,
            &mut archive,
            &CompressionOptions {
                threads: 1,
                zstd_threads: 4,
                zstd_level: 19,
                reorder: ReorderOptions {
                    kmer_size: 3,
                    window_size: 3,
                    ..ReorderOptions::default()
                },
                ..CompressionOptions::default()
            },
            directory.path(),
        )
        .unwrap();
        let decoded = decompress(Cursor::new(archive), &DecompressionOptions::default()).unwrap();
        assert_eq!(decoded.streams[0], (0..32).map(record).collect::<Vec<_>>());
        assert_eq!(decoded.streams[1], (32..64).map(record).collect::<Vec<_>>());
    }

    #[test]
    fn paired_sink_receives_original_order_records() {
        struct Sink {
            ids: Vec<u64>,
        }

        impl PairedRecordSink for Sink {
            fn write_pair(
                &mut self,
                original_id: u64,
                _first: ReadRecord,
                _second: ReadRecord,
            ) -> Result<()> {
                self.ids.push(original_id);
                Ok(())
            }

            fn finish(&mut self) -> Result<()> {
                Ok(())
            }
        }

        let directory = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(directory.path().join("store")).unwrap();
        for index in 0..4 {
            store.append(&record(index), &record(index + 4)).unwrap();
        }
        let mut archive = Vec::new();
        consensus_archive::write_reordered_store_profiled(
            &mut store,
            &mut archive,
            &CompressionOptions::default(),
            directory.path(),
        )
        .unwrap();

        let mut sink = Sink { ids: Vec::new() };
        let metadata = decompress_paired_into(
            Cursor::new(archive),
            &DecompressionOptions::default(),
            &mut sink,
        )
        .unwrap();
        assert_eq!(metadata.pairs, 4);
        assert_eq!(sink.ids, vec![0, 1, 2, 3]);
    }

    #[test]
    fn reordered_file_api_uses_sequence_spools_and_forced_spills() {
        use std::io::Write;

        let directory = tempfile::tempdir().unwrap();
        let input_1 = directory.path().join("reads.1.fastq.gz");
        let input_2 = directory.path().join("reads.2.fastq.gz");
        let output = directory.path().join("reads.spring-rs");
        let mut encoder_1 = gzip_writer(File::create(&input_1).unwrap(), 1);
        let mut encoder_2 = gzip_writer(File::create(&input_2).unwrap(), 1);
        let first = [b"ACGTACGT".as_slice(), b"ACGTACGA", b"ACGTACGG"];
        let second = [b"TGCATGCA".as_slice(), b"TGCATGCT", b"TGCATGCC"];
        for (index, sequence) in first.iter().enumerate() {
            writeln!(encoder_1, "@r{index}/1").unwrap();
            encoder_1.write_all(sequence).unwrap();
            encoder_1.write_all(b"\n+\nIIIIIIII\n").unwrap();
            writeln!(encoder_2, "@r{index}/2").unwrap();
            encoder_2.write_all(second[index]).unwrap();
            encoder_2.write_all(b"\n+\nJJJJJJJJ\n").unwrap();
        }
        encoder_1.finish().unwrap();
        encoder_2.finish().unwrap();
        let options = CompressionOptions {
            quality_mode: QualityMode::IlluminaBinning,
            reorder: ReorderOptions {
                kmer_size: 3,
                window_size: 3,
                sort_chunk_size: 1,
                max_bucket_records: 2,
                temp_dir: Some(directory.path().join("work")),
                ..ReorderOptions::default()
            },
            zstd_level: 3,
            ..CompressionOptions::default()
        };
        compress_paired_fastq_files(&input_1, &input_2, &output, &options).unwrap();
        let work = options.reorder.temp_dir.as_ref().unwrap();
        assert!(!work.join("spring-rs-record-store").exists());
        assert!(!work.join("sequences.0.spool").exists());
        assert!(!work.join("sequences.1.spool").exists());
        assert_eq!(std::fs::read_dir(work).unwrap().count(), 0);
        assert!(!output.with_extension("spring-rs.partial").exists());
        let decoded = decompress(
            BufReader::new(File::open(output).unwrap()),
            &DecompressionOptions::default(),
        )
        .unwrap();
        assert_eq!(
            decoded.streams[0]
                .iter()
                .map(|record| record.sequence.as_slice())
                .collect::<Vec<_>>(),
            first
        );
        assert_eq!(
            decoded.streams[1]
                .iter()
                .map(|record| record.sequence.as_slice())
                .collect::<Vec<_>>(),
            second
        );
        assert!(decoded.streams[0]
            .iter()
            .all(|record| record.quality.as_deref() == Some(b"IIIIIIII")));
        assert!(decoded.streams[1]
            .iter()
            .all(|record| record.quality.as_deref() == Some(b"IIIIIIII")));
    }

    #[test]
    fn failed_file_compression_cleans_workspace_and_partial_output() {
        use std::io::Write;

        let directory = tempfile::tempdir().unwrap();
        let input_1 = directory.path().join("reads.1.fastq.gz");
        let input_2 = directory.path().join("reads.2.fastq.gz");
        let output = directory.path().join("reads.spring-rs");
        let work = directory.path().join("work");
        let mut first = gzip_writer(File::create(&input_1).unwrap(), 1);
        let mut second = gzip_writer(File::create(&input_2).unwrap(), 1);
        first.write_all(b"@r0/1\nACGT\n+\nIIII\n").unwrap();
        second.write_all(b"").unwrap();
        first.finish().unwrap();
        second.finish().unwrap();
        let options = CompressionOptions {
            reorder: ReorderOptions {
                temp_dir: Some(work.clone()),
                ..ReorderOptions::default()
            },
            ..CompressionOptions::default()
        };

        assert!(compress_paired_fastq_files(input_1, input_2, &output, &options).is_err());
        assert!(!output.exists());
        assert!(!output.with_extension("spring-rs.partial").exists());
        assert_eq!(std::fs::read_dir(work).unwrap().count(), 0);
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
