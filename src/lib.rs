//! A Rust-native, library-only sequencing read compressor.
//!
//! The archive format is intentionally independent from the original C++
//! SPRING format. Archives are chunked, self-describing, and support range
//! decoding without a command-line interface.

mod archive;
mod codec;
mod consensus_archive;
mod fastx;
mod record_store;
pub mod reorder;

use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::ops::RangeInclusive;
use std::path::Path;
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
    pub chunk_records: usize,
    pub zstd_level: i32,
    pub reorder: ReorderOptions,
}

#[derive(Clone, Debug)]
pub struct ReorderOptions {
    pub enabled: bool,
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
            enabled: true,
            kmer_size: 21,
            window_size: 31,
            max_bucket_records: 4_096,
            sort_chunk_size: 2_000_000,
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
            chunk_records: 131_072,
            zstd_level: 19,
            reorder: ReorderOptions::default(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct DecompressionOptions {
    pub range: Option<RangeInclusive<u64>>,
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
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ArchiveRecords {
    pub streams: Vec<Vec<ReadRecord>>,
}

/// Compress one logical read stream into a Rust-native archive.
pub fn compress<I, W>(
    records: I,
    output: W,
    options: &CompressionOptions,
) -> Result<CompressionMetadata>
where
    I: IntoIterator<Item = Result<ReadRecord>>,
    W: Write,
{
    compress_streams(vec![records.into_iter()], output, options)
}

/// Compress paired or multiple logical streams. Streams must contain the same number of records.
pub fn compress_streams<I, W>(
    streams: Vec<I>,
    output: W,
    options: &CompressionOptions,
) -> Result<CompressionMetadata>
where
    I: IntoIterator<Item = Result<ReadRecord>>,
    W: Write,
{
    archive::write_archive(streams, output, options)
}

/// Stream paired gzipped FASTQ files into a native archive without buffering
/// the complete dataset in memory.
pub fn compress_paired_fastq_files(
    input_1: impl AsRef<Path>,
    input_2: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &CompressionOptions,
) -> Result<CompressionMetadata> {
    let reader_1 = FastxReader::new(
        BufReader::new(gzip_reader(BufReader::new(File::open(input_1)?))),
        Format::Fastq,
    );
    let reader_2 = FastxReader::new(
        BufReader::new(gzip_reader(BufReader::new(File::open(input_2)?))),
        Format::Fastq,
    );
    let archive = File::create(output)?;
    archive::write_archive_streams(vec![reader_1, reader_2], archive, options)
}

/// Reads paired FASTQ streams and writes an RS3 global-reorder archive.
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
    let ingest_start = std::time::Instant::now();
    let reader_1 = FastxReader::new(
        BufReader::new(gzip_reader(BufReader::new(File::open(input_1)?))),
        Format::Fastq,
    );
    let reader_2 = FastxReader::new(
        BufReader::new(gzip_reader(BufReader::new(File::open(input_2)?))),
        Format::Fastq,
    );
    let store_root = temp_dir.as_ref().join("spring-rs-record-store");
    let mut store = PairStore::create(&store_root)?;
    let mut first = reader_1;
    let mut second = reader_2;
    loop {
        match (first.next(), second.next()) {
            (None, None) => break,
            (Some(Ok(first)), Some(Ok(second))) => {
                first.validate(options.preserve_quality)?;
                second.validate(options.preserve_quality)?;
                store.append(&first, &second)?;
            }
            (None, Some(_)) | (Some(_), None) => {
                return Err(Error::InvalidInput(
                    "paired streams have different record counts".into(),
                ));
            }
            (Some(Err(error)), _) | (_, Some(Err(error))) => return Err(error),
        }
    }
    let ingest = ingest_start.elapsed();
    let archive = File::create(output)?;
    let (metadata, mut profile) =
        consensus_archive::write_reordered_store_profiled(&mut store, archive, options, temp_dir)?;
    profile.ingest = ingest;
    profile.store_bytes = directory_size(store.root())?;
    Ok((metadata, profile))
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
    let mut input = std::io::BufReader::new(input);
    let mut magic = [0; 9];
    input.read_exact(&mut magic)?;
    let mut data = magic.to_vec();
    input.read_to_end(&mut data)?;
    if &magic == b"SPRNG-RS3" {
        consensus_archive::read_reordered(std::io::Cursor::new(data), options)
    } else {
        archive::read_archive(std::io::Cursor::new(data), options)
    }
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
    fn round_trip_single_stream_and_range() {
        let input = vec![Ok(record(1)), Ok(record(2)), Ok(record(3))];
        let mut archive = Vec::new();
        let metadata = compress(input, &mut archive, &CompressionOptions::default()).unwrap();
        assert_eq!(metadata.records, 3);
        let decoded = decompress(
            Cursor::new(archive),
            &DecompressionOptions { range: Some(1..=1) },
        )
        .unwrap();
        assert_eq!(decoded.streams[0], vec![record(2)]);
    }

    #[test]
    fn defaults_use_sixteen_workers_and_large_chunks() {
        let options = CompressionOptions::default();
        assert_eq!(options.threads, 16);
        assert_eq!(options.chunk_records, 131_072);
        assert_eq!(options.zstd_level, 19);
    }

    #[test]
    fn round_trip_paired_streams() {
        let left = vec![Ok(record(1)), Ok(record(2))];
        let right = vec![Ok(record(3)), Ok(record(4))];
        let mut archive = Vec::new();
        compress_streams(
            vec![left, right],
            &mut archive,
            &CompressionOptions::default(),
        )
        .unwrap();
        let decoded = decompress(Cursor::new(archive), &DecompressionOptions::default()).unwrap();
        assert_eq!(decoded.streams.len(), 2);
        assert_eq!(decoded.streams[1][1], record(4));
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
    fn streaming_writer_round_trips_multiple_chunks() {
        let records: Vec<_> = (0..5).map(|n| Ok(record(n))).collect();
        let mates: Vec<_> = (5..10).map(|n| Ok(record(n))).collect();
        let mut archive = std::io::Cursor::new(Vec::new());
        let options = CompressionOptions {
            chunk_records: 2,
            ..CompressionOptions::default()
        };
        let metadata = archive::write_archive_streams(
            vec![records.into_iter(), mates.into_iter()],
            &mut archive,
            &options,
        )
        .unwrap();
        assert_eq!(metadata.records, 5);
        let decoded = decompress(
            std::io::Cursor::new(archive.into_inner()),
            &DecompressionOptions::default(),
        )
        .unwrap();
        assert_eq!(decoded.streams[0].len(), 5);
        assert_eq!(decoded.streams[1][4], record(9));
    }

    #[test]
    fn parallel_and_serial_produce_identical_archives() {
        let encode = |threads: usize| {
            let left: Vec<_> = (0..50).map(|n| Ok(record(n))).collect();
            let right: Vec<_> = (50..100).map(|n| Ok(record(n))).collect();
            let options = CompressionOptions {
                chunk_records: 7,
                threads,
                ..CompressionOptions::default()
            };
            let mut archive = std::io::Cursor::new(Vec::new());
            archive::write_archive_streams(
                vec![left.into_iter(), right.into_iter()],
                &mut archive,
                &options,
            )
            .unwrap();
            archive.into_inner()
        };

        let serial = encode(1);
        let parallel = encode(8);
        assert_eq!(serial, parallel, "archive bytes must be deterministic");

        let decoded = decompress(
            std::io::Cursor::new(parallel),
            &DecompressionOptions::default(),
        )
        .unwrap();
        assert_eq!(decoded.streams[0].len(), 50);
        assert_eq!(decoded.streams[1][49], record(99));
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
            &DecompressionOptions { range: Some(1..=1) },
        )
        .unwrap();
        assert_eq!(decoded.streams[0], vec![first[1].clone()]);
        assert_eq!(decoded.streams[1], vec![second[1].clone()]);
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
                threads: 4,
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
    fn reordered_file_api_uses_store_and_forced_spills() {
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
            reorder: ReorderOptions {
                kmer_size: 3,
                window_size: 3,
                sort_chunk_size: 1,
                max_bucket_records: 2,
                ..ReorderOptions::default()
            },
            zstd_level: 3,
            ..CompressionOptions::default()
        };
        compress_reordered_paired_fastq_files(
            &input_1,
            &input_2,
            &output,
            &options,
            directory.path().join("work"),
        )
        .unwrap();
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
