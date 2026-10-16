use clap::{Args, Parser, Subcommand, ValueEnum};
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressStyle};
use seqarc::{
    build_alignment_index, compress_aligned_paired_fastq_files_with_index_profiled,
    compress_aligned_single_fastq_file_with_index_profiled, decompress_aligned_into,
    decompress_aligned_iter, default_thread_count, open_fastx_reader,
    read_archive_info_with_breakdown, select_id_field, truncate_read_id, verify_aligned_archive,
    AlignedRecordGroup, AlignedRecordSink, AlignerKind, AlignmentCompressionOptions,
    AlignmentLayout, ArchiveInfo, DecompressionOptions, EmbeddedIndexMate, EmbeddedIndexRange,
    Error, FastxReader, Format, GenomeIndexBuildOptions, QualityMode, ReadRecord, Result,
    SectionBreakdown, SequenceRecordGroup,
};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

#[derive(Parser)]
#[command(name = "seqarc", version, about = "Archive sequencing reads")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Build a single-file SeqArc genome index.
    Index(IndexArgs),
    /// Compress single-end or paired plain, gzip, or zstd FASTQ files.
    Compress(CompressArgs),
    /// Decompress an aligned SeqArc archive.
    Decompress(DecompressArgs),
    /// Fully decode and validate an aligned SeqArc archive.
    Validate(ValidateArgs),
    /// Verify chunk checksums and framing without decompressing.
    Verify(VerifyArgs),
    /// Print archive metadata without decoding records.
    Info(InfoArgs),
}

#[derive(Args)]
struct IndexArgs {
    /// Input plain, gzip, or zstd FASTA reference.
    reference: PathBuf,
    /// Output .seqarc-index file.
    #[arg(short, long)]
    output: PathBuf,
    /// Replace an existing output file.
    #[arg(long)]
    force: bool,
    /// Directory in which to create temporary index-build files.
    #[arg(long)]
    temp_dir: Option<PathBuf>,
    /// Alignment backend used by this index.
    #[arg(long, value_enum, default_value_t = AlignerChoice::Minibwa)]
    aligner: AlignerChoice,
    /// GTF annotation used to insert known splice junctions into a rustar index.
    #[arg(long)]
    gtf: Option<PathBuf>,
    /// Splice-junction overhang stored in a rustar index.
    #[arg(long, default_value_t = 100)]
    sjdb_overhang: u32,
    /// Rustar suffix-array sparsity. Larger values reduce index memory.
    #[arg(long, default_value_t = 1)]
    sa_sparse_d: u32,
    /// Rustar suffix-array prefix length.
    #[arg(long, default_value_t = 14)]
    sa_index_nbases: u32,
    /// Rustar chromosome-bin size as a base-2 exponent.
    #[arg(long, default_value_t = 18)]
    chr_bin_nbits: u32,
    /// Number of index construction threads.
    #[arg(long, default_value_t = default_thread_count())]
    threads: usize,
}

#[derive(Clone, Copy, ValueEnum)]
enum AlignerChoice {
    Minibwa,
    Rustar,
}

impl From<AlignerChoice> for AlignerKind {
    fn from(value: AlignerChoice) -> Self {
        match value {
            AlignerChoice::Minibwa => AlignerKind::Minibwa,
            AlignerChoice::Rustar => AlignerKind::Rustar,
        }
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum QualityChoice {
    Lossless,
    Illumina,
    Binary,
    Qvz,
}

#[derive(Args)]
struct CompressArgs {
    /// Input .seqarc-index file.
    #[arg(short, long)]
    index: PathBuf,
    /// Input plain, gzip, or zstd FASTQ; mate 1 for paired input.
    #[arg(short = '1', long)]
    mate_1: PathBuf,
    /// Optional plain, gzip, or zstd mate-2 FASTQ. Omit for single-end input.
    #[arg(short = '2', long)]
    mate_2: Option<PathBuf>,
    /// Extract an embedded index using one-based inclusive MATE:START:END coordinates.
    #[arg(long = "index-range", value_parser = parse_index_range)]
    index_ranges: Vec<EmbeddedIndexRange>,
    /// Optional plain, gzip, or zstd index-read FASTQ (I1).
    #[arg(long = "index-read-1")]
    index_read_1: Option<PathBuf>,
    /// Optional plain, gzip, or zstd index-read FASTQ (I2). Requires --index-read-1.
    #[arg(long = "index-read-2", requires = "index_read_1")]
    index_read_2: Option<PathBuf>,
    /// Output .seqarc archive.
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long, default_value_t = default_thread_count())]
    threads: usize,
    #[arg(long, default_value_t = 250_000)]
    chunk_records: usize,
    #[arg(long, default_value_t = 19)]
    zstd_level: i32,
    #[arg(long)]
    no_ids: bool,
    /// Retain at most the first N identifier characters, excluding FASTQ '@'.
    #[arg(long, conflicts_with = "no_ids")]
    id_prefix_chars: Option<usize>,
    /// Keep only the Nth whitespace-delimited identifier field (zero-based).
    #[arg(long, conflicts_with = "no_ids")]
    id_field: Option<usize>,
    #[arg(long)]
    no_quality: bool,
    #[arg(long, value_enum, default_value_t = QualityChoice::Lossless)]
    quality_mode: QualityChoice,
    #[arg(long, requires = "binary_high", requires = "binary_low")]
    binary_threshold: Option<u8>,
    #[arg(long, requires = "binary_threshold", requires = "binary_low")]
    binary_high: Option<u8>,
    #[arg(long, requires = "binary_threshold", requires = "binary_high")]
    binary_low: Option<u8>,
    /// Target entropy in bits per quality symbol for --quality-mode qvz.
    #[arg(long)]
    qvz_ratio: Option<f64>,
    /// Print detailed component sizes and mapping statistics.
    #[arg(long)]
    profile: bool,
    /// Replace an existing output file.
    #[arg(long)]
    force: bool,
    /// Directory in which to create temporary compression files.
    #[arg(long)]
    temp_dir: Option<PathBuf>,
    /// Validate the first N decoded records against the originals after
    /// compression finishes. 0 disables the check.
    #[arg(long, default_value_t = 1_000_000)]
    post_validate_records: u64,
}

#[derive(Clone, Copy, ValueEnum)]
enum OutputFormat {
    Fastq,
    Fasta,
}

impl From<OutputFormat> for Format {
    fn from(value: OutputFormat) -> Self {
        match value {
            OutputFormat::Fastq => Format::Fastq,
            OutputFormat::Fasta => Format::Fasta,
        }
    }
}

#[derive(Args)]
struct DecompressArgs {
    /// Input aligned .seqarc archive.
    archive: PathBuf,
    /// Primary output; mate 1 for paired archives.
    #[arg(short = '1', long)]
    mate_1: PathBuf,
    /// Optional mate-2 output. Required for paired archives.
    #[arg(short = '2', long)]
    mate_2: Option<PathBuf>,
    /// Optional index-read output for I1. Required when the archive stores index reads.
    #[arg(long = "index-output-1")]
    index_output_1: Option<PathBuf>,
    /// Optional index-read output for I2. Requires --index-output-1.
    #[arg(long = "index-output-2", requires = "index_output_1")]
    index_output_2: Option<PathBuf>,
    #[arg(long, value_enum, default_value_t = OutputFormat::Fastq)]
    format: OutputFormat,
    /// Force gzip output regardless of output suffix.
    #[arg(long, conflicts_with = "plain")]
    gzip: bool,
    /// Force plain output regardless of output suffix.
    #[arg(long, conflicts_with = "gzip")]
    plain: bool,
    #[arg(long, default_value_t = 6, value_parser = clap::value_parser!(u32).range(0..=9))]
    gzip_level: u32,
    /// Inclusive zero-based logical-record range as START:END.
    #[arg(long, value_parser = parse_range)]
    range: Option<RangeInclusive<u64>>,
    /// Number of aligned chunk decoder threads.
    #[arg(long, default_value_t = 1)]
    threads: usize,
    /// Replace existing output files.
    #[arg(long)]
    force: bool,
    /// Skip chunks that fail their checksum instead of aborting, recovering
    /// records from all undamaged chunks.
    #[arg(long)]
    best_effort: bool,
}

#[derive(Args)]
struct VerifyArgs {
    /// Input aligned .seqarc archive.
    archive: PathBuf,
}

#[derive(Args)]
struct InfoArgs {
    /// Input aligned .seqarc archive.
    archive: PathBuf,
    /// Break down on-disk storage by archive section and record subsection.
    #[arg(long)]
    breakdown: bool,
}

#[derive(Args)]
struct ValidateArgs {
    /// Input aligned .seqarc archive.
    archive: PathBuf,
    /// Original plain, gzip, or zstd FASTQ; mate 1 for paired comparison.
    #[arg(short = '1', long)]
    mate_1: Option<PathBuf>,
    /// Original plain, gzip, or zstd mate-2 FASTQ for byte comparison.
    #[arg(short = '2', long, requires = "mate_1")]
    mate_2: Option<PathBuf>,
    /// Optional original index-read FASTQ (I1) for byte comparison.
    #[arg(long = "index-read-1")]
    index_read_1: Option<PathBuf>,
    /// Optional original index-read FASTQ (I2) for byte comparison. Requires --index-read-1.
    #[arg(long = "index-read-2", requires = "index_read_1")]
    index_read_2: Option<PathBuf>,
    /// Apply the compression-time identifier prefix before comparison.
    #[arg(long, requires = "mate_1")]
    id_prefix_chars: Option<usize>,
    /// Apply the compression-time identifier field selection before comparison.
    #[arg(long, requires = "mate_1")]
    id_field: Option<usize>,
    /// Stop after validating this many records instead of the whole archive.
    #[arg(long, value_parser = clap::value_parser!(u64).range(1..))]
    records: Option<u64>,
}

fn parse_range(value: &str) -> std::result::Result<RangeInclusive<u64>, String> {
    let (start, end) = value
        .split_once(':')
        .ok_or_else(|| "range must be START:END".to_owned())?;
    let start = start
        .parse::<u64>()
        .map_err(|_| "invalid range start".to_owned())?;
    let end = end
        .parse::<u64>()
        .map_err(|_| "invalid range end".to_owned())?;
    if start > end {
        return Err("range start must not exceed range end".into());
    }
    Ok(start..=end)
}

fn parse_index_range(value: &str) -> std::result::Result<EmbeddedIndexRange, String> {
    let mut fields = value.split(':');
    let mate = match fields.next().map(str::to_ascii_uppercase).as_deref() {
        Some("R1") => EmbeddedIndexMate::R1,
        Some("R2") => EmbeddedIndexMate::R2,
        _ => return Err("index range mate must be R1 or R2".into()),
    };
    let start = fields
        .next()
        .ok_or_else(|| "index range must be MATE:START:END".to_owned())?
        .parse::<usize>()
        .map_err(|_| "invalid index range start".to_owned())?;
    let end = fields
        .next()
        .ok_or_else(|| "index range must be MATE:START:END".to_owned())?
        .parse::<usize>()
        .map_err(|_| "invalid index range end".to_owned())?;
    if fields.next().is_some() {
        return Err("index range must be MATE:START:END".into());
    }
    if start == 0 || end < start {
        return Err("index range uses one-based inclusive coordinates with START <= END".into());
    }
    Ok(EmbeddedIndexRange {
        mate,
        start: start - 1,
        end,
    })
}

fn validate_index_ranges(ranges: &[EmbeddedIndexRange], paired: bool) -> Result<()> {
    if !paired
        && ranges
            .iter()
            .any(|range| range.mate == EmbeddedIndexMate::R2)
    {
        return Err(Error::InvalidInput(
            "R2 index range requires paired-end input".into(),
        ));
    }
    for mate in [EmbeddedIndexMate::R1, EmbeddedIndexMate::R2] {
        let mut ranges = ranges
            .iter()
            .filter(|range| range.mate == mate)
            .collect::<Vec<_>>();
        ranges.sort_by_key(|range| range.start);
        for pair in ranges.windows(2) {
            if pair[1].start < pair[0].end {
                return Err(Error::InvalidInput(format!(
                    "overlapping index ranges for {mate:?}"
                )));
            }
        }
    }
    Ok(())
}

/// Attach the offending path to an I/O error so messages identify which of
/// the several input/output files failed.
fn io_path_error(path: &Path, error: std::io::Error) -> Error {
    Error::Io(std::io::Error::new(
        error.kind(),
        format!("{}: {error}", path.display()),
    ))
}

fn open_archive(path: &Path) -> Result<BufReader<File>> {
    let file = File::open(path).map_err(|error| io_path_error(path, error))?;
    Ok(BufReader::with_capacity(16 * 1024 * 1024, file))
}

/// Flush a renamed file's parent directory entry to stable storage so the
/// rename itself survives a crash.
fn sync_parent_directory(path: &Path) -> Result<()> {
    // A bare filename such as "out.seqarc-index" has an empty parent path;
    // opening it would fail with ENOENT, so fall back to "." as well.
    let parent = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    let directory = File::open(parent).map_err(|error| io_path_error(parent, error))?;
    directory.sync_all()?;
    Ok(())
}

struct AtomicPath {
    final_path: PathBuf,
    partial_path: PathBuf,
    committed: bool,
}

impl AtomicPath {
    fn new(final_path: PathBuf, force: bool) -> Result<Self> {
        if final_path.exists() && !force {
            return Err(Error::InvalidInput(format!(
                "output already exists: {}",
                final_path.display()
            )));
        }
        let file_name = final_path
            .file_name()
            .ok_or_else(|| Error::InvalidInput("output path has no file name".into()))?
            .to_string_lossy();
        let partial_path =
            final_path.with_file_name(format!(".{file_name}.partial.{}", std::process::id()));
        let _ = std::fs::remove_file(&partial_path);
        Ok(Self {
            final_path,
            partial_path,
            committed: false,
        })
    }

    fn partial(&self) -> &Path {
        &self.partial_path
    }

    fn commit(mut self) -> Result<()> {
        // Flush the partial file to stable storage before the rename so a
        // crash never leaves a zero-length or truncated file at the final
        // path. POSIX rename atomically replaces an existing destination, so
        // no explicit remove is needed; removing first would open a crash
        // window in which both the old and new archives are lost.
        let partial = File::open(&self.partial_path).map_err(|error| {
            io_path_error(&self.partial_path, error)
        })?;
        partial.sync_all()?;
        drop(partial);
        std::fs::rename(&self.partial_path, &self.final_path)?;
        self.committed = true;
        sync_parent_directory(&self.final_path)
    }
}

impl Drop for AtomicPath {
    fn drop(&mut self) {
        if !self.committed {
            let _ = std::fs::remove_file(&self.partial_path);
        }
    }
}

enum RecordWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl Write for RecordWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(bytes),
            Self::Gzip(writer) => writer.write(bytes),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

impl RecordWriter {
    fn create(path: &Path, gzip: bool, level: u32) -> Result<Self> {
        let writer = BufWriter::with_capacity(
            4 * 1024 * 1024,
            File::create(path).map_err(|error| io_path_error(path, error))?,
        );
        Ok(if gzip {
            Self::Gzip(seqarc::gzip_writer(writer, level))
        } else {
            Self::Plain(writer)
        })
    }

    fn finish(self) -> Result<()> {
        match self {
            Self::Plain(mut writer) => writer.flush().map_err(Into::into),
            Self::Gzip(writer) => {
                let mut writer = writer.finish()?;
                writer.flush().map_err(Into::into)
            }
        }
    }
}

struct OutputSink {
    first: Option<RecordWriter>,
    second: Option<RecordWriter>,
    index: Vec<Option<RecordWriter>>,
    format: Format,
}

impl AlignedRecordSink for OutputSink {
    fn write_group(&mut self, _original_id: u64, group: AlignedRecordGroup) -> Result<()> {
        match group.sequence {
            SequenceRecordGroup::Single(record) => {
                if self.second.is_some() {
                    return Err(Error::InvalidInput(
                        "single-end archive received two output paths".into(),
                    ));
                }
                write_record(self.first.as_mut().unwrap(), &record, self.format)?;
            }
            SequenceRecordGroup::Paired(first, second) => {
                let second_writer = self.second.as_mut().ok_or_else(|| {
                    Error::InvalidInput("paired archive requires a second output path".into())
                })?;
                write_record(self.first.as_mut().unwrap(), &first, self.format)?;
                write_record(second_writer, &second, self.format)?;
            }
        }
        if group.index_reads.len() != self.index.len() {
            return Err(Error::InvalidInput(format!(
                "archive contains {} index reads per record but {} index output paths were provided",
                group.index_reads.len(),
                self.index.len()
            )));
        }
        for (writer, index_record) in self.index.iter_mut().zip(&group.index_reads) {
            write_record(
                writer.as_mut().expect("index output writer is required"),
                index_record,
                self.format,
            )?;
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        self.first.take().unwrap().finish()?;
        if let Some(second) = self.second.take() {
            second.finish()?;
        }
        for writer in self.index.iter_mut() {
            writer
                .take()
                .ok_or_else(|| Error::InvalidInput("index output writer is missing".into()))?
                .finish()?;
        }
        Ok(())
    }
}

fn write_record(writer: &mut impl Write, record: &ReadRecord, format: Format) -> Result<()> {
    match format {
        Format::Fastq => {
            let quality = record.quality.as_ref().ok_or_else(|| {
                Error::InvalidInput("FASTQ output requires quality values".into())
            })?;
            if quality.len() != record.sequence.len() {
                return Err(Error::InvalidFormat("quality length mismatch".into()));
            }
            writer.write_all(b"@")?;
            writer.write_all(&record.id)?;
            writer.write_all(b"\n")?;
            writer.write_all(&record.sequence)?;
            writer.write_all(b"\n+\n")?;
            writer.write_all(quality)?;
            writer.write_all(b"\n")?;
        }
        Format::Fasta => {
            writer.write_all(b">")?;
            writer.write_all(&record.id)?;
            writer.write_all(b"\n")?;
            writer.write_all(&record.sequence)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

struct CountingSink {
    records: u64,
}

impl AlignedRecordSink for CountingSink {
    fn write_group(&mut self, _: u64, _: AlignedRecordGroup) -> Result<()> {
        self.records += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        Ok(())
    }
}

type OriginalReader = FastxReader<Box<dyn std::io::BufRead + Send>>;

/// Which fields an archival comparison should check. Archives written with
/// `--no-ids`, lossy quality modes, or embedded index ranges without original
/// index files cannot be compared byte-for-byte on every field.
#[derive(Clone, Copy, Debug)]
struct ComparisonMode {
    ids: bool,
    quality: bool,
    index: bool,
}

impl Default for ComparisonMode {
    fn default() -> Self {
        Self {
            ids: true,
            quality: true,
            index: true,
        }
    }
}

struct ComparisonSink {
    first: OriginalReader,
    second: Option<OriginalReader>,
    index: Vec<OriginalReader>,
    id_prefix_chars: Option<usize>,
    id_field: Option<usize>,
    mode: ComparisonMode,
    /// Record limit applied to the decode, if any. `finish` skips the
    /// leftover-original check once this many records have been compared,
    /// because the remaining originals were intentionally not consumed.
    limit: Option<u64>,
    records: u64,
}

impl ComparisonSink {
    fn expected(
        reader: &mut OriginalReader,
        id_prefix_chars: Option<usize>,
        id_field: Option<usize>,
    ) -> Result<ReadRecord> {
        let mut record = reader
            .next()
            .ok_or_else(|| Error::InvalidInput("archive has extra records".into()))??;
        if let Some(id_field) = id_field {
            select_id_field(&mut record, id_field)?;
        }
        if let Some(id_prefix_chars) = id_prefix_chars {
            truncate_read_id(&mut record, id_prefix_chars)?;
        }
        Ok(record)
    }

    fn compare(
        actual: &ReadRecord,
        expected: ReadRecord,
        label: &str,
        original_id: u64,
        mode: ComparisonMode,
    ) -> Result<()> {
        if mode.ids && actual.id != expected.id {
            return Err(Error::InvalidInput(format!(
                "mismatch at {label} {original_id}: identifier differs"
            )));
        }
        if actual.sequence != expected.sequence {
            let first_difference = actual
                .sequence
                .iter()
                .zip(&expected.sequence)
                .position(|(actual, expected)| actual != expected)
                .unwrap_or(actual.sequence.len().min(expected.sequence.len()));
            return Err(Error::InvalidInput(format!(
                "mismatch at {label} {original_id}: sequence differs at position {first_difference} (decoded length {}, original length {})",
                actual.sequence.len(),
                expected.sequence.len()
            )));
        }
        if mode.quality && actual.quality != expected.quality {
            return Err(Error::InvalidInput(format!(
                "mismatch at {label} {original_id}: quality differs"
            )));
        }
        Ok(())
    }
}

impl AlignedRecordSink for ComparisonSink {
    fn write_group(&mut self, original_id: u64, group: AlignedRecordGroup) -> Result<()> {
        match group.sequence {
            SequenceRecordGroup::Single(record) => {
                if self.second.is_some() {
                    return Err(Error::InvalidInput(
                        "single-end archive received two original files".into(),
                    ));
                }
                Self::compare(
                    &record,
                    Self::expected(&mut self.first, self.id_prefix_chars, self.id_field)?,
                    "record",
                    original_id,
                    self.mode,
                )?;
            }
            SequenceRecordGroup::Paired(first, second) => {
                let second_reader = self.second.as_mut().ok_or_else(|| {
                    Error::InvalidInput("paired archive requires two original files".into())
                })?;
                Self::compare(
                    &first,
                    Self::expected(&mut self.first, self.id_prefix_chars, self.id_field)?,
                    "pair mate 1",
                    original_id,
                    self.mode,
                )?;
                Self::compare(
                    &second,
                    Self::expected(second_reader, self.id_prefix_chars, self.id_field)?,
                    "pair mate 2",
                    original_id,
                    self.mode,
                )?;
            }
        }
        if self.mode.index {
            if group.index_reads.len() != self.index.len() {
                return Err(Error::InvalidInput(format!(
                    "archive contains {} index reads per record but {} index originals were provided",
                    group.index_reads.len(),
                    self.index.len()
                )));
            }
            for (reader, actual) in self.index.iter_mut().zip(&group.index_reads) {
                Self::compare(
                    actual,
                    Self::expected(reader, self.id_prefix_chars, self.id_field)?,
                    "index read",
                    original_id,
                    self.mode,
                )?;
            }
        }
        self.records += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        let truncated = self.limit.is_some_and(|limit| self.records >= limit);
        if truncated {
            return Ok(());
        }
        if self.first.next().is_some()
            || self
                .second
                .as_mut()
                .is_some_and(|reader| reader.next().is_some())
            || self.index.iter_mut().any(|reader| reader.next().is_some())
        {
            return Err(Error::InvalidInput(
                "original FASTQ contains extra records".into(),
            ));
        }
        Ok(())
    }
}

fn quality_mode(args: &CompressArgs) -> Result<QualityMode> {
    if args.no_quality && !matches!(args.quality_mode, QualityChoice::Lossless) {
        return Err(Error::InvalidInput(
            "--quality-mode cannot be used with --no-quality".into(),
        ));
    }
    let binary_values = [args.binary_threshold, args.binary_high, args.binary_low];
    if !matches!(args.quality_mode, QualityChoice::Binary)
        && binary_values.iter().any(Option::is_some)
    {
        return Err(Error::InvalidInput(
            "binary quality options require --quality-mode binary".into(),
        ));
    }
    if !matches!(args.quality_mode, QualityChoice::Qvz) && args.qvz_ratio.is_some() {
        return Err(Error::InvalidInput(
            "--qvz-ratio requires --quality-mode qvz".into(),
        ));
    }
    match args.quality_mode {
        QualityChoice::Lossless => Ok(QualityMode::Lossless),
        QualityChoice::Illumina => Ok(QualityMode::IlluminaBinning),
        QualityChoice::Qvz => {
            let ratio = args.qvz_ratio.ok_or_else(|| {
                Error::InvalidInput("qvz quality mode requires --qvz-ratio".into())
            })?;
            if ratio <= 0.0 || ratio.is_nan() {
                return Err(Error::InvalidInput(
                    "--qvz-ratio must be a positive number".into(),
                ));
            }
            Ok(QualityMode::Qvz { ratio })
        }
        QualityChoice::Binary => {
            let threshold = args.binary_threshold.ok_or_else(|| {
                Error::InvalidInput("binary quality mode requires --binary-threshold".into())
            })?;
            let high = args.binary_high.ok_or_else(|| {
                Error::InvalidInput("binary quality mode requires --binary-high".into())
            })?;
            let low = args.binary_low.ok_or_else(|| {
                Error::InvalidInput("binary quality mode requires --binary-low".into())
            })?;
            if threshold > 93 || high > 93 || low > 93 {
                return Err(Error::InvalidInput(
                    "binary quality values must be in the Phred range 0..=93".into(),
                ));
            }
            Ok(QualityMode::Binary {
                threshold,
                high,
                low,
            })
        }
    }
}

fn percentage(count: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        count as f64 * 100.0 / total as f64
    }
}

fn run(cli: Cli) -> Result<()> {
    match cli.command {
        Command::Index(args) => {
            if args.threads == 0 || args.sa_sparse_d == 0 {
                return Err(Error::InvalidInput(
                    "--threads and --sa-sparse-d must be nonzero".into(),
                ));
            }
            if matches!(args.aligner, AlignerChoice::Minibwa) && args.gtf.is_some() {
                return Err(Error::InvalidInput(
                    "--gtf is only supported with --aligner rustar".into(),
                ));
            }
            let output = AtomicPath::new(args.output, args.force)?;
            build_alignment_index(
                args.reference,
                output.partial(),
                &GenomeIndexBuildOptions {
                    temp_dir: args.temp_dir,
                    aligner: args.aligner.into(),
                    threads: args.threads,
                    rustar_gtf: args.gtf,
                    rustar_sjdb_overhang: args.sjdb_overhang,
                    rustar_sa_sparse_d: args.sa_sparse_d,
                    rustar_sa_index_nbases: args.sa_index_nbases,
                    rustar_chr_bin_nbits: args.chr_bin_nbits,
                },
            )?;
            output.commit()
        }
        Command::Compress(args) => {
            if args.threads == 0 || args.chunk_records == 0 {
                return Err(Error::InvalidInput(
                    "--threads and --chunk-records must be nonzero".into(),
                ));
            }
            // Fail fast on missing inputs instead of after loading the index.
            let mut inputs: Vec<&PathBuf> = vec![&args.index, &args.mate_1];
            inputs.extend(args.mate_2.iter());
            inputs.extend(args.index_read_1.iter());
            inputs.extend(args.index_read_2.iter());
            for input in inputs {
                if !input.exists() {
                    return Err(Error::InvalidInput(format!(
                        "input file not found: {}",
                        input.display()
                    )));
                }
            }
            validate_index_ranges(&args.index_ranges, args.mate_2.is_some())?;
            let quality = quality_mode(&args)?;
            let post_validation = (args.post_validate_records > 0).then(|| PostValidation {
                archive: args.output.clone(),
                mate_1: args.mate_1.clone(),
                mate_2: args.mate_2.clone(),
                index_read_1: args.index_read_1.clone(),
                index_read_2: args.index_read_2.clone(),
                id_prefix_chars: args.id_prefix_chars,
                id_field: args.id_field,
                mode: ComparisonMode {
                    ids: !args.no_ids,
                    quality: !args.no_quality && matches!(quality, QualityMode::Lossless),
                    index: args.index_read_1.is_some() || args.index_ranges.is_empty(),
                },
                limit: args.post_validate_records,
            });
            let options = AlignmentCompressionOptions {
                chunk_records: args.chunk_records,
                threads: args.threads,
                temp_dir: args.temp_dir.clone(),
                preserve_ids: !args.no_ids,
                id_prefix_chars: args.id_prefix_chars,
                id_field: args.id_field,
                preserve_quality: !args.no_quality,
                quality_mode: quality,
                zstd_level: args.zstd_level,
                embedded_index_ranges: args.index_ranges.clone(),
                index_read_files: [args.index_read_1, args.index_read_2]
                    .into_iter()
                    .flatten()
                    .collect(),
                ..AlignmentCompressionOptions::default()
            };
            let output = AtomicPath::new(args.output, args.force)?;
            let (metadata, profile) = if let Some(mate_2) = args.mate_2 {
                compress_aligned_paired_fastq_files_with_index_profiled(
                    args.index,
                    args.mate_1,
                    mate_2,
                    output.partial(),
                    &options,
                )?
            } else {
                compress_aligned_single_fastq_file_with_index_profiled(
                    args.index,
                    args.mate_1,
                    output.partial(),
                    &options,
                )?
            };
            output.commit()?;
            eprintln!(
                "compressed {} records ({})",
                metadata.records,
                if metadata.streams == 1 {
                    "single-end"
                } else {
                    "paired-end"
                }
            );
            assert!(profile.rejected_alignments <= profile.raw_reads);
            let total_reads = profile.mapped_reads + profile.raw_reads;
            let no_hit_reads = profile
                .raw_reads
                .saturating_sub(profile.rejected_alignments);
            eprintln!(
                "alignment: {:.2}% aligned, {:.2}% unaligned ({:.2}% no hit, {:.2}% not compressible)",
                percentage(profile.mapped_reads, total_reads),
                percentage(profile.raw_reads, total_reads),
                percentage(no_hit_reads, total_reads),
                percentage(profile.rejected_alignments, total_reads)
            );
            if args.profile {
                eprintln!(
                    "timings: total {:?}, index {:?}, map/spool {:?}, reference {:?}, chunks {:?}",
                    profile.total,
                    profile.index,
                    profile.map_and_spool,
                    profile.reference,
                    profile.chunks
                );
                eprintln!("archive bytes: {}", profile.archive_bytes);
                eprintln!(
                    "reference: {} bases, {} bytes",
                    profile.reference_bases, profile.reference_compressed_bytes
                );
                eprintln!(
                    "compressed sections: ids {}, qualities {}, records {} bytes",
                    profile.ids_compressed_bytes,
                    profile.qualities_compressed_bytes,
                    profile.records_compressed_bytes
                );
                eprintln!(
                    "index reads: {} raw, {} compressed bytes",
                    profile.indexes_raw_bytes, profile.indexes_compressed_bytes
                );
                eprintln!(
                    "mapping: {} mapped, {} raw ({} no-hit, {} not-compressible), {} candidates, {} non-first",
                    profile.mapped_reads,
                    profile.raw_reads,
                    no_hit_reads,
                    profile.rejected_alignments,
                    profile.candidate_hits,
                    profile.selected_non_first_hits
                );
                eprintln!(
                    "pipeline: {} chunks, {} spool bytes, max {} pending chunks",
                    profile.chunk_count, profile.spool_bytes, profile.max_pending_chunks
                );
            }
            if let Some(post_validation) = post_validation {
                post_compress_validation(&post_validation)?;
            }
            Ok(())
        }
        Command::Decompress(args) => decompress_command(args),
        Command::Validate(args) => validate_command(args),
        Command::Verify(args) => verify_command(args),
        Command::Info(args) => info_command(args),
    }
}

fn decompress_command(args: DecompressArgs) -> Result<()> {
    if args.threads == 0 {
        return Err(Error::InvalidInput("--threads must be nonzero".into()));
    }
    let output_paths = [
        Some(&args.mate_1),
        args.mate_2.as_ref(),
        args.index_output_1.as_ref(),
        args.index_output_2.as_ref(),
    ]
    .into_iter()
    .flatten()
    .collect::<Vec<_>>();
    for (index, path) in output_paths.iter().enumerate() {
        for other in &output_paths[index + 1..] {
            if path == other {
                return Err(Error::InvalidInput("output paths must be different".into()));
            }
        }
    }
    let first_path = AtomicPath::new(args.mate_1, args.force)?;
    let second_path = args
        .mate_2
        .map(|path| AtomicPath::new(path, args.force))
        .transpose()?;
    let index_paths = [args.index_output_1, args.index_output_2]
        .into_iter()
        .flatten()
        .map(|path| AtomicPath::new(path, args.force))
        .collect::<Result<Vec<_>>>()?;
    let gzip_first = output_is_gzip(&first_path.final_path, args.gzip, args.plain);
    let gzip_second = second_path
        .as_ref()
        .map(|path| output_is_gzip(&path.final_path, args.gzip, args.plain));
    let gzip_index = index_paths
        .iter()
        .map(|path| output_is_gzip(&path.final_path, args.gzip, args.plain))
        .collect::<Vec<_>>();
    let mut index_writers = Vec::with_capacity(index_paths.len());
    for (path, gzip) in index_paths.iter().zip(&gzip_index) {
        index_writers.push(Some(RecordWriter::create(
            path.partial(),
            *gzip,
            args.gzip_level,
        )?));
    }
    let mut sink = OutputSink {
        first: Some(RecordWriter::create(
            first_path.partial(),
            gzip_first,
            args.gzip_level,
        )?),
        second: second_path
            .as_ref()
            .zip(gzip_second)
            .map(|(path, gzip)| RecordWriter::create(path.partial(), gzip, args.gzip_level))
            .transpose()?,
        index: index_writers,
        format: args.format.into(),
    };
    let options = DecompressionOptions {
        original_range: args.range,
        threads: args.threads,
        best_effort: args.best_effort,
        ..DecompressionOptions::default()
    };
    let mut records = decompress_aligned_iter(open_archive(&args.archive)?, &options)?;
    let index_streams = records.index_streams();
    if index_streams != sink.index.len() {
        return Err(Error::InvalidInput(format!(
            "archive declares {index_streams} index streams but {} index output paths were provided",
            sink.index.len()
        )));
    }
    let layout = records.layout();
    // Check the layout eagerly so empty (zero-record) archives are validated
    // too; the per-record sink checks never fire for them.
    match layout {
        AlignmentLayout::SingleEnd if sink.second.is_some() => {
            return Err(Error::InvalidInput(
                "single-end archive received two output paths".into(),
            ));
        }
        AlignmentLayout::PairedEnd if sink.second.is_none() => {
            return Err(Error::InvalidInput(
                "paired archive requires a second output path".into(),
            ));
        }
        _ => {}
    }
    let total_chunks = records.total_chunks();
    let progress = decompression_progress(total_chunks);
    if progress.is_hidden() {
        eprintln!("[seqarc] stage=decompress start total_chunks={total_chunks}");
    }
    let mut record_count = 0u64;
    while let Some(record) = records.next() {
        progress.set_position(records.decoded_chunks());
        let (original_id, group) = record?;
        sink.write_group(original_id, group)?;
        record_count += 1;
    }
    progress.set_position(records.decoded_chunks());
    sink.finish()?;
    if progress.is_hidden() {
        eprintln!(
            "[seqarc] stage=decompress done decoded_chunks={} total_chunks={}",
            records.decoded_chunks(),
            total_chunks
        );
    }
    progress.finish_and_clear();
    first_path.commit()?;
    if let Some(second_path) = second_path {
        second_path.commit()?;
    }
    for index_path in index_paths {
        index_path.commit()?;
    }
    if records.skipped_chunks() > 0 {
        eprintln!(
            "best-effort recovery: skipped {} corrupt chunks (about {} records lost)",
            records.skipped_chunks(),
            records.skipped_records()
        );
        return Err(Error::Checksum(format!(
            "recovered {} records from undamaged chunks; {} chunks failed their checksums",
            record_count,
            records.skipped_chunks()
        )));
    }
    eprintln!("decompressed {} records ({:?})", record_count, layout);
    Ok(())
}

fn verify_command(args: VerifyArgs) -> Result<()> {
    let mut archive = open_archive(&args.archive)?;
    let report = verify_aligned_archive(&mut archive)?;
    if report.is_ok() {
        println!(
            "verified {} chunks, {} records: all checksums match",
            report.scanned_chunks, report.records
        );
        return Ok(());
    }
    if report.truncated {
        eprintln!(
            "archive is truncated: {}/{} declared chunks scanned",
            report.scanned_chunks, report.declared_chunks
        );
    }
    if report.bad_chunks > 0 {
        eprintln!(
            "{} of {} chunks failed their checksums (first bad chunk: {})",
            report.bad_chunks,
            report.scanned_chunks,
            report
                .first_bad_chunk
                .map(|chunk| chunk.to_string())
                .unwrap_or_else(|| "unknown".into())
        );
    }
    if report.trailing_bytes {
        eprintln!("archive contains trailing bytes after the final chunk");
    }
    println!(
        "verified {} chunks, {} records: {} damaged",
        report.scanned_chunks, report.records, report.bad_chunks
    );
    Err(Error::Checksum(format!(
        "archive {} failed verification",
        args.archive.display()
    )))
}

fn info_command(args: InfoArgs) -> Result<()> {
    let mut archive = open_archive(&args.archive)?;
    let info = read_archive_info_with_breakdown(&mut archive, args.breakdown)?;
    let file_size = std::fs::metadata(&args.archive)?.len();
    println!("archive: {}", args.archive.display());
    println!(
        "format version: {}.{}",
        info.version_major, info.version_minor
    );
    println!("layout: {:?}", info.layout);
    println!("records: {}", info.records);
    println!(
        "chunks: {}/{} scanned",
        info.scanned_chunks, info.declared_chunks
    );
    println!(
        "ids: {} (mode: {})",
        if info.preserve_ids {
            "preserved"
        } else {
            "not preserved"
        },
        info.id_mode
    );
    println!(
        "quality: {}",
        if info.preserve_quality {
            "preserved"
        } else {
            "not preserved"
        }
    );
    println!("index streams: {}", info.index_streams);
    for range in &info.embedded_indexes {
        println!(
            "embedded index: {:?}:{}:{}",
            range.mate,
            range.start + 1,
            range.end
        );
    }
    println!(
        "position widths: absolute={} delta={}",
        info.absolute_position_width, info.delta_position_width
    );
    println!(
        "reference block: {} bytes (compressed)",
        info.reference_block_bytes
    );
    println!("index metadata: {} bytes", info.index_metadata_bytes);
    println!("file size: {} bytes", file_size);
    if let Some(breakdown) = &info.breakdown {
        print_storage_breakdown(&info, breakdown, file_size);
    }
    println!("note: chunk checksums are not checked; use `seqarc verify` for that");
    if info.truncated {
        eprintln!(
            "warning: archive is truncated: {}/{} declared chunks present",
            info.scanned_chunks, info.declared_chunks
        );
    }
    if info.trailing_bytes {
        eprintln!("warning: archive contains trailing bytes after the final chunk");
    }
    Ok(())
}

/// Print the per-section storage table produced by `info --breakdown`.
///
/// Percentages are relative to the archive file size; sections storing
/// compressed data also show their decompressed size and ratio. The
/// `accounted` line sums every on-disk contribution, so for a complete
/// archive it equals the file size.
fn print_storage_breakdown(info: &ArchiveInfo, breakdown: &SectionBreakdown, file_size: u64) {
    let rows = vec![
        ("preamble:".to_owned(), breakdown.preamble, None),
        ("chunk headers:".to_owned(), breakdown.chunk_headers, None),
        (
            "payload directories:".to_owned(),
            breakdown.payload_directories,
            None,
        ),
        (
            "identifiers:".to_owned(),
            breakdown.ids.compressed,
            Some(breakdown.ids.raw),
        ),
        (
            "qualities:".to_owned(),
            breakdown.qualities.compressed,
            Some(breakdown.qualities.raw),
        ),
        (
            "records:".to_owned(),
            breakdown.records.compressed,
            Some(breakdown.records.raw),
        ),
        (
            "reference block:".to_owned(),
            breakdown.reference_block,
            None,
        ),
        ("index metadata:".to_owned(), breakdown.index_metadata, None),
    ];
    let subsections = subsection_rows(info, breakdown);
    let width = rows
        .iter()
        .map(|(label, _, _)| label.len())
        .chain(subsections.iter().map(|(label, _, _)| label.len()))
        .max()
        .unwrap_or(0);
    println!("storage breakdown:");
    for (label, compressed, raw) in &rows {
        println!(
            "  {label:<width$} {:>13} {}{}",
            display_bytes(*compressed),
            percent_of(*compressed, file_size),
            raw.map(|value| raw_note(value, *compressed))
                .unwrap_or_default()
        );
    }
    for (label, compressed, raw) in &subsections {
        println!(
            "    {label:<width$} {:>13} {}{}",
            display_bytes(*compressed),
            percent_of(*compressed, file_size),
            raw_note(*raw, *compressed)
        );
    }
    println!(
        "  {:<width$} {} of {} bytes",
        "accounted:",
        display_bytes(breakdown.accounted_bytes()),
        display_bytes(file_size)
    );
}

/// Flatten the record subsection breakdown into labeled rows in on-disk
/// order: flags, the per-mate subsections, one row per embedded index range,
/// and the index-read stream subsections.
fn subsection_rows(info: &ArchiveInfo, breakdown: &SectionBreakdown) -> Vec<(String, u64, u64)> {
    let subsections = &breakdown.subsections;
    let mut rows = vec![(
        "flags:".to_owned(),
        subsections.flags.compressed,
        subsections.flags.raw,
    )];
    let names = [
        "lengths",
        "positions",
        "noise positions",
        "noise",
        "reference skips",
        "emits",
        "reverses",
        "unaligned sequences",
    ];
    let paired = info.layout == AlignmentLayout::PairedEnd;
    for (mate, entry) in subsections.mates.iter().enumerate() {
        let prefix = if paired {
            format!("R{} ", mate + 1)
        } else {
            String::new()
        };
        let sizes = [
            entry.lengths,
            entry.positions,
            entry.noise_positions,
            entry.noise,
            entry.reference_skips,
            entry.emits,
            entry.reverses,
            entry.unaligned_sequences,
        ];
        for (name, size) in names.iter().zip(sizes) {
            rows.push((format!("{prefix}{name}:"), size.compressed, size.raw));
        }
    }
    for (range, sizes) in info
        .embedded_indexes
        .iter()
        .zip(&subsections.embedded_index_ranges)
    {
        rows.push((
            format!(
                "index range {:?}:{}:{}:",
                range.mate,
                range.start + 1,
                range.end
            ),
            sizes.compressed,
            sizes.raw,
        ));
    }
    for (stream, entry) in subsections.index_streams.iter().enumerate() {
        for (name, sizes) in [
            ("qualities", entry.qualities),
            ("lengths", entry.lengths),
            ("sequences", entry.sequences),
        ] {
            rows.push((
                format!("I{} {name}:", stream + 1),
                sizes.compressed,
                sizes.raw,
            ));
        }
    }
    rows
}

/// Format a byte count with thousands separators.
fn display_bytes(bytes: u64) -> String {
    let digits = bytes.to_string();
    let mut grouped = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, character) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            grouped.push(',');
        }
        grouped.push(character);
    }
    grouped
}

/// Format a byte count as a percentage of `total`.
fn percent_of(bytes: u64, total: u64) -> String {
    let percent = 100.0 * bytes as f64 / total.max(1) as f64;
    format!("({percent:5.1}%)")
}

/// Suffix showing the decompressed size and compression ratio, if any.
fn raw_note(raw: u64, compressed: u64) -> String {
    if raw == 0 {
        return String::new();
    }
    let ratio = if compressed > 0 {
        format!(" ({:.1}x)", raw as f64 / compressed as f64)
    } else {
        String::new()
    };
    format!("  raw {} bytes{}", display_bytes(raw), ratio)
}

#[allow(clippy::too_many_arguments)]
/// Everything needed to re-check a freshly written archive against its
/// originals. Captured before compression consumes the input paths.
struct PostValidation {
    archive: PathBuf,
    mate_1: PathBuf,
    mate_2: Option<PathBuf>,
    index_read_1: Option<PathBuf>,
    index_read_2: Option<PathBuf>,
    id_prefix_chars: Option<usize>,
    id_field: Option<usize>,
    mode: ComparisonMode,
    limit: u64,
}

fn post_compress_validation(validation: &PostValidation) -> Result<()> {
    // Validation re-opens the original input paths. Streams such as pipes,
    // /dev/stdin, or process substitution cannot be re-opened (the second
    // read yields EOF immediately), so skip the check for them.
    let reopenable = |path: &Path| {
        std::fs::metadata(path)
            .map(|metadata| metadata.is_file())
            .unwrap_or(false)
    };
    if !reopenable(&validation.mate_1)
        || validation
            .mate_2
            .as_deref()
            .is_some_and(|path| !reopenable(path))
        || [
            validation.index_read_1.as_deref(),
            validation.index_read_2.as_deref(),
        ]
        .into_iter()
        .flatten()
        .any(|path| !reopenable(path))
    {
        eprintln!("[seqarc] stage=post_validate skipped input_not_reopenable");
        return Ok(());
    }
    let result = (|| -> Result<()> {
        let started = std::time::Instant::now();
        let archive = open_archive(&validation.archive)?;
        let options = DecompressionOptions {
            max_records: Some(validation.limit),
            ..DecompressionOptions::default()
        };
        let index_files: Vec<PathBuf> = [
            validation.index_read_1.clone(),
            validation.index_read_2.clone(),
        ]
        .into_iter()
        .flatten()
        .collect();
        let mut sink = comparison_sink(
            &validation.mate_1,
            validation.mate_2.as_deref(),
            &index_files,
            validation.id_prefix_chars,
            validation.id_field,
            validation.mode,
            Some(validation.limit),
        )?;
        decompress_aligned_into(archive, &options, &mut sink)?;
        eprintln!(
            "[seqarc] stage=post_validate done elapsed={:.2}s records={}",
            started.elapsed().as_secs_f64(),
            sink.records
        );
        Ok(())
    })();
    result.map_err(|error| {
        // A validation failure right after writing means the archive or the
        // comparison is wrong; remove the archive so a corrupt file cannot be
        // mistaken for a good one, then report.
        let _ = std::fs::remove_file(&validation.archive);
        Error::InvalidInput(format!(
            "post-compression validation failed (archive removed): {error}"
        ))
    })
}

fn comparison_sink(
    mate_1: &Path,
    mate_2: Option<&Path>,
    index_files: &[PathBuf],
    id_prefix_chars: Option<usize>,
    id_field: Option<usize>,
    mode: ComparisonMode,
    limit: Option<u64>,
) -> Result<ComparisonSink> {
    Ok(ComparisonSink {
        first: FastxReader::new(open_fastx_reader(mate_1)?, Format::Fastq),
        second: mate_2
            .map(|second| {
                open_fastx_reader(second).map(|reader| FastxReader::new(reader, Format::Fastq))
            })
            .transpose()?,
        index: index_files
            .iter()
            .map(|path| {
                open_fastx_reader(path).map(|reader| FastxReader::new(reader, Format::Fastq))
            })
            .collect::<Result<Vec<_>>>()?,
        id_prefix_chars,
        id_field,
        mode,
        limit,
        records: 0,
    })
}

fn output_is_gzip(path: &Path, force_gzip: bool, force_plain: bool) -> bool {
    force_gzip
        || (!force_plain
            && path
                .extension()
                .is_some_and(|extension| extension.eq_ignore_ascii_case("gz")))
}

fn decompression_progress(total_chunks: u64) -> ProgressBar {
    let progress = ProgressBar::new(total_chunks);
    progress.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} decoding chunks [{elapsed_precise}] {pos}/{len} [{wide_bar:.cyan/blue}] {percent}% ETA {eta}",
        )
        .expect("valid decompression progress template"),
    );
    progress.enable_steady_tick(std::time::Duration::from_millis(250));
    progress
}

fn validate_command(args: ValidateArgs) -> Result<()> {
    let archive = open_archive(&args.archive)?;
    let options = DecompressionOptions {
        max_records: args.records,
        ..DecompressionOptions::default()
    };
    let limit_note = args
        .records
        .map(|limit| format!(" (record limit {limit})"))
        .unwrap_or_default();
    let index_files: Vec<PathBuf> = [args.index_read_1, args.index_read_2]
        .into_iter()
        .flatten()
        .collect();
    if let Some(first) = args.mate_1 {
        let mut sink = comparison_sink(
            &first,
            args.mate_2.as_deref(),
            &index_files,
            args.id_prefix_chars,
            args.id_field,
            ComparisonMode::default(),
            args.records,
        )?;
        let metadata = decompress_aligned_into(archive, &options, &mut sink)?;
        println!(
            "validated {} records{limit_note} against originals ({:?})",
            sink.records, metadata.layout
        );
    } else {
        let mut sink = CountingSink { records: 0 };
        let metadata = decompress_aligned_into(archive, &options, &mut sink)?;
        println!(
            "validated {} records{limit_note} ({:?})",
            sink.records, metadata.layout
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("seqarc: {error}");
            ExitCode::FAILURE
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_inclusive_range() {
        assert_eq!(parse_range("4:9").unwrap(), 4..=9);
        assert!(parse_range("9:4").is_err());
    }

    #[test]
    fn parses_one_based_embedded_index_ranges() {
        assert_eq!(
            parse_index_range("R2:9:16").unwrap(),
            EmbeddedIndexRange {
                mate: EmbeddedIndexMate::R2,
                start: 8,
                end: 16,
            }
        );
        assert!(parse_index_range("R1:0:8").is_err());
        assert!(parse_index_range("R3:1:8").is_err());
    }

    #[test]
    fn validates_embedded_index_range_layout() {
        let range = |mate, start, end| EmbeddedIndexRange { mate, start, end };
        assert!(validate_index_ranges(
            &[
                range(EmbeddedIndexMate::R1, 0, 4),
                range(EmbeddedIndexMate::R1, 4, 8),
            ],
            false,
        )
        .is_ok());
        assert!(validate_index_ranges(
            &[
                range(EmbeddedIndexMate::R1, 0, 4),
                range(EmbeddedIndexMate::R1, 3, 8),
            ],
            false,
        )
        .is_err());
        assert!(validate_index_ranges(&[range(EmbeddedIndexMate::R2, 0, 4)], false,).is_err());
    }

    #[test]
    fn percentage_handles_counts_and_empty_inputs() {
        assert_eq!(percentage(1, 4), 25.0);
        assert_eq!(percentage(0, 0), 0.0);
    }

    #[test]
    fn compression_requires_binary_parameters() {
        let cli = Cli::try_parse_from([
            "seqarc",
            "compress",
            "--index",
            "genome.seqarc-index",
            "-1",
            "r1.fastq.gz",
            "-2",
            "r2.fastq.gz",
            "-o",
            "reads.seqarc",
            "--quality-mode",
            "binary",
        ])
        .unwrap();
        let Command::Compress(args) = cli.command else {
            panic!("wrong command")
        };
        assert!(quality_mode(&args).is_err());

        assert!(Cli::try_parse_from([
            "seqarc",
            "compress",
            "--index",
            "genome.seqarc-index",
            "-1",
            "r1.fastq.gz",
            "-o",
            "reads.seqarc",
            "--no-ids",
            "--id-prefix-chars",
            "8",
        ])
        .is_err());
    }
}
