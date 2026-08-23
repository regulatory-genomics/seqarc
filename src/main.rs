use clap::{Args, Parser, Subcommand, ValueEnum};
use flate2::write::GzEncoder;
use indicatif::{ProgressBar, ProgressStyle};
use seqarc::{
    build_genome_index_with_options, compress_aligned_paired_fastq_files_with_index_profiled,
    compress_aligned_single_fastq_file_with_index_profiled, decompress_aligned_into,
    decompress_aligned_iter, open_fastx_reader, truncate_read_id, AlignedRecordGroup,
    AlignedRecordSink, AlignmentCompressionOptions, DecompressionOptions, Error, FastxReader,
    Format, GenomeIndexBuildOptions, QualityMode, ReadRecord, Result,
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
}

#[derive(Clone, Copy, ValueEnum)]
enum QualityChoice {
    Lossless,
    Illumina,
    Binary,
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
    /// Output .seqarc archive.
    #[arg(short, long)]
    output: PathBuf,
    #[arg(long, default_value_t = 16)]
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
    /// Print detailed component sizes and mapping statistics.
    #[arg(long)]
    profile: bool,
    /// Replace an existing output file.
    #[arg(long)]
    force: bool,
    /// Directory in which to create temporary compression files.
    #[arg(long)]
    temp_dir: Option<PathBuf>,
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
}

#[derive(Args)]
struct ValidateArgs {
    /// Input aligned .seqarc archive.
    archive: PathBuf,
    /// Original plain, gzip, or zstd FASTQ; mate 1 for paired comparison.
    #[arg(long)]
    mate_1: Option<PathBuf>,
    /// Original plain, gzip, or zstd mate-2 FASTQ for byte comparison.
    #[arg(long, requires = "mate_1")]
    mate_2: Option<PathBuf>,
    /// Apply the compression-time identifier prefix before comparison.
    #[arg(long, requires = "mate_1")]
    id_prefix_chars: Option<usize>,
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

struct AtomicPath {
    final_path: PathBuf,
    partial_path: PathBuf,
    force: bool,
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
            force,
            committed: false,
        })
    }

    fn partial(&self) -> &Path {
        &self.partial_path
    }

    fn commit(mut self) -> Result<()> {
        if self.force && self.final_path.exists() {
            std::fs::remove_file(&self.final_path)?;
        }
        std::fs::rename(&self.partial_path, &self.final_path)?;
        self.committed = true;
        Ok(())
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
        let writer = BufWriter::with_capacity(4 * 1024 * 1024, File::create(path)?);
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
    format: Format,
}

impl AlignedRecordSink for OutputSink {
    fn write_group(&mut self, _original_id: u64, group: AlignedRecordGroup) -> Result<()> {
        match group {
            AlignedRecordGroup::Single(record) => {
                if self.second.is_some() {
                    return Err(Error::InvalidInput(
                        "single-end archive received two output paths".into(),
                    ));
                }
                write_record(self.first.as_mut().unwrap(), &record, self.format)
            }
            AlignedRecordGroup::Paired(first, second) => {
                let second_writer = self.second.as_mut().ok_or_else(|| {
                    Error::InvalidInput("paired archive requires a second output path".into())
                })?;
                write_record(self.first.as_mut().unwrap(), &first, self.format)?;
                write_record(second_writer, &second, self.format)
            }
        }
    }

    fn finish(&mut self) -> Result<()> {
        self.first.take().unwrap().finish()?;
        if let Some(second) = self.second.take() {
            second.finish()?;
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

struct ComparisonSink {
    first: OriginalReader,
    second: Option<OriginalReader>,
    id_prefix_chars: Option<usize>,
    records: u64,
}

impl ComparisonSink {
    fn expected(reader: &mut OriginalReader, id_prefix_chars: Option<usize>) -> Result<ReadRecord> {
        let mut record = reader
            .next()
            .ok_or_else(|| Error::InvalidInput("archive has extra records".into()))??;
        if let Some(id_prefix_chars) = id_prefix_chars {
            truncate_read_id(&mut record, id_prefix_chars)?;
        }
        Ok(record)
    }
}

impl AlignedRecordSink for ComparisonSink {
    fn write_group(&mut self, original_id: u64, group: AlignedRecordGroup) -> Result<()> {
        match group {
            AlignedRecordGroup::Single(record) => {
                if self.second.is_some() {
                    return Err(Error::InvalidInput(
                        "single-end archive received two original files".into(),
                    ));
                }
                if record != Self::expected(&mut self.first, self.id_prefix_chars)? {
                    return Err(Error::InvalidInput(format!(
                        "mismatch at record {original_id}"
                    )));
                }
            }
            AlignedRecordGroup::Paired(first, second) => {
                let second_reader = self.second.as_mut().ok_or_else(|| {
                    Error::InvalidInput("paired archive requires two original files".into())
                })?;
                if first != Self::expected(&mut self.first, self.id_prefix_chars)? {
                    return Err(Error::InvalidInput(format!(
                        "mismatch at pair {original_id}, mate 1"
                    )));
                }
                if second != Self::expected(second_reader, self.id_prefix_chars)? {
                    return Err(Error::InvalidInput(format!(
                        "mismatch at pair {original_id}, mate 2"
                    )));
                }
            }
        }
        self.records += 1;
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if self.first.next().is_some()
            || self
                .second
                .as_mut()
                .is_some_and(|reader| reader.next().is_some())
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
    match args.quality_mode {
        QualityChoice::Lossless => Ok(QualityMode::Lossless),
        QualityChoice::Illumina => Ok(QualityMode::IlluminaBinning),
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
            let output = AtomicPath::new(args.output, args.force)?;
            build_genome_index_with_options(
                args.reference,
                output.partial(),
                &GenomeIndexBuildOptions {
                    temp_dir: args.temp_dir,
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
            let options = AlignmentCompressionOptions {
                chunk_records: args.chunk_records,
                threads: args.threads,
                temp_dir: args.temp_dir.clone(),
                preserve_ids: !args.no_ids,
                id_prefix_chars: args.id_prefix_chars,
                preserve_quality: !args.no_quality,
                quality_mode: quality_mode(&args)?,
                zstd_level: args.zstd_level,
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
            debug_assert!(profile.rejected_alignments <= profile.raw_reads);
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
            Ok(())
        }
        Command::Decompress(args) => decompress_command(args),
        Command::Validate(args) => validate_command(args),
    }
}

fn decompress_command(args: DecompressArgs) -> Result<()> {
    if args.threads == 0 {
        return Err(Error::InvalidInput("--threads must be nonzero".into()));
    }
    if args
        .mate_2
        .as_ref()
        .is_some_and(|path| path == &args.mate_1)
    {
        return Err(Error::InvalidInput(
            "mate output paths must be different".into(),
        ));
    }
    let first_path = AtomicPath::new(args.mate_1, args.force)?;
    let second_path = args
        .mate_2
        .map(|path| AtomicPath::new(path, args.force))
        .transpose()?;
    let gzip_first = output_is_gzip(&first_path.final_path, args.gzip, args.plain);
    let gzip_second = second_path
        .as_ref()
        .map(|path| output_is_gzip(&path.final_path, args.gzip, args.plain));
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
        format: args.format.into(),
    };
    let options = DecompressionOptions {
        original_range: args.range,
        threads: args.threads,
        ..DecompressionOptions::default()
    };
    let mut records = decompress_aligned_iter(
        BufReader::with_capacity(16 * 1024 * 1024, File::open(args.archive)?),
        &options,
    )?;
    let layout = records.layout();
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
    eprintln!("decompressed {} records ({:?})", record_count, layout);
    Ok(())
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
    let archive = BufReader::with_capacity(16 * 1024 * 1024, File::open(args.archive)?);
    if let Some(first) = args.mate_1 {
        let mut sink = ComparisonSink {
            first: FastxReader::new(open_fastx_reader(first)?, Format::Fastq),
            second: args
                .mate_2
                .map(|second| {
                    open_fastx_reader(second).map(|reader| FastxReader::new(reader, Format::Fastq))
                })
                .transpose()?,
            id_prefix_chars: args.id_prefix_chars,
            records: 0,
        };
        let metadata =
            decompress_aligned_into(archive, &DecompressionOptions::default(), &mut sink)?;
        println!(
            "validated {} records against originals ({:?})",
            sink.records, metadata.layout
        );
    } else {
        let mut sink = CountingSink { records: 0 };
        let metadata =
            decompress_aligned_into(archive, &DecompressionOptions::default(), &mut sink)?;
        println!("validated {} records ({:?})", sink.records, metadata.layout);
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
