use crate::minibwa::{self, CigarKind, CigarOp, Index, Strand};
use crate::{
    open_fastx_reader, temporary_directory, AlignmentLayout, Error, FastxReader, Format,
    ReadRecord, Result,
};
use indicatif::ProgressBar;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};

pub(crate) const MAP_TASK_READS: usize = 256;

/// Which read mapper an index or archive uses.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AlignerKind {
    /// BWA-style short-read mapper (default).
    #[default]
    Minibwa,
    /// Spliced RNA-seq mapper with GTF-driven junction support.
    Rustar,
}

/// Options for building a genome alignment index.
#[derive(Clone, Debug)]
pub struct GenomeIndexBuildOptions {
    /// Scratch directory for temporary files; defaults to the system temp dir.
    pub temp_dir: Option<PathBuf>,
    /// Mapper backend for the index.
    pub aligner: AlignerKind,
    /// Worker threads for index construction.
    pub threads: usize,
    /// GTF annotation for rustar splice-junction database construction.
    pub rustar_gtf: Option<PathBuf>,
    /// rustar `sjdbOverhang` for splice-junction database construction.
    pub rustar_sjdb_overhang: u32,
    /// rustar `--sa-sparseD` suffix-array sparsity factor.
    pub rustar_sa_sparse_d: u32,
    /// rustar `--sa-index-nbases` suffix-array seed length.
    pub rustar_sa_index_nbases: u32,
    /// rustar `--chr-bin-nbits` chromosome address width as a base-2 exponent.
    pub rustar_chr_bin_nbits: u32,
}

impl Default for GenomeIndexBuildOptions {
    fn default() -> Self {
        Self {
            temp_dir: None,
            aligner: AlignerKind::Minibwa,
            threads: 4,
            rustar_gtf: None,
            rustar_sjdb_overhang: 100,
            rustar_sa_sparse_d: 1,
            rustar_sa_index_nbases: 14,
            rustar_chr_bin_nbits: 18,
        }
    }
}

/// One selected alignment of a read against the reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GenomeAlignment {
    /// Zero-based reference offset of the first aligned base.
    pub start: i64,
    /// Zero-based read offset where the alignment starts.
    pub query_start: i32,
    /// Zero-based read offset just past the alignment end.
    pub query_end: i32,
    /// Alignment path as (read-length, operation) pairs.
    pub cigar: Vec<CigarOp>,
    /// True when the read maps to the reverse strand.
    pub is_reverse: bool,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct AlignmentSelection {
    pub(crate) candidate_hits: u32,
    pub(crate) selected_non_first: bool,
    pub(crate) rejected_alignment: bool,
}

pub(crate) struct SelectedAlignment {
    pub(crate) hit: Option<GenomeAlignment>,
    pub(crate) selection: AlignmentSelection,
}

struct AlignmentCandidate {
    alignment: GenomeAlignment,
    score: i32,
}

pub(crate) enum AlignmentIndex {
    Minibwa(Index),
    Rustar(Box<crate::rustar_index::RustarIndex>),
}

pub(crate) struct MappingProgressStats {
    layout: AlignmentLayout,
    pub(crate) processed: [u64; 2],
    pub(crate) aligned: [u64; 2],
}

impl MappingProgressStats {
    pub(crate) fn new(layout: AlignmentLayout) -> Self {
        Self {
            layout,
            processed: [0; 2],
            aligned: [0; 2],
        }
    }

    pub(crate) fn observe(
        &mut self,
        mate: usize,
        selected: &[SelectedAlignment],
        progress: &ProgressBar,
    ) {
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

    pub(crate) fn message(&self) -> String {
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

pub(crate) trait ReferenceBases: Sync {
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

impl AlignmentIndex {
    pub(crate) fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let mut magic = [0; 8];
        let is_rustar = File::open(path)?.read_exact(&mut magic).is_ok() && magic == *b"SRSTAR01";
        if is_rustar {
            return crate::rustar_index::load(path)
                .map(Box::new)
                .map(Self::Rustar);
        }
        Index::load(path)
            .map(Self::Minibwa)
            .map_err(|error| Error::InvalidInput(error.to_string()))
    }

    pub(crate) fn reference_len(&self) -> usize {
        match self {
            Self::Minibwa(index) => index.reference_len(),
            Self::Rustar(index) => index.reference_len(),
        }
    }
}

impl ReferenceBases for AlignmentIndex {
    fn len(&self) -> usize {
        self.reference_len()
    }

    #[inline(always)]
    fn base(&self, position: usize) -> u8 {
        match self {
            Self::Minibwa(index) => index.reference_base_at(position),
            Self::Rustar(index) => index.reference_base_at(position),
        }
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

pub(crate) struct ReferenceMask {
    bits: Vec<u64>,
    length: usize,
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct KeptRun {
    pub(crate) original_start: usize,
    pub(crate) compact_start: usize,
    pub(crate) length: usize,
}

/// Build a minibwa genome index with default options.
pub fn build_genome_index(
    genome_fasta: impl AsRef<Path>,
    index_file: impl AsRef<Path>,
) -> Result<Index> {
    build_genome_index_with_options(
        genome_fasta,
        index_file,
        &GenomeIndexBuildOptions::default(),
    )
}

/// Build a minibwa genome index with explicit options.
pub fn build_genome_index_with_options(
    genome_fasta: impl AsRef<Path>,
    index_file: impl AsRef<Path>,
    options: &GenomeIndexBuildOptions,
) -> Result<Index> {
    if options.aligner != AlignerKind::Minibwa {
        return Err(Error::InvalidInput(
            "build_genome_index only supports the minibwa backend".into(),
        ));
    }
    let workspace = temporary_directory(options.temp_dir.as_deref(), "seqarc-index-")?;
    let normalized = normalize_genome(genome_fasta.as_ref(), None, workspace.path())?;

    let index_file = index_file.as_ref();
    let threads = u32::try_from(options.threads.max(1))
        .map_err(|_| Error::InvalidInput("index thread count is too large".into()))?;
    Index::build_from_fasta_in(
        &normalized.fasta,
        index_file,
        threads,
        workspace.path(),
    )
    .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let index =
        Index::load(index_file).map_err(|error| Error::InvalidInput(error.to_string()))?;
    if index.reference_len() == 0 {
        return Err(Error::InvalidInput(
            "built genome index has an empty reference".into(),
        ));
    }
    Ok(index)
}

/// Build a genome index with the configured aligner backend.
pub fn build_alignment_index(
    genome_fasta: impl AsRef<Path>,
    index_file: impl AsRef<Path>,
    options: &GenomeIndexBuildOptions,
) -> Result<()> {
    match options.aligner {
        AlignerKind::Minibwa => {
            build_genome_index_with_options(genome_fasta, index_file, options).map(|_| ())
        }
        AlignerKind::Rustar => crate::rustar_index::build(genome_fasta, index_file, options),
    }
}

/// Normalize a FASTA reference into the single-line form the index builders expect.
pub(crate) struct NormalizedGenome {
    pub(crate) fasta: PathBuf,
    #[allow(dead_code)]
    pub(crate) gtf: Option<PathBuf>,
}

struct GtfFeature {
    start: u64,
    end: u64,
    mapped_start: Option<u64>,
    mapped_end: Option<u64>,
    line: String,
}

const GENOME_CONTIG_NAME: &str = "genome";

pub(crate) fn normalize_genome(
    genome_fasta: impl AsRef<Path>,
    gtf: Option<&Path>,
    workspace: &Path,
) -> Result<NormalizedGenome> {
    let mut gtf_features = match gtf {
        Some(path) => load_gtf(path)?,
        None => BTreeMap::new(),
    };

    let fasta_path = workspace.join("genome.clean.fa");
    let file = File::create(&fasta_path)?;
    let mut writer = std::io::BufWriter::new(file);
    writer.write_all(b">genome\n")?;
    let mut output = Vec::with_capacity(1 << 20);
    let mut line_length = 0usize;
    let mut reference_length = 0u64;
    let mut total_kept = 0u64;

    let reader = open_fastx_reader(genome_fasta)?;
    for record in FastxReader::new(reader, Format::Fasta) {
        let record = record?;
        let contig = (!gtf_features.is_empty()).then(|| contig_name(&record.id));
        let offset = total_kept;
        let events = contig
            .as_ref()
            .and_then(|contig| gtf_features.get(contig))
            .map(|features| {
                let mut events = Vec::with_capacity(features.len() * 2);
                for (index, feature) in features.iter().enumerate() {
                    events.push((feature.start, 0u8, index));
                    events.push((feature.end, 1u8, index));
                }
                events.sort_unstable();
                events
            });
        let mut kept = 0u64;
        let mut position = 0u64;
        let mut event_index = 0usize;
        for base in record.sequence {
            let base = base.to_ascii_uppercase();
            position += 1;
            if let Some(events) = events.as_ref() {
                while event_index < events.len()
                    && events[event_index].0 == position
                    && events[event_index].1 == 0
                {
                    let index = events[event_index].2;
                    gtf_features.get_mut(contig.as_ref().unwrap()).unwrap()[index].mapped_start =
                        Some(offset + kept + 1);
                    event_index += 1;
                }
            }
            if matches!(base, b'A' | b'C' | b'G' | b'T') {
                output.push(base);
                reference_length += 1;
                total_kept += 1;
                kept += 1;
                line_length += 1;
                if line_length == 80 {
                    output.push(b'\n');
                    line_length = 0;
                }
                if output.len() >= 1 << 20 {
                    writer.write_all(&output)?;
                    output.clear();
                }
            }
            if let Some(events) = events.as_ref() {
                while event_index < events.len()
                    && events[event_index].0 == position
                    && events[event_index].1 == 1
                {
                    let index = events[event_index].2;
                    gtf_features.get_mut(contig.as_ref().unwrap()).unwrap()[index].mapped_end =
                        Some(offset + kept);
                    event_index += 1;
                }
            }
        }
    }

    if reference_length == 0 {
        return Err(Error::InvalidInput(
            "genome FASTA does not contain any ACGT bases".into(),
        ));
    }
    if line_length != 0 {
        output.push(b'\n');
    }
    writer.write_all(&output)?;
    writer.flush()?;
    writer.get_ref().sync_all()?;

    let gtf_path = match gtf {
        Some(_) => {
            let path = workspace.join("genome.clean.gtf");
            write_gtf(&path, &gtf_features)?;
            Some(path)
        }
        None => None,
    };

    Ok(NormalizedGenome {
        fasta: fasta_path,
        gtf: gtf_path,
    })
}

fn contig_name(id: &[u8]) -> String {
    std::str::from_utf8(id)
        .ok()
        .and_then(|id| id.split_whitespace().next())
        .unwrap_or("")
        .to_owned()
}

fn load_gtf(path: &Path) -> Result<BTreeMap<String, Vec<GtfFeature>>> {
    let reader = open_fastx_reader(path)?;
    let mut features: BTreeMap<String, Vec<GtfFeature>> = BTreeMap::new();
    for (line_number, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.strip_suffix('\r').unwrap_or(&line);
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            return Err(Error::InvalidFormat(format!(
                "GTF line {} has fewer than 9 columns",
                line_number + 1
            )));
        }
        if fields[0].is_empty() {
            return Err(Error::InvalidFormat(format!(
                "GTF line {} has an empty seqname",
                line_number + 1
            )));
        }
        let start = fields[3].parse::<u64>().map_err(|_| {
            Error::InvalidFormat(format!(
                "GTF line {} has an invalid start coordinate",
                line_number + 1
            ))
        })?;
        let end = fields[4].parse::<u64>().map_err(|_| {
            Error::InvalidFormat(format!(
                "GTF line {} has an invalid end coordinate",
                line_number + 1
            ))
        })?;
        features
            .entry(fields[0].to_owned())
            .or_default()
            .push(GtfFeature {
                start,
                end,
                mapped_start: None,
                mapped_end: None,
                line: line.to_owned(),
            });
    }
    Ok(features)
}

fn write_gtf(path: &Path, features: &BTreeMap<String, Vec<GtfFeature>>) -> Result<()> {
    let mut writer = std::io::BufWriter::new(File::create(path)?);
    for features in features.values() {
        for feature in features {
            let (Some(start), Some(end)) = (feature.mapped_start, feature.mapped_end) else {
                continue;
            };
            if start > end {
                continue;
            }
            let mut fields = feature.line.split('\t');
            let source = fields.nth(1).unwrap_or("");
            let kind = fields.next().unwrap_or("");
            let score = fields.nth(2).unwrap_or("");
            let strand = fields.next().unwrap_or("");
            let frame = fields.next().unwrap_or("");
            let attributes = fields.next().unwrap_or("");
            write!(
                writer,
                "{GENOME_CONTIG_NAME}\t{source}\t{kind}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attributes}"
            )?;
            for extra in fields {
                write!(writer, "\t{extra}")?;
            }
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Align a batch of raw DNA sequences against a minibwa index.
///
/// The returned entries preserve input order. `None` indicates that a sequence
/// could not be mapped. When minibwa reports multiple hits, only its first hit
/// is returned.
pub fn align_dna_batch(
    index: &Index,
    sequences: &[impl AsRef<[u8]>],
) -> Result<Vec<Option<GenomeAlignment>>> {
    let names: Vec<String> = (0..sequences.len())
        .map(|index| format!("read-{index}"))
        .collect();
    let sequences = sequences.iter().map(AsRef::as_ref).collect::<Vec<&[u8]>>();
    let options = minibwa::Opts::new();
    let aligner = minibwa::Aligner::new(index, &options);
    let mut buffer = minibwa::ThreadBuf::new();
    let hits = aligner
        .map_many(&mut buffer, &names, &sequences)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;

    Ok(hits
        .into_iter()
        .map(|hits| hits.into_iter().next().map(GenomeAlignment::from))
        .collect())
}

impl From<minibwa::Hit> for GenomeAlignment {
    fn from(hit: minibwa::Hit) -> Self {
        Self {
            start: hit.ref_start,
            query_start: hit.query_start,
            query_end: hit.query_end,
            cigar: hit.cigar,
            is_reverse: hit.strand == Strand::Reverse,
        }
    }
}

pub(crate) fn map_records<F>(
    index: &AlignmentIndex,
    records: &[&ReadRecord],
    mapping_pool: &rayon::ThreadPool,
    candidate_cost: &F,
) -> Result<Vec<SelectedAlignment>>
where
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &AlignmentIndex) -> Result<usize> + Sync,
{
    // Reads fully consumed by embedded index ranges have empty biological
    // sequences, which aligners reject; treat them as unmapped instead.
    let mappable = records
        .iter()
        .enumerate()
        .filter(|(_, record)| !record.sequence.is_empty())
        .collect::<Vec<_>>();
    if mappable.len() == records.len() {
        return map_index_records(index, records, mapping_pool, candidate_cost);
    }
    let mapped = map_index_records(
        index,
        &mappable.iter().map(|(_, record)| **record).collect::<Vec<_>>(),
        mapping_pool,
        candidate_cost,
    )?;
    let mut selected: Vec<Option<SelectedAlignment>> =
        (0..records.len()).map(|_| None).collect();
    for ((position, _), alignment) in mappable.into_iter().zip(mapped) {
        selected[position] = Some(alignment);
    }
    Ok(selected
        .into_iter()
        .map(|selected| {
            selected.unwrap_or(SelectedAlignment {
                hit: None,
                selection: AlignmentSelection::default(),
            })
        })
        .collect())
}

fn map_index_records<F>(
    index: &AlignmentIndex,
    records: &[&ReadRecord],
    mapping_pool: &rayon::ThreadPool,
    candidate_cost: &F,
) -> Result<Vec<SelectedAlignment>>
where
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &AlignmentIndex) -> Result<usize> + Sync,
{
    match index {
        AlignmentIndex::Minibwa(raw_index) => {
            map_minibwa_records(raw_index, index, records, mapping_pool, candidate_cost)
        }
        AlignmentIndex::Rustar(raw_index) => {
            map_rustar_records(raw_index, index, records, mapping_pool, candidate_cost)
        }
    }
}

fn map_minibwa_records<F>(
    index: &Index,
    reference: &AlignmentIndex,
    records: &[&ReadRecord],
    mapping_pool: &rayon::ThreadPool,
    candidate_cost: &F,
) -> Result<Vec<SelectedAlignment>>
where
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &AlignmentIndex) -> Result<usize> + Sync,
{
    let opts = minibwa::Opts::new();
    let aligner = minibwa::Aligner::new(index, &opts);
    mapping_pool.install(|| {
        use rayon::prelude::*;
        records
            .par_chunks(MAP_TASK_READS)
            .map_init(minibwa::ThreadBuf::new, |buffer, batch| {
                align_batch(
                    &aligner,
                    buffer,
                    reference,
                    batch.iter().copied(),
                    candidate_cost,
                )
            })
            .collect::<Result<Vec<_>>>()
            .map(|batches| batches.into_iter().flatten().collect())
    })
}

fn align_batch<'a, I, F>(
    aligner: &minibwa::Aligner<'_>,
    buffer: &mut minibwa::ThreadBuf,
    reference: &AlignmentIndex,
    records: I,
    candidate_cost: &F,
) -> Result<Vec<SelectedAlignment>>
where
    I: IntoIterator<Item = &'a ReadRecord>,
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &AlignmentIndex) -> Result<usize>,
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
                .map(|(record, hits)| {
                    select_minibwa_alignment(record, hits, reference, candidate_cost)
                })
                .collect()
        })
}

fn map_rustar_records<F>(
    index: &crate::rustar_index::RustarIndex,
    reference: &AlignmentIndex,
    records: &[&ReadRecord],
    mapping_pool: &rayon::ThreadPool,
    candidate_cost: &F,
) -> Result<Vec<SelectedAlignment>>
where
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &AlignmentIndex) -> Result<usize> + Sync,
{
    mapping_pool.install(|| {
        use rayon::prelude::*;
        records
            .par_chunks(MAP_TASK_READS)
            .map(|batch| {
                batch
                    .iter()
                    .map(|record| {
                        let encoded = record
                            .sequence
                            .iter()
                            .copied()
                            .map(rustar_aligner::io::fastq::encode_base)
                            .collect::<Vec<_>>();
                        let name = rustar_read_name(&record.id);
                        let (transcripts, _, _, _) = rustar_aligner::align::align_read(
                            &encoded,
                            &name,
                            &index.inner,
                            &index.params,
                        )
                        .map_err(|error| Error::InvalidInput(error.to_string()))?;
                        let candidates = transcripts
                            .into_iter()
                            .map(|transcript| rustar_candidate(transcript, record.sequence.len()))
                            .collect::<Result<Vec<_>>>()?;
                        select_alignment(record, candidates, reference, candidate_cost)
                    })
                    .collect()
            })
            .collect::<Result<Vec<Vec<_>>>>()
            .map(|batches| batches.into_iter().flatten().collect())
    })
}

fn rustar_read_name(id: &[u8]) -> String {
    if let Ok(id) = std::str::from_utf8(id) {
        return id.to_owned();
    }
    let mut name = String::with_capacity(5 + id.len() * 2);
    name.push_str("read-");
    const HEX: &[u8; 16] = b"0123456789abcdef";
    for &byte in id {
        name.push(HEX[(byte >> 4) as usize] as char);
        name.push(HEX[(byte & 15) as usize] as char);
    }
    name
}

fn rustar_candidate(
    transcript: rustar_aligner::align::Transcript,
    read_len: usize,
) -> Result<AlignmentCandidate> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    if transcript.read_length() != read_len {
        return Err(Error::InvalidInput(
            "rustar CIGAR does not consume the complete read".into(),
        ));
    }
    let cigar = transcript
        .cigar
        .into_iter()
        .map(|op| {
            let kind = match op.kind() {
                Kind::Match => CigarKind::Match,
                Kind::Insertion => CigarKind::Ins,
                Kind::Deletion => CigarKind::Del,
                Kind::Skip => CigarKind::RefSkip,
                Kind::SoftClip => CigarKind::SoftClip,
                Kind::HardClip => CigarKind::HardClip,
                Kind::Pad => CigarKind::Pad,
                Kind::SequenceMatch => CigarKind::Eq,
                Kind::SequenceMismatch => CigarKind::Diff,
            };
            Ok(CigarOp {
                kind,
                len: u32::try_from(op.len()).map_err(|_| {
                    Error::InvalidInput("rustar CIGAR operation is too long".into())
                })?,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(AlignmentCandidate {
        score: transcript.score,
        alignment: GenomeAlignment {
            start: i64::try_from(transcript.genome_start)
                .map_err(|_| Error::InvalidInput("rustar alignment start is too large".into()))?,
            query_start: 0,
            query_end: i32::try_from(read_len)
                .map_err(|_| Error::InvalidInput("read is too long".into()))?,
            cigar,
            is_reverse: transcript.is_reverse,
        },
    })
}

pub(crate) fn select_minibwa_alignment<R, F>(
    record: &ReadRecord,
    hits: Vec<minibwa::Hit>,
    reference: &R,
    candidate_cost: &F,
) -> Result<SelectedAlignment>
where
    R: ReferenceBases + ?Sized,
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &R) -> Result<usize>,
{
    let candidates = hits
        .into_iter()
        .map(|hit| AlignmentCandidate {
            score: hit.score,
            alignment: GenomeAlignment::from(hit),
        })
        .collect();
    select_alignment(record, candidates, reference, candidate_cost)
}

fn select_alignment<R, F>(
    record: &ReadRecord,
    candidates: Vec<AlignmentCandidate>,
    reference: &R,
    candidate_cost: &F,
) -> Result<SelectedAlignment>
where
    R: ReferenceBases + ?Sized,
    F: Fn(&ReadRecord, Option<&GenomeAlignment>, &R) -> Result<usize>,
{
    let raw_cost = candidate_cost(record, None, reference)?;
    let candidate_hits = u32::try_from(candidates.len()).unwrap_or(u32::MAX);
    let mut best: Option<(usize, i32, usize, GenomeAlignment)> = None;
    for (index, candidate) in candidates.into_iter().enumerate() {
        let score = candidate.score;
        let alignment = candidate.alignment;
        let Ok(encoded_cost) = candidate_cost(record, Some(&alignment), reference) else {
            continue;
        };
        let candidate = (encoded_cost, score, index, alignment);
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

impl ReferenceMask {
    pub(crate) fn new(length: usize) -> Self {
        Self {
            bits: vec![0; length.div_ceil(64)],
            length,
        }
    }

    pub(crate) fn mark_hit(&mut self, hit: &GenomeAlignment) -> Result<()> {
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

    pub(crate) fn mark_range(&mut self, start: usize, end: usize) -> Result<()> {
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

pub(crate) fn compact_reference<R: ReferenceBases + ?Sized>(
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

pub(crate) fn compact_boundary(
    runs: &[KeptRun],
    original: usize,
    reference_len: usize,
) -> Result<usize> {
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

pub(crate) fn effective_cigar(hit: &GenomeAlignment, read_len: usize) -> Result<Vec<CigarOp>> {
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

pub(crate) fn normalize_cigar(cigar: impl IntoIterator<Item = CigarOp>) -> Result<Vec<CigarOp>> {
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

pub(crate) fn reverse_complement(sequence: &mut [u8]) {
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

fn nonnegative_usize(value: i64, field: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| Error::InvalidFormat(format!("negative {field}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_index_concatenates_chromosomes() {
        let directory = tempfile::tempdir().unwrap();
        let genome = directory.path().join("genome.fa");
        std::fs::write(&genome, b">chr1\nacnGT\n>chr2\ntgcaRY\n").unwrap();
        let index_file = directory.path().join("genome.seqarc-index");

        let index = build_genome_index(&genome, &index_file).unwrap();

        assert!(index_file.exists());
        assert_eq!(index.reference_len(), 8);
        assert_eq!(
            (0..index.reference_len())
                .map(|position| index.reference_base(position).unwrap())
                .collect::<Vec<_>>(),
            b"ACGTTGCA"
        );
    }

    #[test]
    fn dna_batch_alignment_returns_first_hit_and_strand() {
        let directory = tempfile::tempdir().unwrap();
        let genome = directory.path().join("genome.fa");
        let prefix = directory.path().join("genome.seqarc-index");
        let reference =
            b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
        std::fs::write(&genome, [b">chr1\n".as_slice(), reference, b"\n"].concat()).unwrap();
        let index = build_genome_index(&genome, &prefix).unwrap();

        let forward = reference[20..60].to_vec();
        let reverse = forward
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => *base,
            })
            .collect::<Vec<_>>();
        let alignments = align_dna_batch(
            &index,
            &[
                forward,
                reverse,
                b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN".to_vec(),
            ],
        )
        .unwrap();

        assert_eq!(alignments.len(), 3);
        assert_eq!(alignments[0].as_ref().unwrap().start, 20);
        assert!(!alignments[0].as_ref().unwrap().is_reverse);
        assert!(!alignments[0].as_ref().unwrap().cigar.is_empty());
        assert!(alignments[1].as_ref().unwrap().is_reverse);
        assert!(alignments[2].is_none());
    }

    #[test]
    fn normalize_genome_concatenates_and_remaps_gtf() {
        let directory = tempfile::tempdir().unwrap();
        let genome = directory.path().join("genome.fa");
        // chr1: A C N N G T -> kept A C G T (4)
        // chr2: t g c a N T -> kept T G C A T (5)
        std::fs::write(&genome, b">chr1\nACNNGT\n>chr2\ntgcaNT\n").unwrap();
        let gtf = directory.path().join("genes.gtf");
        std::fs::write(
            &gtf,
            b"chr1\ttest\texon\t2\t6\t.\t+\t.\tgene_id \"g1\";\n\
chr2\ttest\texon\t2\t5\t.\t-\t.\tgene_id \"g2\";\n\
chr1\ttest\texon\t3\t4\t.\t+\t.\tgene_id \"g3\";\n",
        )
        .unwrap();

        let normalized = normalize_genome(&genome, Some(&gtf), directory.path()).unwrap();

        let reference = std::fs::read_to_string(&normalized.fasta).unwrap();
        assert_eq!(reference, ">genome\nACGTTGCAT\n");
        let gtf_lines = std::fs::read_to_string(normalized.gtf.unwrap())
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect::<Vec<_>>();
        assert_eq!(
            gtf_lines,
            vec![
                "genome\ttest\texon\t2\t4\t.\t+\t.\tgene_id \"g1\";".to_owned(),
                "genome\ttest\texon\t6\t8\t.\t-\t.\tgene_id \"g2\";".to_owned(),
            ]
        );
    }
}
