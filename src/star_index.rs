use crate::{temporary_directory, Error, GenomeIndexBuildOptions, Result};
use memmap2::{Mmap, MmapOptions};
use crate::star::error::Error as StarError;
use crate::star::genome::{Genome, GenomeSeq};
use crate::star::index::SectionSink;
use crate::star::index::packed_array::PackedArray;
use crate::star::index::sa_index::SaIndex;
use crate::star::index::suffix_array::SuffixArray;
use crate::star::index::GenomeIndex;
use crate::star::junction::sjdb_insert::PreparedJunction;
use crate::star::junction::SpliceJunctionDb;
use crate::star::params::Parameters;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::Arc;

const INDEX_MAGIC: &[u8; 8] = b"SRSTAR01";
const INDEX_VERSION: u32 = 1;
const INDEX_HEADER_SIZE: u64 = 32;
const MMAP_ALIGNMENT: u64 = 64 * 1024;

#[derive(Clone, Debug)]
struct Section {
    offset: u64,
    length: u64,
}

pub(crate) struct StarIndex {
    pub(crate) inner: GenomeIndex,
    pub(crate) params: Parameters,
}

impl StarIndex {
    pub(crate) fn reference_len(&self) -> usize {
        self.inner.genome.n_genome_real as usize
    }

    #[inline(always)]
    pub(crate) fn reference_base_at(&self, position: usize) -> u8 {
        match self.inner.genome.sequence.base(position) {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        }
    }
}

pub(crate) fn build(
    genome_fasta: impl AsRef<Path>,
    output: impl AsRef<Path>,
    options: &GenomeIndexBuildOptions,
) -> Result<()> {
    if !(1..=14).contains(&options.star_sa_index_nbases) {
        return Err(Error::InvalidInput(
            "star SA index prefix length must be in 1..=14".into(),
        ));
    }
    if options.star_chr_bin_nbits > 30 {
        return Err(Error::InvalidInput(
            "star chromosome bin exponent must be at most 30".into(),
        ));
    }
    let workspace = temporary_directory(options.temp_dir.as_deref(), "seqarc-star-index-")?;
    let mut params = default_parameters()?;
    params.run_thread_n = NonZeroUsize::new(options.threads.max(1)).unwrap();
    // Single scratch space: caps-sa's external-memory chunks and the
    // normalized FASTA/GTF copies share the workspace (there is no
    // separate STAR-style index directory any more — sections stream
    // straight into the container).
    params.temp_dir = Some(workspace.path().to_path_buf());
    let normalized = crate::aligners::normalize_genome(
        genome_fasta.as_ref(),
        options.star_gtf.as_deref(),
        workspace.path(),
    )?;
    params.genome_fasta_files = vec![normalized.fasta];
    params.sjdb_gtf_file = normalized.gtf;
    params.genome_sa_sparse_d = options.star_sa_sparse_d.max(1);
    params.genome_sa_index_nbases = options.star_sa_index_nbases;
    params.genome_chr_bin_nbits = options.star_chr_bin_nbits;
    params.sjdb_overhang = options.star_sjdb_overhang;

    let output = output.as_ref();
    let mut sink = ContainerSink::create(output)?;
    GenomeIndex::generate_streaming_with(&params, &mut sink)
        .and_then(|()| sink.finish())
        .map_err(|error| {
            // Best effort: don't leave a truncated container behind.
            let _ = std::fs::remove_file(output);
            Error::InvalidInput(error.to_string())
        })?;
    Ok(())
}

/// Single-pass container writer: sections stream straight into the
/// `.seqarc-index` file. The table of contents is reserved (as a placeholder)
/// before the first section — the section list is known once the junction
/// preparation has decided whether `sjdbInfo.txt` exists — and back-patched
/// by [`ContainerSink::finish`] once every section's final offset/length is
/// known. Layout matches the old pack-the-directory path: header, TOC, then
/// sections in alphabetical order, each 64 KiB aligned.
struct ContainerSink {
    file: File,
    names: Vec<String>,
    entries: Vec<(String, Section)>,
    next_offset: u64,
    planned: bool,
}

impl ContainerSink {
    fn create(path: &Path) -> std::io::Result<Self> {
        Ok(Self {
            // Read+write: the SAindex pass preads the just-streamed SA
            // section through a cloned handle.
            file: std::fs::OpenOptions::new()
                .read(true)
                .write(true)
                .create(true)
                .truncate(true)
                .open(path)?,
            names: Vec::new(),
            entries: Vec::new(),
            next_offset: 0,
            planned: false,
        })
    }
}

/// Map an I/O failure on the container file.
fn container_io(e: std::io::Error) -> StarError {
    StarError::Index(format!("container write failed: {e}"))
}

impl SectionSink for ContainerSink {
    fn plan(&mut self, names: &[&str]) -> std::result::Result<(), StarError> {
        assert!(!self.planned, "container sink planned twice");
        self.planned = true;
        let toc_len = names.iter().try_fold(0u64, |total, name| {
            total
                .checked_add(20 + name.len() as u64)
                .ok_or_else(|| StarError::Index("star index table overflow".into()))
        })?;
        let writer = &mut self.file;
        writer.write_all(INDEX_MAGIC).map_err(container_io)?;
        writer
            .write_all(&INDEX_VERSION.to_le_bytes())
            .map_err(container_io)?;
        writer
            .write_all(&(names.len() as u32).to_le_bytes())
            .map_err(container_io)?;
        writer
            .write_all(&toc_len.to_le_bytes())
            .map_err(container_io)?;
        writer
            .write_all(&0u64.to_le_bytes())
            .map_err(container_io)?;
        // Placeholder TOC, back-patched by `finish`.
        writer
            .write_all(vec![0u8; toc_len as usize].as_slice())
            .map_err(container_io)?;
        self.names = names.iter().map(|name| (*name).to_owned()).collect();
        self.next_offset = align_u64(INDEX_HEADER_SIZE + toc_len, MMAP_ALIGNMENT);
        Ok(())
    }

    fn write_section(&mut self, name: &str, data: &[u8]) -> std::result::Result<(), StarError> {
        let offset = self.next_offset;
        pad_to(&mut self.file, offset).map_err(container_io)?;
        self.file
            .write_all(data)
            .map_err(|e| StarError::Index(format!("writing section {name}: {e}")))?;
        self.next_offset = align_u64(
            offset
                .checked_add(data.len() as u64)
                .ok_or_else(|| StarError::Index("star index size overflow".into()))?,
            MMAP_ALIGNMENT,
        );
        self.entries
            .push((name.to_owned(), Section { offset, length: data.len() as u64 }));
        Ok(())
    }

    fn stream_section(
        &mut self,
        name: &str,
        stream: &mut dyn FnMut(&mut File) -> std::result::Result<u64, StarError>,
    ) -> std::result::Result<(), StarError> {
        let offset = self.next_offset;
        pad_to(&mut self.file, offset).map_err(container_io)?;
        let length = stream(&mut self.file).map_err(|e| StarError::Index(e.to_string()))?;
        self.next_offset = align_u64(
            offset
                .checked_add(length)
                .ok_or_else(|| StarError::Index("star index size overflow".into()))?,
            MMAP_ALIGNMENT,
        );
        self.entries.push((name.to_owned(), Section { offset, length }));
        Ok(())
    }

    fn section_reader(&self, name: &str) -> std::result::Result<(File, u64), StarError> {
        let entry = self
            .entries
            .iter()
            .find(|(written, _)| written == name)
            .ok_or_else(|| StarError::Index(format!("section {name} not written yet")))?;
        let file = self
            .file
            .try_clone()
            .map_err(|e| StarError::Index(format!("cloning container handle: {e}")))?;
        Ok((file, entry.1.offset))
    }
}

impl ContainerSink {
    /// Back-patch the table of contents with the final section layout and
    /// flush. Sections must appear in `plan()` order — the placeholder was
    /// sized for exactly that list.
    fn finish(mut self) -> std::result::Result<(), StarError> {
        if !self.planned || self.entries.len() != self.names.len() {
            return Err(StarError::Index(
                "container sink finished with missing sections".into(),
            ));
        }
        self.file
            .seek(SeekFrom::Start(INDEX_HEADER_SIZE))
            .map_err(container_io)?;
        for (name, section) in &self.entries {
            let name_len = u16::try_from(name.len())
                .map_err(|_| StarError::Index("star index filename is too long".into()))?;
            self.file
                .write_all(&name_len.to_le_bytes())
                .map_err(container_io)?;
            self.file
                .write_all(&0u16.to_le_bytes())
                .map_err(container_io)?;
            self.file
                .write_all(&section.offset.to_le_bytes())
                .map_err(container_io)?;
            self.file
                .write_all(&section.length.to_le_bytes())
                .map_err(container_io)?;
            self.file
                .write_all(name.as_bytes())
                .map_err(container_io)?;
        }
        self.file.flush().map_err(container_io)?;
        self.file.sync_all().map_err(container_io)?;
        Ok(())
    }
}

pub(crate) fn load(path: impl AsRef<Path>) -> Result<StarIndex> {
    let path = path.as_ref();
    let mut file = File::open(path)?;
    let sections = read_sections(&mut file)?;
    let mut params = default_parameters()?;
    params.genome_dir = path.to_path_buf();
    params.out_filter_multimap_nmax = 1_000;
    params.out_filter_multimap_score_range = 100;

    let chr_name = read_text_lines(&mut file, required(&sections, "chrName.txt")?)?;
    let chr_length = parse_u64_lines(
        &read_section(&mut file, required(&sections, "chrLength.txt")?)?,
        "chrLength.txt",
    )?;
    let chr_start = parse_u64_lines(
        &read_section(&mut file, required(&sections, "chrStart.txt")?)?,
        "chrStart.txt",
    )?;
    if chr_name.len() != chr_length.len() || chr_start.len() != chr_name.len() + 1 {
        return Err(Error::InvalidFormat(
            "star chromosome metadata lengths disagree".into(),
        ));
    }
    let n_chr_real = chr_name.len();
    let n_genome_real = chr_start[n_chr_real];
    let genome_section = required(&sections, "Genome")?;
    // The Genome section is exactly `n_genome` bytes (padded genome, sjdb
    // extension included) — the single source of truth for the genome length.
    let n_genome = genome_section.length;
    if n_genome < n_genome_real {
        return Err(Error::InvalidFormat(format!(
            "star Genome section ({n_genome}) shorter than the real genome ({n_genome_real})"
        )));
    }
    let genome_map = map_section(&file, genome_section, 0)?;
    advise_random(&genome_map);
    let genome = Genome {
        sequence: GenomeSeq::Mapped {
            fwd: Arc::new(genome_map),
            n_genome: usize_from(n_genome, "star genome length")?,
        },
        n_genome,
        n_genome_real,
        n_chr_real,
        chr_name,
        chr_length,
        chr_start,
    };

    let sa_section = required(&sections, "SA")?;
    let sa_bytes = read_section(&mut file, sa_section)?;
    let gstrand_bit = SuffixArray::calculate_gstrand_bit(n_genome);
    let word_length = gstrand_bit + 1;
    let sa_length = packed_length(sa_section.length, word_length)?;
    let suffix_array = SuffixArray {
        // Read the suffix array into owned memory: the bulk of the index is
        // accessed randomly, and demand paging a multi-GB mmap with
        // MADV_RANDOM semantics turns the first pass into I/O stalls.
        data: PackedArray::from_bytes(word_length, sa_length, sa_bytes),
        gstrand_bit,
        gstrand_mask: (1u64 << gstrand_bit) - 1,
    };

    let sai_section = required(&sections, "SAindex")?;
    let sai_header = read_section_prefix(&mut file, sai_section, 8)?;
    let nbases = u64::from_le_bytes(sai_header.try_into().unwrap());
    let nbases = u32::try_from(nbases)
        .map_err(|_| Error::InvalidFormat("star SAindex nbases is too large".into()))?;
    let header_len = 8u64
        .checked_add(8u64.saturating_mul(u64::from(nbases) + 1))
        .ok_or_else(|| Error::InvalidFormat("star SAindex header overflow".into()))?;
    if header_len > sai_section.length {
        return Err(Error::InvalidFormat(
            "truncated star SAindex header".into(),
        ));
    }
    let header = read_section_prefix(
        &mut file,
        sai_section,
        usize_from(header_len, "star SAindex header")?,
    )?;
    let genome_sa_index_start = header[8..]
        .chunks_exact(8)
        .map(|bytes| u64::from_le_bytes(bytes.try_into().unwrap()))
        .collect::<Vec<_>>();
    let sai_data = {
        let mut sai_bytes = read_section(&mut file, sai_section)?;
        sai_bytes.drain(..usize_from(header_len, "star SAindex header")?);
        sai_bytes
    };
    let sai_word_length = gstrand_bit + 3;
    let num_indices = SaIndex::calculate_num_indices(nbases);
    let expected_sai_len = PackedArray::data_byte_len_for(
        sai_word_length,
        usize_from(num_indices, "star SAindex entries")?,
    ) as u64;
    if sai_section.length - header_len != expected_sai_len {
        return Err(Error::InvalidFormat(
            "star SAindex packed length mismatch".into(),
        ));
    }
    let sa_index = SaIndex {
        nbases,
        genome_sa_index_start,
        // Owned memory, like the suffix array above: SAindex lookups are the
        // hottest random accesses during seeding.
        data: PackedArray::from_bytes(
            sai_word_length,
            usize_from(num_indices, "star SAindex entries")?,
            sai_data,
        ),
        gstrand_bit,
    };

    let (prepared_junctions, sjdb_overhang) = sections
        .get("sjdbInfo.txt")
        .map(|section| parse_sjdb_info(&read_section(&mut file, section)?, &genome))
        .transpose()?
        .unwrap_or_default();
    let junction_db = if prepared_junctions.is_empty() {
        SpliceJunctionDb::empty()
    } else {
        let raw = prepared_junctions
            .iter()
            .map(|junction| {
                (
                    junction.chr_idx,
                    junction.stored_start(),
                    junction.stored_end(),
                    junction.strand,
                )
            })
            .collect::<Vec<_>>();
        SpliceJunctionDb::from_raw_junctions(&raw)
    };
    params.sjdb_overhang = sjdb_overhang;

    Ok(StarIndex {
        inner: GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db,
            prepared_junctions,
            sjdb_overhang,
        },
        params,
    })
}

fn default_parameters() -> Result<Parameters> {
    Ok(Parameters::seqarc_defaults())
}

fn read_sections(file: &mut File) -> Result<BTreeMap<String, Section>> {
    let file_len = file.metadata()?.len();
    let mut header = [0; INDEX_HEADER_SIZE as usize];
    file.read_exact(&mut header)?;
    if &header[..8] != INDEX_MAGIC {
        return Err(Error::InvalidFormat("bad star index magic".into()));
    }
    if u32::from_le_bytes(header[8..12].try_into().unwrap()) != INDEX_VERSION {
        return Err(Error::InvalidFormat(
            "unsupported star index version".into(),
        ));
    }
    let count = u32::from_le_bytes(header[12..16].try_into().unwrap()) as usize;
    let toc_len = u64::from_le_bytes(header[16..24].try_into().unwrap());
    if u64::from_le_bytes(header[24..32].try_into().unwrap()) != 0
        || INDEX_HEADER_SIZE.saturating_add(toc_len) > file_len
    {
        return Err(Error::InvalidFormat("invalid star index header".into()));
    }
    let toc = read_at(
        file,
        INDEX_HEADER_SIZE,
        usize_from(toc_len, "star index table")?,
    )?;
    let mut cursor = 0usize;
    let mut sections = BTreeMap::new();
    for _ in 0..count {
        if cursor + 20 > toc.len() {
            return Err(Error::InvalidFormat("truncated star index table".into()));
        }
        let name_len = u16::from_le_bytes(toc[cursor..cursor + 2].try_into().unwrap()) as usize;
        let flags = u16::from_le_bytes(toc[cursor + 2..cursor + 4].try_into().unwrap());
        let offset = u64::from_le_bytes(toc[cursor + 4..cursor + 12].try_into().unwrap());
        let length = u64::from_le_bytes(toc[cursor + 12..cursor + 20].try_into().unwrap());
        cursor += 20;
        if flags != 0 || cursor + name_len > toc.len() {
            return Err(Error::InvalidFormat("invalid star index entry".into()));
        }
        let name = std::str::from_utf8(&toc[cursor..cursor + name_len])
            .map_err(|_| Error::InvalidFormat("non-UTF-8 star index filename".into()))?
            .to_owned();
        cursor += name_len;
        if name.is_empty() || name.contains('/') || name.contains('\\') {
            return Err(Error::InvalidFormat("invalid star index filename".into()));
        }
        let end = offset
            .checked_add(length)
            .ok_or_else(|| Error::InvalidFormat("star index section overflow".into()))?;
        if offset < INDEX_HEADER_SIZE + toc_len || end > file_len {
            return Err(Error::InvalidFormat(
                "star index section is out of bounds".into(),
            ));
        }
        if sections.insert(name, Section { offset, length }).is_some() {
            return Err(Error::InvalidFormat(
                "duplicate star index section".into(),
            ));
        }
    }
    if cursor != toc.len() {
        return Err(Error::InvalidFormat(
            "trailing star index table bytes".into(),
        ));
    }
    let mut ranges = sections.values().collect::<Vec<_>>();
    ranges.sort_by_key(|section| section.offset);
    if ranges
        .windows(2)
        .any(|pair| pair[0].offset + pair[0].length > pair[1].offset)
    {
        return Err(Error::InvalidFormat(
            "overlapping star index sections".into(),
        ));
    }
    Ok(sections)
}

fn required<'a>(sections: &'a BTreeMap<String, Section>, name: &str) -> Result<&'a Section> {
    sections
        .get(name)
        .ok_or_else(|| Error::InvalidFormat(format!("missing star index section {name}")))
}

fn read_text_lines(file: &mut File, section: &Section) -> Result<Vec<String>> {
    let bytes = read_section(file, section)?;
    let text = std::str::from_utf8(&bytes)
        .map_err(|_| Error::InvalidFormat("star metadata is not UTF-8".into()))?;
    Ok(text.lines().map(str::to_owned).collect())
}

fn parse_u64_lines(bytes: &[u8], name: &str) -> Result<Vec<u64>> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| Error::InvalidFormat(format!("{name} is not UTF-8")))?;
    text.lines()
        .map(|line| {
            line.parse::<u64>()
                .map_err(|_| Error::InvalidFormat(format!("invalid integer in {name}")))
        })
        .collect()
}

fn parse_sjdb_info(bytes: &[u8], genome: &Genome) -> Result<(Vec<PreparedJunction>, u32)> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| Error::InvalidFormat("sjdbInfo.txt is not UTF-8".into()))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| Error::InvalidFormat("empty sjdbInfo.txt".into()))?;
    let mut header = header.split('\t');
    let expected = parse_field::<usize>(header.next(), "sjdbInfo count")?;
    let overhang = parse_field::<u32>(header.next(), "sjdbInfo overhang")?;
    let mut junctions = Vec::with_capacity(expected);
    for line in lines.filter(|line| !line.is_empty()) {
        let mut fields = line.split('\t');
        let stored_start = parse_field::<u64>(fields.next(), "sjdb stored start")?;
        let stored_end = parse_field::<u64>(fields.next(), "sjdb stored end")?;
        let motif = parse_field::<u8>(fields.next(), "sjdb motif")?;
        let shift_left = parse_field::<u8>(fields.next(), "sjdb left shift")?;
        let shift_right = parse_field::<u8>(fields.next(), "sjdb right shift")?;
        let strand = parse_field::<u8>(fields.next(), "sjdb strand")?;
        if fields.next().is_some() {
            return Err(Error::InvalidFormat("extra sjdbInfo fields".into()));
        }
        let (start_pos, end_pos) = if motif == 0 {
            (stored_start, stored_end)
        } else {
            (
                stored_start.saturating_sub(u64::from(shift_left)),
                stored_end.saturating_sub(u64::from(shift_left)),
            )
        };
        let original_start = start_pos + u64::from(shift_left);
        let chr_idx = genome
            .chr_start
            .iter()
            .take(genome.n_chr_real)
            .rposition(|&start| start <= original_start)
            .ok_or_else(|| Error::InvalidFormat("sjdb junction precedes genome".into()))?;
        if original_start >= genome.chr_start[chr_idx] + genome.chr_length[chr_idx] {
            return Err(Error::InvalidFormat(
                "sjdb junction falls in chromosome padding".into(),
            ));
        }
        junctions.push(PreparedJunction {
            chr_idx,
            start_pos,
            end_pos,
            motif,
            shift_left,
            shift_right,
            strand,
        });
    }
    if junctions.len() != expected {
        return Err(Error::InvalidFormat("sjdbInfo count mismatch".into()));
    }
    Ok((junctions, overhang))
}

fn parse_field<T: std::str::FromStr>(value: Option<&str>, field: &str) -> Result<T> {
    value
        .and_then(|value| value.parse().ok())
        .ok_or_else(|| Error::InvalidFormat(format!("invalid {field}")))
}

fn read_section(file: &mut File, section: &Section) -> Result<Vec<u8>> {
    read_at(
        file,
        section.offset,
        usize_from(section.length, "star index section")?,
    )
}

fn read_section_prefix(file: &mut File, section: &Section, length: usize) -> Result<Vec<u8>> {
    if length as u64 > section.length {
        return Err(Error::InvalidFormat(
            "truncated star index section".into(),
        ));
    }
    read_at(file, section.offset, length)
}

fn read_at(file: &mut File, offset: u64, length: usize) -> Result<Vec<u8>> {
    file.seek(SeekFrom::Start(offset))?;
    let mut bytes = vec![0; length];
    file.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn map_section(file: &File, section: &Section, relative_offset: u64) -> Result<Mmap> {
    if relative_offset > section.length {
        return Err(Error::InvalidFormat(
            "star mmap offset exceeds section".into(),
        ));
    }
    let length = section.length - relative_offset;
    if length == 0 {
        return Err(Error::InvalidFormat("empty star mmap section".into()));
    }
    unsafe {
        MmapOptions::new()
            .offset(section.offset + relative_offset)
            .len(usize_from(length, "star mmap section")?)
            .map(file)
            .map_err(Error::from)
    }
}

#[cfg(unix)]
fn advise_random(mmap: &Mmap) {
    let _ = mmap.advise(memmap2::Advice::Random);
}

#[cfg(not(unix))]
fn advise_random(_mmap: &Mmap) {}

fn packed_length(byte_length: u64, word_length: u32) -> Result<usize> {
    if byte_length < 8 {
        return Err(Error::InvalidFormat(
            "star packed array is too short".into(),
        ));
    }
    let bits = (byte_length - 8)
        .checked_mul(8)
        .ok_or_else(|| Error::InvalidFormat("star packed length overflow".into()))?;
    usize_from(
        bits.div_ceil(u64::from(word_length)) + 1,
        "star packed entries",
    )
}

fn pad_to(file: &mut File, offset: u64) -> std::io::Result<()> {
    let position = file.stream_position()?;
    if position > offset {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "star index layout overlap",
        ));
    }
    let mut remaining = offset - position;
    let zeros = [0u8; 8192];
    while remaining != 0 {
        let length = usize::try_from(remaining.min(zeros.len() as u64)).unwrap();
        file.write_all(&zeros[..length])?;
        remaining -= length as u64;
    }
    Ok(())
}

fn align_u64(value: u64, alignment: u64) -> u64 {
    value.div_ceil(alignment) * alignment
}

fn usize_from(value: u64, field: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| Error::InvalidFormat(format!("{field} is too large")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_random_dna(length: usize) -> Vec<u8> {
        let mut state = 0x9e37_79b9u32;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 17;
                state ^= state << 5;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    #[test]
    fn index_build_is_deterministic() {
        // No temp paths or timestamps may leak into the container: two builds
        // from identical inputs must be byte-identical.
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        let gtf = directory.path().join("genes.gtf");
        std::fs::write(&fasta, [b">chr1\n".as_slice(), &[b'A'; 600], b"\n"].concat()).unwrap();
        std::fs::write(
            &gtf,
            b"chr1\ttest\texon\t1\t60\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr1\ttest\texon\t301\t360\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
        )
        .unwrap();
        let options = GenomeIndexBuildOptions {
            temp_dir: Some(directory.path().to_path_buf()),
            aligner: crate::AlignerKind::Star,
            threads: 1,
            star_gtf: Some(gtf),
            star_sjdb_overhang: 30,
            star_sa_sparse_d: 1,
            star_sa_index_nbases: 3,
            star_chr_bin_nbits: 8,
        };
        let one = directory.path().join("one.seqarc-index");
        let two = directory.path().join("two.seqarc-index");
        build(&fasta, &one, &options).unwrap();
        build(&fasta, &two, &options).unwrap();
        assert_eq!(
            std::fs::read(&one).unwrap(),
            std::fs::read(&two).unwrap(),
        );
    }

    #[test]
    fn single_file_index_loads_and_aligns_across_a_splice() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        let gtf = directory.path().join("genes.gtf");
        let output = directory.path().join("genome.seqarc-index");
        let exon_1 = pseudo_random_dna(60);
        let intron_body = pseudo_random_dna(100);
        let exon_2 = pseudo_random_dna(60);
        let reference = [
            exon_1.as_slice(),
            b"GT",
            intron_body.as_slice(),
            b"AG",
            exon_2.as_slice(),
        ]
        .concat();
        std::fs::write(
            &fasta,
            [b">chr1\n".as_slice(), reference.as_slice(), b"\n"].concat(),
        )
        .unwrap();
        std::fs::write(
            &gtf,
            b"chr1\ttest\texon\t1\t60\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr1\ttest\texon\t165\t224\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
        )
        .unwrap();
        let options = GenomeIndexBuildOptions {
            temp_dir: Some(directory.path().to_path_buf()),
            aligner: crate::AlignerKind::Star,
            threads: 1,
            star_gtf: Some(gtf),
            star_sjdb_overhang: 30,
            star_sa_sparse_d: 1,
            star_sa_index_nbases: 3,
            star_chr_bin_nbits: 8,
        };

        build(&fasta, &output, &options).unwrap();
        let index = load(&output).unwrap();
        assert_eq!(index.reference_len(), 256);
        assert!(!index.inner.prepared_junctions.is_empty());

        let read = [&exon_1[30..], &exon_2[..30]].concat();
        let encoded = read
            .iter()
            .copied()
            .map(crate::star::io::fastq::encode_base)
            .collect::<Vec<_>>();
        let (transcripts, _, _, _) =
            crate::star::align::align_read(&encoded, "spliced", &index.inner, &index.params)
                .unwrap();
        assert!(transcripts.iter().any(|transcript| {
            transcript
                .cigar
                .iter()
                .any(|op| op.kind() == noodles::sam::alignment::record::cigar::op::Kind::Skip)
        }));

        let fastq = directory.path().join("reads.fastq");
        let archive = directory.path().join("reads.seqarc");
        let mut fastq_bytes = b"@spliced\n".to_vec();
        fastq_bytes.extend_from_slice(&read);
        fastq_bytes.extend_from_slice(b"\n+\n");
        fastq_bytes.extend(std::iter::repeat_n(b'I', read.len()));
        fastq_bytes.push(b'\n');
        std::fs::write(&fastq, fastq_bytes).unwrap();
        let compression_options = crate::AlignmentCompressionOptions {
            threads: 1,
            chunk_records: 1,
            ..crate::AlignmentCompressionOptions::default()
        };
        let (_, profile) = crate::compress_aligned_single_fastq_file_with_index_profiled(
            &output,
            &fastq,
            &archive,
            &compression_options,
        )
        .unwrap();
        assert_eq!(profile.mapped_reads, 1);
        let decoded = crate::decompress_aligned_single(
            File::open(&archive).unwrap(),
            &crate::DecompressionOptions::default(),
        )
        .unwrap();
        assert_eq!(decoded.streams[0][0].sequence, read);

        let bad_version = directory.path().join("bad-version.seqarc-index");
        std::fs::copy(&output, &bad_version).unwrap();
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .open(&bad_version)
            .unwrap();
        file.seek(SeekFrom::Start(8)).unwrap();
        file.write_all(&2u32.to_le_bytes()).unwrap();
        drop(file);
        assert!(load(&bad_version).is_err());

        let truncated = directory.path().join("truncated.seqarc-index");
        std::fs::copy(&output, &truncated).unwrap();
        std::fs::OpenOptions::new()
            .write(true)
            .open(&truncated)
            .unwrap()
            .set_len(64)
            .unwrap();
        assert!(load(&truncated).is_err());
    }

    #[test]
    fn stripped_intron_bases_shift_junction_and_spliced_read_still_aligns() {
        let directory = tempfile::tempdir().unwrap();
        let fasta = directory.path().join("genome.fa");
        let gtf = directory.path().join("genes.gtf");
        let output = directory.path().join("genome.seqarc-index");

        let exon_1 = pseudo_random_dna(60);
        let intron_left = pseudo_random_dna(20);
        let intron_right = pseudo_random_dna(30);
        let exon_2 = pseudo_random_dna(60);
        // Intron = GT + 20 + NNNN + 30 + AG: 58 original bases, 4 of them
        // stripped, so the transformed intron is 54 kept bases.
        let reference = [
            exon_1.as_slice(),
            b"GT",
            intron_left.as_slice(),
            b"NNNN",
            intron_right.as_slice(),
            b"AG",
            exon_2.as_slice(),
        ]
        .concat();
        std::fs::write(
            &fasta,
            [b">chr1\n".as_slice(), reference.as_slice(), b"\n"].concat(),
        )
        .unwrap();
        // Exon 1 covers 1..60; exon 2 starts at original position
        // 60 + 58 + 1 = 119 and ends at 60 + 58 + 60 = 178.
        std::fs::write(
            &gtf,
            b"chr1\ttest\texon\t1\t60\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr1\ttest\texon\t119\t178\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
        )
        .unwrap();
        let options = GenomeIndexBuildOptions {
            temp_dir: Some(directory.path().to_path_buf()),
            aligner: crate::AlignerKind::Star,
            threads: 1,
            star_gtf: Some(gtf),
            star_sjdb_overhang: 30,
            star_sa_sparse_d: 1,
            star_sa_index_nbases: 3,
            star_chr_bin_nbits: 8,
        };

        build(&fasta, &output, &options).unwrap();
        let index = load(&output).unwrap();
        assert_eq!(index.reference_len(), 256);
        assert!(!index.inner.prepared_junctions.is_empty());

        let read = [&exon_1[30..], &exon_2[..30]].concat();
        let encoded = read
            .iter()
            .copied()
            .map(crate::star::io::fastq::encode_base)
            .collect::<Vec<_>>();
        let (transcripts, _, _, _) =
            crate::star::align::align_read(&encoded, "spliced", &index.inner, &index.params)
                .unwrap();
        assert!(
            transcripts.iter().any(|transcript| {
                transcript.cigar.iter().any(|op| {
                    op.kind() == noodles::sam::alignment::record::cigar::op::Kind::Skip
                        && op.len() == 54
                })
            }),
            "expected a spliced alignment with a 54-base intron, got: {transcripts:?}"
        );
    }
}
