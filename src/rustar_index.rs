use crate::{temporary_directory, Error, GenomeIndexBuildOptions, Result};
use memmap2::{Mmap, MmapOptions};
use rustar_aligner::genome::{Genome, GenomeSeq};
use rustar_aligner::index::packed_array::PackedArray;
use rustar_aligner::index::sa_index::SaIndex;
use rustar_aligner::index::suffix_array::SuffixArray;
use rustar_aligner::index::GenomeIndex;
use rustar_aligner::junction::sjdb_insert::PreparedJunction;
use rustar_aligner::junction::SpliceJunctionDb;
use rustar_aligner::params::Parameters;
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

pub(crate) struct RustarIndex {
    pub(crate) inner: GenomeIndex,
    pub(crate) params: Parameters,
}

impl RustarIndex {
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
    if !(1..=14).contains(&options.rustar_sa_index_nbases) {
        return Err(Error::InvalidInput(
            "rustar SA index prefix length must be in 1..=14".into(),
        ));
    }
    if options.rustar_chr_bin_nbits > 30 {
        return Err(Error::InvalidInput(
            "rustar chromosome bin exponent must be at most 30".into(),
        ));
    }
    let workspace = temporary_directory(options.temp_dir.as_deref(), "seqarc-rustar-index-")?;
    let index_dir = workspace.path().join("star");
    std::fs::create_dir(&index_dir)?;
    let mut params = default_parameters()?;
    params.run_mode_in = vec!["genomeGenerate".into()];
    params.run_thread_n = NonZeroUsize::new(options.threads.max(1)).unwrap();
    params.temp_dir = Some(workspace.path().join("work"));
    std::fs::create_dir_all(workspace.path().join("work"))?;
    params.genome_dir = index_dir.clone();
    let normalized = crate::aligners::normalize_genome(
        genome_fasta.as_ref(),
        options.rustar_gtf.as_deref(),
        workspace.path(),
    )?;
    params.genome_fasta_files = vec![normalized.fasta];
    params.sjdb_gtf_file = normalized.gtf;
    params.genome_sa_sparse_d = options.rustar_sa_sparse_d.max(1);
    params.genome_sa_index_nbases = options.rustar_sa_index_nbases;
    params.genome_chr_bin_nbits = options.rustar_chr_bin_nbits;
    params.sjdb_overhang = options.rustar_sjdb_overhang;
    GenomeIndex::generate_streaming(&params)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    pack_directory(&index_dir, output.as_ref())
}

pub(crate) fn load(path: impl AsRef<Path>) -> Result<RustarIndex> {
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
            "rustar chromosome metadata lengths disagree".into(),
        ));
    }
    let n_chr_real = chr_name.len();
    let n_genome_real = chr_start[n_chr_real];
    let n_genome = sections
        .get("genomeParameters.txt")
        .map(|section| read_section(&mut file, section))
        .transpose()?
        .and_then(|bytes| parse_genome_size(&bytes))
        .unwrap_or(n_genome_real);
    let genome_section = required(&sections, "Genome")?;
    if genome_section.length != n_genome {
        return Err(Error::InvalidFormat(format!(
            "rustar Genome section length mismatch: expected {n_genome}, got {}",
            genome_section.length
        )));
    }
    let genome_map = map_section(&file, genome_section, 0)?;
    advise_random(&genome_map);
    let genome = Genome {
        sequence: GenomeSeq::Mapped {
            fwd: Arc::new(genome_map),
            n_genome: usize_from(n_genome, "rustar genome length")?,
        },
        n_genome,
        n_genome_real,
        n_chr_real,
        chr_name,
        chr_length,
        chr_start,
        transform_blocks: None,
    };

    let sa_section = required(&sections, "SA")?;
    let sa_map = map_section(&file, sa_section, 0)?;
    advise_random(&sa_map);
    let gstrand_bit = SuffixArray::calculate_gstrand_bit(n_genome);
    let word_length = gstrand_bit + 1;
    let sa_length = packed_length(sa_section.length, word_length)?;
    let suffix_array = SuffixArray {
        data: PackedArray::from_mmap(word_length, sa_length, sa_map),
        gstrand_bit,
        gstrand_mask: (1u64 << gstrand_bit) - 1,
    };

    let sai_section = required(&sections, "SAindex")?;
    let sai_header = read_section_prefix(&mut file, sai_section, 8)?;
    let nbases = u64::from_le_bytes(sai_header.try_into().unwrap());
    let nbases = u32::try_from(nbases)
        .map_err(|_| Error::InvalidFormat("rustar SAindex nbases is too large".into()))?;
    let header_len = 8u64
        .checked_add(8u64.saturating_mul(u64::from(nbases) + 1))
        .ok_or_else(|| Error::InvalidFormat("rustar SAindex header overflow".into()))?;
    if header_len > sai_section.length {
        return Err(Error::InvalidFormat(
            "truncated rustar SAindex header".into(),
        ));
    }
    let header = read_section_prefix(
        &mut file,
        sai_section,
        usize_from(header_len, "rustar SAindex header")?,
    )?;
    let genome_sa_index_start = header[8..]
        .chunks_exact(8)
        .map(|bytes| u64::from_le_bytes(bytes.try_into().unwrap()))
        .collect::<Vec<_>>();
    let sai_map = map_section(&file, sai_section, header_len)?;
    advise_random(&sai_map);
    let sai_word_length = gstrand_bit + 3;
    let num_indices = SaIndex::calculate_num_indices(nbases);
    let expected_sai_len = PackedArray::data_byte_len_for(
        sai_word_length,
        usize_from(num_indices, "rustar SAindex entries")?,
    ) as u64;
    if sai_section.length - header_len != expected_sai_len {
        return Err(Error::InvalidFormat(
            "rustar SAindex packed length mismatch".into(),
        ));
    }
    let sa_index = SaIndex {
        nbases,
        genome_sa_index_start,
        data: PackedArray::from_mmap(
            sai_word_length,
            usize_from(num_indices, "rustar SAindex entries")?,
            sai_map,
        ),
        word_length: sai_word_length,
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

    Ok(RustarIndex {
        inner: GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db,
            transcriptome: None,
            prepared_junctions,
            sjdb_overhang,
        },
        params,
    })
}

fn default_parameters() -> Result<Parameters> {
    Parameters::try_parse_from(["seqarc", "--readFilesIn", "reads.fastq"])
        .map_err(|error| Error::InvalidInput(format!("rustar parameters: {error}")))
}

fn pack_directory(directory: &Path, output: &Path) -> Result<()> {
    let mut files = std::fs::read_dir(directory)?
        .map(|entry| {
            let entry = entry?;
            let file_type = entry.file_type()?;
            if !file_type.is_file() {
                return Ok(None);
            }
            let name = entry.file_name().into_string().map_err(|_| {
                std::io::Error::new(std::io::ErrorKind::InvalidData, "non-UTF-8 index filename")
            })?;
            let length = entry.metadata()?.len();
            Ok(Some((name, entry.path(), length)))
        })
        .collect::<std::io::Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    files.sort_by(|left, right| left.0.cmp(&right.0));
    if files.is_empty() {
        return Err(Error::InvalidInput(
            "rustar index generation produced no files".into(),
        ));
    }
    let toc_len = files.iter().try_fold(0u64, |total, (name, _, _)| {
        total
            .checked_add(20 + name.len() as u64)
            .ok_or_else(|| Error::InvalidInput("rustar index table overflow".into()))
    })?;
    let mut next_offset = align_u64(INDEX_HEADER_SIZE + toc_len, MMAP_ALIGNMENT);
    let mut entries = Vec::with_capacity(files.len());
    for (name, path, length) in files {
        let offset = next_offset;
        next_offset = align_u64(
            offset
                .checked_add(length)
                .ok_or_else(|| Error::InvalidInput("rustar index size overflow".into()))?,
            MMAP_ALIGNMENT,
        );
        entries.push((name, path, Section { offset, length }));
    }

    let mut writer = File::create(output)?;
    writer.write_all(INDEX_MAGIC)?;
    writer.write_all(&INDEX_VERSION.to_le_bytes())?;
    writer.write_all(&(entries.len() as u32).to_le_bytes())?;
    writer.write_all(&toc_len.to_le_bytes())?;
    writer.write_all(&0u64.to_le_bytes())?;
    for (name, _, section) in &entries {
        let name_len = u16::try_from(name.len())
            .map_err(|_| Error::InvalidInput("rustar index filename is too long".into()))?;
        writer.write_all(&name_len.to_le_bytes())?;
        writer.write_all(&0u16.to_le_bytes())?;
        writer.write_all(&section.offset.to_le_bytes())?;
        writer.write_all(&section.length.to_le_bytes())?;
        writer.write_all(name.as_bytes())?;
    }
    for (_, path, section) in entries {
        pad_to(&mut writer, section.offset)?;
        let mut reader = File::open(path)?;
        let copied = std::io::copy(&mut reader, &mut writer)?;
        if copied != section.length {
            return Err(Error::InvalidInput(
                "rustar index component changed while packing".into(),
            ));
        }
    }
    writer.flush()?;
    Ok(())
}

fn read_sections(file: &mut File) -> Result<BTreeMap<String, Section>> {
    let file_len = file.metadata()?.len();
    let mut header = [0; INDEX_HEADER_SIZE as usize];
    file.read_exact(&mut header)?;
    if &header[..8] != INDEX_MAGIC {
        return Err(Error::InvalidFormat("bad rustar index magic".into()));
    }
    if u32::from_le_bytes(header[8..12].try_into().unwrap()) != INDEX_VERSION {
        return Err(Error::InvalidFormat(
            "unsupported rustar index version".into(),
        ));
    }
    let count = u32::from_le_bytes(header[12..16].try_into().unwrap()) as usize;
    let toc_len = u64::from_le_bytes(header[16..24].try_into().unwrap());
    if u64::from_le_bytes(header[24..32].try_into().unwrap()) != 0
        || INDEX_HEADER_SIZE.saturating_add(toc_len) > file_len
    {
        return Err(Error::InvalidFormat("invalid rustar index header".into()));
    }
    let toc = read_at(
        file,
        INDEX_HEADER_SIZE,
        usize_from(toc_len, "rustar index table")?,
    )?;
    let mut cursor = 0usize;
    let mut sections = BTreeMap::new();
    for _ in 0..count {
        if cursor + 20 > toc.len() {
            return Err(Error::InvalidFormat("truncated rustar index table".into()));
        }
        let name_len = u16::from_le_bytes(toc[cursor..cursor + 2].try_into().unwrap()) as usize;
        let flags = u16::from_le_bytes(toc[cursor + 2..cursor + 4].try_into().unwrap());
        let offset = u64::from_le_bytes(toc[cursor + 4..cursor + 12].try_into().unwrap());
        let length = u64::from_le_bytes(toc[cursor + 12..cursor + 20].try_into().unwrap());
        cursor += 20;
        if flags != 0 || cursor + name_len > toc.len() {
            return Err(Error::InvalidFormat("invalid rustar index entry".into()));
        }
        let name = std::str::from_utf8(&toc[cursor..cursor + name_len])
            .map_err(|_| Error::InvalidFormat("non-UTF-8 rustar index filename".into()))?
            .to_owned();
        cursor += name_len;
        if name.is_empty() || name.contains('/') || name.contains('\\') {
            return Err(Error::InvalidFormat("invalid rustar index filename".into()));
        }
        let end = offset
            .checked_add(length)
            .ok_or_else(|| Error::InvalidFormat("rustar index section overflow".into()))?;
        if offset < INDEX_HEADER_SIZE + toc_len || end > file_len {
            return Err(Error::InvalidFormat(
                "rustar index section is out of bounds".into(),
            ));
        }
        if sections.insert(name, Section { offset, length }).is_some() {
            return Err(Error::InvalidFormat(
                "duplicate rustar index section".into(),
            ));
        }
    }
    if cursor != toc.len() {
        return Err(Error::InvalidFormat(
            "trailing rustar index table bytes".into(),
        ));
    }
    let mut ranges = sections.values().collect::<Vec<_>>();
    ranges.sort_by_key(|section| section.offset);
    if ranges
        .windows(2)
        .any(|pair| pair[0].offset + pair[0].length > pair[1].offset)
    {
        return Err(Error::InvalidFormat(
            "overlapping rustar index sections".into(),
        ));
    }
    Ok(sections)
}

fn required<'a>(sections: &'a BTreeMap<String, Section>, name: &str) -> Result<&'a Section> {
    sections
        .get(name)
        .ok_or_else(|| Error::InvalidFormat(format!("missing rustar index section {name}")))
}

fn read_text_lines(file: &mut File, section: &Section) -> Result<Vec<String>> {
    let bytes = read_section(file, section)?;
    let text = std::str::from_utf8(&bytes)
        .map_err(|_| Error::InvalidFormat("rustar metadata is not UTF-8".into()))?;
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

fn parse_genome_size(bytes: &[u8]) -> Option<u64> {
    let text = std::str::from_utf8(bytes).ok()?;
    text.lines().find_map(|line| {
        line.strip_prefix("genomeFileSizes\t")?
            .split_whitespace()
            .next()?
            .parse()
            .ok()
    })
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
        usize_from(section.length, "rustar index section")?,
    )
}

fn read_section_prefix(file: &mut File, section: &Section, length: usize) -> Result<Vec<u8>> {
    if length as u64 > section.length {
        return Err(Error::InvalidFormat(
            "truncated rustar index section".into(),
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
            "rustar mmap offset exceeds section".into(),
        ));
    }
    let length = section.length - relative_offset;
    if length == 0 {
        return Err(Error::InvalidFormat("empty rustar mmap section".into()));
    }
    unsafe {
        MmapOptions::new()
            .offset(section.offset + relative_offset)
            .len(usize_from(length, "rustar mmap section")?)
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
            "rustar packed array is too short".into(),
        ));
    }
    let bits = (byte_length - 8)
        .checked_mul(8)
        .ok_or_else(|| Error::InvalidFormat("rustar packed length overflow".into()))?;
    usize_from(
        bits.div_ceil(u64::from(word_length)) + 1,
        "rustar packed entries",
    )
}

fn pad_to(file: &mut File, offset: u64) -> Result<()> {
    let position = file.stream_position()?;
    if position > offset {
        return Err(Error::InvalidInput("rustar index layout overlap".into()));
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
            aligner: crate::AlignerKind::Rustar,
            threads: 1,
            rustar_gtf: Some(gtf),
            rustar_sjdb_overhang: 30,
            rustar_sa_sparse_d: 1,
            rustar_sa_index_nbases: 3,
            rustar_chr_bin_nbits: 8,
        };

        build(&fasta, &output, &options).unwrap();
        let index = load(&output).unwrap();
        assert_eq!(index.reference_len(), 256);
        assert!(!index.inner.prepared_junctions.is_empty());

        let read = [&exon_1[30..], &exon_2[..30]].concat();
        let encoded = read
            .iter()
            .copied()
            .map(rustar_aligner::io::fastq::encode_base)
            .collect::<Vec<_>>();
        let (transcripts, _, _, _) =
            rustar_aligner::align::align_read(&encoded, "spliced", &index.inner, &index.params)
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
            aligner: crate::AlignerKind::Rustar,
            threads: 1,
            rustar_gtf: Some(gtf),
            rustar_sjdb_overhang: 30,
            rustar_sa_sparse_d: 1,
            rustar_sa_index_nbases: 3,
            rustar_chr_bin_nbits: 8,
        };

        build(&fasta, &output, &options).unwrap();
        let index = load(&output).unwrap();
        assert_eq!(index.reference_len(), 256);
        assert!(!index.inner.prepared_junctions.is_empty());

        let read = [&exon_1[30..], &exon_2[..30]].concat();
        let encoded = read
            .iter()
            .copied()
            .map(rustar_aligner::io::fastq::encode_base)
            .collect::<Vec<_>>();
        let (transcripts, _, _, _) =
            rustar_aligner::align::align_read(&encoded, "spliced", &index.inner, &index.params)
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
