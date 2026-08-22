use minibwa_rs as raw;
use std::fmt;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;

const INDEX_MAGIC: &[u8; 8] = b"SARCIDX1";
const INDEX_VERSION: u32 = 1;
const INDEX_HEADER_SIZE: u64 = 64;

#[derive(Debug)]
pub struct Error(String);

impl fmt::Display for Error {
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        output.write_str(&self.0)
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

pub struct Index {
    inner: raw::map_algo::mb_idx_t,
}

impl Index {
    pub fn build_from_fasta(
        fasta: impl AsRef<Path>,
        output: impl AsRef<Path>,
        methylation: bool,
        threads: u32,
    ) -> Result<()> {
        if methylation {
            return Err(Error("methylation indexes are not supported".into()));
        }
        let workspace = tempfile::tempdir().map_err(|error| Error(error.to_string()))?;
        let prefix = workspace.path().join("index");
        let mut arguments = vec![
            "index".to_owned(),
            "-t".to_owned(),
            threads.max(1).to_string(),
        ];
        arguments.push(fasta.as_ref().to_string_lossy().into_owned());
        arguments.push(prefix.to_string_lossy().into_owned());
        let (status, message) = raw::index::main_index(&arguments);
        if status != 0 {
            return Err(Error(if message.is_empty() {
                "minibwa index build failed".into()
            } else {
                message
            }));
        }
        pack_components(
            prefix.with_extension("l2b"),
            prefix.with_extension("mbw"),
            output,
        )
    }

    pub fn load(path: impl AsRef<Path>, methylation: bool) -> Result<Self> {
        if methylation {
            return Err(Error("methylation indexes are not supported".into()));
        }
        load_container(path)
    }

    pub fn contig_name(&self, tid: i32) -> Option<&str> {
        raw::map_algo::mb_idx_ctg_name(&self.inner, tid)
    }

    pub fn reference_len(&self) -> usize {
        self.inner.l2b.tot_len as usize
    }

    pub fn reference_base(&self, position: usize) -> Option<u8> {
        (position < self.reference_len()).then(|| self.reference_base_at(position))
    }

    #[inline(always)]
    pub(crate) fn reference_base_at(&self, position: usize) -> u8 {
        b"ACGT"[raw::l2bit::l2b_get0(&self.inner.l2b, position as u64) as usize]
    }
}

fn pack_components(
    l2b_path: impl AsRef<Path>,
    mbw_path: impl AsRef<Path>,
    output: impl AsRef<Path>,
) -> Result<()> {
    let mut l2b = File::open(l2b_path).map_err(|error| Error(error.to_string()))?;
    let mut mbw = File::open(mbw_path).map_err(|error| Error(error.to_string()))?;
    let l2b_len = l2b
        .metadata()
        .map_err(|error| Error(error.to_string()))?
        .len();
    let mbw_len = mbw
        .metadata()
        .map_err(|error| Error(error.to_string()))?
        .len();
    let mut prefix = [0; 24];
    l2b.read_exact(&mut prefix)
        .map_err(|error| Error(error.to_string()))?;
    if &prefix[..4] != b"L2B\x01" {
        return Err(Error("invalid L2B component".into()));
    }
    let reference_len = u64::from_le_bytes(prefix[16..24].try_into().unwrap());
    l2b.seek(SeekFrom::Start(0))
        .map_err(|error| Error(error.to_string()))?;

    let l2b_offset = INDEX_HEADER_SIZE;
    let mbw_offset = align_u64(l2b_offset + l2b_len, 8);
    let mut output = File::create(output).map_err(|error| Error(error.to_string()))?;
    output
        .write_all(&[0; INDEX_HEADER_SIZE as usize])
        .map_err(|error| Error(error.to_string()))?;
    std::io::copy(&mut l2b, &mut output).map_err(|error| Error(error.to_string()))?;
    while output
        .stream_position()
        .map_err(|error| Error(error.to_string()))?
        < mbw_offset
    {
        output
            .write_all(&[0])
            .map_err(|error| Error(error.to_string()))?;
    }
    std::io::copy(&mut mbw, &mut output).map_err(|error| Error(error.to_string()))?;
    output
        .seek(SeekFrom::Start(0))
        .map_err(|error| Error(error.to_string()))?;
    output
        .write_all(INDEX_MAGIC)
        .and_then(|_| output.write_all(&INDEX_VERSION.to_le_bytes()))
        .and_then(|_| output.write_all(&0u32.to_le_bytes()))
        .and_then(|_| output.write_all(&reference_len.to_le_bytes()))
        .and_then(|_| output.write_all(&l2b_offset.to_le_bytes()))
        .and_then(|_| output.write_all(&l2b_len.to_le_bytes()))
        .and_then(|_| output.write_all(&mbw_offset.to_le_bytes()))
        .and_then(|_| output.write_all(&mbw_len.to_le_bytes()))
        .and_then(|_| output.write_all(&0u64.to_le_bytes()))
        .map_err(|error| Error(error.to_string()))?;
    output.flush().map_err(|error| Error(error.to_string()))
}

fn load_container(path: impl AsRef<Path>) -> Result<Index> {
    #[cfg(all(unix, not(target_arch = "wasm32"), target_endian = "little"))]
    {
        let map = raw::mapped_u64::mmap_file(path.as_ref(), 0)
            .ok_or_else(|| Error("failed to mmap minibwa index".into()))?;
        let bytes = map.as_bytes();
        if bytes.len() < INDEX_HEADER_SIZE as usize || &bytes[..8] != INDEX_MAGIC {
            return Err(Error("invalid minibwa container magic".into()));
        }
        let version = read_u32_at(bytes, 8)?;
        if version != INDEX_VERSION {
            return Err(Error("unsupported minibwa container version".into()));
        }
        if read_u32_at(bytes, 12)? != 0 {
            return Err(Error("unsupported minibwa container flags".into()));
        }
        let reference_len = read_u64_at(bytes, 16)?;
        let l2b_offset = usize_from(read_u64_at(bytes, 24)?, "L2B offset")?;
        let l2b_len = usize_from(read_u64_at(bytes, 32)?, "L2B length")?;
        let mbw_offset = usize_from(read_u64_at(bytes, 40)?, "MBW offset")?;
        let mbw_len = usize_from(read_u64_at(bytes, 48)?, "MBW length")?;
        validate_section(bytes.len(), l2b_offset, l2b_len)?;
        validate_section(bytes.len(), mbw_offset, mbw_len)?;
        if l2b_offset % 8 != 0 || mbw_offset % 8 != 0 || l2b_offset + l2b_len > mbw_offset {
            return Err(Error("invalid minibwa container sections".into()));
        }
        let l2b = parse_l2b(Arc::clone(&map), l2b_offset, l2b_len)?;
        if l2b.tot_len != reference_len || l2b.n_ctg != 1 || l2b.n_ambi != 0 || l2b.n_mask != 0 {
            return Err(Error("invalid processed reference metadata".into()));
        }
        let mut bwt = parse_mbw(map, mbw_offset, mbw_len)?;
        if bwt.seq_len != reference_len.saturating_mul(2) {
            return Err(Error("reference and BWT lengths differ".into()));
        }
        raw::bwt::mb_bwt_cache(&mut bwt, 10);
        Ok(Index {
            inner: raw::map_algo::mb_idx_t {
                is_meth: 0,
                l2b,
                bwt,
            },
        })
    }
    #[cfg(not(all(unix, not(target_arch = "wasm32"), target_endian = "little")))]
    {
        let _ = path;
        Err(Error(
            "single-file minibwa indexes require little-endian Unix".into(),
        ))
    }
}

#[cfg(all(unix, not(target_arch = "wasm32"), target_endian = "little"))]
fn parse_l2b(
    map: Arc<raw::mapped_u64::UnixMmap>,
    base: usize,
    length: usize,
) -> Result<raw::l2bit::l2b_t> {
    let bytes = map.as_bytes();
    let section = &bytes[base..base + length];
    if section.len() < 64 || &section[..4] != b"L2B\x01" {
        return Err(Error("invalid L2B section".into()));
    }
    let fields = (0..7)
        .map(|index| read_u64_at(section, 8 + index * 8))
        .collect::<Result<Vec<_>>>()?;
    let (n_ctg, tot_len, n_ambi, n_mask, len_name, len_comm, n_pac) = (
        fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6],
    );
    if n_ctg != 1 || n_ambi != 0 || n_mask != 0 || n_pac != tot_len.div_ceil(32) {
        return Err(Error("invalid processed L2B metadata".into()));
    }
    let mut offset = 64usize;
    let mut contigs = Vec::with_capacity(n_ctg as usize);
    let mut contig_offset = 0u64;
    for _ in 0..n_ctg {
        let len = read_u64_at(section, offset)?;
        offset += 8;
        contigs.push(raw::l2bit::l2b_ctg_t {
            name: String::new(),
            comm: None,
            len,
            off: contig_offset,
        });
        contig_offset = contig_offset
            .checked_add(len)
            .ok_or_else(|| Error("contig length overflow".into()))?;
    }
    if contig_offset != tot_len {
        return Err(Error("invalid contig lengths".into()));
    }
    let ambi = read_intervals(section, &mut offset, n_ambi)?;
    let mask = read_intervals(section, &mut offset, n_mask)?;
    let pac_bytes = usize_from(n_pac, "packed reference length")?
        .checked_mul(8)
        .ok_or_else(|| Error("packed reference overflow".into()))?;
    let pac_end = offset
        .checked_add(pac_bytes)
        .ok_or_else(|| Error("packed reference overflow".into()))?;
    if pac_end > section.len() {
        return Err(Error("truncated packed reference".into()));
    }
    let pac = raw::mapped_u64::mapped_u64_slice(Arc::clone(&map), base + offset, n_pac as usize)
        .ok_or_else(|| Error("invalid packed reference range".into()))?;
    offset = pac_end;
    let names_end = offset
        .checked_add(len_name as usize)
        .ok_or_else(|| Error("reference name overflow".into()))?;
    let comments_end = names_end
        .checked_add(len_comm as usize)
        .ok_or_else(|| Error("reference comment overflow".into()))?;
    if comments_end > section.len() {
        return Err(Error("truncated reference metadata".into()));
    }
    let cat_name = section[offset..names_end].to_vec();
    let cat_comm = section[names_end..comments_end].to_vec();
    populate_contig_strings(&mut contigs, &cat_name, &cat_comm)?;
    Ok(raw::l2bit::l2b_t {
        tot_len,
        n_ctg,
        m_ctg: n_ctg,
        ctg: contigs,
        n_pac,
        m_pac: n_pac,
        n_ambi,
        m_ambi: n_ambi,
        n_mask,
        m_mask: n_mask,
        ambi,
        mask,
        pac,
        cat_name,
        cat_comm,
    })
}

#[cfg(all(unix, not(target_arch = "wasm32"), target_endian = "little"))]
fn parse_mbw(
    map: Arc<raw::mapped_u64::UnixMmap>,
    base: usize,
    length: usize,
) -> Result<raw::bwt::mb_bwt_t> {
    let bytes = map.as_bytes();
    let section = &bytes[base..base + length];
    if section.len() < 56 || &section[..4] != b"MBW\x02" {
        return Err(Error("invalid MBW section".into()));
    }
    let mut bwt = raw::bwt::mb_bwt_init();
    bwt.sa_bit = read_u32_at(section, 4)?;
    bwt.primary = read_u64_at(section, 8)?;
    for index in 1..5 {
        bwt.L2[index] = read_u64_at(section, 8 + index * 8)?;
    }
    bwt.seq_len = bwt.L2[4];
    bwt.data_len = raw::bwt::mb_bwt_data_len(bwt.seq_len);
    let data_bytes = usize_from(bwt.data_len, "BWT data length")?
        .checked_mul(8)
        .ok_or_else(|| Error("BWT size overflow".into()))?;
    let n_sa_offset = 48usize
        .checked_add(data_bytes)
        .ok_or_else(|| Error("BWT size overflow".into()))?;
    if n_sa_offset + 8 > section.len() {
        return Err(Error("truncated BWT data".into()));
    }
    bwt.data =
        raw::mapped_u64::mapped_u64_slice(Arc::clone(&map), base + 48, bwt.data_len as usize)
            .ok_or_else(|| Error("invalid BWT data range".into()))?;
    bwt.n_sa = read_u64_at(section, n_sa_offset)?;
    if bwt.sa_bit != u32::MAX && bwt.n_sa != 0 {
        let step = 1u64
            .checked_shl(bwt.sa_bit)
            .ok_or_else(|| Error("invalid suffix-array sampling".into()))?;
        let expected = (bwt.seq_len + step) >> bwt.sa_bit;
        if bwt.n_sa != expected {
            return Err(Error("invalid suffix-array length".into()));
        }
        let sa_bytes = usize_from(bwt.n_sa, "suffix-array length")?
            .checked_mul(8)
            .ok_or_else(|| Error("suffix-array size overflow".into()))?;
        let sa_offset = n_sa_offset + 8;
        if sa_offset
            .checked_add(sa_bytes)
            .is_none_or(|end| end > section.len())
        {
            return Err(Error("truncated suffix array".into()));
        }
        bwt.sa = raw::mapped_u64::mapped_u64_slice(map, base + sa_offset, bwt.n_sa as usize)
            .ok_or_else(|| Error("invalid suffix-array range".into()))?;
    }
    Ok(bwt)
}

fn read_intervals(
    section: &[u8],
    offset: &mut usize,
    count: u64,
) -> Result<Vec<raw::l2bit::l2b_intv_t>> {
    let mut intervals = Vec::with_capacity(count as usize);
    for _ in 0..count {
        let st = read_u64_at(section, *offset)?;
        let en = read_u64_at(section, *offset + 8)?;
        *offset += 16;
        intervals.push(raw::l2bit::l2b_intv_t { st, en });
    }
    Ok(intervals)
}

fn populate_contig_strings(
    contigs: &mut [raw::l2bit::l2b_ctg_t],
    names: &[u8],
    comments: &[u8],
) -> Result<()> {
    let mut name_offset = 0;
    let mut comment_offset = 0;
    for contig in contigs {
        let remaining_names = names
            .get(name_offset..)
            .ok_or_else(|| Error("invalid contig names".into()))?;
        let name_end = name_offset
            + remaining_names
                .iter()
                .position(|&value| value == 0)
                .ok_or_else(|| Error("invalid contig names".into()))?;
        contig.name = String::from_utf8_lossy(&names[name_offset..name_end]).into_owned();
        name_offset = name_end + 1;
        let remaining_comments = comments
            .get(comment_offset..)
            .ok_or_else(|| Error("invalid contig comments".into()))?;
        let comment_end = comment_offset
            + remaining_comments
                .iter()
                .position(|&value| value == 0)
                .ok_or_else(|| Error("invalid contig comments".into()))?;
        if comment_end != comment_offset {
            contig.comm =
                Some(String::from_utf8_lossy(&comments[comment_offset..comment_end]).into_owned());
        }
        comment_offset = comment_end + 1;
    }
    if name_offset != names.len() || comment_offset != comments.len() {
        return Err(Error("trailing contig metadata".into()));
    }
    Ok(())
}

fn read_u32_at(bytes: &[u8], offset: usize) -> Result<u32> {
    let end = offset
        .checked_add(4)
        .ok_or_else(|| Error("index offset overflow".into()))?;
    bytes
        .get(offset..end)
        .map(|value| u32::from_le_bytes(value.try_into().unwrap()))
        .ok_or_else(|| Error("truncated index".into()))
}

fn read_u64_at(bytes: &[u8], offset: usize) -> Result<u64> {
    let end = offset
        .checked_add(8)
        .ok_or_else(|| Error("index offset overflow".into()))?;
    bytes
        .get(offset..end)
        .map(|value| u64::from_le_bytes(value.try_into().unwrap()))
        .ok_or_else(|| Error("truncated index".into()))
}

fn validate_section(file_len: usize, offset: usize, length: usize) -> Result<()> {
    if offset < INDEX_HEADER_SIZE as usize
        || offset.checked_add(length).is_none_or(|end| end > file_len)
    {
        return Err(Error("index section out of bounds".into()));
    }
    Ok(())
}

fn usize_from(value: u64, what: &str) -> Result<usize> {
    usize::try_from(value).map_err(|_| Error(format!("{what} is too large")))
}

fn align_u64(value: u64, alignment: u64) -> u64 {
    (value + alignment - 1) & !(alignment - 1)
}

#[derive(Clone)]
pub struct Opts {
    inner: raw::options::mb_opt_t,
}

impl Opts {
    pub fn new() -> Self {
        let mut inner = raw::options::mb_opt_t::default();
        raw::options::mb_opt_init(&mut inner);
        Self { inner }
    }

    pub fn set_paired(mut self, paired: bool) -> Self {
        if paired {
            self.inner.flag |= raw::options::MB_F_PE;
        } else {
            self.inner.flag &= !raw::options::MB_F_PE;
        }
        self
    }
}

impl Default for Opts {
    fn default() -> Self {
        Self::new()
    }
}

pub struct ThreadBuf {
    inner: raw::map_algo::mb_tbuf_t,
}

impl ThreadBuf {
    pub fn new() -> Self {
        Self {
            inner: raw::map_algo::mb_tbuf_init(0),
        }
    }
}

impl Default for ThreadBuf {
    fn default() -> Self {
        Self::new()
    }
}

pub struct Aligner<'a> {
    index: &'a Index,
    options: &'a Opts,
}

impl<'a> Aligner<'a> {
    pub fn new(index: &'a Index, options: &'a Opts) -> Self {
        Self { index, options }
    }

    pub fn map_many<N: AsRef<[u8]>, S: AsRef<[u8]>>(
        &self,
        buffer: &mut ThreadBuf,
        names: &[N],
        sequences: &[S],
    ) -> Result<Vec<Vec<Hit>>> {
        if names.len() != sequences.len() {
            return Err(Error("query name and sequence counts differ".into()));
        }
        let names = names
            .iter()
            .map(|name| {
                std::str::from_utf8(name.as_ref())
                    .map_err(|_| Error("query name is not UTF-8".into()))
            })
            .collect::<Result<Vec<_>>>()?;
        let sequences = sequences
            .iter()
            .map(|sequence| {
                let sequence = sequence.as_ref();
                if sequence.is_empty() {
                    return Err(Error("query sequence is empty".into()));
                }
                std::str::from_utf8(sequence)
                    .map_err(|_| Error("query sequence is not UTF-8".into()))
            })
            .collect::<Result<Vec<_>>>()?;
        let lengths = sequences
            .iter()
            .map(|sequence| {
                i32::try_from(sequence.len())
                    .map_err(|_| Error("query sequence is too long".into()))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut hit_counts = vec![0; sequences.len()];
        let mut options = self.options.inner.clone();
        options.flag &= !raw::options::MB_F_PE;
        let hits = raw::map_algo::mb_map_batch(
            &options,
            &self.index.inner,
            sequences.len() as i32,
            &lengths,
            &sequences,
            &mut hit_counts,
            Some(&mut buffer.inner),
            Some(&names),
        );
        hits.into_iter()
            .map(|hits| {
                hits.into_iter()
                    .map(|hit| Hit::from_raw(hit, self.index))
                    .collect()
            })
            .collect()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarKind {
    Match,
    Ins,
    Del,
    RefSkip,
    SoftClip,
    HardClip,
    Pad,
    Eq,
    Diff,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CigarOp {
    pub kind: CigarKind,
    pub len: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hit {
    pub tid: i64,
    pub contig: Option<String>,
    pub ref_start: i64,
    pub ref_end: i64,
    pub query_start: i32,
    pub query_end: i32,
    pub strand: Strand,
    pub mapq: u8,
    pub score: i32,
    pub n_sub: i32,
    pub blen: i32,
    pub mlen: i32,
    pub proper_pair: bool,
    pub is_primary: bool,
    pub is_secondary: bool,
    pub is_supplementary: bool,
    pub cigar: Vec<CigarOp>,
}

impl Hit {
    fn from_raw(hit: raw::pe::mb_hit_t, index: &Index) -> Result<Self> {
        let cigar = hit
            .p
            .as_ref()
            .map(|extra| extra.cigar())
            .unwrap_or_default()
            .iter()
            .map(|&word| {
                Ok(CigarOp {
                    kind: cigar_kind(word & 0x0f)?,
                    len: word >> 4,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let is_secondary = hit.parent != hit.id;
        let is_supplementary = !is_secondary && hit.sam_pri() == 0;
        Ok(Self {
            tid: hit.tid,
            contig: index.contig_name(hit.tid as i32).map(str::to_owned),
            ref_start: hit.ts,
            ref_end: hit.te,
            query_start: hit.qs,
            query_end: hit.qe,
            strand: if hit.rev() != 0 {
                Strand::Reverse
            } else {
                Strand::Forward
            },
            mapq: hit.mapq.clamp(0, 255) as u8,
            score: hit.score,
            n_sub: hit.n_sub,
            blen: hit.blen,
            mlen: hit.mlen,
            proper_pair: hit.proper_pair() != 0,
            is_primary: !is_secondary && !is_supplementary,
            is_secondary,
            is_supplementary,
            cigar,
        })
    }
}

fn cigar_kind(code: u32) -> Result<CigarKind> {
    match code {
        0 => Ok(CigarKind::Match),
        1 => Ok(CigarKind::Ins),
        2 => Ok(CigarKind::Del),
        3 => Ok(CigarKind::RefSkip),
        4 => Ok(CigarKind::SoftClip),
        5 => Ok(CigarKind::HardClip),
        6 => Ok(CigarKind::Pad),
        7 => Ok(CigarKind::Eq),
        8 => Ok(CigarKind::Diff),
        _ => Err(Error("unknown minibwa CIGAR operation".into())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_foreign_index_magic() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("old.index");
        let mut header = vec![0; INDEX_HEADER_SIZE as usize];
        header[8..12].copy_from_slice(&INDEX_VERSION.to_le_bytes());
        std::fs::write(&path, header).unwrap();
        assert!(Index::load(path, false).is_err());
    }
}
