//! Global minimizer ordering and lossless consensus residuals.

use crate::reference_free::record_store::PairStore;
use crate::{Error, ReadRecord, ReorderOptions, Result};
use bed_utils::extsort::ExternalSorterBuilder;
use rkyv::{Archive, Deserialize, Serialize};
use std::collections::VecDeque;
#[cfg(test)]
use std::fs::{File, OpenOptions};
#[cfg(test)]
use std::io::{Read, Write};
use std::path::Path;

#[derive(Archive, Serialize, Deserialize, Clone, Debug, PartialEq, Eq, Ord, PartialOrd)]
#[rkyv(derive(Ord, PartialOrd, Eq, PartialEq))]
pub struct MinimizerEntry {
    pub hash: u64,
    pub read_id: u64,
    pub position: u32,
    pub reverse: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrderedRead {
    pub original_id: u64,
    pub record: ReadRecord,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConsensusCluster {
    pub reference: Vec<u8>,
    pub reads: Vec<ConsensusRead>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConsensusRead {
    pub original_id: u64,
    pub offset: u32,
    pub reverse: bool,
    pub length: u32,
    pub mismatches: Vec<(u32, u8)>,
}

pub fn canonical_minimizers(sequence: &[u8], k: usize, window: usize) -> Vec<(u64, u32, bool)> {
    if k == 0 || window == 0 || sequence.len() < k {
        return Vec::new();
    }
    let kmers = sequence.len() - k + 1;
    let span = window.min(kmers);
    let values = canonical_kmers(sequence, k);
    let mut result = Vec::new();
    let mut previous = None;
    let mut minima: VecDeque<(u64, usize, bool)> = VecDeque::with_capacity(span);
    for (position, value) in values.into_iter().enumerate() {
        if let Some((hash, reverse)) = value {
            while minima.back().is_some_and(|candidate| candidate.0 > hash) {
                minima.pop_back();
            }
            minima.push_back((hash, position, reverse));
        }
        if position + 1 < span {
            continue;
        }
        let window_start = position + 1 - span;
        while minima
            .front()
            .is_some_and(|candidate| candidate.1 < window_start)
        {
            minima.pop_front();
        }
        if let Some(&(hash, position, reverse)) = minima.front() {
            if previous != Some((hash, position)) {
                result.push((hash, position as u32, reverse));
                previous = Some((hash, position));
            }
        }
    }
    result
}

/// Pick the deterministic primary seed used to group reads for consensus.
pub(crate) fn primary_minimizer(
    sequence: &[u8],
    k: usize,
    window: usize,
) -> Option<(u64, u32, bool)> {
    let minimizers = canonical_minimizers(sequence, k, window);
    if minimizers.is_empty() {
        return None;
    }
    minimizers
        .into_iter()
        .min_by_key(|(hash, position, reverse)| (*hash, *position, *reverse))
}

fn canonical_kmers(sequence: &[u8], k: usize) -> Vec<Option<(u64, bool)>> {
    if k > 32 {
        return (0..=sequence.len() - k)
            .map(|start| canonical_kmer(&sequence[start..start + k]))
            .collect();
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let reverse_shift = 2 * (k - 1);
    let mut forward = 0u64;
    let mut reverse = 0u64;
    let mut valid = 0usize;
    let mut output = Vec::with_capacity(sequence.len() - k + 1);
    for (index, &base) in sequence.iter().enumerate() {
        if let Some(value) = base_value(base) {
            forward = ((forward << 2) | value as u64) & mask;
            reverse = (reverse >> 2) | (((3 - value) as u64) << reverse_shift);
            valid += 1;
        } else {
            forward = 0;
            reverse = 0;
            valid = 0;
        }
        if index + 1 < k {
            continue;
        }
        output.push(if valid >= k {
            Some(if forward <= reverse {
                (forward, false)
            } else {
                (reverse, true)
            })
        } else {
            None
        });
    }
    output
}

fn canonical_kmer(kmer: &[u8]) -> Option<(u64, bool)> {
    if kmer.len() > 32 {
        return Some((hash_bytes(kmer), false));
    }
    let mut forward = 0u64;
    let mut reverse = 0u64;
    for &base in kmer {
        let value = base_value(base)?;
        forward = (forward << 2) | value as u64;
        reverse = (reverse >> 2) | (((3 - value) as u64) << (2 * (kmer.len() - 1)));
    }
    if forward <= reverse {
        Some((forward, false))
    } else {
        Some((reverse, true))
    }
}

fn base_value(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn hash_bytes(bytes: &[u8]) -> u64 {
    let mut hash = 0x9e3779b185ebca87u64;
    for &byte in bytes {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
        hash ^= hash >> 29;
    }
    hash
}

pub fn sort_store_minimizers(
    store: &mut PairStore,
    options: &ReorderOptions,
    temp_dir: impl AsRef<Path>,
) -> Result<Vec<MinimizerEntry>> {
    let entries = (0..store.count()).flat_map(|read_id| {
        let read = store.read_first(read_id).ok();
        let minimizers = read
            .as_ref()
            .map(|record| {
                canonical_minimizers(&record.sequence, options.kmer_size, options.window_size)
            })
            .unwrap_or_default();
        minimizers
            .into_iter()
            .map(move |(hash, position, reverse)| MinimizerEntry {
                hash,
                read_id,
                position,
                reverse,
            })
    });
    let sorter = ExternalSorterBuilder::new()
        .with_chunk_size(options.sort_chunk_size.max(1))
        .with_tmp_dir(temp_dir)
        .with_compression(3)
        .num_threads(options.sort_threads.max(1))
        .build()
        .map_err(Error::Io)?;
    sorter
        .sort(entries)
        .map_err(|error| Error::InvalidInput(error.to_string()))?
        .map(|entry| entry.map_err(|error| Error::InvalidInput(error.to_string())))
        .collect()
}

#[cfg(test)]
const ASSIGNMENT_SIZE: usize = 21;

/// Writes one externally sorted primary-minimizer assignment per read.
/// The returned fixed-width file can be replayed without retaining the global
/// permutation or minimizer index in memory.
#[cfg(test)]
pub(crate) fn write_store_assignments(
    store: &mut PairStore,
    options: &ReorderOptions,
    temp_dir: impl AsRef<Path>,
    mate: usize,
) -> Result<std::path::PathBuf> {
    validate_options(options)?;
    if mate > 1 {
        return Err(Error::InvalidInput("mate index must be 0 or 1".into()));
    }
    let temp_dir = temp_dir.as_ref();
    std::fs::create_dir_all(temp_dir)?;
    let sorter = ExternalSorterBuilder::new()
        .with_chunk_size(options.sort_chunk_size.max(1))
        .with_tmp_dir(temp_dir)
        .with_compression(3)
        .num_threads(options.sort_threads.max(1))
        .build()
        .map_err(Error::Io)?;
    let count = store.count();
    let entries = (0..count).map(|read_id| {
        let entry = store
            .read_sequence(read_id, mate)
            .ok()
            .and_then(|sequence| {
                primary_minimizer(&sequence, options.kmer_size, options.window_size)
            });
        match entry {
            Some((hash, position, reverse)) => MinimizerEntry {
                hash,
                read_id,
                position,
                reverse,
            },
            None => MinimizerEntry {
                hash: u64::MAX,
                read_id,
                position: 0,
                reverse: false,
            },
        }
    });
    let sorted = sorter
        .sort(entries)
        .map_err(|error| Error::InvalidInput(error.to_string()))?;
    let path = temp_dir.join(format!("primary-assignments.{mate}.bin"));
    let mut output = File::create(&path)?;
    for entry in sorted {
        let entry = entry.map_err(|error| Error::InvalidInput(error.to_string()))?;
        output.write_all(&entry.hash.to_le_bytes())?;
        output.write_all(&entry.read_id.to_le_bytes())?;
        output.write_all(&entry.position.to_le_bytes())?;
        output.write_all(&[entry.reverse as u8])?;
    }
    output.flush()?;
    Ok(path)
}

#[cfg(test)]
pub(crate) struct AssignmentReader {
    file: File,
}

#[cfg(test)]
impl AssignmentReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self> {
        Ok(Self {
            file: OpenOptions::new().read(true).open(path)?,
        })
    }
}

#[cfg(test)]
impl Iterator for AssignmentReader {
    type Item = Result<MinimizerEntry>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut bytes = [0u8; ASSIGNMENT_SIZE];
        match self.file.read_exact(&mut bytes) {
            Ok(()) => Some(Ok(MinimizerEntry {
                hash: u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
                read_id: u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
                position: u32::from_le_bytes(bytes[16..20].try_into().unwrap()),
                reverse: bytes[20] != 0,
            })),
            Err(error) if error.kind() == std::io::ErrorKind::UnexpectedEof => None,
            Err(error) => Some(Err(Error::Io(error))),
        }
    }
}

pub(crate) fn validate_options(options: &ReorderOptions) -> Result<()> {
    if options.kmer_size == 0 || options.window_size == 0 {
        return Err(Error::InvalidInput(
            "minimizer k and window must be nonzero".into(),
        ));
    }
    if options.kmer_size > 63 {
        return Err(Error::InvalidInput("minimizer k cannot exceed 63".into()));
    }
    Ok(())
}

pub fn decode_consensus_read(cluster: &ConsensusCluster, read: &ConsensusRead) -> Result<Vec<u8>> {
    let start = read.offset as usize;
    let end = start
        .checked_add(read.length as usize)
        .ok_or_else(|| Error::InvalidFormat("consensus range overflow".into()))?;
    if end > cluster.reference.len() {
        return Err(Error::InvalidFormat("consensus range out of bounds".into()));
    }
    let mut sequence = cluster.reference[start..end].to_vec();
    for &(position, base) in &read.mismatches {
        let position = position as usize;
        if position >= sequence.len() {
            return Err(Error::InvalidFormat(
                "consensus mismatch out of bounds".into(),
            ));
        }
        sequence[position] = base;
    }
    if read.reverse {
        reverse_complement_in_place(&mut sequence);
    }
    Ok(sequence)
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizers_are_reverse_complement_canonical() {
        let sequence = b"ACGTTGCATGTCAGTACGAT";
        let k = 5;
        let reverse_sequence = sequence
            .iter()
            .rev()
            .map(|base| match base {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => unreachable!(),
            })
            .collect::<Vec<_>>();

        let forward = canonical_minimizers(sequence, k, 1);
        let mut reverse = canonical_minimizers(&reverse_sequence, k, 1)
            .into_iter()
            .map(|(hash, position, is_reverse)| {
                (
                    hash,
                    (sequence.len() - k - position as usize) as u32,
                    !is_reverse,
                )
            })
            .collect::<Vec<_>>();
        reverse.sort_by_key(|(_, position, _)| *position);
        assert_eq!(forward, reverse);
    }

    #[test]
    fn rolling_minimizers_match_naive_windows() {
        for sequence in [
            b"ACGTACGTACGTACGT".as_slice(),
            b"AAAAACCCCCGGGGGTTTTT",
            b"ACGTNNNNACGTACGTACGT",
        ] {
            for (k, window) in [(3, 4), (5, 7), (8, 3)] {
                assert_eq!(
                    canonical_minimizers(sequence, k, window),
                    naive_minimizers(sequence, k, window)
                );
            }
        }
    }

    fn naive_minimizers(sequence: &[u8], k: usize, window: usize) -> Vec<(u64, u32, bool)> {
        if sequence.len() < k || k == 0 || window == 0 {
            return Vec::new();
        }
        let kmers = sequence.len() - k + 1;
        let span = window.min(kmers);
        let values = (0..kmers)
            .map(|start| canonical_kmer(&sequence[start..start + k]))
            .collect::<Vec<_>>();
        let mut output = Vec::new();
        let mut previous = None;
        for start in 0..=kmers - span {
            let best = values
                .iter()
                .enumerate()
                .skip(start)
                .take(span)
                .filter_map(|(position, value)| {
                    value.map(|(hash, reverse)| (hash, position, reverse))
                })
                .min_by_key(|(hash, position, _)| (*hash, *position));
            if let Some((hash, position, reverse)) = best {
                if previous != Some((hash, position)) {
                    output.push((hash, position as u32, reverse));
                    previous = Some((hash, position));
                }
            }
        }
        output
    }

    #[test]
    fn consensus_residual_round_trip() {
        let cluster = ConsensusCluster {
            reference: b"ACGTAC".to_vec(),
            reads: Vec::new(),
        };
        let read = ConsensusRead {
            original_id: 0,
            offset: 1,
            reverse: false,
            length: 4,
            mismatches: vec![(2, b'N')],
        };
        assert_eq!(decode_consensus_read(&cluster, &read).unwrap(), b"CGNA");
    }
}
