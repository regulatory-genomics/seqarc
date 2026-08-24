//! Lossless encoding for a collection of index reads using a predicted whitelist.
//!
//! This module is intentionally self-contained. It is not registered by the crate
//! yet, so it can be developed and tested independently of the archive format.

use std::collections::BTreeMap;
use std::fmt;

const MAGIC: &[u8; 8] = b"IDXREAD\0";
const VERSION: u8 = 1;
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

const MODE_EXACT: u8 = 0;
const MODE_ONE_MISMATCH: u8 = 1;
const MODE_TWO_MISMATCHES: u8 = 2;
const MODE_RAW: u8 = 3;

/// An error encountered while encoding or decoding index reads.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum IndexReadError {
    InvalidInput(String),
    InvalidFormat(String),
}

impl fmt::Display for IndexReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(message) => write!(f, "invalid index reads: {message}"),
            Self::InvalidFormat(message) => write!(f, "invalid index-read data: {message}"),
        }
    }
}

impl std::error::Error for IndexReadError {}

/// Predict a whitelist and encode all reads against it with up to two mismatches.
///
/// The returned byte string contains the predicted whitelist and can therefore be
/// decoded without any external data. Reads that cannot be represented against a
/// whitelist entry are stored raw. Input order is preserved.
pub fn encode_index_reads<T: AsRef<[u8]>>(reads: &[T]) -> Result<Vec<u8>, IndexReadError> {
    let whitelist = Whitelist::predict(reads)?;
    encode_with_whitelist(reads, &whitelist)
}

/// Decode index reads previously produced by [`encode_index_reads`].
pub fn decode_index_reads(encoded: &[u8]) -> Result<Vec<Vec<u8>>, IndexReadError> {
    let mut input = Input::new(encoded);
    if input.read_exact(MAGIC.len())? != MAGIC {
        return Err(invalid_format("invalid magic"));
    }
    if input.read_byte()? != VERSION {
        return Err(invalid_format("unsupported version"));
    }

    let read_count = input.read_usize_varint("read count")?;
    let whitelist_count = input.read_usize_varint("whitelist count")?;
    if whitelist_count > input.remaining() {
        return Err(invalid_format("whitelist count exceeds remaining data"));
    }

    let mut entries = Vec::with_capacity(whitelist_count);
    let mut seen = BTreeMap::new();
    for ordinal in 0..whitelist_count {
        let length = input.read_usize_varint("whitelist sequence length")?;
        let sequence = input.read_exact(length)?.to_vec();
        if sequence.len() <= 1 || !is_canonical(&sequence) {
            return Err(invalid_format(
                "whitelist entries must be canonical sequences longer than one base",
            ));
        }
        if seen.insert(sequence.clone(), ordinal).is_some() {
            return Err(invalid_format("duplicate whitelist entry"));
        }
        entries.push(WhitelistEntry { sequence, count: 0 });
    }

    if read_count > input.remaining() {
        return Err(invalid_format("read count exceeds remaining data"));
    }
    let mut reads = Vec::with_capacity(read_count);
    for _ in 0..read_count {
        let mode = input.read_byte()?;
        if mode == MODE_RAW {
            let length = input.read_usize_varint("raw sequence length")?;
            reads.push(input.read_exact(length)?.to_vec());
            continue;
        }

        let ordinal = input.read_usize_varint("whitelist ordinal")?;
        let entry = entries
            .get(ordinal)
            .ok_or_else(|| invalid_format("whitelist ordinal is out of bounds"))?;
        let mismatch_count = match mode {
            MODE_EXACT => 0,
            MODE_ONE_MISMATCH => 1,
            MODE_TWO_MISMATCHES => 2,
            _ => return Err(invalid_format("unknown record mode")),
        };
        let mut sequence = entry.sequence.clone();
        let mut previous_position = None;
        for _ in 0..mismatch_count {
            let position = input.read_usize_varint("mismatch position")?;
            if position >= sequence.len() {
                return Err(invalid_format("mismatch position is out of bounds"));
            }
            if previous_position.is_some_and(|previous| position <= previous) {
                return Err(invalid_format(
                    "mismatch positions are not strictly increasing",
                ));
            }
            let original = input.read_byte()?;
            if original == sequence[position] {
                return Err(invalid_format(
                    "mismatch does not change the whitelist base",
                ));
            }
            sequence[position] = original;
            previous_position = Some(position);
        }
        reads.push(sequence);
    }

    if input.remaining() != 0 {
        return Err(invalid_format("trailing bytes"));
    }
    Ok(reads)
}

#[derive(Clone, Debug)]
struct WhitelistEntry {
    sequence: Vec<u8>,
    count: usize,
}

#[derive(Clone, Debug)]
struct Whitelist {
    entries: Vec<WhitelistEntry>,
    ordinals: BTreeMap<Vec<u8>, usize>,
}

impl Whitelist {
    fn predict<T: AsRef<[u8]>>(reads: &[T]) -> Result<Self, IndexReadError> {
        let mut counts_by_length: BTreeMap<usize, BTreeMap<Vec<u8>, usize>> = BTreeMap::new();
        for read in reads {
            let sequence = read.as_ref();
            if sequence.len() <= 1 || !is_canonical(sequence) {
                continue;
            }
            let count = counts_by_length
                .entry(sequence.len())
                .or_default()
                .entry(sequence.to_vec())
                .or_default();
            *count = count
                .checked_add(1)
                .ok_or_else(|| invalid_input("sequence count overflow"))?;
        }

        let mut entries = Vec::new();
        for counts in counts_by_length.into_values() {
            let threshold = ordmag_threshold(counts.values().copied());
            let mut predicted = counts
                .into_iter()
                .filter(|(_, count)| *count >= threshold)
                .map(|(sequence, count)| WhitelistEntry { sequence, count })
                .collect::<Vec<_>>();
            predicted.sort_by(|left, right| {
                right
                    .count
                    .cmp(&left.count)
                    .then_with(|| left.sequence.cmp(&right.sequence))
            });
            entries.extend(predicted);
        }
        entries.sort_by(|left, right| {
            right
                .count
                .cmp(&left.count)
                .then_with(|| left.sequence.cmp(&right.sequence))
        });

        Self::from_entries(entries)
    }

    fn from_entries(entries: Vec<WhitelistEntry>) -> Result<Self, IndexReadError> {
        let mut ordinals = BTreeMap::new();
        for (ordinal, entry) in entries.iter().enumerate() {
            if entry.sequence.len() <= 1 || !is_canonical(&entry.sequence) {
                return Err(invalid_input(
                    "whitelist entries must be canonical sequences longer than one base",
                ));
            }
            if ordinals.insert(entry.sequence.clone(), ordinal).is_some() {
                return Err(invalid_input("duplicate whitelist entry"));
            }
        }
        Ok(Self { entries, ordinals })
    }

    fn candidate(&self, sequence: &[u8]) -> Option<Candidate> {
        if sequence.len() <= 1 || !is_canonical(sequence) {
            return None;
        }
        if let Some(&ordinal) = self.ordinals.get(sequence) {
            return Some(Candidate {
                ordinal,
                mismatches: Vec::new(),
            });
        }

        let mut mutated = sequence.to_vec();
        let mut best = None;
        for first_position in 0..mutated.len() {
            let first_original = mutated[first_position];
            for first_base in BASES {
                if first_base == first_original {
                    continue;
                }
                mutated[first_position] = first_base;
                if let Some(&ordinal) = self.ordinals.get(&mutated) {
                    update_candidate(
                        &mut best,
                        Candidate {
                            ordinal,
                            mismatches: vec![(first_position, first_original)],
                        },
                        &self.entries,
                    );
                }
                for second_position in (first_position + 1)..mutated.len() {
                    let second_original = mutated[second_position];
                    for second_base in BASES {
                        if second_base == second_original {
                            continue;
                        }
                        mutated[second_position] = second_base;
                        if let Some(&ordinal) = self.ordinals.get(&mutated) {
                            update_candidate(
                                &mut best,
                                Candidate {
                                    ordinal,
                                    mismatches: vec![
                                        (first_position, first_original),
                                        (second_position, second_original),
                                    ],
                                },
                                &self.entries,
                            );
                        }
                    }
                    mutated[second_position] = second_original;
                }
            }
            mutated[first_position] = first_original;
        }
        best
    }
}

#[derive(Clone, Debug)]
struct Candidate {
    ordinal: usize,
    mismatches: Vec<(usize, u8)>,
}

fn update_candidate(
    best: &mut Option<Candidate>,
    candidate: Candidate,
    entries: &[WhitelistEntry],
) {
    let replace = best.as_ref().is_none_or(|current| {
        let candidate_entry = &entries[candidate.ordinal];
        let current_entry = &entries[current.ordinal];
        candidate
            .mismatches
            .len()
            .cmp(&current.mismatches.len())
            .then_with(|| current_entry.count.cmp(&candidate_entry.count))
            .then_with(|| candidate_entry.sequence.cmp(&current_entry.sequence))
            .is_lt()
    });
    if replace {
        *best = Some(candidate);
    }
}

fn encode_with_whitelist<T: AsRef<[u8]>>(
    reads: &[T],
    whitelist: &Whitelist,
) -> Result<Vec<u8>, IndexReadError> {
    let mut output = Vec::new();
    output.extend(MAGIC);
    output.push(VERSION);
    write_usize_varint(&mut output, reads.len())?;
    write_usize_varint(&mut output, whitelist.entries.len())?;
    for entry in &whitelist.entries {
        write_usize_varint(&mut output, entry.sequence.len())?;
        output.extend(&entry.sequence);
    }

    for read in reads {
        let sequence = read.as_ref();
        if let Some(candidate) = whitelist.candidate(sequence) {
            output.push(match candidate.mismatches.len() {
                0 => MODE_EXACT,
                1 => MODE_ONE_MISMATCH,
                2 => MODE_TWO_MISMATCHES,
                _ => unreachable!(),
            });
            write_usize_varint(&mut output, candidate.ordinal)?;
            for (position, original) in candidate.mismatches {
                write_usize_varint(&mut output, position)?;
                output.push(original);
            }
        } else {
            output.push(MODE_RAW);
            write_usize_varint(&mut output, sequence.len())?;
            output.extend(sequence);
        }
    }
    Ok(output)
}

fn ordmag_threshold(counts: impl IntoIterator<Item = usize>) -> usize {
    let mut sorted = counts
        .into_iter()
        .filter(|count| *count > 0)
        .collect::<Vec<_>>();
    if sorted.is_empty() {
        return usize::MAX;
    }
    sorted.sort_unstable();
    if sorted.len() <= 2 {
        return sorted[0];
    }

    let recovered = estimate_recovered_ordmag(&sorted).clamp(1, sorted.len());
    let baseline_index =
        ((recovered as f64 * 0.01).round() as usize).min(sorted.len().saturating_sub(1));
    let top_n = find_within_ordmag(&sorted, baseline_index).clamp(1, sorted.len());
    sorted[sorted.len().saturating_sub(top_n.saturating_add(1))]
}

fn estimate_recovered_ordmag(sorted_counts: &[usize]) -> usize {
    let mut expected = 1usize;
    let mut best = (expected, f64::INFINITY);
    for _ in 0..18 {
        let baseline_index =
            ((expected as f64 * 0.01).round() as usize).min(sorted_counts.len().saturating_sub(1));
        let observed = find_within_ordmag(sorted_counts, baseline_index);
        let difference = observed.abs_diff(expected) as f64;
        let loss = difference * difference / expected as f64;
        if loss < best.1 {
            best = (expected, loss);
        }
        expected = expected.saturating_mul(2);
    }
    best.0
}

fn find_within_ordmag(sorted_counts: &[usize], baseline_index: usize) -> usize {
    let baseline = sorted_counts[sorted_counts.len() - baseline_index - 1];
    let cutoff = ((baseline as f64 * 0.1).round() as usize).max(1);
    let first = sorted_counts.partition_point(|count| *count < cutoff);
    sorted_counts.len() - first
}

fn is_canonical(sequence: &[u8]) -> bool {
    sequence.iter().all(|base| BASES.contains(base))
}

fn write_usize_varint(output: &mut Vec<u8>, value: usize) -> Result<(), IndexReadError> {
    let value = u64::try_from(value).map_err(|_| invalid_input("value exceeds u64"))?;
    write_varint(output, value);
    Ok(())
}

fn write_varint(output: &mut Vec<u8>, mut value: u64) {
    loop {
        let mut byte = (value & 0x7f) as u8;
        value >>= 7;
        if value != 0 {
            byte |= 0x80;
        }
        output.push(byte);
        if value == 0 {
            break;
        }
    }
}

struct Input<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Input<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    fn remaining(&self) -> usize {
        self.bytes.len() - self.position
    }

    fn read_byte(&mut self) -> Result<u8, IndexReadError> {
        let byte = self
            .bytes
            .get(self.position)
            .copied()
            .ok_or_else(|| invalid_format("unexpected end of data"))?;
        self.position += 1;
        Ok(byte)
    }

    fn read_exact(&mut self, length: usize) -> Result<&'a [u8], IndexReadError> {
        let end = self
            .position
            .checked_add(length)
            .ok_or_else(|| invalid_format("length overflow"))?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or_else(|| invalid_format("unexpected end of data"))?;
        self.position = end;
        Ok(value)
    }

    fn read_usize_varint(&mut self, name: &str) -> Result<usize, IndexReadError> {
        let value = self.read_varint(name)?;
        usize::try_from(value).map_err(|_| invalid_format(format!("{name} exceeds usize")))
    }

    fn read_varint(&mut self, name: &str) -> Result<u64, IndexReadError> {
        let mut value = 0u64;
        for shift in (0..=63).step_by(7) {
            let byte = self.read_byte()?;
            let payload = (byte & 0x7f) as u64;
            if shift == 63 && payload > 1 {
                return Err(invalid_format(format!("{name} varint overflow")));
            }
            value |= payload << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
        }
        Err(invalid_format(format!("{name} varint is too long")))
    }
}

fn invalid_input(message: impl Into<String>) -> IndexReadError {
    IndexReadError::InvalidInput(message.into())
}

fn invalid_format(message: impl Into<String>) -> IndexReadError {
    IndexReadError::InvalidFormat(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn whitelist(entries: &[(&[u8], usize)]) -> Whitelist {
        Whitelist::from_entries(
            entries
                .iter()
                .map(|(sequence, count)| WhitelistEntry {
                    sequence: sequence.to_vec(),
                    count: *count,
                })
                .collect(),
        )
        .unwrap()
    }

    fn round_trip(reads: &[&[u8]]) {
        let encoded = encode_index_reads(reads).unwrap();
        assert_eq!(decode_index_reads(&encoded).unwrap(), reads);
    }

    #[test]
    fn empty_input_round_trips() {
        round_trip(&[]);
    }

    #[test]
    fn exact_reads_round_trip_in_order() {
        round_trip(&[b"ACGT", b"TGCA", b"ACGT", b"CCCC"]);
    }

    #[test]
    fn one_mismatch_is_stored_losslessly() {
        let reads = [b"ACGA".as_slice()];
        let encoded = encode_with_whitelist(&reads, &whitelist(&[(b"ACGT", 10)])).unwrap();
        assert_eq!(decode_index_reads(&encoded).unwrap(), reads);
        assert_eq!(*encoded.last().unwrap(), b'A');
    }

    #[test]
    fn two_mismatches_are_stored_losslessly() {
        let reads = [b"AAGT".as_slice()];
        let encoded = encode_with_whitelist(&reads, &whitelist(&[(b"CCGT", 10)])).unwrap();
        assert_eq!(decode_index_reads(&encoded).unwrap(), reads);
    }

    #[test]
    fn more_than_two_mismatches_use_raw_mode() {
        let reads = [b"TTTT".as_slice()];
        let encoded = encode_with_whitelist(&reads, &whitelist(&[(b"AAAA", 10)])).unwrap();
        assert_eq!(decode_index_reads(&encoded).unwrap(), reads);
        let records_start = MAGIC.len() + 1 + 1 + 1 + 1 + 4;
        assert_eq!(encoded[records_start], MODE_RAW);
    }

    #[test]
    fn noncanonical_and_short_reads_round_trip_raw() {
        round_trip(&[b"ACNT", b"acgt", b"A", b""]);
    }

    #[test]
    fn mixed_lengths_round_trip() {
        round_trip(&[b"ACGT", b"ACGTA", b"TGCA", b"ACGTA"]);
    }

    #[test]
    fn prediction_and_encoding_are_deterministic() {
        let reads = [b"AAAA", b"CCCC", b"AAAA", b"AAAT", b"CCCC", b"AAAA"];
        let first = encode_index_reads(&reads).unwrap();
        let second = encode_index_reads(&reads).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn candidate_prefers_fewer_mismatches() {
        let whitelist = whitelist(&[(b"AAAT", 1), (b"AACC", 100)]);
        let candidate = whitelist.candidate(b"AAAA").unwrap();
        assert_eq!(whitelist.entries[candidate.ordinal].sequence, b"AAAT");
    }

    #[test]
    fn candidate_prefers_higher_count_then_lexicographic_order() {
        let whitelist = whitelist(&[(b"AAAT", 10), (b"AACA", 20), (b"AAGA", 20)]);
        let candidate = whitelist.candidate(b"AAAA").unwrap();
        assert_eq!(whitelist.entries[candidate.ordinal].sequence, b"AACA");
    }

    #[test]
    fn low_diversity_prediction_does_not_fail() {
        let reads = vec![b"ACGT".as_slice(); 20];
        let predicted = Whitelist::predict(&reads).unwrap();
        assert_eq!(predicted.entries.len(), 1);
        assert_eq!(predicted.entries[0].sequence, b"ACGT");
    }

    #[test]
    fn malformed_data_is_rejected() {
        let encoded = encode_index_reads(&[b"ACGT"]).unwrap();

        let mut bad_magic = encoded.clone();
        bad_magic[0] ^= 1;
        assert!(decode_index_reads(&bad_magic).is_err());

        assert!(decode_index_reads(&encoded[..encoded.len() - 1]).is_err());

        let mut trailing = encoded;
        trailing.push(0);
        assert!(decode_index_reads(&trailing).is_err());
    }

    #[test]
    fn invalid_mismatch_payload_is_rejected() {
        let reads = [b"ACGA".as_slice()];
        let mut encoded = encode_with_whitelist(&reads, &whitelist(&[(b"ACGT", 10)])).unwrap();
        *encoded.last_mut().unwrap() = b'T';
        assert!(decode_index_reads(&encoded).is_err());
    }

    #[test]
    fn invalid_mode_and_ordinal_are_rejected() {
        let reads = [b"ACGT".as_slice()];
        let encoded = encode_with_whitelist(&reads, &whitelist(&[(b"ACGT", 10)])).unwrap();
        let records_start = MAGIC.len() + 1 + 1 + 1 + 1 + 4;

        let mut invalid_mode = encoded.clone();
        invalid_mode[records_start] = 4;
        assert!(decode_index_reads(&invalid_mode).is_err());

        let mut invalid_ordinal = encoded;
        invalid_ordinal[records_start + 1] = 1;
        assert!(decode_index_reads(&invalid_ordinal).is_err());
    }

    #[test]
    fn varint_overflow_is_rejected() {
        let mut encoded = Vec::from(MAGIC.as_slice());
        encoded.push(VERSION);
        encoded.extend([0xff; 10]);
        encoded.push(0x02);
        assert!(decode_index_reads(&encoded).is_err());
    }
}
