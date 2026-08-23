use rkyv::{Archive, Deserialize, Serialize};

use crate::{Error, Result};

#[derive(Archive, Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub(crate) struct PackedSequence {
    length: u32,
    bytes: Vec<u8>,
    packed: bool,
}

#[allow(dead_code)]
impl PackedSequence {
    pub(crate) fn try_from_bytes(sequence: &[u8]) -> Result<Self> {
        let length = u32::try_from(sequence.len())
            .map_err(|_| Error::InvalidInput("sequence exceeds 4 GiB".into()))?;
        let packed = sequence
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'));

        let bytes = if packed {
            let mut bytes = vec![0; sequence.len().div_ceil(4)];
            for (index, base) in sequence.iter().enumerate() {
                let bits = match base {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => unreachable!(),
                };
                bytes[index / 4] |= bits << (6 - (index % 4) * 2);
            }
            bytes
        } else {
            sequence.to_vec()
        };

        Ok(Self {
            length,
            bytes,
            packed,
        })
    }

    #[cfg(test)]
    pub(crate) fn from_bytes(sequence: &[u8]) -> Self {
        Self::try_from_bytes(sequence).unwrap()
    }

    pub(crate) fn from_encoded(length: u32, packed: bool, bytes: Vec<u8>) -> Result<Self> {
        let expected = if packed {
            (length as usize).div_ceil(4)
        } else {
            length as usize
        };
        if bytes.len() != expected {
            return Err(Error::InvalidFormat(
                "packed sequence body length mismatch".into(),
            ));
        }
        Ok(Self {
            length,
            bytes,
            packed,
        })
    }

    pub(crate) fn len(&self) -> u32 {
        self.length
    }

    pub(crate) fn is_packed(&self) -> bool {
        self.packed
    }

    pub(crate) fn encoded_bytes(&self) -> &[u8] {
        &self.bytes
    }

    pub(crate) fn decode(&self) -> Vec<u8> {
        if !self.packed {
            return self.bytes.clone();
        }

        const BASES: &[u8; 4] = b"ACGT";
        (0..self.length as usize)
            .map(|index| BASES[((self.bytes[index / 4] >> (6 - (index % 4) * 2)) & 0b11) as usize])
            .collect()
    }

    pub(crate) fn encoded_size(&self) -> usize {
        self.bytes.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packs_clean_sequences_including_partial_final_bytes() {
        let sequence = b"ACGTACGTACGTA";
        for length in 0..=sequence.len() {
            let packed = PackedSequence::from_bytes(&sequence[..length]);
            assert!(packed.packed);
            assert_eq!(packed.encoded_size(), length.div_ceil(4));
            assert_eq!(packed.decode(), &sequence[..length]);
        }
    }

    #[test]
    fn uses_two_bits_per_base_in_archive_order() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"A", &[0b0000_0000]),
            (b"AC", &[0b0001_0000]),
            (b"ACG", &[0b0001_1000]),
            (b"ACGT", &[0b0001_1011]),
            (b"ACGTA", &[0b0001_1011, 0b0000_0000]),
        ];
        for (sequence, expected) in cases {
            let packed = PackedSequence::from_bytes(sequence);
            assert_eq!(packed.bytes, *expected);
            assert_eq!(packed.decode(), *sequence);
        }
    }

    #[test]
    fn arbitrary_symbols_round_trip_through_raw_fallback() {
        let cases: &[&[u8]] = &[b"acgt", b"ACGTNRYKMSWBDHV", b"ACGT-.*", b"A\0C\xffT"];
        for sequence in cases {
            let packed = PackedSequence::from_bytes(sequence);
            assert!(!packed.packed);
            assert_eq!(packed.encoded_size(), sequence.len());
            assert_eq!(packed.decode(), *sequence);
        }
    }

    #[test]
    fn every_byte_value_round_trips_exactly() {
        let sequence = (0..=u8::MAX).collect::<Vec<_>>();
        let packed = PackedSequence::from_bytes(&sequence);
        assert!(!packed.packed);
        assert_eq!(packed.decode(), sequence);
    }
}
