use crate::star::index::packed_array::PackedArray;

/// Suffix array for genome indexing.
///
/// Stores positions in the genome (forward + reverse complement) sorted
/// by their suffixes, enabling fast exact match search.
#[derive(Clone)]
pub struct SuffixArray {
    /// Packed array of suffix positions (with strand bit)
    pub data: PackedArray,

    /// Strand bit position (position where bit distinguishes forward/reverse)
    pub gstrand_bit: u32,

    /// Mask for extracting position (all bits below gstrand_bit)
    pub gstrand_mask: u64,
}

impl SuffixArray {
    /// Calculate GstrandBit from genome size.
    ///
    /// Formula: max(32, floor(log2(nGenome)) + 1)
    pub fn calculate_gstrand_bit(n_genome: u64) -> u32 {
        if n_genome == 0 {
            return 32;
        }
        let log2_bits = 64 - n_genome.leading_zeros();
        u32::max(32, log2_bits)
    }

    /// Get the number of suffixes in the array.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Decode a suffix array entry into (position, is_reverse).
    pub fn decode(&self, sa_value: u64) -> (u64, bool) {
        let is_reverse = (sa_value >> self.gstrand_bit) != 0;
        let position = sa_value & self.gstrand_mask;
        (position, is_reverse)
    }

    /// Read a suffix array entry.
    pub fn get(&self, index: usize) -> u64 {
        self.data.read(index)
    }

    /// Software-prefetch the SA entry at `index` (hint only; see
    /// [`crate::star::index::packed_array::PackedArray::prefetch`]). Used to overlap the
    /// random SA loads of the binary search with the current probe's genome compare.
    #[inline]
    pub fn prefetch(&self, index: usize) {
        self.data.prefetch(index);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gstrand_bit_calculation() {
        assert_eq!(SuffixArray::calculate_gstrand_bit(0), 32);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1), 32);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1000), 32);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1u64 << 32), 33);
        assert_eq!(SuffixArray::calculate_gstrand_bit(1u64 << 33), 34);
    }
}
