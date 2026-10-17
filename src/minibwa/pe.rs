#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Extra {
    pub dp_score: i32,
    pub dp_max0: i32,
    pub dp_max: i32,
    pub dp_max2: i32,
    pub n_ambi: u32,
    pub n_cigar: i32,
    pub cigar: Vec<u32>,
}

impl Extra {
    #[inline(always)]
    pub(crate) fn add_n_ambi(&mut self, value: u32) {
        self.n_ambi = self.n_ambi.saturating_add(value);
    }

    /// The `n_cigar` leading CIGAR words.
    pub(crate) fn cigar(&self) -> &[u32] {
        let n = (self.n_cigar.max(0) as usize).min(self.cigar.len());
        &self.cigar[..n]
    }

    /// Mutable view of the `n_cigar` leading CIGAR words.
    pub(crate) fn cigar_mut(&mut self) -> &mut [u32] {
        let n = (self.n_cigar.max(0) as usize).min(self.cigar.len());
        &mut self.cigar[..n]
    }

    pub(crate) fn set_cigar_from_vec(&mut self, values: Vec<u32>) {
        self.n_cigar = values.len() as i32;
        self.cigar = values;
    }

    /// Appends CIGAR words after the current `n_cigar` prefix and raises
    /// `n_cigar` accordingly.
    pub(crate) fn extend_cigar_from_slice(&mut self, values: &[u32]) {
        if values.is_empty() {
            return;
        }
        let n = self.n_cigar.max(0) as usize;
        self.cigar.truncate(n);
        self.cigar.extend_from_slice(values);
        self.n_cigar = self.cigar.len() as i32;
    }

    pub(crate) fn remove_cigar(&mut self, index: usize) -> u32 {
        let len = self.n_cigar.max(0) as usize;
        assert!(index < len);
        self.n_cigar -= 1;
        self.cigar.remove(index)
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Hit {
    pub tid: i64,
    pub ts: i64,
    pub te: i64,
    pub id: i32,
    pub cnt: i32,
    pub score: i32,
    pub score0: i32,
    pub as_: i32,
    pub qs: i32,
    pub qe: i32,
    pub parent: i32,
    pub n_sub: i32,
    pub subsc: i32,
    pub mlen: i32,
    pub blen: i32,
    pub mapq: i32,
    pub hash: u32,
    pub flags: u32,
    pub p: Option<Box<Extra>>,
}

impl Hit {
    const REV: u32 = 1 << 0;
    const PROPER_PAIR: u32 = 1 << 1;
    const SAM_PRI: u32 = 1 << 2;
    const FLT: u32 = 1 << 3;
    const INV: u32 = 1 << 4;
    const SPLIT_SHIFT: u32 = 5;
    const SPLIT_MASK: u32 = 0b11 << Self::SPLIT_SHIFT;
    const SPLIT_INV: u32 = 1 << 7;
    const FRAC_HIGH_SHIFT: u32 = 9;
    const FRAC_HIGH_MASK: u32 = 0xff << Self::FRAC_HIGH_SHIFT;

    #[inline(always)]
    fn bit(&self, mask: u32) -> u8 {
        ((self.flags & mask) != 0) as u8
    }

    #[inline(always)]
    fn set_bit(&mut self, mask: u32, value: u8) {
        if value != 0 {
            self.flags |= mask;
        } else {
            self.flags &= !mask;
        }
    }

    #[inline(always)]
    pub(crate) fn rev(&self) -> u8 {
        self.bit(Self::REV)
    }

    #[inline(always)]
    pub(crate) fn set_rev(&mut self, value: u8) {
        self.set_bit(Self::REV, value);
    }

    #[inline(always)]
    pub(crate) fn proper_pair(&self) -> u8 {
        self.bit(Self::PROPER_PAIR)
    }

    #[inline(always)]
    pub(crate) fn sam_pri(&self) -> u8 {
        self.bit(Self::SAM_PRI)
    }

    #[inline(always)]
    pub(crate) fn set_sam_pri(&mut self, value: u8) {
        self.set_bit(Self::SAM_PRI, value);
    }

    #[inline(always)]
    pub(crate) fn flt(&self) -> u8 {
        self.bit(Self::FLT)
    }

    #[inline(always)]
    pub(crate) fn inv(&self) -> u8 {
        self.bit(Self::INV)
    }

    #[inline(always)]
    pub(crate) fn set_inv(&mut self, value: u8) {
        self.set_bit(Self::INV, value);
    }

    #[inline(always)]
    pub(crate) fn split(&self) -> u8 {
        ((self.flags & Self::SPLIT_MASK) >> Self::SPLIT_SHIFT) as u8
    }

    #[inline(always)]
    pub(crate) fn set_split(&mut self, value: u8) {
        self.flags = (self.flags & !Self::SPLIT_MASK)
            | (((value as u32) << Self::SPLIT_SHIFT) & Self::SPLIT_MASK);
    }

    #[inline(always)]
    pub(crate) fn split_inv(&self) -> u8 {
        self.bit(Self::SPLIT_INV)
    }

    #[inline(always)]
    pub(crate) fn set_split_inv(&mut self, value: u8) {
        self.set_bit(Self::SPLIT_INV, value);
    }

    #[inline(always)]
    pub(crate) fn set_frac_high(&mut self, value: u8) {
        self.flags =
            (self.flags & !Self::FRAC_HIGH_MASK) | ((value as u32) << Self::FRAC_HIGH_SHIFT);
    }
}
