#![allow(unused_variables)]

use crate::minibwa::qsufsort::suffix_sort;
use std::io::Seek;

pub(crate) type BgInt = u64;
pub(crate) type SBgInt = i64;

pub(crate) const ALPHABET_SIZE: usize = 4;
pub(crate) const BIT_PER_CHAR: u32 = 2;
pub(crate) const CHAR_PER_WORD: u64 = 16;
pub(crate) const CHAR_PER_BYTE: u32 = 4;
pub(crate) const BITS_IN_WORD: u32 = 32;
pub(crate) const BITS_IN_BYTE: u32 = 8;
pub(crate) const ALL_ONE_MASK: u32 = 0xffff_ffff;
pub(crate) const DNA_OCC_CNT_TABLE_SIZE_IN_WORD: usize = 65536;
pub(crate) const OCC_VALUE_PER_WORD: u64 = 2;
pub(crate) const OCC_INTERVAL: u64 = 256;
pub(crate) const OCC_INTERVAL_MAJOR: u64 = 65536;
pub(crate) const MIN_AVAILABLE_WORD: u64 = 0x10000;
pub(crate) const BWTINC_INSERT_SORT_NUM_ITEM: i64 = 7;

fn c_path(path: &std::path::Path) -> std::path::PathBuf {
    let lossy = path.to_string_lossy();
    let end = lossy
        .as_bytes()
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(lossy.len());
    std::path::PathBuf::from(lossy[..end].to_string())
}

fn c_strerror(errno: i32) -> String {
    std::io::Error::from_raw_os_error(errno)
        .to_string()
        .trim_end_matches(&format!(" (os error {errno})"))
        .to_string()
}

fn einval() -> i32 {
    #[cfg(not(target_arch = "wasm32"))]
    {
        libc::EINVAL
    }
    #[cfg(target_arch = "wasm32")]
    {
        22
    }
}

fn bwt_save_error_and_exit(path: &std::path::Path, err: std::io::Error) -> ! {
    let reason = err
        .raw_os_error()
        .map(c_strerror)
        .unwrap_or_else(|| err.to_string());
    eprintln!(
        "bwt_save_bwt_code_and_occ(): Error writing to {} : {}",
        path.display(),
        reason
    );
    std::process::exit(1);
}

fn bwt_packed_seek_error_and_exit(path: &std::path::Path) -> ! {
    eprintln!(
        "bwt_inc_construct_from_packed() : Can't seek on {} : {}",
        path.display(),
        c_strerror(einval())
    );
    std::process::exit(1);
}

fn bwt_packed_read_error_and_exit(path: &std::path::Path, err: std::io::Error) -> ! {
    let reason = err
        .raw_os_error()
        .map(c_strerror)
        .unwrap_or_else(|| err.to_string());
    eprintln!(
        "bwt_inc_construct_from_packed() : Can't read from {} : {}",
        path.display(),
        reason
    );
    std::process::exit(1);
}

fn bwt_packed_seek_back(fp: &mut std::fs::File, path: &std::path::Path, n: u64) {
    let Ok(delta) = i64::try_from(n) else {
        bwt_packed_seek_error_and_exit(path);
    };
    if fp.seek(std::io::SeekFrom::Current(-delta)).is_err() {
        bwt_packed_seek_error_and_exit(path);
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct GenBwt {
    pub text_length: BgInt,
    pub inverse_sa0: BgInt,
    pub cumulative_freq: Vec<BgInt>,
    pub bwt_code: Vec<u32>,
    pub occ_value: Vec<u32>,
    pub occ_value_major: Vec<BgInt>,
    pub decode_table: Vec<u32>,
    pub bwt_size_in_word: BgInt,
    pub occ_size_in_word: BgInt,
    pub occ_major_size_in_word: BgInt,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct BwtInc {
    pub bwt: GenBwt,
    pub number_of_iteration_done: u32,
    pub cumulative_count_in_current_build: Vec<BgInt>,
    pub available_word: BgInt,
    pub build_size: BgInt,
    pub initial_max_build_size: BgInt,
    pub inc_max_build_size: BgInt,
    pub first_char_in_last_iteration: u32,
    pub working_memory: Vec<u32>,
    pub packed_text: Vec<u32>,
    pub text_buffer: Vec<u8>,
    pub packed_shift: Vec<u32>,
    pub packed_text_offset: usize,
    pub text_buffer_offset: usize,
}

/// Original C static function `text_length_from_byte_packed` from `minibwa/bwtgen.c:97`.
pub(crate) fn text_length_from_byte_packed(
    byte_packed_length: BgInt,
    bit_per_char: u32,
    last_byte_length: u32,
) -> BgInt {
    (byte_packed_length - 1) * (BITS_IN_BYTE / bit_per_char) as u64 + last_byte_length as u64
}

/// Original C static function `initialize_val` from `minibwa/bwtgen.c:103`.
pub(crate) fn initialize_val(start_addr: &mut [u32], length: BgInt, init_value: u32) {
    for x in start_addr.iter_mut().take(length as usize) {
        *x = init_value;
    }
}

/// Original C static function `initialize_val_bg` from `minibwa/bwtgen.c:109`.
pub(crate) fn initialize_val_bg(start_addr: &mut [BgInt], length: BgInt, init_value: BgInt) {
    for x in start_addr.iter_mut().take(length as usize) {
        *x = init_value;
    }
}

/// Reinterprets a `BgInt` (u64) slice as `SBgInt` (i64) in place. The sort
/// helpers in `qsufsort` compare keys as signed integers, while buffers are
/// stored as u64. Both types have identical layout; values are only read
/// through the signed view.
#[inline]
fn as_qsint_slice_mut(values: &mut [BgInt]) -> &mut [SBgInt] {
    // SAFETY: `BgInt` (u64) and `SBgInt` (i64) have the same size and alignment.
    unsafe { std::slice::from_raw_parts_mut(values.as_mut_ptr() as *mut SBgInt, values.len()) }
}

/// Allocates a `Vec<T>` of `len` **uninitialized** elements. Every element
/// must be written before it is read; used for large index-build scratch
/// buffers where zero-initialization would double the build cost.
fn uninit_vec<T>(len: usize) -> Vec<T> {
    let mut values = Vec::<std::mem::MaybeUninit<T>>::with_capacity(len);
    unsafe {
        // SAFETY: callers guarantee every element is initialized before read.
        values.set_len(len);
        let ptr = values.as_mut_ptr() as *mut T;
        let cap = values.capacity();
        std::mem::forget(values);
        Vec::from_raw_parts(ptr, len, cap)
    }
}

/// Original C static function `generate_dna_occ_count_table` from `minibwa/bwtgen.c:115`.
pub(crate) fn generate_dna_occ_count_table(dna_decode_table: &mut [u32]) {
    for i in 0..DNA_OCC_CNT_TABLE_SIZE_IN_WORD {
        let mut c = i as u32;
        dna_decode_table[i] = 0;
        for _ in 0..8 {
            let t = c & 3;
            dna_decode_table[i] = dna_decode_table[i].wrapping_add(1 << (t * 8));
            c >>= 2;
        }
    }
}

/// Original C static function `bwt_occ_value_major_size_in_word` from `minibwa/bwtgen.c:131`.
pub(crate) fn bwt_occ_value_major_size_in_word(num_char: BgInt) -> BgInt {
    let num_of_occ_value = (num_char + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    let num_of_occ_interval_per_major = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
    (num_of_occ_value + num_of_occ_interval_per_major - 1) / num_of_occ_interval_per_major * ALPHABET_SIZE as u64
}

/// Original C static function `bwt_occ_value_minor_size_in_word` from `minibwa/bwtgen.c:140`.
pub(crate) fn bwt_occ_value_minor_size_in_word(num_char: BgInt) -> BgInt {
    let num_of_occ_value = (num_char + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    (num_of_occ_value + OCC_VALUE_PER_WORD - 1) / OCC_VALUE_PER_WORD * ALPHABET_SIZE as u64
}

/// Original C static function `bwt_resident_size_in_word` from `minibwa/bwtgen.c:147`.
pub(crate) fn bwt_resident_size_in_word(num_char: BgInt) -> BgInt {
    let num_char_round_up_to_occ_interval = (num_char + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;
    (num_char_round_up_to_occ_interval + CHAR_PER_WORD - 1) / CHAR_PER_WORD
}

/// Original C static function `bwt_inc_set_build_size_and_text_addr` from `minibwa/bwtgen.c:158`.
pub(crate) fn bwt_inc_set_build_size_and_text_addr(bwt_inc: &mut BwtInc) {
    let word_scale = (std::mem::size_of::<BgInt>() / 4) as u64;
    let max_build_size;
    if bwt_inc.bwt.text_length == 0 {
        max_build_size = (bwt_inc.available_word - (2 + OCC_INTERVAL / CHAR_PER_WORD) * word_scale)
            / (2 * CHAR_PER_WORD + 1)
            * CHAR_PER_WORD
            / word_scale;
        if bwt_inc.initial_max_build_size > 0 {
            bwt_inc.build_size = bwt_inc.initial_max_build_size.min(max_build_size);
        } else {
            bwt_inc.build_size = max_build_size;
        }
    } else {
        max_build_size = (bwt_inc.available_word
            - bwt_inc.bwt.bwt_size_in_word
            - bwt_inc.bwt.occ_size_in_word
            - (3 + bwt_inc.number_of_iteration_done as u64 * OCC_INTERVAL / BIT_PER_CHAR as u64)
                * word_scale)
            / 3
            / word_scale;
        if max_build_size < CHAR_PER_WORD {
            eprintln!("bwt_inc_set_build_size_and_text_addr(): Not enough space allocated to continue construction!");
            std::process::exit(1);
        }
        if bwt_inc.inc_max_build_size > 0 {
            bwt_inc.build_size = bwt_inc.inc_max_build_size.min(max_build_size);
        } else {
            bwt_inc.build_size = max_build_size;
        }
        if bwt_inc.build_size < CHAR_PER_WORD {
            bwt_inc.build_size = CHAR_PER_WORD;
        }
    }
    if bwt_inc.build_size < CHAR_PER_WORD {
        eprintln!(
            "bwt_inc_set_build_size_and_text_addr(): Not enough space allocated to continue construction!"
        );
        std::process::exit(1);
    }
    bwt_inc.build_size = bwt_inc.build_size / CHAR_PER_WORD * CHAR_PER_WORD;
    bwt_inc.packed_text_offset = (2 * (bwt_inc.build_size + 1) * word_scale) as usize;
    bwt_inc.text_buffer_offset = ((bwt_inc.build_size + 1) * word_scale) as usize;
}

/// Original C global function `leading_zero` from `minibwa/bwtgen.c:203`.
pub(crate) fn leading_zero(input: u32) -> u32 {
    const LEADING_ZERO_8BIT: [u32; 256] = [
        8, 7, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    if (input & 0xffff0000) != 0 {
        if (input & 0xff000000) != 0 {
            LEADING_ZERO_8BIT[(input >> 24) as usize]
        } else {
            8 + LEADING_ZERO_8BIT[(input >> 16) as usize]
        }
    } else if (input & 0x0000ff00) != 0 {
        16 + LEADING_ZERO_8BIT[(input >> 8) as usize]
    } else {
        24 + LEADING_ZERO_8BIT[input as usize]
    }
}

/// Original C static function `ceil_log2` from `minibwa/bwtgen.c:232`.
pub(crate) fn ceil_log2(input: u32) -> u32 {
    if input <= 1 {
        0
    } else {
        BITS_IN_WORD - leading_zero(input - 1)
    }
}

/// Original C static function `bit_per_byte_packed_char` from `minibwa/bwtgen.c:239`.
pub(crate) fn bit_per_byte_packed_char(alphabet_size: u32) -> u32 {
    let mut bit_per_char = ceil_log2(alphabet_size);
    if BITS_IN_BYTE / (BITS_IN_BYTE / bit_per_char) > bit_per_char {
        bit_per_char = BITS_IN_BYTE / (BITS_IN_BYTE / bit_per_char);
    }
    bit_per_char
}

/// Original C static function `bit_per_word_packed_char` from `minibwa/bwtgen.c:249`.
pub(crate) fn bit_per_word_packed_char(alphabet_size: u32) -> u32 {
    ceil_log2(alphabet_size)
}

/// Original C static function `convert_byte_packed_to_word_packed` from `minibwa/bwtgen.c:254`.
pub(crate) fn convert_byte_packed_to_word_packed(
    input: &[u8],
    output: &mut [u32],
    alphabet_size: u32,
    text_length: BgInt,
) {
    let bit_per_byte_packed_char = bit_per_byte_packed_char(alphabet_size);
    let bit_per_word_packed_char = bit_per_word_packed_char(alphabet_size);
    let char_per_byte = BITS_IN_BYTE / bit_per_byte_packed_char;
    let char_per_word = BITS_IN_WORD / bit_per_word_packed_char;
    let byte_per_iteration = char_per_word / char_per_byte;
    let mask = ALL_ONE_MASK >> (BITS_IN_WORD - bit_per_word_packed_char)
        << (BITS_IN_WORD - bit_per_word_packed_char);
    let shift = BITS_IN_WORD - BITS_IN_BYTE + bit_per_byte_packed_char - bit_per_word_packed_char;
    let mut byte_processed = 0usize;
    let mut word_processed = 0usize;
    let mut buffer = [0u32; BITS_IN_WORD as usize];
    while ((word_processed + 1) as u64) * (char_per_word as u64) < text_length {
        let mut k = 0usize;
        for _ in 0..byte_per_iteration {
            let mut c = (input[byte_processed] as u32) << shift;
            for _ in 0..char_per_byte {
                buffer[k] = c & mask;
                c <<= bit_per_byte_packed_char;
                k += 1;
            }
            byte_processed += 1;
        }
        let mut c = 0u32;
        for (i, &b) in buffer.iter().enumerate().take(char_per_word as usize) {
            c |= b >> (bit_per_word_packed_char * i as u32);
        }
        output[word_processed] = c;
        word_processed += 1;
    }
    let mut k = 0usize;
    let remaining = text_length - word_processed as u64 * char_per_word as u64;
    for _ in 0..((remaining - 1) / char_per_byte as u64 + 1) {
        let mut c = (input[byte_processed] as u32) << shift;
        for _ in 0..char_per_byte {
            buffer[k] = c & mask;
            c <<= bit_per_byte_packed_char;
            k += 1;
        }
        byte_processed += 1;
    }
    let mut c = 0u32;
    for (i, &b) in buffer.iter().enumerate().take(remaining as usize) {
        c |= b >> (bit_per_word_packed_char * i as u32);
    }
    output[word_processed] = c;
}

/// Original C global function `bwt_create` from `minibwa/bwtgen.c:319`.
pub(crate) fn bwt_create(text_length: BgInt, decode_table: Option<Vec<u32>>) -> GenBwt {
    let mut cumulative_freq = vec![0; ALPHABET_SIZE + 1];
    initialize_val_bg(&mut cumulative_freq, (ALPHABET_SIZE + 1) as u64, 0);
    let decode_table = match decode_table {
        Some(t) => t,
        None => {
            let mut t = vec![0; DNA_OCC_CNT_TABLE_SIZE_IN_WORD];
            generate_dna_occ_count_table(&mut t);
            t
        }
    };
    let occ_major_size_in_word = bwt_occ_value_major_size_in_word(text_length);
    GenBwt {
        text_length: 0,
        inverse_sa0: 0,
        cumulative_freq,
        bwt_code: Vec::new(),
        occ_value: Vec::new(),
        occ_value_major: vec![0; occ_major_size_in_word as usize],
        decode_table,
        bwt_size_in_word: 0,
        occ_size_in_word: 0,
        occ_major_size_in_word,
    }
}

/// Original C global function `bwt_inc_create` from `minibwa/bwtgen.c:350`.
pub(crate) fn bwt_inc_create(
    text_length: BgInt,
    mut initial_max_build_size: u32,
    mut inc_max_build_size: u32,
) -> BwtInc {
    if text_length < inc_max_build_size as u64 {
        inc_max_build_size = text_length as u32;
    }
    if text_length < initial_max_build_size as u64 {
        initial_max_build_size = text_length as u32;
    }
    let mut packed_shift = vec![0; CHAR_PER_WORD as usize];
    for (i, x) in packed_shift.iter_mut().enumerate() {
        *x = BITS_IN_WORD - (i as u32 + 1) * BIT_PER_CHAR;
    }
    let n_iter = (text_length - initial_max_build_size as u64) / inc_max_build_size as u64 + 1;
    let word_scale = (std::mem::size_of::<BgInt>() / 4) as u64;
    let mut available_word = bwt_resident_size_in_word(text_length)
        + bwt_occ_value_minor_size_in_word(text_length)
        + OCC_INTERVAL / BIT_PER_CHAR as u64 * n_iter * 2 * word_scale
        + inc_max_build_size as u64 / 5 * 3 * word_scale;
    if available_word < MIN_AVAILABLE_WORD {
        available_word = MIN_AVAILABLE_WORD;
    }
    eprintln!("[bwt_inc_create] text_length={text_length}, available_word={available_word}");
    BwtInc {
        bwt: bwt_create(text_length, None),
        number_of_iteration_done: 0,
        cumulative_count_in_current_build: vec![0; ALPHABET_SIZE + 1],
        available_word,
        build_size: 0,
        initial_max_build_size: initial_max_build_size as u64,
        inc_max_build_size: inc_max_build_size as u64,
        first_char_in_last_iteration: 0,
        working_memory: Vec::new(),
        packed_text: Vec::new(),
        text_buffer: Vec::new(),
        packed_shift,
        packed_text_offset: 0,
        text_buffer_offset: 0,
    }
}

/// Original C static function `bwt_inc_put_packed_text_to_rank` from `minibwa/bwtgen.c:382`.
pub(crate) fn bwt_inc_put_packed_text_to_rank(
    packed_text: &[u32],
    rank: &mut [BgInt],
    cumulative_count: &mut [BgInt],
    num_char: BgInt,
) {
    let last_word = (num_char - 1) / CHAR_PER_WORD;
    let num_char_in_last_word = num_char - last_word * CHAR_PER_WORD;
    let packed_mask = ALL_ONE_MASK >> (BITS_IN_WORD - BIT_PER_CHAR);
    let mut rank_index = num_char - 1;

    let mut t =
        packed_text[last_word as usize] >> (BITS_IN_WORD - num_char_in_last_word as u32 * BIT_PER_CHAR);
    for _ in 0..num_char_in_last_word {
        let c = t & packed_mask;
        cumulative_count[c as usize + 1] += 1;
        rank[rank_index as usize] = c as BgInt;
        rank_index = rank_index.wrapping_sub(1);
        t >>= BIT_PER_CHAR;
    }

    for i in (0..last_word as usize).rev() {
        t = packed_text[i];
        for _ in 0..CHAR_PER_WORD {
            let c = t & packed_mask;
            cumulative_count[c as usize + 1] += 1;
            rank[rank_index as usize] = c as BgInt;
            rank_index = rank_index.wrapping_sub(1);
            t >>= BIT_PER_CHAR;
        }
    }

    cumulative_count[2] += cumulative_count[1];
    cumulative_count[3] += cumulative_count[2];
    cumulative_count[4] += cumulative_count[3];
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_table_and_size_helpers_match_expected_values() {
        assert_eq!(text_length_from_byte_packed(3, 2, 2), 10);
        assert_eq!(bwt_occ_value_major_size_in_word(0), 4);
        assert_eq!(bwt_occ_value_major_size_in_word(65_536), 8);
        assert_eq!(bwt_occ_value_minor_size_in_word(0), 4);
        assert_eq!(bwt_occ_value_minor_size_in_word(256), 4);
        assert_eq!(bwt_occ_value_minor_size_in_word(257), 8);
        assert_eq!(bwt_resident_size_in_word(1), 16);
        assert_eq!(bwt_resident_size_in_word(256), 16);
        assert_eq!(bwt_resident_size_in_word(257), 32);

        let mut table = vec![0u32; DNA_OCC_CNT_TABLE_SIZE_IN_WORD];
        generate_dna_occ_count_table(&mut table);
        assert_eq!(table[0], 8);
        assert_eq!(table[0xffff], 8 << 24);
        assert_eq!(table[0x1b1b], 0x02020202);
    }

    #[test]
    fn bit_width_helpers_follow_original_rounding() {
        assert_eq!(leading_zero(0), 32);
        assert_eq!(leading_zero(1), 31);
        assert_eq!(leading_zero(0x8000_0000), 0);
        assert_eq!(ceil_log2(1), 0);
        assert_eq!(ceil_log2(2), 1);
        assert_eq!(ceil_log2(3), 2);
        assert_eq!(bit_per_byte_packed_char(4), 2);
        assert_eq!(bit_per_word_packed_char(4), 2);
    }

    #[test]
    fn byte_packed_to_word_packed_preserves_dna_order() {
        let input = [0x1b, 0x1b, 0x1b, 0x1b, 0x10];
        let mut output = [0u32; 2];
        convert_byte_packed_to_word_packed(&input, &mut output, 4, 18);
        assert_eq!(output[0], 0x1b1b1b1b);
        assert_eq!(output[1], 0x10000000);
    }

    #[test]
    fn bwt_create_and_inc_create_initialize_metadata() {
        let bwt = bwt_create(1000, None);
        assert_eq!(bwt.text_length, 0);
        assert_eq!(bwt.cumulative_freq, vec![0; 5]);
        assert_eq!(bwt.occ_major_size_in_word, bwt_occ_value_major_size_in_word(1000));
        assert_eq!(bwt.occ_value_major.len(), bwt.occ_major_size_in_word as usize);
        assert_eq!(bwt.decode_table[0x1b1b], 0x02020202);

        let mut inc = bwt_inc_create(1000, 320, 160);
        assert_eq!(inc.packed_shift[0], 30);
        assert_eq!(inc.packed_shift[15], 0);
        assert!(inc.available_word >= MIN_AVAILABLE_WORD);
        bwt_inc_set_build_size_and_text_addr(&mut inc);
        assert_eq!(inc.build_size % CHAR_PER_WORD, 0);
        assert!(inc.packed_text_offset > inc.text_buffer_offset);
    }

    #[test]
    fn packed_text_to_rank_decodes_reverse_then_cumulative_counts() {
        let packed = [0x1b1b1b1b, 0x10000000];
        let mut rank = [0u64; 18];
        let mut cumulative = [0u64; 5];
        bwt_inc_put_packed_text_to_rank(&packed, &mut rank, &mut cumulative, 18);
        assert_eq!(rank, [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1]);
        assert_eq!(cumulative, [0, 5, 10, 14, 18]);
    }

    #[test]
    fn incremental_sort_and_rank_builders_follow_key_groups() {
        let mut key = [5, 1, 3, 1, 4, 3, 2, 5, 0, 4];
        let mut seq = [50, 10, 30, 11, 40, 31, 20, 51, 0, 41];
        bwt_inc_sort_key(&mut key, &mut seq, 10);
        assert_eq!(key, [0, 1, 1, 2, 3, 3, 4, 4, 5, 5]);
        for pair in key.windows(2) {
            assert!(pair[0] <= pair[1]);
        }

        let mut sorted_rank = [1, 2, 2, 5, 8, 9];
        let mut grouped_seq = [0, 1, 2, 3, 4, 5];
        let mut relative = [0u64; 6];
        bwt_inc_build_relative_rank(
            &mut sorted_rank,
            &mut grouped_seq,
            &mut relative,
            5,
            4,
            &[0, 2, 4, 5, 6],
        );
        assert_eq!(relative[5], 5);
        assert_eq!(relative[1], 1);
        assert_eq!(relative[2], 2);
        assert_eq!(grouped_seq[5], (-1i64) as u64);
    }

    #[test]
    fn trailing_bwt_code_is_cleared_after_text_length() {
        let mut bwt = bwt_create(300, None);
        bwt.text_length = 18;
        bwt.bwt_code = vec![0xffff_ffff; bwt_resident_size_in_word(18) as usize];
        bwt_clear_trailing_bwt_code(&mut bwt);
        assert_eq!(bwt.bwt_code[0], 0xffff_ffff);
        assert_eq!(bwt.bwt_code[1], 0xf000_0000);
        assert!(bwt.bwt_code[2..].iter().all(|&x| x == 0));
    }

    #[test]
    fn bwt_file_size_and_save_emit_original_binary_layout() {
        assert_eq!(bwt_file_size_in_word(0), 0);
        assert_eq!(bwt_file_size_in_word(1), 1);
        assert_eq!(bwt_file_size_in_word(16), 1);
        assert_eq!(bwt_file_size_in_word(17), 2);

        let mut bwt = bwt_create(32, None);
        bwt.text_length = 17;
        bwt.inverse_sa0 = 9;
        bwt.cumulative_freq = vec![0, 3, 7, 11, 17];
        bwt.bwt_code = vec![0x0123_4567, 0x89ab_cdef, 0xffff_ffff];
        let path = std::env::temp_dir().join("minibwa-rs-test.bwt");
        bwt_save_bwt_code_and_occ(&bwt, &path, None).expect("save GenBwt");
        let bytes = std::fs::read(&path).expect("read GenBwt");
        let _ = std::fs::remove_file(&path);

        let mut expected = Vec::new();
        expected.extend_from_slice(&9u64.to_le_bytes());
        for value in [3u64, 7, 11, 17] {
            expected.extend_from_slice(&value.to_le_bytes());
        }
        expected.extend_from_slice(&0x0123_4567u32.to_le_bytes());
        expected.extend_from_slice(&0x89ab_cdefu32.to_le_bytes());
        assert_eq!(bytes, expected);
    }

    #[test]
    fn bwt_inc_free_accepts_null_and_owned_state() {
        bwt_inc_free(None);
        bwt_inc_free(Some(bwt_inc_create(64, 16, 16)));
    }

    #[test]
    fn public_file_helpers_stop_paths_at_embedded_nul_like_c() {
        let bases: Vec<u8> = (0..32).map(|i| ((i * 5 + 1) & 3) as u8).collect();
        let mut packed =
            vec![0u8; (bases.len() + CHAR_PER_BYTE as usize - 1) / CHAR_PER_BYTE as usize];
        for (i, &base) in bases.iter().enumerate() {
            packed[i / CHAR_PER_BYTE as usize] |=
                base << (BITS_IN_BYTE - ((i % CHAR_PER_BYTE as usize) as u32 + 1) * BIT_PER_CHAR);
        }
        let last_len = (bases.len() % CHAR_PER_BYTE as usize) as u8;
        if last_len == 0 {
            packed.push(0);
        }
        packed.push(last_len);

        let dir =
            std::env::temp_dir().join(format!("minibwa-rs-bwtgen-api-nul-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("create temp dir");
        let pac = dir.join("in.pac");
        let raw = dir.join("out.raw");
        std::fs::write(&pac, &packed).expect("write pac");
        let pac_arg = std::path::PathBuf::from(format!("{}\0hidden", pac.to_string_lossy()));
        let raw_arg = std::path::PathBuf::from(format!("{}\0hidden", raw.to_string_lossy()));

        let bwt_inc = bwt_inc_construct_from_packed(&pac_arg, 64, 64).expect("construct GenBwt");
        assert_eq!(bwt_inc.bwt.text_length, bases.len() as u64);
        bwt_save_bwt_code_and_occ(&bwt_inc.bwt, &raw_arg, None).expect("save GenBwt");
        assert!(raw.exists());

        let raw2 = dir.join("out2.raw");
        let raw2_arg = std::path::PathBuf::from(format!("{}\0hidden", raw2.to_string_lossy()));
        mb_bwtgen(&pac_arg, &raw2_arg, 64).expect("generate GenBwt");
        assert!(raw2.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

/// Original C static function `forward_dna_all_occ_count_no_limit` from `minibwa/bwtgen.c:426`.
pub(crate) fn forward_dna_all_occ_count_no_limit(
    dna: &[u32],
    index: BgInt,
    occ_count: &mut [BgInt],
    dna_decode_table: &[u32],
) {
    const TRUNCATE_RIGHT_MASK: [u32; 16] = [
        0x00000000, 0xC0000000, 0xF0000000, 0xFC000000, 0xFF000000, 0xFFC00000, 0xFFF00000,
        0xFFFC0000, 0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00, 0xFFFFFF00, 0xFFFFFFC0,
        0xFFFFFFF0, 0xFFFFFFFC,
    ];

    occ_count[0] = 0;
    occ_count[1] = 0;
    occ_count[2] = 0;
    occ_count[3] = 0;

    let iteration = index / 256;
    let word_to_count = (index - iteration * 256) / 16;
    let char_to_count = index - iteration * 256 - word_to_count * 16;
    let mut dna_i = 0usize;

    for _ in 0..iteration {
        let mut sum = 0u32;
        for _ in 0..16 {
            let word = dna[dna_i];
            sum = sum.wrapping_add(dna_decode_table[(word >> 16) as usize]);
            sum = sum.wrapping_add(dna_decode_table[(word & 0x0000ffff) as usize]);
            dna_i += 1;
        }
        if (sum & 0xfefefeff) != 0 {
            occ_count[0] += (sum & 0x000000ff) as BgInt;
            sum >>= 8;
            occ_count[1] += (sum & 0x000000ff) as BgInt;
            sum >>= 8;
            occ_count[2] += (sum & 0x000000ff) as BgInt;
            sum >>= 8;
            occ_count[3] += sum as BgInt;
        } else if sum == 0x00000100 {
            occ_count[0] += 256;
        } else if sum == 0x00010000 {
            occ_count[1] += 256;
        } else if sum == 0x01000000 {
            occ_count[2] += 256;
        } else if sum == 0x00000000 {
            occ_count[3] += 256;
        } else {
            eprintln!("forward_dna_all_occ_count_no_limit(): DNA occ sum exception!");
            std::process::exit(1);
        }
    }

    let mut sum = 0u32;
    for _ in 0..word_to_count {
        let word = dna[dna_i];
        sum = sum.wrapping_add(dna_decode_table[(word >> 16) as usize]);
        sum = sum.wrapping_add(dna_decode_table[(word & 0x0000ffff) as usize]);
        dna_i += 1;
    }

    if char_to_count > 0 {
        let c = dna[dna_i] & TRUNCATE_RIGHT_MASK[char_to_count as usize];
        sum = sum.wrapping_add(dna_decode_table[(c >> 16) as usize]);
        sum = sum.wrapping_add(dna_decode_table[(c & 0xffff) as usize]);
        sum = sum.wrapping_add(char_to_count as u32).wrapping_sub(16);
    }

    occ_count[0] += (sum & 0x000000ff) as BgInt;
    sum >>= 8;
    occ_count[1] += (sum & 0x000000ff) as BgInt;
    sum >>= 8;
    occ_count[2] += (sum & 0x000000ff) as BgInt;
    sum >>= 8;
    occ_count[3] += sum as BgInt;
}

/// Original C static function `bwt_inc_build_packed_bwt` from `minibwa/bwtgen.c:499`.
pub(crate) fn bwt_inc_build_packed_bwt(
    relative_rank: &[BgInt],
    bwt: &mut [u32],
    num_char: BgInt,
    cumulative_count: &[BgInt],
    packed_shift: &[u32],
) {
    let inverse_sa0 = relative_rank[0];
    let mut previous_rank = relative_rank[0];
    for i in 1..=num_char as usize {
        let current_rank = relative_rank[i];
        let c = ((previous_rank > cumulative_count[1]) as u32)
            + ((previous_rank > cumulative_count[2]) as u32)
            + ((previous_rank > cumulative_count[3]) as u32);
        if c > 0 {
            let mut r = current_rank;
            if r > inverse_sa0 {
                r -= 1;
            }
            let word_index = r / CHAR_PER_WORD;
            let char_index = r - word_index * CHAR_PER_WORD;
            bwt[word_index as usize] |= c << packed_shift[char_index as usize];
        }
        previous_rank = current_rank;
    }
}

#[inline]
fn bwt_occ_value_explicit(bwt: &GenBwt, occ_index_explicit: BgInt, character: u32) -> BgInt {
    let occ_index_major = occ_index_explicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR;
    let major_idx = occ_index_major * ALPHABET_SIZE as u64 + character as u64;
    let minor_idx = occ_index_explicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE as u64 + character as u64;
    let packed = bwt.occ_value[minor_idx as usize];
    bwt.occ_value_major[major_idx as usize]
        + if occ_index_explicit % OCC_VALUE_PER_WORD == 0 {
            (packed >> 16) as u64
        } else {
            (packed & 0xffff) as u64
        }
}

#[inline]
fn forward_dna_occ_count(
    dna: &[u32],
    index: u32,
    character: u32,
    dna_decode_table: &[u32],
) -> u32 {
    const TRUNCATE_RIGHT_MASK: [u32; 16] = [
        0x00000000, 0xC0000000, 0xF0000000, 0xFC000000, 0xFF000000, 0xFFC00000, 0xFFF00000,
        0xFFFC0000, 0xFFFF0000, 0xFFFFC000, 0xFFFFF000, 0xFFFFFC00, 0xFFFFFF00, 0xFFFFFFC0,
        0xFFFFFFF0, 0xFFFFFFFC,
    ];

    let word_to_count = index / 16;
    let char_to_count = index - word_to_count * 16;
    let mut sum = 0u32;

    for i in 0..word_to_count as usize {
        let word = dna[i];
        sum = sum.wrapping_add(dna_decode_table[(word >> 16) as usize]);
        sum = sum.wrapping_add(dna_decode_table[(word & 0x0000ffff) as usize]);
    }

    if char_to_count > 0 {
        let c = dna[word_to_count as usize] & TRUNCATE_RIGHT_MASK[char_to_count as usize];
        sum = sum.wrapping_add(dna_decode_table[(c >> 16) as usize]);
        sum = sum.wrapping_add(dna_decode_table[(c & 0xffff) as usize]);
        sum = sum.wrapping_add(char_to_count).wrapping_sub(16);
    }

    (sum >> (character * 8)) & 0x000000ff
}

#[inline]
fn backward_dna_occ_count(
    dna: &[u32],
    word_index: usize,
    index: u32,
    character: u32,
    dna_decode_table: &[u32],
) -> u32 {
    const TRUNCATE_LEFT_MASK: [u32; 16] = [
        0x00000000, 0x00000003, 0x0000000f, 0x0000003f, 0x000000ff, 0x000003ff, 0x00000fff,
        0x00003fff, 0x0000ffff, 0x0003ffff, 0x000fffff, 0x003fffff, 0x00ffffff, 0x03ffffff,
        0x0fffffff, 0x3fffffff,
    ];

    let word_to_count = index / 16;
    let char_to_count = index - word_to_count * 16;
    let mut sum = 0u32;
    let mut pos = word_index - word_to_count as usize - 1;

    if char_to_count > 0 {
        let c = dna[pos] & TRUNCATE_LEFT_MASK[char_to_count as usize];
        sum = sum.wrapping_add(dna_decode_table[(c >> 16) as usize]);
        sum = sum.wrapping_add(dna_decode_table[(c & 0xffff) as usize]);
        sum = sum.wrapping_add(char_to_count).wrapping_sub(16);
    }

    for _ in 0..word_to_count {
        pos += 1;
        let word = dna[pos];
        sum = sum.wrapping_add(dna_decode_table[(word >> 16) as usize]);
        sum = sum.wrapping_add(dna_decode_table[(word & 0x0000ffff) as usize]);
    }

    (sum >> (character * 8)) & 0x000000ff
}

#[inline]
fn bwt_occ_value(bwt: &GenBwt, mut index: BgInt, character: u32) -> BgInt {
    if index > bwt.inverse_sa0 {
        index -= 1;
    }

    let occ_explicit_index = (index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
    let occ_index = occ_explicit_index * OCC_INTERVAL;
    let occ_value = bwt_occ_value_explicit(bwt, occ_explicit_index, character);

    if occ_index == index {
        return occ_value;
    }

    let word_index = (occ_index / CHAR_PER_WORD) as usize;
    if occ_index < index {
        occ_value
            + forward_dna_occ_count(
                &bwt.bwt_code[word_index..],
                (index - occ_index) as u32,
                character,
                &bwt.decode_table,
            ) as u64
    } else {
        occ_value
            - backward_dna_occ_count(
                &bwt.bwt_code,
                word_index,
                (occ_index - index) as u32,
                character,
                &bwt.decode_table,
            ) as u64
    }
}

/// Original C static function `bwt_inc_get_absolute_rank` from `minibwa/bwtgen.c:638`.
pub(crate) fn bwt_inc_get_absolute_rank(
    bwt: &GenBwt,
    absolute_rank: &mut [BgInt],
    seq: &mut [BgInt],
    packed_text: &[u32],
    num_char: BgInt,
    cumulative_count: &[BgInt],
    first_char_in_last_iteration: u32,
) -> BgInt {
    let mut seq_index_from_start = [0u64; ALPHABET_SIZE];
    let mut seq_index_from_end = [0u64; ALPHABET_SIZE];
    for i in 0..ALPHABET_SIZE {
        seq_index_from_start[i] = cumulative_count[i];
        seq_index_from_end[i] = cumulative_count[i + 1] - 1;
    }

    let shift = BITS_IN_WORD - BIT_PER_CHAR;
    let packed_mask = ALL_ONE_MASK >> shift;
    let mut sa_index = bwt.inverse_sa0;
    let mut rank_index = num_char - 1;
    let last_word = num_char / CHAR_PER_WORD;

    for i in (0..last_word as usize).rev() {
        let mut t = packed_text[i];
        for _ in 0..CHAR_PER_WORD {
            let c = t & packed_mask;
            sa_index = bwt.cumulative_freq[c as usize] + bwt_occ_value(bwt, sa_index, c) + 1;
            if sa_index > bwt.inverse_sa0 {
                let idx = seq_index_from_end[c as usize] as usize;
                seq[idx] = rank_index;
                absolute_rank[idx] = sa_index;
                seq_index_from_end[c as usize] -= 1;
            } else {
                let idx = seq_index_from_start[c as usize] as usize;
                seq[idx] = rank_index;
                absolute_rank[idx] = sa_index;
                seq_index_from_start[c as usize] += 1;
            }
            rank_index = rank_index.wrapping_sub(1);
            t >>= BIT_PER_CHAR;
        }
    }

    let idx = seq_index_from_start[first_char_in_last_iteration as usize] as usize;
    absolute_rank[idx] = bwt.inverse_sa0;
    seq[idx] = num_char;
    seq_index_from_start[first_char_in_last_iteration as usize]
}

/// Original C static function `bwt_inc_sort_key` from `minibwa/bwtgen.c:690`.
pub(crate) fn bwt_inc_sort_key(key: &mut [BgInt], seq: &mut [BgInt], num_item: BgInt) {
    const INSERT_SORT_NUM_ITEM: usize = BWTINC_INSERT_SORT_NUM_ITEM as usize;
    const EQUAL_KEY_THRESHOLD: usize = 4;
    if num_item < 2 {
        return;
    }

    let mut low_index = 0usize;
    let mut high_index = num_item as usize - 1;
    let mut low_stack = [0usize; 32];
    let mut high_stack = [0usize; 32];
    let mut stack_depth = 0usize;

    loop {
        loop {
            if high_index - low_index < INSERT_SORT_NUM_ITEM {
                for i in low_index + 1..=high_index {
                    let temp_seq = seq[i];
                    let temp_key = key[i];
                    let mut j = i;
                    while j > low_index && key[j - 1] > temp_key {
                        seq[j] = seq[j - 1];
                        key[j] = key[j - 1];
                        j -= 1;
                    }
                    if j != i {
                        seq[j] = temp_seq;
                        key[j] = temp_key;
                    }
                }
                break;
            }

            let mut mid_index = (low_index & high_index) + ((low_index ^ high_index) / 2);
            if key[low_index] > key[mid_index] {
                swap_key_seq(key, seq, low_index, mid_index);
            }
            if key[low_index] > key[high_index] {
                swap_key_seq(key, seq, low_index, high_index);
            }
            if key[mid_index] > key[high_index] {
                swap_key_seq(key, seq, mid_index, high_index);
            }

            let mut number_of_equal_key = 0usize;
            let mut low_partition_index = low_index + 1;
            let mut high_partition_index = high_index - 1;

            loop {
                while low_partition_index <= high_partition_index
                    && key[low_partition_index] <= key[mid_index]
                {
                    number_of_equal_key += (key[low_partition_index] == key[mid_index]) as usize;
                    low_partition_index += 1;
                }
                while low_partition_index < high_partition_index {
                    if key[mid_index] >= key[high_partition_index] {
                        number_of_equal_key +=
                            (key[mid_index] == key[high_partition_index]) as usize;
                        break;
                    }
                    high_partition_index -= 1;
                }
                if low_partition_index >= high_partition_index {
                    break;
                }
                swap_key_seq(key, seq, low_partition_index, high_partition_index);
                if high_partition_index == mid_index {
                    mid_index = low_partition_index;
                }
                low_partition_index += 1;
                high_partition_index -= 1;
            }

            high_partition_index = low_partition_index;
            low_partition_index -= 1;

            swap_key_seq(key, seq, mid_index, low_partition_index);

            if high_index - low_index + INSERT_SORT_NUM_ITEM
                <= EQUAL_KEY_THRESHOLD * number_of_equal_key
            {
                mid_index = low_index;
                loop {
                    while mid_index < low_partition_index
                        && key[mid_index] < key[low_partition_index]
                    {
                        mid_index += 1;
                    }
                    while mid_index < low_partition_index
                        && key[low_partition_index] == key[low_partition_index - 1]
                    {
                        low_partition_index -= 1;
                    }
                    if mid_index >= low_partition_index {
                        break;
                    }
                    swap_key_seq(key, seq, mid_index, low_partition_index - 1);
                    mid_index += 1;
                    low_partition_index -= 1;
                }
            }

            if low_partition_index - low_index > high_index - high_partition_index {
                low_stack[stack_depth] = low_index;
                high_stack[stack_depth] = low_partition_index - 1;
                stack_depth += 1;
                low_index = high_partition_index;
            } else {
                low_stack[stack_depth] = high_partition_index;
                high_stack[stack_depth] = high_index;
                stack_depth += 1;
                if low_partition_index > low_index {
                    high_index = low_partition_index - 1;
                } else {
                    break;
                }
            }
        }

        if stack_depth > 0 {
            stack_depth -= 1;
            low_index = low_stack[stack_depth];
            high_index = high_stack[stack_depth];
        } else {
            return;
        }
    }
}

#[inline]
fn swap_key_seq(key: &mut [BgInt], seq: &mut [BgInt], a: usize, b: usize) {
    seq.swap(a, b);
    key.swap(a, b);
}

/// Original C static function `bwt_inc_build_relative_rank` from `minibwa/bwtgen.c:868`.
pub(crate) fn bwt_inc_build_relative_rank(
    sorted_rank: &mut [BgInt],
    seq: &mut [BgInt],
    relative_rank: &mut [BgInt],
    num_item: BgInt,
    mut old_inverse_sa0: BgInt,
    cumulative_count: &[BgInt],
) {
    let mut last_index = num_item;
    let mut last_rank = sorted_rank[num_item as usize];
    if last_rank > old_inverse_sa0 {
        sorted_rank[num_item as usize] -= 1;
    }
    let mut s = seq[num_item as usize];
    relative_rank[s as usize] = num_item;
    if last_rank == old_inverse_sa0 {
        old_inverse_sa0 += 1;
        last_rank += 1;
    }

    let mut c = ALPHABET_SIZE as u64 - 1;
    let mut freq = cumulative_count[c as usize];
    for i in (0..num_item).rev() {
        let r = sorted_rank[i as usize];
        if r > old_inverse_sa0 {
            sorted_rank[i as usize] -= 1;
        }
        s = seq[i as usize];
        if i < freq {
            if last_index >= freq {
                last_rank += 1;
            }
            c -= 1;
            freq = cumulative_count[c as usize];
        }
        if r == last_rank {
            relative_rank[s as usize] = last_index;
        } else {
            if i == last_index - 1 {
                if last_index < num_item && (seq[(last_index + 1) as usize] as SBgInt) < 0 {
                    seq[last_index as usize] = seq[(last_index + 1) as usize].wrapping_sub(1);
                } else {
                    seq[last_index as usize] = (-1i64) as BgInt;
                }
            }
            last_index = i;
            last_rank = r;
            relative_rank[s as usize] = i;
            if r == old_inverse_sa0 {
                old_inverse_sa0 += 1;
                last_rank += 1;
            }
        }
    }
}

/// Original C static function `bwt_inc_build_bwt` from `minibwa/bwtgen.c:925`.
pub(crate) fn bwt_inc_build_bwt(
    insert_bwt: &mut [u32],
    relative_rank: &[BgInt],
    num_char: BgInt,
    cumulative_count: &[BgInt],
) {
    let mut previous_rank = relative_rank[0];
    for i in 1..=num_char as usize {
        let current_rank = relative_rank[i];
        let c = ((previous_rank >= cumulative_count[1]) as u32)
            + ((previous_rank >= cumulative_count[2]) as u32)
            + ((previous_rank >= cumulative_count[3]) as u32);
        insert_bwt[current_rank as usize] = c;
        previous_rank = current_rank;
    }
}

/// Original C static function `bwt_inc_merge_bwt` from `minibwa/bwtgen.c:943`.
pub(crate) fn bwt_inc_merge_bwt(
    sorted_rank: &[BgInt],
    old_bwt: &[u32],
    insert_bwt: &[u32],
    merged_bwt: &mut [u32],
    num_old_bwt: BgInt,
    num_insert_bwt: BgInt,
) {
    let mut o_index = 0u64;
    let mut i_index = 0u64;
    let mut m_index = 0u64;
    let mut m_word = 0u64;
    let mut m_char = 0u64;
    merged_bwt[0] = 0;

    while o_index < num_old_bwt {
        while i_index <= num_insert_bwt && sorted_rank[i_index as usize] <= o_index {
            let rank = sorted_rank[i_index as usize];
            if rank != 0 {
                merged_bwt[m_word as usize] |= insert_bwt[i_index as usize]
                    << (BITS_IN_WORD - (m_char as u32 + 1) * BIT_PER_CHAR);
                m_index += 1;
                m_char += 1;
                if m_char == CHAR_PER_WORD {
                    m_char = 0;
                    m_word += 1;
                    merged_bwt[m_word as usize] = 0;
                }
            }
            i_index += 1;
        }

        let o = if i_index <= num_insert_bwt {
            sorted_rank[i_index as usize]
        } else {
            num_old_bwt
        };
        let num_insert = o - o_index;
        let mut o_word = o_index / CHAR_PER_WORD;
        let o_char = o_index - o_word * CHAR_PER_WORD;

        if o_char > m_char {
            let left_shift = ((o_char - m_char) as u32) * BIT_PER_CHAR;
            let right_shift = ((CHAR_PER_WORD + m_char - o_char) as u32) * BIT_PER_CHAR;
            merged_bwt[m_word as usize] |= (old_bwt[o_word as usize]
                << (o_char as u32 * BIT_PER_CHAR)
                >> (m_char as u32 * BIT_PER_CHAR))
                | (old_bwt[o_word as usize + 1] >> right_shift);
            o_index += num_insert.min(CHAR_PER_WORD - m_char);
            while o > o_index {
                o_word += 1;
                m_word += 1;
                merged_bwt[m_word as usize] = (old_bwt[o_word as usize] << left_shift)
                    | (old_bwt[o_word as usize + 1] >> right_shift);
                o_index += CHAR_PER_WORD;
            }
        } else if o_char < m_char {
            let right_shift = ((m_char - o_char) as u32) * BIT_PER_CHAR;
            let left_shift = ((CHAR_PER_WORD + o_char - m_char) as u32) * BIT_PER_CHAR;
            merged_bwt[m_word as usize] |= old_bwt[o_word as usize]
                << (o_char as u32 * BIT_PER_CHAR)
                >> (m_char as u32 * BIT_PER_CHAR);
            o_index += num_insert.min(CHAR_PER_WORD - m_char);
            while o > o_index {
                o_word += 1;
                m_word += 1;
                merged_bwt[m_word as usize] = (old_bwt[o_word as usize - 1] << left_shift)
                    | (old_bwt[o_word as usize] >> right_shift);
                o_index += CHAR_PER_WORD;
            }
        } else {
            merged_bwt[m_word as usize] |= old_bwt[o_word as usize]
                << (m_char as u32 * BIT_PER_CHAR)
                >> (m_char as u32 * BIT_PER_CHAR);
            o_index += num_insert.min(CHAR_PER_WORD - m_char);
            while o > o_index {
                o_word += 1;
                m_word += 1;
                merged_bwt[m_word as usize] = old_bwt[o_word as usize];
                o_index += CHAR_PER_WORD;
            }
        }

        o_index = o;
        m_index += num_insert;
        m_word = m_index / CHAR_PER_WORD;
        m_char = m_index - m_word * CHAR_PER_WORD;
        if m_char == 0 {
            merged_bwt[m_word as usize] = 0;
        } else {
            let offset = BITS_IN_WORD - m_char as u32 * BIT_PER_CHAR;
            merged_bwt[m_word as usize] = merged_bwt[m_word as usize] >> offset << offset;
        }
    }

    while i_index <= num_insert_bwt {
        if sorted_rank[i_index as usize] != 0 {
            merged_bwt[m_word as usize] |= insert_bwt[i_index as usize]
                << (BITS_IN_WORD - (m_char as u32 + 1) * BIT_PER_CHAR);
            m_char += 1;
            if m_char == CHAR_PER_WORD {
                m_char = 0;
                m_word += 1;
                merged_bwt[m_word as usize] = 0;
            }
        }
        i_index += 1;
    }
}

/// Original C global function `bwt_clear_trailing_bwt_code` from `minibwa/bwtgen.c:1053`.
pub(crate) fn bwt_clear_trailing_bwt_code(bwt: &mut GenBwt) {
    let bwt_resident_size_in_word = bwt_resident_size_in_word(bwt.text_length);
    let word_index = bwt.text_length / CHAR_PER_WORD;
    let offset = (bwt.text_length - word_index * CHAR_PER_WORD) as u32 * BIT_PER_CHAR;
    if offset > 0 {
        let rshift = BITS_IN_WORD - offset;
        bwt.bwt_code[word_index as usize] = bwt.bwt_code[word_index as usize] >> rshift << rshift;
    } else if word_index < bwt_resident_size_in_word {
        bwt.bwt_code[word_index as usize] = 0;
    }
    for i in word_index + 1..bwt_resident_size_in_word {
        bwt.bwt_code[i as usize] = 0;
    }
}

/// Original C global function `bwt_generate_occ_value_from_bwt` from `minibwa/bwtgen.c:1077`.
pub(crate) fn bwt_generate_occ_value_from_bwt(
    bwt: &[u32],
    occ_value: &mut [u32],
    occ_value_major: &mut [BgInt],
    text_length: BgInt,
    decode_table: &[u32],
) {
    let word_between_occ_value = OCC_INTERVAL / CHAR_PER_WORD;
    let number_of_occ_value = (text_length + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    let number_of_occ_interval_per_major = OCC_INTERVAL_MAJOR / OCC_INTERVAL;
    let number_of_occ_value_major =
        (number_of_occ_value + number_of_occ_interval_per_major - 1) / number_of_occ_interval_per_major;

    let mut temp_occ_value0 = [0u64; ALPHABET_SIZE];
    let mut temp_occ_value1 = [0u64; ALPHABET_SIZE];
    occ_value_major[0] = 0;
    occ_value_major[1] = 0;
    occ_value_major[2] = 0;
    occ_value_major[3] = 0;

    let mut occ_index = 0u64;
    let mut bwt_index = 0usize;
    for occ_major_index in 1..number_of_occ_value_major {
        for _ in 0..number_of_occ_interval_per_major / 2 {
            let mut sum = 0u64;
            temp_occ_value1.copy_from_slice(&temp_occ_value0);
            for _ in 0..word_between_occ_value {
                let c = bwt[bwt_index];
                sum += decode_table[(c >> 16) as usize] as u64;
                sum += decode_table[(c & 0x0000ffff) as usize] as u64;
                bwt_index += 1;
            }
            if (sum & 0xfefefeff) != 0 {
                temp_occ_value1[0] += sum & 0x000000ff;
                sum >>= 8;
                temp_occ_value1[1] += sum & 0x000000ff;
                sum >>= 8;
                temp_occ_value1[2] += sum & 0x000000ff;
                sum >>= 8;
                temp_occ_value1[3] += sum;
            } else if sum == 0x00000100 {
                temp_occ_value1[0] += 256;
            } else if sum == 0x00010000 {
                temp_occ_value1[1] += 256;
            } else if sum == 0x01000000 {
                temp_occ_value1[2] += 256;
            } else {
                temp_occ_value1[3] += 256;
            }
            for c in 0..ALPHABET_SIZE {
                occ_value[occ_index as usize * 4 + c] =
                    ((temp_occ_value0[c] << 16) | temp_occ_value1[c]) as u32;
                temp_occ_value0[c] = temp_occ_value1[c];
            }
            occ_index += 1;

            sum = 0;
            for _ in 0..word_between_occ_value {
                let c = bwt[bwt_index];
                sum += decode_table[(c >> 16) as usize] as u64;
                sum += decode_table[(c & 0x0000ffff) as usize] as u64;
                bwt_index += 1;
            }
            if (sum & 0xfefefeff) != 0 {
                temp_occ_value0[0] += sum & 0x000000ff;
                sum >>= 8;
                temp_occ_value0[1] += sum & 0x000000ff;
                sum >>= 8;
                temp_occ_value0[2] += sum & 0x000000ff;
                sum >>= 8;
                temp_occ_value0[3] += sum;
            } else if sum == 0x00000100 {
                temp_occ_value0[0] += 256;
            } else if sum == 0x00010000 {
                temp_occ_value0[1] += 256;
            } else if sum == 0x01000000 {
                temp_occ_value0[2] += 256;
            } else {
                temp_occ_value0[3] += 256;
            }
        }

        for c in 0..ALPHABET_SIZE {
            occ_value_major[occ_major_index as usize * 4 + c] =
                occ_value_major[(occ_major_index - 1) as usize * 4 + c] + temp_occ_value0[c];
            temp_occ_value0[c] = 0;
        }
    }

    while occ_index < (number_of_occ_value - 1) / 2 {
        let mut sum = 0u64;
        temp_occ_value1.copy_from_slice(&temp_occ_value0);
        for _ in 0..word_between_occ_value {
            let c = bwt[bwt_index];
            sum += decode_table[(c >> 16) as usize] as u64;
            sum += decode_table[(c & 0x0000ffff) as usize] as u64;
            bwt_index += 1;
        }
        if (sum & 0xfefefeff) != 0 {
            temp_occ_value1[0] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value1[1] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value1[2] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value1[3] += sum;
        } else if sum == 0x00000100 {
            temp_occ_value1[0] += 256;
        } else if sum == 0x00010000 {
            temp_occ_value1[1] += 256;
        } else if sum == 0x01000000 {
            temp_occ_value1[2] += 256;
        } else {
            temp_occ_value1[3] += 256;
        }
        for c in 0..ALPHABET_SIZE {
            occ_value[occ_index as usize * 4 + c] =
                ((temp_occ_value0[c] << 16) | temp_occ_value1[c]) as u32;
            temp_occ_value0[c] = temp_occ_value1[c];
        }
        occ_index += 1;

        sum = 0;
        for _ in 0..word_between_occ_value {
            let c = bwt[bwt_index];
            sum += decode_table[(c >> 16) as usize] as u64;
            sum += decode_table[(c & 0x0000ffff) as usize] as u64;
            bwt_index += 1;
        }
        if (sum & 0xfefefeff) != 0 {
            temp_occ_value0[0] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value0[1] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value0[2] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value0[3] += sum;
        } else if sum == 0x00000100 {
            temp_occ_value0[0] += 256;
        } else if sum == 0x00010000 {
            temp_occ_value0[1] += 256;
        } else if sum == 0x01000000 {
            temp_occ_value0[2] += 256;
        } else {
            temp_occ_value0[3] += 256;
        }
    }

    let mut sum = 0u64;
    temp_occ_value1.copy_from_slice(&temp_occ_value0);
    if occ_index * 2 < number_of_occ_value - 1 {
        for _ in 0..word_between_occ_value {
            let c = bwt[bwt_index];
            sum += decode_table[(c >> 16) as usize] as u64;
            sum += decode_table[(c & 0x0000ffff) as usize] as u64;
            bwt_index += 1;
        }
        if (sum & 0xfefefeff) != 0 {
            temp_occ_value1[0] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value1[1] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value1[2] += sum & 0x000000ff;
            sum >>= 8;
            temp_occ_value1[3] += sum;
        } else if sum == 0x00000100 {
            temp_occ_value1[0] += 256;
        } else if sum == 0x00010000 {
            temp_occ_value1[1] += 256;
        } else if sum == 0x01000000 {
            temp_occ_value1[2] += 256;
        } else {
            temp_occ_value1[3] += 256;
        }
    }

    for c in 0..ALPHABET_SIZE {
        occ_value[occ_index as usize * 4 + c] = ((temp_occ_value0[c] << 16) | temp_occ_value1[c]) as u32;
    }
}

/// Original C static function `bwt_inc_construct` from `minibwa/bwtgen.c:1289`.
pub(crate) fn bwt_inc_construct(bwt_inc: &mut BwtInc, num_char: BgInt) {
    let merged_bwt_size_in_word = bwt_resident_size_in_word(bwt_inc.bwt.text_length + num_char);
    let merged_occ_size_in_word = bwt_occ_value_minor_size_in_word(bwt_inc.bwt.text_length + num_char);
    initialize_val_bg(
        &mut bwt_inc.cumulative_count_in_current_build,
        (ALPHABET_SIZE + 1) as u64,
        0,
    );

    let first_char_in_this_iteration;
    let new_inverse_sa0;
    let mut merged_bwt;

    if bwt_inc.bwt.text_length == 0 {
        let num_items = num_char as usize + 1;
        let mut seq = uninit_vec::<u64>(num_items);
        let mut relative_rank = uninit_vec::<u64>(num_items);

        bwt_inc_put_packed_text_to_rank(
            &bwt_inc.packed_text,
            &mut relative_rank,
            &mut bwt_inc.cumulative_count_in_current_build,
            num_char,
        );

        first_char_in_this_iteration = relative_rank[0] as u32;
        relative_rank[num_char as usize] = 0;

        let qs_v = as_qsint_slice_mut(&mut relative_rank[..num_char as usize + 1]);
        let qs_i = as_qsint_slice_mut(&mut seq[..num_char as usize + 1]);
        suffix_sort(qs_v, qs_i, num_char as i64, ALPHABET_SIZE as i64 - 1, 0, 0);
        new_inverse_sa0 = relative_rank[0];
        drop(seq);

        merged_bwt = uninit_vec::<u32>(merged_bwt_size_in_word as usize + 1);
        initialize_val(&mut merged_bwt, merged_bwt_size_in_word, 0);
        bwt_inc_build_packed_bwt(
            &relative_rank,
            &mut merged_bwt,
            num_char,
            &bwt_inc.cumulative_count_in_current_build,
            &bwt_inc.packed_shift,
        );

        bwt_inc.first_char_in_last_iteration = ALPHABET_SIZE as u32;
    } else {
        let num_items = num_char as usize + 1;
        let mut sorted_rank = uninit_vec::<u64>(num_items);
        let mut seq = uninit_vec::<u64>(num_items);
        let mut relative_rank = uninit_vec::<u64>(num_items);

        first_char_in_this_iteration = bwt_inc.packed_text[0] >> (BITS_IN_WORD - BIT_PER_CHAR);

        forward_dna_all_occ_count_no_limit(
            &bwt_inc.packed_text,
            num_char,
            &mut bwt_inc.cumulative_count_in_current_build[1..],
            &bwt_inc.bwt.decode_table,
        );
        bwt_inc.cumulative_count_in_current_build[bwt_inc.first_char_in_last_iteration as usize + 1] += 1;
        bwt_inc.cumulative_count_in_current_build[2] += bwt_inc.cumulative_count_in_current_build[1];
        bwt_inc.cumulative_count_in_current_build[3] += bwt_inc.cumulative_count_in_current_build[2];
        bwt_inc.cumulative_count_in_current_build[4] += bwt_inc.cumulative_count_in_current_build[3];

        let old_inverse_sa0_relative_rank = bwt_inc_get_absolute_rank(
            &bwt_inc.bwt,
            &mut sorted_rank,
            &mut seq,
            &bwt_inc.packed_text,
            num_char,
            &bwt_inc.cumulative_count_in_current_build,
            bwt_inc.first_char_in_last_iteration,
        );

        for i in 0..ALPHABET_SIZE {
            let lo = bwt_inc.cumulative_count_in_current_build[i];
            let hi = bwt_inc.cumulative_count_in_current_build[i + 1];
            if lo > old_inverse_sa0_relative_rank || hi <= old_inverse_sa0_relative_rank {
                bwt_inc_sort_key(
                    &mut sorted_rank[lo as usize..],
                    &mut seq[lo as usize..],
                    hi - lo,
                );
            } else {
                if lo < old_inverse_sa0_relative_rank {
                    bwt_inc_sort_key(
                        &mut sorted_rank[lo as usize..],
                        &mut seq[lo as usize..],
                        old_inverse_sa0_relative_rank - lo,
                    );
                }
                if hi > old_inverse_sa0_relative_rank + 1 {
                    bwt_inc_sort_key(
                        &mut sorted_rank[(old_inverse_sa0_relative_rank + 1) as usize..],
                        &mut seq[(old_inverse_sa0_relative_rank + 1) as usize..],
                        hi - old_inverse_sa0_relative_rank - 1,
                    );
                }
            }
        }

        bwt_inc_build_relative_rank(
            &mut sorted_rank,
            &mut seq,
            &mut relative_rank,
            num_char,
            bwt_inc.bwt.inverse_sa0,
            &bwt_inc.cumulative_count_in_current_build,
        );
        assert_eq!(relative_rank[num_char as usize], old_inverse_sa0_relative_rank);

        let qs_v = as_qsint_slice_mut(&mut relative_rank[..num_char as usize + 1]);
        let qs_i = as_qsint_slice_mut(&mut seq[..num_char as usize + 1]);
        suffix_sort(qs_v, qs_i, num_char as i64, num_char as i64, 1, 1);
        drop(seq);

        let new_inverse_sa0_relative_rank = relative_rank[0];
        new_inverse_sa0 = sorted_rank[new_inverse_sa0_relative_rank as usize] + new_inverse_sa0_relative_rank;
        sorted_rank[new_inverse_sa0_relative_rank as usize] = 0;

        let mut insert_bwt = uninit_vec::<u32>(num_char as usize + 1);
        bwt_inc_build_bwt(
            &mut insert_bwt,
            &relative_rank,
            num_char,
            &bwt_inc.cumulative_count_in_current_build,
        );
        drop(relative_rank);

        merged_bwt = uninit_vec::<u32>(merged_bwt_size_in_word as usize + 2);
        let mut old_bwt = std::mem::take(&mut bwt_inc.bwt.bwt_code);
        old_bwt.resize(old_bwt.len() + 2, 0);
        bwt_inc_merge_bwt(
            &sorted_rank,
            &old_bwt,
            &insert_bwt,
            &mut merged_bwt,
            bwt_inc.bwt.text_length,
            num_char,
        );
    }

    bwt_inc.bwt.text_length += num_char;
    bwt_inc.bwt.bwt_code = merged_bwt;
    bwt_inc.bwt.bwt_size_in_word = merged_bwt_size_in_word;
    bwt_inc.bwt.occ_size_in_word = merged_occ_size_in_word;
    bwt_inc.bwt.occ_value.clear();
    bwt_inc.bwt.occ_value.resize(merged_occ_size_in_word as usize, 0);

    bwt_clear_trailing_bwt_code(&mut bwt_inc.bwt);
    bwt_generate_occ_value_from_bwt(
        &bwt_inc.bwt.bwt_code,
        &mut bwt_inc.bwt.occ_value,
        &mut bwt_inc.bwt.occ_value_major,
        bwt_inc.bwt.text_length,
        &bwt_inc.bwt.decode_table,
    );

    bwt_inc.bwt.inverse_sa0 = new_inverse_sa0;

    bwt_inc.bwt.cumulative_freq[1] +=
        bwt_inc.cumulative_count_in_current_build[1] - (bwt_inc.first_char_in_last_iteration <= 0) as u64;
    bwt_inc.bwt.cumulative_freq[2] +=
        bwt_inc.cumulative_count_in_current_build[2] - (bwt_inc.first_char_in_last_iteration <= 1) as u64;
    bwt_inc.bwt.cumulative_freq[3] +=
        bwt_inc.cumulative_count_in_current_build[3] - (bwt_inc.first_char_in_last_iteration <= 2) as u64;
    bwt_inc.bwt.cumulative_freq[4] +=
        bwt_inc.cumulative_count_in_current_build[4] - (bwt_inc.first_char_in_last_iteration <= 3) as u64;

    bwt_inc.first_char_in_last_iteration = first_char_in_this_iteration;
    bwt_inc_set_build_size_and_text_addr(bwt_inc);
    bwt_inc.number_of_iteration_done += 1;
}

/// Original C global function `bwt_inc_construct_from_packed` from `minibwa/bwtgen.c:1427`.
pub(crate) fn bwt_inc_construct_from_packed(
    input_file_name: &std::path::Path,
    initial_max_build_size: BgInt,
    inc_max_build_size: BgInt,
) -> std::io::Result<BwtInc> {
    let input_file_name = c_path(input_file_name);
    let mut fp = match std::fs::File::open(&input_file_name) {
        Ok(fp) => fp,
        Err(err) => {
            let reason = err
                .raw_os_error()
                .map(c_strerror)
                .unwrap_or_else(|| err.to_string());
            eprintln!(
                "bwt_inc_construct_from_packed() : Cannot open {} : {}",
                input_file_name.display(),
                reason
            );
            std::process::exit(1);
        }
    };
    use std::io::{Read, Seek};

    if fp.seek(std::io::SeekFrom::End(-1)).is_err() {
        eprintln!(
            "bwt_inc_construct_from_packed() : Can't seek on {} : {}",
            input_file_name.display(),
            c_strerror(einval())
        );
        std::process::exit(1);
    }
    let packed_file_len = match fp.stream_position() {
        Ok(pos) => pos,
        Err(_) => bwt_packed_seek_error_and_exit(&input_file_name),
    };
    let mut last = [0u8; 1];
    if let Err(err) = fp.read_exact(&mut last) {
        bwt_packed_read_error_and_exit(&input_file_name, err);
    }
    let last_byte_length = last[0] as u32;
    let total_text_length = text_length_from_byte_packed(packed_file_len, BIT_PER_CHAR, last_byte_length);

    let mut bwt_inc = bwt_inc_create(
        total_text_length,
        initial_max_build_size as u32,
        inc_max_build_size as u32,
    );
    bwt_inc_set_build_size_and_text_addr(&mut bwt_inc);

    let mut text_to_load = if bwt_inc.build_size > total_text_length {
        total_text_length
    } else {
        total_text_length
            - ((total_text_length - bwt_inc.build_size + CHAR_PER_WORD - 1) / CHAR_PER_WORD
                * CHAR_PER_WORD)
    };
    let mut text_size_in_byte = text_to_load / CHAR_PER_BYTE as u64;
    bwt_inc.packed_text.clear();
    bwt_inc.packed_text.resize(
        (text_to_load as usize + CHAR_PER_WORD as usize - 1) / CHAR_PER_WORD as usize + 1,
        0,
    );
    bwt_inc.text_buffer.resize(text_size_in_byte as usize + 1, 0);
    bwt_packed_seek_back(&mut fp, &input_file_name, text_size_in_byte + 2);
    if let Err(err) = fp.read_exact(&mut bwt_inc.text_buffer) {
        bwt_packed_read_error_and_exit(&input_file_name, err);
    }
    bwt_packed_seek_back(&mut fp, &input_file_name, text_size_in_byte + 1);
    convert_byte_packed_to_word_packed(
        &bwt_inc.text_buffer,
        &mut bwt_inc.packed_text,
        ALPHABET_SIZE as u32,
        text_to_load,
    );
    bwt_inc_construct(&mut bwt_inc, text_to_load);

    let mut processed_text_length = text_to_load;
    while processed_text_length < total_text_length {
        text_to_load = bwt_inc.build_size / CHAR_PER_WORD * CHAR_PER_WORD;
        if text_to_load > total_text_length - processed_text_length {
            text_to_load = total_text_length - processed_text_length;
        }
        text_size_in_byte = text_to_load / CHAR_PER_BYTE as u64;
        bwt_inc.packed_text.clear();
        bwt_inc.packed_text.resize(
            (text_to_load as usize + CHAR_PER_WORD as usize - 1) / CHAR_PER_WORD as usize + 1,
            0,
        );
        bwt_inc.text_buffer.resize(text_size_in_byte as usize, 0);
        bwt_packed_seek_back(&mut fp, &input_file_name, text_size_in_byte);
        if let Err(err) = fp.read_exact(&mut bwt_inc.text_buffer) {
            bwt_packed_read_error_and_exit(&input_file_name, err);
        }
        bwt_packed_seek_back(&mut fp, &input_file_name, text_size_in_byte);
        convert_byte_packed_to_word_packed(
            &bwt_inc.text_buffer,
            &mut bwt_inc.packed_text,
            ALPHABET_SIZE as u32,
            text_to_load,
        );
        bwt_inc_construct(&mut bwt_inc, text_to_load);
        processed_text_length += text_to_load;
        if bwt_inc.number_of_iteration_done % 10 == 0 {
            eprintln!(
                "[bwt_inc_construct_from_packed] {} iterations done. {} characters processed.",
                bwt_inc.number_of_iteration_done, processed_text_length
            );
        }
    }

    Ok(bwt_inc)
}

/// Original C global function `bwt_inc_free` from `minibwa/bwtgen.c:1537`.
pub(crate) fn bwt_inc_free(bwt_inc: Option<BwtInc>) {
    if bwt_inc.is_none() {
        return;
    }
}

/// Original C static function `bwt_file_size_in_word` from `minibwa/bwtgen.c:1550`.
pub(crate) fn bwt_file_size_in_word(num_char: BgInt) -> BgInt {
    (num_char + CHAR_PER_WORD - 1) / CHAR_PER_WORD
}

/// Original C global function `bwt_save_bwt_code_and_occ` from `minibwa/bwtgen.c:1556`.
pub(crate) fn bwt_save_bwt_code_and_occ(
    bwt: &GenBwt,
    bwt_file_name: &std::path::Path,
    occ_value_file_name: Option<&std::path::Path>,
) -> std::io::Result<()> {
    let bwt_file_name = c_path(bwt_file_name);
    let bwt_file = match std::fs::File::create(&bwt_file_name) {
        Ok(file) => file,
        Err(err) => {
            let reason = err
                .raw_os_error()
                .map(c_strerror)
                .unwrap_or_else(|| err.to_string());
            eprintln!(
                "bwt_save_bwt_code_and_occ(): Cannot open {} for writing: {}",
                bwt_file_name.display(),
                reason
            );
            std::process::exit(1);
        }
    };
    let mut bwt_file = std::io::BufWriter::with_capacity(64 * 1024, bwt_file);
    let bwt_length = bwt_file_size_in_word(bwt.text_length);
    use std::io::Write;
    if let Err(err) = bwt_file.write_all(&bwt.inverse_sa0.to_le_bytes()) {
        bwt_save_error_and_exit(&bwt_file_name, err);
    }
    for i in 1..=ALPHABET_SIZE {
        if let Err(err) = bwt_file.write_all(&bwt.cumulative_freq[i].to_le_bytes()) {
            bwt_save_error_and_exit(&bwt_file_name, err);
        }
    }
    for &word in bwt.bwt_code.iter().take(bwt_length as usize) {
        if let Err(err) = bwt_file.write_all(&word.to_le_bytes()) {
            bwt_save_error_and_exit(&bwt_file_name, err);
        }
    }
    if let Err(err) = bwt_file.flush() {
        bwt_save_error_and_exit(&bwt_file_name, err);
    }
    let _ = occ_value_file_name;
    Ok(())
}

/// Original C global function `mb_bwtgen` from `minibwa/bwtgen.c:1588`.
pub(crate) fn mb_bwtgen(
    fn_pac: &std::path::Path,
    fn_bwt: &std::path::Path,
    block_size: i32,
) -> std::io::Result<()> {
    let fn_pac = c_path(fn_pac);
    let fn_bwt = c_path(fn_bwt);
    let bwt_inc = bwt_inc_construct_from_packed(&fn_pac, block_size as u64, block_size as u64)?;
    eprintln!(
        "[bwt_gen] Finished constructing GenBwt in {} iterations.",
        bwt_inc.number_of_iteration_done
    );
    bwt_save_bwt_code_and_occ(&bwt_inc.bwt, &fn_bwt, None)?;
    bwt_inc_free(Some(bwt_inc));
    Ok(())
}
