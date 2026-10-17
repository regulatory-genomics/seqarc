#![allow(unused_variables)]

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use crate::minibwa::mapped_u64::U64Storage;

const MB_MAGIC: &[u8; 4] = b"MBW\x02";
const BWT_CNT_SHIFT: u32 = 56;
const BWT_CNT_MASK: u64 = (1u64 << BWT_CNT_SHIFT) - 1;

fn c_path<P: AsRef<Path>>(path: P) -> PathBuf {
    let lossy = path.as_ref().to_string_lossy();
    let end = lossy
        .as_bytes()
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(lossy.len());
    PathBuf::from(lossy[..end].to_string())
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct SeedHit {
    pub x: [u64; 2],
    pub size: u64,
    pub info: u64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Bwt {
    pub primary: u64,
    pub l2: [u64; 5],
    pub seq_len: u64,
    pub data_len: u64,
    pub data: U64Storage,
    pub cnt_table: [u32; 256],
    pub pre_len: u32,
    pub pre: Vec<SeedHit>,
    pub sa_bit: u32,
    pub n_sa: u64,
    pub sa: U64Storage,
}

#[repr(C)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct SeedHitVec {
    pub n: usize,
    pub m: usize,
    pub a: Vec<SeedHit>,
}

#[inline(always)]
pub(crate) fn mb_sai_v_clear(v: &mut SeedHitVec) {
    v.n = 0;
    v.a.clear();
    v.m = v.a.capacity();
}

#[inline(always)]
fn mb_sai_v_grow(v: &mut SeedHitVec) {
    v.m = v.n + 1;
    v.m += (v.m >> 1) + 16;
    v.a.reserve(v.m.saturating_sub(v.a.capacity()));
    v.m = v.a.capacity();
}

#[inline(always)]
fn mb_sai_v_push_reserved(v: &mut SeedHitVec, x: SeedHit) {
    if v.n >= v.a.capacity() {
        mb_sai_v_grow(v);
    }
    let n = v.n;
    unsafe {
        v.a.as_mut_ptr().add(n).write(x);
        if n == v.a.len() {
            v.a.set_len(n + 1);
        }
    }
    v.n = n + 1;
}

#[repr(C)]
#[derive(Debug)]
pub(crate) struct SmemEntryRef {
    pub min_len: i32,
    pub min_occ: i32,
    pub st: i32,
    pub en: i32,
    pub q: *const u8,
    pub v: *mut SeedHitVec,
    pub stage: i32,
    pub x: i32,
    pub i: i32,
    pub kmer: u32,
    pub p: SeedHit,
}

#[inline(always)]
unsafe fn smem_q_ref(s: &SmemEntryRef, i: i32) -> u8 {
    unsafe { *s.q.add(i as usize) }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct TinyQueue {
    pub front: i32,
    pub count: i32,
    pub cap: i32,
    pub a: Vec<i32>,
}

/// Original C static function `bwt_gen_cnt_table` from `minibwa/bwt.c:14`.
pub(crate) fn bwt_gen_cnt_table(cnt: &mut [u32; 256]) {
    for i in 0..256usize {
        let mut x = 0u32;
        for j in 0..4u32 {
            x |= (((i & 3) as u32 == j) as u32
                + (((i >> 2) & 3) as u32 == j) as u32
                + (((i >> 4) & 3) as u32 == j) as u32
                + ((i >> 6) as u32 == j) as u32)
                << (j << 3);
        }
        cnt[i] = x;
    }
}

/// Original C global function `mb_bwt_init` from `minibwa/bwt.c:25`.
pub(crate) fn mb_bwt_init() -> Bwt {
    let mut bwt = Bwt {
        primary: 0,
        l2: [0; 5],
        seq_len: 0,
        data_len: 0,
        data: U64Storage::new(),
        cnt_table: [0; 256],
        pre_len: 0,
        pre: Vec::new(),
        sa_bit: u32::MAX,
        n_sa: 0,
        sa: U64Storage::new(),
    };
    bwt_gen_cnt_table(&mut bwt.cnt_table);
    bwt
}

/// Original C static function `mb_bwt_data_len` from `minibwa/bwt.c:50`.
pub(crate) fn mb_bwt_data_len(len: u64) -> u64 {
    let bwt_len = (len + 127) / 128 * 4;
    let occ_len = ((len + 127) / 128 + 1) * 4;
    bwt_len + occ_len
}

/// GenBwt layout. Each block consists of u64[4]+u32[8], 64 bytes in total. The
/// lower 56 bits of each u64[4] (see BWT_CNT_SHIFT) store the accumulative
/// count of A/C/G/T bases. The higher 8 bits store the count of A/C/G/T in the
/// next 64nt. u32[8] keeps a GenBwt substring of 128nt in length. Because it
/// follows little endian, it can also be considered as u64[4] etc.
///
/// Original C global function `mb_bwt_init_from_raw` from `minibwa/bwt.c:64`.
pub(crate) fn mb_bwt_init_from_raw(is_byte: i32, raw: &[u8], len: u64, primary: u64) -> Bwt {
    let mut c = [0u64; 4];
    let mut x = [0u64; 4];
    let mut bwt = mb_bwt_init();
    bwt.primary = primary;
    bwt.seq_len = len;
    bwt.data_len = mb_bwt_data_len(len);
    bwt.data = U64Storage::with_capacity(bwt.data_len as usize);
    if len == 0 {
        bwt.data.resize(bwt.data_len as usize, 0);
        return bwt;
    }

    let mut last_c: Option<usize> = None;
    for i in 0..len {
        let a = if is_byte != 0 {
            raw[i as usize] & 3
        } else {
            let word_start = ((i >> 4) as usize) * 4;
            let b = u32::from_ne_bytes([
                raw[word_start],
                raw[word_start + 1],
                raw[word_start + 2],
                raw[word_start + 3],
            ]);
            ((b >> ((!i & 0xf) << 1)) & 3) as u8
        };
        if (i & 0x7f) == 0 {
            if i > 0 {
                bwt.data.extend_from_slice(&x);
            }
            last_c = Some(bwt.data.len());
            bwt.data.extend_from_slice(&c);
            x = [0; 4];
        } else if (i & 0x3f) == 0 {
            if let Some(last_c) = last_c {
                for j in 0..4usize {
                    bwt.data[last_c + j] |= (c[j] - bwt.data[last_c + j]) << BWT_CNT_SHIFT;
                }
            }
        }
        c[a as usize] += 1;
        x[((i & 0x7f) >> 5) as usize] |= (a as u64) << ((i & 0x1f) << 1);
    }
    bwt.data.extend_from_slice(&x);
    bwt.data.extend_from_slice(&c);
    if bwt.data.len() != bwt.data_len as usize {
        c_assert_fail_bwt_init_from_raw(101, "k == bwt->data_len");
    }
    bwt.l2[0] = 0;
    for i in 0..4usize {
        bwt.l2[i + 1] = bwt.l2[i] + c[i];
    }
    if bwt.l2[4] != len {
        c_assert_fail_bwt_init_from_raw(104, "bwt->L2[4] == len");
    }
    bwt
}

fn c_assert_fail_bwt_init_from_raw(line: u32, assertion: &str) -> ! {
    eprintln!("minibwa: bwt.c:{line}: mb_bwt_init_from_raw: Assertion `{assertion}' failed.");
    std::process::abort();
}

/// Original C static function `mb_bwt_set_intv` from `minibwa/bwt.h:71`.
#[inline(always)]
pub(crate) fn mb_bwt_set_intv(bwt: &Bwt, c: i32, ik: &mut SeedHit) {
    let c = c as usize;
    ik.x[0] = unsafe { *bwt.l2.get_unchecked(c) } + 1;
    ik.x[1] = unsafe { *bwt.l2.get_unchecked(3 - c) } + 1;
    ik.size = unsafe { *bwt.l2.get_unchecked(c + 1) } - unsafe { *bwt.l2.get_unchecked(c) };
    ik.info = 0;
}

/// Original C static function `mb_bwt_block_prefetch` from `minibwa/bwt.c:118`.
#[inline(always)]
pub(crate) fn mb_bwt_block_prefetch(bwt: &Bwt, k: u64) {
    if k > 0 && !bwt.data.is_empty() {
        let block = ((k - 1 - ((k - 1 >= bwt.primary) as u64)) >> 7) << 3;
        let cell = bwt.data.as_ptr().wrapping_add(block as usize);
        crate::minibwa::s2n_lite::_mm_prefetch(cell as *const u8, 3);
    }
}

#[inline(always)]
pub(crate) fn mb_bwt_pre_prefetch(bwt: &Bwt, kmer: u32) {
    if !bwt.pre.is_empty() {
        let p = unsafe { bwt.pre.as_ptr().add(kmer as usize) };
        crate::minibwa::s2n_lite::_mm_prefetch(p as *const u8, 3);
    }
}

/// Original C static function `rank_aux1` from `minibwa/bwt.c:123`.
#[inline(always)]
pub(crate) fn rank_aux1(mut y: u64, c: u8) -> i32 {
    y = ((if (c & 2) != 0 { y } else { !y }) >> 1)
        & (if (c & 1) != 0 { y } else { !y })
        & 0x5555_5555_5555_5555u64;
    y.count_ones() as i32
}

/// Original C global function `mb_bwt_rank11` from `minibwa/bwt.c:136`.
#[inline(always)]
pub(crate) fn mb_bwt_rank11(bwt: &Bwt, mut k: u64, c: u8) -> u64 {
    if k == 0 {
        return 0;
    }
    if k == bwt.seq_len + 1 {
        return unsafe { *bwt.l2.get_unchecked(c as usize + 1) }
            - unsafe { *bwt.l2.get_unchecked(c as usize) };
    }
    k -= 1;
    k -= (k >= bwt.primary) as u64;
    let mask = if (k & 0x7f) >= 64 {
        (1u64 << (64 - BWT_CNT_SHIFT)) - 1
    } else {
        0
    };
    let block = ((k >> 7) << 3) as usize;
    let word = unsafe { *bwt.data.get_unchecked(block + c as usize) };
    let mut n = (word & BWT_CNT_MASK) + ((word >> BWT_CNT_SHIFT) & mask);
    let mut p = block + 4;
    let end = p + ((k & 0x7f) >> 5) as usize;
    if (k & 0x7f) >= 64 {
        p += 2;
    }
    if p < end {
        n += rank_aux1(unsafe { *bwt.data.get_unchecked(p) }, c) as u64;
        p += 1;
    }
    n += rank_aux1(
        unsafe { *bwt.data.get_unchecked(p) } << ((!k & 0x1f) << 1),
        c,
    ) as u64;
    if c == 0 {
        n -= !k & 0x1f;
    }
    n
}

/// Original C static function `seek_block` from `minibwa/bwt.c:157`.
#[inline(always)]
pub(crate) fn seek_block(bwt: &Bwt, k: u64, cnt: &mut [u64; 4]) -> usize {
    let p = ((k >> 7) << 3) as usize;
    let mask = if (k & 0x7f) >= 64 {
        (1u64 << (64 - BWT_CNT_SHIFT)) - 1
    } else {
        0
    };
    let data = bwt.data.as_ptr();
    let x0 = unsafe { *data.add(p) };
    let x1 = unsafe { *data.add(p + 1) };
    let x2 = unsafe { *data.add(p + 2) };
    let x3 = unsafe { *data.add(p + 3) };
    cnt[0] = (x0 & BWT_CNT_MASK) + ((x0 >> BWT_CNT_SHIFT) & mask);
    cnt[1] = (x1 & BWT_CNT_MASK) + ((x1 >> BWT_CNT_SHIFT) & mask);
    cnt[2] = (x2 & BWT_CNT_MASK) + ((x2 >> BWT_CNT_SHIFT) & mask);
    cnt[3] = (x3 & BWT_CNT_MASK) + ((x3 >> BWT_CNT_SHIFT) & mask);
    (p + 4) * 2
}

/// Original C static function `rank_aux4` from `minibwa/bwt.c:168`.
#[inline(always)]
pub(crate) fn rank_aux4(bwt: &Bwt, x: u32) -> u32 {
    unsafe {
        *bwt.cnt_table.get_unchecked((x & 0xff) as usize)
            + *bwt.cnt_table.get_unchecked(((x >> 8) & 0xff) as usize)
            + *bwt.cnt_table.get_unchecked(((x >> 16) & 0xff) as usize)
            + *bwt.cnt_table.get_unchecked((x >> 24) as usize)
    }
}

#[inline(always)]
unsafe fn bwt_word32(bwt: &Bwt, q: usize) -> u32 {
    unsafe { *(bwt.data.as_ptr() as *const u32).add(q) }
}

/// Original C global function `mb_bwt_rank1a` from `minibwa/bwt.c:173`.
#[inline(always)]
pub(crate) fn mb_bwt_rank1a(bwt: &Bwt, mut k: u64, cnt: &mut [u64; 4]) {
    if k == 0 {
        *cnt = [0; 4];
        return;
    }
    k -= 1;
    k -= (k >= bwt.primary) as u64;
    let mut q = seek_block(bwt, k, cnt);
    let end = q + ((k & 0x7f) >> 4) as usize;
    if (k & 0x7f) >= 64 {
        q += 4;
    }
    let mut x = 0u32;
    while q < end {
        x = x.wrapping_add(rank_aux4(bwt, unsafe { bwt_word32(bwt, q) }));
        q += 1;
    }
    let val = unsafe { bwt_word32(bwt, q) };
    let tmp = val << ((!k & 0xf) << 1);
    x = x.wrapping_add(rank_aux4(bwt, tmp).wrapping_sub((!k & 0xf) as u32));
    cnt[0] += (x & 0xff) as u64;
    cnt[1] += ((x >> 8) & 0xff) as u64;
    cnt[2] += ((x >> 16) & 0xff) as u64;
    cnt[3] += (x >> 24) as u64;
}

/// Original C global function `mb_bwt_rank2a` from `minibwa/bwt.c:192`.
#[inline(always)]
pub(crate) fn mb_bwt_rank2a(
    bwt: &Bwt,
    mut k: u64,
    mut l: u64,
    cntk: &mut [u64; 4],
    cntl: &mut [u64; 4],
) {
    let mut k1 = k.wrapping_sub(1);
    let mut l1 = l.wrapping_sub(1);
    k1 -= (k1 >= bwt.primary) as u64;
    l1 -= (l1 >= bwt.primary) as u64;
    mb_bwt_block_prefetch(bwt, k);
    if (k1 >> 7) != (l1 >> 7) || k == 0 || l == 0 {
        mb_bwt_block_prefetch(bwt, l);
        mb_bwt_rank1a(bwt, k, cntk);
        mb_bwt_rank1a(bwt, l, cntl);
    } else if l - k == 1 {
        let z = k - (k > bwt.primary) as u64;
        mb_bwt_rank1a(bwt, k, cntk);
        *cntl = *cntk;
        let base = (((z >> 7) << 3) + 4 + (((z & 127) >> 5) as u64)) as usize;
        let c = ((unsafe { *bwt.data.get_unchecked(base) } >> ((z & 31) << 1)) & 3) as usize;
        cntl[c] += 1;
    } else {
        k = k1;
        l = l1;
        let mut q = seek_block(bwt, k, cntk);
        let endk = q + ((k & 0x7f) >> 4) as usize;
        let endl = q + ((l & 0x7f) >> 4) as usize;
        if (k & 0x7f) >= 64 {
            q += 4;
        }
        let mut x = 0u32;
        while q < endk {
            x = x.wrapping_add(rank_aux4(bwt, unsafe { bwt_word32(bwt, q) }));
            q += 1;
        }
        let mut y = x;
        let val = unsafe { bwt_word32(bwt, q) };
        let tmp = val << ((!k & 0xf) << 1);
        x = x.wrapping_add(rank_aux4(bwt, tmp).wrapping_sub((!k & 0xf) as u32));
        while q < endl {
            y = y.wrapping_add(rank_aux4(bwt, unsafe { bwt_word32(bwt, q) }));
            q += 1;
        }
        let val = unsafe { bwt_word32(bwt, q) };
        let tmp = val << ((!l & 0xf) << 1);
        y = y.wrapping_add(rank_aux4(bwt, tmp).wrapping_sub((!l & 0xf) as u32));
        *cntl = *cntk;
        cntk[0] += (x & 0xff) as u64;
        cntk[1] += ((x >> 8) & 0xff) as u64;
        cntk[2] += ((x >> 16) & 0xff) as u64;
        cntk[3] += (x >> 24) as u64;
        cntl[0] += (y & 0xff) as u64;
        cntl[1] += ((y >> 8) & 0xff) as u64;
        cntl[2] += ((y >> 16) & 0xff) as u64;
        cntl[3] += (y >> 24) as u64;
    }
}

#[inline(always)]
pub(crate) fn mb_bwt_extend_forward_one(bwt: &Bwt, ik: &SeedHit, c: usize) -> SeedHit {
    debug_assert!(c < 4);
    let mut tk = [0u64; 4];
    let mut tl = [0u64; 4];
    mb_bwt_rank2a(bwt, ik.x[1], ik.x[1] + ik.size, &mut tk, &mut tl);
    tl[0] -= tk[0];
    tl[1] -= tk[1];
    tl[2] -= tk[2];
    tl[3] -= tk[3];
    let x0 = ik.x[0] + ((ik.x[1] <= bwt.primary && ik.x[1] + ik.size > bwt.primary) as u64);
    let x0 = match c {
        0 => x0 + tl[3] + tl[2] + tl[1],
        1 => x0 + tl[3] + tl[2],
        2 => x0 + tl[3],
        _ => x0,
    };
    SeedHit {
        x: [
            x0,
            unsafe { *bwt.l2.get_unchecked(c) } + 1 + unsafe { *tk.get_unchecked(c) },
        ],
        size: unsafe { *tl.get_unchecked(c) },
        info: 0,
    }
}

#[inline(always)]
pub(crate) fn mb_bwt_extend_back(bwt: &Bwt, ik: &SeedHit, ok: &mut [SeedHit; 4]) {
    let mut tk = [0u64; 4];
    let mut tl = [0u64; 4];
    mb_bwt_rank2a(bwt, ik.x[0], ik.x[0] + ik.size, &mut tk, &mut tl);
    ok[0].x[0] = unsafe { *bwt.l2.get_unchecked(0) } + 1 + tk[0];
    ok[1].x[0] = unsafe { *bwt.l2.get_unchecked(1) } + 1 + tk[1];
    ok[2].x[0] = unsafe { *bwt.l2.get_unchecked(2) } + 1 + tk[2];
    ok[3].x[0] = unsafe { *bwt.l2.get_unchecked(3) } + 1 + tk[3];
    tl[0] -= tk[0];
    tl[1] -= tk[1];
    tl[2] -= tk[2];
    tl[3] -= tk[3];
    ok[0].size = tl[0];
    ok[1].size = tl[1];
    ok[2].size = tl[2];
    ok[3].size = tl[3];
    ok[3].x[1] = ik.x[1] + ((ik.x[0] <= bwt.primary && ik.x[0] + ik.size > bwt.primary) as u64);
    ok[2].x[1] = ok[3].x[1] + tl[3];
    ok[1].x[1] = ok[2].x[1] + tl[2];
    ok[0].x[1] = ok[1].x[1] + tl[1];
}

#[inline(always)]
pub(crate) fn mb_bwt_extend_back_one(bwt: &Bwt, ik: &SeedHit, c: usize) -> SeedHit {
    debug_assert!(c < 4);
    let mut tk = [0u64; 4];
    let mut tl = [0u64; 4];
    mb_bwt_rank2a(bwt, ik.x[0], ik.x[0] + ik.size, &mut tk, &mut tl);
    tl[0] -= tk[0];
    tl[1] -= tk[1];
    tl[2] -= tk[2];
    tl[3] -= tk[3];
    let x1 = ik.x[1] + ((ik.x[0] <= bwt.primary && ik.x[0] + ik.size > bwt.primary) as u64);
    let x1 = match c {
        0 => x1 + tl[3] + tl[2] + tl[1],
        1 => x1 + tl[3] + tl[2],
        2 => x1 + tl[3],
        _ => x1,
    };
    SeedHit {
        x: [
            unsafe { *bwt.l2.get_unchecked(c) } + 1 + unsafe { *tk.get_unchecked(c) },
            x1,
        ],
        size: unsafe { *tl.get_unchecked(c) },
        info: 0,
    }
}

#[inline(always)]
pub(crate) fn tq_reset(q: &mut TinyQueue, n: i32) {
    let mut cap = n;
    if cap <= 0 {
        q.cap = 0;
        q.front = 0;
        q.count = 0;
        return;
    }
    cap -= 1;
    cap |= cap >> 1;
    cap |= cap >> 2;
    cap |= cap >> 4;
    cap |= cap >> 8;
    cap |= cap >> 16;
    cap += 1;
    q.cap = cap;
    if q.a.len() < q.cap as usize {
        q.a.resize(q.cap as usize, 0);
    }
    q.front = 0;
    q.count = 0;
}

/// Original C static function `tq_push` from `minibwa/bwt.c:321`.
#[inline(always)]
pub(crate) fn tq_push(q: &mut TinyQueue, x: i32) {
    debug_assert!(q.cap > 0);
    debug_assert!(q.count < q.cap);
    let idx = ((q.count + q.front) & (q.cap - 1)) as usize;
    debug_assert!(idx < q.a.len());
    unsafe {
        *q.a.get_unchecked_mut(idx) = x;
    }
    q.count += 1;
}

/// Original C static function `tq_shift` from `minibwa/bwt.c:326`.
#[inline(always)]
pub(crate) fn tq_shift(q: &mut TinyQueue) -> i32 {
    if q.count == 0 {
        return -1;
    }
    debug_assert!(q.cap > 0);
    let idx = q.front as usize;
    debug_assert!(idx < q.a.len());
    let x = unsafe { *q.a.get_unchecked(idx) };
    q.front += 1;
    q.front &= q.cap - 1;
    q.count -= 1;
    x
}

#[inline(always)]
pub(crate) fn se_one_step_back_ref(bwt: &Bwt, s: &mut SmemEntryRef) {
    let c = unsafe { smem_q_ref(s, s.i) } as usize;
    debug_assert!(c < 4);
    let child = mb_bwt_extend_back_one(bwt, &s.p, c);
    if child.size < s.min_occ as u64 {
        s.x = s.i + 1;
        s.stage = 1;
    } else {
        s.p = child;
        s.i -= 1;
        mb_bwt_block_prefetch(bwt, s.p.x[0]);
        mb_bwt_block_prefetch(bwt, s.p.x[0] + s.p.size);
    }
}

pub(crate) fn mb_bwt_smem_batch_ref_with_queue(
    bwt: &Bwt,
    n: i32,
    a: &mut [SmemEntryRef],
    _v: &mut [SeedHitVec],
    tq: &mut TinyQueue,
) {
    tq_reset(tq, n);
    for i in 0..n {
        let s = unsafe { a.get_unchecked_mut(i as usize) };
        tq_push(tq, i);
        s.stage = 1;
        s.x = s.st;
        let sv = unsafe { &mut *s.v };
        if sv.m < 64 {
            sv.m = 64;
            sv.a.reserve(sv.m.saturating_sub(sv.a.capacity()));
            sv.m = sv.a.capacity();
        }
    }

    while tq.count > 0 {
        let idx = tq_shift(tq);
        let s = unsafe { a.get_unchecked_mut(idx as usize) };
        if s.stage == 1 {
            if s.en - s.x < s.min_len {
                continue;
            }
            let mut xn = -1;
            let mut i = s.x;
            while i < s.x + s.min_len {
                if unsafe { smem_q_ref(s, i) } > 3 {
                    xn = i;
                }
                i += 1;
            }
            if xn >= 0 {
                s.x = xn + 1;
            } else {
                s.i = s.x + s.min_len - 1;
                if !bwt.pre.is_empty() && s.min_len >= bwt.pre_len as i32 {
                    s.kmer = 0;
                    let mut i = 0;
                    while i < bwt.pre_len as i32 {
                        s.kmer = s.kmer << 2 | unsafe { smem_q_ref(s, s.i) } as u32;
                        s.i -= 1;
                        i += 1;
                    }
                    mb_bwt_pre_prefetch(bwt, s.kmer);
                    s.stage = 2;
                } else {
                    mb_bwt_set_intv(bwt, unsafe { smem_q_ref(s, s.i) } as i32, &mut s.p);
                    s.i -= 1;
                    s.stage = 3;
                }
            }
        } else if s.stage == 2 || s.stage == 5 {
            s.p = unsafe { *bwt.pre.get_unchecked(s.kmer as usize) };
            if s.p.size < s.min_occ as u64 {
                s.i += bwt.pre_len as i32;
                mb_bwt_set_intv(bwt, unsafe { smem_q_ref(s, s.i) } as i32, &mut s.p);
                s.i -= 1;
            }
            s.stage += 1;
        } else if s.stage == 3 {
            if s.i < s.x {
                mb_bwt_block_prefetch(bwt, s.p.x[1]);
                mb_bwt_block_prefetch(bwt, s.p.x[1] + s.p.size);
                s.i = s.x + s.min_len;
                s.stage = 4;
            } else {
                se_one_step_back_ref(bwt, s);
            }
        } else if s.stage == 4 {
            if s.i == s.en {
                s.p.info = (s.x as u64) << 32 | s.i as u64;
                let sv = unsafe { &mut *s.v };
                mb_sai_v_push_reserved(sv, s.p);
                continue;
            } else {
                let c = 3 - unsafe { smem_q_ref(s, s.i) } as i32;
                if c >= 0 && {
                    let child = mb_bwt_extend_forward_one(bwt, &s.p, c as usize);
                    if child.size >= s.min_occ as u64 {
                        s.p = child;
                        true
                    } else {
                        false
                    }
                } {
                    s.i += 1;
                    mb_bwt_block_prefetch(bwt, s.p.x[1]);
                    mb_bwt_block_prefetch(bwt, s.p.x[1] + s.p.size);
                } else {
                    s.p.info = (s.x as u64) << 32 | s.i as u64;
                    let sv = unsafe { &mut *s.v };
                    mb_sai_v_push_reserved(sv, s.p);
                    if c < 0 {
                        s.x = s.i + 1;
                        s.stage = 1;
                    } else if !bwt.pre.is_empty() && s.i - s.x - 1 >= bwt.pre_len as i32 {
                        s.kmer = 0;
                        let mut i = 0;
                        while i < bwt.pre_len as i32 {
                            s.kmer = s.kmer << 2 | unsafe { smem_q_ref(s, s.i) } as u32;
                            s.i -= 1;
                            i += 1;
                        }
                        mb_bwt_pre_prefetch(bwt, s.kmer);
                        s.stage = 5;
                    } else {
                        mb_bwt_set_intv(bwt, unsafe { smem_q_ref(s, s.i) } as i32, &mut s.p);
                        s.i -= 1;
                        s.stage = 6;
                    }
                }
            }
        } else if s.stage == 6 {
            if s.i < s.x + 1 {
                s.x = s.i + 1;
                s.stage = 1;
            } else {
                se_one_step_back_ref(bwt, s);
            }
        }
        tq_push(tq, idx);
    }
}

/// Compute inverse CSA.
///
/// Original C static function `bwt_inv_psi` from `minibwa/bwt.c:460`.
#[inline(always)]
pub(crate) fn bwt_inv_psi(bwt: &Bwt, k: u64) -> u64 {
    let x = k - (k > bwt.primary) as u64;
    let block = (((x >> 7) << 3) + 4 + ((x & 127) >> 5)) as usize;
    let c = ((unsafe { *bwt.data.get_unchecked(block) } >> ((x & 31) << 1)) & 3) as u8;
    let next = unsafe { *bwt.l2.get_unchecked(c as usize) } + 1 + mb_bwt_rank11(bwt, k, c);
    if k == bwt.primary {
        0
    } else {
        next
    }
}

/// `bwt->bwt` and `bwt->occ` must be precalculated.
///
/// Original C global function `mb_bwt_gen_sa` from `minibwa/bwt.c:469`.
pub(crate) fn mb_bwt_gen_sa(bwt: &mut Bwt, sa_bit: u32) {
    assert!(!bwt.data.is_empty());
    drop(std::mem::take(&mut bwt.sa));
    bwt.sa_bit = sa_bit;
    bwt.n_sa = (bwt.seq_len + (1u64 << sa_bit)) >> sa_bit;
    bwt.sa = uninit_u64_vec(bwt.n_sa as usize).into();

    let mut isa = 0u64;
    let mut sa = bwt.seq_len;
    let mask = (1u64 << sa_bit) - 1;
    for _ in 0..bwt.seq_len {
        if (isa & mask) == 0 {
            bwt.sa[(isa >> bwt.sa_bit) as usize] = sa;
        }
        sa = sa.wrapping_sub(1);
        isa = bwt_inv_psi(bwt, isa);
    }
    if (isa & mask) == 0 {
        bwt.sa[(isa >> bwt.sa_bit) as usize] = sa;
    }
    bwt.sa[0] = u64::MAX;
}

#[inline]
pub(crate) fn mb_bwt_sa_batch_with_scratch(
    bwt: &Bwt,
    n: i64,
    x: &mut [u64],
    z: &mut Vec<(u64, u64)>,
) {
    let mask = (1u64 << bwt.sa_bit) - 1;
    if n <= 0 {
        return;
    }
    z.clear();
    z.resize(n as usize, (0, 0));
    let n = n as usize;
    let z_ptr = z.as_mut_ptr();
    let x_ptr = x.as_mut_ptr();
    let sa_ptr = bwt.sa.as_ptr();
    for i in 0..n {
        let x_i = unsafe { *x_ptr.add(i) };
        unsafe {
            (*z_ptr.add(i)).0 = x_i;
            (*z_ptr.add(i)).1 = i as u64;
        }
        if (x_i & mask) == 0 {
            let sa = unsafe { sa_ptr.add((x_i >> bwt.sa_bit) as usize) };
            crate::minibwa::s2n_lite::_mm_prefetch(sa as *const u8, 3);
        } else {
            mb_bwt_block_prefetch(bwt, x_i);
        }
    }
    let mut step = 0u64;
    let mut r = n;
    while r > 0 {
        let r0 = r;
        r = 0;
        for i in 0..r0 {
            let zi = unsafe { *z_ptr.add(i) };
            if (zi.0 & mask) == 0 {
                unsafe {
                    *x_ptr.add(zi.1 as usize) =
                        step.wrapping_add(*sa_ptr.add((zi.0 >> bwt.sa_bit) as usize));
                }
            } else {
                unsafe { *z_ptr.add(r) = zi };
                r += 1;
            }
        }
        for i in 0..r {
            let x_i = bwt_inv_psi(bwt, unsafe { (*z_ptr.add(i)).0 });
            unsafe {
                (*z_ptr.add(i)).0 = x_i;
            }
            if (x_i & mask) == 0 {
                let sa = unsafe { sa_ptr.add((x_i >> bwt.sa_bit) as usize) };
                crate::minibwa::s2n_lite::_mm_prefetch(sa as *const u8, 3);
            } else {
                mb_bwt_block_prefetch(bwt, x_i);
            }
        }
        step += 1;
    }
}

fn mb_bwt_count_kmer_uninit(bwt: &Bwt, depth: i32, s: &mut [std::mem::MaybeUninit<SeedHit>]) {
    #[derive(Clone, Copy, Debug, Default)]
    struct CountPair64 {
        p: SeedHit,
        d: i32,
        c: i32,
    }

    let mut stack = [CountPair64::default(); 64];
    let mut s_top = 0usize;
    let mut str_ = [0u8; 16];
    assert!(depth <= 15);
    for a in 0..4i32 {
        mb_bwt_set_intv(bwt, a, &mut stack[s_top].p);
        stack[s_top].d = 1;
        stack[s_top].c = a;
        s_top += 1;
    }
    while s_top > 0 {
        s_top -= 1;
        let top = stack[s_top];
        let mut ok = [SeedHit::default(); 4];
        if top.d > 0 {
            str_[(depth - top.d) as usize] = top.c as u8;
        }
        mb_bwt_extend_back(bwt, &top.p, &mut ok);
        for a in 0..4i32 {
            str_[(depth - top.d - 1) as usize] = a as u8;
            if top.d != depth - 1 {
                stack[s_top].p = ok[a as usize];
                stack[s_top].d = top.d + 1;
                stack[s_top].c = a;
                s_top += 1;
            } else {
                let mut x = 0u64;
                for i in 0..depth {
                    x |= (str_[i as usize] as u64) << (i * 2);
                }
                s[x as usize].write(ok[a as usize]);
            }
        }
    }
}

/// Original C global function `mb_bwt_cache` from `minibwa/bwt.c:576`.
pub(crate) fn mb_bwt_cache(bwt: &mut Bwt, len: i32) {
    bwt.pre.clear();
    bwt.pre_len = len as u32;
    let n = 1usize << (len * 2);
    let mut pre = Vec::<std::mem::MaybeUninit<SeedHit>>::with_capacity(n);
    let pre_slice = unsafe { std::slice::from_raw_parts_mut(pre.as_mut_ptr(), n) };
    mb_bwt_count_kmer_uninit(bwt, len, pre_slice);
    let ptr = pre.as_mut_ptr() as *mut SeedHit;
    let cap = pre.capacity();
    std::mem::forget(pre);
    let pre = unsafe { Vec::from_raw_parts(ptr, n, cap) };
    bwt.pre = pre;
}

/// Original C static function `read_huge` from `minibwa/bwt.c:588`.
pub(crate) fn read_huge<R: Read>(fp: &mut R, size: u64, a: &mut [u8]) -> u64 {
    let bufsize = 0x1000000usize;
    let mut offset = 0usize;
    let mut remaining = size as usize;
    while remaining > 0 {
        let x = bufsize.min(remaining);
        match fp.read(&mut a[offset..offset + x]) {
            Ok(0) | Err(_) => break,
            Ok(n) => {
                remaining -= n;
                offset += n;
            }
        }
    }
    offset as u64
}

fn uninit_u64_vec(len: usize) -> Vec<u64> {
    let mut a = Vec::<std::mem::MaybeUninit<u64>>::with_capacity(len);
    unsafe {
        a.set_len(len);
        let ptr = a.as_mut_ptr() as *mut u64;
        let cap = a.capacity();
        std::mem::forget(a);
        Vec::from_raw_parts(ptr, len, cap)
    }
}

/// Original C global function `mb_bwt_load_raw` from `minibwa/bwt.c:600`.
pub(crate) fn mb_bwt_load_raw<P: AsRef<Path>>(fn_: P) -> Option<Bwt> {
    let mut fp = File::open(c_path(fn_)).ok()?;
    let file_len = fp.metadata().ok()?.len();
    let raw_size = file_len
        .checked_sub(std::mem::size_of::<u64>() as u64 * 5)
        .map(|n| n >> 2)
        .unwrap_or(0);
    let mut hdr = [0u8; 40];
    let _ = fp.read(&mut hdr);
    let primary = u64::from_le_bytes(hdr[0..8].try_into().unwrap());
    let mut l2 = [0u64; 5];
    for i in 1..5usize {
        let start = i * 8;
        l2[i] = u64::from_le_bytes(hdr[start..start + 8].try_into().unwrap());
    }
    let mut raw = vec![0u8; (raw_size << 2) as usize];
    read_huge(&mut fp, raw_size << 2, &mut raw);
    Some(mb_bwt_init_from_raw(0, &raw, l2[4], primary))
}

/// Original C global function `mb_bwt_save` from `minibwa/bwt.c:622`.
pub(crate) fn mb_bwt_save<P: AsRef<Path>>(fn_: P, bwt: &Bwt) -> i32 {
    let fp = match File::create(c_path(fn_)) {
        Ok(fp) => fp,
        Err(_) => return -1,
    };
    let mut fp = fp;
    let _ = fp.write_all(MB_MAGIC);
    let _ = fp.write_all(&bwt.sa_bit.to_le_bytes());
    let _ = fp.write_all(&bwt.primary.to_le_bytes());
    for i in 1..5usize {
        let _ = fp.write_all(&bwt.l2[i].to_le_bytes());
    }
    let _ = write_u64_slice_le(&mut fp, &bwt.data[..bwt.data_len as usize]);
    let _ = fp.write_all(&bwt.n_sa.to_le_bytes());
    if bwt.sa_bit != u32::MAX && bwt.n_sa > 0 && !bwt.sa.is_empty() {
        let _ = write_u64_slice_le(&mut fp, &bwt.sa[..bwt.n_sa as usize]);
    }
    let _ = fp.flush();
    0
}

fn write_u64_slice_le<W: Write>(writer: &mut W, words: &[u64]) -> std::io::Result<()> {
    #[cfg(target_endian = "little")]
    {
        let bytes = unsafe {
            std::slice::from_raw_parts(words.as_ptr().cast::<u8>(), std::mem::size_of_val(words))
        };
        writer.write_all(bytes)
    }
    #[cfg(not(target_endian = "little"))]
    {
        for &word in words {
            writer.write_all(&word.to_le_bytes())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_matches_bruteforce_for_constructed_byte_bwt() {
        let raw: Vec<u8> = (0..257).map(|i| ((i * 7 + i / 3) & 3) as u8).collect();
        let primary = 73;
        let bwt = mb_bwt_init_from_raw(1, &raw, raw.len() as u64, primary);

        for k in 0..=raw.len() as u64 + 1 {
            let mut cnt = [0u64; 4];
            mb_bwt_rank1a(&bwt, k, &mut cnt);
            for c in 0..4u8 {
                let mut brute = 0u64;
                for j in 0..k {
                    if j == primary {
                        continue;
                    }
                    let z = j - (j > primary) as u64;
                    if raw[z as usize] == c {
                        brute += 1;
                    }
                }
                assert_eq!(mb_bwt_rank11(&bwt, k, c), brute, "rank11 k={k} c={c}");
                assert_eq!(cnt[c as usize], brute, "rank1a k={k} c={c}");
            }
        }
    }

    #[test]
    fn rank2_matches_two_rank1_queries() {
        let raw: Vec<u8> = (0..513).map(|i| ((i * 5 + i / 11) & 3) as u8).collect();
        let bwt = mb_bwt_init_from_raw(1, &raw, raw.len() as u64, 211);

        for &(k, l) in &[
            (0, 0),
            (0, 1),
            (1, 2),
            (31, 32),
            (63, 65),
            (128, 191),
            (212, 320),
            (500, 514),
        ] {
            let mut cntk = [0u64; 4];
            let mut cntl = [0u64; 4];
            let mut expect_k = [0u64; 4];
            let mut expect_l = [0u64; 4];
            mb_bwt_rank2a(&bwt, k, l, &mut cntk, &mut cntl);
            mb_bwt_rank1a(&bwt, k, &mut expect_k);
            mb_bwt_rank1a(&bwt, l, &mut expect_l);
            assert_eq!(cntk, expect_k, "k={k} l={l}");
            assert_eq!(cntl, expect_l, "k={k} l={l}");
        }
    }

    #[cfg(unix)]
    #[test]
    fn save_ignores_post_open_write_errors_like_original() {
        if !std::path::Path::new("/dev/full").exists() {
            return;
        }
        let mut bwt = mb_bwt_init_from_raw(1, &[0, 1, 2, 3], 4, 2);
        bwt.sa_bit = u32::MAX;
        assert_eq!(mb_bwt_save("/dev/full", &bwt), 0);
    }

    #[test]
    fn init_from_raw_empty_input_returns_zero_filled_bwt() {
        let bwt = mb_bwt_init_from_raw(1, &[], 0, 0);
        assert_eq!(bwt.seq_len, 0);
        assert_eq!(bwt.data_len, mb_bwt_data_len(0));
        assert_eq!(bwt.data, vec![0; bwt.data_len as usize]);
        assert_eq!(bwt.l2, [0; 5]);
    }

    #[test]
    fn load_raw_matches_init_from_word_packed_raw() {
        let raw_bases: Vec<u8> = (0..201).map(|i| ((i * 13 + 1) & 3) as u8).collect();
        let mut packed = vec![0u8; ((raw_bases.len() + 15) / 16) * 4];
        for (i, &base) in raw_bases.iter().enumerate() {
            let word = i >> 4;
            let shift = ((!i & 0xf) << 1) as u32;
            let mut val = u32::from_le_bytes(packed[word * 4..word * 4 + 4].try_into().unwrap());
            val |= (base as u32) << shift;
            packed[word * 4..word * 4 + 4].copy_from_slice(&val.to_le_bytes());
        }
        let expected = mb_bwt_init_from_raw(0, &packed, raw_bases.len() as u64, 37);
        let path = std::env::temp_dir().join("minibwa-rs-raw.mbw");
        {
            let mut fp = File::create(&path).expect("create raw GenBwt");
            fp.write_all(&37u64.to_le_bytes()).unwrap();
            fp.write_all(&expected.l2[1].to_le_bytes()).unwrap();
            fp.write_all(&expected.l2[2].to_le_bytes()).unwrap();
            fp.write_all(&expected.l2[3].to_le_bytes()).unwrap();
            fp.write_all(&expected.l2[4].to_le_bytes()).unwrap();
            fp.write_all(&packed).unwrap();
        }
        let loaded = mb_bwt_load_raw(&path).expect("load raw GenBwt");
        let _ = std::fs::remove_file(&path);
        assert_eq!(loaded, expected);
    }
}
