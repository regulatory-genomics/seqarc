#![allow(unused_variables)]

pub(crate) const KSW_NEG_INF: i32 = -0x40000000;
pub(crate) const KSW_EZ_SCORE_ONLY: i32 = 0x01;
pub(crate) const KSW_EZ_RIGHT: i32 = 0x02;
pub(crate) const KSW_EZ_GENERIC_SC: i32 = 0x04;
pub(crate) const KSW_EZ_APPROX_MAX: i32 = 0x08;
pub(crate) const KSW_EZ_APPROX_DROP: i32 = 0x10;
pub(crate) const KSW_EZ_EXTZ_ONLY: i32 = 0x40;
pub(crate) const KSW_EZ_REV_CIGAR: i32 = 0x80;
pub(crate) const KSW_LL_STOP: i32 = 0x20000;
pub(crate) const KSW_LL_SUBO: i32 = 0x40000;
pub(crate) const KSW_CIGAR_MATCH: u32 = 0;
pub(crate) const KSW_CIGAR_INS: u32 = 1;
pub(crate) const KSW_CIGAR_DEL: u32 = 2;
pub(crate) const KSW_CIGAR_N_SKIP: u32 = 3;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct KswExtz {
    pub max: u32,
    pub zdropped: u32,
    pub max_q: i32,
    pub max_t: i32,
    pub mqe: i32,
    pub mqe_t: i32,
    pub mte: i32,
    pub mte_q: i32,
    pub score: i32,
    pub m_cigar: i32,
    pub n_cigar: i32,
    pub reach_end: i32,
    pub cigar: Vec<u32>,
}

impl Default for KswExtz {
    fn default() -> Self {
        Self {
            max: 0,
            zdropped: 0,
            max_q: -1,
            max_t: -1,
            mqe: KSW_NEG_INF,
            mqe_t: -1,
            mte: KSW_NEG_INF,
            mte_q: -1,
            score: KSW_NEG_INF,
            m_cigar: 0,
            n_cigar: 0,
            reach_end: 0,
            cigar: Vec::new(),
        }
    }
}

/// Per-thread scratch buffers for the SSE gapped-extension kernels. Reusing
/// one instance per worker thread avoids repeated large allocations.
#[derive(Default)]
pub(crate) struct Ksw2Scratch {
    pub u: Vec<u8>,
    pub v: Vec<u8>,
    pub x: Vec<u8>,
    pub y: Vec<u8>,
    pub x2: Vec<u8>,
    pub y2: Vec<u8>,
    pub s: Vec<u8>,
    pub sf_qr: Vec<u8>,
    pub sf: Vec<u8>,
    pub qr: Vec<u8>,
    pub h: Vec<i32>,
    pub p: Vec<u8>,
    pub off: Vec<i32>,
    pub off_end: Vec<i32>,
}


/// Original C static function `ksw_push_cigar` from `minibwa/ksw2.h:124`.
pub(crate) fn ksw_push_cigar(
    n_cigar: &mut i32,
    m_cigar: &mut i32,
    cigar: &mut Vec<u32>,
    op: u32,
    len: i32,
) {
    if *n_cigar == 0 || op != (cigar[*n_cigar as usize - 1] & 0xf) {
        if *n_cigar == *m_cigar {
            *m_cigar = if *m_cigar != 0 { *m_cigar << 1 } else { 4 };
            cigar.reserve((*m_cigar as usize).saturating_sub(cigar.capacity()));
        }
        if *n_cigar as usize == cigar.len() {
            cigar.push((len as u32) << 4 | op);
        } else {
            cigar[*n_cigar as usize] = (len as u32) << 4 | op;
        }
        *n_cigar += 1;
    } else {
        cigar[*n_cigar as usize - 1] += (len as u32) << 4;
    }
}

/// Original C static function `ksw_backtrack` from `minibwa/ksw2.h:140`.
pub(crate) fn ksw_backtrack(
    is_rot: i32,
    is_rev: i32,
    min_intron_len: i32,
    p: &[u8],
    off: &[i32],
    off_end: Option<&[i32]>,
    n_col: i32,
    i0: i32,
    j0: i32,
    m_cigar_: &mut i32,
    n_cigar_: &mut i32,
    cigar_: &mut Vec<u32>,
) {
    let mut n_cigar = 0;
    let mut m_cigar = *m_cigar_;
    let cigar = cigar_;
    let mut i = i0;
    let mut j = j0;
    let mut state = 0i32;
    while i >= 0 && j >= 0 {
        let mut force_state = -1;
        let tmp = if is_rot != 0 {
            let r = i + j;
            if i < off[r as usize] {
                force_state = 2;
            }
            if let Some(off_end) = off_end {
                if i > off_end[r as usize] {
                    force_state = 1;
                }
            }
            if force_state < 0 {
                p[(r * n_col + i - off[r as usize]) as usize]
            } else {
                0
            }
        } else {
            if j < off[i as usize] {
                force_state = 2;
            }
            if let Some(off_end) = off_end {
                if j > off_end[i as usize] {
                    force_state = 1;
                }
            }
            if force_state < 0 {
                p[(i * n_col + j - off[i as usize]) as usize]
            } else {
                0
            }
        };
        if state == 0 {
            state = (tmp & 7) as i32;
        } else if ((tmp >> (state + 2)) & 1) == 0 {
            state = 0;
        }
        if state == 0 {
            state = (tmp & 7) as i32;
        }
        if force_state >= 0 {
            state = force_state;
        }
        if state == 0 {
            ksw_push_cigar(&mut n_cigar, &mut m_cigar, cigar, KSW_CIGAR_MATCH, 1);
            i -= 1;
            j -= 1;
        } else if state == 1 || (state == 3 && min_intron_len <= 0) {
            ksw_push_cigar(&mut n_cigar, &mut m_cigar, cigar, KSW_CIGAR_DEL, 1);
            i -= 1;
        } else if state == 3 && min_intron_len > 0 {
            ksw_push_cigar(&mut n_cigar, &mut m_cigar, cigar, KSW_CIGAR_N_SKIP, 1);
            i -= 1;
        } else {
            ksw_push_cigar(&mut n_cigar, &mut m_cigar, cigar, KSW_CIGAR_INS, 1);
            j -= 1;
        }
    }
    if i >= 0 {
        let op = if min_intron_len > 0 && i >= min_intron_len {
            KSW_CIGAR_N_SKIP
        } else {
            KSW_CIGAR_DEL
        };
        ksw_push_cigar(&mut n_cigar, &mut m_cigar, cigar, op, i + 1);
    }
    if j >= 0 {
        ksw_push_cigar(&mut n_cigar, &mut m_cigar, cigar, KSW_CIGAR_INS, j + 1);
    }
    if is_rev == 0 {
        cigar[..n_cigar as usize].reverse();
    }
    cigar.truncate(n_cigar as usize);
    *m_cigar_ = m_cigar;
    *n_cigar_ = n_cigar;
}

/// Original C static function `ksw_reset_extz` from `minibwa/ksw2.h:174`.
pub(crate) fn ksw_reset_extz(ez: &mut KswExtz) {
    ez.max_q = -1;
    ez.max_t = -1;
    ez.mqe_t = -1;
    ez.mte_q = -1;
    ez.max = 0;
    ez.score = KSW_NEG_INF;
    ez.mqe = KSW_NEG_INF;
    ez.mte = KSW_NEG_INF;
    ez.n_cigar = 0;
    ez.zdropped = 0;
    ez.reach_end = 0;
}

/// Original C static function `ksw_apply_zdrop` from `minibwa/ksw2.h:181`.
pub(crate) fn ksw_apply_zdrop(
    ez: &mut KswExtz,
    is_rot: i32,
    h: i32,
    a: i32,
    b: i32,
    zdrop: i32,
    e: i8,
) -> i32 {
    let (r, t) = if is_rot != 0 { (a, b) } else { (a + b, a) };
    if h > ez.max as i32 {
        ez.max = h as u32;
        ez.max_t = t;
        ez.max_q = r - t;
    } else if t >= ez.max_t && r - t >= ez.max_q {
        let tl = t - ez.max_t;
        let ql = (r - t) - ez.max_q;
        let l = if tl > ql { tl - ql } else { ql - tl };
        if zdrop >= 0 && ez.max as i32 - h > zdrop + l * e as i32 {
            ez.zdropped = 1;
            return 1;
        }
    }
    0
}

/// Original C static function `ksw_gen_nt4_mat` from `minibwa/ksw2.h:199`.
pub(crate) fn ksw_gen_nt4_mat(mat: &mut [i8; 25], a: i8, b: i8, b_ts: i8, b_ambi: i8, mt: i32) {
    let a = if a < 0 { a.wrapping_neg() } else { a };
    let b = if b > 0 { -b } else { b };
    let b_ambi = if b_ambi > 0 { -b_ambi } else { b_ambi };
    for i in 0..4usize {
        for j in 0..4usize {
            mat[i * 5 + j] = if i == j { a } else { b };
        }
        mat[i * 5 + 4] = b_ambi;
    }
    for j in 0..5usize {
        mat[4 * 5 + j] = b_ambi;
    }
    if b_ts != 0 && b_ts != b {
        let b_ts = if b_ts > 0 { -b_ts } else { b_ts };
        mat[2] = b_ts;
        mat[8] = b_ts;
        mat[10] = b_ts;
        mat[16] = b_ts;
    }
    if mt == 1 {
        mat[1 * 5 + 3] = a;
    } else if mt == 2 {
        mat[2 * 5] = a;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cigar_push_merges_adjacent_ops() {
        let mut n = 0;
        let mut m = 0;
        let mut cigar = Vec::new();
        ksw_push_cigar(&mut n, &mut m, &mut cigar, KSW_CIGAR_MATCH, 3);
        ksw_push_cigar(&mut n, &mut m, &mut cigar, KSW_CIGAR_MATCH, 2);
        ksw_push_cigar(&mut n, &mut m, &mut cigar, KSW_CIGAR_INS, 1);
        assert_eq!(n, 2);
        assert_eq!(
            cigar[..n as usize],
            [5 << 4 | KSW_CIGAR_MATCH, 1 << 4 | KSW_CIGAR_INS]
        );
    }

    #[test]
    fn backtrack_builds_and_reverses_cigar() {
        let p = [0u8; 4];
        let off = [0i32; 2];
        let off_end = [1i32; 2];
        let mut m = 0;
        let mut n = 0;
        let mut cigar = Vec::new();
        ksw_backtrack(0,
            0,
            0,
            &p,
            &off,
            Some(&off_end),
            2,
            1,
            1,
            &mut m,
            &mut n,
            &mut cigar,
        );
        assert_eq!(n, 1);
        assert_eq!(cigar[..n as usize], [2 << 4 | KSW_CIGAR_MATCH]);

        let p = [1u8, 0, 0, 0];
        ksw_backtrack(0,
            0,
            5,
            &p,
            &off,
            Some(&off_end),
            2,
            1,
            1,
            &mut m,
            &mut n,
            &mut cigar,
        );
        assert_eq!(
            cigar[..n as usize],
            [
                1 << 4 | KSW_CIGAR_INS,
                1 << 4 | KSW_CIGAR_DEL,
                1 << 4 | KSW_CIGAR_MATCH
            ]
        );
    }
}
