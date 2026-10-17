#![allow(unused_variables)]

use crate::minibwa::align::{mb_append_cigar, mb_update_extra};
use crate::minibwa::kommon::KOM_NT4_TABLE;
use crate::minibwa::ksw2::{
    KswExtz, ksw_gen_nt4_mat, KSW_EZ_EXTZ_ONLY, KSW_EZ_GENERIC_SC, KSW_EZ_REV_CIGAR,
    KSW_EZ_RIGHT, KSW_LL_SUBO, Ksw2Scratch,
};
use crate::minibwa::ksw2_extz2_sse::ksw_extz2_sse;
use crate::minibwa::ksw2_ll_sse::{ksw_ll_i16_core, ksw_ll_qinit, ksw_ll_u8_core};
use crate::minibwa::l2bit::{l2b_getseq, l2b_meth_rev, L2bMeth, L2b};
use crate::minibwa::lchain::{Mb128, radix_sort_mb128x};
use crate::minibwa::map_algo::{
    mb_hit_sort, mb_set_mapq, mb_set_parent, mb_set_sam_pri, mb_sync_high_cov, MB_PARENT_UNSET,
};
use crate::minibwa::mbpriv::{mb_is_sr_mode, mb_seq_rev, KOM_DBG_FLAG, MB_DBG_ALN_PE};
use crate::minibwa::options::{MapOpts, MB_F_METH, MB_F_PRIMARY5};
use std::sync::atomic::Ordering;

/// Per-hit dynamic data: DP scores plus the packed CIGAR / tag words.
/// `cigar` holds the `n_cigar` CIGAR words, optionally followed by the packed
/// tag words (the CS tag) when `cs()` is set.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Extra {
    pub dp_score: i32,
    pub dp_max0: i32,
    pub dp_max: i32,
    pub dp_max2: i32,
    pub n_ambi_cs: u32,
    pub n_cigar: i32,
    pub cigar: Vec<u32>,
}

impl Extra {
    pub(crate) const CS_FLAG: u32 = 1 << 31;
    const N_AMBI_MASK: u32 = !Self::CS_FLAG;

    #[inline(always)]
    pub(crate) fn n_ambi(&self) -> u32 {
        self.n_ambi_cs & Self::N_AMBI_MASK
    }

    #[inline(always)]
    pub(crate) fn set_n_ambi(&mut self, value: u32) {
        self.n_ambi_cs = (self.n_ambi_cs & Self::CS_FLAG) | (value & Self::N_AMBI_MASK);
    }

    #[inline(always)]
    pub(crate) fn add_n_ambi(&mut self, value: u32) {
        self.set_n_ambi(self.n_ambi().saturating_add(value));
    }

    #[inline(always)]
    pub(crate) fn set_cs(&mut self, value: u32) {
        if value != 0 {
            self.n_ambi_cs |= Self::CS_FLAG;
        } else {
            self.n_ambi_cs &= !Self::CS_FLAG;
        }
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
        self.set_cs(0);
        self.cigar = values;
    }

    pub(crate) fn truncate_cigar(&mut self, len: usize) {
        self.n_cigar = self.n_cigar.min(len as i32);
        self.set_cs(0);
    }

    /// Replaces the words after `n_cigar` with `values` and raises the CS flag
    /// when non-empty.
    pub(crate) fn set_tag_words_from_slice(&mut self, values: &[u32]) {
        let n = self.n_cigar.max(0) as usize;
        self.cigar.truncate(n);
        self.cigar.extend_from_slice(values);
        self.set_cs(!values.is_empty() as u32);
    }

    /// Appends CIGAR words after the current `n_cigar` prefix (dropping any
    /// packed tag words) and raises `n_cigar` accordingly.
    pub(crate) fn extend_cigar_from_slice(&mut self, values: &[u32]) {
        if values.is_empty() {
            return;
        }
        let n = self.n_cigar.max(0) as usize;
        self.cigar.truncate(n);
        self.cigar.extend_from_slice(values);
        self.n_cigar = self.cigar.len() as i32;
        self.set_cs(0);
    }

    pub(crate) fn remove_cigar(&mut self, index: usize) -> u32 {
        let len = self.n_cigar.max(0) as usize;
        assert!(index < len);
        self.n_cigar -= 1;
        self.set_cs(0);
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
    const RESCUED: u32 = 1 << 8;
    const FRAC_HIGH_SHIFT: u32 = 9;
    const FRAC_HIGH_MASK: u32 = 0xff << Self::FRAC_HIGH_SHIFT;
    const SEED_RATIO_SHIFT: u32 = 17;
    const SEED_RATIO_MASK: u32 = 0xff << Self::SEED_RATIO_SHIFT;

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
    pub(crate) fn set_proper_pair(&mut self, value: u8) {
        self.set_bit(Self::PROPER_PAIR, value);
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
    pub(crate) fn rescued(&self) -> u8 {
        self.bit(Self::RESCUED)
    }

    #[inline(always)]
    pub(crate) fn set_rescued(&mut self, value: u8) {
        self.set_bit(Self::RESCUED, value);
    }

    #[inline(always)]
    pub(crate) fn frac_high(&self) -> u8 {
        ((self.flags & Self::FRAC_HIGH_MASK) >> Self::FRAC_HIGH_SHIFT) as u8
    }

    #[inline(always)]
    pub(crate) fn set_frac_high(&mut self, value: u8) {
        self.flags =
            (self.flags & !Self::FRAC_HIGH_MASK) | ((value as u32) << Self::FRAC_HIGH_SHIFT);
    }

    #[inline(always)]
    pub(crate) fn seed_ratio(&self) -> u8 {
        ((self.flags & Self::SEED_RATIO_MASK) >> Self::SEED_RATIO_SHIFT) as u8
    }

    #[inline(always)]
    pub(crate) fn set_seed_ratio(&mut self, value: u8) {
        self.flags =
            (self.flags & !Self::SEED_RATIO_MASK) | ((value as u32) << Self::SEED_RATIO_SHIFT);
    }

}


#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct Pestat {
    pub lo: i32,
    pub hi: i32,
    pub failed: i32,
    pub avg: f64,
    pub std: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct PairAux {
    pub score: i32,
    pub sub_sc: i32,
    pub n_sub: i32,
    pub n_pp: i32,
    pub i: [i32; 2],
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct HitVec {
    pub n: i32,
    pub m: i32,
    pub a: Vec<Hit>,
}

/// Original C static function `mb_insert_dir` from `minibwa/pe.c:10`.
pub(crate) fn mb_insert_dir(h0: &Hit, h1: &Hit, dist: &mut i64) -> i32 {
    let p0 = if h0.rev() != 0 { h0.te } else { h0.ts };
    let p1 = if h1.rev() != 0 { h1.te } else { h1.ts };
    *dist = if p0 > p1 { p0 - p1 } else { p1 - p0 };
    (((h0.rev() as i32) << 1 | h1.rev() as i32) ^ if p0 < p1 { 0 } else { 3 }) as i32
}

/// Original C static function `mb_pair_score` from `minibwa/pe.c:19`.
pub(crate) fn mb_pair_score(h0: &Hit, h1: &Hit, pes: &[Pestat; 4], match_sc: i32) -> f64 {
    const MB_SQRT1_2: f64 = 0.707106781186547524401;
    let mut dist = 0i64;
    let dir = mb_insert_dir(h0, h1, &mut dist) as usize;
    if pes[dir].failed == 0 && dist >= pes[dir].lo as i64 && dist <= pes[dir].hi as i64 {
        let p0 = h0.p.as_ref().expect("mb_pair_score requires h0->p");
        let p1 = h1.p.as_ref().expect("mb_pair_score requires h1->p");
        let ns = (dist as f64 - pes[dir].avg) / pes[dir].std;
        return p0.dp_max as f64
            + p1.dp_max as f64
            + 0.721 * libm::erfc(ns.abs() * MB_SQRT1_2).mul_add(2.0, 0.0).ln() * match_sc as f64;
    }
    -1.0
}

/// Original C static function `mb_hit_sum_score` from `minibwa/pe.c:41`.
pub(crate) fn mb_hit_sum_score(n_hit: i32, hit: &[Hit]) -> i32 {
    let mut n_pri = 0;
    let mut sc = 0;
    for h in hit.iter().take(n_hit as usize) {
        if h.id == h.parent && h.p.is_some() {
            n_pri += 1;
            sc = h.p.as_ref().unwrap().dp_max;
        }
    }
    if n_pri == 0 {
        return 0;
    }
    if n_pri == 1 {
        return sc;
    }
    let mut a = Vec::with_capacity(n_pri as usize);
    for (i, h) in hit.iter().take(n_hit as usize).enumerate() {
        if h.id == h.parent && h.p.is_some() {
            a.push(((h.qs as u64) << 32) | i as u64);
        }
    }
    a.sort_unstable();
    sc = 0;
    let mut qe = 0;
    for x in a {
        let h = &hit[x as u32 as usize];
        if h.qe <= qe {
            continue;
        }
        let dp_max = h.p.as_ref().unwrap().dp_max;
        if h.qs < qe {
            sc += ((h.qe - qe) as f64 / (h.qe - h.qs) as f64 * dp_max as f64 + 0.499) as i32;
        } else {
            sc += dp_max;
        }
        qe = h.qe;
    }
    sc
}

/// Original C static function `mb_pair_hits` from `minibwa/pe.c:140`.
pub(crate) fn mb_pair_hits(
    opt: &MapOpts,
    l2b: &L2b,
    n_hit: [i32; 2],
    hit: &mut [Vec<Hit>; 2],
    pes: &[Pestat; 4],
    ret: &mut PairAux,
) {
    *ret = PairAux {
        score: -1,
        sub_sc: -1,
        i: [-1, -1],
        ..Default::default()
    };
    if n_hit[0] == 0 || n_hit[1] == 0 {
        return;
    }
    let mut pa = Vec::<Mb128>::with_capacity((n_hit[0] + n_hit[1]) as usize);
    for r in 0..2usize {
        for i in 0..n_hit[r] as usize {
            let h = &mut hit[r][i];
            h.set_proper_pair(0);
            let p = l2b.ctg[h.tid as usize].off + if h.rev() != 0 { h.te } else { h.ts } as u64;
            pa.push(Mb128 {
                x: p,
                y: ((i as u64) << 2) | ((h.rev() as u64) << 1) | r as u64,
            });
        }
    }
    radix_sort_mb128x(&mut pa);
    let mut y = [-1i32; 4];
    let mut pp = Vec::<(u64, u64)>::new();
    for i in 0..pa.len() {
        let pi = pa[i];
        let pi_read = (pi.y & 1) as usize;
        let pi_idx = (pi.y >> 2) as usize;
        let hi_tid = hit[pi_read][pi_idx].tid;
        for r in 0..2usize {
            let dir = (r << 1) | ((pi.y >> 1) as usize & 1);
            if pes[dir].failed != 0 {
                continue;
            }
            let which = (r << 1) | ((pi.y as usize & 1) ^ 1);
            if y[which] < 0 {
                continue;
            }
            let mut k = y[which] as isize;
            while k >= 0 {
                let pk = pa[k as usize];
                if (pk.y & 3) as usize != which {
                    k -= 1;
                    continue;
                }
                let pk_read = (pk.y & 1) as usize;
                let pk_idx = (pk.y >> 2) as usize;
                let hk_tid = hit[pk_read][pk_idx].tid;
                if hi_tid != hk_tid {
                    break;
                }
                let dist = pi.x as i64 - pk.x as i64;
                if dist > pes[dir].hi as i64 {
                    break;
                }
                if dist >= pes[dir].lo as i64 {
                    hit[pk_read][pk_idx].set_proper_pair(1);
                    hit[pi_read][pi_idx].set_proper_pair(1);
                    let mut s = if pk_read == pi_read {
                        -1.0
                    } else {
                        let (left, right) = hit.split_at(1);
                        if pk_read == 0 {
                            mb_pair_score(&left[0][pk_idx], &right[0][pi_idx], pes, opt.a)
                        } else {
                            mb_pair_score(&left[0][pi_idx], &right[0][pk_idx], pes, opt.a)
                        }
                    };
                    if s < 0.0 {
                        s = 0.0;
                    }
                    let yv = if (pk.y & 1) == 0 {
                        ((pk_idx as u64) << 32) | pi_idx as u64
                    } else {
                        ((pi_idx as u64) << 32) | pk_idx as u64
                    };
                    let hash = hit[pk_read][pk_idx].hash ^ hit[pi_read][pi_idx].hash;
                    pp.push((((s + 0.499) as u64) << 32 | hash as u64, yv));
                }
                k -= 1;
            }
        }
        y[(pi.y & 3) as usize] = i as i32;
    }
    ret.n_pp = pp.len() as i32;
    if !pp.is_empty() {
        let mut max = 0u64;
        let mut max2 = 0u64;
        let tmp = (opt.a + opt.b).max(opt.q + opt.e);
        for &(x, yv) in &pp {
            if x > max {
                max2 = max;
                max = x;
                ret.i[0] = (yv >> 32) as i32;
                ret.i[1] = yv as u32 as i32;
            } else if x > max2 {
                max2 = x;
            }
        }
        ret.score = (max >> 32) as i32;
        ret.sub_sc = (max2 >> 32) as i32;
        for &(x, _) in &pp {
            if (x >> 32) as i32 >= ret.score - tmp {
                ret.n_sub += 1;
            }
        }
    }
}

/// A linear algorithm to find ungapped alignment.
///
/// Original C static function `mb_ungap` from `minibwa/pe.c:220`.
pub(crate) fn mb_ungap(
    qlen: i32,
    qseq: &[u8],
    tlen: i32,
    tseq: &[u8],
    kmer: i32,
    mt: L2bMeth,
    max_i: &mut i32,
    n_good: &mut i32,
    n_kmer: &mut i32,
) -> i32 {
    const C2T: [u8; 4] = [0, 3, 2, 3];
    const G2A: [u8; 4] = [0, 1, 0, 3];
    *max_i = -1;
    *n_good = 0;
    *n_kmer = 0;
    if qlen >= u16::MAX as i32 || kmer <= 0 {
        return 0;
    }
    let cap = 1usize << (kmer * 2);
    let mask = cap - 1;
    let mut a = vec![0i32; tlen.max(0) as usize];
    let mut h = vec![0u16; cap];
    let mut l = 0i32;
    let mut x = 0usize;
    for i in 0..qlen as usize {
        if qseq[i] < 4 {
            let c = if mt == L2bMeth::C2t {
                C2T[qseq[i] as usize]
            } else if mt == L2bMeth::G2a {
                G2A[qseq[i] as usize]
            } else {
                qseq[i]
            };
            x = ((x << 2) | c as usize) & mask;
            l += 1;
            if l >= kmer {
                if h[x] == 0 {
                    *n_kmer += 1;
                }
                h[x] = i as u16;
            }
        } else {
            x = 0;
            l = 0;
        }
    }
    l = 0;
    x = 0;
    for i in 0..tlen as usize {
        if tseq[i] < 4 {
            let c = if mt == L2bMeth::C2t {
                C2T[tseq[i] as usize]
            } else if mt == L2bMeth::G2a {
                G2A[tseq[i] as usize]
            } else {
                tseq[i]
            };
            x = ((x << 2) | c as usize) & mask;
            l += 1;
            if l >= kmer && h[x] > 0 && i >= h[x] as usize {
                a[i - h[x] as usize] += 1;
            }
        } else {
            x = 0;
            l = 0;
        }
    }
    let mut max = 0;
    for (i, &v) in a.iter().enumerate() {
        if max < v {
            max = v;
            *max_i = i as i32;
        }
    }
    for &v in &a {
        if v > max >> 1 {
            *n_good += 1;
        }
    }
    max
}

/// Original C static function `mb_matesw_align` from `minibwa/pe.c:245`.
pub(crate) fn mb_matesw_align(
    scratch: &mut Ksw2Scratch,
    opt: &MapOpts,
    qlen: i32,
    qseq: &mut [u8],
    tlen: i32,
    tseq: &mut [u8],
    h: &mut Hit,
    min_sc: i32,
    mt: L2bMeth,
    ez: &mut KswExtz,
) {
    *h = Hit::default();
    let max_sc = qlen.min(tlen);
    let b_mm = (opt.b + opt.a - 1) / opt.a;
    let b_ts = (opt.b_ts + opt.a - 1) / opt.a;
    let b_ambi = (opt.b_ambi + opt.a - 1) / opt.a;
    let gapo = (opt.q + opt.a - 1) / opt.a;
    let gape = (opt.e + opt.a - 1) / opt.a;
    if max_sc >= 32767 || qlen <= 0 || tlen <= 0 {
        return;
    }
    let mut mat = [0i8; 25];
    ksw_gen_nt4_mat(&mut mat, 1, b_mm as i8, b_ts as i8, b_ambi as i8, mt as i32);
    let sz = if max_sc < 255 - b_mm { 1 } else { 2 };
    let xtra = KSW_LL_SUBO | opt.min_len;
    let qp = ksw_ll_qinit(sz, qlen, qseq, 5, &mat);
    let rst = if sz == 1 {
        ksw_ll_u8_core(&qp, tlen, tseq, gapo, gape, xtra)
    } else {
        ksw_ll_i16_core(&qp, tlen, tseq, gapo, gape, xtra)
    };
    if (KOM_DBG_FLAG.load(Ordering::Relaxed) & MB_DBG_ALN_PE) != 0 {
        eprintln!(
            "===> qlen={}; tlen={}; score={}; qe={}; te={} <===",
            qlen,
            tlen,
            rst.score,
            rst.qe + 1,
            rst.te + 1
        );
        let alphabet = b"ACGTN";
        eprintln!(
            "{}",
            qseq.iter()
                .take(qlen.max(0) as usize)
                .map(|&c| alphabet[c.min(4) as usize] as char)
                .collect::<String>()
        );
        eprintln!(
            "{}",
            tseq.iter()
                .take(tlen.max(0) as usize)
                .map(|&c| alphabet[c.min(4) as usize] as char)
                .collect::<String>()
        );
    }
    if rst.score >= opt.min_dp_max && rst.score >= min_sc {
        let te = rst.te + 1;
        let qe = rst.qe + 1;
        if te <= 0 || qe <= 0 {
            return;
        }
        mb_seq_rev(qe as u32, qseq);
        mb_seq_rev(te as u32, tseq);
        ksw_gen_nt4_mat(
            &mut mat,
            opt.a as i8,
            opt.b as i8,
            opt.b_ts as i8,
            opt.b_ambi as i8,
            mt as i32,
        );
        let mut ksw_flag = KSW_EZ_EXTZ_ONLY | KSW_EZ_RIGHT | KSW_EZ_REV_CIGAR;
        if mt != L2bMeth::None {
            ksw_flag |= KSW_EZ_GENERIC_SC;
        }
        ksw_extz2_sse(scratch,
            qe,
            qseq,
            te,
            tseq,
            5,
            &mat,
            opt.q as i8,
            opt.e as i8,
            opt.bw,
            opt.zdrop,
            opt.end_bonus,
            ksw_flag,
            ez,
        );
        mb_seq_rev(qe as u32, qseq);
        mb_seq_rev(te as u32, tseq);
        if ez.n_cigar > 0 && ez.max as i32 >= opt.min_dp_max * opt.a {
            h.p = Some(Box::new(Extra::default()));
            let cigar = ez.cigar.clone();
            mb_append_cigar(h, ez.n_cigar as u32, &cigar);
            h.set_rescued(1);
            h.qe = qe;
            h.te = te as i64;
            h.ts = (te
                - if ez.reach_end != 0 {
                    ez.mqe_t + 1
                } else {
                    ez.max_t + 1
                }) as i64;
            h.qs = qe - if ez.reach_end != 0 { qe } else { ez.max_q + 1 };
            if let Some(p) = &mut h.p {
                p.dp_max = ez.max as i32;
                p.dp_score = ez.max as i32;
                p.dp_max2 = ((ez.max as f64 / rst.score as f64) * rst.score2 as f64 + 0.499) as i32;
                if p.dp_max2 < 0 {
                    p.dp_max2 = 0;
                }
            }
            h.score = rst.score;
            h.score0 = rst.score;
            h.subsc = rst.score2;
            h.cnt = 0;
            h.as_ = -1;
            h.parent = MB_PARENT_UNSET;
            if (KOM_DBG_FLAG.load(Ordering::Relaxed) & MB_DBG_ALN_PE) != 0 {
                let cigar = ez
                    .cigar
                    .iter()
                    .take(ez.n_cigar as usize)
                    .map(|&c| {
                        format!(
                            "{}{}",
                            c >> 4,
                            crate::minibwa::align::MB_CIGAR_STR.as_bytes()[(c & 0xf) as usize] as char
                        )
                    })
                    .collect::<String>();
                eprintln!(
                    "max={}; ts={}; qs={}; reach_end={}; cigar={}",
                    ez.max, h.ts, h.qs, ez.reach_end, cigar
                );
            }
            mb_update_extra(h,
                &qseq[h.qs.max(0) as usize..qe as usize],
                &tseq[h.ts.max(0) as usize..te as usize],
                &mat,
                opt.q as i8,
                opt.e as i8,
                opt.flag,
                0,
            );
        }
    }
}

/// Original C static function `mb_matesw_core` from `minibwa/pe.c:314`.
pub(crate) fn mb_matesw_core(
    scratch: &mut Ksw2Scratch,
    opt: &MapOpts,
    l2b: &L2b,
    pes: &[Pestat; 4],
    h0: &Hit,
    r0: i32,
    len: i32,
    seq: &mut [Vec<u8>; 2],
    mt0: L2bMeth,
    h1: &mut HitVec,
    min_sc: i32,
    ez: &mut KswExtz,
) -> Option<f64> {
    let mut skip = [false; 4];
    for dir in 0..4usize {
        skip[dir] = pes[dir].failed != 0;
    }
    if skip.iter().all(|&x| x) {
        return None;
    }
    let pos5 = if h0.rev() != 0 { h0.te } else { h0.ts };
    let mut ret = None;
    for dir in 0..4usize {
        if skip[dir] {
            continue;
        }
        let is_rev = (((dir >> 1) != (dir & 1)) as u8 ^ h0.rev()) != 0;
        let is_larger = if (dir >> 1) != (dir & 1) {
            (((dir >> 1) as i32) ^ is_rev as i32) != 0
        } else {
            (((dir >> 1) as i32) ^ r0 ^ h0.rev() as i32) != 0
        };
        let mut ts = (if is_larger {
            pos5 + pes[dir].lo as i64
        } else {
            pos5 - pes[dir].hi as i64
        }) - if !is_rev { 0 } else { len as i64 };
        let mut te = (if is_larger {
            pos5 + pes[dir].hi as i64
        } else {
            pos5 - pes[dir].lo as i64
        }) + if !is_rev { len as i64 } else { 0 };
        if ts < 0 {
            ts = 0;
        }
        let ctg_len = l2b.ctg[h0.tid as usize].len as i64;
        if te > ctg_len {
            te = ctg_len;
        }
        if te - ts <= len as i64 {
            continue;
        }
        let mut ts2 = ts;
        let mut te2 = te;
        let mut refseq = vec![0u8; (te - ts) as usize];
        let mt = if is_rev { l2b_meth_rev(mt0) } else { mt0 };
        l2b_getseq(l2b, h0.tid, ts, te, &mut refseq);
        let seq_idx = is_rev as usize;
        let mut max_i = 0;
        let mut n_good = 0;
        let mut n_kmer = 0;
        let max_ug = mb_ungap(len,
            &seq[seq_idx],
            (te - ts) as i32,
            &refseq,
            7,
            mt,
            &mut max_i,
            &mut n_good,
            &mut n_kmer,
        );
        if max_ug >= 10 && max_ug >= len >> 1 && n_good == 1 {
            ts2 = ts + max_i as i64 - (len / 2) as i64;
            te2 = ts2 + (len * 2) as i64;
            if ts2 < ts {
                ts2 = ts;
            }
            if te2 > te {
                te2 = te;
            }
        }
        if max_ug >= 10 || (max_ug as f64) >= n_kmer as f64 * 0.33 {
            let mut ht = Hit {
                p: None,
                ..Default::default()
            };
            let offset = (ts2 - ts) as usize;
            mb_matesw_align(                scratch,
                opt,
                len,
                &mut seq[seq_idx],
                (te2 - ts2) as i32,
                &mut refseq[offset..offset + (te2 - ts2) as usize],
                &mut ht,
                min_sc,
                mt,
                ez,
            );
            if ht.p.is_some() {
                ht.tid = h0.tid;
                ht.ts += ts2;
                ht.te += ts2;
                ht.set_rev(is_rev as u8);
                if is_rev {
                    let qt = ht.qs;
                    ht.qs = len - ht.qe;
                    ht.qe = len - qt;
                }
                let score = mb_pair_score(h0, &ht, pes, opt.a);
                h1.a.push(ht);
                h1.n = h1.a.len() as i32;
                h1.m = h1.n;
                ret = Some(score);
            }
        }
    }
    ret
}

/// Original C static function `mb_matesw` from `minibwa/pe.c:375`.
pub(crate) fn mb_matesw(
    scratch: &mut Ksw2Scratch,
    opt: &MapOpts,
    l2b: &L2b,
    n_hit: &mut [i32; 2],
    hit: &mut [Vec<Hit>; 2],
    pes: &[Pestat; 4],
    paux0: &PairAux,
    qlen: [i32; 2],
    qseq: [&str; 2],
    is_meth: i32,
) -> i32 {
    if opt.max_rescue == 0 {
        return 0;
    }
    let mut n_res = 0i32;
    for r in 0..2usize {
        let mut m = 0;
        if n_hit[r] == 0 {
            continue;
        }
        let best = hit[r][0]
            .p
            .as_ref()
            .expect("mb_matesw requires h0[0]->p")
            .dp_max;
        for i in 0..n_hit[r] as usize {
            if m >= opt.max_rescue {
                break;
            }
            let dp = hit[r][i]
                .p
                .as_ref()
                .expect("mb_matesw requires h0[i]->p")
                .dp_max;
            if hit[r][i].proper_pair() == 0 && dp >= best - opt.pen_unpair * opt.a {
                m += 1;
                n_res += 1;
            }
        }
    }
    if n_res == 0 {
        return 0;
    }
    let mut a = Vec::<(u64, u64)>::with_capacity(n_res as usize);
    for r in 0..2usize {
        let mut m = 0;
        if n_hit[r] == 0 {
            continue;
        }
        let best = hit[r][0]
            .p
            .as_ref()
            .expect("mb_matesw requires h0[0]->p")
            .dp_max;
        for i in 0..n_hit[r] as usize {
            if m >= opt.max_rescue {
                break;
            }
            let dp = hit[r][i]
                .p
                .as_ref()
                .expect("mb_matesw requires h0[i]->p")
                .dp_max;
            if hit[r][i].proper_pair() == 0 && dp >= best - opt.pen_unpair * opt.a {
                a.push((
                    ((dp as u64) << 32) | hit[r][i].hash as u64,
                    ((i as u64) << 1) | r as u64,
                ));
                m += 1;
            }
        }
    }
    a.sort_by_key(|&(x, _)| x);
    let mut qs = [
        [
            Vec::<u8>::with_capacity(qlen[0] as usize),
            Vec::<u8>::with_capacity(qlen[0] as usize),
        ],
        [
            Vec::<u8>::with_capacity(qlen[1] as usize),
            Vec::<u8>::with_capacity(qlen[1] as usize),
        ],
    ];
    for r in 0..2usize {
        let bases = qseq[r].as_bytes();
        let len = qlen[r].max(0) as usize;
        qs[r][0].resize(len, 0);
        qs[r][1].resize(len, 0);
        let fwd = qs[r][0].as_mut_ptr();
        let rev = qs[r][1].as_mut_ptr();
        for (i, &base) in bases.iter().take(len).enumerate() {
            let c = KOM_NT4_TABLE[base as usize];
            unsafe {
                *fwd.add(i) = c;
                *rev.add(len - 1 - i) = if c < 4 { 3 - c } else { 4 };
            }
        }
    }
    let min_sc = [
        mb_hit_sum_score(n_hit[0], &hit[0]) / opt.a - opt.pen_unpair,
        mb_hit_sum_score(n_hit[1], &hit[1]) / opt.a - opt.pen_unpair,
    ];
    let mut ha = [
        HitVec {
            n: n_hit[0],
            m: n_hit[0],
            a: std::mem::take(&mut hit[0]),
        },
        HitVec {
            n: n_hit[1],
            m: n_hit[1],
            a: std::mem::take(&mut hit[1]),
        },
    ];
    let mut max = [paux0.score, paux0.score];
    let mut max2 = [paux0.sub_sc, paux0.sub_sc];
    let mut skip = [false, false];
    let mut ez = KswExtz {
        m_cigar: 16,
        cigar: Vec::with_capacity(16),
        ..Default::default()
    };
    for &(_, yv) in a.iter().rev() {
        let r = (yv & 1) as usize;
        let j = (yv >> 1) as usize;
        if skip[r] {
            continue;
        }
        let mt = if is_meth == 0 {
            L2bMeth::None
        } else if r == 0 {
            L2bMeth::G2a
        } else {
            L2bMeth::C2t
        };
        let sc = if r == 0 {
            let (left, right) = ha.split_at_mut(1);
            mb_matesw_core(
                scratch,
                opt,
                l2b,
                pes,
                &left[0].a[j],
                r as i32,
                qlen[1 - r],
                &mut qs[1 - r],
                mt,
                &mut right[0],
                min_sc[1 - r],
                &mut ez,
            )
        } else {
            let (left, right) = ha.split_at_mut(1);
            mb_matesw_core(
                scratch,
                opt,
                l2b,
                pes,
                &right[0].a[j],
                r as i32,
                qlen[1 - r],
                &mut qs[1 - r],
                mt,
                &mut left[0],
                min_sc[1 - r],
                &mut ez,
            )
        };
        if let Some(sc) = sc {
            let sc = sc as i32;
            if sc > max[r] {
                max2[r] = max[r];
                max[r] = sc;
            } else if sc > max2[r] {
                max2[r] = sc;
            }
            if max[r] == max2[r] {
                skip[r] = true;
            }
        }
    }
    let n_add = (ha[0].n - n_hit[0]) + (ha[1].n - n_hit[1]);
    for r in 0..2usize {
        n_hit[r] = ha[r].n;
        hit[r] = std::mem::take(&mut ha[r].a);
    }
    n_add
}

/// Original C global function `mb_pair` from `minibwa/pe.c:456`.
pub(crate) fn mb_pair(
    scratch: &mut Ksw2Scratch,
    opt: &MapOpts,
    l2b: &L2b,
    n_hit: &mut [i32; 2],
    hit: &mut [Vec<Hit>; 2],
    pes: &[Pestat; 4],
    qlen: [i32; 2],
    qseq: [&str; 2],
) {
    const PE_BONUS: i32 = 4;
    let is_meth = ((opt.flag & MB_F_METH) != 0) as i32;
    let mut paux = PairAux::default();
    if n_hit[0] == 0 && n_hit[1] == 0 {
        return;
    }
    let seed_ratio = [
        if n_hit[0] > 0 {
            hit[0][0].seed_ratio()
        } else {
            255
        },
        if n_hit[1] > 0 {
            hit[1][0].seed_ratio()
        } else {
            255
        },
    ];
    let min_seed_ratio = seed_ratio[0].min(seed_ratio[1]);
    mb_pair_hits(opt, l2b, *n_hit, hit, pes, &mut paux);
    let do_matesw = !(paux.n_pp > 0 && paux.score == paux.sub_sc);
    if do_matesw && opt.max_rescue > 0 {
        let sub_diff = (opt.a + opt.b).max(opt.q + opt.e);
        if mb_matesw(scratch, opt, l2b, n_hit, hit, pes, &paux, qlen, qseq, is_meth) > 0 {
            for r in 0..2usize {
                for h in hit[r].iter_mut().take(n_hit[r] as usize) {
                    if h.rescued() == 0 {
                        h.n_sub = 0;
                        h.subsc = 0;
                        if let Some(p) = &mut h.p {
                            p.dp_max2 = 0;
                        }
                    }
                }
                mb_hit_sort(&mut n_hit[r], &mut hit[r]);
                mb_set_parent(opt.mask_level,
                    opt.mask_len,
                    n_hit[r],
                    &mut hit[r],
                    sub_diff,
                    0,
                );
                mb_set_mapq(qlen[r],
                    n_hit[r],
                    &mut hit[r],
                    opt.min_chain_score,
                    opt.a,
                    mb_is_sr_mode(opt, qlen[r]),
                    opt.max_sr_len,
                );
            }
            mb_pair_hits(opt, l2b, *n_hit, hit, pes, &mut paux);
        }
    }
    if paux.n_pp != 0 {
        let mut dp_max_se = [0i32; 2];
        let mut dp_max_se2 = [0i32; 2];
        for r in 0..2usize {
            for h in hit[r].iter().take(n_hit[r] as usize) {
                let dp = h.p.as_ref().expect("mb_pair requires hit[r][i]->p").dp_max;
                if dp_max_se[r] < dp {
                    dp_max_se2[r] = dp_max_se[r];
                    dp_max_se[r] = dp;
                } else if dp_max_se2[r] < dp {
                    dp_max_se2[r] = dp;
                }
            }
        }
        let score_se = dp_max_se[0] + dp_max_se[1];
        let score_se2 = dp_max_se2[0] + dp_max_se2[1];
        mb_sync_high_cov(n_hit[0], &mut hit[0]);
        mb_sync_high_cov(n_hit[1], &mut hit[1]);
        if paux.score >= score_se - opt.pen_unpair * opt.a {
            let i0 = paux.i[0] as usize;
            let i1 = paux.i[1] as usize;
            let mut score2 = paux.sub_sc;
            let identity = (hit[0][i0].mlen + hit[1][i1].mlen) as f64
                / (hit[0][i0].blen + hit[1][i1].blen) as f64;
            if (hit[0][i0].id != hit[0][i0].parent || hit[1][i1].id != hit[1][i1].parent)
                && score2 < score_se - opt.pen_unpair * opt.a
            {
                score2 = score_se - opt.pen_unpair * opt.a;
            }
            let frac_high =
                hit[0][i0].frac_high() as f64 / 255.0 + hit[1][i1].frac_high() as f64 / 255.0;
            let mut diff = paux.score - score2;
            let cap = paux.score + PE_BONUS * opt.a - score_se2;
            if diff > cap {
                diff = cap;
            }
            let mut mapq_pe = (6.02 * identity * identity * diff as f64 / opt.a as f64
                - 4.343 * ((paux.n_sub + 1) as f64).ln()
                + 0.499) as i32;
            if mapq_pe < 0 {
                mapq_pe = 0;
            }
            mapq_pe = (mapq_pe as f64 * (1.0 - 0.5 * frac_high) + 0.499) as i32;
            if min_seed_ratio < 50 {
                let sr = min_seed_ratio as f64;
                mapq_pe = (mapq_pe as f64 * sr * sr / 2500.0) as i32;
            }
            if mapq_pe > 60 {
                mapq_pe = 60;
            }
            if mapq_pe <= 0 && paux.score > score2 {
                mapq_pe = 1;
            }
            for (r, &idx) in [i0, i1].iter().enumerate() {
                if n_hit[r] == 1 {
                    if hit[r][idx].mapq < mapq_pe {
                        hit[r][idx].mapq = mapq_pe;
                    }
                } else if hit[r][idx].mapq > mapq_pe {
                    hit[r][idx].mapq = mapq_pe;
                } else {
                    hit[r][idx].mapq =
                        (0.2 * hit[r][idx].mapq as f64 + 0.8 * mapq_pe as f64 + 0.499) as i32;
                }
                if hit[r][idx].id != hit[r][idx].parent {
                    let old_parent = hit[r][idx].parent;
                    let new_parent = hit[r][idx].id;
                    let parent_id = hit[r][old_parent as usize].id;
                    for child in hit[r].iter_mut().take(n_hit[r] as usize) {
                        if child.parent == parent_id {
                            child.parent = new_parent;
                        }
                    }
                    hit[r][old_parent as usize].mapq = 0;
                }
                let q = hit[r][idx].clone();
                for i in 0..n_hit[r] as usize {
                    if i == idx || hit[r][i].id != hit[r][i].parent {
                        continue;
                    }
                    let ol = if hit[r][i].qe <= q.qs || hit[r][i].qs >= q.qe {
                        0
                    } else {
                        hit[r][i].qe.min(q.qe) - hit[r][i].qs.max(q.qs)
                    };
                    if ol as f32 > opt.mask_level * (hit[r][i].qe - hit[r][i].qs) as f32 {
                        let old_id = hit[r][i].id;
                        for j in 0..n_hit[r] as usize {
                            if hit[r][j].parent == old_id {
                                hit[r][j].parent = q.id;
                            }
                        }
                        hit[r][i].mapq = 0;
                    }
                }
            }
        } else {
            let diff = score_se - opt.pen_unpair * opt.a - paux.score;
            let mapq_pe = 6 * diff / opt.a;
            for r in 0..2usize {
                for h in hit[r].iter_mut().take(n_hit[r] as usize) {
                    if h.mapq > mapq_pe {
                        h.mapq = mapq_pe;
                    }
                }
            }
        }
    }
    mb_set_sam_pri(
        n_hit[0],
        &mut hit[0],
        ((opt.flag & MB_F_PRIMARY5) != 0) as i32,
    );
    mb_set_sam_pri(
        n_hit[1],
        &mut hit[1],
        ((opt.flag & MB_F_PRIMARY5) != 0) as i32,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_dir_matches_orientation_and_distance() {
        let h0 = Hit {
            ts: 100,
            te: 150,
            ..Default::default()
        };
        let h1 = Hit {
            ts: 300,
            te: 350,
            flags: 1,
            ..Default::default()
        };
        let mut dist = 0;
        assert_eq!(mb_insert_dir(&h0, &h1, &mut dist), 1);
        assert_eq!(dist, 250);
        assert_eq!(mb_insert_dir(&h1, &h0, &mut dist), 1);
        assert_eq!(dist, 250);
    }

    #[test]
    fn pair_score_uses_proper_orientation_distribution() {
        let h0 = Hit {
            ts: 100,
            te: 150,
            p: Some(Box::new(Extra {
                    dp_max: 80,
                    ..Default::default()
                })),
            ..Default::default()
        };
        let h1 = Hit {
            ts: 300,
            te: 350,
            flags: 1,
            p: Some(Box::new(Extra {
                    dp_max: 70,
                    ..Default::default()
                })),
            ..Default::default()
        };
        let mut pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        pes[1] = Pestat {
            lo: 200,
            hi: 300,
            failed: 0,
            avg: 250.0,
            std: 25.0,
        };
        let s = mb_pair_score(&h0, &h1, &pes, 2);
        assert!(s > 140.0 && s < 151.0);
        pes[1].failed = 1;
        assert_eq!(mb_pair_score(&h0, &h1, &pes, 2), -1.0);
    }

    #[test]
    #[should_panic(expected = "mb_pair_score requires h0->p")]
    fn pair_score_requires_alignment_extra_like_c_deref() {
        let h0 = Hit {
            ts: 100,
            te: 150,
            ..Default::default()
        };
        let h1 = Hit {
            ts: 300,
            te: 350,
            flags: 1,
            p: Some(Box::new(Extra {
                    dp_max: 70,
                    ..Default::default()
                })),
            ..Default::default()
        };
        let mut pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        pes[1] = Pestat {
            lo: 200,
            hi: 300,
            failed: 0,
            avg: 250.0,
            std: 25.0,
        };
        let _ = mb_pair_score(&h0, &h1, &pes, 2);
    }

    #[test]
    fn hit_sum_score_merges_overlapping_primary_query_ranges() {
        let hits = vec![
            Hit {
                id: 0,
                parent: 0,
                qs: 0,
                qe: 100,
                p: Some(Box::new(Extra {
                        dp_max: 100,
                        ..Default::default()
                    })),
                ..Default::default()
            },
            Hit {
                id: 1,
                parent: 1,
                qs: 50,
                qe: 150,
                p: Some(Box::new(Extra {
                        dp_max: 80,
                        ..Default::default()
                    })),
                ..Default::default()
            },
            Hit {
                id: 2,
                parent: 0,
                qs: 0,
                qe: 150,
                p: Some(Box::new(Extra {
                        dp_max: 1000,
                        ..Default::default()
                    })),
                ..Default::default()
            },
        ];
        assert_eq!(mb_hit_sum_score(hits.len() as i32, &hits), 140);
    }

    #[test]
    fn pair_hits_finds_best_proper_pair_and_marks_hits() {
        let l2b = L2b {
            n_ctg: 1,
            ctg: vec![crate::minibwa::l2bit::L2bCtg {
                name: "ctg".into(),
                len: 2000,
                off: 0,
                comm: None,
            }],
            ..Default::default()
        };
        let mut opt = MapOpts::default();
        crate::minibwa::options::mb_opt_init(&mut opt);
        let mut hit = [
            vec![Hit {
                tid: 0,
                ts: 100,
                te: 150,
                id: 0,
                parent: 0,
                hash: 1,
                p: Some(Box::new(Extra {
                        dp_max: 80,
                        ..Default::default()
                    })),
                ..Default::default()
            }],
            vec![Hit {
                tid: 0,
                ts: 300,
                te: 350,
                flags: 1,
                id: 0,
                parent: 0,
                hash: 2,
                p: Some(Box::new(Extra {
                        dp_max: 70,
                        ..Default::default()
                    })),
                ..Default::default()
            }],
        ];
        let mut pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        pes[1] = Pestat {
            lo: 200,
            hi: 300,
            failed: 0,
            avg: 250.0,
            std: 25.0,
        };
        let mut ret = PairAux::default();
        mb_pair_hits(&opt, &l2b, [1, 1], &mut hit, &pes, &mut ret);
        assert_eq!(ret.n_pp, 1);
        assert_eq!(ret.i, [0, 0]);
        assert!(ret.score > 140);
        assert_eq!((hit[0][0].proper_pair(), hit[1][0].proper_pair()), (1, 1));
    }

    #[test]
    #[should_panic(expected = "mb_pair requires hit[r][i]->p")]
    fn pair_final_single_end_score_requires_alignment_extra_like_c_deref() {
        let l2b = L2b {
            n_ctg: 1,
            ctg: vec![crate::minibwa::l2bit::L2bCtg {
                name: "ctg".into(),
                len: 2000,
                off: 0,
                comm: None,
            }],
            ..Default::default()
        };
        let mut opt = MapOpts::default();
        crate::minibwa::options::mb_opt_init(&mut opt);
        opt.max_rescue = 0;
        let mut n_hit = [2, 1];
        let mut hit = [
            vec![
                Hit {
                    tid: 0,
                    ts: 100,
                    te: 150,
                    id: 0,
                    parent: 0,
                    hash: 1,
                    p: Some(Box::new(Extra {
                            dp_max: 80,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
                Hit {
                    tid: 0,
                    ts: 1000,
                    te: 1050,
                    id: 1,
                    parent: 1,
                    hash: 3,
                    p: None,
                    ..Default::default()
                },
            ],
            vec![Hit {
                tid: 0,
                ts: 300,
                te: 350,
                flags: 1,
                id: 0,
                parent: 0,
                hash: 2,
                p: Some(Box::new(Extra {
                        dp_max: 70,
                        ..Default::default()
                    })),
                ..Default::default()
            }],
        ];
        let mut pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        pes[1] = Pestat {
            lo: 200,
            hi: 300,
            failed: 0,
            avg: 250.0,
            std: 25.0,
        };

        mb_pair(&mut Ksw2Scratch::default(), &opt,
            &l2b,
            &mut n_hit,
            &mut hit,
            &pes,
            [12, 12],
            ["ACGTACGTACGT", "TGCATGCATGCA"],
        );
    }

    #[test]
    #[should_panic(expected = "n_pri > 0")]
    fn pair_lift_uses_parent_hit_id_even_when_that_leaves_no_primary_like_c_assert() {
        let l2b = L2b {
            n_ctg: 1,
            ctg: vec![crate::minibwa::l2bit::L2bCtg {
                name: "ctg".into(),
                len: 2000,
                off: 0,
                comm: None,
            }],
            ..Default::default()
        };
        let mut opt = MapOpts::default();
        crate::minibwa::options::mb_opt_init(&mut opt);
        opt.max_rescue = 0;
        let mut n_hit = [3, 1];
        let mut hit = [
            vec![
                Hit {
                    tid: 0,
                    ts: 800,
                    te: 850,
                    id: 7,
                    parent: 7,
                    mapq: 50,
                    hash: 7,
                    p: Some(Box::new(Extra {
                            dp_max: 40,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
                Hit {
                    tid: 0,
                    ts: 100,
                    te: 150,
                    id: 42,
                    parent: 0,
                    mlen: 10,
                    blen: 10,
                    hash: 11,
                    p: Some(Box::new(Extra {
                            dp_max: 100,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
                Hit {
                    tid: 0,
                    ts: 900,
                    te: 950,
                    id: 8,
                    parent: 7,
                    hash: 8,
                    p: Some(Box::new(Extra {
                            dp_max: 30,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
            ],
            vec![Hit {
                tid: 0,
                ts: 300,
                te: 350,
                flags: 1,
                id: 0,
                parent: 0,
                mlen: 10,
                blen: 10,
                hash: 12,
                p: Some(Box::new(Extra {
                        dp_max: 100,
                        ..Default::default()
                    })),
                ..Default::default()
            }],
        ];
        let mut pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        pes[1] = Pestat {
            lo: 200,
            hi: 300,
            failed: 0,
            avg: 250.0,
            std: 25.0,
        };

        mb_pair(&mut Ksw2Scratch::default(), &opt,
            &l2b,
            &mut n_hit,
            &mut hit,
            &pes,
            [12, 12],
            ["ACGTACGTACGT", "TGCATGCATGCA"],
        );

        let _ = hit;
    }

    #[test]
    fn pair_caps_unpaired_mapq_against_best_pair() {
        let l2b = L2b {
            n_ctg: 1,
            ctg: vec![crate::minibwa::l2bit::L2bCtg {
                name: "ctg".into(),
                len: 2000,
                off: 0,
                comm: None,
            }],
            ..Default::default()
        };
        let mut opt = MapOpts::default();
        crate::minibwa::options::mb_opt_init(&mut opt);
        opt.max_rescue = 0;
        let mut n_hit = [2, 2];
        let mut hit = [
            vec![
                Hit {
                    tid: 0,
                    ts: 1000,
                    te: 1050,
                    id: 0,
                    parent: 0,
                    mlen: 50,
                    blen: 50,
                    mapq: 60,
                    hash: 1,
                    p: Some(Box::new(Extra {
                            dp_max: 100,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
                Hit {
                    tid: 0,
                    ts: 100,
                    te: 150,
                    id: 1,
                    parent: 1,
                    mlen: 50,
                    blen: 50,
                    mapq: 60,
                    hash: 2,
                    p: Some(Box::new(Extra {
                            dp_max: 80,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
            ],
            vec![
                Hit {
                    tid: 0,
                    ts: 1500,
                    te: 1550,
                    flags: 1,
                    id: 0,
                    parent: 0,
                    mlen: 50,
                    blen: 50,
                    mapq: 60,
                    hash: 4,
                    p: Some(Box::new(Extra {
                            dp_max: 100,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
                Hit {
                    tid: 0,
                    ts: 300,
                    te: 350,
                    flags: 1,
                    id: 1,
                    parent: 1,
                    mlen: 50,
                    blen: 50,
                    mapq: 60,
                    hash: 8,
                    p: Some(Box::new(Extra {
                            dp_max: 80,
                            ..Default::default()
                        })),
                    ..Default::default()
                },
            ],
        ];
        let mut pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        pes[1] = Pestat {
            lo: 200,
            hi: 300,
            failed: 0,
            avg: 200.0,
            std: 25.0,
        };

        mb_pair(&mut Ksw2Scratch::default(), &opt,
            &l2b,
            &mut n_hit,
            &mut hit,
            &pes,
            [50, 50],
            ["ACGT", "TGCA"],
        );

        assert!(hit[0][0].mapq < 60);
        assert_eq!(hit[0][0].mapq, hit[1][0].mapq);
    }

    #[test]
    fn ungap_detects_single_diagonal_kmer_support() {
        let qseq = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3];
        let tseq = [4, 4, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 4];
        let mut max_i = 0;
        let mut n_good = 0;
        let mut n_kmer = 0;
        let max = mb_ungap(qseq.len() as i32,
            &qseq,
            tseq.len() as i32,
            &tseq,
            3,
            L2bMeth::None,
            &mut max_i,
            &mut n_good,
            &mut n_kmer,
        );
        assert!(max > 0);
        assert_eq!(max_i, 2);
        assert!(n_good >= 1);
        assert!(n_kmer > 0);
    }

    #[test]
    #[should_panic(expected = "mb_matesw requires h0[0]->p")]
    fn matesw_candidate_selection_requires_alignment_extra_like_c_deref() {
        let l2b = L2b {
            n_ctg: 1,
            ctg: vec![crate::minibwa::l2bit::L2bCtg {
                name: "ctg".into(),
                len: 1000,
                off: 0,
                comm: None,
            }],
            ..Default::default()
        };
        let mut opt = MapOpts::default();
        crate::minibwa::options::mb_opt_init(&mut opt);
        opt.max_rescue = 10;
        opt.pen_unpair = 5;
        let mut n_hit = [1, 1];
        let mut hit = [
            vec![Hit {
                tid: 0,
                ts: 100,
                te: 150,
                id: 0,
                parent: 0,
                p: None,
                ..Default::default()
            }],
            vec![Hit {
                tid: 0,
                ts: 300,
                te: 350,
                id: 0,
                parent: 0,
                p: Some(Box::new(Extra {
                        dp_max: 80,
                        ..Default::default()
                    })),
                ..Default::default()
            }],
        ];
        let pes = [Pestat {
            failed: 1,
            ..Default::default()
        }; 4];
        let paux = PairAux {
            score: -1,
            sub_sc: -1,
            ..Default::default()
        };

        let n_add = mb_matesw(&mut Ksw2Scratch::default(), &opt,
            &l2b,
            &mut n_hit,
            &mut hit,
            &pes,
            &paux,
            [12, 12],
            ["ACGTACGTACGT", "TGCATGCATGCA"],
            0,
        );
    }
}
