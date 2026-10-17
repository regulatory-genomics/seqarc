#![allow(unused_variables)]

use crate::minibwa::l2bit::L2b;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Anchor {
    pub sid: i32,
    pub len: i32,
    pub qpos: i32,
    pub flag: u32,
    pub flt: u32,
    pub tpos: i64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Mb128 {
    pub x: u64,
    pub y: u64,
}

fn radix_insert_sort_mb128x(a: &mut [Mb128]) {
    for i in 1..a.len() {
        if a[i].x < a[i - 1].x {
            let tmp = a[i];
            let mut j = i;
            while j > 0 && tmp.x < a[j - 1].x {
                a[j] = a[j - 1];
                j -= 1;
            }
            a[j] = tmp;
        }
    }
}

fn radix_sort_mb128x_rec(a: &mut [Mb128], n_bits: u32, s: u32) {
    const RS_MIN_SIZE: usize = 64;
    let size = 1usize << n_bits;
    let mask = size as u64 - 1;
    let mut counts = [0usize; 256];
    for x in a.iter() {
        counts[((x.x >> s) & mask) as usize] += 1;
    }
    let mut starts = [0usize; 256];
    let mut ends = [0usize; 256];
    let mut sum = 0usize;
    for k in 0..size {
        starts[k] = sum;
        sum += counts[k];
        ends[k] = sum;
    }
    let mut b = starts;
    let mut k = 0usize;
    while k < size {
        if b[k] != ends[k] {
            let mut l = ((a[b[k]].x >> s) & mask) as usize;
            if l != k {
                let mut tmp = a[b[k]];
                loop {
                    std::mem::swap(&mut tmp, &mut a[b[l]]);
                    b[l] += 1;
                    l = ((tmp.x >> s) & mask) as usize;
                    if l == k {
                        break;
                    }
                }
                a[b[k]] = tmp;
                b[k] += 1;
            } else {
                b[k] += 1;
            }
        } else {
            k += 1;
        }
    }
    if s != 0 {
        let next_s = s.saturating_sub(n_bits);
        for k in 0..size {
            let st = starts[k];
            let en = ends[k];
            if en - st > RS_MIN_SIZE {
                radix_sort_mb128x_rec(&mut a[st..en], n_bits, next_s);
            } else if en - st > 1 {
                radix_insert_sort_mb128x(&mut a[st..en]);
            }
        }
    }
}

pub(crate) fn radix_sort_mb128x(a: &mut [Mb128]) {
    const RS_MIN_SIZE: usize = 64;
    if a.len() <= RS_MIN_SIZE {
        radix_insert_sort_mb128x(a);
    } else {
        radix_sort_mb128x_rec(a, 8, 56);
    }
}

pub(crate) fn mb_log2(x: f32) -> f32 {
    let mut i = x.to_bits();
    let mut log_2 = ((i >> 23) & 255) as f32 - 128.0;
    i &= !(255 << 23);
    i += 127 << 23;
    let f = f32::from_bits(i);
    log_2 += (-0.34484843f32 * f + 2.02466578f32) * f - 0.67487759f32;
    log_2
}

#[inline(always)]
fn penalty_to_i64(x: f32) -> i64 {
    debug_assert!(x.is_finite());
    debug_assert!(x >= i64::MIN as f32 && x <= i64::MAX as f32);
    unsafe { x.to_int_unchecked::<i64>() }
}

/// Anchor format in minibwa (`Anchor`):
///   - `tid2`: target (reference) sequence ID
///   - `len`:  length of the seed (q_span)
///   - `qpos`: query coordinate of the last base
///   - `tpos`: target coordinate of the last base
///   - `qocc`: query occurrence
///   - `tocc`: target occurrence
///   - `flt`:  filtered flag
///
/// For chaining, we need:
///   - Target position (`tid2`, `tpos`) for gap calculation on reference
///   - Query position (`qpos`) for gap calculation on query
///   - Seed length (`len`) as q_span for scoring
///
/// Original C static function `mb_chain_bk_end` from `minibwa/lchain.c:25`.
pub(crate) fn mb_chain_bk_end(
    max_drop: i32,
    z: &[Mb128],
    f: &[i32],
    p: &[i64],
    t: &mut [i32],
    k: i64,
) -> i64 {
    let mut i = z[k as usize].y as i64;
    let mut end_i: i64;
    let mut max_i = i;
    let mut max_s = 0i32;
    if i < 0 || t[i as usize] != 0 {
        return i;
    }
    loop {
        t[i as usize] = 2;
        i = p[i as usize];
        end_i = i;
        let s = if i < 0 {
            z[k as usize].x as i32
        } else {
            z[k as usize].x as i32 - f[i as usize]
        };
        if s > max_s {
            max_s = s;
            max_i = i;
        } else if max_s - s > max_drop {
            break;
        }
        if !(i >= 0 && t[i as usize] == 0) {
            break;
        }
    }
    i = z[k as usize].y as i64;
    while i >= 0 && i != end_i {
        t[i as usize] = 0;
        i = p[i as usize];
    }
    max_i
}

/// Original C static function `mb_chain_backtrack` from `minibwa/lchain.c:43`.
pub(crate) fn mb_chain_backtrack(
    n: i64,
    f: &[i32],
    p: &[i64],
    v: &mut [i32],
    t: &mut [i32],
    min_sc: i32,
    max_drop: i32,
    n_u_: &mut i32,
    n_v_: &mut i32,
    u: &mut Vec<u64>,
) {
    *n_u_ = 0;
    *n_v_ = 0;
    u.clear();
    let n_z = (0..n as usize).filter(|&i| f[i] >= min_sc).count();
    if n_z == 0 {
        return;
    }
    let mut z = Vec::with_capacity(n_z);
    for i in 0..n as usize {
        if f[i] >= min_sc {
            z.push(Mb128 {
                x: f[i] as u64,
                y: i as u64,
            });
        }
    }
    radix_sort_mb128x(&mut z);

    t[..n as usize].fill(0);
    let mut n_v = 0i64;
    let mut n_u = 0i32;
    for k in (0..n_z as i64).rev() {
        if t[z[k as usize].y as usize] == 0 {
            let n_v0 = n_v;
            let end_i = mb_chain_bk_end(max_drop, &z, f, p, t, k);
            let mut i = z[k as usize].y as i64;
            while i != end_i {
                n_v += 1;
                t[i as usize] = 1;
                i = p[i as usize];
            }
            let sc = if i < 0 {
                z[k as usize].x as i32
            } else {
                z[k as usize].x as i32 - f[i as usize]
            };
            if sc >= min_sc && n_v > n_v0 {
                n_u += 1;
            } else {
                n_v = n_v0;
            }
        }
    }

    u.resize(n_u as usize, 0);
    t[..n as usize].fill(0);
    n_v = 0;
    n_u = 0;
    for k in (0..n_z as i64).rev() {
        if t[z[k as usize].y as usize] == 0 {
            let n_v0 = n_v;
            let end_i = mb_chain_bk_end(max_drop, &z, f, p, t, k);
            let mut i = z[k as usize].y as i64;
            while i != end_i {
                v[n_v as usize] = i as i32;
                n_v += 1;
                t[i as usize] = 1;
                i = p[i as usize];
            }
            let sc = if i < 0 {
                z[k as usize].x as i32
            } else {
                z[k as usize].x as i32 - f[i as usize]
            };
            if sc >= min_sc && n_v > n_v0 {
                u[n_u as usize] = (sc as u64) << 32 | (n_v - n_v0) as u64;
                n_u += 1;
            } else {
                n_v = n_v0;
            }
        }
    }
    debug_assert!(n_v < i32::MAX as i64);
    *n_u_ = n_u;
    *n_v_ = n_v as i32;
}

/// Original C static function `compact_a` from `minibwa/lchain.c:94`.
///
/// Unlike the C original (which returned a fresh allocation to the caller),
/// this compacts the chains back into the input `a` and returns it, so the
/// caller's reserved anchor-buffer capacity survives across reads. `perm` is
/// a caller-owned permutation scratch buffer, reused across calls.
pub(crate) fn compact_a(
    l2b: &L2b,
    n_u: i32,
    u: &mut [u64],
    n_v: i32,
    v: Vec<i32>,
    mut a: Vec<Anchor>,
    perm: &mut Vec<Anchor>,
) -> Vec<Anchor> {
    perm.clear();
    perm.resize(n_v as usize, Anchor::default());
    let b: &mut Vec<Anchor> = perm;
    let mut k = 0usize;
    for i in 0..n_u as usize {
        let k0 = k;
        let ni = u[i] as i32;
        for j in 0..ni as usize {
            b[k] = a[v[k0 + (ni as usize - j - 1)] as usize];
            k += 1;
        }
    }

    let mut w = vec![Mb128::default(); n_u as usize];
    k = 0;
    for i in 0..n_u as usize {
        let ctg = &l2b.ctg[(b[k].sid >> 1) as usize];
        w[i].x = b[k].tpos as u64 + ctg.off * 2 + ctg.len * (b[k].sid as u64 & 1);
        w[i].y = (k as u64) << 32 | i as u64;
        k += u[i] as u32 as usize;
    }
    radix_sort_mb128x(&mut w);
    let mut u2 = vec![0u64; n_u as usize];
    k = 0;
    for i in 0..n_u as usize {
        let j = w[i].y as u32 as usize;
        let n = u[j] as u32 as usize;
        u2[i] = u[j];
        let src = (w[i].y >> 32) as usize;
        a[k..k + n].copy_from_slice(&b[src..src + n]);
        k += n;
    }
    u[..n_u as usize].copy_from_slice(&u2);
    // a[..k] now holds the compacted chains; returning it (not a fresh copy
    // of b) preserves the caller's allocation.
    a.truncate(k);
    a
}

/// Original C global function `mb_lchain_dp` from `minibwa/lchain.c:163`.
pub(crate) fn mb_lchain_dp(
    l2b: &L2b,
    mut max_dist_x: i32,
    mut max_dist_y: i32,
    bw: i32,
    max_skip: i32,
    max_iter: i32,
    min_sc: i32,
    chn_pen_gap: f32,
    n: i64,
    a: Vec<Anchor>,
    perm: &mut Vec<Anchor>,
    n_u_: &mut i32,
    _u: &mut Vec<u64>,
) -> Vec<Anchor> {
    _u.clear();
    *n_u_ = 0;
    if n == 0 || a.is_empty() {
        return Vec::new();
    }
    if max_dist_x < bw {
        max_dist_x = bw;
    }
    if max_dist_y < bw {
        max_dist_y = bw;
    }
    let n_usize = n as usize;
    let mut v = vec![0i32; n_usize];
    let mut p = vec![0i64; n_usize];
    let mut f = vec![0i32; n_usize];
    let mut t = vec![0i32; n_usize];

    let mut mmax_f = 0i32;
    let max_drop = bw;
    let mut max_ii = -1i64;
    let a_ptr = a.as_ptr();
    let f_ptr = f.as_mut_ptr();
    let p_ptr = p.as_mut_ptr();
    let v_ptr = v.as_mut_ptr();
    let t_ptr = t.as_mut_ptr();
    for i in 0..n {
        let i_usize = i as usize;
        let ai = unsafe { *a_ptr.add(i_usize) };
        let mut max_j = -1i64;
        let mut max_f = ai.len;
        let mut n_skip = 0i32;
        let mut j = i - 1;
        while j >= 0 && j >= i - max_iter as i64 {
            let j_usize = j as usize;
            let aj = unsafe { *a_ptr.add(j_usize) };
            if ai.tpos - aj.tpos >= (max_dist_x + ai.len) as i64 {
                break;
            }
            let mut sc = i32::MIN;
            let dq = (ai.qpos - aj.qpos) as i64;
            if dq > 0 && dq <= (max_dist_y + ai.len) as i64 && ai.sid == aj.sid {
                let dr = ai.tpos - aj.tpos;
                if dr > 0 && dq <= (max_dist_x + ai.len) as i64 {
                    let dd = if dr > dq { dr - dq } else { dq - dr };
                    if dd <= bw as i64 {
                        let dg = if dr < dq { dr } else { dq };
                        let mut sc0 = if (ai.len as i64) < dg {
                            ai.len as i64
                        } else {
                            dg
                        };
                        if dd != 0 {
                            let lin_pen = chn_pen_gap * dd as f32;
                            let log_pen = if dd >= 1 {
                                mb_log2((dd + 1) as f32)
                            } else {
                                0.0
                            };
                            sc0 -= penalty_to_i64(lin_pen + 0.5 * log_pen);
                        }
                        sc = sc0 as i32;
                    }
                }
            }
            if sc != i32::MIN {
                sc += unsafe { *f_ptr.add(j_usize) };
                if sc > max_f {
                    max_f = sc;
                    max_j = j;
                    if n_skip > 0 {
                        n_skip -= 1;
                    }
                } else if unsafe { *t_ptr.add(j_usize) } == i as i32 {
                    n_skip += 1;
                    if n_skip > max_skip {
                        break;
                    }
                }
                let pj = unsafe { *p_ptr.add(j_usize) };
                if pj >= 0 {
                    unsafe {
                        *t_ptr.add(pj as usize) = i as i32;
                    }
                }
            }
            j -= 1;
        }
        let end_j = j;
        if max_ii < 0
            || ai.tpos - unsafe { (*a_ptr.add(max_ii as usize)).tpos }
                > (max_dist_x + ai.len) as i64
        {
            let mut max = i32::MIN;
            max_ii = -1;
            j = i - 1;
            while j >= end_j && j >= 0 {
                let fj = unsafe { *f_ptr.add(j as usize) };
                if max < fj {
                    max = fj;
                    max_ii = j;
                }
                j -= 1;
            }
        }
        if max_ii >= 0 && max_ii < end_j {
            let aj = unsafe { *a_ptr.add(max_ii as usize) };
            let mut tmp = i32::MIN;
            let dq = (ai.qpos - aj.qpos) as i64;
            if dq > 0 && dq <= (max_dist_y + ai.len) as i64 && ai.sid == aj.sid {
                let dr = ai.tpos - aj.tpos;
                if dr > 0 && dq <= (max_dist_x + ai.len) as i64 {
                    let dd = if dr > dq { dr - dq } else { dq - dr };
                    if dd <= bw as i64 {
                        let dg = if dr < dq { dr } else { dq };
                        let mut sc0 = if (ai.len as i64) < dg {
                            ai.len as i64
                        } else {
                            dg
                        };
                        if dd != 0 {
                            let lin_pen = chn_pen_gap * dd as f32;
                            let log_pen = if dd >= 1 {
                                mb_log2((dd + 1) as f32)
                            } else {
                                0.0
                            };
                            sc0 -= penalty_to_i64(lin_pen + 0.5 * log_pen);
                        }
                        tmp = sc0 as i32;
                    }
                }
            }
            let f_max_ii = unsafe { *f_ptr.add(max_ii as usize) };
            if tmp != i32::MIN && max_f < tmp + f_max_ii {
                max_f = tmp + f_max_ii;
                max_j = max_ii;
            }
        }
        unsafe {
            *f_ptr.add(i_usize) = max_f;
            *p_ptr.add(i_usize) = max_j;
        }
        let vi = if max_j >= 0 && unsafe { *v_ptr.add(max_j as usize) } > max_f {
            unsafe { *v_ptr.add(max_j as usize) }
        } else {
            max_f
        };
        unsafe {
            *v_ptr.add(i_usize) = vi;
        }
        if max_ii < 0
            || (ai.tpos - unsafe { (*a_ptr.add(max_ii as usize)).tpos }
                <= (max_dist_x + ai.len) as i64
                && unsafe { *f_ptr.add(max_ii as usize) } < max_f)
        {
            max_ii = i;
        }
        if mmax_f < max_f {
            mmax_f = max_f;
        }
    }

    let mut n_u = 0i32;
    let mut n_v = 0i32;
    mb_chain_backtrack(n, &f, &p, &mut v, &mut t, min_sc, max_drop, &mut n_u, &mut n_v, _u,
    );
    *n_u_ = n_u;
    if n_u == 0 {
        return Vec::new();
    }
    compact_a(l2b, n_u, _u, n_v, v, a, perm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::minibwa::l2bit::{L2bCtg, L2b};

    #[test]
    fn lchain_dp_chains_simple_colinear_anchors() {
        let l2b = L2b {
            tot_len: 1000,
            n_ctg: 1,
            m_ctg: 1,
            ctg: vec![L2bCtg {
                name: "ctg".to_string(),
                comm: None,
                len: 1000,
                off: 0,
            }],
            ..Default::default()
        };
        let anchors = vec![
            Anchor {
                sid: 0,
                len: 10,
                qpos: 9,
                tpos: 109,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 10,
                qpos: 29,
                tpos: 129,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 10,
                qpos: 49,
                tpos: 149,
                ..Default::default()
            },
        ];
        let mut n_u = 0;
        let mut u = Vec::new();
        let out = mb_lchain_dp(&l2b,
            100,
            100,
            20,
            25,
            5000,
            1,
            0.01,
            anchors.len() as i64,
            anchors,
            &mut Vec::new(),
            &mut n_u,
            &mut u,
        );
        assert_eq!(n_u, 1);
        assert_eq!(u.len(), 1);
        assert_eq!(u[0] as u32, 3);
        assert_eq!((u[0] >> 32) as i32, 30);
        assert_eq!(out.len(), 3);
        assert_eq!(
            out.iter().map(|a| a.qpos).collect::<Vec<_>>(),
            vec![9, 29, 49]
        );
    }

    #[test]
    fn lchain_dp_splits_different_targets() {
        let mut l2b = L2b {
            tot_len: 1000,
            n_ctg: 1,
            m_ctg: 1,
            ctg: vec![L2bCtg {
                name: "ctg".to_string(),
                comm: None,
                len: 1000,
                off: 0,
            }],
            ..Default::default()
        };
        l2b.n_ctg = 2;
        l2b.tot_len = 2000;
        l2b.ctg.push(L2bCtg {
            name: "ctg2".to_string(),
            comm: None,
            len: 1000,
            off: 1000,
        });
        let anchors = vec![
            Anchor {
                sid: 0,
                len: 10,
                qpos: 9,
                tpos: 109,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 10,
                qpos: 29,
                tpos: 129,
                ..Default::default()
            },
            Anchor {
                sid: 2,
                len: 10,
                qpos: 9,
                tpos: 209,
                ..Default::default()
            },
            Anchor {
                sid: 2,
                len: 10,
                qpos: 29,
                tpos: 229,
                ..Default::default()
            },
        ];
        let mut n_u = 0;
        let mut u = Vec::new();
        let out = mb_lchain_dp(&l2b,
            100,
            100,
            20,
            25,
            5000,
            1,
            0.01,
            anchors.len() as i64,
            anchors,
            &mut Vec::new(),
            &mut n_u,
            &mut u,
        );
        assert_eq!(n_u, 2);
        assert_eq!(u.iter().map(|x| *x as u32).collect::<Vec<_>>(), vec![2, 2]);
        assert_eq!(out.len(), 4);
    }
}
