#![allow(unused_variables)]

use crate::minibwa::align::mb_align_skeleton_with_scratch;
use crate::minibwa::bwt::{Bwt, SeedHit};
use crate::minibwa::bwt::{SeedHitVec, mb_sai_v_clear};
use crate::minibwa::kommon::KOM_NT4_TABLE;
use crate::minibwa::l2bit::L2b;
use crate::minibwa::lchain::{Mb128, Anchor, mb_lchain_dp, radix_sort_mb128x};
use crate::minibwa::mbpriv::{
    mb_hash64, mb_hash_str, mb_is_sr_mode,
};
use crate::minibwa::options::{mb_opt_adap, MapOpts, MB_F_NO_ALN, MB_F_PRIMARY5};
use crate::minibwa::pe::Hit;
use crate::minibwa::seed::{
    mb_anchor_sort, AnchorVec, mb_anchor_with_scratch, mb_seed_intv_batch,
};

pub(crate) const MB_PARENT_UNSET: i32 = -1;
pub(crate) const MB_PARENT_TMP_PRI: i32 = -2;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MapIndex {
    pub l2b: L2b,
    pub bwt: Bwt,
}

#[derive(Default)]
pub(crate) struct MapBuffer {
    pub anchor_v: AnchorVec,
    pub anchor_aux: Vec<(i64, i64)>,
    pub anchor_sa: Vec<u64>,
    pub anchor_sa_batch: Vec<(u64, u64)>,
    pub anchor_batch: Vec<(i64, i64)>,
    pub chain_w: Vec<u64>,
    /// Permutation scratch for `compact_a`, so anchor compaction does not
    /// allocate per read.
    pub chain_perm: Vec<Anchor>,
    pub align_tseq: Vec<u8>,
    pub ksw: crate::minibwa::ksw2::Ksw2Scratch,
    pub align_qseq0: Vec<u8>,
}

/// Original C global function `mb_idx_ctg_name` from `minibwa/map-algo.c:52`.
pub(crate) fn mb_idx_ctg_name(idx: &MapIndex, tid: i32) -> Option<&str> {
    if tid >= 0 && tid < idx.l2b.n_ctg as i32 {
        Some(&idx.l2b.ctg[tid as usize].name)
    } else {
        None
    }
}

/// Original C global function `mb_tbuf_init` from `minibwa/map-algo.c:59`.
pub(crate) fn mb_tbuf_init() -> MapBuffer {
    MapBuffer {
        ksw: crate::minibwa::ksw2::Ksw2Scratch::default(),
        anchor_v: AnchorVec::default(),
        anchor_aux: Vec::new(),
        anchor_sa: Vec::new(),
        anchor_sa_batch: Vec::new(),
        anchor_batch: Vec::new(),
        chain_w: Vec::new(),
        chain_perm: Vec::new(),
        align_tseq: Vec::new(),
        align_qseq0: Vec::new(),
    }
}

/// Original C global function `mb_cal_mblen` from `minibwa/map-algo.c:97`.
pub(crate) fn mb_cal_mblen(n: i32, a: &[Anchor], blen_: &mut i32) -> i32 {
    *blen_ = 0;
    if n <= 0 {
        return 0;
    }
    let mut mlen = a[0].len as i64;
    let mut blen = a[0].len as i64;
    for i in 1..n as usize {
        let span = a[i].len;
        let tl = a[i].tpos as i32 - a[i - 1].tpos as i32;
        let ql = a[i].qpos - a[i - 1].qpos;
        blen += tl.max(ql) as i64;
        mlen += if tl > span && ql > span {
            span
        } else {
            tl.min(ql)
        } as i64;
    }
    *blen_ = blen as i32;
    mlen as i32
}

/// Original C static function `mb_cal_fuzzy_len` from `minibwa/map-algo.c:115`.
pub(crate) fn mb_cal_fuzzy_len(r: &mut Hit, a: &[Anchor]) {
    r.mlen = mb_cal_mblen(r.cnt, &a[r.as_ as usize..], &mut r.blen);
}

/// Original C static function `mb_hit_set_coor` from `minibwa/map-algo.c:120`.
pub(crate) fn mb_hit_set_coor(r: &mut Hit, qlen: i32, l2b: &L2b, a: &[Anchor]) {
    let k = r.as_ as usize;
    let ak0 = &a[k];
    let ak1 = &a[k + r.cnt as usize - 1];
    r.tid = (ak0.sid >> 1) as i64;
    r.set_rev((ak0.sid & 1) as u8);
    r.ts = ak0.tpos + 1 - ak0.len as i64;
    r.te = ak1.tpos + 1;
    if r.rev() == 0 {
        r.qs = ak0.qpos + 1 - ak0.len;
        r.qe = ak1.qpos + 1;
    } else {
        r.qs = qlen - (ak1.qpos + 1);
        r.qe = qlen - (ak0.qpos + 1 - ak0.len);
    }
    mb_cal_fuzzy_len(r, a);
}

/// Original C global function `mb_cal_high_cov` from `minibwa/map-algo.c:139`.
pub(crate) fn mb_cal_high_cov(n: i32, sai: &[SeedHit], max_occ: i32) -> i32 {
    let mut b = Vec::new();
    for x in sai.iter().take(n as usize) {
        if x.size > max_occ as u64 {
            b.push(x.info);
        }
    }
    if b.is_empty() {
        return 0;
    }
    b.sort_unstable();
    let mut hi_st = (b[0] >> 32) as i32;
    let mut hi_en = b[0] as u32 as i32;
    let mut hi_cov = 0;
    for &x in b.iter().skip(1) {
        let st = (x >> 32) as i32;
        let en = x as u32 as i32;
        if st > hi_en {
            hi_cov += hi_en - hi_st;
            hi_st = st;
            hi_en = en;
        } else {
            hi_en = hi_en.max(en);
        }
    }
    hi_cov + hi_en - hi_st
}

/// Original C global function `mb_sync_high_cov` from `minibwa/map-algo.c:165`.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::minibwa::l2bit::{L2bCtg, L2b};
    use crate::minibwa::options::{mb_opt_init, MapOpts};
    use crate::minibwa::pe::Extra;

    #[test]
    fn mblen_and_hit_coordinates_match_anchor_chain() {
        let anchors = vec![
            Anchor {
                sid: 0,
                len: 10,
                qpos: 19,
                tpos: 109,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 8,
                qpos: 39,
                tpos: 129,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 12,
                qpos: 59,
                tpos: 149,
                ..Default::default()
            },
        ];
        let mut blen = 0;
        assert_eq!(mb_cal_mblen(anchors.len() as i32, &anchors, &mut blen), 30);
        assert_eq!(blen, 50);

        let l2b = L2b {
            tot_len: 1000,
            n_ctg: 1,
            m_ctg: 1,
            ctg: vec![L2bCtg {
                name: "ctg".into(),
                comm: None,
                len: 1000,
                off: 0,
            }],
            ..Default::default()
        };
        let mut hit = Hit {
            cnt: anchors.len() as i32,
            as_: 0,
            ..Default::default()
        };
        mb_hit_set_coor(&mut hit, 100, &l2b, &anchors);
        assert_eq!((hit.tid, hit.rev(), hit.ts, hit.te), (0, 0, 100, 150));
        assert_eq!((hit.qs, hit.qe, hit.mlen, hit.blen), (10, 60, 30, 50));

        let mut rev = hit.clone();
        let mut rev_anchors = anchors.clone();
        for a in &mut rev_anchors {
            a.sid = 1;
        }
        mb_hit_set_coor(&mut rev, 100, &l2b, &rev_anchors);
        assert_eq!((rev.rev(), rev.qs, rev.qe), (1, 40, 90));
    }

    #[test]
    fn generated_hits_are_sorted_and_coordinate_filled() {
        let l2b = L2b {
            tot_len: 1000,
            n_ctg: 1,
            m_ctg: 1,
            ctg: vec![L2bCtg {
                name: "ctg".into(),
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
                qpos: 19,
                tpos: 109,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 8,
                qpos: 39,
                tpos: 129,
                ..Default::default()
            },
            Anchor {
                sid: 0,
                len: 12,
                qpos: 79,
                tpos: 179,
                ..Default::default()
            },
        ];
        let u = vec![(50u64 << 32) | 2, (80u64 << 32) | 1];
        let hits = mb_gen_hit(17, 100, &l2b, u.len() as i32, &u, &anchors);
        assert_eq!(hits.len(), 2);
        assert!(hits[0].score >= hits[1].score);
        assert_eq!((hits[0].id, hits[0].parent), (0, MB_PARENT_UNSET));
        assert_eq!((hits[0].score, hits[0].score0, hits[0].cnt), (80, 80, 1));
        assert_eq!(
            (hits[0].as_, hits[0].qs, hits[0].qe, hits[0].ts, hits[0].te),
            (2, 68, 80, 168, 180)
        );
        assert_eq!((hits[1].score, hits[1].cnt, hits[1].as_), (50, 2, 0));
    }

    #[test]
    fn parent_and_sync_logic_matches_primary_selection() {
        let mut hits = vec![
            Hit {
                id: 5,
                parent: 5,
                qs: 0,
                qe: 100,
                score: 100,
                ..Default::default()
            },
            Hit {
                id: 8,
                parent: 5,
                qs: 10,
                qe: 90,
                score: 80,
                ..Default::default()
            },
            Hit {
                id: 9,
                parent: MB_PARENT_TMP_PRI,
                qs: 120,
                qe: 160,
                score: 40,
                ..Default::default()
            },
        ];
        mb_sync_hits(hits.len() as i32, &mut hits);
        assert_eq!(
            hits.iter().map(|h| (h.id, h.parent)).collect::<Vec<_>>(),
            vec![(0, 0), (1, 0), (2, 2)]
        );
        assert_eq!(
            hits.iter().map(|h| h.sam_pri()).collect::<Vec<_>>(),
            vec![1, 0, 0]
        );

        mb_set_sam_pri(hits.len() as i32, &mut hits, 1);
        assert_eq!(
            hits.iter().map(|h| h.sam_pri()).collect::<Vec<_>>(),
            vec![1, 0, 0]
        );

        let mut parent_hits = vec![
            Hit {
                qs: 0,
                qe: 100,
                score: 100,
                parent: MB_PARENT_UNSET,
                ..Default::default()
            },
            Hit {
                qs: 20,
                qe: 80,
                score: 60,
                parent: MB_PARENT_UNSET,
                ..Default::default()
            },
            Hit {
                qs: 130,
                qe: 170,
                score: 30,
                parent: MB_PARENT_UNSET,
                ..Default::default()
            },
        ];
        mb_set_parent(0.5,
            10,
            parent_hits.len() as i32,
            &mut parent_hits,
            0,
            0,
        );
        assert_eq!(
            parent_hits
                .iter()
                .map(|h| (h.id, h.parent))
                .collect::<Vec<_>>(),
            vec![(0, 0), (1, 0), (2, 2)]
        );
    }

    #[test]
    fn select_sub_sort_and_filter_compact_hits() {
        let mut hits = vec![
            Hit {
                id: 0,
                parent: 0,
                qs: 0,
                qe: 100,
                tid: 0,
                ts: 0,
                te: 100,
                score: 100,
                hash: 3,
                ..Default::default()
            },
            Hit {
                id: 1,
                parent: 0,
                qs: 10,
                qe: 90,
                tid: 0,
                ts: 10,
                te: 90,
                score: 80,
                hash: 1,
                p: Some(Box::new(Extra {
                        dp_max: 90,
                        ..Default::default()
                    }),
                ),
                ..Default::default()
            },
            Hit {
                id: 2,
                parent: 0,
                qs: 0,
                qe: 100,
                tid: 0,
                ts: 0,
                te: 100,
                score: 70,
                hash: 9,
                ..Default::default()
            },
        ];
        let mut n = hits.len() as i32;
        mb_select_sub(0.5, 0, 1, &mut n, &mut hits);
        assert_eq!(n, 2);
        assert_eq!(
            hits.iter().map(|h| (h.id, h.parent)).collect::<Vec<_>>(),
            vec![(0, 0), (1, 0)]
        );

        hits[0].score = 30;
        hits[1].score = 20;
        hits[1].p = Some(Box::new(Extra {
            dp_max: 120,
            ..Default::default()
        }));
        mb_hit_sort(&mut n, &mut hits);
        assert_eq!(hits.iter().map(|h| h.hash).collect::<Vec<_>>(), vec![1, 3]);

        let mut opt = MapOpts::default();
        mb_opt_init(&mut opt);
        opt.min_chain_score = 40;
        opt.min_dp_max = 50;
        opt.a = 2;
        hits[0].mlen = 100;
        hits[1].mlen = 10;
        hits[1].p = Some(Box::new(Extra {
            dp_max: 200,
            ..Default::default()
        }));
        mb_filter_hits(&opt, 100, &mut n, &mut hits);
        assert_eq!(n, 1);
        assert_eq!(hits[0].hash, 1);
    }

    #[test]
    fn squeeze_anchors_updates_hit_offsets_in_offset_order() {
        let mut anchors = vec![
            Anchor {
                qpos: 0,
                ..Default::default()
            },
            Anchor {
                qpos: 1,
                ..Default::default()
            },
            Anchor {
                qpos: 20,
                ..Default::default()
            },
            Anchor {
                qpos: 21,
                ..Default::default()
            },
            Anchor {
                qpos: 40,
                ..Default::default()
            },
        ];
        let mut hits = vec![
            Hit {
                as_: 4,
                cnt: 1,
                ..Default::default()
            },
            Hit {
                as_: 2,
                cnt: 2,
                ..Default::default()
            },
        ];
        let n = mb_squeeze_a(hits.len() as i32, &mut hits, &mut anchors);
        assert_eq!(n, 3);
        assert_eq!(hits.iter().map(|h| h.as_).collect::<Vec<_>>(), vec![2, 0]);
        assert_eq!(
            anchors.iter().take(3).map(|a| a.qpos).collect::<Vec<_>>(),
            vec![20, 21, 40]
        );
    }

}

/// Original C global function `mb_gen_hit` from `minibwa/map-algo.c:174`.
pub(crate) fn mb_gen_hit(
    hash: u32,
    qlen: i32,
    l2b: &L2b,
    n_u: i32,
    u: &[u64],
    a: &[Anchor],
) -> Vec<Hit> {
    if n_u <= 0 {
        return Vec::new();
    }
    let mut z = Vec::with_capacity(n_u as usize);
    let mut k = 0usize;
    for &ui in u.iter().take(n_u as usize) {
        let h = (mb_hash64(
            (mb_hash64(a[k].tpos as u64).wrapping_add(mb_hash64(a[k].qpos as u64))) ^ hash as u64,
        ) & 0xffff_ffff) as u32;
        let x = ui ^ h as u64;
        let y = ((k as u64) << 32) | (ui as u32 as u64);
        z.push(Mb128 { x, y });
        k += ui as u32 as usize;
    }
    radix_sort_mb128x(&mut z);
    z.reverse();

    let mut r = Vec::with_capacity(n_u as usize);
    for (i, zy) in z.iter().enumerate() {
        let mut ri = Hit {
            id: i as i32,
            parent: MB_PARENT_UNSET,
            score: (zy.x >> 32) as i32,
            score0: (zy.x >> 32) as i32,
            hash: zy.x as u32,
            cnt: zy.y as u32 as i32,
            as_: (zy.y >> 32) as i32,
            ..Default::default()
        };
        mb_hit_set_coor(&mut ri, qlen, l2b, a);
        r.push(ri);
    }
    r
}

/// Original C global function `mb_split_hit` from `minibwa/map-algo.c:211`.
pub(crate) fn mb_split_hit(
    r: &mut Hit,
    r2: &mut Hit,
    n: i32,
    qlen: i32,
    a: &[Anchor],
    l2b: &L2b,
) {
    if n <= 0 || n >= r.cnt {
        return;
    }
    *r2 = r.clone();
    r2.id = -1;
    r2.set_sam_pri(0);
    r2.p = None;
    r2.set_split_inv(0);
    r2.cnt = r.cnt - n;
    r2.score = ((r.score as f32 * (r2.cnt as f32 / r.cnt as f32)) + 0.499) as i32;
    r2.as_ = r.as_ + n;
    if r.parent == r.id {
        r2.parent = MB_PARENT_TMP_PRI;
    }
    mb_hit_set_coor(r2, qlen, l2b, a);
    r.cnt -= r2.cnt;
    r.score -= r2.score;
    mb_hit_set_coor(r, qlen, l2b, a);
    r.set_split(r.split() | 1);
    r2.set_split(r2.split() | 2);
}

/// Original C global function `mb_sync_hits` from `minibwa/map-algo.c:230`.
pub(crate) fn mb_sync_hits(n_regs: i32, regs: &mut [Hit]) {
    if n_regs <= 0 {
        return;
    }
    let n = n_regs as usize;
    let mut max_id = -1;
    for r in regs.iter().take(n) {
        max_id = max_id.max(r.id);
    }
    let n_tmp = max_id + 1;
    let mut tmp = vec![-1i32; n_tmp.max(0) as usize];
    for (i, r) in regs.iter().take(n).enumerate() {
        if r.id >= 0 {
            tmp[r.id as usize] = i as i32;
        }
    }
    for i in 0..n {
        let old_parent = regs[i].parent;
        regs[i].id = i as i32;
        if old_parent == MB_PARENT_TMP_PRI {
            regs[i].parent = i as i32;
        } else if old_parent >= 0
            && (old_parent as usize) < tmp.len()
            && tmp[old_parent as usize] >= 0
        {
            regs[i].parent = tmp[old_parent as usize];
        } else {
            regs[i].parent = MB_PARENT_UNSET;
        }
    }
    mb_set_sam_pri(n_regs, regs, 0);
}

/// Original C static function `update_sub` from `minibwa/map-algo.c:258`.
pub(crate) fn update_sub(
    ri: &mut Hit,
    rp: &mut Hit,
    mask_level: f32,
    mask_len: i32,
    sub_diff: i32,
    uncov_len: i32,
) -> i32 {
    let si = ri.qs;
    let ei = ri.qe;
    let sj = rp.qs;
    let ej = rp.qe;
    if ej <= si || sj >= ei {
        return 0;
    }
    let min = (ej - sj).min(ei - si);
    let max = (ej - sj).max(ei - si);
    let ol = ei.min(ej) - si.max(sj);
    if ol as f64 / min as f64 - uncov_len as f64 / max as f64 > mask_level as f64
        && uncov_len <= mask_len
    {
        let mut cnt_sub = 0;
        let mut sci = ri.score;
        ri.parent = rp.parent;
        rp.subsc = rp.subsc.max(sci);
        if rp.p.is_some()
            && ri.p.is_some()
            && (rp.tid != ri.tid || rp.ts != ri.ts || rp.te != ri.te || ol != min)
        {
            let ri_dp_max = ri.p.as_ref().unwrap().dp_max;
            let rp_extra = rp.p.as_mut().unwrap();
            sci = ri_dp_max;
            rp_extra.dp_max2 = rp_extra.dp_max2.max(sci);
            if rp_extra.dp_max - ri_dp_max <= sub_diff {
                cnt_sub = 1;
            }
        }
        if cnt_sub != 0 {
            rp.n_sub += 1;
        }
        1
    } else {
        0
    }
}

/// Original C global function `mb_set_parent` from `minibwa/map-algo.c:279`.
pub(crate) fn mb_set_parent(
    mask_level: f32,
    mask_len: i32,
    n: i32,
    r: &mut [Hit],
    sub_diff: i32,
    hard_mask_level: i32,
) {
    if n <= 0 {
        return;
    }
    let n = n as usize;
    for (i, ri) in r.iter_mut().take(n).enumerate() {
        ri.id = i as i32;
    }
    let mut w = vec![0i32; n];
    let mut k = 1usize;
    r[0].parent = 0;
    for i in 1..n {
        let si = r[i].qs;
        let ei = r[i].qe;
        let mut cov = Vec::new();
        let mut uncov_len = 0i32;
        if hard_mask_level == 0 {
            for &wj in w.iter().take(k) {
                let rp = &r[wj as usize];
                let mut sj = rp.qs;
                let mut ej = rp.qe;
                if ej <= si || sj >= ei {
                    continue;
                }
                if sj < si {
                    sj = si;
                }
                if ej > ei {
                    ej = ei;
                }
                cov.push((sj as u64) << 32 | ej as u32 as u64);
            }
            if cov.is_empty() {
                w[k] = i as i32;
                k += 1;
                r[i].parent = i as i32;
                r[i].n_sub = 0;
                continue;
            }
            cov.sort_unstable();
            let mut x = si;
            for &cv in &cov {
                let st = (cv >> 32) as i32;
                let en = cv as u32 as i32;
                if st > x {
                    uncov_len += st - x;
                }
                x = x.max(en);
            }
            if ei > x {
                uncov_len += ei - x;
            }
        }

        let mut max_ol = 0;
        let mut max_j = -1i32;
        for (j, &wj) in w.iter().take(k).enumerate() {
            let rp = &r[wj as usize];
            let sj = rp.qs;
            let ej = rp.qe;
            let ol = if ej <= si || sj >= ei {
                0
            } else {
                ei.min(ej) - si.max(sj)
            };
            if max_ol < ol {
                max_ol = ol;
                max_j = j as i32;
            }
        }
        let mut n_par = 0;
        if max_j >= 0 {
            let wi = w[max_j as usize] as usize;
            if wi < i {
                let (left, right) = r.split_at_mut(i);
                n_par += update_sub(
                    &mut right[0],
                    &mut left[wi],
                    mask_level,
                    mask_len,
                    sub_diff,
                    uncov_len,
                );
            } else {
                let (left, right) = r.split_at_mut(wi);
                n_par += update_sub(
                    &mut left[i],
                    &mut right[0],
                    mask_level,
                    mask_len,
                    sub_diff,
                    uncov_len,
                );
            }
            let mut j = 0usize;
            while j < k && n_par == 0 {
                let wi = w[j] as usize;
                if wi < i {
                    let (left, right) = r.split_at_mut(i);
                    n_par += update_sub(
                        &mut right[0],
                        &mut left[wi],
                        mask_level,
                        mask_len,
                        sub_diff,
                        uncov_len,
                    );
                } else {
                    let (left, right) = r.split_at_mut(wi);
                    n_par += update_sub(
                        &mut left[i],
                        &mut right[0],
                        mask_level,
                        mask_len,
                        sub_diff,
                        uncov_len,
                    );
                }
                j += 1;
            }
        }
        if n_par == 0 {
            w[k] = i as i32;
            k += 1;
            r[i].parent = i as i32;
            r[i].n_sub = 0;
        }
    }
}

/// Original C global function `mb_set_sam_pri` from `minibwa/map-algo.c:330`.
pub(crate) fn mb_set_sam_pri(n: i32, r: &mut [Hit], is_primary5: i32) {
    if n <= 0 {
        return;
    }
    let mut n_pri = 0;
    let mut min_i = -1;
    let mut min_qs = -1;
    let mut first_i = -1;
    for (i, ri) in r.iter_mut().take(n as usize).enumerate() {
        ri.set_sam_pri(0);
        if ri.id != ri.parent {
            continue;
        }
        n_pri += 1;
        if n_pri == 1 {
            first_i = i as i32;
        }
        if min_qs < 0 || ri.qs < min_qs {
            min_i = i as i32;
            min_qs = ri.qs;
        }
    }
    debug_assert!(n_pri > 0);
    if is_primary5 != 0 {
        r[min_i as usize].set_sam_pri(1);
    } else {
        r[first_i as usize].set_sam_pri(1);
    }
}

/// Original C global function `mb_select_sub` from `minibwa/map-algo.c:346`.
pub(crate) fn mb_select_sub(
    pri_ratio: f32,
    min_diff: i32,
    best_n: i32,
    n_: &mut i32,
    r: &mut Vec<Hit>,
) {
    if pri_ratio > 0.0 && *n_ > 0 {
        let n = *n_ as usize;
        let mut keep = vec![0u8; n];
        let mut n_2nd = 0;
        for i in 0..n {
            let p = r[i].parent;
            if p == i as i32 || r[i].inv() != 0 {
                keep[i] = 1;
            } else if p >= 0
                && (p as usize) < r.len()
                && ((r[i].score as f32) >= r[p as usize].score as f32 * pri_ratio
                    || r[i].score + min_diff >= r[p as usize].score)
                && n_2nd < best_n
                && !(r[i].qs == r[p as usize].qs
                    && r[i].qe == r[p as usize].qe
                    && r[i].tid == r[p as usize].tid
                    && r[i].ts == r[p as usize].ts
                    && r[i].te == r[p as usize].te)
            {
                keep[i] = 1;
                n_2nd += 1;
            }
        }
        let mut k = 0usize;
        for i in 0..n {
            if keep[i] != 0 {
                if k < i {
                    r[k] = std::mem::take(&mut r[i]);
                }
                k += 1;
            } else if r[i].p.is_some() {
                r[i].p = None;
            }
        }
        if k != n {
            r.truncate(k);
            mb_sync_hits(k as i32, r);
        }
        *n_ = k as i32;
    }
}

/// Original C global function `mb_hit_sort` from `minibwa/map-algo.c:366`.
pub(crate) fn mb_hit_sort(n_regs: &mut i32, r: &mut Vec<Hit>) {
    let n = *n_regs as usize;
    if n <= 1 {
        return;
    }
    let mut aux = Vec::new();
    for i in 0..n {
        if r[i].inv() != 0 || r[i].cnt >= 0 {
            let score = r[i].p.as_ref().map(|p| p.dp_max).unwrap_or(r[i].score);
            aux.push(Mb128 {
                x: ((score as u64) << 32) | r[i].hash as u64,
                y: i as u64,
            });
        } else if r[i].p.is_some() {
            r[i].p = None;
        }
    }
    radix_sort_mb128x(&mut aux);
    let mut old = std::mem::take(r);
    let mut t = Vec::with_capacity(aux.len());
    for zy in aux.iter().rev() {
        t.push(std::mem::take(&mut old[zy.y as usize]));
    }
    *n_regs = t.len() as i32;
    *r = t;
}

/// Original C global function `mb_filter_hits` from `minibwa/map-algo.c:394`.
pub(crate) fn mb_filter_hits(opt: &MapOpts, qlen: i32, n_regs: &mut i32, regs: &mut Vec<Hit>) {
    let mut k = 0usize;
    for i in 0..*n_regs as usize {
        let mut flt = regs[i].flt() != 0;
        if regs[i].p.is_some() {
            let dp_max = regs[i].p.as_ref().unwrap().dp_max;
            if regs[i].mlen < opt.min_chain_score {
                flt = true;
            } else if dp_max < opt.min_dp_max * opt.a {
                flt = true;
            }
            if flt {
                regs[i].p = None;
            }
        }
        if !flt {
            if k < i {
                regs[k] = std::mem::take(&mut regs[i]);
            }
            k += 1;
        }
    }
    regs.truncate(k);
    *n_regs = k as i32;
}

/// Original C global function `mb_squeeze_a` from `minibwa/map-algo.c:413`.
pub(crate) fn mb_squeeze_a(n_regs: i32, regs: &mut [Hit], a: &mut [Anchor]) -> i32 {
    let mut aux = Vec::with_capacity(n_regs as usize);
    for (i, r) in regs.iter().take(n_regs as usize).enumerate() {
        aux.push(((r.as_ as u64) << 32) | i as u32 as u64);
    }
    aux.sort_unstable();
    let mut as_ = 0usize;
    for &x in &aux {
        let i = x as u32 as usize;
        let r = &mut regs[i];
        if r.as_ as usize != as_ {
            let src = r.as_ as usize;
            let cnt = r.cnt as usize;
            a.copy_within(src..src + cnt, as_);
            r.as_ = as_ as i32;
        }
        as_ += r.cnt as usize;
    }
    as_ as i32
}

/// Original C static function `mb_set_inv_mapq` from `minibwa/map-algo.c:437`.
pub(crate) fn mb_set_inv_mapq(n_regs: i32, regs: &mut [Hit]) {
    let n = n_regs as usize;
    if n_regs < 3 {
        return;
    }
    if !regs.iter().take(n).any(|r| r.inv() != 0) {
        return;
    }
    let mut aux = Vec::new();
    for (i, r) in regs.iter().take(n).enumerate() {
        if r.parent == i as i32 || r.parent < 0 {
            aux.push(Mb128 {
                x: ((r.tid as u64) << 32) | r.ts as u64,
                y: i as u64,
            });
        }
    }
    radix_sort_mb128x(&mut aux);
    for i in 1..aux.len().saturating_sub(1) {
        let inv_i = aux[i].y as usize;
        if regs[inv_i].inv() != 0 {
            let l_mapq = regs[aux[i - 1].y as usize].mapq;
            let r_mapq = regs[aux[i + 1].y as usize].mapq;
            regs[inv_i].mapq = l_mapq.min(r_mapq);
        }
    }
}

/// Original C global function `mb_set_mapq` from `minibwa/map-algo.c:463`.
pub(crate) fn mb_set_mapq(
    qlen: i32,
    n_regs: i32,
    regs: &mut [Hit],
    min_chain_sc: i32,
    match_sc: i32,
    is_sr: i32,
    max_sr_len: i32,
) {
    const MAPQ_COEF_LEN: i32 = 50;
    const MAPQ_COEF_FAC: f64 = 3.0;
    const Q_COEF: f64 = 40.0;
    if n_regs == 0 {
        return;
    }
    for r in regs.iter_mut().take(n_regs as usize) {
        if r.inv() != 0 {
            r.mapq = 0;
        } else if r.parent == r.id {
            let subsc = r.subsc.max(min_chain_sc);
            let pen_chn = if (r.score as f64) > qlen as f64 * 0.1 {
                1.0
            } else {
                10.0 * r.score as f64 / qlen as f64
            };
            let mut mapq;
            if let Some(p) = &r.p {
                if p.dp_max2 > 0 && p.dp_max > 0 {
                    let identity = r.mlen as f64 / r.blen as f64;
                    let mut x = if r.blen < MAPQ_COEF_LEN {
                        1.0
                    } else {
                        MAPQ_COEF_FAC / (r.blen as f64).ln()
                    };
                    x *= identity * identity;
                    let mapq_sr = (6.02 * x * x * (p.dp_max - p.dp_max2) as f64 / match_sc as f64
                        + 0.499) as i32;
                    x = p.dp_max2 as f64 / p.dp_max as f64;
                    if subsc > r.score0 {
                        x *= subsc as f64 / r.score0 as f64;
                    }
                    let mapq_lr = (pen_chn
                        * identity
                        * Q_COEF
                        * (1.0 - x * x)
                        * (p.dp_max as f64 / match_sc as f64).ln())
                        as i32;
                    if is_sr != 0 {
                        mapq = mapq_sr;
                    } else if max_sr_len < 0 {
                        mapq = mapq_lr;
                    } else {
                        mapq = if qlen < max_sr_len {
                            mapq_sr
                        } else {
                            (mapq_lr as f64
                                - (mapq_lr - mapq_sr) as f64
                                    * 2.0f64.powf(1.0 - qlen as f64 / max_sr_len as f64)
                                + 0.499) as i32
                        };
                    }
                } else {
                    let x = subsc as f64 / r.score0 as f64;
                    let identity = r.mlen as f64 / r.blen as f64;
                    mapq = (pen_chn
                        * identity
                        * Q_COEF
                        * (1.0 - x)
                        * (p.dp_max as f64 / match_sc as f64).ln())
                        as i32;
                }
            } else {
                let x = subsc as f64 / r.score0 as f64;
                mapq = (pen_chn * Q_COEF * (1.0 - x) * (r.score as f64).ln()) as i32;
            }
            mapq -= (4.343f64 * ((r.n_sub + 1) as f64).ln() + 0.499) as i32;
            mapq = mapq.max(0);
            r.mapq = mapq.min(60);
            if let Some(p) = &r.p {
                if p.dp_max > p.dp_max2 && r.mapq == 0 {
                    r.mapq = 1;
                }
            }
        } else {
            r.mapq = 0;
        }
    }
    mb_set_inv_mapq(n_regs, regs);
}

/// Original C global function `mb_map_sai` from `minibwa/map-algo.c:535`.
///
/// `seq` is the NT4-encoded query (same encoding as the seeding path in
/// `mb_map_batch`), so the caller encodes each read once instead of twice.
pub(crate) fn mb_map_sai(
    opt: &MapOpts,
    idx: &MapIndex,
    qlen: i64,
    seq: &[u8],
    u: &mut SeedHitVec,
    n_hit_: &mut i32,
    b: &mut MapBuffer,
    qname: Option<&str>,
) -> Vec<Hit> {
    const MIN_RECHAIN_LEN: i32 = 1000;
    const MIN_RECHAIN_RATIO: f64 = 0.1;
    *n_hit_ = 0;
    if u.n == 0 {
        mb_sai_v_clear(u);
        return Vec::new();
    }
    // Seeds the per-hit tie-break hash in `mb_gen_hit`; the query name is the
    // only functional use of `qname` (the debug printers were removed).
    let mut hash = qname.map(mb_hash_str).unwrap_or(0);
    hash ^= (mb_hash64(qlen as u64).wrapping_add(mb_hash64(opt.seed as u64))) as u32;
    hash = mb_hash64(hash as u64) as u32;
    let hi_cov = mb_cal_high_cov(u.n as i32, &u.a, opt.max_occ);
    let is_sr = mb_is_sr_mode(opt, qlen as i32);
    let sub_diff = (opt.a + opt.b).max(opt.q + opt.e);
    let chn_pen_gap = opt.chain_gap_scale * 0.01 * opt.min_len as f32;

    let mut v = std::mem::take(&mut b.anchor_v);
    mb_anchor_with_scratch(idx,
        u,
        qlen as i32,
        opt.max_occ,
        &mut v,
        &mut b.anchor_aux,
        &mut b.anchor_sa,
        &mut b.anchor_sa_batch,
        &mut b.anchor_batch,
    );
    mb_sai_v_clear(u);

    let mut n_hit = 0;
    let mut w = std::mem::take(&mut b.chain_w);
    w.clear();
    let anchors = std::mem::take(&mut v.a);
    let mut a =  {
        mb_lchain_dp(&idx.l2b,
            opt.max_gap,
            opt.max_gap,
            opt.bw,
            opt.max_chain_skip,
            opt.max_chain_iter,
            opt.min_chain_score,
            chn_pen_gap,
            v.n,
            anchors,
            &mut b.chain_perm,
            &mut n_hit,
            &mut w,
        )
    };

    if opt.bw_long > opt.bw * 2 && is_sr == 0 && n_hit > 0 {
        let mut best = 0usize;
        for i in 1..n_hit as usize {
            if (w[i] >> 32) > (w[best] >> 32) {
                best = i;
            }
        }
        let mut as_ = 0usize;
        for x in w.iter().take(best) {
            as_ += *x as u32 as usize;
        }
        let cnt = w[best] as u32 as usize;
        let st = a[as_].qpos + 1 - a[as_].len;
        let en = a[as_ + cnt - 1].qpos + 1;
        if qlen as i32 - (en - st) > MIN_RECHAIN_LEN
            && (en - st) as f64 > qlen as f64 * MIN_RECHAIN_RATIO
        {
            let n_a = w
                .iter()
                .take(n_hit as usize)
                .map(|x| *x as u32 as i64)
                .sum::<i64>();
            mb_anchor_sort(&idx.l2b, n_a, &mut a);
            w.clear();
            a =  {
                mb_lchain_dp(&idx.l2b,
                    opt.max_gap,
                    opt.max_gap,
                    opt.bw_long,
                    opt.max_chain_skip,
                    opt.max_chain_iter,
                    opt.min_chain_score,
                    chn_pen_gap,
                    n_a,
                    a,
                    &mut b.chain_perm,
                    &mut n_hit,
                    &mut w,
                )
            };
        }
    }

    let mut hit = mb_gen_hit(hash, qlen as i32, &idx.l2b, n_hit, &w, &a);
    mb_set_parent(opt.mask_level,
        opt.mask_len,
        n_hit,
        &mut hit,
        sub_diff,
        0,
    );
    mb_select_sub(opt.pri_ratio,
        opt.min_len * 2,
        opt.best_n,
        &mut n_hit,
        &mut hit,
    );

    if (opt.flag & MB_F_NO_ALN) == 0 {
            mb_align_skeleton_with_scratch(&mut b.ksw,
                opt,
                idx,
                qlen as i32,
                &seq,
                &mut n_hit,
                &mut hit,
                &mut a,
                &mut b.align_tseq,
                &mut b.align_qseq0,
            );
        mb_set_parent(opt.mask_level,
            opt.mask_len,
            n_hit,
            &mut hit,
            sub_diff,
            0,
        );
        mb_select_sub(opt.pri_ratio,
            opt.min_len * 2,
            opt.best_n,
            &mut n_hit,
            &mut hit,
        );
        mb_set_sam_pri(n_hit, &mut hit, ((opt.flag & MB_F_PRIMARY5) != 0) as i32);
    }
        for h in hit.iter_mut().take(n_hit as usize) {
            h.set_frac_high((255.0 * hi_cov as f64 / qlen as f64) as u8);
        }
        mb_set_mapq(qlen as i32,
            n_hit,
            &mut hit,
            opt.min_chain_score,
            opt.a,
            is_sr,
            opt.max_sr_len,
        );
    *n_hit_ = n_hit;
    v.n = 0;
    v.m = a.capacity() as i64;
    v.a = a;
    b.anchor_v = v;
    b.chain_w = w;
    hit
}

/// Align a set of sequences in batch.
///
/// # Arguments
/// - `opt`: options, typically initialized by `mb_opt_init()`
/// - `idx`: index
/// - `n_seq`: number of sequences
/// - `qlen`: query lengths, of size `n_seq`
/// - `seq`: query sequences, ASCII or 0/1/2/3 encoded, of size `n_seq`
/// - `n_hit`: (out) number of hits, of size `n_seq`
/// - `b`: thread buffer; can be `NULL`
/// - `qname`: query name, of size `n_seq`
///
/// # Returns
/// hits, of size `n_seq`
///
/// Original C global function `mb_map_batch` from `minibwa/map-algo.c:656`.
pub(crate) fn mb_map_batch(
    opt: &MapOpts,
    idx: &MapIndex,
    n_seq: i32,
    qlen: &[i32],
    seq: &[&str],
    n_hit: &mut [i32],
    b0: Option<&mut MapBuffer>,
    qname: Option<&[&str]>,
) -> Vec<Vec<Hit>> {
    if n_seq <= 0 {
        return Vec::new();
    }
    let mut owned;
    let b = if let Some(b) = b0 {
        b
    } else {
        owned = mb_tbuf_init();
        &mut owned
    };
    let mut hit = vec![Vec::new(); n_seq as usize];
    let sb_max = opt.sb_seq.min(n_seq).max(0) as usize;
    let mut sai = vec![SeedHitVec::default(); sb_max.max(1)];
    let mut seq4 = vec![Vec::<u8>::new(); sb_max.max(1)];
    let mut i = 0i32;
    let mut sb_st = 0i32;
    let mut sb_len = 0i32;
    while i <= n_seq {
        if i == n_seq || sb_len >= opt.sb_len || i - sb_st >= opt.sb_seq {
            let sb_n = i - sb_st;
            if sb_n == 0 {
                sb_st = i;
                sb_len = 0;
                i += 1;
                continue;
            }
            for k in 0..sb_n as usize {
                let idx_k = sb_st as usize + k;
                seq4[k].clear();
                seq4[k].extend(
                    seq[idx_k]
                        .bytes()
                        .take(qlen[idx_k] as usize)
                        .map(|c| KOM_NT4_TABLE[c as usize]),
                );
                // sai[k] is cleared (capacity-preserving) inside
                // mb_seed_intv_batch; do not reset it here or the grown
                // SeedHitVec capacity is dropped every sub-batch.
            }
            let mut seq4_refs = seq4[..sb_n as usize]
                .iter()
                .map(|s| s.as_ptr())
                .collect::<Vec<_>>();
            mb_seed_intv_batch(&idx.bwt,
                sb_n,
                &qlen[sb_st as usize..],
                &seq4_refs,
                opt.min_len,
                opt.max_sub_occ,
                &mut sai[..sb_n as usize],
            );
            seq4_refs.clear();
            for k in 0..sb_n as usize {
                let idx_k = sb_st as usize + k;
                let mut opt_adap = MapOpts::default();
                mb_opt_adap(opt, qlen[idx_k], &mut opt_adap);
                let name = qname.and_then(|names| names.get(idx_k).copied());
                hit[idx_k] = mb_map_sai(
                    &opt_adap,
                    idx,
                    qlen[idx_k] as i64,
                    &seq4[k],
                    &mut sai[k],
                    &mut n_hit[idx_k],
                    b,
                    name,
                );
            }
            sb_st = i;
            sb_len = 0;
        }
        if i < n_seq {
            sb_len += qlen[i as usize];
        }
        i += 1;
    }
    hit
}
