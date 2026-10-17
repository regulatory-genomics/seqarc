use std::fmt::Write as _;

use crate::minibwa::align::{
    MB_CIGAR_DEL, MB_CIGAR_EQ_MATCH, MB_CIGAR_INS, MB_CIGAR_MATCH, MB_CIGAR_N_SKIP,
    MB_CIGAR_X_MISMATCH,
};
use crate::minibwa::pe::Hit;

/// Original C static function `write_indel_ds` from `minibwa/cs.c:18`.
pub(crate) fn write_indel_ds(s: &mut String, len: i64, seq: &[u8], ll: i64, lr: i64) {
    const BASES: &[u8; 5] = b"acgtn";
    if ll + lr >= len {
        s.push('[');
        for i in 0..len as usize {
            s.push(BASES[seq[i] as usize] as char);
        }
        s.push(']');
    } else {
        let mut k = 0usize;
        if ll > 0 {
            s.push('[');
            for i in 0..ll as usize {
                s.push(BASES[seq[k + i] as usize] as char);
            }
            s.push(']');
            k += ll as usize;
        }
        for i in 0..(len - lr - ll) as usize {
            s.push(BASES[seq[k + i] as usize] as char);
        }
        k += (len - lr - ll) as usize;
        if lr > 0 {
            s.push('[');
            for i in 0..lr as usize {
                s.push(BASES[seq[k + i] as usize] as char);
            }
            s.push(']');
        }
    }
}

/// Original C global function `mb_write_cs_ds` from `minibwa/cs.c:47`.
pub(crate) fn mb_write_cs_ds(s: &mut String, tseq: &[u8], qseq: &[u8], r: &Hit, is_ds: i32) {
    const BASES: &[u8; 5] = b"acgtn";
    let Some(p) = &r.p else {
        return;
    };
    let mut q_len = 0usize;
    let mut t_len = 0usize;
    s.push(if is_ds != 0 { 'd' } else { 'c' });
    s.push_str("s:Z:");
    for &cg in p.cigar().iter().take(p.n_cigar as usize) {
        let op = cg & 0xf;
        let len = (cg >> 4) as usize;
        if op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH {
            q_len += len;
            t_len += len;
        } else if op == MB_CIGAR_INS {
            q_len += len;
        } else if op == MB_CIGAR_DEL || op == MB_CIGAR_N_SKIP {
            t_len += len;
        }
    }
    let mut q_off = 0usize;
    let mut t_off = 0usize;
    for &cg in p.cigar().iter().take(p.n_cigar as usize) {
        let op = cg & 0xf;
        let len = (cg >> 4) as usize;
        assert!(
            (MB_CIGAR_MATCH..=MB_CIGAR_N_SKIP).contains(&op)
                || op == MB_CIGAR_EQ_MATCH
                || op == MB_CIGAR_X_MISMATCH
        );
        if op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH {
            let mut l_tmp = 0i32;
            for j in 0..len {
                if qseq[q_off + j] != tseq[t_off + j] {
                    if l_tmp > 0 {
                        write!(s, ":{l_tmp}").unwrap();
                        l_tmp = 0;
                    }
                    s.push('*');
                    s.push(BASES[tseq[t_off + j] as usize] as char);
                    s.push(BASES[qseq[q_off + j] as usize] as char);
                } else {
                    l_tmp += 1;
                }
            }
            if l_tmp > 0 {
                write!(s, ":{l_tmp}").unwrap();
            }
            q_off += len;
            t_off += len;
        } else if op == MB_CIGAR_INS {
            if is_ds != 0 {
                let y = q_off;
                let mut z = 1usize;
                while z <= len {
                    if y < z || qseq[y + len - z] != qseq[y - z] {
                        break;
                    }
                    z += 1;
                }
                let lr = z - 1;
                z = 0;
                while z < len {
                    if y + len + z >= q_len || qseq[y + len + z] != qseq[y + z] {
                        break;
                    }
                    z += 1;
                }
                let ll = z;
                s.push('+');
                write_indel_ds(s, len as i64, &qseq[y..], ll as i64, lr as i64);
            } else {
                s.push('+');
                for &b in &qseq[q_off..q_off + len] {
                    s.push(BASES[b as usize] as char);
                }
            }
            q_off += len;
        } else if op == MB_CIGAR_DEL {
            if is_ds != 0 {
                let x = t_off;
                let mut z = 1usize;
                while z <= len {
                    if x < z || tseq[x + len - z] != tseq[x - z] {
                        break;
                    }
                    z += 1;
                }
                let lr = z - 1;
                z = 0;
                while z < len {
                    if x + len + z >= t_len || tseq[x + z] != tseq[x + len + z] {
                        break;
                    }
                    z += 1;
                }
                let ll = z;
                s.push('-');
                write_indel_ds(s, len as i64, &tseq[x..], ll as i64, lr as i64);
            } else {
                s.push('-');
                for &b in &tseq[t_off..t_off + len] {
                    s.push(BASES[b as usize] as char);
                }
            }
            t_off += len;
        } else {
            assert!(len >= 2);
            write!(
                s,
                "~{}{}{}{}{}",
                BASES[tseq[t_off] as usize] as char,
                BASES[tseq[t_off + 1] as usize] as char,
                len,
                BASES[tseq[t_off + len - 2] as usize] as char,
                BASES[tseq[t_off + len - 1] as usize] as char,
            )
            .unwrap();
            t_off += len;
        }
    }
    assert_eq!(t_off as i64, r.te - r.ts);
    assert_eq!(q_off as i32, r.qe - r.qs);
}

/// Original C global function `mb_write_md` from `minibwa/cs.c:129`.
pub(crate) fn mb_write_md(s: &mut String, tseq: &[u8], qseq: &[u8], r: &Hit) {
    const BASES: &[u8; 5] = b"ACGTN";
    let Some(p) = &r.p else {
        return;
    };
    s.push_str("MD:Z:");
    let mut q_off = 0usize;
    let mut t_off = 0usize;
    let mut l_md = 0i32;
    for &cg in p.cigar().iter().take(p.n_cigar as usize) {
        let op = cg & 0xf;
        let len = (cg >> 4) as usize;
        assert!(
            (MB_CIGAR_MATCH..=MB_CIGAR_N_SKIP).contains(&op)
                || op == MB_CIGAR_EQ_MATCH
                || op == MB_CIGAR_X_MISMATCH
        );
        if op == MB_CIGAR_MATCH || op == MB_CIGAR_EQ_MATCH || op == MB_CIGAR_X_MISMATCH {
            for j in 0..len {
                if qseq[q_off + j] != tseq[t_off + j] {
                    write!(s, "{l_md}{}", BASES[tseq[t_off + j] as usize] as char).unwrap();
                    l_md = 0;
                } else {
                    l_md += 1;
                }
            }
            q_off += len;
            t_off += len;
        } else if op == MB_CIGAR_INS {
            q_off += len;
        } else if op == MB_CIGAR_DEL {
            write!(s, "{l_md}^").unwrap();
            for &b in &tseq[t_off..t_off + len] {
                s.push(BASES[b as usize] as char);
            }
            l_md = 0;
            t_off += len;
        } else if op == MB_CIGAR_N_SKIP {
            t_off += len;
        }
    }
    if l_md > 0 {
        write!(s, "{l_md}").unwrap();
    }
    assert_eq!(t_off as i64, r.te - r.ts);
    assert_eq!(q_off as i32, r.qe - r.qs);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ds_marks_repeat_context() {
        let mut s = String::new();
        write_indel_ds(&mut s, 3, &[0, 0, 0], 1, 1);
        assert_eq!(s, "[a]a[a]");
    }
}
