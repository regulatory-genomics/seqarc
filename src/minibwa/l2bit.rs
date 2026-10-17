#![allow(unused_variables)]

use crate::minibwa::kommon::KOM_NT4_TABLE;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::minibwa::mapped_u64::U64Storage;

const L2B_MAGIC: &[u8; 4] = b"L2B\x01";

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum L2bMeth {
    #[default]
    None = 0,
    C2t = 1,
    G2a = 2,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct L2bCtg {
    pub name: String,
    pub comm: Option<String>,
    pub len: u64,
    pub off: u64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct L2bIntv {
    pub st: u64,
    pub en: u64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct L2b {
    pub tot_len: u64,
    pub n_ctg: u64,
    pub m_ctg: u64,
    pub ctg: Vec<L2bCtg>,
    pub n_pac: u64,
    pub m_pac: u64,
    pub n_ambi: u64,
    pub m_ambi: u64,
    pub n_mask: u64,
    pub m_mask: u64,
    pub ambi: Vec<L2bIntv>,
    pub mask: Vec<L2bIntv>,
    pub pac: U64Storage,
    pub cat_name: Vec<u8>,
    pub cat_comm: Vec<u8>,
}

/// Original C static function `l2b_pos2cid` from `minibwa/l2bit.c:9`.
pub(crate) fn l2b_pos2cid(l2b: &L2b, s: i64, len: i64, cst: &mut i64) -> i64 {
    let mut lo = 0i64;
    let mut hi = l2b.n_ctg as i64;
    while lo < hi {
        let mid = (lo + hi) / 2;
        let ctg = &l2b.ctg[mid as usize];
        if ctg.off as i64 <= s && s < (ctg.off + ctg.len) as i64 {
            *cst = s - ctg.off as i64;
            return if s + len <= (ctg.off + ctg.len) as i64 {
                mid
            } else {
                -1
            };
        } else if s < ctg.off as i64 {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    -1
}

/// Original C global function `l2b_intv2cid` from `minibwa/l2bit.c:9`.
pub(crate) fn l2b_intv2cid(l2b: &L2b, st: u64, en: u64, cst: &mut i64, rev: &mut i32) -> i64 {
    assert!(st < en);
    if en > l2b.tot_len * 2 {
        return -3;
    }
    if st < l2b.tot_len && l2b.tot_len < en {
        return -2;
    }
    *rev = (st >= l2b.tot_len) as i32;
    let s = if st < l2b.tot_len {
        st
    } else {
        l2b.tot_len * 2 - en
    } as i64;
    l2b_pos2cid(l2b, s, (en - st) as i64, cst)
}

/// Original C global function `l2b_intv2cid_meth` from `minibwa/l2bit.c:36`.
pub(crate) fn l2b_intv2cid_meth(
    l2b: &L2b,
    st: u64,
    en: u64,
    mt: &mut L2bMeth,
    cst: &mut i64,
    rev: &mut i32,
) -> i64 {
    assert!(st < en);
    let len = (en - st) as i64;
    let tot_len = l2b.tot_len;
    if en > tot_len * 4 {
        return -3;
    }
    let copy = (st / tot_len) as i32;
    *mt = if copy == 0 || copy == 3 {
        L2bMeth::C2t
    } else {
        L2bMeth::G2a
    };
    *rev = (copy >= 2) as i32;
    let mut s = (st - tot_len * copy as u64) as i64;
    if copy >= 2 {
        s = tot_len as i64 - len - s;
    }
    if s < 0 {
        -2
    } else {
        l2b_pos2cid(l2b, s, len, cst)
    }
}

/// Original C global function `l2b_getseq` from `minibwa/l2bit.c:29`.
pub(crate) fn l2b_getseq(l2b: &L2b, tid: i64, mut st: i64, mut en: i64, seq: &mut [u8]) -> i64 {
    if tid < 0 || tid >= l2b.n_ctg as i64 {
        return -1;
    }
    let ctg = &l2b.ctg[tid as usize];
    if st < 0 {
        st = 0;
    }
    if en > ctg.len as i64 {
        en = ctg.len as i64;
    }
    st += ctg.off as i64;
    en += ctg.off as i64;
    let len = (en - st).max(0) as usize;
    let pac = l2b.pac.as_ptr();
    let seq_ptr = seq.as_mut_ptr();
    let mut pos = st as u64;
    let mut out = 0usize;
    while out < len {
        let in_word = (pos & 31) as usize;
        let take = (32 - in_word).min(len - out);
        let mut word = unsafe { *pac.add((pos >> 5) as usize) >> (in_word * 2) };
        for j in 0..take {
            unsafe {
                *seq_ptr.add(out + j) = (word & 3) as u8;
            }
            word >>= 2;
        }
        out += take;
        pos += take as u64;
    }
    let mut n_ambi = 0i32;
    let aid = l2b_getambi(
        l2b,
        tid,
        st - ctg.off as i64,
        en - ctg.off as i64,
        &mut n_ambi,
    );
    if aid >= 0 {
        for i in 0..n_ambi as i64 {
            let iv = l2b.ambi[(aid + i) as usize];
            let mut s = iv.st as i64;
            let mut e = iv.en as i64;
            if s < st {
                s = st;
            }
            if e > en {
                e = en;
            }
            if s < e {
                let off = (s - st) as usize;
                let len = (e - s) as usize;
                unsafe {
                    std::ptr::write_bytes(seq_ptr.add(off), 4, len);
                }
            }
        }
    }
    en - st
}

/// Original C global function `l2b_meth_convert` from `minibwa/l2bit.c:93`.
pub(crate) fn l2b_meth_convert(mt: L2bMeth, len: i64, seq: &mut [u8]) {
    if mt == L2bMeth::C2t {
        for i in 0..len as usize {
            if seq[i] == 1 {
                seq[i] = 3;
            }
        }
    } else if mt == L2bMeth::G2a {
        for i in 0..len as usize {
            if seq[i] == 2 {
                seq[i] = 0;
            }
        }
    }
}

/// Original C static function `l2b_get0` from `minibwa/l2bit.h:43`.
#[inline(always)]
pub(crate) fn l2b_get0(l2b: &L2b, i: u64) -> i32 {
    ((l2b.pac[(i >> 5) as usize] >> ((i & 31) << 1)) & 3) as i32
}

/// Original C static function `l2b_meth_rev` from `minibwa/l2bit.h:54`.
pub(crate) fn l2b_meth_rev(mt: L2bMeth) -> L2bMeth {
    if mt == L2bMeth::None {
        L2bMeth::None
    } else if mt == L2bMeth::C2t {
        L2bMeth::G2a
    } else {
        L2bMeth::C2t
    }
}

/// Original C global function `l2b_getambi` from `minibwa/l2bit.c:56`.
pub(crate) fn l2b_getambi(l2b: &L2b, tid: i64, mut st: i64, mut en: i64, n_ambi: &mut i32) -> i64 {
    *n_ambi = 0;
    if tid < 0 || tid >= l2b.n_ctg as i64 {
        return -1;
    }
    if st < 0 {
        st = 0;
    }
    if en > l2b.ctg[tid as usize].len as i64 {
        en = l2b.ctg[tid as usize].len as i64;
    }
    if st >= en {
        return -1;
    }
    let g_beg = l2b.ctg[tid as usize].off as i64 + st;
    let g_end = l2b.ctg[tid as usize].off as i64 + en;

    let mut lo = 0i64;
    let mut hi = l2b.n_ambi as i64;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if l2b.ambi[mid as usize].en as i64 > g_beg {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    let i_st = lo;
    lo = i_st;
    hi = l2b.n_ambi as i64;
    while lo < hi {
        let mid = (lo + hi) / 2;
        if l2b.ambi[mid as usize].st as i64 >= g_end {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    let i_en = lo;
    *n_ambi = (i_en - i_st) as i32;
    if *n_ambi == 0 {
        -1
    } else {
        i_st
    }
}

/// Original C static function `l2b_format_seq` from `minibwa/l2bit.c:88`.
pub(crate) fn l2b_format_seq(len: u64, seq: &mut [u8], rng: &mut u64) {
    for i in 0..len as usize {
        let b = seq[i];
        let mut c = KOM_NT4_TABLE[b as usize];
        if c > 4 {
            c = 4;
        }
        if c == 4 {
            *rng = rng.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = *rng;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            c |= ((z ^ (z >> 31)) & 3) as u8;
        }
        if b < b'A' || b > b'Z' {
            c |= 1 << 3;
        }
        seq[i] = c;
    }
}

/// Original C static function `l2b_add_seq` from `minibwa/l2bit.c:101`.
pub(crate) fn l2b_add_seq(
    l2b: &mut L2b,
    len: u64,
    seq: &[u8],
    name: &str,
    comm: Option<&str>,
    rng: &mut u64,
) {
    let off = l2b.tot_len;
    let name_end = name
        .as_bytes()
        .iter()
        .position(|&b| b == 0)
        .unwrap_or(name.len());
    let comm_end = comm.map(|s| s.as_bytes().iter().position(|&b| b == 0).unwrap_or(s.len()));
    l2b.ctg.push(L2bCtg {
        name: name[..name_end].to_string(),
        comm: comm.zip(comm_end).map(|(s, end)| s[..end].to_string()),
        len,
        off: l2b.tot_len,
    });
    l2b.n_ctg += 1;
    l2b.m_ctg = l2b.ctg.capacity() as u64;
    l2b.tot_len += len;

    l2b.n_pac = (l2b.tot_len + 31) / 32;
    if l2b.pac.len() < l2b.n_pac as usize {
        l2b.pac.resize(l2b.n_pac as usize, 0);
    }
    l2b.m_pac = l2b.pac.capacity() as u64;

    let mut ambi_len = 0u64;
    let mut mask_len = 0u64;
    for i in 0..len {
        let c = seq[i as usize] as u64;
        let x = off + i;
        if (c & (1 << 3)) != 0 {
            mask_len += 1;
        } else if mask_len > 0 {
            l2b.mask.push(L2bIntv {
                st: x - mask_len,
                en: x,
            });
            l2b.n_mask += 1;
            mask_len = 0;
        }
        if (c & (1 << 2)) != 0 {
            ambi_len += 1;
        } else if ambi_len > 0 {
            l2b.ambi.push(L2bIntv {
                st: x - ambi_len,
                en: x,
            });
            l2b.n_ambi += 1;
            ambi_len = 0;
        }
        l2b.pac[(x >> 5) as usize] |= (c & 3) << ((x & 0x1f) * 2);
    }
    l2b.m_mask = l2b.mask.capacity() as u64;
    l2b.m_ambi = l2b.ambi.capacity() as u64;
}

fn l2b_add_seq_bytes(
    l2b: &mut L2b,
    len: u64,
    seq: &[u8],
    name: &[u8],
    comm: Option<&[u8]>,
    rng: &mut u64,
) {
    let name_lossy = String::from_utf8_lossy(name);
    let comm_lossy = comm.map(String::from_utf8_lossy);
    l2b_add_seq(
        l2b,
        len,
        seq,
        name_lossy.as_ref(),
        comm_lossy.as_deref(),
        rng,
    );
}

fn l2b_collate_raw(l2b: &mut L2b, names: &[Vec<u8>], comments: &[Option<Vec<u8>>]) {
    if !l2b.cat_name.is_empty() || !l2b.cat_comm.is_empty() {
        return;
    }
    for name in names {
        l2b.cat_name.extend_from_slice(name);
        l2b.cat_name.push(0);
    }
    for comm in comments {
        if let Some(comm) = comm {
            l2b.cat_comm.extend_from_slice(comm);
        }
        l2b.cat_comm.push(0);
    }
}

struct KseqRecord {
    name: Vec<u8>,
    comment: Option<Vec<u8>>,
    seq: Vec<u8>,
}

fn is_c_isspace_byte(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

fn c_string_prefix(bytes: &[u8]) -> &[u8] {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());
    &bytes[..end]
}

fn c_path<P: AsRef<Path>>(path: P) -> PathBuf {
    PathBuf::from(
        String::from_utf8_lossy(c_string_prefix(path.as_ref().to_string_lossy().as_bytes()))
            .into_owned(),
    )
}

fn c_str_len(s: &str) -> usize {
    c_string_prefix(s.as_bytes()).len()
}

fn strip_line_ending(mut line: Vec<u8>) -> Vec<u8> {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
    line
}

fn read_one<R: BufRead + ?Sized>(reader: &mut R) -> io::Result<Option<u8>> {
    let buf = reader.fill_buf()?;
    if buf.is_empty() {
        return Ok(None);
    }
    let c = buf[0];
    reader.consume(1);
    Ok(Some(c))
}

fn read_kseq_record<R: BufRead + ?Sized>(
    reader: &mut R,
    last_char: &mut Option<u8>,
) -> io::Result<Option<KseqRecord>> {
    match last_char.take() {
        Some(_) => {}
        None => loop {
            match read_one(reader)? {
                Some(b'>' | b'@') => break,
                Some(_) => continue,
                None => return Ok(None),
            }
        },
    }

    let mut header = Vec::new();
    reader.read_until(b'\n', &mut header)?;
    let header = strip_line_ending(header);
    let split = header
        .iter()
        .position(|&b| is_c_isspace_byte(b))
        .unwrap_or(header.len());
    let name = c_string_prefix(&header[..split]).to_vec();
    let comment = if split < header.len() {
        let comm = c_string_prefix(&header[split + 1..]).to_vec();
        if comm.is_empty() {
            None
        } else {
            Some(comm)
        }
    } else {
        None
    };

    let mut seq = Vec::new();
    let mut c = loop {
        match read_one(reader)? {
            Some(b'\n' | b'\r') => continue,
            other => break other,
        }
    };
    while let Some(ch) = c {
        match ch {
            b'>' | b'@' => {
                *last_char = Some(ch);
                return Ok(Some(KseqRecord { name, comment, seq }));
            }
            b'+' => {
                let mut discard = Vec::new();
                reader.read_until(b'\n', &mut discard)?;
                let mut qual_len = 0usize;
                while qual_len < seq.len() {
                    let mut qline = Vec::new();
                    if reader.read_until(b'\n', &mut qline)? == 0 {
                        return Ok(None);
                    }
                    qual_len += strip_line_ending(qline).len();
                }
                if qual_len != seq.len() {
                    return Ok(None);
                }
                return Ok(Some(KseqRecord { name, comment, seq }));
            }
            b'\n' | b'\r' => {}
            _ => {
                seq.push(ch);
                let mut rest = Vec::new();
                reader.read_until(b'\n', &mut rest)?;
                seq.extend_from_slice(&strip_line_ending(rest));
            }
        }
        c = read_one(reader)?;
    }
    Ok(Some(KseqRecord { name, comment, seq }))
}

/// Original C global function `l2b_import` from `minibwa/l2bit.c:175`.
pub(crate) fn l2b_import<P: AsRef<Path>>(fn_: P, seed: u64) -> Option<L2b> {
    let fn_ = c_path(fn_);
    let path = fn_.to_string_lossy();
    let raw: Box<dyn Read> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&fn_).ok()?)
    };
    let mut buffered = BufReader::with_capacity(1 << 20, raw);
    let is_gzip = buffered.fill_buf().ok()?.starts_with(&[0x1f, 0x8b]);
    let mut reader: Box<dyn BufRead> = if is_gzip {
        Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(buffered),
        ))
    } else {
        Box::new(buffered)
    };
    let mut l2b = L2b::default();
    let mut rng = seed;
    let mut last_char = None;
    let mut names = Vec::new();
    let mut comments = Vec::new();
    while let Some(mut rec) = read_kseq_record(&mut *reader, &mut last_char).ok()? {
        l2b_format_seq(rec.seq.len() as u64, &mut rec.seq, &mut rng);
        l2b_add_seq_bytes(
            &mut l2b,
            rec.seq.len() as u64,
            &rec.seq,
            &rec.name,
            rec.comment.as_deref(),
            &mut rng,
        );
        names.push(rec.name);
        comments.push(rec.comment);
    }
    l2b_collate_raw(&mut l2b, &names, &comments);
    Some(l2b)
}

/// Original C global function `l2b_save` from `minibwa/l2bit.c:202`.
pub(crate) fn l2b_save<P: AsRef<Path>>(fn_: P, l2b: &L2b) -> i32 {
    let fn_ = c_path(fn_);
    let path = fn_.to_string_lossy();
    let mut fp: Box<dyn Write> = if path == "-" {
        Box::new(BufWriter::with_capacity(1 << 20, io::stdout()))
    } else {
        match File::create(&fn_) {
            Ok(fp) => Box::new(BufWriter::with_capacity(1 << 20, fp)),
            Err(_) => return -1,
        }
    };
    let use_cat = l2b_cat_strs_valid(l2b);
    let len_name: u64 = if use_cat {
        l2b.cat_name.len() as u64
    } else {
        l2b.ctg.iter().map(|c| c_str_len(&c.name) as u64 + 1).sum()
    };
    let len_comm: u64 = if use_cat {
        l2b.cat_comm.len() as u64
    } else {
        l2b.ctg
            .iter()
            .map(|c| {
                c.comm
                    .as_ref()
                    .map(|s| c_str_len(s) as u64 + 1)
                    .unwrap_or(1)
            })
            .sum()
    };
    let _ = fp.write_all(L2B_MAGIC);
    let _ = fp.write_all(&0u32.to_le_bytes());
    for x in [
        l2b.n_ctg,
        l2b.tot_len,
        l2b.n_ambi,
        l2b.n_mask,
        len_name,
        len_comm,
        l2b.n_pac,
    ] {
        let _ = fp.write_all(&x.to_le_bytes());
    }
    for ctg in &l2b.ctg {
        let _ = fp.write_all(&ctg.len.to_le_bytes());
    }
    for iv in &l2b.ambi {
        let _ = fp.write_all(&iv.st.to_le_bytes());
        let _ = fp.write_all(&iv.en.to_le_bytes());
    }
    for iv in &l2b.mask {
        let _ = fp.write_all(&iv.st.to_le_bytes());
        let _ = fp.write_all(&iv.en.to_le_bytes());
    }
    let _ = write_u64_slice_le(&mut fp, &l2b.pac);
    if use_cat {
        let _ = fp.write_all(&l2b.cat_name);
        let _ = fp.write_all(&l2b.cat_comm);
    } else {
        for ctg in &l2b.ctg {
            let _ = fp.write_all(c_string_prefix(ctg.name.as_bytes()));
            let _ = fp.write_all(&[0]);
        }
        for ctg in &l2b.ctg {
            if let Some(comm) = &ctg.comm {
                let _ = fp.write_all(c_string_prefix(comm.as_bytes()));
            }
            let _ = fp.write_all(&[0]);
        }
    }
    let _ = fp.flush();
    0
}

fn l2b_cat_strs_valid(l2b: &L2b) -> bool {
    if l2b.cat_name.is_empty() || l2b.cat_comm.is_empty() {
        return false;
    }
    let mut p_name = 0usize;
    let mut p_comm = 0usize;
    for _ in 0..l2b.n_ctg as usize {
        if p_name >= l2b.cat_name.len() {
            return false;
        }
        let Some(name_rel_end) = l2b.cat_name[p_name..].iter().position(|&b| b == 0) else {
            return false;
        };
        p_name += name_rel_end + 1;
        if p_comm >= l2b.cat_comm.len() {
            return false;
        }
        if l2b.cat_comm[p_comm] == 0 {
            p_comm += 1;
        } else {
            let Some(comm_rel_end) = l2b.cat_comm[p_comm..].iter().position(|&b| b == 0) else {
                return false;
            };
            p_comm += comm_rel_end + 1;
        }
    }
    p_name == l2b.cat_name.len() && p_comm == l2b.cat_comm.len()
}

fn write_u64_slice_le<W: Write + ?Sized>(writer: &mut W, words: &[u64]) -> std::io::Result<()> {
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

/// Original C global function `l2b_save_pac` from `minibwa/l2bit.c:286`.
pub(crate) fn l2b_save_pac<P: AsRef<Path>>(fn_: P, l2b: &L2b, both_strand: i32) -> i32 {
    let fn_ = c_path(fn_);
    let path = fn_.to_string_lossy();
    let mut fp: Box<dyn Write> = if path == "-" {
        Box::new(io::stdout())
    } else {
        match File::create(&fn_) {
            Ok(fp) => Box::new(fp),
            Err(_) => return -1,
        }
    };
    let n_pac = ((if both_strand != 0 {
        l2b.tot_len * 2
    } else {
        l2b.tot_len
    }) + 3)
        / 4;
    let mut pac = vec![0u8; n_pac as usize];
    let mut x = 0u64;
    for i in 0..l2b.tot_len {
        pac[(x >> 2) as usize] |= (l2b_get0(l2b, i) as u8) << ((!x & 3) * 2);
        x += 1;
    }
    if both_strand != 0 {
        for i in (0..l2b.tot_len).rev() {
            pac[(x >> 2) as usize] |= (3 - l2b_get0(l2b, i) as u8) << ((!x & 3) * 2);
            x += 1;
        }
    }
    let n_write = (x >> 2) + if (x & 3) == 0 { 0 } else { 1 };
    let _ = fp.write_all(&pac[..n_write as usize]);
    if x % 4 == 0 {
        let _ = fp.write_all(&[0]);
    }
    let _ = fp.write_all(&[(x % 4) as u8]);
    0
}

/// Original C static function `l2b_c2t` from `minibwa/l2bit.c:321`.
pub(crate) fn l2b_c2t(b: u8) -> u8 {
    if b == 1 {
        3
    } else {
        b
    }
}

/// Original C static function `l2b_g2a` from `minibwa/l2bit.c:322`.
pub(crate) fn l2b_g2a(b: u8) -> u8 {
    if b == 2 {
        0
    } else {
        b
    }
}

/// Original C global function `l2b_save_pac_meth` from `minibwa/l2bit.c:324`.
pub(crate) fn l2b_save_pac_meth<P: AsRef<Path>>(fn_: P, l2b: &L2b, both_strand: i32) -> i32 {
    let fn_ = c_path(fn_);
    let path = fn_.to_string_lossy();
    let mut fp: Box<dyn Write> = if path == "-" {
        Box::new(io::stdout())
    } else {
        match File::create(&fn_) {
            Ok(fp) => Box::new(fp),
            Err(_) => return -1,
        }
    };
    let mut len = l2b.tot_len * 2;
    if both_strand != 0 {
        len *= 2;
    }
    let n_pac = (len + 3) / 4;
    let mut pac = vec![0u8; n_pac as usize];
    let mut x = 0u64;
    for i in 0..l2b.tot_len {
        let b = l2b_c2t(l2b_get0(l2b, i) as u8);
        pac[(x >> 2) as usize] |= b << ((!x & 3) * 2);
        x += 1;
    }
    for i in 0..l2b.tot_len {
        let b = l2b_g2a(l2b_get0(l2b, i) as u8);
        pac[(x >> 2) as usize] |= b << ((!x & 3) * 2);
        x += 1;
    }
    if both_strand != 0 {
        for i in (0..l2b.tot_len).rev() {
            let b = 3 - l2b_g2a(l2b_get0(l2b, i) as u8);
            pac[(x >> 2) as usize] |= b << ((!x & 3) * 2);
            x += 1;
        }
        for i in (0..l2b.tot_len).rev() {
            let b = 3 - l2b_c2t(l2b_get0(l2b, i) as u8);
            pac[(x >> 2) as usize] |= b << ((!x & 3) * 2);
            x += 1;
        }
    }
    let n_write = (x >> 2) + if (x & 3) == 0 { 0 } else { 1 };
    let _ = fp.write_all(&pac[..n_write as usize]);
    if x % 4 == 0 {
        let _ = fp.write_all(&[0]);
    }
    let _ = fp.write_all(&[(x % 4) as u8]);
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pos2cid_returns_minus_one_for_positions_outside_contigs() {
        let l2b = L2b {
            n_ctg: 2,
            ctg: vec![
                L2bCtg {
                    name: "a".into(),
                    comm: None,
                    len: 10,
                    off: 0,
                },
                L2bCtg {
                    name: "b".into(),
                    comm: None,
                    len: 5,
                    off: 10,
                },
            ],
            ..Default::default()
        };
        let mut cst = -1;
        assert_eq!(l2b_pos2cid(&l2b, -1, 1, &mut cst), -1);
        assert_eq!(l2b_pos2cid(&l2b, 15, 1, &mut cst), -1);
        assert_eq!(l2b_pos2cid(&l2b, 8, 4, &mut cst), -1);
        assert_eq!(l2b_pos2cid(&l2b, 10, 5, &mut cst), 1);
        assert_eq!(cst, 0);
    }

    #[test]
    fn pos2cid_empty_contig_list_returns_minus_one() {
        let l2b = L2b::default();
        let mut cst = -1;
        assert_eq!(l2b_pos2cid(&l2b, 0, 1, &mut cst), -1);
        assert_eq!(cst, -1);
    }

    #[test]
    fn import_accepts_fastq_and_preserves_non_newline_sequence_bytes() {
        let path =
            std::env::temp_dir().join(format!("minibwa-rs-l2b-fastq-{}.fq", std::process::id()));
        std::fs::write(&path, b"@ctg comment\nAC GT\n+\n!!!!!\n").expect("write FASTQ");
        let l2b = l2b_import(&path, 11).expect("import FASTQ");
        let _ = std::fs::remove_file(&path);

        assert_eq!(l2b.n_ctg, 1);
        assert_eq!(l2b.ctg[0].name, "ctg");
        assert_eq!(l2b.ctg[0].comm.as_deref(), Some("comment"));
        assert_eq!(l2b.ctg[0].len, 5);
        assert_eq!(l2b.cat_name, b"ctg\0");
        assert_eq!(l2b.cat_comm, b"comment\0");

        let mut seq = vec![0; 5];
        assert_eq!(l2b_getseq(&l2b, 0, 0, 5, &mut seq), 5);
        assert_eq!(seq[0], 0);
        assert_eq!(seq[1], 1);
        assert_eq!(seq[3], 2);
        assert_eq!(seq[4], 3);
        assert_eq!(seq[2], 4);
    }

    #[test]
    fn import_metadata_uses_c_isspace_and_nul_truncation() {
        let path =
            std::env::temp_dir().join(format!("minibwa-rs-l2b-nul-meta-{}", std::process::id()));
        std::fs::write(&path, b">ctg\x01raw\0hidden comment\0tail\nACGT\n").expect("write FASTA");

        let l2b = l2b_import(&path, 3).expect("import raw metadata FASTA");
        let _ = std::fs::remove_file(&path);

        assert_eq!(l2b.n_ctg, 1);
        assert_eq!(l2b.cat_name, b"ctg\x01raw\0");
        assert_eq!(l2b.cat_comm, b"comment\0");
        assert_eq!(l2b.ctg[0].name, "ctg\u{1}raw");
        assert_eq!(l2b.ctg[0].comm.as_deref(), Some("comment"));
    }

    #[test]
    fn save_pac_meth_emits_expected_small_packed_sequence() {
        let mut l2b = L2b::default();
        let mut seq = b"ACGT".to_vec();
        let mut rng = 3;
        l2b_format_seq(seq.len() as u64, &mut seq, &mut rng);
        l2b_add_seq(&mut l2b, seq.len() as u64, &seq, "s", None, &mut rng);
        let path = std::env::temp_dir().join("minibwa-rs-small-meth.pac");
        assert_eq!(l2b_save_pac_meth(&path, &l2b, 1), 0);
        let bytes = std::fs::read(&path).expect("read meth pac");
        let _ = std::fs::remove_file(&path);
        assert_eq!(bytes, vec![0x3b, 0x13, 0x3b, 0x13, 0, 0]);
    }
}
