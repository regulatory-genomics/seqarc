#[cfg(not(unix))]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Timeval {
    pub tv_sec: i64,
    pub tv_usec: i64,
}


const fn make_kom_nt4_table() -> [u8; 256] {
    let mut table = [4u8; 256];
    table[0] = 0;
    table[1] = 1;
    table[2] = 2;
    table[3] = 3;
    table[b'A' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'a' as usize] = 0;
    table[b'c' as usize] = 1;
    table[b'g' as usize] = 2;
    table[b't' as usize] = 3;
    table
}

pub(crate) const KOM_NT4_TABLE: [u8; 256] = make_kom_nt4_table();


/// Original C global function `kom_parse_num` from `minibwa/kommon.c:19`.
pub(crate) fn kom_parse_num(str_: &str) -> (i64, usize) {
    let bytes = str_.as_bytes();
    let nul = bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len());
    let (mut x, mut p) = parse_c_float_prefix(&str_[..nul]);
    if p < bytes.len() {
        if bytes[p] == b'G' || bytes[p] == b'g' {
            x *= 1e9;
            p += 1;
        } else if bytes[p] == b'M' || bytes[p] == b'm' {
            x *= 1e6;
            p += 1;
        } else if bytes[p] == b'K' || bytes[p] == b'k' {
            x *= 1e3;
            p += 1;
        }
    }
    ((x + 0.499) as i64, p)
}

fn parse_c_float_prefix(s: &str) -> (f64, usize) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let bytes = s.as_bytes();
        let nul = bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len());
        let mut c_buf = Vec::with_capacity(nul + 1);
        c_buf.extend_from_slice(&bytes[..nul]);
        c_buf.push(0);
        let mut end: *mut std::ffi::c_char = std::ptr::null_mut();
        let x = unsafe { libc::strtod(c_buf.as_ptr() as *const std::ffi::c_char, &mut end) };
        let consumed = if end.is_null() {
            0
        } else {
            unsafe { end.offset_from(c_buf.as_ptr() as *const std::ffi::c_char) as usize }
        };
        (x, consumed)
    }
    #[cfg(target_arch = "wasm32")]
    {
        parse_decimal_float_prefix(s)
    }
}

#[cfg(target_arch = "wasm32")]
fn parse_decimal_float_prefix(s: &str) -> (f64, usize) {
    let bytes = s.as_bytes();
    let nul = bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len());
    let bytes = &bytes[..nul];
    let mut i = 0usize;
    while i < bytes.len() && bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    if i < bytes.len() && matches!(bytes[i], b'+' | b'-') {
        i += 1;
    }
    let mut saw_digit = false;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        saw_digit = true;
        i += 1;
    }
    if i < bytes.len() && bytes[i] == b'.' {
        i += 1;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            saw_digit = true;
            i += 1;
        }
    }
    if !saw_digit {
        return (0.0, 0);
    }
    if i < bytes.len() && matches!(bytes[i], b'e' | b'E') {
        let exp = i;
        i += 1;
        if i < bytes.len() && matches!(bytes[i], b'+' | b'-') {
            i += 1;
        }
        let exp_digits = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i == exp_digits {
            i = exp;
        }
    }
    let parsed = std::str::from_utf8(&bytes[start..i])
        .ok()
        .and_then(|x| x.parse::<f64>().ok())
        .unwrap_or(0.0);
    (parsed, i)
}

pub(crate) fn kom_strtol(str_: &str) -> (i64, usize) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let bytes = str_.as_bytes();
        let nul = bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len());
        let mut c_buf = Vec::with_capacity(nul + 1);
        c_buf.extend_from_slice(&bytes[..nul]);
        c_buf.push(0);
        let mut end: *mut std::ffi::c_char = std::ptr::null_mut();
        let x = unsafe { libc::strtol(c_buf.as_ptr() as *const std::ffi::c_char, &mut end, 10) };
        let consumed = if end.is_null() {
            0
        } else {
            unsafe { end.offset_from(c_buf.as_ptr() as *const std::ffi::c_char) as usize }
        };
        (x as i64, consumed)
    }
    #[cfg(target_arch = "wasm32")]
    {
        parse_decimal_int_prefix(str_)
    }
}

#[cfg(target_arch = "wasm32")]
fn parse_decimal_int_prefix(s: &str) -> (i64, usize) {
    let bytes = s.as_bytes();
    let nul = bytes.iter().position(|&c| c == 0).unwrap_or(bytes.len());
    let bytes = &bytes[..nul];
    let mut i = 0usize;
    while i < bytes.len() && bytes[i].is_ascii_whitespace() {
        i += 1;
    }
    let start = i;
    let neg = if i < bytes.len() && matches!(bytes[i], b'+' | b'-') {
        let neg = bytes[i] == b'-';
        i += 1;
        neg
    } else {
        false
    };
    let digits = i;
    let mut value = 0i64;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        value = value
            .saturating_mul(10)
            .saturating_add((bytes[i] - b'0') as i64);
        i += 1;
    }
    if i == digits {
        return (0, 0);
    }
    if neg {
        value = value.saturating_neg();
    }
    (value, i.max(start))
}

pub(crate) fn kom_atoi(str_: &str) -> i32 {
    kom_strtol(str_).0 as i32
}

pub(crate) fn kom_atol(str_: &str) -> i64 {
    kom_strtol(str_).0
}

/// Original C global function `kom_panic` from `minibwa/kommon.c:31`.
pub(crate) fn kom_panic(func: &str, msg: &str) -> ! {
    eprintln!("[E::{func}] {msg} ABORT!");
    std::process::abort()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_num_matches_suffix_rounding_rules() {
        assert_eq!(kom_parse_num("12.4k rest"), (12400, 5));
        assert_eq!(kom_parse_num("1.5M"), (1500000, 4));
        assert_eq!(kom_parse_num("2Gx"), (2000000000, 2));
        assert_eq!(kom_parse_num("17"), (17, 2));
        assert_eq!(kom_parse_num("  .5k"), (500, 5));
        assert_eq!(kom_parse_num("abc"), (0, 0));
    }

    #[test]
    fn integer_helpers_use_c_strtol_boundaries_and_overflow() {
        let neg = "-999999999999999999999999999999999999999tail";
        assert_eq!(kom_strtol("  -42xyz"), (-42, 5));
        assert_eq!(kom_atol("17\0hidden"), 17);
        assert_eq!(kom_strtol("abc"), (0, 0));
        assert_eq!(kom_atol(neg), libc::c_long::MIN as i64);
        assert_eq!(kom_atoi("4294967297"), 1);
    }

}
