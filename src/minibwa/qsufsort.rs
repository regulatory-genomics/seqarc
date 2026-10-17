#![allow(unused_variables)]

pub(crate) const INSERT_SORT_NUM_ITEM: i64 = 16;

/// Makes suffix array p of x. x becomes inverse of p. p and x are both of size
/// n+1. Contents of x[0..n-1] are integers in the range l..k-1. Original
/// contents of x[n] is disregarded, the n-th symbol being regarded as
/// end-of-string smaller than all other symbols.
///
/// Original C global function `suffix_sort` from `minibwa/qsufsort.c:56`.
pub(crate) fn suffix_sort(
    v: &mut [i64],
    idx: &mut [i64],
    num_char: i64,
    largest_input_symbol: i64,
    smallest_input_symbol: i64,
    skip_transform: i32,
) {
    let mut num_sorted_pos = 1i64;
    if skip_transform == 0 {
        let mut num_symbol_aggregated = 0;
        let new_alphabet_size = transform(
            v,
            idx,
            num_char,
            largest_input_symbol,
            smallest_input_symbol,
            num_char,
            &mut num_symbol_aggregated,
        );
        bucket_sort(v, idx, num_char, new_alphabet_size);
        idx[0] = -1;
        v[num_char as usize] = 0;
        num_sorted_pos = num_symbol_aggregated;
    }

    while idx[0] >= -num_char {
        let mut i = 0i64;
        let mut negated_sorted_group_length = 0i64;
        loop {
            let s = idx[i as usize];
            if s < 0 {
                i -= s;
                negated_sorted_group_length += s;
            } else {
                if negated_sorted_group_length != 0 {
                    idx[(i + negated_sorted_group_length) as usize] = negated_sorted_group_length;
                    negated_sorted_group_length = 0;
                }
                let j = v[s as usize] + 1;
                sort_split(v, idx, i, j - 1, num_sorted_pos);
                i = j;
            }
            if i > num_char {
                break;
            }
        }
        if negated_sorted_group_length != 0 {
            idx[(i + negated_sorted_group_length) as usize] = negated_sorted_group_length;
        }
        num_sorted_pos *= 2;
    }
}

/// Sorting routine called for each unsorted group. Sorts the array of integers
/// (suffix numbers) of length n starting at p. The algorithm is a ternary-split
/// quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
/// Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
/// function is based on Program 7.
///
/// Original C static function `sort_split` from `minibwa/qsufsort.c:113`.
pub(crate) fn sort_split(
    v_arr: &mut [i64],
    idx: &mut [i64],
    lowest_pos: i64,
    highest_pos: i64,
    num_sorted_char: i64,
) {
    let num_item = highest_pos - lowest_pos + 1;
    if num_item <= INSERT_SORT_NUM_ITEM {
        insert_sort_split(v_arr, idx, lowest_pos, highest_pos, num_sorted_char);
        return;
    }

    let v = choose_pivot(v_arr, idx, lowest_pos, highest_pos, num_sorted_char);
    let mut a = lowest_pos;
    let mut b = lowest_pos;
    let mut c = highest_pos;
    let mut d = highest_pos;

    loop {
        while c >= b && v_arr[(idx[b as usize] + num_sorted_char) as usize] <= v {
            if v_arr[(idx[b as usize] + num_sorted_char) as usize] == v {
                idx.swap(a as usize, b as usize);
                a += 1;
            }
            b += 1;
        }
        while c >= b && v_arr[(idx[c as usize] + num_sorted_char) as usize] >= v {
            if v_arr[(idx[c as usize] + num_sorted_char) as usize] == v {
                idx.swap(c as usize, d as usize);
                d -= 1;
            }
            c -= 1;
        }
        if b > c {
            break;
        }
        idx.swap(b as usize, c as usize);
        b += 1;
        c -= 1;
    }

    let mut s = a - lowest_pos;
    let mut t = b - a;
    s = s.min(t);
    let mut l = lowest_pos;
    let mut m = b - s;
    while m < b {
        idx.swap(l as usize, m as usize);
        l += 1;
        m += 1;
    }

    s = d - c;
    t = highest_pos - d;
    s = s.min(t);
    l = b;
    m = highest_pos - s + 1;
    while m <= highest_pos {
        idx.swap(l as usize, m as usize);
        l += 1;
        m += 1;
    }

    s = b - a;
    t = d - c;
    if s > 0 {
        sort_split(v_arr, idx, lowest_pos, lowest_pos + s - 1, num_sorted_char);
    }

    a = lowest_pos + s;
    b = highest_pos - t;
    if a == b {
        v_arr[idx[a as usize] as usize] = a;
        idx[a as usize] = -1;
    } else {
        c = a;
        while c <= b {
            v_arr[idx[c as usize] as usize] = b;
            c += 1;
        }
    }

    if t > 0 {
        sort_split(v_arr, idx, highest_pos - t + 1, highest_pos, num_sorted_char);
    }
}

/// Algorithm by Bentley & McIlroy.
///
/// Original C static function `choose_pivot` from `minibwa/qsufsort.c:194`.
pub(crate) fn choose_pivot(
    v: &[i64],
    idx: &[i64],
    lowest_pos: i64,
    highest_pos: i64,
    num_sorted_char: i64,
) -> i64 {
    let num_item = highest_pos - lowest_pos + 1;
    let m = lowest_pos + num_item / 2;
    let s = num_item / 8;
    let key1 = v[(idx[lowest_pos as usize] + num_sorted_char) as usize];
    let key2 = v[(idx[(lowest_pos + s) as usize] + num_sorted_char) as usize];
    let key3 = v[(idx[(lowest_pos + 2 * s) as usize] + num_sorted_char) as usize];
    let keyl = if key1 < key2 {
        if key2 < key3 {
            key2
        } else if key1 < key3 {
            key3
        } else {
            key1
        }
    } else if key2 > key3 {
        key2
    } else if key1 > key3 {
        key3
    } else {
        key1
    };
    let key1 = v[(idx[(m - s) as usize] + num_sorted_char) as usize];
    let key2 = v[(idx[m as usize] + num_sorted_char) as usize];
    let key3 = v[(idx[(m + s) as usize] + num_sorted_char) as usize];
    let keym = if key1 < key2 {
        if key2 < key3 {
            key2
        } else if key1 < key3 {
            key3
        } else {
            key1
        }
    } else if key2 > key3 {
        key2
    } else if key1 > key3 {
        key3
    } else {
        key1
    };
    let key1 = v[(idx[(highest_pos - 2 * s) as usize] + num_sorted_char) as usize];
    let key2 = v[(idx[(highest_pos - s) as usize] + num_sorted_char) as usize];
    let key3 = v[(idx[highest_pos as usize] + num_sorted_char) as usize];
    let keyn = if key1 < key2 {
        if key2 < key3 {
            key2
        } else if key1 < key3 {
            key3
        } else {
            key1
        }
    } else if key2 > key3 {
        key2
    } else if key1 > key3 {
        key3
    } else {
        key1
    };
    if keyl < keym {
        if keym < keyn {
            keym
        } else if keyl < keyn {
            keyn
        } else {
            keyl
        }
    } else if keym > keyn {
        keym
    } else if keyl > keyn {
        keyn
    } else {
        keyl
    }
}

/// Quadratic sorting method to use for small subarrays.
///
/// Original C static function `insert_sort_split` from `minibwa/qsufsort.c:227`.
pub(crate) fn insert_sort_split(
    v: &mut [i64],
    idx: &mut [i64],
    lowest_pos: i64,
    highest_pos: i64,
    num_sorted_char: i64,
) {
    let num_item = highest_pos - lowest_pos + 1;
    let mut key = vec![0i64; num_item as usize];
    let mut pos = vec![0i64; num_item as usize];
    for i in 0..num_item as usize {
        pos[i] = idx[lowest_pos as usize + i];
        key[i] = v[(pos[i] + num_sorted_char) as usize];
    }
    for i in 1..num_item as usize {
        let tmp_key = key[i];
        let tmp_pos = pos[i];
        let mut j = i;
        while j > 0 && key[j - 1] > tmp_key {
            key[j] = key[j - 1];
            pos[j] = pos[j - 1];
            j -= 1;
        }
        key[j] = tmp_key;
        pos[j] = tmp_pos;
    }

    let mut negative_sorted_length = -1i64;
    let mut i = num_item - 1;
    let mut group_num = highest_pos;
    while i > 0 {
        idx[(i + lowest_pos) as usize] = pos[i as usize];
        v[idx[(i + lowest_pos) as usize] as usize] = group_num;
        if key[(i - 1) as usize] == key[i as usize] {
            negative_sorted_length = 0;
        } else {
            if negative_sorted_length < 0 {
                idx[(i + lowest_pos) as usize] = negative_sorted_length;
            }
            group_num = i + lowest_pos - 1;
            negative_sorted_length -= 1;
        }
        i -= 1;
    }
    idx[lowest_pos as usize] = pos[0];
    v[idx[lowest_pos as usize] as usize] = group_num;
    if negative_sorted_length < 0 {
        idx[lowest_pos as usize] = negative_sorted_length;
    }
}

/// Bucketsort for first iteration.
///
/// Input: x[0..n-1] holds integers in the range 1..k-1, all of which appear
/// at least once. x[n] is 0. (This is the corresponding output of transform.) k
/// must be at most n+1. p is array of size n+1 whose contents are disregarded.
///
/// Output: x is V and p is I after the initial sorting stage of the refined
/// suffix sorting algorithm.
///
/// Original C static function `bucket_sort` from `minibwa/qsufsort.c:288`.
pub(crate) fn bucket_sort(v: &mut [i64], idx: &mut [i64], num_char: i64, alphabet_size: i64) {
    for i in 0..alphabet_size as usize {
        idx[i] = -1;
    }
    for i in 0..=num_char as usize {
        let c = v[i] as usize;
        v[i] = idx[c];
        idx[c] = i as i64;
    }
    let mut current_index = num_char;
    let mut i = alphabet_size;
    while i > 0 {
        let mut c = idx[(i - 1) as usize];
        let mut d = v[c as usize];
        let group_num = current_index;
        v[c as usize] = group_num;
        if d >= 0 {
            idx[current_index as usize] = c;
            while d >= 0 {
                c = d;
                d = v[c as usize];
                v[c as usize] = group_num;
                current_index -= 1;
                idx[current_index as usize] = c;
            }
        } else {
            idx[current_index as usize] = -1;
        }
        current_index -= 1;
        i -= 1;
    }
}

/// Transforms the alphabet of x by attempting to aggregate several symbols into
/// one, while preserving the suffix order of x. The alphabet may also be
/// compacted, so that x on output comprises all integers of the new alphabet
/// with no skipped numbers.
///
/// Input: x is an array of size n+1 whose first n elements are positive
/// integers in the range l..k-1. p is array of size n+1, used for temporary
/// storage. q controls aggregation and compaction by defining the maximum intue
/// for any symbol during transformation: q must be at least k-l; if q<=n,
/// compaction is guaranteed; if k-l>n, compaction is never done; if q is
/// INT_MAX, the maximum number of symbols are aggregated into one.
///
/// Output: Returns an integer j in the range 1..q representing the size of the
/// new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
/// set to the number of old symbols grouped into one. Only x[n] is 0.
///
/// Original C static function `transform` from `minibwa/qsufsort.c:344`.
pub(crate) fn transform(
    v: &mut [i64],
    idx: &mut [i64],
    num_char: i64,
    largest_input_symbol: i64,
    smallest_input_symbol: i64,
    max_new_alphabet_size: i64,
    num_symbol_aggregated: &mut i64,
) -> i64 {
    let max_num_input_symbol = largest_input_symbol - smallest_input_symbol + 1;
    let mut max_num_bit = 0i64;
    let mut i = max_num_input_symbol;
    while i != 0 {
        max_num_bit += 1;
        i >>= 1;
    }
    let max_symbol = i64::MAX >> max_num_bit;
    let mut c = max_num_input_symbol;
    let mut a = 0i64;
    let mut min_symbol_in_chunk = 0i64;
    let mut max_symbol_in_chunk = 0i64;
    while a < num_char && max_symbol_in_chunk <= max_symbol && c <= max_new_alphabet_size {
        min_symbol_in_chunk =
            (min_symbol_in_chunk << max_num_bit) | (v[a as usize] - smallest_input_symbol + 1);
        max_symbol_in_chunk = c;
        c = (max_symbol_in_chunk << max_num_bit) | max_num_input_symbol;
        a += 1;
    }

    let mask = (1i64 << ((a - 1) * max_num_bit)) - 1;
    v[num_char as usize] = smallest_input_symbol - 1;
    for i in 0..=max_symbol_in_chunk as usize {
        idx[i] = 0;
    }
    c = min_symbol_in_chunk;
    i = a;
    while i <= num_char {
        idx[c as usize] = 1;
        c = ((c & mask) << max_num_bit) | (v[i as usize] - smallest_input_symbol + 1);
        i += 1;
    }
    i = 1;
    while i < a {
        idx[c as usize] = 1;
        c = (c & mask) << max_num_bit;
        i += 1;
    }
    let mut new_alphabet_size = 1i64;
    for i in 0..=max_symbol_in_chunk as usize {
        if idx[i] != 0 {
            idx[i] = new_alphabet_size;
            new_alphabet_size += 1;
        }
    }
    c = min_symbol_in_chunk;
    i = 0;
    let mut j = a;
    while j <= num_char {
        v[i as usize] = idx[c as usize];
        c = ((c & mask) << max_num_bit) | (v[j as usize] - smallest_input_symbol + 1);
        i += 1;
        j += 1;
    }
    while i < num_char {
        v[i as usize] = idx[c as usize];
        c = (c & mask) << max_num_bit;
        i += 1;
    }
    v[num_char as usize] = 0;
    *num_symbol_aggregated = a;
    new_alphabet_size
}

#[cfg(test)]
mod tests {
    

}
