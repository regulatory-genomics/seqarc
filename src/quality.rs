use crate::{Error, Result};

const HEADER_LEN: usize = 8;

pub(crate) fn encode_stream_bytes(rows: &[&[u8]]) -> Result<Vec<u8>> {
    let mut encoded = Vec::new();
    let mut start = 0;
    while start < rows.len() {
        let width = rows[start].len();
        let end = rows[start..]
            .iter()
            .position(|row| row.len() != width)
            .map_or(rows.len(), |offset| start + offset);
        let body_len = width
            .checked_mul(end - start)
            .ok_or_else(|| Error::InvalidInput("quality run body size overflow".into()))?;
        let run_len = HEADER_LEN
            .checked_add(body_len)
            .ok_or_else(|| Error::InvalidInput("quality run encoded size overflow".into()))?;
        encoded.reserve(run_len);
        encoded.extend_from_slice(
            &u32::try_from(width)
                .map_err(|_| Error::InvalidInput("quality row width exceeds u32".into()))?
                .to_le_bytes(),
        );
        encoded.extend_from_slice(
            &u32::try_from(end - start)
                .map_err(|_| Error::InvalidInput("quality row count exceeds u32".into()))?
                .to_le_bytes(),
        );
        for row in &rows[start..end] {
            encoded.extend_from_slice(row);
        }
        start = end;
    }
    Ok(encoded)
}

pub(crate) fn decode_stream_bytes(input: &[u8], expected_records: usize) -> Result<Vec<Vec<u8>>> {
    let mut cursor = 0;
    let mut decoded = Vec::with_capacity(expected_records);
    while decoded.len() < expected_records {
        let width = read_u32(input, &mut cursor, "quality run width")? as usize;
        let records = read_u32(input, &mut cursor, "quality run record count")? as usize;
        if records == 0 {
            return Err(Error::InvalidFormat("quality run has zero records".into()));
        }
        let remaining = expected_records - decoded.len();
        if records > remaining {
            return Err(Error::InvalidFormat(
                "quality run exceeds expected record count".into(),
            ));
        }
        let body_len = width
            .checked_mul(records)
            .ok_or_else(|| Error::InvalidFormat("quality run body size overflow".into()))?;
        let body = read_bytes(input, &mut cursor, body_len)?;
        if width == 0 {
            decoded.extend((0..records).map(|_| Vec::new()));
        } else {
            for row in body.chunks_exact(width) {
                decoded.push(row.to_vec());
            }
        }
    }
    if cursor != input.len() {
        return Err(Error::InvalidFormat(
            "trailing bytes after quality stream".into(),
        ));
    }
    Ok(decoded)
}

fn read_u32(input: &[u8], cursor: &mut usize, what: &str) -> Result<u32> {
    let bytes = read_bytes(input, cursor, 4)?;
    Ok(u32::from_le_bytes(bytes.try_into().map_err(|_| {
        Error::InvalidFormat(format!("truncated {what}"))
    })?))
}

fn read_bytes<'a>(input: &'a [u8], cursor: &mut usize, length: usize) -> Result<&'a [u8]> {
    let end = cursor
        .checked_add(length)
        .ok_or_else(|| Error::InvalidFormat("quality stream offset overflow".into()))?;
    if end > input.len() {
        return Err(Error::InvalidFormat("truncated quality stream".into()));
    }
    let bytes = &input[*cursor..end];
    *cursor = end;
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_and_decodes_empty_fixed_variable_and_zero_width_rows() {
        assert_eq!(encode_stream_bytes(&[]).unwrap(), Vec::<u8>::new());
        let rows = vec![
            b"AB".as_slice(),
            b"CD".as_slice(),
            b"E".as_slice(),
            b"".as_slice(),
        ];
        let encoded = encode_stream_bytes(&rows).unwrap();
        assert_eq!(decode_stream_bytes(&encoded, rows.len()).unwrap(), rows);
        assert_eq!(decode_stream_bytes(&[], 0).unwrap(), Vec::<Vec<u8>>::new());
    }

    #[test]
    fn rejects_truncation_zero_count_count_mismatch_overflow_and_trailing() {
        let valid = encode_stream_bytes(&[b"AB".as_slice()]).unwrap();
        for end in 0..valid.len() {
            assert!(decode_stream_bytes(&valid[..end], 1).is_err());
        }
        assert!(decode_stream_bytes(&[0; 8], 1).is_err());
        assert!(decode_stream_bytes(&valid, 2).is_err());
        let mut too_many = valid.clone();
        too_many[4..8].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(decode_stream_bytes(&too_many, 1).is_err());
        assert!(decode_stream_bytes(&[u32::MAX as u8; 8], 1).is_err());
        assert!(decode_stream_bytes(&[0; 1], 0).is_err());
        assert!(decode_stream_bytes(&[valid.as_slice(), &[0]].concat(), 1).is_err());
    }
}
