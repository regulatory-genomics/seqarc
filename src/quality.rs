use crate::{Error, Result};
use std::io::Read;

const HEADER_LEN: usize = 8;

pub(crate) fn encode_run_bytes(rows: &[Vec<u8>]) -> Result<Vec<u8>> {
    let first = rows
        .first()
        .ok_or_else(|| Error::InvalidInput("quality run cannot be empty".into()))?;
    let width = first.len();
    let width = u32::try_from(width)
        .map_err(|_| Error::InvalidInput("quality row width exceeds u32".into()))?;
    let records = u32::try_from(rows.len())
        .map_err(|_| Error::InvalidInput("quality row count exceeds u32".into()))?;
    let body_len = (width as usize)
        .checked_mul(rows.len())
        .ok_or_else(|| Error::InvalidInput("quality run body size overflow".into()))?;

    for row in rows {
        if row.len() != width as usize {
            return Err(Error::InvalidInput(
                "quality rows in a run have different widths".into(),
            ));
        }
    }

    let mut encoded = Vec::with_capacity(
        HEADER_LEN
            .checked_add(body_len)
            .ok_or_else(|| Error::InvalidInput("quality run encoded size overflow".into()))?,
    );
    encoded.extend_from_slice(&width.to_le_bytes());
    encoded.extend_from_slice(&records.to_le_bytes());
    for row in rows {
        encoded.extend_from_slice(row);
    }
    Ok(encoded)
}

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

pub(crate) struct QualityStreamReader<R> {
    input: R,
    expected_records: usize,
    records_read: usize,
}

impl<R: Read> QualityStreamReader<R> {
    pub(crate) fn new(input: R, expected_records: usize) -> Self {
        Self {
            input,
            expected_records,
            records_read: 0,
        }
    }

    pub(crate) fn next(&mut self) -> Result<Option<Vec<Vec<u8>>>> {
        if self.records_read == self.expected_records {
            return Ok(None);
        }

        let width = read_u32_reader(&mut self.input, "quality run width")? as usize;
        let records = read_u32_reader(&mut self.input, "quality run record count")? as usize;
        if records == 0 {
            return Err(Error::InvalidFormat("quality run has zero records".into()));
        }
        let remaining = self.expected_records - self.records_read;
        if records > remaining {
            return Err(Error::InvalidFormat(
                "quality run exceeds expected record count".into(),
            ));
        }
        let body_len = width
            .checked_mul(records)
            .ok_or_else(|| Error::InvalidFormat("quality run body size overflow".into()))?;
        let mut body = vec![0; body_len];
        read_exact(&mut self.input, &mut body, "quality run body")?;
        self.records_read += records;
        let rows = if width == 0 {
            (0..records).map(|_| Vec::new()).collect()
        } else {
            body.chunks_exact(width).map(|row| row.to_vec()).collect()
        };
        Ok(Some(rows))
    }

    pub(crate) fn finish(self) -> Result<R> {
        if self.records_read != self.expected_records {
            return Err(Error::InvalidFormat(
                "quality stream record count mismatch".into(),
            ));
        }
        Ok(self.input)
    }
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

fn read_u32_reader<R: Read>(input: &mut R, what: &str) -> Result<u32> {
    let mut bytes = [0; 4];
    read_exact(input, &mut bytes, what)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_exact<R: Read>(input: &mut R, bytes: &mut [u8], what: &str) -> Result<()> {
    input.read_exact(bytes).map_err(|error| match error.kind() {
        std::io::ErrorKind::UnexpectedEof => Error::InvalidFormat(format!("truncated {what}")),
        _ => Error::Io(error),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

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
    fn encodes_one_run_and_rejects_invalid_input() {
        assert_eq!(
            encode_run_bytes(&[b"A".to_vec(), b"B".to_vec()]).unwrap(),
            [
                1u32.to_le_bytes().as_slice(),
                2u32.to_le_bytes().as_slice(),
                b"AB"
            ]
            .concat()
        );
        assert!(encode_run_bytes(&[]).is_err());
        assert!(encode_run_bytes(&[b"A".to_vec(), b"BC".to_vec()]).is_err());
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

    #[test]
    fn streaming_reader_preserves_following_bytes() {
        let encoded = encode_stream_bytes(&[b"AA".as_slice(), b"B".as_slice()]).unwrap();
        let mut input = encoded.clone();
        input.extend_from_slice(b"following");
        let mut reader = QualityStreamReader::new(Cursor::new(input), 2);
        assert_eq!(reader.next().unwrap().unwrap(), vec![b"AA".to_vec()]);
        assert_eq!(reader.next().unwrap().unwrap(), vec![b"B".to_vec()]);
        assert!(reader.next().unwrap().is_none());
        let input = reader.finish().unwrap().into_inner();
        assert_eq!(&input[encoded.len()..], b"following");
    }
}
