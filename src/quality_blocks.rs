use crate::{Error, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BlockMode {
    Raw = 0,
    DensePacked = 1,
}

struct Block {
    mode: BlockMode,
    width: u32,
    records: u32,
    alphabet: Vec<u8>,
    body: Vec<u8>,
}

pub(crate) struct QualityBlocks {
    #[cfg(test)]
    blocks: Vec<Block>,
    #[cfg(test)]
    block_records: u32,
}

impl QualityBlocks {
    #[cfg(test)]
    pub(crate) fn new(block_records: usize) -> Self {
        Self {
            blocks: Vec::new(),
            block_records: block_records.max(1) as u32,
        }
    }

    #[cfg(test)]
    pub(crate) fn push(&mut self, qualities: &[Vec<u8>]) -> Result<()> {
        let mut index = 0;
        while index < qualities.len() {
            let width = qualities[index].len();
            let end = qualities[index..]
                .iter()
                .take(self.block_records as usize)
                .position(|quality| quality.len() != width)
                .map_or_else(
                    || (index + self.block_records as usize).min(qualities.len()),
                    |offset| index + offset,
                );
            let block = encode_block(&qualities[index..end])?;
            self.blocks.push(block);
            index = end;
        }
        Ok(())
    }

    #[cfg(test)]
    pub(crate) fn write_to<W: std::io::Write>(&self, payload: &mut W) -> Result<()> {
        payload.write_all(&self.block_records.to_le_bytes())?;
        payload.write_all(&(self.blocks.len() as u32).to_le_bytes())?;
        for block in &self.blocks {
            write_block(payload, block)?;
        }
        Ok(())
    }

    pub(crate) fn decode(payload: &[u8]) -> Result<Vec<Vec<u8>>> {
        let mut cursor = 0usize;
        let _block_records = read_u32(payload, &mut cursor)? as usize;
        let block_count = read_u32(payload, &mut cursor)? as usize;
        let mut qualities = Vec::new();
        for _ in 0..block_count {
            let mode = match read_u8(payload, &mut cursor)? {
                0 => BlockMode::Raw,
                1 => BlockMode::DensePacked,
                other => {
                    return Err(Error::InvalidFormat(format!(
                        "unknown quality block mode {other}"
                    )))
                }
            };
            let width = read_u32(payload, &mut cursor)?;
            let records = read_u32(payload, &mut cursor)? as usize;
            let alphabet_len = read_u16(payload, &mut cursor)? as usize;
            let alphabet = read_bytes(payload, &mut cursor, alphabet_len)?.to_vec();
            let body_len = usize::try_from(read_u64(payload, &mut cursor)?)
                .map_err(|_| Error::InvalidFormat("quality block body too large".into()))?;
            let body = read_bytes(payload, &mut cursor, body_len)?.to_vec();
            let decoded = decode_block(mode, width, records, &alphabet, &body)?;
            qualities.extend(decoded);
        }
        if cursor != payload.len() {
            return Err(Error::InvalidFormat(
                "trailing data after quality blocks".into(),
            ));
        }
        Ok(qualities)
    }
}

#[cfg(test)]
pub(crate) struct QualityBlockWriter<W> {
    output: W,
    block_records: usize,
    block_count_position: u64,
    block_count: u32,
}

#[cfg(test)]
impl<W: std::io::Write + std::io::Seek> QualityBlockWriter<W> {
    pub(crate) fn new(mut output: W, block_records: usize) -> Result<Self> {
        let block_records = block_records.max(1);
        let encoded_block_records = u32::try_from(block_records)
            .map_err(|_| Error::InvalidInput("quality block size exceeds u32".into()))?;
        output.write_all(&encoded_block_records.to_le_bytes())?;
        let block_count_position = output.stream_position()?;
        output.write_all(&0u32.to_le_bytes())?;
        Ok(Self {
            output,
            block_records,
            block_count_position,
            block_count: 0,
        })
    }

    pub(crate) fn push(&mut self, qualities: &[Vec<u8>]) -> Result<()> {
        let mut index = 0;
        while index < qualities.len() {
            let width = qualities[index].len();
            let end = qualities[index..]
                .iter()
                .take(self.block_records)
                .position(|quality| quality.len() != width)
                .map_or_else(
                    || (index + self.block_records).min(qualities.len()),
                    |offset| index + offset,
                );
            write_block(&mut self.output, &encode_block(&qualities[index..end])?)?;
            self.block_count = self
                .block_count
                .checked_add(1)
                .ok_or_else(|| Error::InvalidInput("too many quality blocks".into()))?;
            index = end;
        }
        Ok(())
    }

    pub(crate) fn finish(mut self) -> Result<W> {
        let end = self.output.stream_position()?;
        self.output
            .seek(std::io::SeekFrom::Start(self.block_count_position))?;
        self.output.write_all(&self.block_count.to_le_bytes())?;
        self.output.seek(std::io::SeekFrom::Start(end))?;
        Ok(self.output)
    }
}

fn write_block<W: std::io::Write>(output: &mut W, block: &Block) -> Result<()> {
    output.write_all(&[block.mode as u8])?;
    output.write_all(&block.width.to_le_bytes())?;
    output.write_all(&block.records.to_le_bytes())?;
    output.write_all(&(block.alphabet.len() as u16).to_le_bytes())?;
    output.write_all(&block.alphabet)?;
    output.write_all(&(block.body.len() as u64).to_le_bytes())?;
    output.write_all(&block.body)?;
    Ok(())
}

fn read_u8(buf: &[u8], cursor: &mut usize) -> Result<u8> {
    if *cursor + 1 > buf.len() {
        return Err(Error::InvalidFormat(
            "truncated u8 in quality blocks".into(),
        ));
    }
    let v = buf[*cursor];
    *cursor += 1;
    Ok(v)
}

fn read_u16(buf: &[u8], cursor: &mut usize) -> Result<u16> {
    if *cursor + 2 > buf.len() {
        return Err(Error::InvalidFormat(
            "truncated u16 in quality blocks".into(),
        ));
    }
    let v = u16::from_le_bytes(buf[*cursor..*cursor + 2].try_into().unwrap());
    *cursor += 2;
    Ok(v)
}

fn read_u32(buf: &[u8], cursor: &mut usize) -> Result<u32> {
    if *cursor + 4 > buf.len() {
        return Err(Error::InvalidFormat(
            "truncated u32 in quality blocks".into(),
        ));
    }
    let v = u32::from_le_bytes(buf[*cursor..*cursor + 4].try_into().unwrap());
    *cursor += 4;
    Ok(v)
}

fn read_u64(buf: &[u8], cursor: &mut usize) -> Result<u64> {
    if *cursor + 8 > buf.len() {
        return Err(Error::InvalidFormat(
            "truncated u64 in quality blocks".into(),
        ));
    }
    let v = u64::from_le_bytes(buf[*cursor..*cursor + 8].try_into().unwrap());
    *cursor += 8;
    Ok(v)
}

fn read_bytes<'a>(buf: &'a [u8], cursor: &mut usize, length: usize) -> Result<&'a [u8]> {
    if buf.len().saturating_sub(*cursor) < length {
        return Err(Error::InvalidFormat(
            "truncated bytes in quality blocks".into(),
        ));
    }
    let value = &buf[*cursor..*cursor + length];
    *cursor += length;
    Ok(value)
}

fn encode_block(block: &[Vec<u8>]) -> Result<Block> {
    let width = block.first().map(|q| q.len()).unwrap_or(0);
    if width == 0 || block.is_empty() {
        return Ok(Block {
            mode: BlockMode::Raw,
            width: width as u32,
            records: block.len() as u32,
            alphabet: Vec::new(),
            body: Vec::new(),
        });
    }
    for q in block {
        if q.len() != width {
            return Err(Error::InvalidInput(
                "variable-length quality values are not supported".into(),
            ));
        }
    }

    let mut counts = [0usize; 256];
    for q in block {
        for &b in q {
            counts[b as usize] += 1;
        }
    }
    let mut alphabet: Vec<u8> = Vec::new();
    for byte in 0u16..=255 {
        let byte = byte as u8;
        if counts[byte as usize] > 0 {
            alphabet.push(byte);
        }
    }

    let raw: Vec<u8> = {
        let mut buf = Vec::with_capacity(block.len() * width);
        for q in block {
            buf.extend_from_slice(q);
        }
        buf
    };

    let mut symbol_map = [0u8; 256];
    for (i, &byte) in alphabet.iter().enumerate() {
        symbol_map[byte as usize] = i as u8;
    }

    let total = raw.len();
    let mut chosen_mode = BlockMode::Raw;
    let mut chosen_alphabet: Vec<u8> = Vec::new();
    let mut chosen_body = raw.clone();

    if alphabet.len() >= 2 && alphabet.len() <= 4 {
        let bits = alphabet.len().next_power_of_two().trailing_zeros() as usize;
        let mask = (1u8 << bits) - 1;
        let group = 8 / bits;
        let mut symbols: Vec<u8> = Vec::with_capacity(total);
        for &byte in &raw {
            symbols.push(symbol_map[byte as usize]);
        }
        let mut packed: Vec<u8> = Vec::new();
        for chunk in symbols.chunks(group) {
            let mut byte = 0u8;
            for (offset, &value) in chunk.iter().enumerate() {
                byte |= (value & mask) << (offset * bits);
            }
            packed.push(byte);
        }
        if packed.len() < chosen_body.len() {
            chosen_mode = BlockMode::DensePacked;
            chosen_alphabet = alphabet.clone();
            chosen_body = packed;
        }
    }

    Ok(Block {
        mode: chosen_mode,
        width: width as u32,
        records: block.len() as u32,
        alphabet: chosen_alphabet,
        body: chosen_body,
    })
}

pub(crate) fn encode_block_bytes(block: &[Vec<u8>]) -> Result<Vec<u8>> {
    let block = encode_block(block)?;
    let mut bytes = Vec::new();
    write_block(&mut bytes, &block)?;
    Ok(bytes)
}

fn decode_block(
    mode: BlockMode,
    width: u32,
    records: usize,
    alphabet: &[u8],
    body: &[u8],
) -> Result<Vec<Vec<u8>>> {
    let expected = (width as usize)
        .checked_mul(records)
        .ok_or_else(|| Error::InvalidFormat("quality block size overflow".into()))?;
    let raw = match mode {
        BlockMode::Raw => {
            if !alphabet.is_empty() || body.len() != expected {
                return Err(Error::InvalidFormat(
                    "quality block body size mismatch".into(),
                ));
            }
            body.to_vec()
        }
        BlockMode::DensePacked => {
            if !(2..=4).contains(&alphabet.len()) {
                return Err(Error::InvalidFormat(
                    "invalid packed quality alphabet".into(),
                ));
            }
            let mut raw = unpack_packed(body, alphabet)?;
            if raw.len() < expected {
                return Err(Error::InvalidFormat(
                    "quality block body size mismatch".into(),
                ));
            }
            raw.truncate(expected);
            raw
        }
    };
    let mut out = Vec::with_capacity(records);
    for row in 0..records {
        let start = row * width as usize;
        out.push(raw[start..start + width as usize].to_vec());
    }
    Ok(out)
}

fn unpack_packed(body: &[u8], alphabet: &[u8]) -> Result<Vec<u8>> {
    if alphabet.is_empty() {
        return Ok(Vec::new());
    }
    let bits = alphabet.len().next_power_of_two().trailing_zeros() as usize;
    let mask = (1u8 << bits) - 1;
    let group = 8 / bits;
    let mut out = Vec::with_capacity(body.len() * group);
    for &byte in body {
        for offset in 0..group {
            let value = (byte >> (offset * bits)) & mask;
            let symbol = alphabet
                .get(value as usize)
                .copied()
                .ok_or_else(|| Error::InvalidFormat("packed symbol out of range".into()))?;
            out.push(symbol);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_block_raw() {
        let block: Vec<Vec<u8>> = vec![
            b"DDDDDDDD".to_vec(),
            b"DDDD9DDD".to_vec(),
            b"DD,DD9DD".to_vec(),
        ];
        let encoded = encode_block(&block).unwrap();
        assert_eq!(encoded.mode, BlockMode::DensePacked);
        let decoded = decode_block(
            encoded.mode,
            encoded.width,
            encoded.records as usize,
            &encoded.alphabet,
            &encoded.body,
        )
        .unwrap();
        assert_eq!(decoded, block);
    }

    #[test]
    fn round_trips_quality_blocks_with_packing() {
        let qualities: Vec<Vec<u8>> = (0..1000)
            .map(|i| {
                let mut q = vec![b'D'; 150];
                for (idx, byte) in q.iter_mut().enumerate() {
                    let value = match (i + idx) % 17 {
                        0 => b'!',
                        1 => b',',
                        2..=5 => b'9',
                        _ => b'D',
                    };
                    *byte = value;
                }
                q
            })
            .collect();
        let mut blocks = QualityBlocks::new(256);
        blocks.push(&qualities).unwrap();
        assert!(blocks
            .blocks
            .iter()
            .any(|b| b.mode == BlockMode::DensePacked));
        let mut buf = Vec::new();
        blocks.write_to(&mut buf).unwrap();
        let decoded = QualityBlocks::decode(&buf).unwrap();
        assert_eq!(decoded, qualities);
    }

    #[test]
    fn falls_back_to_raw_for_large_alphabet() {
        let block: Vec<Vec<u8>> = (0..16)
            .map(|i| (0..32).map(|j| (i + j) as u8).collect())
            .collect();
        let encoded = encode_block(&block).unwrap();
        assert_eq!(encoded.mode, BlockMode::Raw);
        let decoded = decode_block(
            encoded.mode,
            encoded.width,
            encoded.records as usize,
            &encoded.alphabet,
            &encoded.body,
        )
        .unwrap();
        assert_eq!(decoded, block);
    }

    #[test]
    fn round_trips_partial_packed_byte_and_variable_widths() {
        let qualities = vec![
            b"DD9".to_vec(),
            b"9DD".to_vec(),
            b"D9".to_vec(),
            b"99".to_vec(),
        ];
        let mut blocks = QualityBlocks::new(256);
        blocks.push(&qualities).unwrap();
        assert_eq!(blocks.blocks.len(), 2);
        assert!(blocks
            .blocks
            .iter()
            .all(|block| block.mode == BlockMode::DensePacked));

        let mut encoded = Vec::new();
        blocks.write_to(&mut encoded).unwrap();
        assert_eq!(QualityBlocks::decode(&encoded).unwrap(), qualities);
    }

    #[test]
    fn streaming_writer_matches_buffered_encoding() {
        let qualities = vec![
            b"DD9".to_vec(),
            b"9DD".to_vec(),
            b"D9".to_vec(),
            b"99".to_vec(),
            b"!,9D".to_vec(),
        ];
        let mut buffered = QualityBlocks::new(2);
        buffered.push(&qualities).unwrap();
        let mut expected = Vec::new();
        buffered.write_to(&mut expected).unwrap();

        let cursor = std::io::Cursor::new(Vec::new());
        let mut streaming = QualityBlockWriter::new(cursor, 2).unwrap();
        streaming.push(&qualities).unwrap();
        let actual = streaming.finish().unwrap().into_inner();

        assert_eq!(actual, expected);
        assert_eq!(QualityBlocks::decode(&actual).unwrap(), qualities);
    }
}
