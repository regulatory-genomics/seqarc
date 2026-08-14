use crate::{Error, Result};

pub const MIN_LEVEL: i32 = 1;
pub const MAX_LEVEL: i32 = 22;

pub struct Encoder {
    compressor: zstd::bulk::Compressor<'static>,
}

impl Encoder {
    /// Creates a single-threaded bulk encoder.
    ///
    /// Parallelism is provided by the archive pipeline, which runs one encoder
    /// per worker thread across independent chunks. zstd's internal
    /// multithreading is ineffective here because individual sections are
    /// smaller than a single compression job at high levels.
    pub fn new(level: i32) -> Result<Self> {
        if !(MIN_LEVEL..=MAX_LEVEL).contains(&level) {
            return Err(Error::InvalidInput(format!(
                "Zstd level must be between {MIN_LEVEL} and {MAX_LEVEL}"
            )));
        }
        Ok(Self {
            compressor: zstd::bulk::Compressor::new(level)?,
        })
    }

    pub fn encode(&mut self, input: &[u8]) -> Result<Vec<u8>> {
        self.compressor.compress(input).map_err(Into::into)
    }
}

pub struct Decoder {
    decompressor: zstd::bulk::Decompressor<'static>,
}

impl Decoder {
    pub fn new() -> Result<Self> {
        Ok(Self {
            decompressor: zstd::bulk::Decompressor::new()?,
        })
    }

    pub fn decode(&mut self, input: &[u8], raw_len: usize) -> Result<Vec<u8>> {
        self.decompressor
            .decompress(input, raw_len)
            .map_err(Into::into)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zstd_round_trip() {
        let input = b"ACGTACGTACGTIIIIIIIIIIIIIIII".repeat(1024);
        let mut encoder = Encoder::new(19).unwrap();
        let encoded = encoder.encode(&input).unwrap();
        assert!(encoded.len() < input.len());
        let mut decoder = Decoder::new().unwrap();
        assert_eq!(decoder.decode(&encoded, input.len()).unwrap(), input);
    }

    #[test]
    fn rejects_invalid_level() {
        assert!(Encoder::new(23).is_err());
    }
}
