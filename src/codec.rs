use crate::Result;

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
