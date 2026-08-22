use crate::{Error, Result};
use bed_utils::extsort::{ExternalChunk, ExternalChunkBuilder};
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

const CHUNK_SIZE: usize = 8 * 1024 * 1024;

pub(crate) struct Lz4ByteWriter {
    path: PathBuf,
    builder: ExternalChunkBuilder<Vec<u8>>,
    buffer: Vec<u8>,
}

impl Lz4ByteWriter {
    pub(crate) fn create(path: impl AsRef<Path>, level: u32) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::create(&path)?;
        let builder = ExternalChunkBuilder::new(file, level)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        Ok(Self {
            path,
            builder,
            buffer: Vec::with_capacity(CHUNK_SIZE),
        })
    }

    fn write_chunk(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let chunk = std::mem::replace(&mut self.buffer, Vec::with_capacity(CHUNK_SIZE));
        self.builder
            .add(chunk)
            .map_err(|error| Error::InvalidInput(error.to_string()))
    }

    pub(crate) fn finish(mut self) -> Result<PathBuf> {
        self.flush()?;
        let chunk = self
            .builder
            .finish()
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        drop(chunk);
        Ok(self.path)
    }
}

impl Write for Lz4ByteWriter {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut written = 0;
        while written < bytes.len() {
            let available = CHUNK_SIZE - self.buffer.len();
            let amount = available.min(bytes.len() - written);
            self.buffer
                .extend_from_slice(&bytes[written..written + amount]);
            written += amount;
            if self.buffer.len() == CHUNK_SIZE {
                self.write_chunk().map_err(as_io_error)?;
            }
        }
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.write_chunk().map_err(as_io_error)
    }
}

fn as_io_error(error: Error) -> std::io::Error {
    match error {
        Error::Io(error) => error,
        other => std::io::Error::other(other.to_string()),
    }
}

pub(crate) fn copy_to<W: Write>(path: impl AsRef<Path>, output: &mut W) -> Result<u64> {
    let chunk = ExternalChunk::<Vec<u8>>::open(File::open(path)?)
        .map_err(|error| Error::InvalidFormat(error.to_string()))?;
    let mut copied = 0u64;
    for bytes in chunk {
        let bytes = bytes.map_err(|error| Error::InvalidFormat(error.to_string()))?;
        output.write_all(&bytes)?;
        copied = copied
            .checked_add(bytes.len() as u64)
            .ok_or_else(|| Error::InvalidInput("external byte stream too large".into()))?;
    }
    Ok(copied)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_chunked_bytes() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("bytes.lz4");
        let expected = vec![b'A'; CHUNK_SIZE + 123];
        let mut writer = Lz4ByteWriter::create(&path, 1).unwrap();
        writer.write_all(&expected).unwrap();
        writer.finish().unwrap();

        let mut actual = Vec::new();
        assert_eq!(copy_to(path, &mut actual).unwrap(), expected.len() as u64);
        assert_eq!(actual, expected);
    }
}
