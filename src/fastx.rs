use crate::{Error, ReadRecord, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

const INPUT_BUFFER_SIZE: usize = 1024 * 1024;
const DECODED_BUFFER_SIZE: usize = 8 * 1024 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum InputCompression {
    Plain,
    Gzip,
    Zstd,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Fastq,
    Fasta,
}

pub fn gzip_reader<R: std::io::Read>(reader: R) -> MultiGzDecoder<R> {
    MultiGzDecoder::new(reader)
}

pub fn gzip_writer<W: Write>(writer: W, level: u32) -> GzEncoder<W> {
    GzEncoder::new(writer, Compression::new(level.min(9)))
}

/// Open a plain, gzip, or zstd-compressed input detected from its magic bytes.
pub fn open_fastx_reader(path: impl AsRef<Path>) -> Result<Box<dyn BufRead + Send>> {
    let mut input = BufReader::with_capacity(INPUT_BUFFER_SIZE, File::open(path)?);
    let compression = detect_compression(input.fill_buf()?);
    let reader: Box<dyn BufRead + Send> = match compression {
        InputCompression::Plain => Box::new(input),
        InputCompression::Gzip => Box::new(BufReader::with_capacity(
            DECODED_BUFFER_SIZE,
            gzip_reader(input),
        )),
        InputCompression::Zstd => Box::new(BufReader::with_capacity(
            DECODED_BUFFER_SIZE,
            zstd::stream::read::Decoder::with_buffer(input)?,
        )),
    };
    Ok(reader)
}

fn detect_compression(prefix: &[u8]) -> InputCompression {
    if prefix.starts_with(&[0x1f, 0x8b]) {
        return InputCompression::Gzip;
    }
    if prefix.len() >= 4 {
        let magic = u32::from_le_bytes(prefix[..4].try_into().unwrap());
        if magic == 0xfd2f_b528 || magic & 0xffff_fff0 == 0x184d_2a50 {
            return InputCompression::Zstd;
        }
    }
    InputCompression::Plain
}

pub struct FastxReader<R> {
    reader: R,
    format: Format,
    pending: Option<String>,
}

impl<R: BufRead> FastxReader<R> {
    pub fn new(reader: R, format: Format) -> Self {
        Self {
            reader,
            format,
            pending: None,
        }
    }
}

impl<R: BufRead> Iterator for FastxReader<R> {
    type Item = Result<ReadRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.format {
            Format::Fastq => self.next_fastq(),
            Format::Fasta => self.next_fasta(),
        }
    }
}

impl<R: BufRead> FastxReader<R> {
    fn line(&mut self) -> Result<Option<String>> {
        let mut line = String::new();
        let n = self.reader.read_line(&mut line)?;
        if n == 0 {
            Ok(None)
        } else {
            Ok(Some(line.trim_end_matches(['\n', '\r']).to_owned()))
        }
    }
    fn next_fastq(&mut self) -> Option<Result<ReadRecord>> {
        let id = match self.line() {
            Ok(Some(v)) => v,
            Ok(None) => return None,
            Err(e) => return Some(Err(e)),
        };
        let seq = match self.line() {
            Ok(Some(v)) => v,
            Ok(None) => return Some(Err(Error::InvalidInput("truncated FASTQ".into()))),
            Err(e) => return Some(Err(e)),
        };
        let plus = match self.line() {
            Ok(Some(v)) => v,
            Ok(None) => return Some(Err(Error::InvalidInput("truncated FASTQ".into()))),
            Err(e) => return Some(Err(e)),
        };
        let qual = match self.line() {
            Ok(Some(v)) => v,
            Ok(None) => return Some(Err(Error::InvalidInput("truncated FASTQ".into()))),
            Err(e) => return Some(Err(e)),
        };
        if !id.starts_with('@') || !plus.starts_with('+') {
            return Some(Err(Error::InvalidInput("invalid FASTQ record".into())));
        }
        Some(Ok(ReadRecord::new(
            id.as_bytes()[1..].to_vec(),
            seq.into_bytes(),
            Some(qual.into_bytes()),
        )))
    }
    fn next_fasta(&mut self) -> Option<Result<ReadRecord>> {
        let id = if let Some(pending) = self.pending.take() {
            pending
        } else {
            match self.line() {
                Ok(Some(line)) => line,
                Ok(None) => return None,
                Err(error) => return Some(Err(error)),
            }
        };
        if !id.starts_with('>') {
            return Some(Err(Error::InvalidInput("invalid FASTA header".into())));
        }
        let mut sequence = String::new();
        loop {
            match self.line() {
                Ok(Some(line)) if line.starts_with('>') => {
                    self.pending = Some(line);
                    break;
                }
                Ok(Some(line)) => sequence.push_str(&line),
                Ok(None) => break,
                Err(e) => return Some(Err(e)),
            }
        }
        Some(Ok(ReadRecord::new(
            id.as_bytes()[1..].to_vec(),
            sequence.into_bytes(),
            None,
        )))
    }
}

pub struct FastxWriter<W> {
    writer: W,
    format: Format,
}
impl<W: Write> FastxWriter<W> {
    pub fn new(writer: W, format: Format) -> Self {
        Self { writer, format }
    }
    pub fn write_record(&mut self, record: &ReadRecord) -> Result<()> {
        match self.format {
            Format::Fastq => {
                let q = record
                    .quality
                    .as_ref()
                    .ok_or_else(|| Error::InvalidInput("FASTQ output requires quality".into()))?;
                if q.len() != record.sequence.len() {
                    return Err(Error::InvalidInput("quality length mismatch".into()));
                }
                writeln!(self.writer, "@{}", String::from_utf8_lossy(&record.id))?;
                self.writer.write_all(&record.sequence)?;
                self.writer.write_all(b"\n+\n")?;
                self.writer.write_all(q)?;
                self.writer.write_all(b"\n")?;
            }
            Format::Fasta => {
                writeln!(self.writer, ">{}", String::from_utf8_lossy(&record.id))?;
                self.writer.write_all(&record.sequence)?;
                self.writer.write_all(b"\n")?;
            }
        }
        Ok(())
    }
    pub fn finish(mut self) -> Result<()> {
        self.writer.flush().map_err(Into::into)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Read, Write};

    #[test]
    fn opens_plain_gzip_and_concatenated_zstd_with_misleading_extensions() {
        let directory = tempfile::tempdir().unwrap();
        let parts = [
            b"@read0\nACGT\n+\nIIII\n".as_slice(),
            b"@read1\nTGCA\n+\nJJJJ\n",
        ];
        let expected = parts.concat();

        let plain = directory.path().join("plain.fastq.zst");
        std::fs::write(&plain, &expected).unwrap();

        let gzip = directory.path().join("gzip.fastq");
        let mut gzip_writer = gzip_writer(File::create(&gzip).unwrap(), 1);
        gzip_writer.write_all(&expected).unwrap();
        gzip_writer.finish().unwrap();

        let zstd = directory.path().join("zstd.fastq.gz");
        let mut encoded = Vec::from([0x50, 0x2a, 0x4d, 0x18, 0, 0, 0, 0]);
        for part in parts {
            let mut encoder = zstd::stream::write::Encoder::new(&mut encoded, 1).unwrap();
            encoder.write_all(part).unwrap();
            encoder.finish().unwrap();
        }
        std::fs::write(&zstd, encoded).unwrap();

        for path in [plain, gzip, zstd] {
            let mut decoded = Vec::new();
            open_fastx_reader(path)
                .unwrap()
                .read_to_end(&mut decoded)
                .unwrap();
            assert_eq!(decoded, expected);
        }
    }

    #[test]
    fn malformed_zstd_is_reported_while_reading() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("truncated.fastq");
        std::fs::write(&path, [0x28, 0xb5, 0x2f, 0xfd, 0x00]).unwrap();
        let mut reader = open_fastx_reader(path).unwrap();
        let mut decoded = Vec::new();
        assert!(reader.read_to_end(&mut decoded).is_err());
    }
}
