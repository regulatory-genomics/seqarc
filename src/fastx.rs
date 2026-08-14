use crate::{Error, ReadRecord, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::{BufRead, Write};

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
