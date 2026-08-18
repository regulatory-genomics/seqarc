use crate::packed_sequence::PackedSequence;
use crate::{Error, Result};
use rkyv::{Archive, Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

const MAGIC: &[u8; 8] = b"SPRSP001";
const HEADER_LEN: u64 = 16;

#[derive(Archive, Serialize, Deserialize, Clone, Debug, PartialEq, Eq)]
pub(crate) struct SortableRead {
    pub hash: u64,
    pub original_id: u64,
    pub position: u32,
    pub reverse: bool,
    pub sequence: PackedSequence,
}

impl SortableRead {
    pub(crate) fn compare(left: &Self, right: &Self) -> std::cmp::Ordering {
        (left.hash, left.original_id, left.position, left.reverse).cmp(&(
            right.hash,
            right.original_id,
            right.position,
            right.reverse,
        ))
    }
}

pub(crate) struct SequenceSpoolWriter {
    path: PathBuf,
    file: BufWriter<File>,
    count: u64,
}

impl SequenceSpoolWriter {
    pub(crate) fn create(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let mut file = BufWriter::with_capacity(8 * 1024 * 1024, File::create(&path)?);
        file.write_all(MAGIC)?;
        file.write_all(&0u64.to_le_bytes())?;
        Ok(Self {
            path,
            file,
            count: 0,
        })
    }

    pub(crate) fn push(&mut self, read: &SortableRead) -> Result<()> {
        self.file.write_all(&read.hash.to_le_bytes())?;
        self.file.write_all(&read.original_id.to_le_bytes())?;
        self.file.write_all(&read.position.to_le_bytes())?;
        let flags = (read.reverse as u8) | ((read.sequence.is_packed() as u8) << 1);
        self.file.write_all(&[flags])?;
        self.file.write_all(&read.sequence.len().to_le_bytes())?;
        self.file.write_all(read.sequence.encoded_bytes())?;
        self.count = self
            .count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidInput("too many sequence spool records".into()))?;
        Ok(())
    }

    pub(crate) fn finish(mut self) -> Result<PathBuf> {
        self.file.flush()?;
        self.file
            .get_mut()
            .seek(SeekFrom::Start(MAGIC.len() as u64))?;
        self.file.get_mut().write_all(&self.count.to_le_bytes())?;
        self.file.flush()?;
        Ok(self.path)
    }
}

pub(crate) struct SequenceSpoolReader {
    file: BufReader<File>,
    file_len: u64,
    count: u64,
    read: u64,
    position: u64,
    failed: bool,
}

impl SequenceSpoolReader {
    pub(crate) fn open(path: impl AsRef<Path>) -> Result<Self> {
        let file = File::open(path)?;
        let file_len = file.metadata()?.len();
        if file_len < HEADER_LEN {
            return Err(Error::InvalidFormat(
                "truncated sequence spool header".into(),
            ));
        }
        let mut file = BufReader::with_capacity(8 * 1024 * 1024, file);
        let mut magic = [0u8; 8];
        file.read_exact(&mut magic)?;
        if magic != *MAGIC {
            return Err(Error::InvalidFormat("invalid sequence spool magic".into()));
        }
        let mut count = [0u8; 8];
        file.read_exact(&mut count)?;
        Ok(Self {
            file,
            file_len,
            count: u64::from_le_bytes(count),
            read: 0,
            position: HEADER_LEN,
            failed: false,
        })
    }

    pub(crate) fn record_count(&self) -> u64 {
        self.count
    }

    fn next_record(&mut self) -> Result<Option<SortableRead>> {
        if self.read == self.count {
            if self.position != self.file_len {
                return Err(Error::InvalidFormat(
                    "trailing data in sequence spool".into(),
                ));
            }
            return Ok(None);
        }
        let mut header = [0u8; 25];
        self.file.read_exact(&mut header).map_err(|error| {
            if error.kind() == std::io::ErrorKind::UnexpectedEof {
                Error::InvalidFormat("truncated sequence spool row header".into())
            } else {
                Error::Io(error)
            }
        })?;
        self.position += header.len() as u64;
        let hash = u64::from_le_bytes(header[0..8].try_into().unwrap());
        let original_id = u64::from_le_bytes(header[8..16].try_into().unwrap());
        let position = u32::from_le_bytes(header[16..20].try_into().unwrap());
        let flags = header[20];
        if flags & !3 != 0 {
            return Err(Error::InvalidFormat("invalid sequence spool flags".into()));
        }
        let reverse = flags & 1 != 0;
        let packed = flags & 2 != 0;
        let sequence_len = u32::from_le_bytes(header[21..25].try_into().unwrap());
        let length = if packed {
            (sequence_len as usize).div_ceil(4)
        } else {
            sequence_len as usize
        };
        let remaining = self.file_len.saturating_sub(self.position);
        if length as u64 > remaining {
            return Err(Error::InvalidFormat("truncated sequence spool row".into()));
        }
        let mut bytes = vec![0; length];
        self.file.read_exact(&mut bytes)?;
        self.position += length as u64;
        let value = SortableRead {
            hash,
            original_id,
            position,
            reverse,
            sequence: PackedSequence::from_encoded(sequence_len, packed, bytes)?,
        };
        self.read += 1;
        Ok(Some(value))
    }
}

impl Iterator for SequenceSpoolReader {
    type Item = Result<SortableRead>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.failed {
            return None;
        }
        match self.next_record() {
            Ok(Some(read)) => Some(Ok(read)),
            Ok(None) => None,
            Err(error) => {
                self.failed = true;
                Some(Err(error))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(id: u64, sequence: &[u8]) -> SortableRead {
        SortableRead {
            hash: id % 3,
            original_id: id,
            position: id as u32,
            reverse: !id.is_multiple_of(2),
            sequence: PackedSequence::from_bytes(sequence),
        }
    }

    #[test]
    fn spool_round_trips_packed_and_raw_sequences() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.spool");
        let expected = vec![read(0, b"ACGTACG"), read(1, b"ACNTac")];
        let mut writer = SequenceSpoolWriter::create(&path).unwrap();
        for value in &expected {
            writer.push(value).unwrap();
        }
        writer.finish().unwrap();

        let reader = SequenceSpoolReader::open(path).unwrap();
        assert_eq!(reader.record_count(), expected.len() as u64);
        assert_eq!(reader.collect::<Result<Vec<_>>>().unwrap(), expected);
    }

    #[test]
    fn spool_rejects_truncated_rows() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.spool");
        let mut writer = SequenceSpoolWriter::create(&path).unwrap();
        writer.push(&read(0, b"ACGT")).unwrap();
        writer.finish().unwrap();
        let file = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
        file.set_len(std::fs::metadata(&path).unwrap().len() - 1)
            .unwrap();

        let error = SequenceSpoolReader::open(path)
            .unwrap()
            .next()
            .unwrap()
            .unwrap_err();
        assert!(matches!(error, Error::InvalidFormat(_)));
    }
}
