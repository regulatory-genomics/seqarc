use crate::{Error, ReadRecord, Result};
use std::fs::{create_dir_all, File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug)]
struct FieldRef {
    offset: u64,
    length: u32,
}

#[derive(Clone, Copy, Debug)]
struct PairIndex {
    id_1: FieldRef,
    seq_1: FieldRef,
    qual_1: FieldRef,
    id_2: FieldRef,
    seq_2: FieldRef,
    qual_2: FieldRef,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StoreStats {
    pub index_reads: u64,
    pub blob_reads: u64,
    pub bytes_read: u64,
}

pub(crate) struct PairMetadata {
    pub id: Vec<u8>,
    pub quality: Vec<u8>,
}

pub struct PairStore {
    root: PathBuf,
    ids_1: File,
    seqs_1: File,
    quals_1: File,
    ids_2: File,
    seqs_2: File,
    quals_2: File,
    index: File,
    count: u64,
    stats: StoreStats,
}

impl PairStore {
    pub fn create(root: impl AsRef<Path>) -> Result<Self> {
        let root = root.as_ref().to_path_buf();
        create_dir_all(&root)?;
        Ok(Self {
            ids_1: create_blob(&root, "ids.1")?,
            seqs_1: create_blob(&root, "seqs.1")?,
            quals_1: create_blob(&root, "quals.1")?,
            ids_2: create_blob(&root, "ids.2")?,
            seqs_2: create_blob(&root, "seqs.2")?,
            quals_2: create_blob(&root, "quals.2")?,
            index: create_blob(&root, "pairs.index")?,
            root,
            count: 0,
            stats: StoreStats::default(),
        })
    }

    pub fn append(&mut self, first: &ReadRecord, second: &ReadRecord) -> Result<u64> {
        let entry = PairIndex {
            id_1: append_field(&mut self.ids_1, &first.id)?,
            seq_1: append_field(&mut self.seqs_1, &first.sequence)?,
            qual_1: append_field(&mut self.quals_1, first.quality.as_deref().unwrap_or(&[]))?,
            id_2: append_field(&mut self.ids_2, &second.id)?,
            seq_2: append_field(&mut self.seqs_2, &second.sequence)?,
            qual_2: append_field(&mut self.quals_2, second.quality.as_deref().unwrap_or(&[]))?,
        };
        write_index(&mut self.index, entry)?;
        let id = self.count;
        self.count += 1;
        Ok(id)
    }

    pub(crate) fn open(root: impl AsRef<Path>) -> Result<Self> {
        let root = root.as_ref().to_path_buf();
        let index = open_blob(&root, "pairs.index")?;
        let index_len = index.metadata()?.len();
        if index_len % INDEX_SIZE != 0 {
            return Err(Error::InvalidFormat(
                "record store index has an invalid length".into(),
            ));
        }
        Ok(Self {
            ids_1: open_blob(&root, "ids.1")?,
            seqs_1: open_blob(&root, "seqs.1")?,
            quals_1: open_blob(&root, "quals.1")?,
            ids_2: open_blob(&root, "ids.2")?,
            seqs_2: open_blob(&root, "seqs.2")?,
            quals_2: open_blob(&root, "quals.2")?,
            index,
            root,
            count: index_len / INDEX_SIZE,
            stats: StoreStats::default(),
        })
    }

    pub fn count(&self) -> u64 {
        self.count
    }

    pub fn read_pair(&mut self, pair_id: u64) -> Result<(ReadRecord, ReadRecord)> {
        let entry = self.read_index_entry(pair_id)?;
        let first = ReadRecord {
            id: read_field_counted(&mut self.ids_1, entry.id_1, &mut self.stats)?,
            sequence: read_field_counted(&mut self.seqs_1, entry.seq_1, &mut self.stats)?,
            quality: optional_field_counted(&mut self.quals_1, entry.qual_1, &mut self.stats)?,
        };
        let second = ReadRecord {
            id: read_field_counted(&mut self.ids_2, entry.id_2, &mut self.stats)?,
            sequence: read_field_counted(&mut self.seqs_2, entry.seq_2, &mut self.stats)?,
            quality: optional_field_counted(&mut self.quals_2, entry.qual_2, &mut self.stats)?,
        };
        Ok((first, second))
    }

    pub fn read_first(&mut self, pair_id: u64) -> Result<ReadRecord> {
        Ok(self.read_pair(pair_id)?.0)
    }

    pub fn read_sequence(&mut self, pair_id: u64, mate: usize) -> Result<Vec<u8>> {
        let entry = self.read_index_entry(pair_id)?;
        match mate {
            0 => read_field_counted(&mut self.seqs_1, entry.seq_1, &mut self.stats),
            1 => read_field_counted(&mut self.seqs_2, entry.seq_2, &mut self.stats),
            _ => Err(Error::InvalidInput("mate index must be 0 or 1".into())),
        }
    }

    pub fn read_mate(&mut self, pair_id: u64, mate: usize) -> Result<ReadRecord> {
        let entry = self.read_index_entry(pair_id)?;
        match mate {
            0 => Ok(ReadRecord {
                id: read_field_counted(&mut self.ids_1, entry.id_1, &mut self.stats)?,
                sequence: read_field_counted(&mut self.seqs_1, entry.seq_1, &mut self.stats)?,
                quality: optional_field_counted(&mut self.quals_1, entry.qual_1, &mut self.stats)?,
            }),
            1 => Ok(ReadRecord {
                id: read_field_counted(&mut self.ids_2, entry.id_2, &mut self.stats)?,
                sequence: read_field_counted(&mut self.seqs_2, entry.seq_2, &mut self.stats)?,
                quality: optional_field_counted(&mut self.quals_2, entry.qual_2, &mut self.stats)?,
            }),
            _ => Err(Error::InvalidInput("mate index must be 0 or 1".into())),
        }
    }

    pub(crate) fn read_pair_metadata(
        &mut self,
        pair_id: u64,
        ids: bool,
        quality: bool,
    ) -> Result<(PairMetadata, PairMetadata)> {
        let entry = self.read_index_entry(pair_id)?;
        let first_id = if ids {
            read_field_counted(&mut self.ids_1, entry.id_1, &mut self.stats)?
        } else {
            Vec::new()
        };
        let second_id = if ids {
            read_field_counted(&mut self.ids_2, entry.id_2, &mut self.stats)?
        } else {
            Vec::new()
        };
        let first_quality = if quality {
            read_field_counted(&mut self.quals_1, entry.qual_1, &mut self.stats)?
        } else {
            Vec::new()
        };
        let second_quality = if quality {
            read_field_counted(&mut self.quals_2, entry.qual_2, &mut self.stats)?
        } else {
            Vec::new()
        };
        Ok((
            PairMetadata {
                id: first_id,
                quality: first_quality,
            },
            PairMetadata {
                id: second_id,
                quality: second_quality,
            },
        ))
    }

    pub fn stats(&self) -> StoreStats {
        self.stats
    }
    pub fn reset_stats(&mut self) {
        self.stats = StoreStats::default();
    }

    pub(crate) fn add_stats(&mut self, stats: StoreStats) {
        self.stats.index_reads += stats.index_reads;
        self.stats.blob_reads += stats.blob_reads;
        self.stats.bytes_read += stats.bytes_read;
    }

    fn read_index_entry(&mut self, pair_id: u64) -> Result<PairIndex> {
        if pair_id >= self.count {
            return Err(Error::InvalidInput("pair id out of bounds".into()));
        }
        self.index.seek(SeekFrom::Start(pair_id * INDEX_SIZE))?;
        self.stats.index_reads += 1;
        self.stats.bytes_read += INDEX_SIZE;
        read_index(&mut self.index)
    }

    pub fn root(&self) -> &Path {
        &self.root
    }
}

const FIELD_SIZE: usize = 12;
const INDEX_SIZE: u64 = (FIELD_SIZE * 6) as u64;

fn create_blob(root: &Path, name: &str) -> Result<File> {
    Ok(OpenOptions::new()
        .create(true)
        .read(true)
        .write(true)
        .truncate(true)
        .open(root.join(name))?)
}

fn open_blob(root: &Path, name: &str) -> Result<File> {
    Ok(OpenOptions::new().read(true).open(root.join(name))?)
}

fn append_field(file: &mut File, bytes: &[u8]) -> Result<FieldRef> {
    let offset = file.seek(SeekFrom::End(0))?;
    let length = u32::try_from(bytes.len())
        .map_err(|_| Error::InvalidInput("record field exceeds 4 GiB".into()))?;
    file.write_all(bytes)?;
    Ok(FieldRef { offset, length })
}

fn write_index(file: &mut File, entry: PairIndex) -> Result<()> {
    for field in [
        entry.id_1,
        entry.seq_1,
        entry.qual_1,
        entry.id_2,
        entry.seq_2,
        entry.qual_2,
    ] {
        file.write_all(&field.offset.to_le_bytes())?;
        file.write_all(&field.length.to_le_bytes())?;
    }
    Ok(())
}

fn read_index(file: &mut File) -> Result<PairIndex> {
    let mut fields = [FieldRef {
        offset: 0,
        length: 0,
    }; 6];
    for field in &mut fields {
        let mut offset = [0; 8];
        let mut length = [0; 4];
        file.read_exact(&mut offset)?;
        file.read_exact(&mut length)?;
        *field = FieldRef {
            offset: u64::from_le_bytes(offset),
            length: u32::from_le_bytes(length),
        };
    }
    Ok(PairIndex {
        id_1: fields[0],
        seq_1: fields[1],
        qual_1: fields[2],
        id_2: fields[3],
        seq_2: fields[4],
        qual_2: fields[5],
    })
}

fn read_field(file: &mut File, field: FieldRef) -> Result<Vec<u8>> {
    file.seek(SeekFrom::Start(field.offset))?;
    let mut bytes = vec![0; field.length as usize];
    file.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn read_field_counted(file: &mut File, field: FieldRef, stats: &mut StoreStats) -> Result<Vec<u8>> {
    let bytes = read_field(file, field)?;
    stats.blob_reads += 1;
    stats.bytes_read += field.length as u64;
    Ok(bytes)
}

fn optional_field_counted(
    file: &mut File,
    field: FieldRef,
    stats: &mut StoreStats,
) -> Result<Option<Vec<u8>>> {
    let bytes = read_field_counted(file, field, stats)?;
    Ok((!bytes.is_empty()).then_some(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stores_and_reads_pairs_without_retaining_records() {
        let temp = tempfile::tempdir().unwrap();
        let mut store = PairStore::create(temp.path()).unwrap();
        let first = ReadRecord::new("r1", b"ACGT".to_vec(), Some(b"IIII".to_vec()));
        let second = ReadRecord::new("r2", b"TGCA".to_vec(), Some(b"JJJJ".to_vec()));
        assert_eq!(store.append(&first, &second).unwrap(), 0);
        assert_eq!(store.read_pair(0).unwrap(), (first.clone(), second.clone()));

        let mut reader = PairStore::open(temp.path()).unwrap();
        assert_eq!(reader.count(), 1);
        assert_eq!(reader.read_pair(0).unwrap(), (first, second));
    }
}
