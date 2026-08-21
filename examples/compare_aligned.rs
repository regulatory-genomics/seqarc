use flate2::read::MultiGzDecoder;
use spring_rs::{
    decompress_aligned_paired_into, gzip_reader, DecompressionOptions, Error, FastxReader, Format,
    PairedRecordSink, ReadRecord, Result,
};
use std::fs::File;
use std::io::BufReader;

type FastqReader = FastxReader<BufReader<MultiGzDecoder<File>>>;

struct CompareSink {
    first: FastqReader,
    second: FastqReader,
    compared: u64,
}

impl CompareSink {
    fn next_expected(reader: &mut FastqReader) -> Result<ReadRecord> {
        match reader.next() {
            Some(record) => record,
            None => Err(Error::InvalidInput("archive has extra records".into())),
        }
    }

    fn compare_record(
        actual: &ReadRecord,
        expected: &ReadRecord,
        original_id: u64,
        mate: usize,
    ) -> Result<()> {
        if actual != expected {
            return Err(Error::InvalidInput(format!(
                "mismatch at pair {original_id}, mate {mate}"
            )));
        }
        Ok(())
    }
}

impl PairedRecordSink for CompareSink {
    fn write_pair(
        &mut self,
        original_id: u64,
        first: ReadRecord,
        second: ReadRecord,
    ) -> Result<()> {
        let expected_first = Self::next_expected(&mut self.first)?;
        let expected_second = Self::next_expected(&mut self.second)?;
        Self::compare_record(&first, &expected_first, original_id, 1)?;
        Self::compare_record(&second, &expected_second, original_id, 2)?;
        self.compared += 1;
        if self.compared % 100_000 == 0 {
            eprintln!("compared {} pairs", self.compared);
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if self.first.next().is_some() || self.second.next().is_some() {
            return Err(Error::InvalidInput(
                "original FASTQ contains extra records".into(),
            ));
        }
        println!(
            "byte-for-byte comparison passed for {} pairs",
            self.compared
        );
        Ok(())
    }
}

fn main() -> Result<()> {
    let mut args = std::env::args_os().skip(1);
    let archive = args
        .next()
        .ok_or_else(|| Error::InvalidInput("missing archive".into()))?;
    let input_1 = args
        .next()
        .ok_or_else(|| Error::InvalidInput("missing mate-1 FASTQ".into()))?;
    let input_2 = args
        .next()
        .ok_or_else(|| Error::InvalidInput("missing mate-2 FASTQ".into()))?;

    let mut sink = CompareSink {
        first: FastxReader::new(
            BufReader::new(gzip_reader(File::open(input_1)?)),
            Format::Fastq,
        ),
        second: FastxReader::new(
            BufReader::new(gzip_reader(File::open(input_2)?)),
            Format::Fastq,
        ),
        compared: 0,
    };
    decompress_aligned_paired_into(
        BufReader::new(File::open(archive)?),
        &DecompressionOptions::default(),
        &mut sink,
    )?;
    Ok(())
}
