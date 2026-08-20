use spring_rs::{
    decompress_paired_into, DecompressionOptions, OutputOrder, PairedRecordSink, ReadRecord,
};
use std::fs::File;
use std::io::BufReader;

#[derive(Default)]
struct CountSink {
    pairs: u64,
    bases: u64,
}

impl PairedRecordSink for CountSink {
    fn write_pair(
        &mut self,
        _original_id: u64,
        first: ReadRecord,
        second: ReadRecord,
    ) -> spring_rs::Result<()> {
        self.pairs += 1;
        self.bases += (first.sequence.len() + second.sequence.len()) as u64;
        Ok(())
    }

    fn finish(&mut self) -> spring_rs::Result<()> {
        Ok(())
    }
}

fn main() -> spring_rs::Result<()> {
    let archive = std::env::args()
        .nth(1)
        .ok_or_else(|| spring_rs::Error::InvalidInput("archive path is required".into()))?;
    let output_order = match std::env::args().nth(2).as_deref() {
        Some("compression") => OutputOrder::Compression,
        Some("original") | None => OutputOrder::Original,
        Some(other) => {
            return Err(spring_rs::Error::InvalidInput(format!(
                "unknown output order {other}"
            )))
        }
    };
    let options = DecompressionOptions {
        output_order,
        temp_dir: Some("/tmp/opencode/spring-rs-decode-sort".into()),
        ..DecompressionOptions::default()
    };
    let mut sink = CountSink::default();
    let started = std::time::Instant::now();
    let metadata = decompress_paired_into(
        BufReader::with_capacity(16 * 1024 * 1024, File::open(archive)?),
        &options,
        &mut sink,
    )?;
    println!(
        "{output_order:?}: {} pairs, {} bases, {:.3}s, {metadata:?}",
        sink.pairs,
        sink.bases,
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}
