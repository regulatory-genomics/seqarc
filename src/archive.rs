use crate::{
    codec, ArchiveRecords, CompressionMetadata, CompressionOptions, DecompressionOptions, Error,
    QualityMode, ReadRecord, Result,
};
use crossbeam_channel::{bounded, Receiver};
use std::io::{Read, Seek, SeekFrom, Write};

const MAGIC: &[u8; 9] = b"SPRNG-RS2";
const VERSION: u32 = 2;
const MAX_SECTION_SIZE: usize = 1 << 30;

/// One prepared chunk: the per-stream, per-field record data ready to encode.
type ChunkStreams = Vec<Vec<ReadRecord>>;

pub fn write_archive<I, W>(
    streams: Vec<I>,
    mut output: W,
    options: &CompressionOptions,
) -> Result<CompressionMetadata>
where
    I: IntoIterator<Item = Result<ReadRecord>>,
    W: Write,
{
    if streams.is_empty() || streams.len() > u16::MAX as usize {
        return Err(Error::InvalidInput(
            "archive must contain one or more streams".into(),
        ));
    }
    let mut data: Vec<Vec<ReadRecord>> = Vec::with_capacity(streams.len());
    for stream in streams {
        let mut records = Vec::new();
        for item in stream {
            let mut record = item?;
            record.validate(options.preserve_quality)?;
            if !options.preserve_ids {
                record.id.clear();
            }
            if !options.preserve_quality {
                record.quality = None;
            } else {
                apply_quality(&mut record, &options.quality_mode);
            }
            records.push(record);
        }
        data.push(records);
    }
    let count = data[0].len();
    if data.iter().any(|s| s.len() != count) {
        return Err(Error::InvalidInput(
            "paired streams have different record counts".into(),
        ));
    }
    let chunk_records = options.chunk_records.max(1).min(u32::MAX as usize);
    output.write_all(MAGIC)?;
    write_u32(&mut output, VERSION)?;
    write_u16(&mut output, data.len() as u16)?;
    let mut flags = 0u16;
    if options.preserve_ids {
        flags |= 1;
    }
    if options.preserve_quality {
        flags |= 2;
    }
    if options.preserve_order {
        flags |= 4;
    }
    write_u16(&mut output, flags)?;
    write_u64(&mut output, count as u64)?;
    write_u32(&mut output, chunk_records as u32)?;
    write_quality_mode(&mut output, &options.quality_mode)?;
    let chunks = count.div_ceil(chunk_records);
    write_u64(&mut output, chunks as u64)?;

    // Slice the fully-materialized streams into owned per-chunk views, then run
    // the same parallel pipeline used for streaming input.
    let ids = options.preserve_ids;
    let quality = options.preserve_quality;
    let chunk_iter = (0..count).step_by(chunk_records).map(|start| {
        let end = (start + chunk_records).min(count);
        Ok(data
            .iter()
            .map(|stream| stream[start..end].to_vec())
            .collect::<ChunkStreams>())
    });
    encode_pipeline(chunk_iter, &mut output, options, ids, quality)?;

    Ok(CompressionMetadata {
        streams: data.len(),
        records: count as u64,
        preserved_ids: options.preserve_ids,
        preserved_quality: options.preserve_quality,
        preserved_order: options.preserve_order,
    })
}

pub fn write_archive_streams<I, W>(
    mut streams: Vec<I>,
    mut output: W,
    options: &CompressionOptions,
) -> Result<CompressionMetadata>
where
    I: Iterator<Item = Result<ReadRecord>> + Send,
    W: Write + Seek,
{
    if streams.is_empty() || streams.len() > u16::MAX as usize {
        return Err(Error::InvalidInput(
            "archive must contain one or more streams".into(),
        ));
    }
    let num_streams = streams.len();
    let chunk_records = options.chunk_records.max(1).min(u32::MAX as usize);
    let (record_count_position, chunk_count_position) =
        write_header(&mut output, num_streams, options, chunk_records)?;

    let ids = options.preserve_ids;
    let quality = options.preserve_quality;
    let preserve_quality = options.preserve_quality;
    let preserve_ids = options.preserve_ids;
    let quality_mode = options.quality_mode.clone();

    // Producer side: pull records from each stream and assemble prepared chunks.
    let chunk_iter = std::iter::from_fn(move || {
        let mut chunk: ChunkStreams = Vec::with_capacity(num_streams);
        let mut any = false;
        for stream in &mut streams {
            let mut records_for_stream = Vec::with_capacity(chunk_records);
            for _ in 0..chunk_records {
                match stream.next() {
                    Some(Ok(mut record)) => {
                        any = true;
                        if let Err(error) = record.validate(preserve_quality) {
                            return Some(Err(error));
                        }
                        if !preserve_ids {
                            record.id.clear();
                        }
                        if !preserve_quality {
                            record.quality = None;
                        } else {
                            apply_quality(&mut record, &quality_mode);
                        }
                        records_for_stream.push(record);
                    }
                    Some(Err(error)) => return Some(Err(error)),
                    None => break,
                }
            }
            chunk.push(records_for_stream);
        }
        if !any {
            return None;
        }
        let count = chunk[0].len();
        if chunk.iter().any(|stream| stream.len() != count) {
            return Some(Err(Error::InvalidInput(
                "paired streams have different record counts".into(),
            )));
        }
        Some(Ok(chunk))
    });

    let (records, chunks) = encode_pipeline(chunk_iter, &mut output, options, ids, quality)?;

    let end = output.stream_position()?;
    output.seek(SeekFrom::Start(record_count_position))?;
    write_u64(&mut output, records)?;
    output.seek(SeekFrom::Start(chunk_count_position))?;
    write_u64(&mut output, chunks)?;
    output.seek(SeekFrom::Start(end))?;

    Ok(CompressionMetadata {
        streams: num_streams,
        records,
        preserved_ids: options.preserve_ids,
        preserved_quality: options.preserve_quality,
        preserved_order: options.preserve_order,
    })
}

/// Producer/consumer compression pipeline.
///
/// A single producer (the calling thread) yields prepared chunks in order. A
/// pool of worker threads compresses chunks independently, each with its own
/// single-threaded zstd encoder. A bounded ordering buffer reassembles encoded
/// chunks into input order before writing so the archive is deterministic.
fn encode_pipeline<I, W>(
    chunks: I,
    output: &mut W,
    options: &CompressionOptions,
    ids: bool,
    quality: bool,
) -> Result<(u64, u64)>
where
    I: Iterator<Item = Result<ChunkStreams>> + Send,
    W: Write,
{
    let level = options.zstd_level;
    codec::Encoder::new(level)?;
    let workers = options.threads.max(1);
    // Active workers already retain one raw chunk each. Keep only half a
    // worker-pool of queued jobs/results to bound memory without starving the
    // compressors while the producer parses the next chunk.
    let capacity = workers.div_ceil(2).max(1);
    let (job_tx, job_rx) = bounded::<(usize, ChunkStreams)>(capacity);
    let (done_tx, done_rx) = bounded::<Result<(usize, usize, Vec<u8>)>>(capacity);

    std::thread::scope(|scope| -> Result<(u64, u64)> {
        let producer_errors = done_tx.clone();
        scope.spawn(move || {
            for (index, chunk) in chunks.enumerate() {
                match chunk {
                    Ok(chunk) => {
                        if job_tx.send((index, chunk)).is_err() {
                            return;
                        }
                    }
                    Err(error) => {
                        let _ = producer_errors.send(Err(error));
                        return;
                    }
                }
            }
        });

        for _ in 0..workers {
            let job_rx = job_rx.clone();
            let done_tx = done_tx.clone();
            scope.spawn(move || {
                let mut encoder = match codec::Encoder::new(level) {
                    Ok(encoder) => encoder,
                    Err(error) => {
                        let _ = done_tx.send(Err(error));
                        return;
                    }
                };
                while let Ok((index, chunk)) = job_rx.recv() {
                    let count = chunk.first().map_or(0, Vec::len);
                    let encoded = encode_chunk_to_vec(&mut encoder, &chunk, ids, quality);
                    if done_tx
                        .send(encoded.map(|bytes| (index, count, bytes)))
                        .is_err()
                    {
                        break;
                    }
                }
            });
        }
        drop(job_rx);
        drop(done_tx);

        write_ordered(done_rx, output)
    })
}

fn write_ordered<W>(
    done_rx: Receiver<Result<(usize, usize, Vec<u8>)>>,
    output: &mut W,
) -> Result<(u64, u64)>
where
    W: Write,
{
    use std::collections::BTreeMap;

    let mut pending: BTreeMap<usize, (usize, Vec<u8>)> = BTreeMap::new();
    let mut next_to_write = 0usize;
    let mut total_records = 0u64;
    let mut total_chunks = 0u64;

    while let Ok(result) = done_rx.recv() {
        let (index, count, bytes) = result?;
        pending.insert(index, (count, bytes));
        while let Some((count, bytes)) = pending.remove(&next_to_write) {
            output.write_all(&bytes)?;
            total_records += count as u64;
            total_chunks += 1;
            next_to_write += 1;
        }
    }

    if !pending.is_empty() {
        return Err(Error::InvalidInput(
            "compression workers exited before all chunks were written".into(),
        ));
    }

    Ok((total_records, total_chunks))
}

fn write_header<W: Write + Seek>(
    output: &mut W,
    streams: usize,
    options: &CompressionOptions,
    chunk_records: usize,
) -> Result<(u64, u64)> {
    output.write_all(MAGIC)?;
    write_u32(output, VERSION)?;
    write_u16(output, streams as u16)?;
    let mut flags = 0u16;
    if options.preserve_ids {
        flags |= 1;
    }
    if options.preserve_quality {
        flags |= 2;
    }
    if options.preserve_order {
        flags |= 4;
    }
    write_u16(output, flags)?;
    let record_count_position = output.stream_position()?;
    write_u64(output, 0)?;
    write_u32(output, chunk_records as u32)?;
    write_quality_mode(output, &options.quality_mode)?;
    let chunk_count_position = output.stream_position()?;
    write_u64(output, 0)?;
    Ok((record_count_position, chunk_count_position))
}

pub fn read_archive<R: Read>(
    mut input: R,
    options: &DecompressionOptions,
) -> Result<ArchiveRecords> {
    let mut magic = [0; 9];
    input.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(Error::InvalidFormat("bad magic".into()));
    }
    if read_u32(&mut input)? != VERSION {
        return Err(Error::InvalidFormat("unsupported version".into()));
    }
    let streams = read_u16(&mut input)? as usize;
    let flags = read_u16(&mut input)?;
    let records = read_u64(&mut input)? as usize;
    let chunk_records = read_u32(&mut input)? as usize;
    let quality_mode = read_quality_mode(&mut input)?;
    let chunks = read_u64(&mut input)? as usize;
    if streams == 0 || chunk_records == 0 {
        return Err(Error::InvalidFormat("invalid archive dimensions".into()));
    }
    let mut output = ArchiveRecords {
        streams: (0..streams).map(|_| Vec::new()).collect(),
    };
    let mut decoder = codec::Decoder::new()?;
    for chunk in 0..chunks {
        let chunk_len = read_u32(&mut input)? as usize;
        let start = chunk * chunk_records;
        let end = start + chunk_len;
        let decoded = read_chunk(
            &mut input,
            &mut decoder,
            streams,
            chunk_len,
            flags & 1 != 0,
            flags & 2 != 0,
            chunk,
        )?;
        for (stream_index, stream) in decoded.into_iter().enumerate() {
            for (offset, record) in stream.into_iter().enumerate() {
                let index = start + offset;
                if index >= records {
                    return Err(Error::InvalidFormat("chunk exceeds record count".into()));
                }
                if selected(options, index as u64) {
                    output.streams[stream_index].push(record);
                }
            }
        }
        if end > records {
            return Err(Error::InvalidFormat("chunk exceeds record count".into()));
        }
    }
    let _ = quality_mode;
    Ok(output)
}

fn selected(options: &DecompressionOptions, index: u64) -> bool {
    options
        .range
        .as_ref()
        .is_none_or(|range| range.contains(&index))
}

fn apply_quality(record: &mut ReadRecord, mode: &QualityMode) {
    if let Some(quality) = record.quality.as_mut() {
        match mode {
            QualityMode::Lossless => {}
            QualityMode::IlluminaBinning => {
                for q in quality {
                    *q = match *q {
                        b'!'..=b'$' => b'!',
                        b'%'..=b'(' => b'%',
                        b')'..=b',' => b')',
                        b'-'..=b'1' => b'-',
                        b'2'..=b'6' => b'2',
                        b'7'..=b';' => b'7',
                        _ => b'@',
                    };
                }
            }
            QualityMode::Binary {
                threshold,
                high,
                low,
            } => {
                for q in quality {
                    *q = if q.saturating_sub(33) >= *threshold {
                        *high + 33
                    } else {
                        *low + 33
                    };
                }
            }
        }
    }
}

/// Serializes and compresses one prepared chunk into a self-contained buffer.
///
/// The output is byte-identical regardless of how many worker threads are used,
/// which keeps archives deterministic.
fn encode_chunk_to_vec(
    encoder: &mut codec::Encoder,
    streams: &[Vec<ReadRecord>],
    ids: bool,
    quality: bool,
) -> Result<Vec<u8>> {
    let count = streams.first().map_or(0, Vec::len);
    let mut output = Vec::new();
    write_u32(&mut output, count as u32)?;
    for stream in streams {
        let mut lengths = Vec::with_capacity(count * if quality { 12 } else { 8 });
        let mut id_data = Vec::new();
        let mut sequence_data = Vec::new();
        let mut quality_data = Vec::new();
        for record in stream {
            let id = if ids { record.id.as_slice() } else { &[] };
            push_len(&mut lengths, id.len())?;
            push_len(&mut lengths, record.sequence.len())?;
            if quality {
                let qualities = record.quality.as_deref().unwrap_or(&[]);
                push_len(&mut lengths, qualities.len())?;
                quality_data.extend_from_slice(qualities);
            }
            id_data.extend_from_slice(id);
            sequence_data.extend_from_slice(&record.sequence);
        }
        write_section(&mut output, encoder, &lengths)?;
        write_section(&mut output, encoder, &id_data)?;
        write_section(&mut output, encoder, &sequence_data)?;
        if quality {
            write_section(&mut output, encoder, &quality_data)?;
        }
    }
    Ok(output)
}

fn read_chunk<R: Read>(
    input: &mut R,
    decoder: &mut codec::Decoder,
    streams: usize,
    count: usize,
    ids: bool,
    quality: bool,
    chunk: usize,
) -> Result<Vec<Vec<ReadRecord>>> {
    let mut output = Vec::with_capacity(streams);
    for _ in 0..streams {
        let lengths = read_section(input, decoder, chunk)?;
        let id_data = read_section(input, decoder, chunk)?;
        let sequence_data = read_section(input, decoder, chunk)?;
        let quality_data = if quality {
            Some(read_section(input, decoder, chunk)?)
        } else {
            None
        };
        let fields = if quality { 3 } else { 2 };
        if lengths.len() != count.saturating_mul(fields).saturating_mul(4) {
            return Err(Error::InvalidFormat("invalid chunk lengths table".into()));
        }
        let mut lengths_pos = 0;
        let mut id_pos = 0;
        let mut sequence_pos = 0;
        let mut quality_pos = 0;
        let mut records = Vec::with_capacity(count);
        for _ in 0..count {
            let id_len = take_len(&lengths, &mut lengths_pos)?;
            let sequence_len = take_len(&lengths, &mut lengths_pos)?;
            let quality_len = if quality {
                take_len(&lengths, &mut lengths_pos)?
            } else {
                0
            };
            let id = take_field(&id_data, &mut id_pos, id_len)?;
            let sequence = take_field(&sequence_data, &mut sequence_pos, sequence_len)?;
            let qualities = match quality_data.as_deref() {
                Some(data) => Some(take_field(data, &mut quality_pos, quality_len)?),
                None => None,
            };
            records.push(ReadRecord {
                id: if ids { id } else { Vec::new() },
                sequence,
                quality: qualities,
            });
        }
        if lengths_pos != lengths.len()
            || id_pos != id_data.len()
            || sequence_pos != sequence_data.len()
            || quality_data
                .as_ref()
                .is_some_and(|data| quality_pos != data.len())
        {
            return Err(Error::InvalidFormat("trailing chunk section data".into()));
        }
        output.push(records);
    }
    Ok(output)
}

fn write_section<W: Write>(output: &mut W, encoder: &mut codec::Encoder, raw: &[u8]) -> Result<()> {
    let encoded = encoder.encode(raw)?;
    write_u64(output, raw.len() as u64)?;
    write_u64(output, encoded.len() as u64)?;
    write_u32(output, checksum(raw))?;
    output.write_all(&encoded)?;
    Ok(())
}

fn read_section<R: Read>(
    input: &mut R,
    decoder: &mut codec::Decoder,
    chunk: usize,
) -> Result<Vec<u8>> {
    let raw_len = usize_from_u64(read_u64(input)?, "raw section")?;
    let encoded_len = usize_from_u64(read_u64(input)?, "encoded section")?;
    if raw_len > MAX_SECTION_SIZE || encoded_len > MAX_SECTION_SIZE {
        return Err(Error::InvalidFormat("section exceeds size limit".into()));
    }
    let expected = read_u32(input)?;
    let mut encoded = vec![0; encoded_len];
    input.read_exact(&mut encoded)?;
    let raw = decoder.decode(&encoded, raw_len)?;
    if raw.len() != raw_len || checksum(&raw) != expected {
        return Err(Error::InvalidFormat(format!(
            "checksum failure in chunk {chunk}"
        )));
    }
    Ok(raw)
}

fn push_len(output: &mut Vec<u8>, len: usize) -> Result<()> {
    let len =
        u32::try_from(len).map_err(|_| Error::InvalidInput("record field exceeds 4 GiB".into()))?;
    output.extend_from_slice(&len.to_le_bytes());
    Ok(())
}

fn take_len(input: &[u8], pos: &mut usize) -> Result<usize> {
    if input.len().saturating_sub(*pos) < 4 {
        return Err(Error::InvalidFormat("truncated lengths table".into()));
    }
    let len = u32::from_le_bytes(input[*pos..*pos + 4].try_into().unwrap()) as usize;
    *pos += 4;
    Ok(len)
}

fn take_field(input: &[u8], pos: &mut usize, len: usize) -> Result<Vec<u8>> {
    if input.len().saturating_sub(*pos) < len {
        return Err(Error::InvalidFormat("truncated chunk field".into()));
    }
    let value = input[*pos..*pos + len].to_vec();
    *pos += len;
    Ok(value)
}

fn usize_from_u64(value: u64, field: &str) -> Result<usize> {
    usize::try_from(value)
        .map_err(|_| Error::InvalidFormat(format!("{field} length does not fit usize")))
}
fn checksum(bytes: &[u8]) -> u32 {
    bytes.iter().fold(2166136261u32, |hash, byte| {
        (hash ^ *byte as u32).wrapping_mul(16777619)
    })
}
fn write_u16<W: Write>(w: &mut W, v: u16) -> Result<()> {
    w.write_all(&v.to_le_bytes()).map_err(Into::into)
}
fn write_u32<W: Write>(w: &mut W, v: u32) -> Result<()> {
    w.write_all(&v.to_le_bytes()).map_err(Into::into)
}
fn write_u64<W: Write>(w: &mut W, v: u64) -> Result<()> {
    w.write_all(&v.to_le_bytes()).map_err(Into::into)
}
fn read_u16<R: Read>(r: &mut R) -> Result<u16> {
    let mut b = [0; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}
fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
fn read_u64<R: Read>(r: &mut R) -> Result<u64> {
    let mut b = [0; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
fn write_quality_mode<W: Write>(w: &mut W, mode: &QualityMode) -> Result<()> {
    match mode {
        QualityMode::Lossless => {
            w.write_all(&[0])?;
        }
        QualityMode::IlluminaBinning => {
            w.write_all(&[1])?;
        }
        QualityMode::Binary {
            threshold,
            high,
            low,
        } => {
            w.write_all(&[2, *threshold, *high, *low])?;
        }
    }
    Ok(())
}
fn read_quality_mode<R: Read>(r: &mut R) -> Result<QualityMode> {
    let mut tag = [0; 1];
    r.read_exact(&mut tag)?;
    match tag[0] {
        0 => Ok(QualityMode::Lossless),
        1 => Ok(QualityMode::IlluminaBinning),
        2 => {
            let mut b = [0; 3];
            r.read_exact(&mut b)?;
            Ok(QualityMode::Binary {
                threshold: b[0],
                high: b[1],
                low: b[2],
            })
        }
        _ => Err(Error::InvalidFormat("unknown quality mode".into())),
    }
}
