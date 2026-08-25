# SeqArc Aligned Archive Format

This document specifies the current aligned SeqArc archive format written by
the aligned single-end and paired-end compression APIs.

## Conventions

- The archive stores either single-end reads or paired-end reads. Counts and
  ordinals refer to logical records: one read for single-end archives and one
  read pair for paired-end archives.
- All fixed-width integers are unsigned little-endian unless explicitly stated.
- Signed coordinate differences are encoded with ZigZag followed by an unsigned
  varint.
- ZigZag maps signed value `x` to `((x << 1) XOR (x >> 63))` interpreted as
  `u64`.
- Varints use unsigned LEB128: seven payload bits per byte, least-significant
  group first, with bit 7 set when another byte follows.
- Bitmap record `i` is stored in bit `i % 8` of byte `i / 8`; the first record
  uses the least-significant bit.
- Zstandard frames use the compression level selected by the encoder. The level
  is not stored in the archive and is not needed for decompression.
- Version 20 is the only supported aligned archive version.

## Top-Level Layout

```text
global header                         39 bytes
embedded compact reference block     reference_block_length bytes
embedded index metadata block        index_metadata_block_length bytes
chunk 0
chunk 1
...
chunk N-1
```

## Global Header

The global header is exactly 39 bytes.

| Offset | Size | Type | Field |
|---:|---:|---|---|
| 0 | 9 | bytes | Magic: `SEQARC-AL` |
| 9 | 4 | `u32` | Version: `20` |
| 13 | 2 | `u16` | Flags |
| 15 | 8 | `u64` | Reference block length |
| 23 | 8 | `u64` | Index metadata block length |
| 31 | 8 | `u64` | Chunk count |

### Flags

| Bits | Meaning |
|---|---|
| 0 | Read identifiers are preserved |
| 1 | Quality values are preserved |
| 2 | Single-end layout when set; paired-end layout when clear |
| 3-4 | Identifier mode |
| 5-15 | Reserved; must be zero |

Identifier mode values are:

| Value | Mode | Meaning |
|---:|---|---|
| 0 | None | No identifiers are stored |
| 1 | Independent | Each identifier stream is stored independently |
| 2 | Identical mates | Only mate 1 is stored; mate 2 is identical |
| 3 | Slash mates | A shared stem is stored; `/1` and `/2` are restored |

The identifier mode is selected for the complete archive, not independently
for each chunk. Single-end archives permit only `None` and `Independent`.

The stream count is derived from bit 2:

| Layout | Stream count |
|---|---:|
| Single-end | 1 |
| Paired-end | 2 |

## Embedded Compact Reference

The archive is self-contained. It embeds the subset of reference bases used to
reconstruct selected mapped alignments.

The encoder marks every reference position consumed by `M`, `=`, or `X` CIGAR
operations. `D` and `N` spans do not retain reference bases by themselves.
Unmarked positions are removed, and the remaining bases are concatenated in
original reference-coordinate order. Alignment starts and reference-consuming
CIGAR lengths are translated into this compact coordinate space before chunk
encoding. A stored `D` or `N` length is therefore the number of globally
retained bases in its original interval; zero-length translated operations are
omitted.

The run mapping between the original reference and the compact reference is not
stored because archive records already use compact coordinates.

### Reference Block

```text
zstd_frame(raw_reference_payload)
```

The complete Zstandard frame is counted by the global header's
`reference_block_length` field. Exactly one frame is permitted.

### Raw Reference Payload

```text
u64 base_count
u8  packed_bases[ceil(base_count / 4)]
```

The compact reference contains only uppercase A/C/G/T. Genome FASTA processing
uppercases canonical bases and drops all other input bases before index creation.

Four canonical bases are encoded with two bits:

| Base | Code |
|---|---:|
| A | 0 |
| C | 1 |
| G | 2 |
| T | 3 |

Within each packed byte, the first base occupies bits 7-6, the second bits 5-4,
the third bits 3-2, and the fourth bits 1-0. Unused low bits in the final packed
byte must be zero.

## Embedded Index Metadata

The metadata block is stored immediately after the compact reference. Archives
without embedded indexes store a zero range count.

The metadata block starts with:

```text
u64 range_count
```

Each range descriptor follows in command-line declaration order:

```text
u8  source_mate                 0 for R1, 1 for R2
u64 zero_based_start
u64 length
```

Ranges use coordinates in the original biological read. Ranges on the same mate
must be nonempty and non-overlapping. The CLI accepts one-based inclusive
coordinates and converts them to the stored zero-based start and length.

After all descriptors, the metadata stores one merged whitelist:

```text
u64 global_whitelist_length
u8  global_whitelist[global_whitelist_length]
```

Each range first predicts its own local whitelist across every logical record.
Local whitelists retain their counts for correction decisions during compression,
but are not serialized. Their entries are deduplicated into one global whitelist;
counts are dropped and unique sequences are sorted lexicographically.

The global whitelist payload contains a varint entry count followed by entries:

```text
repeat entry_count times:
    varint sequence_length
    u8     sequence[sequence_length]
```

Whitelist entries contain uppercase A/C/G/T and are longer than one base.

## Chunks

Chunks preserve original logical-record order and are independently
decompressible. The default encoder writes 250,000 records per chunk, but chunk
size is not fixed by the format.

### Chunk Header

Each chunk starts with a 44-byte header.

| Size | Type | Field |
|---:|---|---|
| 4 | `u32` | Logical record count |
| 8 | `u64` | Compressed payload length |
| 8 | `u64` | Total raw section length |
| 8 | `u64` | Raw identifier section length |
| 8 | `u64` | Raw quality section length |
| 8 | `u64` | Raw record section length |

The total raw length must equal the sum of all declared raw section lengths.

### Compressed Chunk Payload

The compressed payload is a section directory followed by three independent
Zstandard frames:

```text
u64 compressed_identifier_length
u64 compressed_quality_length
u64 compressed_record_length
u8  compressed_identifiers[compressed_identifier_length]
u8  compressed_qualities[compressed_quality_length]
u8  compressed_records[compressed_record_length]
```

The chunk header's `compressed_payload_length` includes the 24-byte section
directory. Each frame must decompress to the corresponding raw length from the
chunk header.

## Identifier Section

An encoder may preprocess identifiers before archive encoding, including
retaining only a Unicode-character prefix. The leading FASTQ `@` is syntax and
is not part of the identifier. The archive stores and reconstructs the resulting
identifier bytes exactly; discarded suffixes and the preprocessing length are
not stored and cannot be recovered during decompression.

The raw identifier section contains one length-delimited block per stream:

```text
repeat stream_count times:
    u64 block_length
    u8  block[block_length]
```

For paired `Identical mates` or `Slash mates` mode, the second block is present
with zero length.

Identifiers inside a non-empty block are encoded in logical-record order. For each
identifier:

```text
varint common_prefix_length_with_previous_identifier
varint suffix_length
u8     suffix[suffix_length]
```

The previous identifier is initially empty. The common prefix must not exceed
the previous identifier length.

For slash-mate mode, the encoder removes the final `/1` or `/2` before prefix
encoding. The decoder appends `/1` to mate 1 and `/2` to mate 2.

## Quality Section

The raw quality section contains one length-delimited quality stream per archive
stream:

```text
repeat stream_count times:
    u64 stream_length
    u8  stream[stream_length]
```

When qualities are not preserved, every stream length is zero.

A stream is a sequence of raw quality runs. There is no run count; runs are read
until exactly the declared logical-record count has been reconstructed and the
enclosing stream boundary is reached. Consecutive records with the same quality
width are grouped into one run.

### Quality Run

```text
u32 width
u32 record_count
u8  qualities[width * record_count]
```

`record_count` must be nonzero and must not exceed the number of records remaining
in the stream. The quality body is row-major and contains exactly one byte per
quality value. Zero-width rows are valid and therefore have an empty body.

Quality preprocessing may be lossless, Illumina-binned, or binary depending on
encoder options. The archive format records only the resulting quality bytes;
it does not record which preprocessing mode produced them.

When embedded index ranges are configured, extraction is performed before
alignment. Their quality fragments are retained temporarily, reinserted before
quality compression, and stored as part of the ordinary complete R1/R2 quality
rows. There is no separate index-quality stream.

## Record Section

The record section stores one section per stream. For paired archives all
stream-0 records are followed by all stream-1 records. Each stream section has
the following layout:

```text
u8  mapped_bitmap[ceil(record_count / 8)]
u8  reverse_bitmap[ceil(record_count / 8)]
u64 positions_length
u8  positions[positions_length]
record_column_directory
u8  programs[programs_length]
u8  raw_bitmap[raw_bitmap_byte_length]
u8  mismatch_bitmap[mismatch_bitmap_byte_length]
u8  ordinals[ordinal_byte_length]
u8  mismatch_details[mismatch_detail_byte_length]
u8  literals[literals_length]
```

There is one body for every logical-record ordinal in the stream. Whether a body
is mapped or unmapped is selected by the mapped bitmap. Position varints exist
only for mapped records.

### Position Prediction

Position values are signed deltas encoded as ZigZag varints.

For stream 0, including every single-end archive:

- The predictor is the previous mapped mate-1 start in the chunk.
- The initial predictor is zero.

For paired stream 1:

- If mate 1 of the same pair is mapped, the predictor is that mate-1 start.
- Otherwise, the predictor is the previous mapped mate-2 start in the chunk.
- The previous mate-2 start is updated for every mapped mate-2 record.

The positions substream must be consumed exactly after all mapped records for
the mate are decoded.

### Strand Prediction

For stream 0, a set reverse bit means reverse strand.

For paired stream 1:

- If mate 1 of the same pair is mapped, the stored bit is
  `mate_2_reverse XOR mate_1_reverse`.
- Otherwise, the stored bit is the absolute mate-2 reverse value.

### Record Column Directory

```text
u64 programs_length
u64 raw_bitmap_bit_count
u64 raw_bitmap_byte_length
u64 mismatch_bitmap_bit_count
u64 mismatch_bitmap_byte_length
u64 ordinal_byte_length
u64 mismatch_detail_byte_length
u64 literals_length
```

For the raw and mismatch bitmap columns, byte length must equal
`ceil(bit_count / 8)`, and unused low bits in the final byte must be zero. The
directory and all six columns are stored once per mate, not once per record.

### Program Column

The program column concatenates exactly `record_count` programs in logical-record
order. Each program is parsed from its operation count and operation-derived
payload lengths:

```text
varint operation_count
u8     packed_operation_codes[ceil(operation_count / 4)]
varint lengths[non_whitelist_operation_count]
payload match_payloads
```

Archive operation codes are two-bit values. The first operation occupies bits
7-6, and four operations are stored per byte. Unused low bits in the final byte
must be zero.

| Code | Operation | Meaning |
|---:|---|---|
| 0 | `M` Match | Copy reference bases and apply mismatch noise |
| 1 | `L` Literal | Emit adaptive packed/raw/sparse query bases |
| 2 | `R` Reference advance | Advance reference without emitting query bases |
| 3 | `W` Whitelist | Emit one grouped embedded-index boundary |

The mapper's aligned-region CIGAR does not include unaligned query prefixes or
suffixes. SeqArc inserts synthetic soft-clip operations so the effective CIGAR
consumes the complete reduced read. Archive normalization maps `M`, `=`, and `X`
to `M`; maps `I` and `S` to `L`; maps `D` and `N` to `R`; removes `H`, `P`, and
zero-length operations; and merges adjacent operations of the same kind.

Whitelist operations are inserted at the reduced-query boundaries represented
by global index descriptors. A `W` operation has no stored length, consumes
neither reduced query nor reference, and emits one or more descriptor fragments.
Descriptors sharing the same reduced-query boundary are grouped into one `W`
operation, with one token per descriptor. Core `M` or `L` operations are split
when a whitelist boundary falls inside them. Each descriptor owned by the mate
must be represented exactly once, and each unique boundary must have exactly one
`W` operation.

Mapped bodies may contain `M`, `L`, `R`, and `W`. Unmapped bodies use the same
format but may contain only `L` and `W`.

### Operation Payloads

#### Match

The aligned query is compared against the compact reference:

```text
varint mismatch_count
repeat mismatch_count times:
    varint position_delta
    u8     replacement_base
```

Mismatch positions are offsets within the current `M` operation. The first
position is relative to zero; subsequent positions are relative to the previous
mismatch position.

#### Literal And Raw Whitelist Sources

The literal column contains only concatenated sequence bytes. It has no per-run
mode, length, packing, or sparse-exception framing. After whitelist columns are
decoded, the encoder and decoder derive maximal literal runs from consecutive
output sources:

- ordinary `L` operation bases;
- raw descriptors within a grouped `W` operation;
- sources separated only by `R`, which emits no output.

An `M` operation or an exact/one-mismatch whitelist descriptor ends the current
literal run. The literal column stores each run's oriented bytes directly in
derived order. No literal-run count, separator, or length is stored because all
boundaries are derived from operation lengths, descriptor lengths, and the raw
bitmap. The outer record-section Zstandard frame compresses these byte-aligned
DNA symbols.

#### Reference Advance

No payload is stored; only the reference position advances.

#### Whitelist

Whitelist columns follow descriptor occurrences in record and operation order.
Within a grouped `W`, descriptor occurrences follow original-coordinate order for
forward/unmapped reads and reverse original-coordinate order for reverse mapped
reads.

| Raw bit | Meaning |
|---:|---|
| 0 | Descriptor is encoded from the global whitelist |
| 1 | Descriptor sequence is supplied by a shared literal run |

The raw bitmap column contains exactly one bit per descriptor occurrence across
all records in record order and grouped-W descriptor order.

After raw membership is known, the mismatch bitmap column contains exactly one
mode bit per non-raw descriptor:

| Mode bit | Encoding |
|---:|---|
| 0 | Exact global-whitelist ordinal |
| 1 | Global ordinal, one mismatch position, and original base |

Raw descriptors store no mismatch bit, ordinal, or mismatch details.

The global whitelist may contain at most 65,536 entries. The ordinal column
stores one little-endian `u16` ordinal for every non-raw descriptor. Byte
alignment preserves repeated ordinal symbols for the outer Zstandard codec.

The mismatch-detail column stores one byte for every set mismatch bit:

```text
bits 7-2  mismatch position, 0..63
bits 1-0  original A/C/G/T base
```

One-mismatch candidates for descriptors longer than 64 bases are encoded as raw
literal sources. Reads requiring two or more whitelist mismatches, and
noncanonical reads, also set their raw bitmap bit and use the literal column.

Candidate selection uses only the descriptor's local predicted whitelist, while
the selected sequence is represented by its archive-global ordinal. This
prevents another descriptor's entries from affecting correction.

For reverse mapped reads, boundary groups and descriptors within each group occur
in descending original-coordinate order. Decoded whitelist fragments are
reverse-complemented while constructing the oriented read, and the complete
reconstructed read is reverse-complemented back to FASTQ orientation.

## Integrity and Validation

The aligned format has no archive-level checksum. Decoders rely on:

- exact magic and version checks;
- reserved flag checks;
- fixed and aggregate length checks;
- Zstandard frame decoding;
- exact section and substream consumption;
- end-of-file immediately after the declared final chunk;
- bounds checks for coordinates, archive operations, whitelist fields, mismatch
  positions, sparse exceptions, reference payload lengths, and packed-reference
  padding.

Writers should reject values that cannot fit their specified integer widths.
Readers must use checked arithmetic when combining stored lengths and offsets.

## Reconstruction Guarantees

An aligned archive preserves logical-record order and its single/paired
cardinality. Depending on encoder options it can also preserve identifiers and
qualities. Sequence reconstruction is lossless:

- mapped reads are rebuilt from the embedded compact reference, archive
  operations, strand, and noise;
- unmapped reads are rebuilt from literal and whitelist operations;
- embedded indexes are emitted inline from the global whitelist or lossless raw
  fallbacks at their original R1 or R2 coordinates.

If a lossy quality preprocessing mode was selected, the stored quality bytes
are reproduced exactly, but they are not the original input qualities.
