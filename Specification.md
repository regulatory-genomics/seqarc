# SeqArc Aligned Archive Format

This document specifies the current aligned SeqArc archive format written by
the aligned single-end and paired-end compression APIs.

## Conventions

- The archive stores either single-end reads or paired-end reads, optionally
  with up to two index reads (I1 and I2) taken from separate FASTQ files.
  Counts and ordinals refer to logical records: one read for single-end
  archives, one read pair for paired-end archives. Each logical record also
  carries one stored index read per declared index stream. All input FASTQ
  files share the same identifier proper (the first whitespace-delimited token)
  per logical record; trailing descriptions may differ between files. Index
  reads are never stored with their own identifiers and are reconstructed using
  the biological mate-1 identifier, so after any compression-time identifier
  preprocessing the index-read identifiers must equal the transformed mate-1
  identifier.
- All fixed-width integers are unsigned little-endian unless explicitly stated.
- Varints use unsigned LEB128: seven payload bits per byte, least-significant
  group first, with bit 7 set when another byte follows.
- Zstandard frames use the compression level selected by the encoder. The level
  is not stored in the archive and is not needed for decompression.
- Version 0.28 is the only supported aligned archive version.

## Top-Level Layout

```text
global header                         46 bytes
embedded compact reference block     reference_block_length bytes
embedded index metadata block        index_metadata_block_length bytes
chunk 0
chunk 1
...
chunk N-1
```

## Global Header

The global header is exactly 46 bytes.

| Offset | Size | Type | Field |
|---:|---:|---|---|
| 0 | 9 | bytes | Magic: `SEQARC-AL` |
| 9 | 2 | `u16` | Version major: `0` |
| 11 | 2 | `u16` | Version minor: `28` |
| 13 | 2 | `u16` | Flags |
| 15 | 1 | `u8` | Index read stream count: `0`, `1`, or `2` |
| 16 | 1 | `u8` | Absolute position width: `1`, `2`, `3`, `4`, or `8` |
| 17 | 1 | `u8` | Delta position width: `1`, `2`, `3`, `4`, or `8` |
| 18 | 8 | `u64` | Reference block length |
| 26 | 8 | `u64` | Index metadata block length |
| 34 | 8 | `u64` | Chunk count |
| 42 | 4 | `u32` | CRC32C checksum |

The absolute and delta position widths select the fixed-width encoding used by
the R1/R2 position subsections; both must be one of `1`, `2`, `3`, `4`, or `8`.

The CRC32C checksum is computed over the first 42 header bytes (magic through
chunk count) followed by the complete embedded reference block and the
complete index metadata block. It is verified before the reference block is
decoded, so corruption anywhere in the archive preamble is detected before any
alignment data is reconstructed. Because the chunk count is covered, writers
must patch this checksum after writing the final chunk count.

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
without embedded indexes store a zero range count. Declared index read streams
do not occupy descriptor entries; their count lives in the global header.

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

The metadata block ends immediately after the last descriptor.

## Chunks

Chunks preserve original logical-record order and are independently
decompressible. The default encoder writes 250,000 records per chunk, but chunk
size is not fixed by the format.

### Chunk Header

Each chunk starts with a 48-byte header.

| Size | Type | Field |
|---:|---|---|
| 4 | `u32` | Logical record count |
| 8 | `u64` | Compressed payload length |
| 8 | `u64` | Total raw section length |
| 8 | `u64` | Raw identifier section length |
| 8 | `u64` | Raw quality section length |
| 8 | `u64` | Record section length |
| 4 | `u32` | CRC32C checksum |

The total raw length must equal the sum of the raw identifier section length,
the raw quality section length, and the record section length.

The record section length is the byte length of the record section's compressed
subsection directory (see Record Section); unlike the identifier and quality
lengths, it is not a decompressed size.

The CRC32C checksum is computed over the first 44 header bytes (the six fields
above, in order) followed by the complete compressed payload. It detects
corruption in the chunk framing as well as in the payload. Readers must verify
the checksum before decompressing any chunk payload.

### Compressed Chunk Payload

The compressed payload is a 24-byte section directory followed by two
independent Zstandard frames and the record section's subsection directory:

```text
u64 compressed_identifier_length
u64 compressed_quality_length
u64 compressed_record_length
u8  compressed_identifiers[compressed_identifier_length]
u8  compressed_qualities[compressed_quality_length]
u8  record_section[compressed_record_length]
```

The identifier and quality frames must each decompress to the corresponding raw
length from the chunk header. The record section is consumed verbatim; it is
not wrapped in an outer frame, and its byte length must match
`compressed_record_length` and the chunk header's record section length.

## Identifier Section

An encoder may preprocess identifiers before archive encoding, including
retaining only a Unicode-character prefix or retaining only one
whitespace-delimited field of the identifier (zero-based; runs of consecutive
ASCII whitespace act as one separator). The leading FASTQ `@` is syntax and
is not part of the identifier. The archive stores and reconstructs the resulting
identifier bytes exactly; discarded identifier bytes and the preprocessing
parameters are not stored and cannot be recovered during decompression.

The raw identifier section contains one length-delimited block per biological
stream. Index-read identifiers are stored separately, as record-section
subsections:

```text
repeat biological_stream_count times:
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

The raw quality section contains one length-delimited quality stream per
biological archive stream. Index-read qualities are stored separately, as
record-section subsections:

```text
repeat biological_stream_count times:
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

The record section is stored verbatim inside the chunk payload — it is not
wrapped in an outer Zstandard frame. It is itself a directory of independently
compressed Zstandard frames, one per subsection, so each subsection can be
decompressed on its own:

```text
u32 subsection_count
repeat subsection_count times:
    u64 compressed_length
    u64 raw_length
    u8  zstd_frame[compressed_length]
```

The subsection count is:

```text
1 + 8 * biological_stream_count + embedded_index_range_count
    + 3 * index_stream_count
```

Subsections occur in this order:

1. Alignment flags.
2. Eight R1 subsections.
3. Eight R2 subsections for paired archives.
4. One raw subsection for each embedded index range in declaration order.
5. Four subsections for each declared index read stream in I1, I2 order.

Every frame must decode to its declared raw length, and every subsection must be
consumed exactly.

### Alignment Flags

The flag subsection contains one byte per logical record:

| Value | Meaning |
|---:|---|
| `0` | All biological reads are aligned |
| `1` | Only R1 is aligned |
| `2` | Only R2 is aligned |
| `3` | No biological read is aligned |

Single-end archives permit only `0` and `3`.

### Biological Subsections

Each biological stream has eight independently compressed subsections:

1. Read lengths: one little-endian `u16` per logical record.
2. Positions: one compact-reference position per aligned read, written
   little-endian at the width declared by the global header. R1 positions and
   R2-only positions are absolute values using `absolute_position_width`; R2
   positions for records with both mates aligned are signed R2-minus-R1 deltas
   encoded with ZigZag and written using `delta_position_width`. Encoders must
   reject values that exceed the selected width.
3. Noise positions: one little-endian `u16` position delta per noise symbol.
4. Noise symbols: one ASCII row per aligned read, terminated by newline.
5. Reference skips: one little-endian `u64` for every `4` noise symbol.
6. Emitted sequences: one nonempty raw sequence per `5` symbol, terminated by newline.
7. Reverse flags: one byte, `0` or `1`, per aligned read. R1 and R2-only flags
   are absolute orientation values. For a record with both mates aligned, the
   R2 flag is `1` when R2 has the same orientation as R1 and `0` otherwise.
8. Unaligned sequences: concatenated raw bytes; read lengths determine boundaries.

Noise-position deltas are measured from the query cursor after the preceding
event. Exact reference bases between events are copied implicitly. For paired
records with both mates aligned, R2 position and orientation are decoded
relative to R1 as described above. Symbols are:

| Symbol | Meaning |
|---:|---|
| `0`-`3` | Substitute one base using its ordinal among the four bases other than the reference base |
| `4` | Advance the reference by the next value from the reference-skip subsection |
| `5` | Emit the next newline-delimited sequence without advancing the reference |

The substitution alphabet is A/C/G/T/N. Removing the current reference base
from that ordered alphabet leaves exactly four replacement choices represented
by `0` through `3`.

The mapper CIGAR is made query-complete by synthetic soft clips. `M`, `=`, and
`X` produce exact copies and substitution events; `D` and `N` produce `4`
events; `I` and `S` produce `5` events; `H`, `P`, and zero-length operations are
discarded. Reverse reads are encoded in alignment orientation and reverse
complemented after reconstruction.

### Embedded Index Subsections

Each configured index range has one independently compressed raw-byte frame.
Sequences occur in logical-record order and have the fixed length declared by
the metadata descriptor, so no per-record delimiter is stored. Decoders reinsert
the ranges into R1 or R2 in original-coordinate order. Index-range qualities
remain part of the ordinary complete quality rows.

### Index-Read File Subsections

Each declared index read stream has three independently compressed subsections,
in this order:

1. Quality stream: a sequence of raw quality runs as described by the quality
   run format, containing one run-delimited body per logical record. Empty
   when qualities are not preserved.
2. Read lengths: one little-endian `u16` per logical record. Index read
   sequences are nonempty, and writers must reject lengths that cannot fit
   `u16`.
3. Sequences: concatenated raw bytes; the read lengths determine boundaries.

Index reads are never stored with identifiers. All input FASTQ files must
share the same identifier proper per logical record, which compression
validates; trailing descriptions may differ between files. Decoders reconstruct
index reads using the biological mate-1 identifier. Index reads are never
mapped and never contribute mapped flags, positions, CIGAR data, or reference
coverage.

## Integrity and Validation

Every chunk carries a CRC32C checksum over its 44-byte field block and
compressed payload, verified before payload decompression. This localizes
corruption to individual chunks and enables:

- `verify`: a checksum-only scan that validates the global header, chunk
  framing, and every chunk checksum without decompressing;
- best-effort decompression: when enabled, chunks that fail their checksum are
  skipped and records are recovered from all undamaged chunks; the operation
  reports the skipped chunk count and estimated lost records, and exits with
  an error.

Beyond the preamble and per-chunk checksums, decoders rely on:

- exact magic and version checks;
- reserved flag checks;
- fixed and aggregate length checks;
- Zstandard frame decoding;
- exact section and substream consumption;
- end-of-file immediately after the declared final chunk;
- bounds checks for coordinates, read lengths, noise deltas and symbols,
  reference skips, emitted-sequence counts, index ranges, index read lengths,
  reference payload lengths, and packed-reference padding.

Writers should reject values that cannot fit their specified integer widths.
Readers must use checked arithmetic when combining stored lengths and offsets.

## Reconstruction Guarantees

An aligned archive preserves logical-record order and its single/paired
cardinality. Depending on encoder options it can also preserve identifiers and
qualities. Sequence reconstruction is lossless:

- mapped reads are rebuilt from the embedded compact reference, strand, noise,
  reference skips, and emitted sequences;
- unmapped reads are rebuilt from their raw sequence subsections;
- embedded indexes are read from independent raw frames and reinserted at their
  original R1 or R2 coordinates;
- index reads are restored from their independent quality, length, and
  sequence subsections in I1, I2 order, reusing the biological mate-1
  identifier.

If a lossy quality preprocessing mode was selected, the stored quality bytes
are reproduced exactly, but they are not the original input qualities.
