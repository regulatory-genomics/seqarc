# SeqArc

SeqArc is a Rust archiver for sequencing reads. Instead of compressing raw
FASTQ bytes, it aligns reads against a reference genome and stores a compact
alignment-based representation of the data, compressed with Zstandard. The
resulting archive reconstructs the original reads exactly, while ordinary
compression tools only see redundant sequence text. For typical sequencing
data, a SeqArc archive is roughly **30% of the size of the corresponding
gzip-compressed FASTQ**, and decompression is **much faster than gzip**.

## Design goals

SeqArc is built around **asymmetric compression**: compression is allowed to
take longer, but decompression is very fast.

Traditional sequencing archivers such as Spring achieve very high compression
ratios, but decompression takes a long time. Generic tools such as gzip or
zstd decompress quickly, but their compression ratio on FASTQ data is poor.
SeqArc bridges this gap: its compression ratio is close to Spring's, while
decompression runs at super-fast speeds similar to zstd. This makes it well
suited for storing FASTQ files that may need to be read again from time to
time.

SeqArc also supports single-cell sequencing data, including UMI/barcode
extraction (embedded index ranges) and index reads (I1/I2).

## How it works

- Reads are aligned to a reference genome with a built-in aligner
  (BWA-style via `minibwa`, or `rustar`). Mapped reads are stored as compact
  positions, orientation flags, and a per-read noise stream; unmapped reads
  are stored as raw bytes.
- Every archive embeds the subset of reference bases it needs (2-bit packed),
  so archives are fully self-contained and need no external reference.
- Records are grouped into independently decompressible chunks (default
  250,000 records per chunk), enabling metadata inspection, checksum-only
  verification, and partial decoding by record range.
- The global header and every chunk carry a CRC32C checksum, so corruption is
  detected and localized; damaged chunks can be skipped and records recovered
  from undamaged ones.

## Building

Build from source with a Rust toolchain:

```bash
cargo build --release
```

The binary is then available at `target/release/seqarc`.

## Usage

A typical workflow is to build a genome index once, compress FASTQ files
against it, and decompress the archive when the reads are needed again:

```bash
seqarc index genome.fa -o genome.seqarc-index
seqarc compress --index genome.seqarc-index \
    -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -o reads.seqarc
seqarc decompress reads.seqarc -1 out_R1.fastq.gz -2 out_R2.fastq.gz
```

Plain, gzip, and zstd FASTQ inputs are all supported; output compression is
selected by the output file suffix (or forced with `--gzip`/`--plain`).

### Commands

| Command | Description |
|---|---|
| `index` | Build a single-file SeqArc genome index from a FASTA reference |
| `compress` | Compress single-end or paired FASTQ files into an archive |
| `decompress` | Decompress an archive back into FASTQ (or FASTA) files |
| `validate` | Fully decode an archive and optionally byte-compare it against the original FASTQ files |
| `verify` | Scan chunk checksums and framing without decompressing |
| `info` | Print archive metadata (version, layout, options, chunk and record counts) without decoding records |

### Useful options

- `--no-ids` / `--no-quality` — drop read identifiers or quality values.
- `--id-prefix-chars N` / `--id-field N` — retain only part of each identifier.
- `--quality-mode` — lossless (default), `illumina` binning, or `binary`
  thresholding (configured with `--binary-threshold/high/low`).
- `--index-read-1` / `--index-read-2` — store separate index reads (I1/I2).
- `--index-range R1:START:END` — extract a read region (e.g., a UMI) into an
  embedded index before alignment and reinsert it on decompression.
- `--range START:END` — decompress only a range of logical records.
- `--best-effort` — skip chunks that fail their checksum and recover records
  from all undamaged chunks.

Run `seqarc <command> --help` for the complete option list of each command.

## Archive format

The on-disk aligned archive format is specified in detail in
[Specification.md](Specification.md), including the global header, embedded
compact reference, chunk layout, integrity checksums, and reconstruction
guarantees.
