pub mod fasta;

use crate::star::index::SectionSink;

use crate::star::error::Error;
use crate::star::params::Parameters;

use fasta::parse_fasta_files;

/// STAR's genome spacing character (used for inter-chromosome padding).
const GENOME_SPACING_CHAR: u8 = 5;

/// Backing storage for a genome's `[forward | reverse-complement]` sequence.
///
/// `Owned` is the full `2*n_genome` byte buffer built at genomeGenerate time
/// (it is the only variant that supports slicing/mutation). `Mapped` is a
/// read-only memory map of the on-disk `Genome` file, which holds **only the
/// forward strand** (`n_genome` bytes): the reverse-complement half is computed
/// on access in [`GenomeSeq::base`], so loading never materializes the ~`n`-byte
/// RC buffer and the forward bytes are reclaimable file-backed pages rather than
/// an anonymous `Vec`. `Arc<Mmap>` keeps `Genome::clone` (two-pass) cheap.
#[derive(Clone)]
pub enum GenomeSeq {
    Owned(Vec<u8>),
    Mapped {
        fwd: std::sync::Arc<memmap2::Mmap>,
        n_genome: usize,
    },
}

impl GenomeSeq {
    /// Base at absolute position `i` — forward `[0, n_genome)` or
    /// reverse-complement `[n_genome, 2*n_genome)`. For the `Mapped` RC half,
    /// `base(i) = complement(forward[2*n_genome - 1 - i])`, exactly the bytes
    /// the owned builder writes into the second half.
    #[inline]
    pub fn base(&self, i: usize) -> u8 {
        match self {
            GenomeSeq::Owned(v) => v[i],
            GenomeSeq::Mapped { fwd, n_genome } => {
                let n = *n_genome;
                if i < n {
                    fwd[i]
                } else {
                    let f = fwd[2 * n - 1 - i];
                    if f < 4 { 3 - f } else { f }
                }
            }
        }
    }

    /// Total sequence length (`2*n_genome` — forward + reverse complement).
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            GenomeSeq::Owned(v) => v.len(),
            GenomeSeq::Mapped { n_genome, .. } => 2 * n_genome,
        }
    }

    /// Bounds-checked [`base`](Self::base): the base at `i`, or `None` if out of
    /// range.
    #[inline]
    pub fn get(&self, i: usize) -> Option<u8> {
        if i < self.len() {
            Some(self.base(i))
        } else {
            None
        }
    }

    /// The contiguous byte buffer. For `Owned` this is the full
    /// `[forward | RC]`; for `Mapped` it is the forward strand only — callers
    /// that may touch the RC half must use [`base`](Self::base). Used at build
    /// time (always `Owned`) for SA construction and the on-disk write.
    pub fn as_slice(&self) -> &[u8] {
        match self {
            GenomeSeq::Owned(v) => v,
            GenomeSeq::Mapped { fwd, .. } => fwd,
        }
    }
}

impl From<Vec<u8>> for GenomeSeq {
    fn from(v: Vec<u8>) -> Self {
        GenomeSeq::Owned(v)
    }
}

// `memmap2::Mmap` is neither `Debug` nor `PartialEq`, so derive them by hand via
// the byte view. `as_slice()` is the full buffer for `Owned` (the only variant
// tests construct), so equality/printing behave like the old `Vec<u8>` field.
impl std::fmt::Debug for GenomeSeq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GenomeSeq({} bytes)", self.len())
    }
}

impl PartialEq for GenomeSeq {
    fn eq(&self, other: &Self) -> bool {
        self.as_slice() == other.as_slice()
    }
}

/// Packed genome with chromosome metadata.
///
/// The genome sequence is stored as one byte per base:
/// - A=0, C=1, G=2, T=3, N=4, padding=5
/// - Chromosomes are concatenated with padding to bin boundaries
/// - The reverse complement occupies the second half of the `sequence` buffer
#[derive(Clone)]
pub struct Genome {
    /// Forward genome (0..n_genome) + reverse complement (n_genome..2*n_genome).
    /// Owned at build time; a memory map of the forward strand (RC computed on
    /// access) when loaded from disk. Access bases via [`GenomeSeq::base`].
    pub sequence: GenomeSeq,

    /// Total length of the forward (padded) genome.
    pub n_genome: u64,

    /// Length of the real-genome forward strand, excluding any appended
    /// Gsj splice-junction buffer. Equals `n_genome` before `append_sjdb`
    /// and stays pinned at that value afterwards; matches STAR's
    /// `mapGen.nGenomeReal`.
    pub n_genome_real: u64,

    /// Number of real chromosomes (not including scaffold/contigs if excluded).
    pub n_chr_real: usize,

    /// Chromosome names.
    pub chr_name: Vec<String>,

    /// True (unpadded) chromosome lengths.
    pub chr_length: Vec<u64>,

    /// Padded start positions of each chromosome in the genome.
    /// Length = n_chr_real + 1; the last entry is n_genome (total size).
    pub chr_start: Vec<u64>,

}

/// Compute STAR's bin-padded chromosome start offsets: `chr_start[i]` is chromosome
/// `i`'s 0-based start in the padded forward genome; the trailing `chr_start[n]` is
/// the total padded (forward) genome size. Shared by [`Genome::from_fasta`] and
/// the sjdb insertion pipeline.
pub(crate) fn compute_chr_starts(chr_lengths: &[u64], chr_bin_nbits: u32) -> Vec<u64> {
    let bin_size = 1u64 << chr_bin_nbits;
    let mut chr_start = Vec::with_capacity(chr_lengths.len() + 1);
    let mut n: u64 = 0;
    for &len in chr_lengths {
        if n > 0 {
            n = ((n + 1) / bin_size + 1) * bin_size;
        }
        chr_start.push(n);
        n += len;
    }
    n = ((n + 1) / bin_size + 1) * bin_size;
    chr_start.push(n);
    chr_start
}

impl Genome {
    /// Build a genome from FASTA files, matching STAR's layout.
    ///
    /// # Arguments
    /// - `params`: CLI parameters (genomeFastaFiles, genomeChrBinNbits)
    ///
    /// # Returns
    /// A `Genome` with forward + reverse complement sequences and metadata.
    pub fn from_fasta(params: &Parameters) -> Result<Self, Error> {
        let chromosomes = parse_fasta_files(&params.genome_fasta_files)?;
        let bin_nbits = params.genome_chr_bin_nbits;

        // First pass: chromosome names/lengths, validating non-zero length.
        let mut chr_name = Vec::new();
        let mut chr_length = Vec::new();
        for chrom in &chromosomes {
            let len = chrom.sequence.len() as u64;
            if len == 0 {
                return Err(Error::Fasta(format!(
                    "chromosome '{}' has zero length",
                    chrom.name
                )));
            }
            chr_name.push(chrom.name.clone());
            chr_length.push(len);
        }

        let chr_start = compute_chr_starts(&chr_length, bin_nbits);
        let n_genome = *chr_start.last().unwrap();
        let n_chr_real = chromosomes.len();

        // Allocate buffer: forward (0..n_genome) + reverse (n_genome..2*n_genome)
        let total_len = (n_genome * 2) as usize;
        let mut sequence = vec![GENOME_SPACING_CHAR; total_len];

        // Second pass: copy actual chromosome sequences into the buffer
        for (i, chrom) in chromosomes.iter().enumerate() {
            let start = chr_start[i] as usize;
            let len = chrom.sequence.len();
            sequence[start..start + len].copy_from_slice(&chrom.sequence);
        }

        // Build reverse complement in the second half
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        Ok(Genome {
            sequence: sequence.into(),
            n_genome,
            n_genome_real: n_genome,
            n_chr_real,
            chr_name,
            chr_length,
            chr_start,
        })
    }

    /// Append a splice-junction flanking-sequence buffer (`Gsj`) to the
    /// forward genome and rebuild the reverse complement over the extended
    /// forward range. Matches STAR's in-memory layout after
    /// `sjdbBuildIndex.cpp:293` (`memcpy(G+chrStart[nChrReal], Gsj, nGsj)`):
    /// the Gsj bytes live immediately after the forward real-genome bytes
    /// (`chr_start[n_chr_real]` stays pinned at the pre-sjdb forward total,
    /// matching STAR's `chrStart.txt` — verified against STAR docker output).
    ///
    /// `n_genome` grows to include Gsj. `chr_start` / `chr_length` /
    /// `chr_name` / `n_chr_real` are NOT updated — Gsj lives outside the
    /// chromosome accounting.
    pub fn append_sjdb(&mut self, gsj: &[u8]) {
        let old_n = self.n_genome;
        let new_n = old_n + gsj.len() as u64;

        let mut new_seq = vec![GENOME_SPACING_CHAR; (new_n * 2) as usize];
        new_seq[..old_n as usize].copy_from_slice(&self.sequence.as_slice()[..old_n as usize]);
        new_seq[old_n as usize..new_n as usize].copy_from_slice(gsj);

        // Rebuild RC over the extended forward range (STAR stores Gsj_RC
        // implicitly — rustar-aligner keeps the explicit `[fwd | RC]` layout).
        for i in 0..new_n as usize {
            let base = new_seq[i];
            let complement = if base < 4 { 3 - base } else { base };
            new_seq[2 * new_n as usize - 1 - i] = complement;
        }

        self.sequence = new_seq.into();
        self.n_genome = new_n;
    }

    /// Access a base from the genome (forward or reverse strand).
    ///
    /// # Arguments
    /// - `pos`: Position in the genome (0..2*n_genome)
    ///
    /// # Returns
    /// The base value (0-3 for ACGT, 4 for N, 5 for padding), or None if out of bounds.
    ///
    /// `#[inline]`: called per-byte in the hottest loops (`extend_alignment`'s DP
    /// scan, `cluster_seeds`) — uninlined, this is a real function-call boundary on
    /// every base, unlike STAR's plain `G[i]` array indexing it stands in for.
    #[inline]
    pub fn get_base(&self, pos: u64) -> Option<u8> {
        if pos < self.sequence.len() as u64 {
            Some(self.sequence.base(pos as usize))
        } else {
            None
        }
    }

    /// Get the chromosome containing a given genomic position.
    ///
    /// Uses binary search on chr_start for O(log n) lookup.
    ///
    /// # Returns
    /// `(chr_index, offset_within_chr)` or None if position is in padding.
    pub fn position_to_chr(&self, pos: u64) -> Option<(usize, u64)> {
        // Binary search: find the last chr_start that is <= pos
        let idx = self.chr_start[..self.n_chr_real].partition_point(|&start| start <= pos);
        if idx == 0 {
            return None;
        }
        let i = idx - 1;
        let start = self.chr_start[i];
        if pos < start + self.chr_length[i] {
            Some((i, pos - start))
        } else {
            None // Position is in padding between chromosomes
        }
    }

    /// Write the `Genome` section (forward-strand padded sequence, exactly
    /// `n_genome` bytes) through `sink`.
    pub fn write_genome_section(&self, sink: &mut dyn SectionSink) -> Result<(), Error> {
        sink.write_section("Genome", &self.sequence.as_slice()[..self.n_genome as usize])
    }

    /// Write the chromosome metadata text sections through `sink`:
    /// `chrName.txt` (names, one per line), `chrLength.txt` (lengths, one per
    /// line), and `chrStart.txt` (bin-padded start offsets plus the trailing
    /// `n_genome` entry).
    pub fn write_metadata_sections(&self, sink: &mut dyn SectionSink) -> Result<(), Error> {
        use std::io::Write;

        let mut buf = Vec::new();
        for &len in &self.chr_length {
            writeln!(buf, "{len}")?;
        }
        sink.write_section("chrLength.txt", &buf)?;

        buf.clear();
        for name in &self.chr_name {
            writeln!(buf, "{name}")?;
        }
        sink.write_section("chrName.txt", &buf)?;

        buf.clear();
        for &start in &self.chr_start {
            writeln!(buf, "{start}")?;
        }
        sink.write_section("chrStart.txt", &buf)?;
        Ok(())
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_params(fasta_paths: &[std::path::PathBuf], bin_nbits: u32) -> Parameters {
        let mut args = vec!["rustar-aligner", "--runMode", "genomeGenerate"];

        for path in fasta_paths {
            args.push("--genomeFastaFiles");
            args.push(path.to_str().unwrap());
        }

        let bin_str = bin_nbits.to_string();
        args.extend(["--genomeChrBinNbits", &bin_str]);

        Parameters::parse_from(args)
    }

    #[test]
    fn single_chromosome_padding() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap(); // 4 bases

        let params = make_params(&[file.path().to_path_buf()], 3); // bin_size = 8
        let genome = Genome::from_fasta(&params).unwrap();

        // Padding formula: n=4, then ((4+1)/8 + 1)*8 = (0+1)*8 = 8
        assert_eq!(genome.n_genome, 8);
        assert_eq!(genome.n_chr_real, 1);
        assert_eq!(genome.chr_start, vec![0, 8]);
        assert_eq!(genome.chr_length, vec![4]);

        // Check bases
        assert_eq!(genome.get_base(0), Some(0)); // A
        assert_eq!(genome.get_base(1), Some(1)); // C
        assert_eq!(genome.get_base(2), Some(2)); // G
        assert_eq!(genome.get_base(3), Some(3)); // T
        assert_eq!(genome.get_base(4), Some(5)); // padding
        assert_eq!(genome.get_base(7), Some(5)); // padding

        // Check reverse complement
        // Formula: pos(2*n - 1 - i) = complement(pos i)
        // For n=8: pos 15 = complement(pos 0), pos 12 = complement(pos 3), etc.
        assert_eq!(genome.get_base(15), Some(3)); // T (complement of A at pos 0)
        assert_eq!(genome.get_base(14), Some(2)); // G (complement of C at pos 1)
        assert_eq!(genome.get_base(13), Some(1)); // C (complement of G at pos 2)
        assert_eq!(genome.get_base(12), Some(0)); // A (complement of T at pos 3)
        assert_eq!(genome.get_base(11), Some(5)); // padding (complement of padding at pos 4)
        assert_eq!(genome.get_base(8), Some(5)); // padding (complement of padding at pos 7)
    }

    #[test]
    fn two_chromosomes_padding() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "AA").unwrap(); // 2 bases
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TT").unwrap(); // 2 bases

        let params = make_params(&[file.path().to_path_buf()], 2); // bin_size = 4
        let genome = Genome::from_fasta(&params).unwrap();

        // chr1 starts at 0, length 2
        // After chr1: n=2, padding ((2+1)/4 + 1)*4 = (0+1)*4 = 4
        // chr2 starts at 4, length 2
        // After chr2: n=6, padding ((6+1)/4 + 1)*4 = (1+1)*4 = 8
        assert_eq!(genome.n_genome, 8);
        assert_eq!(genome.chr_start, vec![0, 4, 8]);
        assert_eq!(genome.chr_length, vec![2, 2]);

        // chr1 bases
        assert_eq!(genome.get_base(0), Some(0)); // A
        assert_eq!(genome.get_base(1), Some(0)); // A
        assert_eq!(genome.get_base(2), Some(5)); // padding
        assert_eq!(genome.get_base(3), Some(5)); // padding

        // chr2 bases
        assert_eq!(genome.get_base(4), Some(3)); // T
        assert_eq!(genome.get_base(5), Some(3)); // T
        assert_eq!(genome.get_base(6), Some(5)); // padding
        assert_eq!(genome.get_base(7), Some(5)); // padding
    }

    #[test]
    fn reverse_complement_correctness() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">test").unwrap();
        writeln!(file, "ACGTN").unwrap(); // A=0, C=1, G=2, T=3, N=4

        let params = make_params(&[file.path().to_path_buf()], 3); // bin_size = 8
        let genome = Genome::from_fasta(&params).unwrap();

        let n = genome.n_genome as usize;

        // Forward: A C G T N (then padding)
        assert_eq!(genome.sequence.base(0), 0); // A
        assert_eq!(genome.sequence.base(1), 1); // C
        assert_eq!(genome.sequence.base(2), 2); // G
        assert_eq!(genome.sequence.base(3), 3); // T
        assert_eq!(genome.sequence.base(4), 4); // N

        // Reverse complement should be at positions [2n-1, 2n-2, 2n-3, 2n-4, 2n-5]
        // which maps to the reverse of [0,1,2,3,4]
        assert_eq!(genome.sequence.base(2 * n - 1), 3); // T (complement of A at pos 0)
        assert_eq!(genome.sequence.base(2 * n - 1 - 1), 2); // G (complement of C at pos 1)
        assert_eq!(genome.sequence.base(2 * n - 1 - 2), 1); // C (complement of G at pos 2)
        assert_eq!(genome.sequence.base(2 * n - 1 - 3), 0); // A (complement of T at pos 3)
        assert_eq!(genome.sequence.base(2 * n - 1 - 4), 4); // N (complement of N at pos 4)
    }

    #[test]
    fn position_to_chr_mapping() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "AAA").unwrap();
        writeln!(file, ">chr2").unwrap();
        writeln!(file, "TTT").unwrap();

        let params = make_params(&[file.path().to_path_buf()], 2); // bin_size = 4
        let genome = Genome::from_fasta(&params).unwrap();

        // chr1: positions 0-2 (starts at 0, length 3)
        // After chr1: n=3, padding ((3+1)/4 + 1)*4 = (1+1)*4 = 8
        // chr2: positions 8-10 (starts at 8, length 3)
        // After chr2: n=11, padding ((11+1)/4 + 1)*4 = (3+1)*4 = 16

        assert_eq!(genome.n_genome, 16);
        assert_eq!(genome.chr_start, vec![0, 8, 16]);

        assert_eq!(genome.position_to_chr(0), Some((0, 0)));
        assert_eq!(genome.position_to_chr(1), Some((0, 1)));
        assert_eq!(genome.position_to_chr(2), Some((0, 2)));
        assert_eq!(genome.position_to_chr(3), None); // padding
        assert_eq!(genome.position_to_chr(7), None); // padding
        assert_eq!(genome.position_to_chr(8), Some((1, 0)));
        assert_eq!(genome.position_to_chr(9), Some((1, 1)));
        assert_eq!(genome.position_to_chr(10), Some((1, 2)));
        assert_eq!(genome.position_to_chr(11), None); // padding
    }

    #[test]
    fn append_sjdb_extends_forward_and_rebuilds_rc() {
        // Forward "ACGT" + spacer padding (bin 3 → 8 bytes padded).
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        let params = make_params(&[file.path().to_path_buf()], 3);
        let mut genome = Genome::from_fasta(&params).unwrap();
        assert_eq!(genome.n_genome, 8);

        // Gsj bytes: donor "AC" + acceptor "GT" + spacer 5 (five total).
        let gsj: Vec<u8> = vec![0, 1, 2, 3, 5];
        let gsj_len = gsj.len();
        genome.append_sjdb(&gsj);

        // n_genome grows by gsj_len; chr_start pinned at pre-sjdb boundary.
        assert_eq!(genome.n_genome, 8 + gsj_len as u64);
        assert_eq!(genome.chr_start, vec![0, 8]);
        assert_eq!(genome.n_chr_real, 1);

        // Forward is [real 0..8 | gsj 8..13].
        assert_eq!(&genome.sequence.as_slice()[..4], &[0, 1, 2, 3]);
        assert_eq!(&genome.sequence.as_slice()[8..13], gsj.as_slice());

        // RC over the extended forward range. sequence[2n-1-i] = complement(sequence[i]).
        let new_n = genome.n_genome as usize;
        assert_eq!(genome.sequence.base(2 * new_n - 1 - 8), 3); // complement of A at fwd[8]=0
        assert_eq!(genome.sequence.base(2 * new_n - 1 - 12), 5); // spacer stays 5
        assert_eq!(genome.sequence.len(), 2 * new_n);
    }

    #[test]
    fn append_sjdb_with_empty_gsj_is_noop() {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "ACGT").unwrap();
        let params = make_params(&[file.path().to_path_buf()], 3);
        let mut genome = Genome::from_fasta(&params).unwrap();
        let before = genome.sequence.clone();
        let before_n = genome.n_genome;
        genome.append_sjdb(&[]);
        assert_eq!(genome.n_genome, before_n);
        assert_eq!(genome.sequence, before);
    }
}
