pub mod packed_array;
pub mod packed_stream;
pub mod sa_build;
pub mod sa_index;
pub mod suffix_array;

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::star::error::Error;
use crate::star::genome::Genome;
use crate::star::index::packed_stream::PackedStreamWriter;
use crate::star::junction::SpliceJunctionDb;
use crate::star::junction::sjdb_insert::{self, PreparedJunction};
use crate::star::params::Parameters;
use sa_index::SaIndex;
use suffix_array::SuffixArray;

/// Destination for the index sections produced by
/// [`GenomeIndex::generate_streaming_with`]. Implementations decide the
/// physical layout: seqarc streams the sections straight into its packed
/// single-file container (alignment padding between sections, table of
/// contents back-patched at the end).
pub trait SectionSink {
    /// Announce the sections that will be written, before any data.
    /// `sjdbInfo.txt` is only present when junctions were prepared.
    fn plan(&mut self, names: &[&str]) -> Result<(), Error>;

    /// Write a section from an in-memory buffer.
    fn write_section(&mut self, name: &str, data: &[u8]) -> Result<(), Error>;

    /// Stream a section through the sink's file. The closure receives the
    /// file positioned at the section's first byte, must leave the cursor
    /// at the section's last written byte, and returns the section's data
    /// length. The sink records `(offset, length)` and pads to its
    /// alignment before the next section.
    fn stream_section(
        &mut self,
        name: &str,
        stream: &mut dyn FnMut(&mut File) -> Result<u64, Error>,
    ) -> Result<(), Error>;

    /// A file handle plus the absolute offset of a written section's first
    /// byte, for chunked `pread` passes over large streamed sections (the
    /// SAindex build reads the just-written SA this way).
    fn section_reader(&self, name: &str) -> Result<(File, u64), Error>;
}

/// Map an I/O failure inside a section writer.
fn section_io(e: std::io::Error) -> Error {
    Error::Index(format!("index section write failed: {e}"))
}

/// Complete genome index (genome + suffix array + SA index + junction database).
#[derive(Clone)]
pub struct GenomeIndex {
    pub genome: Genome,
    pub suffix_array: SuffixArray,
    pub sa_index: SaIndex,
    pub junction_db: SpliceJunctionDb,
    /// Prepared splice junctions in their post-dedup order (the same
    /// order they occupy in the Gsj buffer appended to the genome).
    /// Populated on both build and load paths when sjdb is present;
    /// empty for indices built without a GTF. Used at align time to
    /// decode Gsj-region SA hits back to real-genome `(donor, acceptor)`
    /// pairs.
    pub prepared_junctions: Vec<PreparedJunction>,
    /// `sjdbOverhang` recorded in `sjdbInfo.txt`. Zero when no sjdb
    /// junctions are present.
    pub sjdb_overhang: u32,
}

/// Output of [`GenomeIndex::build_prep`] — the shared setup
/// stage for the streaming [`GenomeIndex::generate_streaming_with`].
/// Holds everything needed to drive the SA construction and the
/// per-file writes.
struct BuildPrep {
    genome: Genome,
    prepared_junctions: Vec<PreparedJunction>,
}

impl GenomeIndex {
    /// Streaming genome-index generation: hands every index section to
    /// `sink` without materialising the 25 GB-class SA `PackedArray` in
    /// RAM. The flow:
    ///
    /// 1. [`build_prep`][Self::build_prep] — load FASTA, parse GTF,
    ///    build junction database, append `gsj` to the genome.
    /// 2. Announce the section plan to the sink (`sjdbInfo.txt` only when
    ///    junctions exist) and write the `Genome` section — no dependency
    ///    on the SA.
    /// 3. Stream the SA through [`SectionSink::stream_section`]. The
    ///    caps-sa emit callback bit-packs each entry straight into the
    ///    sink's file, so it grows as construction progresses. Total
    ///    peak RSS during this phase ≈
    ///    `genome.sequence` + caps-sa scratch —
    ///    ~5 GB on the human genome vs the ~47 GB an in-memory
    ///    SA build would peak at.
    /// 4. Finalise the SA writer (flush partial-byte + padding).
    /// 5. **Build the SAindex in parallel** from the just-written SA via
    ///    chunked pread + [`SaIndex::build_parallel`]. caps-sa's phase-4
    ///    emit is single-threaded; doing the SAindex k-mer work
    ///    there sat the entire ~16 min of serial extraction
    ///    on top of caps-sa's parallel SA build. Deferring lets
    ///    the SAindex use all available rayon workers.
    /// 6. Stream the SAindex section, then drop it — the sjdb writers
    ///    below don't read it. On the human genome this frees ~1.5 GB
    ///    of resident `PackedArray` before the last writes.
    /// 7. Write the chromosome metadata sections and the sjdb output.
    pub fn generate_streaming_with(
        params: &Parameters,
        sink: &mut dyn SectionSink,
    ) -> Result<(), Error> {
        let BuildPrep {
            genome,
            prepared_junctions,
        } = Self::build_prep(params)?;

        // The sink must know the full section list before any data is
        // written (seqarc reserves the container table of contents up
        // front). Order matches the container layout.
        let mut names = vec![
            "Genome",
            "SA",
            "SAindex",
            "chrLength.txt",
            "chrName.txt",
            "chrStart.txt",
        ];
        if !prepared_junctions.is_empty() {
            names.push("sjdbInfo.txt");
        }
        sink.plan(&names)?;

        log::info!("Writing genome section...");
        genome.write_genome_section(sink)?;

        let gstrand_bit = SuffixArray::calculate_gstrand_bit(genome.n_genome);
        let gstrand_mask = (1u64 << gstrand_bit) - 1;
        let word_length = gstrand_bit + 1;
        let nbases = params.genome_sa_index_nbases;

        log::info!(
            "Streaming SA (gstrand_bit={gstrand_bit}, word_length={word_length}, nbases={nbases})"
        );

        let mut n_entries = 0usize;
        sink.stream_section("SA", &mut |file| {
            // Large buffer — at 33 bit / entry × 5.9 B entries the SA
            // is ~24 GB; an 8 MB BufWriter keeps the write-side syscall
            // rate to ~3000/s, dominated by sequential bandwidth.
            let mut sa_buf = BufWriter::with_capacity(8 * 1024 * 1024, file);
            {
                let mut sa_writer = PackedStreamWriter::new(&mut sa_buf, word_length);
                log::info!("Building suffix array...");
                let (got_gbit, got_gmask, got_entries) = sa_build::build_streaming(
                    &genome,
                    params.temp_dir.as_deref(),
                    params.genome_sa_sparse_d as u64,
                    |packed_value| {
                        // Emit is now lightweight: just bit-pack into the
                        // SA. caps-sa's phase-4 emit loop is
                        // single-threaded, so anything we do here
                        // serialises the whole build. The SAindex is built
                        // afterwards via a parallel pass over the on-disk
                        // SA.
                        sa_writer
                            .write_one(packed_value)
                            .map_err(section_io)?;
                        Ok(())
                    },
                )?;
                debug_assert_eq!(got_gbit, gstrand_bit);
                debug_assert_eq!(got_gmask, gstrand_mask);
                n_entries = got_entries;
                log::info!("Suffix array streamed to disk: {n_entries} entries");
                sa_writer.finish().map_err(section_io)?;
            }
            sa_buf.flush().map_err(section_io)?;
            Ok(crate::star::index::packed_array::PackedArray::data_byte_len_for(
                word_length, n_entries,
            ) as u64)
        })?;

        // Parallel SAindex build via chunked pread on the just-written
        // SA. Avoiding `memmap2::Mmap` keeps the touched SA pages
        // out of process RSS (`read_at` hits the kernel page cache,
        // which is kernel-side memory — not counted in
        // `Maximum resident set size`). The only big new allocation
        // for this phase is the ~2.86 GB `Vec<AtomicU64>` inside
        // `build_parallel`.
        log::info!("Building SAindex in parallel from the streamed SA...");
        let (sa_file, sa_base) = sink.section_reader("SA")?;
        let sa_index = SaIndex::build_parallel(
            &genome,
            &sa_file,
            sa_base,
            word_length,
            gstrand_bit,
            gstrand_mask,
            n_entries,
            nbases,
        )?;
        drop(sa_file);
        log::info!(
            "SA index built: nbases={}, {} indices",
            sa_index.nbases,
            sa_index.data.len()
        );

        // Stream the SAindex section, then drop it — the sjdb writers
        // below don't read it. On the human genome this frees ~1.5 GB of
        // resident `PackedArray` before the last writes.
        sink.stream_section("SAindex", &mut |file| {
            let mut w = BufWriter::with_capacity(8 * 1024 * 1024, file);
            w.write_all(&(sa_index.nbases as u64).to_le_bytes())
                .map_err(section_io)?;
            for &val in &sa_index.genome_sa_index_start {
                w.write_all(&val.to_le_bytes()).map_err(section_io)?;
            }
            w.write_all(sa_index.data.data()).map_err(section_io)?;
            w.flush().map_err(section_io)?;
            Ok(
                8 + 8 * (sa_index.nbases as u64 + 1)
                    + sa_index.data.data().len() as u64,
            )
        })?;
        drop(sa_index);

        // Chromosome metadata sections (container order is alphabetical,
        // matching the layout the old directory-packing produced).
        genome.write_metadata_sections(sink)?;

        // Sjdb section.
        if !prepared_junctions.is_empty() {
            let bytes =
                sjdb_insert::sjdb_info_bytes(&prepared_junctions, params.sjdb_overhang)?;
            sink.write_section("sjdbInfo.txt", &bytes)?;
            log::info!("Wrote sjdb section: {} junctions", prepared_junctions.len());
        }

        Ok(())
    }

    /// Shared setup: load FASTA, parse GTF, build junctions, append
    /// `gsj` to the genome. Used by
    /// [`generate_streaming_with`][Self::generate_streaming_with] so the
    /// SA-input shape is built in one place.
    fn build_prep(params: &Parameters) -> Result<BuildPrep, Error> {
        log::info!("Loading FASTA files...");
        let mut genome = Genome::from_fasta(params)?;

        log::info!(
            "Loaded {} chromosomes, total padded genome size: {} bytes",
            genome.n_chr_real,
            genome.n_genome
        );

        // Parse GTF once and share the result between the junction database
        // and the sjdb insertion pipeline. Junction preparation + Gsj append
        // must happen BEFORE the suffix array is built, because STAR indexes
        // the flanking-sequence buffer alongside the real genome in a single
        // SA (`sjdbBuildIndex.cpp:293`).
        let n_genome_real = genome.n_genome;

        let mut raw: Vec<(usize, u64, u64, u8)> = Vec::new();

        if let Some(ref gtf_path) = params.sjdb_gtf_file {
            let exons = crate::star::junction::gtf::parse_gtf_configured(
                gtf_path,
                &params.sjdb_gtf_feature_exon,
                &params.sjdb_gtf_chr_prefix,
            )?;
            log::debug!("Parsed {} exon features from GTF", exons.len());

            let gtf_raw = crate::star::junction::gtf::extract_junctions_configured(
                exons,
                &genome,
                &params.sjdb_gtf_tag_exon_parent_transcript,
            )?;
            log::info!("Extracted {} annotated junctions from GTF", gtf_raw.len());
            raw.extend(gtf_raw);
        }

        if !params.sjdb_file_chr_start_end.is_empty() {
            let extra = crate::star::junction::chr_start_end::parse_sjdb_chr_start_end(
                &params.sjdb_file_chr_start_end,
                &genome,
            )?;
            log::info!(
                "Parsed {} junctions from --sjdbFileChrStartEnd",
                extra.len()
            );
            raw.extend(extra);
        }

        let (junction_db, prepared_junctions) = if !raw.is_empty() {
            let jdb = SpliceJunctionDb::from_raw_junctions(&raw);

            let prepared: Vec<PreparedJunction> = raw
                .iter()
                .map(|&(chr_idx, intron_start, intron_end, strand)| {
                    sjdb_insert::prepare_junction(
                        chr_idx,
                        intron_start,
                        intron_end,
                        strand,
                        &genome,
                        n_genome_real,
                    )
                })
                .collect();
            let prepared = sjdb_insert::sort_and_dedup(prepared);

            let gsj =
                sjdb_insert::build_gsj(&prepared, &genome, n_genome_real, params.sjdb_overhang)?;
            log::info!(
                "Built Gsj buffer: {} junctions × {} bytes = {} bytes",
                prepared.len(),
                2 * params.sjdb_overhang + 1,
                gsj.len()
            );
            genome.append_sjdb(&gsj);
            log::info!(
                "Extended genome with sjdb: n_genome = {} (pre-sjdb {})",
                genome.n_genome,
                n_genome_real
            );

            (jdb, prepared)
        } else {
            log::info!("No GTF or --sjdbFileChrStartEnd provided, all junctions will be novel");
            (SpliceJunctionDb::empty(), Vec::new())
        };

        log::info!(
            "Junction database initialized: {} annotated junctions",
            junction_db.len()
        );

        Ok(BuildPrep {
            genome,
            prepared_junctions,
        })
    }

    /// Convert a raw SA position for a reverse-strand match to forward genome coordinates.
    ///
    /// The SA stores reverse-strand positions as offsets within the RC genome region.
    /// For chromosome identification and SAM output, we need the corresponding position
    /// in the forward genome (leftmost aligned base in forward coordinates).
    ///
    /// For forward strand: returns the position unchanged.
    /// For reverse strand: `forward_pos = n_genome - 1 - sa_pos - (match_length - 1)`
    ///
    /// The raw SA position is still needed for genome base access (add n_genome offset).
    pub fn sa_pos_to_forward(&self, sa_pos: u64, is_reverse: bool, match_length: usize) -> u64 {
        if is_reverse {
            self.genome
                .n_genome
                .saturating_sub(sa_pos)
                .saturating_sub(match_length as u64)
        } else {
            sa_pos
        }
    }
}
