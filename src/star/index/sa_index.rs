use std::fs::File;
use std::sync::atomic::{AtomicU64, Ordering};

use rayon::prelude::*;

use crate::star::error::Error;
use crate::star::genome::Genome;
use crate::star::index::packed_array::PackedArray;

/// Thread-safe positioned read (no shared cursor), portable across
/// platforms: `pread` on Unix via [`std::os::unix::fs::FileExt::read_at`],
/// the equivalent overlapped read on Windows via
/// [`std::os::windows::fs::FileExt::seek_read`]. Like both primitives it
/// may return a short count at EOF; callers size their buffer with tail
/// padding so a short final read just leaves the padding bytes zero.
#[cfg(unix)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
    use std::os::unix::fs::FileExt;
    file.read_at(buf, offset)
}

#[cfg(windows)]
fn read_at(file: &File, buf: &mut [u8], offset: u64) -> std::io::Result<usize> {
    use std::os::windows::fs::FileExt;
    file.seek_read(buf, offset)
}

/// SA index for fast k-mer lookup during binary search.
///
/// For each k-mer prefix (up to genomeSAindexNbases length), stores the
/// SA range where suffixes starting with that prefix can be found.
#[derive(Clone)]
pub struct SaIndex {
    /// Length of indexed k-mers (typically 14)
    pub nbases: u32,

    /// Cumulative counts for each k-mer level
    /// Length = nbases + 1
    /// genomeSAindexStart[k] = total number of k-mers of length <= k
    pub genome_sa_index_start: Vec<u64>,

    /// Packed array storing SA index entries
    /// Each entry is (gstrand_bit + 3) bits:
    /// - Bits 0..gstrand_bit: SA position (or max value if prefix absent)
    /// - Bit gstrand_bit+1: "contains N" flag
    /// - Bit gstrand_bit+2: "prefix absent" flag
    pub data: PackedArray,

    /// Strand bit position (from SA)
    pub gstrand_bit: u32,
}

impl SaIndex {
    /// Calculate the total number of k-mer prefixes to index.
    ///
    /// Returns sum of 4^1 + 4^2 + ... + 4^nbases = (4^(nbases+1) - 4) / 3
    pub fn calculate_num_indices(nbases: u32) -> u64 {
        if nbases == 0 {
            return 0;
        }
        // Sum of geometric series: (4^(n+1) - 4) / 3
        let power = 4u64.pow(nbases + 1);
        (power - 4) / 3
    }

    /// Parallel SAindex build from the on-disk SA file, using STAR's
    /// `isaStep + binary search` skip algorithm inside each worker.
    ///
    /// ## Algorithm — STAR's `genomeSAindex.cpp::genomeSAindexChunk`
    ///
    /// The SA is sorted lex, so consecutive `SA[isa]` entries share
    /// monotonically non-decreasing k-mer prefixes. STAR walks the
    /// SA with two ideas:
    ///
    /// 1. **Boundaries only**: keep `ind0[iL]` = the last-written
    ///    k-mer index at each level. When the k-mer at `SA[isa]`
    ///    has `indPref[iL] > ind0[iL]`, record `isa` at slot
    ///    `start[iL] + indPref[iL]` and update `ind0[iL]`.
    /// 2. **Skip identical runs**: `isaStep = nSA / 4^nbases` (≈ 22
    ///    on the human genome). Jump forward by `isaStep`; if the
    ///    full k-mer didn't change, keep jumping; if it did,
    ///    binary-search the boundary inside the last `isaStep`
    ///    window. Total visited entries ≈ number of distinct
    ///    k-mers × `O(log isaStep)`, not `nSA × nbases` — typically
    ///    **20-30× less work** than scanning every entry.
    ///
    /// ## Parallelisation
    ///
    /// Each worker owns a chunk `[chunk_start, chunk_end)` of SA
    /// indices and runs STAR's algorithm *locally* over its range
    /// with its own `ind0_local[iL]`. Cross-chunk merge uses
    /// `AtomicU64::fetch_min` on a shared `Vec<AtomicU64>` of length
    /// `num_indices`: the smallest `sa_idx` seen across all chunks
    /// wins.
    ///
    /// Chunk boundary handling: a chunk's `ind0_local[iL]` starts at
    /// `None` (no kmer written *in this chunk*). The first iteration
    /// always writes — the resulting slot may also be written by the
    /// previous chunk's last iteration, but `fetch_min` picks the
    /// earlier `sa_idx` so the answer is the same as the serial
    /// algorithm.
    ///
    /// ## Phase 2: gap-fill with `next_isa | absent_mask`
    ///
    /// After the parallel pass, `firsts[]` has present slots set to
    /// the first-occurrence `sa_idx` and absent slots still at
    /// `u64::MAX`. A sequential **backward** pass per SAindex level
    /// (level = `[start[iL], start[iL] + 4^(iL+1))`) replaces each
    /// absent slot with `next_present_sa_idx | absent_mask`,
    /// matching STAR's encoding. Tail-gaps (after the last present
    /// slot of a level) get `n_entries | absent_mask`.
    ///
    /// ## Memory
    ///
    /// - `firsts: Vec<AtomicU64>`: ~2.86 GB transient (358 M × 8 B
    ///   on the human genome). Reused in-place during gap-fill, then
    ///   dropped after encoding into the output `PackedArray`.
    /// - Per-worker pread buffer: ~4 MB.
    /// - Output `PackedArray`: ~1.5 GB (kept).
    /// - SA bytes: kernel page cache (chunked `pread`), not process RSS.
    pub fn build_parallel(
        genome: &Genome,
        sa_file: &File,
        sa_base: u64,
        sa_word_length: u32,
        gstrand_bit: u32,
        gstrand_mask: u64,
        n_entries: usize,
        nbases: u32,
    ) -> Result<Self, Error> {
        let sai_word_length = gstrand_bit + 3;
        let mut genome_sa_index_start = vec![0u64; (nbases + 1) as usize];
        for k in 1..=nbases {
            genome_sa_index_start[k as usize] =
                genome_sa_index_start[(k - 1) as usize] + 4u64.pow(k);
        }
        let num_indices = Self::calculate_num_indices(nbases) as usize;
        let absent_mask: u64 = 1u64 << (gstrand_bit + 2);

        // `isaStep` from STAR: `nSA / 4^nbases`. With 5.9 B SA entries
        // and `nbases = 14`, this is ~22 — i.e. on average we expect
        // ~22 consecutive SA entries to share their full 14-mer (and
        // can be skipped over). Clamp to ≥ 1 to keep the loop
        // well-formed when `nSA < 4^nbases`.
        let isa_step = (n_entries / (1usize << (2 * nbases as usize))).max(1);
        log::info!(
            "Building SA index in parallel: nbases={nbases}, num_indices={num_indices}, \
             n_entries={n_entries}, isa_step={isa_step}, threads={}",
            rayon::current_num_threads()
        );

        // Sentinel `u64::MAX` for "no `sa_idx` has reached this slot
        // yet". `fetch_min` makes any real `sa_idx` smaller than the
        // sentinel, so the first chunk to write wins — and across
        // chunks, the smallest `sa_idx` wins.
        let firsts: Vec<AtomicU64> = (0..num_indices)
            .into_par_iter()
            .map(|_| AtomicU64::new(u64::MAX))
            .collect();

        let sa_mask: u64 = if sa_word_length == 64 {
            u64::MAX
        } else {
            (1u64 << sa_word_length) - 1
        };
        let n_genome = genome.n_genome as usize;
        let genome_seq: &[u8] = genome.sequence.as_slice();

        // Chunk size: 1 M entries per worker. STAR's algorithm
        // visits at most ~chunk_size / isa_step boundaries per chunk
        // (plus `log isa_step` binary-search probes each), and each
        // boundary touches `genome.sequence` once. Larger chunks
        // grow the per-worker pread buffer (~4 MB at 33-bit packed
        // words for 1 M entries); smaller chunks amortise the pread
        // setup over fewer boundaries.
        const ENTRIES_PER_CHUNK: usize = 1 << 20;
        let n_chunks = n_entries.div_ceil(ENTRIES_PER_CHUNK);

        (0..n_chunks)
            .into_par_iter()
            .try_for_each(|chunk_idx| -> std::io::Result<()> {
                let chunk_start = chunk_idx * ENTRIES_PER_CHUNK;
                let chunk_end = (chunk_start + ENTRIES_PER_CHUNK).min(n_entries);
                let chunk_n = chunk_end - chunk_start;
                if chunk_n == 0 {
                    return Ok(());
                }

                // Pre-read the chunk's SA byte range. `read_at` is
                // thread-safe pread; the +8-byte tail padding makes
                // every entry's 8-byte LE load safe.
                let start_bit = chunk_start as u64 * sa_word_length as u64;
                let end_bit = chunk_end as u64 * sa_word_length as u64;
                let start_byte = start_bit / 8;
                let start_bit_shift = (start_bit % 8) as u32;
                let end_byte_excl = end_bit.div_ceil(8);
                let bytes_to_read = (end_byte_excl - start_byte) as usize + 8;
                let mut buf = vec![0u8; bytes_to_read];
                let _ = read_at(sa_file, &mut buf, sa_base + start_byte)?;

                // Read packed value at chunk-local index `i`
                // (`0..chunk_n`).
                let read_packed = |i: usize| -> u64 {
                    let local_bit = i as u64 * sa_word_length as u64 + start_bit_shift as u64;
                    let local_byte = (local_bit / 8) as usize;
                    let local_shift = (local_bit % 8) as u32;
                    let bytes: &[u8; 8] = (&buf[local_byte..local_byte + 8])
                        .try_into()
                        .expect("padded buf must have 8 bytes at every entry's start");
                    let word = u64::from_le_bytes(*bytes);
                    (word >> local_shift) & sa_mask
                };

                // STAR's `funCalcSAiFromSA`: compute the full
                // `nbases`-long k-mer at chunk-local index `i`, plus
                // the level `iL4` where the first N (if any) appears.
                // If no N: `iL4 = -1`, the k-mer is fully defined.
                // If N at level `iL4`: the returned k-mer has zeros
                // in positions `iL4..nbases` (i.e. only the
                // `iL4`-long prefix is valid).
                let calc_kmer = |i: usize| -> (u64, i32) {
                    let packed = read_packed(i);
                    let pos = packed & gstrand_mask;
                    let is_reverse = (packed >> gstrand_bit) != 0;
                    let genome_pos = if is_reverse {
                        pos as usize + n_genome
                    } else {
                        pos as usize
                    };
                    let mut kmer: u64 = 0;
                    for ii in 0..nbases as usize {
                        if genome_pos + ii >= genome_seq.len() {
                            // Treat past-end as an N for early-break
                            // purposes. Suffix has fewer than nbases
                            // comparable bases.
                            return (kmer << (2 * (nbases as usize - ii)), ii as i32);
                        }
                        let g = genome_seq[genome_pos + ii];
                        if g >= 4 {
                            return (kmer << (2 * (nbases as usize - ii)), ii as i32);
                        }
                        kmer = (kmer << 2) | (g as u64);
                    }
                    (kmer, -1)
                };

                // STAR's `funSAiFindNextIndex`: jump forward by
                // `isa_step` while `(indFull, iL4)` is unchanged;
                // when it changes, binary-search the boundary inside
                // the last `isa_step` window. Returns the chunk-local
                // index of the first entry where `(indFull, iL4)`
                // differs from the input, plus its `(indFull, iL4)`.
                // Returns `chunk_n` if no change is found in this
                // chunk.
                let find_next =
                    |i: usize, ind_full_prev: u64, il4_prev: i32| -> (usize, u64, i32) {
                        let mut next_i = i + isa_step;
                        let mut next_kmer = 0u64;
                        let mut next_il4 = -1i32;
                        while next_i < chunk_n {
                            let (k, l) = calc_kmer(next_i);
                            if k == ind_full_prev && l == il4_prev {
                                next_i += isa_step;
                                continue;
                            }
                            next_kmer = k;
                            next_il4 = l;
                            break;
                        }
                        if next_i >= chunk_n {
                            // Past chunk end. Check the chunk's last
                            // entry; if it still matches, no boundary
                            // exists in this chunk.
                            if chunk_n > 0 {
                                let last_i = chunk_n - 1;
                                let (k, l) = calc_kmer(last_i);
                                if k == ind_full_prev && l == il4_prev {
                                    return (chunk_n, 0, -1);
                                }
                                next_kmer = k;
                                next_il4 = l;
                                next_i = last_i;
                            } else {
                                return (chunk_n, 0, -1);
                            }
                        }
                        // Binary search in `(i .. next_i]` for the first
                        // index where `(indFull, iL4)` differs from
                        // `(ind_full_prev, il4_prev)`.
                        let mut lo = next_i.saturating_sub(isa_step).max(i);
                        let mut hi = next_i.min(chunk_n - 1);
                        while lo + 1 < hi {
                            let mid = lo + (hi - lo) / 2;
                            let (k, l) = calc_kmer(mid);
                            if k == ind_full_prev && l == il4_prev {
                                lo = mid;
                            } else {
                                hi = mid;
                                next_kmer = k;
                                next_il4 = l;
                            }
                        }
                        (hi, next_kmer, next_il4)
                    };

                // Per-chunk last-written kmer index at each level.
                // `None` means "nothing written in this chunk yet";
                // the first iteration always writes (matches STAR's
                // `isa == 0` special-case).
                let mut ind0_local: [Option<u64>; 32] = [None; 32];
                let mut i: usize = 0;
                let (mut ind_full, mut il4) = calc_kmer(i);

                while i < chunk_n {
                    let sa_idx = (chunk_start + i) as u64;
                    for il in 0..nbases as usize {
                        if il as i32 == il4 {
                            // N at level `il`. STAR sets the N flag
                            // on `ind0[il1]` for `il1 >= il4`; we
                            // don't track the N flag in this version
                            // (our `hierarchical_lookup` doesn't
                            // consult it), so just break out of the
                            // level loop. This means our SAindex
                            // file is not byte-identical to STAR's
                            // in the N-bit positions — documented
                            // limitation.
                            break;
                        }
                        let ind_pref = ind_full >> (2 * (nbases as usize - 1 - il));
                        let is_new = match ind0_local[il] {
                            None => true,
                            Some(prev) => ind_pref > prev,
                        };
                        if is_new {
                            let slot = (genome_sa_index_start[il] + ind_pref) as usize;
                            firsts[slot].fetch_min(sa_idx, Ordering::Relaxed);
                            ind0_local[il] = Some(ind_pref);
                        }
                        // `ind_pref < ind0_local[il]` would mean the
                        // SA isn't sorted — we don't error here
                        // (per-chunk locality + sorted SA → can't
                        // happen unless caps-sa produced a broken
                        // SA, which our differential tests guard
                        // against).
                    }
                    let (next_i, next_full, next_il4) = find_next(i, ind_full, il4);
                    i = next_i;
                    ind_full = next_full;
                    il4 = next_il4;
                }
                Ok(())
            })
            .map_err(|e| Error::Index(format!("pread during SAindex build: {e}")))?;

        // Phase 2: backward gap-fill **in place** in `firsts[]`.
        // For each SAindex level, walk slots from high to low. Track
        // the most recently seen present `sa_idx`; replace each
        // `u64::MAX` slot with `next_present_sa_idx | absent_mask`.
        // Levels are processed independently (no dependency between
        // them); parallel across levels is overkill for ~358 M
        // total slots × ~10 ns per slot = ~3.6 s.
        for (il, &level_start_raw) in genome_sa_index_start
            .iter()
            .enumerate()
            .take(nbases as usize)
        {
            let level_start = level_start_raw as usize;
            let level_size = 4u64.pow(il as u32 + 1) as usize;
            // Tail-gap sentinel: STAR uses `nSA | absent_mask` for
            // slots after the last-written present k-mer at this
            // level (i.e. the trailing run that has no further
            // present slot to point at).
            let mut next_present: u64 = n_entries as u64;
            for off in (0..level_size).rev() {
                let slot = level_start + off;
                let v = firsts[slot].load(Ordering::Relaxed);
                if v == u64::MAX {
                    firsts[slot].store(next_present | absent_mask, Ordering::Relaxed);
                } else {
                    next_present = v;
                }
            }
        }

        // Final sequential pack into the output `PackedArray`.
        // Every slot in `firsts` is now valid (either the
        // first-occurrence `sa_idx` or `next | absent_mask`).
        let mut data = PackedArray::new(sai_word_length, num_indices);
        for (i, slot) in firsts.iter().enumerate() {
            data.write(i, slot.load(Ordering::Relaxed));
        }
        drop(firsts);

        Ok(SaIndex {
            nbases,
            genome_sa_index_start,
            data,
            gstrand_bit,
        })
    }
}

impl SaIndex {
    /// Hierarchical SAindex lookup (STAR's maxMappableLength2strands approach).
    ///
    /// Starts with full k-mer at level min(k, nbases), progressively shortens
    /// until a present entry is found. Returns narrowed SA range + matched level.
    ///
    /// Returns: Some((sa_start, sa_end_exclusive, matched_level, bounds_tight))
    ///   - sa_start: first SA index in range
    ///   - sa_end_exclusive: past-the-end SA index
    ///   - matched_level: how many bases the SAindex resolved
    ///   - bounds_tight: both bounds came from present SAindex entries
    ///     (safe to skip first matched_level bases in binary search)
    ///
    /// Returns None if no prefix exists in the index (all levels absent).
    pub fn hierarchical_lookup(
        &self,
        kmer_idx: u64,
        k: u32,
        n_sa: usize,
    ) -> Option<(usize, usize, usize, bool)> {
        let mut lind = k.min(self.nbases);
        let mut ind = kmer_idx;

        // If k > nbases, truncate to nbases by removing trailing bases
        if k > self.nbases {
            ind >>= 2 * (k - self.nbases);
        }

        // Walk down levels until we find a present entry
        while lind > 0 {
            let sai_pos = self.genome_sa_index_start[(lind - 1) as usize] + ind;
            if sai_pos >= self.data.len() as u64 {
                lind -= 1;
                ind >>= 2;
                continue;
            }

            let entry = self.data.read(sai_pos as usize);
            let is_absent = (entry >> (self.gstrand_bit + 2)) & 1 != 0;

            if !is_absent {
                // Found present entry — extract SA position
                let sa_pos_mask = (1u64 << self.gstrand_bit) - 1;
                let sa_start = (entry & sa_pos_mask) as usize;

                // Get upper bound from next k-mer at same level
                let level_end = self.genome_sa_index_start[lind as usize];
                let next_pos = self.genome_sa_index_start[(lind - 1) as usize] + ind + 1;

                let (sa_end, bounds_tight) = if next_pos < level_end {
                    let next_entry = self.data.read(next_pos as usize);
                    let next_absent = (next_entry >> (self.gstrand_bit + 2)) & 1 != 0;
                    if !next_absent {
                        ((next_entry & sa_pos_mask) as usize, true)
                    } else {
                        (n_sa, false)
                    }
                } else {
                    (n_sa, false)
                };

                return Some((sa_start, sa_end, lind as usize, bounds_tight));
            }

            lind -= 1;
            ind >>= 2;
        }

        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calculate_num_indices() {
        // For nbases=1: 4^1 = 4
        assert_eq!(SaIndex::calculate_num_indices(1), 4);

        // For nbases=2: 4^1 + 4^2 = 4 + 16 = 20
        assert_eq!(SaIndex::calculate_num_indices(2), 20);

        // For nbases=3: 4 + 16 + 64 = 84
        assert_eq!(SaIndex::calculate_num_indices(3), 84);

        // For nbases=14 (STAR default): ~357 million
        let n14 = SaIndex::calculate_num_indices(14);
        assert_eq!(n14, 357_913_940);
    }
}
