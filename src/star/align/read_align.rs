/// Read alignment driver function
use crate::star::align::score::{AlignmentScorer, SpliceMotif};
use crate::star::align::seed::Seed;
use crate::star::align::stitch::{cluster_seeds, stitch_seeds_with_jdb_debug};
use crate::star::align::transcript::Transcript;
use crate::star::error::Error;
use crate::star::index::GenomeIndex;
use crate::star::params::{IntronMotifFilter, IntronStrandFilter, MultimapperOrder, Parameters};
use crate::star::stats::UnmappedReason;
use std::hash::{DefaultHasher, Hash, Hasher};

/// Derive a deterministic per-read RNG seed from `run_rng_seed` + the read name.
///
/// STAR seeds `std::mt19937` once per chunk/thread (`runRNGseed*(iChunk+1)`),
/// then advances the state sequentially per read. rustar-aligner parallelises per-read
/// via rayon, so we instead fold the read name into the seed — this keeps tie
/// breaks reproducible regardless of thread count while still honoring the
/// user's `--runRNGseed` value.
pub(crate) fn per_read_seed(run_rng_seed: u64, read_name: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    read_name.hash(&mut hasher);
    run_rng_seed.wrapping_mul(hasher.finish().wrapping_add(1))
}

/// Shuffle the prefix of `items` whose `score_fn` equals the first element's score.
///
/// Mirrors STAR's `ReadAlign_multMapSelect` / `funPrimaryAlignMark`: best-scoring
/// alignments are randomized so primary selection (index 0) is not biased by
/// upstream sort order. Non-tied elements are left alone.
fn shuffle_tied_prefix<T>(items: &mut [T], score_fn: impl Fn(&T) -> i32, seed: u64) {
    let Some(first) = items.first() else {
        return;
    };
    let best = score_fn(first);
    let tied = items.iter().take_while(|t| score_fn(t) == best).count();
    if tied < 2 {
        return;
    }
    crate::star::rng::shuffle_deterministic(&mut items[..tied], seed);
}

/// Result of aligning a single read: (transcripts, chimeric_alignments, n_for_mapq, unmapped_reason)
pub type AlignReadResult = (
    Vec<Transcript>,
    Vec<crate::star::chimeric::ChimericAlignment>,
    usize,
    Option<UnmappedReason>,
);


/// Align a read to the genome.
///
/// # Algorithm
/// 1. Find seeds (exact matches) using MMP search
/// 2. Cluster seeds by genomic proximity
/// 3. Stitch seeds within each cluster using DP
/// 4. Filter transcripts by quality thresholds
/// 5. Sort by score and limit to top N
/// 6. Detect chimeric alignments if enabled
///
/// # Arguments
/// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
/// * `read_name` - Read name (needed for chimeric output)
/// * `index` - Genome index
/// * `params` - User parameters
///
/// # Returns
/// Tuple of (transcripts, chimeric alignments, n_for_mapq, unmapped_reason):
/// - transcripts: sorted by score (best first)
/// - chimeric alignments: sorted by score (best first)
/// - n_for_mapq: effective alignment count for MAPQ calculation (max of transcript count
///   and valid cluster count, to avoid undercounting from coordinate dedup on tandem repeats)
/// - unmapped_reason: `Some(reason)` if no alignments produced, `None` if mapped
pub fn align_read(
    read_seq: &[u8],
    read_name: &str,
    index: &GenomeIndex,
    params: &Parameters,
) -> Result<AlignReadResult, Error> {
    let debug_read = !params.read_name_filter.is_empty() && read_name == params.read_name_filter;

    // Step 1: Find seeds (seedMapMin from params)
    let min_seed_length = params.seed_map_min;
    let seeds = Seed::find_seeds(
        read_seq,
        index,
        min_seed_length,
        params,
        if debug_read { read_name } else { "" },
    )?;

    if debug_read {
        let total_positions: usize = seeds.iter().map(|s| s.sa_end - s.sa_start).sum();
        eprintln!(
            "[DEBUG {}] Seeds: {} seeds, {} total SA positions, read_len={}",
            read_name,
            seeds.len(),
            total_positions,
            read_seq.len()
        );
        let n_lr = seeds.iter().filter(|s| !s.search_rc).count();
        let n_rl = seeds.iter().filter(|s| s.search_rc).count();
        eprintln!(
            "  {} seeds total: {} L→R (sparse), {} R→L (sparse)",
            seeds.len(),
            n_lr,
            n_rl
        );
        for (i, s) in seeds.iter().enumerate() {
            let n_loci = s.sa_end - s.sa_start;
            eprintln!(
                "  seed[{}]: read_pos={}, len={}, n_loci={}, search_rc={}, sa=[{},{})",
                i, s.read_pos, s.length, n_loci, s.search_rc, s.sa_start, s.sa_end
            );
        }
    }

    if seeds.is_empty() {
        if debug_read {
            eprintln!("[DEBUG {read_name}] No seeds found — unmapped");
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Step 2: Cluster seeds (STAR's bin-based windowing)
    // seed_per_window_nmax capacity eviction is handled inside cluster_seeds()
    let clusters = cluster_seeds(&seeds, index, params, read_seq.len(), debug_read);

    if debug_read {
        eprintln!(
            "[DEBUG {}] Clusters: {} clusters",
            read_name,
            clusters.len()
        );
        for (i, cluster) in clusters.iter().enumerate() {
            let chr_name = if cluster.chr_idx < index.genome.chr_name.len() {
                &index.genome.chr_name[cluster.chr_idx]
            } else {
                "unknown"
            };
            eprintln!(
                "  cluster[{}]: chr={}, is_reverse={}, seeds={}, anchor_bin={}",
                i,
                chr_name,
                cluster.is_reverse,
                cluster.alignments.len(),
                cluster.anchor_bin,
            );
            for (j, wa) in cluster.alignments.iter().enumerate() {
                let chr_pos = wa.genome_pos.saturating_sub(
                    if cluster.chr_idx < index.genome.chr_start.len() {
                        index.genome.chr_start[cluster.chr_idx]
                    } else {
                        0
                    },
                ) + 1; // 1-based
                eprintln!(
                    "    wa[{}]: read_pos={}, len={}, genome_pos={} ({}:{}), n_rep={}, is_anchor={}",
                    j,
                    wa.read_pos,
                    wa.length,
                    wa.genome_pos,
                    chr_name,
                    chr_pos,
                    wa.n_rep,
                    wa.is_anchor
                );
                if j >= 5 {
                    eprintln!(
                        "    ... ({} more WA entries)",
                        cluster.alignments.len() - j - 1
                    );
                    break;
                }
            }
        }
    }

    if clusters.is_empty() {
        if debug_read {
            eprintln!("[DEBUG {read_name}] No clusters — unmapped");
        }
        return Ok((Vec::new(), Vec::new(), 0, Some(UnmappedReason::Other)));
    }

    // Cap total clusters (alignWindowsPerReadNmax)
    let mut clusters = clusters;
    clusters.truncate(params.align_windows_per_read_nmax);

    // NOTE: STAR's winReadCoverageRelativeMin filter is long-reads-only
    // (#ifdef COMPILE_FOR_LONG_READS in stitchPieces.cpp). Standard STAR
    // does NOT filter clusters by seed coverage. Removed to match STAR.

    // Step 2b: Detect chimeric alignments from multi-cluster seeds (Tier 2)
    let mut chimeric_alignments = Vec::new();
    if params.chim_segment_min > 0 && clusters.len() > 1 {
        use crate::star::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        chimeric_alignments
            .extend(detector.detect_from_multi_clusters(&clusters, read_seq, read_name, index)?);
    }

    // Step 3: Stitch seeds within each cluster
    let scorer = AlignmentScorer::from_params(params);
    let mut transcripts = Vec::new();
    // Collect all raw (pre-dedup) transcripts for chimericDetectionOld (Tier 1).
    let mut all_raw_transcripts: Vec<crate::star::align::transcript::Transcript> = Vec::new();

    // Use junction DB for annotation-aware scoring if available
    let junction_db = if index.junction_db.is_empty() {
        None
    } else {
        Some(&index.junction_db)
    };

    for (ci, cluster) in clusters.iter().enumerate() {
        let debug_name = if debug_read { read_name } else { "" };
        let cluster_transcripts = stitch_seeds_with_jdb_debug(
            cluster,
            read_seq,
            index,
            &scorer,
            junction_db,
            params.align_transcripts_per_window_nmax,
            debug_name,
        );
        if debug_read {
            eprintln!(
                "[DEBUG {}] Cluster[{}]: {} transcripts from DP",
                read_name,
                ci,
                cluster_transcripts.len()
            );
            for (ti, t) in cluster_transcripts.iter().enumerate().take(5) {
                let chr_name = if t.chr_idx < index.genome.chr_name.len() {
                    &index.genome.chr_name[t.chr_idx]
                } else {
                    "unknown"
                };
                eprintln!(
                    "  transcript[{}]: chr={}:{}-{} ({}) score={} mm={} junctions={} cigar={}",
                    ti,
                    chr_name,
                    t.genome_start,
                    t.genome_end,
                    if t.is_reverse { "-" } else { "+" },
                    t.score,
                    t.n_mismatch,
                    t.n_junction,
                    t.cigar_string()
                );
            }
        }
        if params.chim_segment_min > 0 {
            all_raw_transcripts.extend(cluster_transcripts.iter().cloned());
        }
        transcripts.extend(cluster_transcripts);
    }

    // Step 4a: Deduplicate and score-range filter — BEFORE quality filters.
    // STAR order: multMapSelect (score-range) → mappedFilter (quality gates).
    // Doing quality filters first is wrong: it can remove the high-scoring primary,
    // leaving a lower-scoring secondary that then passes as the "best" alignment.

    // Deduplicate transcripts with identical genomic coordinates AND CIGAR,
    // PRESERVING discovery (window) order. STAR's primary tie-break falls back
    // to the earliest-discovered alignment, so we must not reorder here (the
    // old code sorted by coordinate, destroying that order).
    {
        let mut seen = std::collections::HashSet::new();
        transcripts.retain(|t| {
            seen.insert((
                t.chr_idx,
                t.genome_start,
                t.genome_end,
                t.is_reverse,
                t.cigar_string(),
            ))
        });
    }

    // Deterministic primary tie-break (score, then a fixed positional order).
    transcripts.sort_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then_with(|| a.n_junction.cmp(&b.n_junction))
            .then_with(|| a.chr_idx.cmp(&b.chr_idx))
            .then_with(|| a.genome_start.cmp(&b.genome_start))
            .then_with(|| a.is_reverse.cmp(&b.is_reverse))
    });

    // Primary selection — STAR's multMapSelect (ReadAlign_multMapSelect.cpp).
    //
    // STAR's DEFAULT (`--outMultimapperOrder Old_2.4`) does NOT consult the RNG
    // for primary selection; it marks the deterministic best alignment primary.
    // Only `--outMultimapperOrder Random` shuffles. Previously rustar-aligner
    // shuffled unconditionally, which randomised the primary among equal-score
    // loci and diverged from STAR's deterministic choice. Gate the shuffle on
    // the Random mode so the default is deterministic and STAR-faithful.
    //
    // Under Random, we shuffle the tied top-score prefix with a per-read seed
    // (deterministic per read → thread-count invariant).
    if params.out_multimapper_order == MultimapperOrder::Random {
        shuffle_tied_prefix(
            &mut transcripts,
            |t| t.score,
            per_read_seed(params.run_rng_seed, read_name),
        );
    }

    // Score-range filter: keep only alignments within outFilterMultimapScoreRange of the best.
    // (STAR's multMapSelect step — must run before quality filters.)
    if !transcripts.is_empty() {
        let max_score = transcripts[0].score;
        let score_threshold = max_score - params.out_filter_multimap_score_range;
        transcripts.retain(|t| t.score >= score_threshold);
    }

    // Step 4: Quality filters (STAR's mappedFilter — runs after score-range selection).
    // STAR uses (Lread-1) for relative thresholds and casts to integer
    // (ReadAlign_mappedFilter.cpp lines 8-9)
    let read_length = read_seq.len() as f64;
    let lread_m1 = (read_seq.len() as f64) - 1.0;

    // Log filtering statistics
    let pre_filter_count = transcripts.len();
    let mut filter_reasons = std::collections::HashMap::new();

    transcripts.retain(|t| {
        // Absolute score threshold
        if t.score < params.out_filter_score_min {
            *filter_reasons.entry("score_min").or_insert(0) += 1;
            return false;
        }

        // Relative score threshold: STAR casts to intScore (i32)
        if t.score < (params.out_filter_score_min_over_lread * lread_m1) as i32 {
            *filter_reasons.entry("score_min_relative").or_insert(0) += 1;
            return false;
        }

        // Absolute mismatch count
        if t.n_mismatch > params.out_filter_mismatch_nmax {
            *filter_reasons.entry("mismatch_max").or_insert(0) += 1;
            log::debug!(
                "Filtered {}: {} mismatches > {} max (read_len={}, score={})",
                read_name,
                t.n_mismatch,
                params.out_filter_mismatch_nmax,
                read_length,
                t.score
            );
            return false;
        }

        // Relative mismatch count (mismatches / read_length)
        let mismatch_rate = t.n_mismatch as f64 / read_length;
        if mismatch_rate > params.out_filter_mismatch_nover_lmax {
            *filter_reasons.entry("mismatch_rate").or_insert(0) += 1;
            log::debug!(
                "Filtered {}: {:.1}% mismatch rate > {:.1}% max ({}/{} bases, score={})",
                read_name,
                mismatch_rate * 100.0,
                params.out_filter_mismatch_nover_lmax * 100.0,
                t.n_mismatch,
                read_length,
                t.score
            );
            return false;
        }

        // Absolute matched bases
        let n_matched = t.n_matched();
        if n_matched < params.out_filter_match_nmin as usize {
            *filter_reasons.entry("match_min").or_insert(0) += 1;
            return false;
        }

        // Relative matched bases: STAR casts to uint (u32)
        if (n_matched as f64) < params.out_filter_match_nmin_over_lread * lread_m1 {
            *filter_reasons.entry("match_min_relative").or_insert(0) += 1;
            return false;
        }

        // Junction motif filtering
        match params.out_filter_intron_motifs {
            IntronMotifFilter::None => {
                // Accept all motifs
            }
            IntronMotifFilter::RemoveNoncanonical => {
                // Reject if any junction is non-canonical
                if t.junction_motifs.contains(&SpliceMotif::NonCanonical) {
                    *filter_reasons.entry("noncanonical_junction").or_insert(0) += 1;
                    return false;
                }
            }
            IntronMotifFilter::RemoveNoncanonicalUnannotated => {
                // Only reject if a non-canonical junction is NOT annotated in GTF
                if t.junction_motifs
                    .iter()
                    .zip(t.junction_annotated.iter())
                    .any(|(m, annotated)| *m == SpliceMotif::NonCanonical && !annotated)
                {
                    *filter_reasons
                        .entry("noncanonical_unannotated_junction")
                        .or_insert(0) += 1;
                    return false;
                }
            }
        }

        // Intron strand consistency filtering (outFilterIntronStrands)
        // STAR's RemoveInconsistentStrands removes transcripts that have junctions
        // with mixed intron strand (some imply + strand, some imply - strand).
        // This handles chimeric/impossible transcripts spanning both strands.
        // Note: a reverse-strand read CAN have + strand motifs (antisense reads
        // from + strand genes in unstranded RNA-seq) — this is valid and STAR
        // keeps such reads. Only mixed-strand within one transcript is filtered.
        if params.out_filter_intron_strands == IntronStrandFilter::RemoveInconsistentStrands {
            let mut has_plus = false;
            let mut has_minus = false;
            for motif in &t.junction_motifs {
                match motif.implied_strand() {
                    Some('+') => has_plus = true,
                    Some('-') => has_minus = true,
                    _ => {}
                }
            }
            if has_plus && has_minus {
                *filter_reasons.entry("inconsistent_strand").or_insert(0) += 1;
                return false;
            }
        }

        true
    });

    // Log filtering summary if anything was filtered
    if pre_filter_count > transcripts.len() {
        let filtered = pre_filter_count - transcripts.len();
        log::debug!(
            "Read {read_name}: Filtered {filtered}/{pre_filter_count} transcripts: {filter_reasons:?}"
        );
    }

    if debug_read {
        eprintln!(
            "[DEBUG {}] After quality filters: {}/{} transcripts remain (reasons: {:?})",
            read_name,
            transcripts.len(),
            pre_filter_count,
            filter_reasons
        );
    }

    // Step 3b: chimericDetectionOld (Tier 1) — STAR-faithful post-stitching transcript-pair search.
    // Uses the best post-dedup transcript as the primary segment and searches all raw transcripts
    // (pre-dedup, from all clusters) for the best complementary partner.
    // (Vendored port: rewritten let-chains — edition 2021.)
    if params.chim_segment_min > 0
        && !all_raw_transcripts.is_empty()
    {
        if let Some(tr_best) = transcripts.first() {
            use crate::star::chimeric::detect_chimeric_old;
            let chims = detect_chimeric_old(
                &all_raw_transcripts,
                tr_best,
                read_seq,
                read_name,
                params,
                index,
            )?;
            chimeric_alignments.extend(chims);
        }
    }

    // Step 3c: Soft-clip re-mapping (Phase 12.2) — try to align soft-clipped bases when
    // detect_chimeric_old found no chimeric partner in the existing transcript pool.
    if params.chim_segment_min > 0
        && chimeric_alignments.is_empty()
    {
        if let Some(tr_best) = transcripts.first() {
            use crate::star::chimeric::ChimericDetector;
            let detector = ChimericDetector::new(params);
            if let Some(chim) = detector.detect_from_soft_clips(tr_best, read_seq, read_name, index)? {
                chimeric_alignments.push(chim);
            }
        }
    }

    // Step 3d: Tier 3 — re-seed outer uncovered read regions of each chimeric pair (Phase 17.10).
    // Extends 2-segment chimeras toward multi-junction fusions by seeding the read bases
    // that lie outside both chimeric segments.
    if params.chim_segment_min > 0 && !chimeric_alignments.is_empty() {
        use crate::star::chimeric::ChimericDetector;
        let detector = ChimericDetector::new(params);
        let mut tier3 = Vec::new();
        for chim in &chimeric_alignments {
            let extras =
                detector.detect_from_chimeric_residuals(chim, read_seq, read_name, index)?;
            tier3.extend(extras);
        }
        chimeric_alignments.extend(tier3);
    }

    // Note: STAR sometimes finds 2 equivalent indel placements in homopolymer runs
    // via its recursive stitcher's seed exploration (NH=2 instead of NH=1 for ~5 reads).
    // Generating equivalents post-hoc causes more harm than good (41 false NH=2 vs 5 fixed).
    // The root cause is jR scanning placing insertions at different positions — fixing that
    // would be a better approach than post-hoc enumeration.

    // Step 6: Filter chimeric alignments
    if params.chim_segment_min > 0 {
        chimeric_alignments.retain(|chim| {
            chim.meets_min_segment_length(params.chim_segment_min)
                && chim.meets_min_score(params.chim_score_min)
        });
    }

    // n_for_mapq = transcripts.len() after dedup and filtering.
    // Multi-transcript DP (Phase 16.10) produces multiple transcripts per window
    // for tandem repeats (e.g. rDNA), yielding correct NH → correct MAPQ.
    let n_for_mapq = transcripts.len();

    if debug_read {
        eprintln!(
            "[DEBUG {}] Final: {} transcripts, n_for_mapq={}",
            read_name,
            transcripts.len(),
            n_for_mapq
        );
        for (i, t) in transcripts.iter().enumerate() {
            let chr_name = if t.chr_idx < index.genome.chr_name.len() {
                &index.genome.chr_name[t.chr_idx]
            } else {
                "unknown"
            };
            eprintln!(
                "  FINAL[{}]: chr={}:{}-{} ({}) score={} mm={} junctions={} cigar={}",
                i,
                chr_name,
                t.genome_start,
                t.genome_end,
                if t.is_reverse { "-" } else { "+" },
                t.score,
                t.n_mismatch,
                t.n_junction,
                t.cigar_string()
            );
        }
    }

    // STAR's mappedFilter (ReadAlign_mappedFilter.cpp) classifies the read by its
    // BEST transcript, in strict order: too-short (score/match) → too-many-mismatches
    // → too-many-loci → mapped. Crucially, `multi` (too many loci) is reached ONLY
    // when the best alignment PASSES the score/match/mismatch gates. We mirror that:
    // if the quality filter emptied the set, the best transcript failed → classify
    // short/mismatch; otherwise, if too many loci survive, it's multi; else mapped.
    let unmapped_reason = if transcripts.is_empty() {
        let has_mismatch = filter_reasons.contains_key("mismatch_max")
            || filter_reasons.contains_key("mismatch_rate");
        let has_other = filter_reasons
            .keys()
            .any(|k| *k != "mismatch_max" && *k != "mismatch_rate");
        if has_mismatch && !has_other {
            Some(UnmappedReason::TooManyMismatches)
        } else {
            Some(UnmappedReason::TooShort)
        }
    } else if transcripts.len() > params.out_filter_multimap_nmax as usize {
        // Best alignment passed quality, but too many loci → STAR unmapType=3.
        // `n_for_mapq` already holds this count (computed above); drop the alignments.
        transcripts.clear();
        Some(UnmappedReason::TooManyLoci)
    } else {
        None
    };

    Ok((
        transcripts,
        chimeric_alignments,
        n_for_mapq,
        unmapped_reason,
    ))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::star::align::transcript::Exon;
    use crate::star::genome::Genome;
    use crate::star::index::packed_array::PackedArray;
    use crate::star::index::sa_index::SaIndex;
    use crate::star::index::suffix_array::SuffixArray;
    use noodles::sam::alignment::record::cigar;

    fn default_params() -> Parameters {
        // Parse empty args to get default parameters
        Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"])
    }

    fn make_test_index() -> GenomeIndex {
        // Simple genome: ACGTACGTNN (10 bases)
        let seq = vec![0, 1, 2, 3, 0, 1, 2, 3, 4, 4];
        let n_genome = 64u64; // Padded
        let mut sequence = vec![5u8; (n_genome * 2) as usize];
        sequence[0..seq.len()].copy_from_slice(&seq);

        // Build reverse complement
        for i in 0..n_genome as usize {
            let base = sequence[i];
            let complement = if base < 4 { 3 - base } else { base };
            sequence[2 * n_genome as usize - 1 - i] = complement;
        }

        let genome = Genome {

            sequence: sequence.into(),
            n_genome,
            n_genome_real: n_genome,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![10],
            chr_start: vec![0, n_genome],
        };

        // Create dummy SA and SAindex (would need real index for actual alignment)
        let gstrand_bit = 33;
        let suffix_array = SuffixArray {
            data: PackedArray::new(gstrand_bit, 0),
            gstrand_bit,
            gstrand_mask: (1u64 << gstrand_bit) - 1,
        };

        let word_length = gstrand_bit + 3;
        let sa_index = SaIndex {
            data: PackedArray::new(word_length, 0),
            nbases: 14,
            genome_sa_index_start: vec![0],
            gstrand_bit,
        };

        GenomeIndex {
            genome,
            suffix_array,
            sa_index,
            junction_db: crate::star::junction::SpliceJunctionDb::empty(),
            prepared_junctions: Vec::new(),
            sjdb_overhang: 0,
        }
    }

    #[test]
    fn test_align_read_no_seeds() {
        let index = make_test_index();
        let params = default_params();

        // Read with all N's (no seeds possible)
        let read_seq = vec![4, 4, 4, 4, 4, 4, 4, 4, 4, 4];

        let result = align_read(&read_seq, "READ_001", &index, &params);
        assert!(result.is_ok());

        let (transcripts, chimeras, n_for_mapq, unmapped_reason) = result.unwrap();
        assert_eq!(transcripts.len(), 0); // No alignment
        assert_eq!(chimeras.len(), 0); // No chimeric alignments
        assert_eq!(n_for_mapq, 0);
        assert_eq!(unmapped_reason, Some(UnmappedReason::Other));
    }

    #[test]
    fn test_transcript_filtering_score() {
        let index = make_test_index();
        let mut params = default_params();
        params.out_filter_score_min = 50;

        // Would need actual seeds and alignment to test this properly
        // This test just verifies the function doesn't crash
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_002", &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_transcript_filtering_mismatch() {
        let index = make_test_index();
        let mut params = default_params();
        params.out_filter_mismatch_nmax = 2;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_003", &index, &params);
        assert!(result.is_ok());
    }

    #[test]
    fn test_unmapped_reason_mismatch_classification() {
        // Verify the filter_reasons logic used to derive unmapped_reason.
        // mismatch-only → TooManyMismatches; anything else (or mixed) → TooShort.
        let classify = |reasons: &[&str]| -> UnmappedReason {
            let filter_reasons: std::collections::HashMap<&str, i32> =
                reasons.iter().map(|k| (*k, 1)).collect();
            let has_mismatch = filter_reasons.contains_key("mismatch_max")
                || filter_reasons.contains_key("mismatch_rate");
            let has_other = filter_reasons
                .keys()
                .any(|k| *k != "mismatch_max" && *k != "mismatch_rate");
            if has_mismatch && !has_other {
                UnmappedReason::TooManyMismatches
            } else {
                UnmappedReason::TooShort
            }
        };

        assert_eq!(
            classify(&["mismatch_max"]),
            UnmappedReason::TooManyMismatches
        );
        assert_eq!(
            classify(&["mismatch_rate"]),
            UnmappedReason::TooManyMismatches
        );
        assert_eq!(
            classify(&["mismatch_max", "mismatch_rate"]),
            UnmappedReason::TooManyMismatches
        );
        // Mixed with score filter → TooShort
        assert_eq!(
            classify(&["mismatch_max", "score_min"]),
            UnmappedReason::TooShort
        );
        // Score-only → TooShort
        assert_eq!(classify(&["score_min"]), UnmappedReason::TooShort);
        assert_eq!(classify(&["match_min"]), UnmappedReason::TooShort);
    }

    #[test]
    fn test_transcript_multimap_limit() {
        let index = make_test_index();
        let mut params = default_params();
        params.out_filter_multimap_nmax = 5;

        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let result = align_read(&read_seq, "READ_004", &index, &params);
        assert!(result.is_ok());

        let (transcripts, _chimeras, _n_for_mapq, _reason) = result.unwrap();
        assert!(transcripts.len() <= 5);
    }

    #[test]
    fn test_strand_consistency_filter() {
        use crate::star::align::transcript::Transcript;
        use crate::star::params::IntronStrandFilter;
        use cigar::op::{Kind, Op};

        // Create a transcript with conflicting strand motifs (mixed + and - within one transcript)
        let t_inconsistent = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::CtAc], // +strand and -strand
            junction_annotated: vec![],
        };

        // Create a transcript with consistent strand motifs (all + strand)
        let t_consistent = Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1300,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1300,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_junction: 2,
            junction_motifs: vec![SpliceMotif::GtAg, SpliceMotif::GcAg], // both + strand
            junction_annotated: vec![],
        };

        // Note: STAR's RemoveInconsistentStrands filters transcripts where
        // junctions have MIXED implied strand (both + and - within one transcript).
        // It does NOT compare junction strand vs alignment strand — a reverse-strand
        // read at a GT/AG junction (antisense of + strand gene) is valid and kept.

        // Verify mixed-strand is detected
        let mut has_plus = false;
        let mut has_minus = false;
        for motif in &t_inconsistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && has_minus); // Inconsistent (mixed strands)

        // Verify consistent transcript has no conflict
        has_plus = false;
        has_minus = false;
        for motif in &t_consistent.junction_motifs {
            match motif.implied_strand() {
                Some('+') => has_plus = true,
                Some('-') => has_minus = true,
                _ => {}
            }
        }
        assert!(has_plus && !has_minus); // Consistent (all +)

        // Also verify: single CT/AC on forward-strand or single GT/AG on reverse-strand
        // are NOT filtered (STAR keeps these — antisense reads are valid)
        let ctac_only = vec![SpliceMotif::CtAc];
        let (mut hp, mut hm) = (false, false);
        for m in &ctac_only {
            match m.implied_strand() {
                Some('+') => hp = true,
                Some('-') => hm = true,
                _ => {}
            }
        }
        assert!(!hp && hm); // Only minus → NOT mixed → NOT filtered

        // Verify the filter enum
        assert_ne!(
            IntronStrandFilter::None,
            IntronStrandFilter::RemoveInconsistentStrands
        );
    }

    #[test]
    fn test_noncanonical_unannotated_filter() {
        use crate::star::align::score::SpliceMotif;
        use crate::star::align::transcript::Transcript;
        use cigar::op::{Kind, Op};

        // Helper: check if a transcript would be filtered by RemoveNoncanonicalUnannotated
        // (mirrors the logic in the retain closure)
        let would_filter = |t: &Transcript| -> bool {
            t.junction_motifs
                .iter()
                .zip(t.junction_annotated.iter())
                .any(|(m, annotated)| *m == SpliceMotif::NonCanonical && !annotated)
        };

        let base_transcript = || Transcript {
            chr_idx: 0,
            genome_start: 1000,
            genome_end: 1200,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 1000,
                genome_end: 1200,
                read_start: 0,
                read_end: 100,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 100)],
            score: 100,
            n_mismatch: 0,
            n_junction: 1,
            junction_motifs: vec![],
            junction_annotated: vec![],
        };

        // Case 1: NonCanonical + unannotated → should be filtered
        let mut t1 = base_transcript();
        t1.junction_motifs = vec![SpliceMotif::NonCanonical];
        t1.junction_annotated = vec![false];
        assert!(
            would_filter(&t1),
            "NonCanonical + unannotated should be filtered"
        );

        // Case 2: NonCanonical + annotated → should be KEPT
        let mut t2 = base_transcript();
        t2.junction_motifs = vec![SpliceMotif::NonCanonical];
        t2.junction_annotated = vec![true];
        assert!(
            !would_filter(&t2),
            "NonCanonical + annotated should be kept"
        );

        // Case 3: Canonical + unannotated → should be KEPT
        let mut t3 = base_transcript();
        t3.junction_motifs = vec![SpliceMotif::GtAg];
        t3.junction_annotated = vec![false];
        assert!(!would_filter(&t3), "Canonical + unannotated should be kept");

        // Case 4: Mixed — one canonical + one non-canonical unannotated → filtered
        let mut t4 = base_transcript();
        t4.junction_motifs = vec![SpliceMotif::GtAg, SpliceMotif::NonCanonical];
        t4.junction_annotated = vec![true, false];
        assert!(
            would_filter(&t4),
            "Mixed with unannotated non-canonical should be filtered"
        );

        // Case 5: Mixed — one canonical + one non-canonical annotated → kept
        let mut t5 = base_transcript();
        t5.junction_motifs = vec![SpliceMotif::GtAg, SpliceMotif::NonCanonical];
        t5.junction_annotated = vec![false, true];
        assert!(
            !would_filter(&t5),
            "Mixed with annotated non-canonical should be kept"
        );
    }

    #[test]
    fn shuffle_tied_prefix_is_deterministic() {
        // Same seed + same input → same permutation on reruns.
        let items: Vec<(i32, u32)> = (0..8).map(|i| (100, i)).collect();
        let mut a = items.clone();
        let mut b = items.clone();
        shuffle_tied_prefix(&mut a, |t| t.0, 12345);
        shuffle_tied_prefix(&mut b, |t| t.0, 12345);
        assert_eq!(a, b);
    }

    #[test]
    fn shuffle_tied_prefix_respects_ties() {
        // Only the top-score prefix gets shuffled; lower-scored tail is left alone.
        let mut items = vec![(100, 0u32), (100, 1), (100, 2), (50, 3), (40, 4)];
        shuffle_tied_prefix(&mut items, |t| t.0, 777);
        // Last two elements (non-tied) stay in place.
        assert_eq!(items[3], (50, 3));
        assert_eq!(items[4], (40, 4));
        // Tied prefix contains the original three items in some order.
        let mut top: Vec<u32> = items[..3].iter().map(|t| t.1).collect();
        top.sort_unstable();
        assert_eq!(top, vec![0, 1, 2]);
    }

    #[test]
    fn shuffle_tied_prefix_different_seeds_can_diverge() {
        // Probabilistic: for a tied set of 8, at least two seeds should disagree
        // on the chosen primary. (Exhaustive over a small seed range is fine.)
        let base: Vec<(i32, u32)> = (0..8).map(|i| (100, i)).collect();
        let mut firsts = std::collections::HashSet::new();
        for seed in 0..32u64 {
            let mut v = base.clone();
            shuffle_tied_prefix(&mut v, |t| t.0, seed);
            firsts.insert(v[0].1);
        }
        assert!(
            firsts.len() >= 2,
            "expected different seeds to pick different primaries, got {firsts:?}"
        );
    }

    #[test]
    fn shuffle_tied_prefix_noop_when_no_ties() {
        let mut items = vec![(100, 0u32), (90, 1), (80, 2)];
        let before = items.clone();
        shuffle_tied_prefix(&mut items, |t| t.0, 42);
        assert_eq!(items, before);
    }
}
