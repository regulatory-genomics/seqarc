// Vendored from scverse/rustar-aligner (rev 207773b); trimmed to the one item
// seqarc's mapping path uses — the upstream run-stats machinery was dropped
// with the CLI.

/// Reason a read could not be mapped
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnmappedReason {
    /// No seeds or no clusters found
    Other,
    /// Transcripts generated but all filtered out (score, match count, length)
    TooShort,
    /// All transcripts filtered by mismatch count/rate
    TooManyMismatches,
    /// Too many multi-mapping loci (nTr > outFilterMultimapNmax)
    TooManyLoci,
}
