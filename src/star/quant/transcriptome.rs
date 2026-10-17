// Vendored from scverse/rustar-aligner rev 207773b and trimmed to the
// `--quantTranscriptomeSAMoutput` output-mode enum: the transcriptome
// index machinery that surrounded it was unreachable in seqarc (the
// only production caller fed a `GenomeIndex::transcriptome` field that
// was never read). The enum stays because
// `src/star/params/mod.rs` holds a field of this type.

/// STAR's `--quantTranscriptomeSAMoutput` value.
///
/// Controls how alignments are filtered / adjusted before being projected
/// onto the transcriptome:
///   * `BanSingleEndBanIndelsExtendSoftclip` (STAR default, RSEM-compatible) —
///     reject alignments containing indels; reject PE alignments where both
///     sides of the read come from a single mate; extend soft-clipped bases
///     back into matched bases (re-counting mismatches).
///   * `BanSingleEnd` — only reject single-mate-only PE alignments; keep
///     indels and soft-clips as-is.
///   * `BanSingleEndExtendSoftclip` — reject single-mate-only PE alignments;
///     keep indels; extend soft-clips.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QuantTranscriptomeSAMoutput {
    /// Default: ban indels, ban single-mate PE hits, extend soft-clips.
    #[default]
    BanSingleEndBanIndelsExtendSoftclip,
    /// Ban single-mate PE hits only; keep indels and soft-clips.
    BanSingleEnd,
    /// Ban single-mate PE hits; keep indels; extend soft-clips.
    BanSingleEndExtendSoftclip,
}

impl std::str::FromStr for QuantTranscriptomeSAMoutput {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "BanSingleEnd_BanIndels_ExtendSoftclip" => {
                Ok(Self::BanSingleEndBanIndelsExtendSoftclip)
            }
            "BanSingleEnd" => Ok(Self::BanSingleEnd),
            "BanSingleEnd_ExtendSoftclip" => Ok(Self::BanSingleEndExtendSoftclip),
            _ => Err(format!(
                "unknown --quantTranscriptomeSAMoutput '{s}'; expected \
                 'BanSingleEnd_BanIndels_ExtendSoftclip', 'BanSingleEnd', or \
                 'BanSingleEnd_ExtendSoftclip'"
            )),
        }
    }
}

// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mode_from_str_all_three() {
        use std::str::FromStr;
        assert_eq!(
            QuantTranscriptomeSAMoutput::from_str("BanSingleEnd_BanIndels_ExtendSoftclip").unwrap(),
            QuantTranscriptomeSAMoutput::BanSingleEndBanIndelsExtendSoftclip,
        );
        assert_eq!(
            QuantTranscriptomeSAMoutput::from_str("BanSingleEnd").unwrap(),
            QuantTranscriptomeSAMoutput::BanSingleEnd,
        );
        assert_eq!(
            QuantTranscriptomeSAMoutput::from_str("BanSingleEnd_ExtendSoftclip").unwrap(),
            QuantTranscriptomeSAMoutput::BanSingleEndExtendSoftclip,
        );
        assert!(QuantTranscriptomeSAMoutput::from_str("garbage").is_err());
    }
}
