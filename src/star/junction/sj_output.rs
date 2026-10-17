// Provenance: rest of this module (`SpliceJunctionStats` SJ.out.tab
// accumulation/filtering/writing) was dead in seqarc and removed; only
// `encode_motif` remains, still used by the live sjdb pipeline.

use crate::star::align::score::SpliceMotif;

/// Encode splice motif for SJ.out.tab
///
/// Direct mapping from detected motif to STAR's numeric encoding.
/// The motif is already detected on the forward genome, so no strand
/// transformation is needed — strand is derived separately from the motif.
///
/// STAR convention:
/// 0 = non-canonical
/// 1 = GT/AG
/// 2 = CT/AC
/// 3 = GC/AG
/// 4 = CT/GC
/// 5 = AT/AC
/// 6 = GT/AT
pub(crate) fn encode_motif(motif: SpliceMotif) -> u8 {
    match motif {
        SpliceMotif::GtAg => 1,
        SpliceMotif::CtAc => 2,
        SpliceMotif::GcAg => 3,
        SpliceMotif::CtGc => 4,
        SpliceMotif::AtAc => 5,
        SpliceMotif::GtAt => 6,
        SpliceMotif::NonCanonical => 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_motif_gtag() {
        assert_eq!(encode_motif(SpliceMotif::GtAg), 1);
    }

    #[test]
    fn test_encode_motif_ctac() {
        assert_eq!(encode_motif(SpliceMotif::CtAc), 2);
    }

    #[test]
    fn test_encode_motif_gcag() {
        assert_eq!(encode_motif(SpliceMotif::GcAg), 3);
    }

    #[test]
    fn test_encode_motif_ctgc() {
        assert_eq!(encode_motif(SpliceMotif::CtGc), 4);
    }

    #[test]
    fn test_encode_motif_atac() {
        assert_eq!(encode_motif(SpliceMotif::AtAc), 5);
    }

    #[test]
    fn test_encode_motif_gtat() {
        assert_eq!(encode_motif(SpliceMotif::GtAt), 6);
    }

    #[test]
    fn test_encode_motif_noncanonical() {
        assert_eq!(encode_motif(SpliceMotif::NonCanonical), 0);
    }
}
