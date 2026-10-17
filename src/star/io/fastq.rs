// Vendored from scverse/rustar-aligner (rev 207773b, module io/fastq.rs); trimmed to the two base-encoding helpers seqarc's mapping
// path uses — the upstream FASTQ readers/writers were dropped with the CLI.

/// Encode an ASCII base to genome encoding (A=0, C=1, G=2, T=3, N=4).
pub fn encode_base(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4, // N or any ambiguous base (R, Y, S, W, K, M, etc.)
    }
}

/// Complement an encoded base (A=0↔T=3, C=1↔G=2, N=4→N=4).
pub fn complement_base(encoded: u8) -> u8 {
    match encoded {
        0 => 3,       // A -> T
        1 => 2,       // C -> G
        2 => 1,       // G -> C
        3 => 0,       // T -> A
        _ => encoded, // N -> N
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_base() {
        assert_eq!(encode_base(b'A'), 0);
        assert_eq!(encode_base(b'a'), 0);
        assert_eq!(encode_base(b'C'), 1);
        assert_eq!(encode_base(b'G'), 2);
        assert_eq!(encode_base(b'T'), 3);
        assert_eq!(encode_base(b'N'), 4);
        assert_eq!(encode_base(b'R'), 4);
    }

    #[test]
    fn test_complement_base() {
        assert_eq!(complement_base(0), 3);
        assert_eq!(complement_base(1), 2);
        assert_eq!(complement_base(2), 1);
        assert_eq!(complement_base(3), 0);
        assert_eq!(complement_base(4), 4);
    }
}
