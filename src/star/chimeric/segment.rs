//! ChimericSegment and ChimericAlignment data structures

/// A single segment of a chimeric alignment
#[derive(Debug, Clone)]
pub struct ChimericSegment {
    pub chr_idx: usize,
    pub genome_start: u64,
    pub genome_end: u64,
    pub is_reverse: bool,
    pub read_start: usize,
    pub read_end: usize,
    pub score: i32,
}

impl ChimericSegment {
    /// Get segment length in read coordinates
    pub fn read_length(&self) -> usize {
        self.read_end - self.read_start
    }

    /// Check if segment meets minimum length requirement
    pub fn meets_min_length(&self, min_len: u32) -> bool {
        self.read_length() >= min_len as usize
    }

}

/// A chimeric alignment consisting of two segments
#[derive(Debug, Clone)]
pub struct ChimericAlignment {
    pub donor: ChimericSegment,
    pub acceptor: ChimericSegment,
    pub total_score: i32,
}

impl ChimericAlignment {
    /// Create a new chimeric alignment
    pub fn new(donor: ChimericSegment, acceptor: ChimericSegment) -> Self {
        let total_score = donor.score + acceptor.score;
        Self {
            donor,
            acceptor,
            total_score,
        }
    }

    /// Check if both segments meet minimum length requirement
    pub fn meets_min_segment_length(&self, min_len: u32) -> bool {
        self.donor.meets_min_length(min_len) && self.acceptor.meets_min_length(min_len)
    }

    /// Check if total score meets minimum threshold
    pub fn meets_min_score(&self, min_score: i32) -> bool {
        self.total_score >= min_score
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    fn mock_segment(
        chr_idx: usize,
        genome_start: u64,
        genome_end: u64,
        read_start: usize,
        read_end: usize,
    ) -> ChimericSegment {
        ChimericSegment {
            chr_idx,
            genome_start,
            genome_end,
            is_reverse: false,
            read_start,
            read_end,
            score: 100,
        }
    }

    #[test]
    fn test_segment_lengths() {
        let seg = mock_segment(0, 1000, 1050, 0, 50);
        assert_eq!(seg.read_length(), 50);
    }

    #[test]
    fn test_segment_min_length() {
        let seg = mock_segment(0, 1000, 1050, 0, 50);
        assert!(seg.meets_min_length(20));
        assert!(seg.meets_min_length(50));
        assert!(!seg.meets_min_length(51));
    }

    #[test]
    fn test_chimeric_alignment_creation() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(donor, acceptor);

        assert_eq!(chim.total_score, 200);
    }

    #[test]
    fn test_chimeric_meets_min_segment_length() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(donor, acceptor);

        assert!(chim.meets_min_segment_length(20));
        assert!(chim.meets_min_segment_length(30));
        assert!(!chim.meets_min_segment_length(51));
    }

    #[test]
    fn test_chimeric_meets_min_score() {
        let donor = mock_segment(0, 1000, 1050, 0, 50);
        let acceptor = mock_segment(1, 2000, 2030, 50, 80);
        let chim = ChimericAlignment::new(donor, acceptor);

        assert!(chim.meets_min_score(100));
        assert!(chim.meets_min_score(200));
        assert!(!chim.meets_min_score(201));
    }


}
