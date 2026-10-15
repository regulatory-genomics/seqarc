//! Port of the QVZ lossy quality-value quantizer from Spring 1.1.1
//! (`Spring-1.1.1/src/qvz`).
//!
//! QVZ is an encoder-only transform: it trains per-column conditional
//! quantizers on a block of quality strings and replaces each quality value
//! with a reconstruction point chosen to stay near a target entropy (bits per
//! symbol). The decoder never sees the codebooks — the quantized values are
//! stored losslessly by the rest of the archive — so only the encoder side is
//! ported here. Dead Spring code (clustering, codebook serialization, the
//! arithmetic coder, ratio/fixed-MSE modes, non-MSE distortions) is omitted.
//!
//! Spring leaves its WELL-1024a RNG state uninitialized, so its output is not
//! reproducible; this port seeds the generator deterministically instead.

use crate::{Error, Result};

/// Number of representable quality symbols (Phred 0..=71 as Phred+33 bytes
/// 33..=104), matching Spring's `ALPHABET_SIZE`.
const ALPHABET_SIZE: usize = 72;

/// Lloyd-Max iteration cap, matching Spring's `QUANTIZER_MAX_ITER`.
const QUANTIZER_MAX_ITER: u32 = 100;

/// Maximum supported read length, matching Spring's `MAX_READS_PER_LINE`.
/// Training memory grows linearly with the column count (72 PMF tables per
/// column), so longer reads are rejected rather than silently degrading.
pub(crate) const MAX_READ_LENGTH: usize = 1022;

/// Sentinel used in alphabet reverse-index tables.
const SYMBOL_NOT_FOUND: u32 = u32::MAX;

// ---------------------------------------------------------------------------
// WELL-1024a pseudo-random generator (Spring's well.cpp)
// ---------------------------------------------------------------------------

/// WELL-1024a generator used for randomized lo/hi quantizer selection.
struct Well1024a {
    state: [u32; 32],
    n: usize,
    bit_output: u32,
    bits_left: u32,
}

impl Well1024a {
    /// Creates a generator with a fixed nonzero seed. Spring reads this state
    /// from uninitialized stack memory; a deterministic xorshift walk keeps
    /// the statistical behavior while making output reproducible.
    fn new() -> Self {
        let mut state = [0u32; 32];
        let mut x = 0x9E37_79B9u32;
        for slot in state.iter_mut() {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            *slot = x;
        }
        Self {
            state,
            n: 0,
            bit_output: 0,
            bits_left: 0,
        }
    }

    fn next_u32(&mut self) -> u32 {
        let s = &mut self.state;
        let n = self.n;
        let z0 = s[(n + 31) & 31];
        let v_m1 = s[(n + 3) & 31];
        let v_m2 = s[(n + 24) & 31];
        let v_m3 = s[(n + 10) & 31];
        let z1 = s[n] ^ (v_m1 ^ (v_m1 >> 8));
        let z2 = (v_m2 ^ (v_m2 << 19)) ^ (v_m3 ^ (v_m3 << 14));

        s[n] = z1 ^ z2;
        let n = (n + 31) & 31;
        s[n] = (z0 ^ (z0 << 11)) ^ (z1 ^ (z1 << 7)) ^ (z2 ^ (z2 << 13));
        self.n = n;
        s[n]
    }

    /// Returns `bits` random bits, amortizing generator calls (Spring's
    /// `well_1024a_bits`).
    fn bits(&mut self, bits: u32) -> u32 {
        let mask = (1u32 << bits) - 1;
        if self.bits_left < bits {
            self.bit_output = self.next_u32();
            self.bits_left = 32;
        }
        let value = self.bit_output & mask;
        self.bit_output >>= bits;
        self.bits_left -= bits;
        value
    }
}

// ---------------------------------------------------------------------------
// Alphabets (Spring's pmf.cpp alphabet helpers)
// ---------------------------------------------------------------------------

/// A sorted set of symbols with a reverse index for lookups.
#[derive(Clone, Debug)]
struct Alphabet {
    symbols: Vec<u8>,
    indexes: Vec<u32>,
}

impl Alphabet {
    /// Identity alphabet containing symbols `0..size`.
    fn new(size: usize) -> Self {
        Self::from_symbols((0..size as u8).collect())
    }

    /// Builds an alphabet from the given symbols, keeping them in order and
    /// resolving duplicate lookups to the last position (Spring's
    /// `alphabet_compute_index` overwrites indexes in symbol order).
    fn from_symbols(symbols: Vec<u8>) -> Self {
        let mut indexes = vec![SYMBOL_NOT_FOUND; 256];
        for (index, &symbol) in symbols.iter().enumerate() {
            indexes[symbol as usize] = index as u32;
        }
        Self { symbols, indexes }
    }

    fn size(&self) -> usize {
        self.symbols.len()
    }

    fn index_of(&self, symbol: u8) -> Option<usize> {
        let index = self.indexes[symbol as usize];
        (index != SYMBOL_NOT_FOUND).then_some(index as usize)
    }

    /// Sorted union of two sorted alphabets (Spring's `alphabet_union`).
    fn union(&self, other: &Alphabet) -> Alphabet {
        let mut symbols = Vec::with_capacity(self.size() + other.size());
        let mut i = 0;
        let mut j = 0;
        while i < self.size() && j < other.size() {
            if self.symbols[i] < other.symbols[j] {
                symbols.push(self.symbols[i]);
                i += 1;
            } else if self.symbols[i] == other.symbols[j] {
                symbols.push(self.symbols[i]);
                i += 1;
                j += 1;
            } else {
                symbols.push(other.symbols[j]);
                j += 1;
            }
        }
        symbols.extend_from_slice(&self.symbols[i..]);
        symbols.extend_from_slice(&other.symbols[j..]);
        Alphabet::from_symbols(symbols)
    }
}

// ---------------------------------------------------------------------------
// Probability mass functions (Spring's pmf.cpp)
// ---------------------------------------------------------------------------

/// A PMF that can either accumulate empirical counts or hold explicit
/// probabilities, mirroring Spring's `pmf_t` with its lazy recalculation.
#[derive(Clone, Debug)]
struct Pmf {
    probs: Vec<f64>,
    counts: Vec<u32>,
    total: u32,
    ready: bool,
}

impl Pmf {
    fn new(size: usize) -> Self {
        Self {
            probs: vec![0.0; size],
            counts: vec![0; size],
            total: 0,
            ready: false,
        }
    }

    fn size(&self) -> usize {
        self.probs.len()
    }

    fn increment(&mut self, index: usize) {
        self.counts[index] += 1;
        self.total += 1;
    }

    /// Converts accumulated counts into probabilities (Spring's
    /// `recalculate_pmf`).
    fn recalculate(&mut self) {
        self.ready = true;
        if self.total == 0 {
            return;
        }
        let total = f64::from(self.total);
        for (prob, &count) in self.probs.iter_mut().zip(&self.counts) {
            *prob = f64::from(count) / total;
        }
    }

    /// Shannon entropy in bits (Spring's `get_entropy`).
    fn entropy(&mut self) -> f64 {
        if !self.ready {
            self.recalculate();
        }
        -self
            .probs
            .iter()
            .filter(|probability| **probability > 0.0)
            .map(|probability| probability * probability.log2())
            .sum::<f64>()
    }

    /// Rescales nonzero probabilities to sum to one (Spring's
    /// `renormalize_pmf`).
    fn renormalize(&mut self) {
        let total: f64 = self.probs.iter().sum();
        if total > 0.0 {
            for prob in self.probs.iter_mut() {
                *prob /= total;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Distortion metric (Spring's distortion.cpp, MSE only)
// ---------------------------------------------------------------------------

/// Squared-error distortion matrix over the quality alphabet.
struct DistortionMatrix {
    symbols: usize,
    data: Vec<f64>,
}

impl DistortionMatrix {
    fn mse(symbols: usize) -> Self {
        let mut data = vec![0.0; symbols * symbols];
        for x in 0..symbols {
            for y in 0..symbols {
                let difference = (x as i32 - y as i32) as f64;
                data[x + y * symbols] = difference * difference;
            }
        }
        Self { symbols, data }
    }

    fn get(&self, x: u8, y: u8) -> f64 {
        self.data[x as usize + self.symbols * y as usize]
    }
}

// ---------------------------------------------------------------------------
// Quantizers (Spring's quantizer.cpp)
// ---------------------------------------------------------------------------

/// Maps input symbols to reconstruction points; mirrors Spring's
/// `quantizer_t`. `ratio` stores the probability that this quantizer is the
/// one selected during random lo/hi mixing.
#[derive(Clone, Debug)]
struct Quantizer {
    q: Vec<u8>,
    output_alphabet: Alphabet,
    ratio: f64,
    /// Expected distortion of this quantizer; surfaced by tests.
    #[allow(dead_code)]
    mse: f64,
}

impl Quantizer {
    /// Lloyd-Max quantizer design with the given number of output states
    /// (Spring's `generate_quantizer`).
    fn generate(pmf: &Pmf, dist: &DistortionMatrix, states: usize) -> Quantizer {
        debug_assert!(pmf.ready, "quantizer generation requires a ready PMF");
        let size = pmf.size();
        let mut bounds = vec![0usize; states + 1];
        bounds[states] = size;
        for (j, bound) in bounds.iter_mut().enumerate().skip(1).take(states) {
            *bound = (j * size) / states;
        }
        let mut reconstruction = vec![0u8; states];
        for j in 0..states {
            reconstruction[j] = ((bounds[j] + bounds[j + 1] - 1) / 2) as u8;
        }

        // Alternate between adjusting reconstruction points for fixed bounds
        // and bounds for fixed reconstruction points until convergence.
        let mut changed = true;
        let mut iter = 0;
        while changed && iter < QUANTIZER_MAX_ITER {
            changed = false;
            iter += 1;

            for j in 0..states {
                let mut min_mse = f64::MAX;
                let mut min_r = bounds[j];
                for r in bounds[j]..bounds[j + 1] {
                    let mut mse = 0.0;
                    for i in bounds[j]..bounds[j + 1] {
                        mse += pmf.probs[i] * dist.get(i as u8, r as u8);
                    }
                    if mse < min_mse {
                        min_r = r;
                        min_mse = mse;
                    }
                }
                if min_r as u8 != reconstruction[j] {
                    changed = true;
                    reconstruction[j] = min_r as u8;
                }
            }

            let mut r = 0usize;
            let mut j = 1usize;
            while j + 1 < size && r + 1 < states {
                let mse = dist.get(j as u8, reconstruction[r]);
                let next_mse = dist.get(j as u8, reconstruction[r + 1]);
                if next_mse < mse {
                    r += 1;
                    bounds[r] = j;
                }
                j += 1;
            }
        }

        let mut q = vec![0u8; size];
        let mut mse = 0.0;
        for (j, &reconstruction_point) in reconstruction.iter().enumerate() {
            let (start, end) = (bounds[j], bounds[j + 1]);
            q[start..end].fill(reconstruction_point);
            for (i, &probability) in pmf.probs[start..end].iter().enumerate() {
                mse += dist.get((start + i) as u8, reconstruction_point) * probability;
            }
        }

        Quantizer {
            q,
            output_alphabet: Alphabet::from_symbols(reconstruction),
            ratio: 0.0,
            mse,
        }
    }

    /// Output PMF when this quantizer is applied to the given input PMF
    /// (Spring's `apply_quantizer`). The output is over the input alphabet,
    /// indexed by reconstruction point.
    fn apply(&self, pmf: &Pmf) -> Pmf {
        debug_assert!(pmf.ready, "quantizer application requires a ready PMF");
        let mut probs = vec![0.0; pmf.size()];
        for (index, &probability) in pmf.probs.iter().enumerate() {
            probs[self.q[index] as usize] += probability;
        }
        Pmf {
            probs,
            counts: vec![0; pmf.size()],
            total: 0,
            ready: true,
        }
    }
}

/// Brackets the target entropy with a low and a high quantizer and returns
/// them together with the mixing ratio (Spring's `optimize_for_entropy`).
fn optimize_for_entropy(
    pmf: &mut Pmf,
    dist: &DistortionMatrix,
    target: f64,
) -> (Quantizer, Quantizer, f64) {
    if target == 0.0 {
        let quantizer = Quantizer::generate(pmf, dist, 1);
        return (quantizer.clone(), quantizer, 1.0);
    }

    let size = pmf.size();

    // Applying a quantizer only groups input probabilities, so the output
    // entropy can never exceed the input entropy. When even the input entropy
    // is below the target, the state search is guaranteed to run through every
    // state count and settle on the last pair with a zero ratio; jump straight
    // there instead of generating all intermediate quantizers.
    if size >= 2 && pmf.entropy() < target {
        let lo = Quantizer::generate(pmf, dist, size - 1);
        let hi = Quantizer::generate(pmf, dist, size);
        return (lo, hi, 0.0);
    }

    let mut states = 1usize;
    let mut hi = Quantizer::generate(pmf, dist, states);
    let mut hi_entropy = hi.apply(pmf).entropy();

    let mut lo;
    let mut lo_entropy;
    loop {
        lo = hi;
        lo_entropy = hi_entropy;
        states += 1;
        hi = Quantizer::generate(pmf, dist, states);
        hi_entropy = hi.apply(pmf).entropy();
        if !(hi_entropy < target && states < pmf.size()) {
            break;
        }
    }

    let ratio = if hi_entropy < target {
        0.0
    } else if lo_entropy >= target || hi_entropy == lo_entropy {
        1.0
    } else {
        (target - hi_entropy) / (lo_entropy - hi_entropy)
    };
    (lo, hi, ratio)
}

// ---------------------------------------------------------------------------
// Conditional statistics (Spring's codebook.cpp calculate_statistics)
// ---------------------------------------------------------------------------

/// Conditional PMFs for each column given the previous column's value, plus
/// the marginal PMF per column. Mirrors Spring's `cond_pmf_list_t`.
struct CondPmfList {
    columns: usize,
    pmfs: Vec<Pmf>,
    marginal: Vec<Pmf>,
}

impl CondPmfList {
    fn new(columns: usize) -> Self {
        let count = 1 + ALPHABET_SIZE * (columns - 1);
        Self {
            columns,
            pmfs: (0..count).map(|_| Pmf::new(ALPHABET_SIZE)).collect(),
            marginal: Vec::new(),
        }
    }

    fn conditional(&self, column: usize, prev: u8) -> &Pmf {
        &self.pmfs[self.conditional_index(column, prev)]
    }

    /// The 72 conditional PMFs for `column`, one per possible previous symbol.
    fn conditionals(&self, column: usize) -> &[Pmf] {
        if column == 0 {
            &self.pmfs[..1]
        } else {
            let start = 1 + (column - 1) * ALPHABET_SIZE;
            &self.pmfs[start..start + ALPHABET_SIZE]
        }
    }

    fn conditional_mut(&mut self, column: usize, prev: u8) -> &mut Pmf {
        let index = self.conditional_index(column, prev);
        &mut self.pmfs[index]
    }

    fn conditional_index(&self, column: usize, prev: u8) -> usize {
        if column == 0 {
            0
        } else {
            1 + (column - 1) * ALPHABET_SIZE + prev as usize
        }
    }

    /// Computes the unconditional PMF per column from the conditional ones
    /// via the law of total probability, and finalizes every PMF so that
    /// later stages can read probabilities without mutation.
    fn compute_marginals(&mut self) {
        for pmf in self.pmfs.iter_mut() {
            pmf.recalculate();
        }

        let first = {
            let mut pmf = Pmf::new(ALPHABET_SIZE);
            pmf.probs.copy_from_slice(&self.pmfs[0].probs);
            pmf.ready = true;
            pmf
        };
        self.marginal.push(first);

        for column in 1..self.columns {
            let previous = self.marginal[column - 1].probs.clone();
            let mut probs = vec![0.0; ALPHABET_SIZE];
            for (j, &weight) in previous.iter().enumerate() {
                let conditional = &self.pmfs[self.conditional_index(column, j as u8)];
                for (k, slot) in probs.iter_mut().enumerate() {
                    *slot += weight * conditional.probs[k];
                }
            }
            self.marginal.push(Pmf {
                probs,
                counts: vec![0; ALPHABET_SIZE],
                total: 0,
                ready: true,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// Conditional quantizer list (Spring's codebook.cpp generate_codebooks)
// ---------------------------------------------------------------------------

/// Per-column quantizers indexed by the previous column's reconstruction
/// symbol, mirroring Spring's `cond_quantizer_list_t`. For each column and
/// each left-context symbol there is a low and a high quantizer mixed
/// according to a stored quantized ratio.
struct CondQuantizerList {
    #[allow(dead_code)]
    columns: usize,
    input_alphabets: Vec<Alphabet>,
    quantizers: Vec<Vec<Quantizer>>,
    qratio: Vec<Vec<u8>>,
}

impl CondQuantizerList {
    /// Selects the quantizer for the given column and left context using a
    /// random draw against the quantized ratio (Spring's `choose_quantizer`).
    fn choose_quantizer<'a>(&'a self, well: &mut Well1024a, column: usize, prev: u8) -> &'a Quantizer {
        let index = self.input_alphabets[column]
            .index_of(prev)
            .expect("previous quantizer value is present in the input alphabet");
        if well.bits(7) >= u32::from(self.qratio[column][index]) {
            &self.quantizers[column][2 * index + 1]
        } else {
            &self.quantizers[column][2 * index]
        }
    }

    fn push_column(&mut self, input_alphabet: Alphabet, quantizers: Vec<Quantizer>, ratios: &[f64]) {
        debug_assert_eq!(quantizers.len(), 2 * ratios.len());
        let qratio = ratios
            .iter()
            .map(|&ratio| (ratio * 128.0) as u8)
            .collect();
        self.input_alphabets.push(input_alphabet);
        self.quantizers.push(quantizers);
        self.qratio.push(qratio);
    }
}

/// Trains per-column conditional quantizers for one cluster of reads
/// (Spring's `generate_codebooks` with a single cluster and `MODE_FIXED`).
fn generate_codebooks(
    in_pmfs: &CondPmfList,
    dist: &DistortionMatrix,
    target: f64,
    columns: usize,
) -> CondQuantizerList {
    let mut list = CondQuantizerList {
        columns,
        input_alphabets: Vec::with_capacity(columns),
        quantizers: Vec::with_capacity(columns),
        qratio: Vec::with_capacity(columns),
    };

    // Column 0 has no left context: one lo/hi pair.
    let (q_lo, q_hi, ratio) = {
        let mut pmf = in_pmfs.pmfs[0].clone();
        optimize_for_entropy(&mut pmf, dist, target)
    };
    let mut q_lo = q_lo;
    let mut q_hi = q_hi;
    q_lo.ratio = ratio;
    q_hi.ratio = 1.0 - ratio;
    list.push_column(Alphabet::new(1), vec![q_lo.clone(), q_hi.clone()], &[ratio]);

    // P(Q | X) for the previous column; column 1 derives it directly from
    // the column-0 pair, later columns from the previous column's PMFs.
    let mut prev_qpmf: Vec<Pmf> = Vec::new();

    for column in 1..columns {
        // Union of every output symbol produced by the previous column.
        let union = {
            let previous = &list.quantizers[column - 1];
            let mut union = previous[0].output_alphabet.clone();
            for quantizer in &previous[1..] {
                union = union.union(&quantizer.output_alphabet);
            }
            union
        };

        // P(Q_i | X_i) over the new output union.
        let qpmf = if column == 1 {
            qpmf_from_quantizer_pair(&q_lo, &q_hi, ratio, &union)
        } else {
            compute_qpmf_list(in_pmfs, column, &prev_qpmf, &union, &list)
        };

        // P(X_{i+1} | Q_i) for each output symbol of this column.
        let mut xpmf = compute_xpmf_list(in_pmfs, column, &qpmf, &union);

        // One lo/hi quantizer pair per possible left context.
        let mut quantizers = Vec::with_capacity(2 * union.size());
        let mut ratios = Vec::with_capacity(union.size());
        for pmf in xpmf.iter_mut() {
            let (mut lo, mut hi, ratio) = optimize_for_entropy(pmf, dist, target);
            lo.ratio = ratio;
            hi.ratio = 1.0 - ratio;
            ratios.push(ratio);
            quantizers.push(lo);
            quantizers.push(hi);
        }
        list.push_column(union.clone(), quantizers, &ratios);

        prev_qpmf = qpmf;
    }

    list
}

/// Builds P(Q_1 | X_1) from the single column-0 quantizer pair (Spring's
/// `compute_qpmf_quan_list`). The raw lo/hi mixture weights already sum to
/// one per row, so no renormalization is applied, matching Spring.
fn qpmf_from_quantizer_pair(
    q_lo: &Quantizer,
    q_hi: &Quantizer,
    ratio: f64,
    union: &Alphabet,
) -> Vec<Pmf> {
    let mut qpmf = (0..ALPHABET_SIZE).map(|_| Pmf::new(union.size())).collect::<Vec<_>>();
    for (k, pmf) in qpmf.iter_mut().enumerate() {
        for (idx, &symbol) in union.symbols.iter().enumerate() {
            if q_lo.q[k] == symbol {
                pmf.probs[idx] += ratio;
            }
            if q_hi.q[k] == symbol {
                pmf.probs[idx] += 1.0 - ratio;
            }
        }
        pmf.ready = true;
    }
    qpmf
}

/// Computes P(Q_i | X_i) for column `column >= 2` (Spring's
/// `compute_qpmf_list`). The result holds one PMF per input symbol over the
/// union of this column's outputs.
///
/// Spring recomputes the `p_temp` sum inside the union-index loop, but it
/// only depends on the input symbol and the previous context index, so it is
/// hoisted into a table here. Accumulation order is unchanged.
fn compute_qpmf_list(
    in_pmfs: &CondPmfList,
    column: usize,
    prev_qpmf: &[Pmf],
    union: &Alphabet,
    list: &CondQuantizerList,
) -> Vec<Pmf> {
    let previous = &list.quantizers[column - 1];
    let contexts = &list.input_alphabets[column - 1];

    // p_temp[k][j] = sum_x P(Q_{c-1}=j | X_{c-1}=x) * P(X_c=k | X_{c-1}=x) * P(X_{c-1}=x)
    let previous_conditionals = in_pmfs.conditionals(column - 1);
    let marginal = &in_pmfs.marginal[column - 2].probs;
    let mut p_temp = vec![0.0f64; ALPHABET_SIZE * contexts.size()];
    for (k, row) in p_temp.chunks_mut(contexts.size()).enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            let mut accumulated = 0.0;
            for ((prev_pmf, conditional), &marginal_weight) in prev_qpmf
                .iter()
                .zip(previous_conditionals)
                .zip(marginal)
            {
                accumulated += prev_pmf.probs[j] * marginal_weight * conditional.probs[k];
            }
            *slot = accumulated;
        }
    }

    let mut qpmf = (0..ALPHABET_SIZE).map(|_| Pmf::new(union.size())).collect::<Vec<_>>();
    for (k, pmf) in qpmf.iter_mut().enumerate() {
        for (idx, &symbol) in union.symbols.iter().enumerate() {
            let mut accumulated = 0.0;
            for (j, pair) in previous.as_chunks::<2>().0.iter().enumerate() {
                let q_lo = &pair[0];
                let q_hi = &pair[1];
                let mut p_q_xq = 0.0;
                if q_lo.q[k] == symbol {
                    p_q_xq += q_lo.ratio;
                }
                if q_hi.q[k] == symbol {
                    p_q_xq += q_hi.ratio;
                }
                accumulated += p_q_xq * p_temp[k * contexts.size() + j];
            }
            pmf.probs[idx] = accumulated;
        }
        pmf.ready = true;
        pmf.renormalize();
    }
    qpmf
}

/// Computes P(X_{i+1} | Q_i) for each output symbol of the union (Spring's
/// `compute_xpmf_list`). The result holds one PMF per union symbol over the
/// full input alphabet.
fn compute_xpmf_list(
    in_pmfs: &CondPmfList,
    column: usize,
    qpmf: &[Pmf],
    union: &Alphabet,
) -> Vec<Pmf> {
    let mut xpmf = (0..union.size()).map(|_| Pmf::new(ALPHABET_SIZE)).collect::<Vec<_>>();

    for (idx, pmf) in xpmf.iter_mut().enumerate() {
        for (k, slot) in pmf.probs.iter_mut().enumerate() {
            let mut accumulated = 0.0;
            for (x, q_pmf) in qpmf.iter().enumerate() {
                accumulated += q_pmf.probs[idx]
                    * in_pmfs.conditional(column, x as u8).probs[k]
                    * in_pmfs.marginal[column - 1].probs[x];
            }
            *slot = accumulated;
        }
        pmf.ready = true;
        pmf.renormalize();
    }
    xpmf
}

// ---------------------------------------------------------------------------
// Public entry point (Spring's quantize_quality_qvz + start_qv_quantization)
// ---------------------------------------------------------------------------

/// Quantizes a block of Phred+33 quality rows in place toward `ratio` bits
/// per symbol, mirroring Spring's `quantize_quality_qvz` (single cluster,
/// fixed-rate mode, MSE distortion).
pub(crate) fn quantize(rows: &mut [&mut [u8]], ratio: f64) -> Result<()> {
    if rows.is_empty() {
        return Ok(());
    }
    if ratio <= 0.0 || ratio.is_nan() {
        return Err(Error::InvalidInput(
            "qvz quality mode requires a positive ratio".into(),
        ));
    }

    let mut columns = 0usize;
    for row in rows.iter() {
        columns = columns.max(row.len());
        for &byte in row.iter() {
            if !(33..=104).contains(&byte) {
                return Err(Error::InvalidInput(format!(
                    "qvz quality mode requires Phred+33 scores in the range 33..=104 \
                     (Phred 0..=71); found byte {byte}"
                )));
            }
        }
    }
    if columns == 0 {
        return Ok(());
    }
    if columns > MAX_READ_LENGTH {
        return Err(Error::InvalidInput(format!(
            "qvz quality mode supports read lengths up to {MAX_READ_LENGTH} bases; \
             found {columns}"
        )));
    }

    let dist = DistortionMatrix::mse(ALPHABET_SIZE);

    let mut stats = CondPmfList::new(columns);
    for row in rows.iter() {
        if row.is_empty() {
            continue;
        }
        stats
            .conditional_mut(0, 0)
            .increment((row[0] - 33) as usize);
        for column in 1..row.len() {
            stats
                .conditional_mut(column, row[column - 1] - 33)
                .increment((row[column] - 33) as usize);
        }
    }
    stats.compute_marginals();

    let quantizers = generate_codebooks(&stats, &dist, ratio, columns);

    let mut well = Well1024a::new();
    for row in rows.iter_mut() {
        let mut previous = 0u8;
        for (column, byte) in row.iter_mut().enumerate() {
            let quantizer = quantizers.choose_quantizer(&mut well, column, previous);
            let value = quantizer.q[(*byte - 33) as usize];
            *byte = value + 33;
            previous = value;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform_pmf(size: usize) -> Pmf {
        let mut pmf = Pmf::new(size);
        for index in 0..size {
            pmf.increment(index);
        }
        pmf.recalculate();
        pmf
    }

    #[test]
    fn well_generator_is_deterministic_and_bounded() {
        let mut a = Well1024a::new();
        let mut b = Well1024a::new();
        for _ in 0..1000 {
            assert_eq!(a.next_u32(), b.next_u32());
        }
        let mut well = Well1024a::new();
        for _ in 0..1000 {
            assert!(well.bits(7) < 128);
        }
    }

    #[test]
    fn alphabet_union_merges_sorted_symbols() {
        let a = Alphabet::from_symbols(vec![1, 3, 5]);
        let b = Alphabet::from_symbols(vec![2, 3, 4]);
        let union = a.union(&b);
        assert_eq!(union.symbols, vec![1, 2, 3, 4, 5]);
        assert_eq!(union.index_of(4), Some(3));
        assert_eq!(union.index_of(9), None);
        assert_eq!(Alphabet::new(1).symbols, vec![0]);
    }

    #[test]
    fn pmf_entropy_matches_hand_computed_values() {
        let mut pmf = Pmf::new(4);
        for _ in 0..4 {
            pmf.increment(0);
            pmf.increment(1);
        }
        pmf.increment(2);
        // 4/9, 4/9, 1/9, 0
        let entropy = pmf.entropy();
        let expected = -(2.0 * (4.0 / 9.0) * (4.0f64 / 9.0).log2()
            + (1.0 / 9.0) * (1.0f64 / 9.0).log2());
        assert!((entropy - expected).abs() < 1e-12);
    }

    #[test]
    fn quantizer_of_uniform_pmf_has_monotone_distortion() {
        let dist = DistortionMatrix::mse(ALPHABET_SIZE);
        let pmf = uniform_pmf(ALPHABET_SIZE);
        let mut previous_mse = f64::MAX;
        for states in 1..=ALPHABET_SIZE {
            let quantizer = Quantizer::generate(&pmf, &dist, states);
            assert!(quantizer.mse <= previous_mse + 1e-9);
            assert_eq!(quantizer.output_alphabet.size(), states);
            previous_mse = quantizer.mse;
        }
        // A uniform distribution quantized to every state collapses to identity.
        let identity = Quantizer::generate(&pmf, &dist, ALPHABET_SIZE);
        assert_eq!(identity.q, (0..ALPHABET_SIZE as u8).collect::<Vec<_>>());
    }

    #[test]
    fn apply_quantizer_groups_input_probabilities() {
        let dist = DistortionMatrix::mse(4);
        let mut pmf = Pmf::new(4);
        for index in 0..4 {
            pmf.increment(index);
        }
        pmf.recalculate();
        let quantizer = Quantizer::generate(&pmf, &dist, 2);
        let output = quantizer.apply(&pmf);
        let total: f64 = output.probs.iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
        // Only the two reconstruction points carry probability.
        let used = output.probs.iter().filter(|p| **p > 0.0).count();
        assert_eq!(used, 2);
    }

    #[test]
    fn optimize_for_entropy_brackets_target() {
        let dist = DistortionMatrix::mse(ALPHABET_SIZE);
        let mut pmf = uniform_pmf(ALPHABET_SIZE);
        for target in [0.5f64, 1.0, 2.0, 4.0] {
            let (lo, hi, ratio) = optimize_for_entropy(&mut pmf, &dist, target);
            assert!((0.0..=1.0).contains(&ratio));
            let lo_entropy = lo.apply(&pmf).entropy();
            let hi_entropy = hi.apply(&pmf).entropy();
            // The hi quantizer is the first to reach or exceed the target.
            assert!(hi_entropy >= target - 1e-9 || hi.output_alphabet.size() == ALPHABET_SIZE);
            assert!(lo_entropy <= hi_entropy + 1e-9);
        }
        // A zero target collapses to a single shared one-state quantizer.
        let (lo, hi, ratio) = optimize_for_entropy(&mut pmf, &dist, 0.0);
        assert_eq!(ratio, 1.0);
        assert_eq!(lo.q, hi.q);
    }

    /// Small deterministic LCG so test data has noise like real quality
    /// scores while staying reproducible.
    fn lcg(state: &mut u64) -> u64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *state >> 33
    }

    fn synthetic_rows(rows: usize, columns: usize, seed: u64) -> Vec<Vec<u8>> {
        let mut state = seed;
        (0..rows)
            .map(|_| {
                (0..columns)
                    .map(|_| 33 + 2 + (lcg(&mut state) % 38) as u8)
                    .collect()
            })
            .collect()
    }

    #[test]
    fn quantize_reduces_entropy_and_stays_in_alphabet() {
        let mut rows = synthetic_rows(96, 48, 0x5EED);
        let original = rows.clone();

        let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        quantize(&mut slices, 1.0).unwrap();
        for row in rows.iter() {
            assert_eq!(row.len(), 48);
            assert!(row.iter().all(|byte| (33..=104).contains(byte)));
        }

        // QVZ quantizes conditionally on the previous column's reconstruction
        // point, so the global per-column symbol set need not be tiny; the
        // gain shows up when the stream is entropy coded.
        let flat = |rows: &[Vec<u8>]| rows.concat();
        let before = zstd::stream::encode_all(flat(&original).as_slice(), 3).unwrap();
        let after = zstd::stream::encode_all(flat(&rows).as_slice(), 3).unwrap();
        assert!(
            after.len() * 10 < before.len() * 9,
            "quantized stream should compress better: {} -> {}",
            before.len(),
            after.len()
        );
    }

    #[test]
    fn quantize_is_deterministic() {
        let mut first = synthetic_rows(48, 40, 0xC0FFEE);
        let mut second = synthetic_rows(48, 40, 0xC0FFEE);
        let mut a = first.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        let mut b = second.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        quantize(&mut a, 1.0).unwrap();
        quantize(&mut b, 1.0).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn quantize_rejects_invalid_ratio_scores_and_lengths() {
        let mut rows = [vec![b'I'; 10]];
        let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        assert!(quantize(&mut slices, 0.0).is_err());
        assert!(quantize(&mut slices, f64::NAN).is_err());

        let mut rows = [vec![32; 10]];
        let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        assert!(quantize(&mut slices, 1.0).is_err());
        let mut rows = [vec![105; 10]];
        let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        assert!(quantize(&mut slices, 1.0).is_err());

        let long = MAX_READ_LENGTH + 1;
        let mut rows = [vec![b'I'; long]];
        let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        let error = quantize(&mut slices, 1.0).unwrap_err();
        assert!(error.to_string().contains("1022"));
    }

    #[test]
    fn quantize_handles_empty_blocks_and_rows() {
        let mut empty: Vec<&mut [u8]> = Vec::new();
        quantize(&mut empty, 1.0).unwrap();
        let mut rows: Vec<Vec<u8>> = vec![Vec::new(), vec![b'I'; 3]];
        let mut slices = rows.iter_mut().map(Vec::as_mut_slice).collect::<Vec<_>>();
        quantize(&mut slices, 1.0).unwrap();
        assert!(rows[0].is_empty());
        assert_eq!(rows[1].len(), 3);
    }
}
