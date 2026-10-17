//! Vendored subset of [rustar-aligner](https://github.com/scverse/rustar-aligner)
//! (STAR-like spliced aligner), copied at revision `207773b11df5fb9c229bbbcf5385f6daa7f2fe7c`
//! (crate version 0.2.0) and inlined to drop the git dependency. The module is
//! named `star` after the aligner this reimplements (and the `SRSTAR01` index
//! container magic); upstream's crate is named `rustar-aligner`.
//!
//! Only the modules reachable from seqarc's index build (`star_index.rs`)
//! and read mapping (`aligners.rs`) paths are vendored; the CLI, solo/
//! cell-counting, SAM/BAM output, and other unrelated upstream modules were
//! dropped. Upstream `params/mod.rs` was trimmed accordingly (solo CLI
//! validation), and `io/mod.rs` was trimmed to `pub mod fastq`.
//!
//! Copyright (c) James Ferguson and rustar-aligner contributors.
//! This subset is distributed under the MIT license:
//! Permission is hereby granted, free of charge, to any person obtaining a
//! copy of this software and associated documentation files (the "Software"),
//! to deal in the Software without restriction, including without limitation
//! the rights to use, copy, modify, merge, publish, distribute, sublicense,
//! and/or sell copies of the Software, and to permit persons to whom the
//! Software is furnished to do so, subject to the following conditions: The
//! above copyright notice and this permission notice shall be included in all
//! copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED
//! "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
//! NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
//! PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
//! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//! ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

pub(crate) mod align;
pub(crate) mod chimeric;
pub(crate) mod cpu;
pub(crate) mod error;
pub(crate) mod genome;
pub(crate) mod index;
pub(crate) mod io;
pub(crate) mod junction;
pub(crate) mod params;
pub(crate) mod quant;
pub(crate) mod rng;
pub(crate) mod stats;
