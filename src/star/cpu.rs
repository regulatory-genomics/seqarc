//! Runtime CPU helpers. The compile-time build flavour (OS, arch, SIMD level,
//! git hash, build timestamp) is assembled into `VERSION_BODY` by `build.rs`
//! and read via `env!("VERSION_BODY")`.
//!
//! (Vendored trim: the runtime feature-detection / compatibility-guard helpers
//! — `binary_target`, `detected_features`, `upgrade_hint`, `check_cpu_compat`,
//! `cpu_detected_line` — are not used by seqarc and were removed.)

/// Software-prefetch a cache line for reading, hinting the CPU to pull `ptr`'s
/// line into L1 ahead of a later dependent load. Purely a performance hint —
/// architecturally a no-op, so it never affects results and is safe on any
/// pointer (the CPU ignores faults on prefetch). Used to overlap the random
/// mmap loads of the SA binary search with the current probe's genome compare.
#[inline]
pub fn prefetch_read(ptr: *const u8) {
    #[cfg(target_arch = "x86_64")]
    {
        // SAFETY: _mm_prefetch is a hint; it never dereferences/faults.
        unsafe {
            core::arch::x86_64::_mm_prefetch(ptr.cast::<i8>(), core::arch::x86_64::_MM_HINT_T0);
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        // SAFETY: `prfm` is a hint instruction; it never faults on a bad address.
        unsafe {
            core::arch::asm!("prfm pldl1keep, [{0}]", in(reg) ptr, options(nostack, preserves_flags));
        }
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        let _ = ptr;
    }
}
