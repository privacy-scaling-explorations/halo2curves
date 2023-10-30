#![cfg_attr(feature = "asm", feature(asm_const))]

mod arithmetic;

pub mod bn256;
pub mod hash_to_curve;
pub mod pairing;
pub mod pasta;
pub mod pluto_eris;
pub mod secp256k1;
pub mod serde;
#[macro_use]
pub mod legendre;

#[macro_use]
mod derive;
pub use arithmetic::CurveAffineExt;
pub use pasta_curves::arithmetic::{Coordinates, CurveAffine, CurveExt};

// Re-export ff and group to simplify down stream dependencies
#[cfg(feature = "reexport")]
pub use ff;
#[cfg(not(feature = "reexport"))]
use ff;
#[cfg(feature = "reexport")]
pub use group;
#[cfg(not(feature = "reexport"))]
use group;

#[cfg(test)]
pub mod tests;

#[cfg(all(feature = "prefetch", target_arch = "x86_64"))]
#[inline(always)]
pub fn prefetch<T>(data: &[T], offset: usize) {
    use core::arch::x86_64::_mm_prefetch;
    unsafe {
        _mm_prefetch(
            data.as_ptr().add(offset) as *const i8,
            core::arch::x86_64::_MM_HINT_T0,
        );
    }
}
