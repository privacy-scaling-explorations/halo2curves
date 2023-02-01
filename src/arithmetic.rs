//! This module provides common utilities, traits and structures for group and
//! field arithmetic.
//!
//! This module is temporary, and the extension traits defined here are expected to be
//! upstreamed into the `ff` and `group` crates after some refactoring.

pub trait CurveAffineExt: pasta_curves::arithmetic::CurveAffine {
    fn batch_add<const COMPLETE: bool, const LOAD_POINTS: bool>(
        points: &mut [Self],
        output_indices: &[u32],
        num_points: usize,
        offset: usize,
        bases: &[Self],
        base_positions: &[u32],
    );

    /// Unlike the `Coordinates` trait, this just returns the raw affine coordinates without checking `is_on_curve`
    fn into_coordinates(self) -> (Self::Base, Self::Base) {
        // fallback implementation
        let coordinates = self.coordinates().unwrap();
        (*coordinates.x(), *coordinates.y())
    }
}

/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub(crate) const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub(crate) const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub(crate) const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}
