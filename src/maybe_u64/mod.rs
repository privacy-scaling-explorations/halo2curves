use ff::PrimeField;
use rand_core::RngCore;
#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};
use subtle::CtOption;

mod arithmetic;
mod conversion;
mod field;
// mod field_ext;
mod misc;
mod prime_field;
mod serdes;

#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, Eq)]
#[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
/// Struct for a 256-bits field elements.
/// - U64: if the element is less than 1<<64
/// - Full: otherwise
///
/// This allows for faster computation with u64 ops.
pub enum MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    U64(u64),
    Full(F),
}

impl<F> MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    /// Converts an element of `Fr` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 32] {
        <Self as ff::PrimeField>::to_repr(self)
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fr`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<Self> {
        <Self as ff::PrimeField>::from_repr(*bytes)
    }
}

pub trait MaybeU64Coversion {
    /// Convert a u64 to a full field element
    fn to_full(&self) -> Self;
    /// Convert a full field element into u64 format.
    fn to_u64(&self) -> Self;
    /// random u64
    fn random_u64(rng: impl RngCore) -> Self;
    /// random field element
    fn random_field(rng: impl RngCore) -> Self;
}
