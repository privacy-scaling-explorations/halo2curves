use std::convert::TryInto;

use ff::PrimeField;
use rand_core::RngCore;

use crate::MaybeU64;

use super::MaybeU64Coversion;

impl<F> MaybeU64Coversion for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    /// Convert a u64 to a full field element
    fn to_full(&self) -> Self {
        match *self {
            MaybeU64::Full(p) => panic!("{:?} is already a full element", p),
            MaybeU64::U64(p) => Self::Full(F::from(p)),
        }
    }

    /// Convert a full field element into u64 format.
    fn to_u64(&self) -> Self {
        match *self {
            MaybeU64::Full(p) => {
                let repr = p.to_repr();
                for e in repr.iter().skip(8) {
                    assert_eq!(*e, 0, "{:?} is not an u64", p)
                }
                MaybeU64::U64(u64::from_le_bytes(repr[0..8].try_into().unwrap()))
            }
            MaybeU64::U64(p) => panic!("{} is already a u64", p),
        }
    }

    /// random u64
    fn random_u64(mut rng: impl RngCore) -> Self {
        Self::U64(rng.next_u64())
    }

    /// random field element
    fn random_field(rng: impl RngCore) -> Self {
        Self::Full(F::random(rng))
    }
}

impl<F> From<u64> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    fn from(value: u64) -> Self {
        Self::U64(value)
    }
}

impl<F> From<bool> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    fn from(value: bool) -> Self {
        Self::U64(value as u64)
    }
}
