//! Implements     
//! - ConditionallySelectable
//! - ConstantTimeEq
//! - Default

use ff::PrimeField;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

use crate::MaybeU64;

impl<F> Default for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    fn default() -> Self {
        MaybeU64::U64(0)
    }
}

impl<F> ConditionallySelectable for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    fn conditional_select(a: &Self, b: &Self, c: Choice) -> Self {
        // WARNING: not constant time
        if c.unwrap_u8() == 0 {
            *a
        } else {
            *b
        }
    }
}

impl<F> ConstantTimeEq for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    fn ct_eq(&self, rhs: &Self) -> Choice {
        // WARNING: not constant time
        match (self, rhs) {
            (MaybeU64::U64(a), MaybeU64::U64(b)) => a.ct_eq(b),
            (MaybeU64::U64(a), MaybeU64::Full(b)) => F::from(*a).ct_eq(b),
            (MaybeU64::Full(a), MaybeU64::U64(b)) => F::from(*b).ct_eq(a),
            (MaybeU64::Full(a), MaybeU64::Full(b)) => a.ct_eq(b),
        }
    }
}

impl<F> PartialEq<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    fn eq(&self, rhs: &MaybeU64<F>) -> bool {
        self.ct_eq(rhs).into()
    }
}
