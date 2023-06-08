use ff::{Field, PrimeField};
use rand_core::RngCore;
use subtle::{Choice, CtOption};

use crate::MaybeU64;

impl<F> Field for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    const ZERO: Self = Self::U64(0);
    const ONE: Self = Self::U64(1);

    fn random(rng: impl RngCore) -> Self {
        Self::Full(F::random(rng))
    }

    fn square(&self) -> Self {
        self * self
    }

    fn double(&self) -> Self {
        self + self
    }

    fn invert(&self) -> CtOption<Self> {
        match *self {
            MaybeU64::U64(a) => F::from(a).invert().map(|x| MaybeU64::Full(x)),
            MaybeU64::Full(a) => a.invert().map(|x| MaybeU64::Full(x)),
        }
    }

    fn sqrt(&self) -> CtOption<Self> {
        match *self {
            MaybeU64::U64(a) => F::from(a).sqrt().map(|x| MaybeU64::Full(x)),
            MaybeU64::Full(a) => a.sqrt().map(|x| MaybeU64::Full(x)),
        }
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        let (num, div) = match (num, div) {
            (MaybeU64::U64(a), MaybeU64::U64(b)) => (F::from(*a), F::from(*b)),
            (MaybeU64::U64(a), MaybeU64::Full(b)) => (F::from(*a), *b),
            (MaybeU64::Full(a), MaybeU64::U64(b)) => (*a, F::from(*b)),
            (MaybeU64::Full(a), MaybeU64::Full(b)) => (*a, *b),
        };
        let (c, f) = F::sqrt_ratio(&num, &div);

        (c, Self::Full(f))
    }
}
