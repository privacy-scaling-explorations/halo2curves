use std::cmp::Ordering;
use std::convert::TryInto;
use std::iter::{Product, Sum};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use ff::{Field, PrimeField};

use crate::arithmetic::{adc, mac, sbb};
use crate::MaybeU64;

impl<'a, 'b, F> Add<&'b MaybeU64<F>> for &'a MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn add(self, rhs: &'b MaybeU64<F>) -> Self::Output {
        match (self, rhs) {
            (MaybeU64::U64(a), MaybeU64::U64(b)) => {
                let (c, carry) = adc(*a, *b, 0);
                match carry {
                    0 => MaybeU64::U64(c),
                    1 => MaybeU64::Full(
                        F::from_repr(
                            [
                                c.to_le_bytes().as_ref(),
                                &[carry as u8],
                                vec![0; 23].as_ref(),
                            ]
                            .concat()
                            .try_into()
                            .unwrap(),
                        )
                        .unwrap(),
                    ),
                    _ => panic!("invalid carry: {}", carry),
                }
            }
            (MaybeU64::U64(a), MaybeU64::Full(b)) => MaybeU64::Full(F::from(*a) + b),
            (MaybeU64::Full(a), MaybeU64::U64(b)) => MaybeU64::Full(*a + F::from(*b)),
            (MaybeU64::Full(a), MaybeU64::Full(b)) => MaybeU64::Full(*a + *b),
        }
    }
}

impl<'a, F> Add<MaybeU64<F>> for &'a MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn add(self, rhs: MaybeU64<F>) -> Self::Output {
        self + &rhs
    }
}

impl<'b, F> Add<&'b MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn add(self, rhs: &'b MaybeU64<F>) -> Self::Output {
        &self + rhs
    }
}

impl<F> Add<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn add(self, rhs: MaybeU64<F>) -> Self::Output {
        &self + &rhs
    }
}

impl<F> AddAssign<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    #[inline]
    fn add_assign(&mut self, rhs: MaybeU64<F>) {
        *self = &*self + &rhs;
    }
}

impl<'b, F> AddAssign<&'b MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    #[inline]
    fn add_assign(&mut self, rhs: &'b MaybeU64<F>) {
        *self = &*self + rhs;
    }
}

impl<F, T> Sum<T> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
    T: core::borrow::Borrow<Self>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::ZERO, |acc, item| acc + item.borrow())
    }
}

impl<'a, 'b, F> Sub<&'b MaybeU64<F>> for &'a MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn sub(self, rhs: &'b MaybeU64<F>) -> Self::Output {
        match (self, rhs) {
            (MaybeU64::U64(a), MaybeU64::U64(b)) => {
                let (c, carry) = sbb(*a, *b, 0);
                match carry {
                    0 => MaybeU64::U64(c),
                    u64::MAX => MaybeU64::Full(F::from(*a) - F::from(*b)),
                    _ => panic!("invalid carry: {}", carry),
                }
            }
            (MaybeU64::U64(a), MaybeU64::Full(b)) => MaybeU64::Full(F::from(*a) - b),
            (MaybeU64::Full(a), MaybeU64::U64(b)) => MaybeU64::Full(*a - F::from(*b)),
            (MaybeU64::Full(a), MaybeU64::Full(b)) => MaybeU64::Full(*a - *b),
        }
    }
}

impl<'a, F> Sub<MaybeU64<F>> for &'a MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn sub(self, rhs: MaybeU64<F>) -> Self::Output {
        self - &rhs
    }
}

impl<'b, F> Sub<&'b MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn sub(self, rhs: &'b MaybeU64<F>) -> Self::Output {
        &self - rhs
    }
}

impl<F> Sub<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn sub(self, rhs: MaybeU64<F>) -> Self::Output {
        &self - &rhs
    }
}

impl<F> SubAssign<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    #[inline]
    fn sub_assign(&mut self, rhs: MaybeU64<F>) {
        *self = &*self - &rhs;
    }
}

impl<'b, F> SubAssign<&'b MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    #[inline]
    fn sub_assign(&mut self, rhs: &'b MaybeU64<F>) {
        *self = &*self - rhs;
    }
}

impl<'a, 'b, F> Mul<&'b MaybeU64<F>> for &'a MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn mul(self, rhs: &'b MaybeU64<F>) -> Self::Output {
        match (self, rhs) {
            (MaybeU64::U64(a), MaybeU64::U64(b)) => {
                let (c, carry) = mac(0, *a, *b, 0);
                match carry {
                    0 => MaybeU64::U64(c),
                    _ => MaybeU64::Full(
                        F::from_repr(
                            [
                                c.to_le_bytes().as_ref(),
                                carry.to_le_bytes().as_ref(),
                                vec![0u8; 16].as_ref(),
                            ]
                            .concat()
                            .try_into()
                            .unwrap(),
                        )
                        .unwrap(),
                    ),
                }
            }
            (MaybeU64::U64(a), MaybeU64::Full(b)) => MaybeU64::Full(F::from(*a) * b),
            (MaybeU64::Full(a), MaybeU64::U64(b)) => MaybeU64::Full(*a * F::from(*b)),
            (MaybeU64::Full(a), MaybeU64::Full(b)) => MaybeU64::Full(*a * *b),
        }
    }
}

impl<'a, F> Mul<MaybeU64<F>> for &'a MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn mul(self, rhs: MaybeU64<F>) -> Self::Output {
        self * &rhs
    }
}

impl<'b, F> Mul<&'b MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn mul(self, rhs: &'b MaybeU64<F>) -> Self::Output {
        &self * rhs
    }
}

impl<F> Mul<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn mul(self, rhs: MaybeU64<F>) -> Self::Output {
        &self * &rhs
    }
}

impl<F> MulAssign<MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: MaybeU64<F>) {
        *self = &*self * &rhs;
    }
}

impl<'b, F> MulAssign<&'b MaybeU64<F>> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    #[inline]
    fn mul_assign(&mut self, rhs: &'b MaybeU64<F>) {
        *self = &*self * rhs;
    }
}

impl<F> Neg for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
{
    type Output = MaybeU64<F>;

    fn neg(self) -> <Self as Neg>::Output {
        match self {
            MaybeU64::U64(a) => MaybeU64::Full(-F::from(a)),
            MaybeU64::Full(a) => MaybeU64::Full(-a),
        }
    }
}

impl<F, T> Product<T> for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]>,
    T: core::borrow::Borrow<Self>,
{
    fn product<I: Iterator<Item = T>>(iter: I) -> Self {
        iter.fold(Self::ONE, |acc, item| acc * item.borrow())
    }
}


impl<F> PartialOrd for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]> + Ord,
{
    fn partial_cmp(&self, other: &MaybeU64<F>) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl<F> Ord for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]> + Ord,
{
    fn cmp(&self, other: &MaybeU64<F>) -> Ordering {
        match (self, other) {
            (MaybeU64::U64(a), MaybeU64::U64(b)) => a.cmp(b),
            (MaybeU64::U64(a), MaybeU64::Full(b)) => F::from(*a).cmp(b),
            (MaybeU64::Full(a), MaybeU64::U64(b)) => a.cmp(&F::from(*b)),
            (MaybeU64::Full(a), MaybeU64::Full(b)) => a.cmp(b),
        }
    }
}