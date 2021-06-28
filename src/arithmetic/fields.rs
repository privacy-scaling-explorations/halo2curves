//! This module contains the `Field` abstraction that allows us to write
//! code that generalizes over a pair of fields.

use core::mem::size_of;
use static_assertions::const_assert;
use subtle::{Choice, ConstantTimeEq};

use super::Group;

use std::io::{self, Read, Write};

const_assert!(size_of::<usize>() >= 4);

/// This trait is a common interface for dealing with elements of a finite
/// field.
pub trait BaseExt: ff::Field + Ord + ConstantTimeEq {
    /// This computes a random element of the field using system randomness.
    fn rand() -> Self {
        Self::random(rand::rngs::OsRng)
    }

    /// Returns whether or not this element is zero.
    fn ct_is_zero(&self) -> Choice {
        self.ct_eq(&Self::zero())
    }

    /// Writes this element in its normalized, little endian form into a buffer.
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()>;

    /// Reads a normalized, little endian represented field element from a
    /// buffer.
    fn read<R: Read>(reader: &mut R) -> io::Result<Self>;

    /// Exponentiates `self` by `by`, where `by` is a little-endian order
    /// integer exponent.
    fn pow(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 0x1) as u8).into());
            }
        }
        res
    }

    /// Performs a batch inversion using Montgomery's trick, returns the product
    /// of every inverse. Zero inputs are ignored.
    fn batch_invert(v: &mut [Self]) -> Self {
        let mut tmp = Vec::with_capacity(v.len());

        let mut acc = Self::one();
        for p in v.iter() {
            tmp.push(acc);
            acc = Self::conditional_select(&(acc * p), &acc, p.ct_is_zero());
        }

        acc = acc.invert().unwrap();
        let allinv = acc;

        for (p, tmp) in v.iter_mut().rev().zip(tmp.into_iter().rev()) {
            let skip = p.ct_is_zero();

            let tmp = tmp * acc;
            acc = Self::conditional_select(&(acc * *p), &acc, skip);
            *p = Self::conditional_select(&tmp, p, skip);
        }

        allinv
    }
}

pub trait FieldExt: ff::PrimeField + BaseExt + Group<Scalar = Self> {
    /// Modulus of the field written as a string for display purposes
    const MODULUS: &'static str;

    /// Generator of the $2^S$ multiplicative subgroup
    const ROOT_OF_UNITY: Self;

    /// Inverse of $2$ in the field.
    const TWO_INV: Self;

    /// Inverse of `ROOT_OF_UNITY`
    const ROOT_OF_UNITY_INV: Self;

    /// The value $(T-1)/2$ such that $2^S \cdot T = p - 1$ with $T$ odd.
    const T_MINUS1_OVER2: [u64; 4];

    /// Generator of the $t-order$ multiplicative subgroup
    const DELTA: Self;

    /// Ideally the smallest prime $\alpha$ such that gcd($p - 1$, $\alpha$) = $1$
    const RESCUE_ALPHA: u64;

    /// $RESCUE_INVALPHA \cdot RESCUE_ALPHA = 1 \mod p - 1$ such that
    /// `(a^RESCUE_ALPHA)^RESCUE_INVALPHA = a`.
    const RESCUE_INVALPHA: [u64; 4];

    /// Element of multiplicative order $3$.
    const ZETA: Self;

    /// Obtains a field element that is congruent to the provided little endian
    /// byte representation of an integer.
    fn from_bytes_wide(bytes: &[u8; 64]) -> Self;
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
