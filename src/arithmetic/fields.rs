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
    /// Modulus of the field written as a string for display purposes
    const MODULUS: &'static str;

    /// This computes a random element of the field using system randomness.
    fn rand() -> Self {
        Self::random(rand::rngs::OsRng)
    }

    /// Returns whether or not this element is zero.
    fn ct_is_zero(&self) -> Choice {
        self.ct_eq(&Self::zero())
    }

    /// Obtains a field element that is congruent to the provided little endian
    /// byte representation of an integer.
    fn from_bytes_wide(bytes: &[u8; 64]) -> Self;

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
}

pub trait FieldExt: ff::PrimeField + BaseExt + Group<Scalar = Self> {
    /// Inverse of $2$ in the field.
    const TWO_INV: Self;

    /// Inverse of `ROOT_OF_UNITY`
    const ROOT_OF_UNITY_INV: Self;

    /// Generator of the $t-order$ multiplicative subgroup
    const DELTA: Self;

    /// Element of multiplicative order $3$.
    const ZETA: Self;

    /// Obtains a field element congruent to the integer `v`.
    fn from_u128(v: u128) -> Self;

    // /// Converts this field element to its normalized, little endian byte
    // /// representation.
    // fn to_bytes(&self) -> [u8; 32];

    // /// Attempts to obtain a field element from its normalized, little endian
    // /// byte representation.
    // fn from_bytes(bytes: &[u8; 32]) -> CtOption<Self>;

    /// Gets the lower 128 bits of this field element when expressed
    /// canonically.
    fn get_lower_128(&self) -> u128;
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

/// Compute a + (b * c), returning the result and the new carry over.
#[inline(always)]
pub(crate) const fn macx(a: u64, b: u64, c: u64) -> (u64, u64) {
    let res = (a as u128) + ((b as u128) * (c as u128));
    (res as u64, (res >> 64) as u64)
}

/// Compute a * b, returning the result.
#[inline(always)]
pub(crate) fn mul_512(a: [u64; 4], b: [u64; 4]) -> [u64; 8] {
    let (r0, carry) = macx(0, a[0], b[0]);
    let (r1, carry) = macx(carry, a[0], b[1]);
    let (r2, carry) = macx(carry, a[0], b[2]);
    let (r3, carry_out) = macx(carry, a[0], b[3]);

    let (r1, carry) = macx(r1, a[1], b[0]);
    let (r2, carry) = mac(r2, a[1], b[1], carry);
    let (r3, carry) = mac(r3, a[1], b[2], carry);
    let (r4, carry_out) = mac(carry_out, a[1], b[3], carry);

    let (r2, carry) = macx(r2, a[2], b[0]);
    let (r3, carry) = mac(r3, a[2], b[1], carry);
    let (r4, carry) = mac(r4, a[2], b[2], carry);
    let (r5, carry_out) = mac(carry_out, a[2], b[3], carry);

    let (r3, carry) = macx(r3, a[3], b[0]);
    let (r4, carry) = mac(r4, a[3], b[1], carry);
    let (r5, carry) = mac(r5, a[3], b[2], carry);
    let (r6, carry_out) = mac(carry_out, a[3], b[3], carry);

    [r0, r1, r2, r3, r4, r5, r6, carry_out]
}
