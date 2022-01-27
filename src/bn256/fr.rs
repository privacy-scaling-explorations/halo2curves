use super::LegendreSymbol;
use core::convert::TryInto;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::io::{self, Read, Write};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::arithmetic::{adc, mac, sbb, BaseExt, FieldExt, Group};

#[cfg(all(feature = "asm", target_arch = "x86_64"))]
use super::assembly::assembly_field;

#[derive(Clone, Copy, Eq, Hash)]
pub struct Fr(pub(crate) [u64; 4]);

/// Constant representing the modulus
/// q = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
pub const MODULUS: Fr = Fr([
    0x43e1f593f0000001,
    0x2833e84879b97091,
    0xb85045b68181585d,
    0x30644e72e131a029,
]);

/// INV = -(q^{-1} mod 2^64) mod 2^64
const INV: u64 = 0xc2e1f593efffffff;

/// R = 2^256 mod q
/// 0xe0a77c19a07df2f666ea36f7879462e36fc76959f60cd29ac96341c4ffffffb
const R: Fr = Fr([
    0xac96341c4ffffffb,
    0x36fc76959f60cd29,
    0x666ea36f7879462e,
    0x0e0a77c19a07df2f,
]);

/// R^2 = 2^512 mod q
/// 0x216d0b17f4e44a58c49833d53bb808553fe3ab1e35c59e31bb8e645ae216da7
const R2: Fr = Fr([
    0x1bb8e645ae216da7,
    0x53fe3ab1e35c59e3,
    0x8c49833d53bb8085,
    0x0216d0b17f4e44a5,
]);

/// R^3 = 2^768 mod q
/// 0xcf8594b7fcc657c893cc664a19fcfed2a489cbe1cfbb6b85e94d8e1b4bf0040
const R3: Fr = Fr([
    0x5e94d8e1b4bf0040,
    0x2a489cbe1cfbb6b8,
    0x893cc664a19fcfed,
    0x0cf8594b7fcc657c,
]);

const GENERATOR: Fr = Fr::from_raw([0x07, 0x00, 0x00, 0x00]);

const S: u32 = 28;

const ROOT_OF_UNITY: Fr = Fr::from_raw([
    0xd34f1ed960c37c9c,
    0x3215cf6dd39329c8,
    0x98865ea93dd31f74,
    0x03ddb9f5166d18b7,
]);

impl Group for Fr {
    type Scalar = Fr;

    fn group_zero() -> Self {
        Self::zero()
    }
    fn group_add(&mut self, rhs: &Self) {
        *self += *rhs;
    }
    fn group_sub(&mut self, rhs: &Self) {
        *self -= *rhs;
    }
    fn group_scale(&mut self, by: &Self::Scalar) {
        *self *= *by;
    }
}

impl ::std::fmt::Display for Fr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl fmt::Debug for Fr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let tmp = self.to_bytes();
        write!(f, "0x")?;
        for &b in tmp.iter().rev() {
            write!(f, "{:02x}", b)?;
        }
        Ok(())
    }
}

impl Default for Fr {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl From<bool> for Fr {
    fn from(bit: bool) -> Fr {
        if bit {
            Fr::one()
        } else {
            Fr::zero()
        }
    }
}

impl From<u64> for Fr {
    fn from(val: u64) -> Fr {
        Fr([val, 0, 0, 0]) * R2
    }
}

impl ConstantTimeEq for Fr {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.0[0].ct_eq(&other.0[0])
            & self.0[1].ct_eq(&other.0[1])
            & self.0[2].ct_eq(&other.0[2])
            & self.0[3].ct_eq(&other.0[3])
    }
}

impl PartialEq for Fr {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).unwrap_u8() == 1
    }
}

impl std::cmp::Ord for Fr {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let left = self.to_bytes();
        let right = other.to_bytes();
        left.iter()
            .zip(right.iter())
            .rev()
            .find_map(|(left_byte, right_byte)| match left_byte.cmp(right_byte) {
                std::cmp::Ordering::Equal => None,
                res => Some(res),
            })
            .unwrap_or(std::cmp::Ordering::Equal)
    }
}

impl std::cmp::PartialOrd for Fr {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl ConditionallySelectable for Fr {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fr([
            u64::conditional_select(&a.0[0], &b.0[0], choice),
            u64::conditional_select(&a.0[1], &b.0[1], choice),
            u64::conditional_select(&a.0[2], &b.0[2], choice),
            u64::conditional_select(&a.0[3], &b.0[3], choice),
        ])
    }
}

impl Neg for Fr {
    type Output = Fr;

    #[inline]
    fn neg(self) -> Fr {
        -&self
    }
}

impl<'a> Neg for &'a Fr {
    type Output = Fr;

    #[inline]
    fn neg(self) -> Fr {
        self.neg()
    }
}

impl<'a, 'b> Sub<&'b Fr> for &'a Fr {
    type Output = Fr;

    #[inline]
    fn sub(self, rhs: &'b Fr) -> Fr {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fr> for &'a Fr {
    type Output = Fr;

    #[inline]
    fn add(self, rhs: &'b Fr) -> Fr {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fr> for &'a Fr {
    type Output = Fr;

    #[inline]
    fn mul(self, rhs: &'b Fr) -> Fr {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fr, Fr);
impl_binops_multiplicative!(Fr, Fr);

impl Fr {
    pub fn legendre(&self) -> LegendreSymbol {
        unimplemented!()
    }

    /// Returns zero, the additive identity.
    #[inline]
    pub const fn zero() -> Fr {
        Fr([0, 0, 0, 0])
    }

    /// Returns one, the multiplicative identity.
    #[inline]
    pub const fn one() -> Fr {
        R
    }

    fn from_u512(limbs: [u64; 8]) -> Fr {
        // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
        // with the higher bits multiplied by 2^256. Thus, we perform two reductions
        //
        // 1. the lower bits are multiplied by R^2, as normal
        // 2. the upper bits are multiplied by R^2 * 2^256 = R^3
        //
        // and computing their sum in the field. It remains to see that arbitrary 256-bit
        // numbers can be placed into Montgomery form safely using the reduction. The
        // reduction works so long as the product is less than R=2^256 multiplied by
        // the modulus. This holds because for any `c` smaller than the modulus, we have
        // that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
        // reduction always works so long as `c` is in the field; in this case it is either the
        // constant `R2` or `R3`.
        let d0 = Fr([limbs[0], limbs[1], limbs[2], limbs[3]]);
        let d1 = Fr([limbs[4], limbs[5], limbs[6], limbs[7]]);
        // Convert to Montgomery form
        d0 * R2 + d1 * R3
    }

    /// Converts from an integer represented in little endian
    /// into its (congruent) `Fr` representation.
    pub const fn from_raw(val: [u64; 4]) -> Self {
        let (r0, carry) = mac(0, val[0], R2.0[0], 0);
        let (r1, carry) = mac(0, val[0], R2.0[1], carry);
        let (r2, carry) = mac(0, val[0], R2.0[2], carry);
        let (r3, r4) = mac(0, val[0], R2.0[3], carry);

        let (r1, carry) = mac(r1, val[1], R2.0[0], 0);
        let (r2, carry) = mac(r2, val[1], R2.0[1], carry);
        let (r3, carry) = mac(r3, val[1], R2.0[2], carry);
        let (r4, r5) = mac(r4, val[1], R2.0[3], carry);

        let (r2, carry) = mac(r2, val[2], R2.0[0], 0);
        let (r3, carry) = mac(r3, val[2], R2.0[1], carry);
        let (r4, carry) = mac(r4, val[2], R2.0[2], carry);
        let (r5, r6) = mac(r5, val[2], R2.0[3], carry);

        let (r3, carry) = mac(r3, val[3], R2.0[0], 0);
        let (r4, carry) = mac(r4, val[3], R2.0[1], carry);
        let (r5, carry) = mac(r5, val[3], R2.0[2], carry);
        let (r6, r7) = mac(r6, val[3], R2.0[3], carry);

        let k = r0.wrapping_mul(INV);
        let (_, carry) = mac(r0, k, MODULUS.0[0], 0);
        let (r1, carry) = mac(r1, k, MODULUS.0[1], carry);
        let (r2, carry) = mac(r2, k, MODULUS.0[2], carry);
        let (r3, carry) = mac(r3, k, MODULUS.0[3], carry);
        let (r4, carry2) = adc(r4, 0, carry);

        let k = r1.wrapping_mul(INV);
        let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
        let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
        let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
        let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
        let (r5, carry2) = adc(r5, carry2, carry);

        let k = r2.wrapping_mul(INV);
        let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
        let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
        let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
        let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
        let (r6, carry2) = adc(r6, carry2, carry);

        let k = r3.wrapping_mul(INV);
        let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
        let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
        let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
        let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
        let (r7, _) = adc(r7, carry2, carry);

        let (d0, borrow) = sbb(r4, MODULUS.0[0], 0);
        let (d1, borrow) = sbb(r5, MODULUS.0[1], borrow);
        let (d2, borrow) = sbb(r6, MODULUS.0[2], borrow);
        let (d3, borrow) = sbb(r7, MODULUS.0[3], borrow);

        // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
        // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
        let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
        let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
        let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
        let (d3, _) = adc(d3, MODULUS.0[3] & borrow, carry);

        Fr([d0, d1, d2, d3])
    }
}

#[cfg(all(feature = "asm", target_arch = "x86_64"))]
assembly_field!(Fr, MODULUS, INV);

#[cfg(any(not(feature = "asm"), not(target_arch = "x86_64")))]
impl Fr {
    /// Doubles this field element.
    #[inline]
    pub const fn double(&self) -> Fr {
        // TODO: This can be achieved more efficiently with a bitshift.
        self.add(self)
    }

    /// Squares this element.
    #[inline]
    pub const fn square(&self) -> Fr {
        let (r1, carry) = mac(0, self.0[0], self.0[1], 0);
        let (r2, carry) = mac(0, self.0[0], self.0[2], carry);
        let (r3, r4) = mac(0, self.0[0], self.0[3], carry);

        let (r3, carry) = mac(r3, self.0[1], self.0[2], 0);
        let (r4, r5) = mac(r4, self.0[1], self.0[3], carry);

        let (r5, r6) = mac(r5, self.0[2], self.0[3], 0);

        let r7 = r6 >> 63;
        let r6 = (r6 << 1) | (r5 >> 63);
        let r5 = (r5 << 1) | (r4 >> 63);
        let r4 = (r4 << 1) | (r3 >> 63);
        let r3 = (r3 << 1) | (r2 >> 63);
        let r2 = (r2 << 1) | (r1 >> 63);
        let r1 = r1 << 1;

        let (r0, carry) = mac(0, self.0[0], self.0[0], 0);
        let (r1, carry) = adc(0, r1, carry);
        let (r2, carry) = mac(r2, self.0[1], self.0[1], carry);
        let (r3, carry) = adc(0, r3, carry);
        let (r4, carry) = mac(r4, self.0[2], self.0[2], carry);
        let (r5, carry) = adc(0, r5, carry);
        let (r6, carry) = mac(r6, self.0[3], self.0[3], carry);
        let (r7, _) = adc(0, r7, carry);

        Fr::montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7)
    }

    #[allow(clippy::too_many_arguments)]
    #[inline(always)]
    const fn montgomery_reduce(
        r0: u64,
        r1: u64,
        r2: u64,
        r3: u64,
        r4: u64,
        r5: u64,
        r6: u64,
        r7: u64,
    ) -> Self {
        // The Montgomery reduction here is based on Algorithm 14.32 in
        // Handbook of Applied Cryptography
        // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

        let k = r0.wrapping_mul(INV);
        let (_, carry) = mac(r0, k, MODULUS.0[0], 0);
        let (r1, carry) = mac(r1, k, MODULUS.0[1], carry);
        let (r2, carry) = mac(r2, k, MODULUS.0[2], carry);
        let (r3, carry) = mac(r3, k, MODULUS.0[3], carry);
        let (r4, carry2) = adc(r4, 0, carry);

        let k = r1.wrapping_mul(INV);
        let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
        let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
        let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
        let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
        let (r5, carry2) = adc(r5, carry2, carry);

        let k = r2.wrapping_mul(INV);
        let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
        let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
        let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
        let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
        let (r6, carry2) = adc(r6, carry2, carry);

        let k = r3.wrapping_mul(INV);
        let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
        let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
        let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
        let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
        let (r7, _) = adc(r7, carry2, carry);

        // Result may be within MODULUS of the correct value
        (&Fr([r4, r5, r6, r7])).sub(&MODULUS)
    }

    /// Multiplies `rhs` by `self`, returning the result.
    #[inline]
    pub const fn mul(&self, rhs: &Self) -> Self {
        // Schoolbook multiplication

        let (r0, carry) = mac(0, self.0[0], rhs.0[0], 0);
        let (r1, carry) = mac(0, self.0[0], rhs.0[1], carry);
        let (r2, carry) = mac(0, self.0[0], rhs.0[2], carry);
        let (r3, r4) = mac(0, self.0[0], rhs.0[3], carry);

        let (r1, carry) = mac(r1, self.0[1], rhs.0[0], 0);
        let (r2, carry) = mac(r2, self.0[1], rhs.0[1], carry);
        let (r3, carry) = mac(r3, self.0[1], rhs.0[2], carry);
        let (r4, r5) = mac(r4, self.0[1], rhs.0[3], carry);

        let (r2, carry) = mac(r2, self.0[2], rhs.0[0], 0);
        let (r3, carry) = mac(r3, self.0[2], rhs.0[1], carry);
        let (r4, carry) = mac(r4, self.0[2], rhs.0[2], carry);
        let (r5, r6) = mac(r5, self.0[2], rhs.0[3], carry);

        let (r3, carry) = mac(r3, self.0[3], rhs.0[0], 0);
        let (r4, carry) = mac(r4, self.0[3], rhs.0[1], carry);
        let (r5, carry) = mac(r5, self.0[3], rhs.0[2], carry);
        let (r6, r7) = mac(r6, self.0[3], rhs.0[3], carry);

        Fr::montgomery_reduce(r0, r1, r2, r3, r4, r5, r6, r7)
    }

    /// Subtracts `rhs` from `self`, returning the result.
    #[inline]
    pub const fn sub(&self, rhs: &Self) -> Self {
        let (d0, borrow) = sbb(self.0[0], rhs.0[0], 0);
        let (d1, borrow) = sbb(self.0[1], rhs.0[1], borrow);
        let (d2, borrow) = sbb(self.0[2], rhs.0[2], borrow);
        let (d3, borrow) = sbb(self.0[3], rhs.0[3], borrow);

        // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
        // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
        let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
        let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
        let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
        let (d3, _) = adc(d3, MODULUS.0[3] & borrow, carry);

        Fr([d0, d1, d2, d3])
    }

    /// Adds `rhs` to `self`, returning the result.
    #[inline]
    pub const fn add(&self, rhs: &Self) -> Self {
        let (d0, carry) = adc(self.0[0], rhs.0[0], 0);
        let (d1, carry) = adc(self.0[1], rhs.0[1], carry);
        let (d2, carry) = adc(self.0[2], rhs.0[2], carry);
        let (d3, _) = adc(self.0[3], rhs.0[3], carry);

        // Attempt to subtract the modulus, to ensure the value
        // is smaller than the modulus.
        (&Fr([d0, d1, d2, d3])).sub(&MODULUS)
    }

    /// Negates `self`.
    #[inline]
    pub const fn neg(&self) -> Self {
        // Subtract `self` from `MODULUS` to negate. Ignore the final
        // borrow because it cannot underflow; self is guaranteed to
        // be in the field.
        let (d0, borrow) = sbb(MODULUS.0[0], self.0[0], 0);
        let (d1, borrow) = sbb(MODULUS.0[1], self.0[1], borrow);
        let (d2, borrow) = sbb(MODULUS.0[2], self.0[2], borrow);
        let (d3, _) = sbb(MODULUS.0[3], self.0[3], borrow);

        // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
        // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
        let mask = (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

        Fr([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
    }
}

impl From<Fr> for [u8; 32] {
    fn from(value: Fr) -> [u8; 32] {
        value.to_bytes()
    }
}

impl<'a> From<&'a Fr> for [u8; 32] {
    fn from(value: &'a Fr) -> [u8; 32] {
        value.to_bytes()
    }
}

impl ff::Field for Fr {
    fn random(mut rng: impl RngCore) -> Self {
        Self::from_u512([
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
            rng.next_u64(),
        ])
    }

    fn zero() -> Self {
        Self::zero()
    }

    fn one() -> Self {
        Self::one()
    }

    fn is_zero(&self) -> Choice {
        self.ct_is_zero()
    }

    fn double(&self) -> Self {
        self.double()
    }

    #[inline(always)]
    fn square(&self) -> Self {
        self.square()
    }

    /// Computes the square root of this element, if it exists.
    fn sqrt(&self) -> CtOption<Self> {
        unimplemented!()
    }

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    fn invert(&self) -> CtOption<Self> {
        let tmp = self.pow(&[
            0x43e1f593efffffff,
            0x2833e84879b97091,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ]);

        CtOption::new(tmp, !self.ct_eq(&Self::zero()))
    }
}

impl ff::PrimeField for Fr {
    type Repr = [u8; 32];

    const NUM_BITS: u32 = 254;
    const CAPACITY: u32 = 253;
    const S: u32 = S;

    fn from_repr(repr: Self::Repr) -> CtOption<Self> {
        let mut tmp = Fr([0, 0, 0, 0]);

        tmp.0[0] = u64::from_le_bytes(repr[0..8].try_into().unwrap());
        tmp.0[1] = u64::from_le_bytes(repr[8..16].try_into().unwrap());
        tmp.0[2] = u64::from_le_bytes(repr[16..24].try_into().unwrap());
        tmp.0[3] = u64::from_le_bytes(repr[24..32].try_into().unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], MODULUS.0[0], 0);
        let (_, borrow) = sbb(tmp.0[1], MODULUS.0[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], MODULUS.0[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], MODULUS.0[3], borrow);

        // If the element is smaller than MODULUS then the
        // subtraction will underflow, producing a borrow value
        // of 0xffff...ffff. Otherwise, it'll be zero.
        let is_some = (borrow as u8) & 1;

        // Convert to Montgomery form by computing
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    fn to_repr(&self) -> Self::Repr {
        // Turn into canonical form by computing
        // (a.R) / R = a
        #[cfg(all(feature = "asm", target_arch = "x86_64"))]
        let tmp = Fr::montgomery_reduce(&[self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0]);

        #[cfg(any(not(feature = "asm"), not(target_arch = "x86_64")))]
        let tmp = Fr::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);

        let mut res = [0; 32];
        res[0..8].copy_from_slice(&tmp.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp.0[3].to_le_bytes());

        res
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr()[0] & 1)
    }

    fn multiplicative_generator() -> Self {
        GENERATOR
    }

    fn root_of_unity() -> Self {
        ROOT_OF_UNITY
    }
}

impl BaseExt for Fr {
    const MODULUS: &'static str =
        "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001";

    /// Converts a 512-bit little endian integer into
    /// a `Fr` by reducing by the modulus.
    fn from_bytes_wide(bytes: &[u8; 64]) -> Self {
        Self::from_u512([
            u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
            u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
            u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
            u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
            u64::from_le_bytes(bytes[32..40].try_into().unwrap()),
            u64::from_le_bytes(bytes[40..48].try_into().unwrap()),
            u64::from_le_bytes(bytes[48..56].try_into().unwrap()),
            u64::from_le_bytes(bytes[56..64].try_into().unwrap()),
        ])
    }

    fn ct_is_zero(&self) -> Choice {
        self.ct_eq(&Self::zero())
    }

    /// Writes this element in its normalized, little endian form into a buffer.
    fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        let compressed = self.to_bytes();
        writer.write_all(&compressed[..])
    }

    /// Reads a normalized, little endian represented field element from a
    /// buffer.
    fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut compressed = [0u8; 32];
        reader.read_exact(&mut compressed[..])?;
        Option::from(Self::from_bytes(&compressed))
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "invalid point encoding in proof"))
    }
}

impl FieldExt for Fr {
    // const ROOT_OF_UNITY: Self = ROOT_OF_UNITY;

    const TWO_INV: Self = Fr::from_raw([
        0xa1f0fac9f8000001,
        0x9419f4243cdcb848,
        0xdc2822db40c0ac2e,
        0x183227397098d014,
    ]);

    const ROOT_OF_UNITY_INV: Self = Self::from_raw([
        0x0ed3e50a414e6dba,
        0xb22625f59115aba7,
        0x1bbe587180f34361,
        0x048127174daabc26,
    ]);

    // 0x09226b6e22c6f0ca64ec26aad4c86e715b5f898e5e963f25870e56bbe533e9a2
    const DELTA: Self = Self::from_raw([
        0x870e56bbe533e9a2,
        0x5b5f898e5e963f25,
        0x64ec26aad4c86e71,
        0x09226b6e22c6f0ca,
    ]);

    const RESCUE_ALPHA: u64 = 5;

    const RESCUE_INVALPHA: [u64; 4] = [
        0xcfe7f7a98ccccccd,
        0x535cb9d394945a0d,
        0x93736af8679aad17,
        0x26b6a528b427b354,
    ];

    const T_MINUS1_OVER2: [u64; 4] = [
        0xcdcb848a1f0fac9f,
        0x0c0ac2e9419f4243,
        0x098d014dc2822db4,
        0x0000000183227397,
    ];

    const ZETA: Self = Self::from_raw([
        0xb8ca0b2d36636f23,
        0xcc37a73fec2bc5e9,
        0x048b6e193fd84104,
        0x30644e72e131a029,
    ]);

    fn from_u128(v: u128) -> Self {
        Fr::from_raw([v as u64, (v >> 64) as u64, 0, 0])
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fr`, failing if the input is not canonical.
    fn from_bytes(bytes: &[u8; 32]) -> CtOption<Fr> {
        <Self as ff::PrimeField>::from_repr(*bytes)
    }

    /// Converts an element of `Fr` into a byte representation in
    /// little-endian byte order.
    fn to_bytes(&self) -> [u8; 32] {
        <Self as ff::PrimeField>::to_repr(self)
    }

    /// Gets the lower 128 bits of this field element when expressed
    /// canonically.
    fn get_lower_128(&self) -> u128 {
        #[cfg(all(feature = "asm", target_arch = "x86_64"))]
        let tmp = Fr::montgomery_reduce(&[self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0]);

        #[cfg(any(not(feature = "asm"), not(target_arch = "x86_64")))]
        let tmp = Fr::montgomery_reduce(self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0);

        u128::from(tmp.0[0]) | (u128::from(tmp.0[1]) << 64)
    }
}

#[cfg(test)]
use ff::{Field, PrimeField};

#[test]
fn test_zeta() {
    let a = Fr::ZETA;
    assert!(a != Fr::one());
    let b = a * a;
    assert!(b != Fr::one());
    let c = b * a;
    println!("{:?}", c);
    assert!(c == Fr::one());
}

#[test]
fn test_root_of_unity() {
    assert_eq!(
        Fr::root_of_unity().pow_vartime(&[1 << Fr::S, 0, 0, 0]),
        Fr::one()
    );
}

#[test]
fn test_inv_root_of_unity() {
    assert_eq!(Fr::ROOT_OF_UNITY_INV, Fr::root_of_unity().invert().unwrap());
}

#[test]
fn test_inv() {
    // Compute -(r^{-1} mod 2^64) mod 2^64 by exponentiating
    // by totient(2**64) - 1

    let mut inv = 1u64;
    for _ in 0..63 {
        inv = inv.wrapping_mul(inv);
        inv = inv.wrapping_mul(MODULUS.0[0]);
    }
    inv = inv.wrapping_neg();

    assert_eq!(inv, INV);
}

#[test]
fn test_inv_2() {
    assert_eq!(Fr::TWO_INV, Fr::from(2).invert().unwrap());
}

#[test]
fn test_from_u512() {
    assert_eq!(
        Fr::from_raw([
            0x7e7140b5196b9e6f,
            0x9abac9e4157b6172,
            0xf04bc41062fd7322,
            0x1185fa9c9fef6326,
        ]),
        Fr::from_u512([
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa,
            0xaaaaaaaaaaaaaaaa
        ])
    );
}

#[test]
fn test_field() {
    crate::tests::field::random_field_tests::<Fr>("fr".to_string());
}
