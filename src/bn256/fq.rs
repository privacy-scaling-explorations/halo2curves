#[cfg(feature = "asm")]
use crate::bn256::assembly::field_arithmetic_asm;
#[cfg(not(feature = "asm"))]
use crate::{
    arithmetic::{bigint_geq, macx},
    field_arithmetic_4, field_specific_4,
};

use crate::arithmetic::{adc, mac, sbb};
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::serde::SerdeObject;
use core::convert::TryInto;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// Number of 64 bit limbs to represent the field element
pub(super) const NUM_BITS: u32 = 254;

// Inverter constant
const BYIL: usize = 6;

// Jabobi constant
const JACOBI_L: usize = 5;

/// Constant representing the modulus
/// q = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
const MODULUS: Fq = Fq([
    0x3c208c16d87cfd47,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
]);

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 8] = [
    0xd87c_fd47,
    0x3c20_8c16,
    0x6871_ca8d,
    0x9781_6a91,
    0x8181_585d,
    0xb850_45b6,
    0xe131_a029,
    0x3064_4e72,
];

/// INV = -(q^{-1} mod 2^64) mod 2^64
const INV: u64 = 0x87d20782e4866389;

/// R = 2^256 mod q
const R: Fq = Fq([
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f,
]);

/// R^2 = 2^512 mod q
const R2: Fq = Fq([
    0xf32cfc5b538afa89,
    0xb5e71911d44501fb,
    0x47ab1eff0a417ff6,
    0x06d89f71cab8351f,
]);

/// R^3 = 2^768 mod q
const R3: Fq = Fq([
    0xb1cd6dafda1530df,
    0x62f210e6a7283db6,
    0xef7f0b0c0ada0afb,
    0x20fd6e902d592544,
]);

pub const NEGATIVE_ONE: Fq = Fq([
    0x68c3488912edefaa,
    0x8d087f6872aabf4f,
    0x51e1a24709081231,
    0x2259d6b14729c0fa,
]);

pub(crate) const MODULUS_STR: &str =
    "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47";

/// Obtained with:
/// `sage: GF(0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47).primitive_element()`
const MULTIPLICATIVE_GENERATOR: Fq = Fq::from_raw([0x03, 0x0, 0x0, 0x0]);

const TWO_INV: Fq = Fq::from_raw([
    0x9e10460b6c3e7ea4,
    0xcbc0b548b438e546,
    0xdc2822db40c0ac2e,
    0x183227397098d014,
]);

/// `0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46`
const ROOT_OF_UNITY: Fq = Fq::from_raw([
    0x3c208c16d87cfd46,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
]);

/// `0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd46`
const ROOT_OF_UNITY_INV: Fq = Fq::from_raw([
    0x3c208c16d87cfd46,
    0x97816a916871ca8d,
    0xb85045b68181585d,
    0x30644e72e131a029,
]);

// `0x9`
const DELTA: Fq = Fq::from_raw([0x9, 0, 0, 0]);

/// `ZETA^3 = 1 mod r` where `ZETA^2 != 1 mod r`
const ZETA: Fq = Fq::from_raw([
    0xe4bd44e5607cfd48,
    0xc28f069fbb966e3d,
    0x5e6dd9e7e0acccb0,
    0x30644e72e131a029,
]);

/// Size of the 2-adic sub-group of the field.
const S: u32 = 0;

use crate::{
    const_montgomery_4, extend_field_legendre, field_bits, impl_add_binop_specify_impl,
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_field, impl_from_u64,
    impl_from_uniform_bytes, impl_prime_field, impl_serde_object, impl_sub_binop_specify_output,
    impl_sum_prod, pow_vartime,
};
impl_binops_additive!(Fq, Fq);
impl_binops_multiplicative!(Fq, Fq);
impl_add_binop_specify_impl!(Fq);
impl_field!(Fq, sparse);
impl_serde_object!(Fq);
impl_prime_field!(Fq, [u8; 32], le);
impl_sum_prod!(Fq);
extend_field_legendre!(Fq);
impl_from_uniform_bytes!(Fq, 64);
impl_from_uniform_bytes!(Fq, 48);
impl_from_u64!(Fq);

const_montgomery_4!(Fq);
#[cfg(not(feature = "asm"))]
field_arithmetic_4!(Fq, sparse);
#[cfg(feature = "asm")]
field_arithmetic_asm!(Fq, MODULUS, INV);

#[cfg(target_pointer_width = "64")]
field_bits!(Fq);
#[cfg(not(target_pointer_width = "64"))]
field_bits!(Fq);

#[cfg(feature = "derive_serde")]
crate::serialize_deserialize_primefield!(Fq, [u8; 32]);

impl Fq {
    /// Computes the square root of this element, if it exists.
    fn sqrt(&self) -> CtOption<Self> {
        let tmp = self.pow([
            0x4f082305b61f3f52,
            0x65e05aa45a1c72a3,
            0x6e14116da0605617,
            0x0c19139cb84c680a,
        ]);

        CtOption::new(tmp, tmp.square().ct_eq(self))
    }
}

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fq, "field_arithmetic");
    crate::field_testing_suite!(Fq, "conversion");
    crate::field_testing_suite!(Fq, "serialization");
    crate::field_testing_suite!(Fq, "quadratic_residue");
    crate::field_testing_suite!(Fq, "bits");
    crate::field_testing_suite!(Fq, "serialization_check");
    crate::field_testing_suite!(Fq, "constants", MODULUS_STR);
    crate::field_testing_suite!(Fq, "sqrt");
    crate::field_testing_suite!(Fq, "zeta");
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 64);
}
