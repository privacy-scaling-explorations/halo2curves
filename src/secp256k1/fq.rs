use crate::arithmetic::{adc, bigint_geq, mac, macx, sbb};
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::serde::SerdeObject;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// Number of 64 bit limbs to represent the field element
pub(crate) const NUM_BITS: u32 = 256;

// Inverter constant
const BYIL: usize = 6;

// Jabobi constant
const JACOBI_L: usize = 5;

/// Constant representing the modulus
/// q = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141
const MODULUS: Fq = Fq([
    0xbfd25e8cd0364141,
    0xbaaedce6af48a03b,
    0xfffffffffffffffe,
    0xffffffffffffffff,
]);

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 8] = [
    0xd036_4141,
    0xbfd2_5e8c,
    0xaf48_a03b,
    0xbaae_dce6,
    0xffff_fffe,
    0xffff_ffff,
    0xffff_ffff,
    0xffff_ffff,
];

///Constant representing the modulus as static str
const MODULUS_STR: &str = "0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141";

/// INV = -(q^{-1} mod 2^64) mod 2^64
const INV: u64 = 0x4b0dff665588b13f;

/// R = 2^256 mod q
/// 0x14551231950b75fc4402da1732fc9bebf
const R: Fq = Fq([0x402da1732fc9bebf, 0x4551231950b75fc4, 0x1, 0]);

/// R^2 = 2^512 mod q
/// 0x9d671cd581c69bc5e697f5e45bcd07c6741496c20e7cf878896cf21467d7d140
const R2: Fq = Fq([
    0x896cf21467d7d140,
    0x741496c20e7cf878,
    0xe697f5e45bcd07c6,
    0x9d671cd581c69bc5,
]);

/// R^3 = 2^768 mod q
/// 0x555d800c18ef116db1b31347f1d0b2da0017648444d4322c7bc0cfe0e9ff41ed
const R3: Fq = Fq([
    0x7bc0cfe0e9ff41ed,
    0x0017648444d4322c,
    0xb1b31347f1d0b2da,
    0x555d800c18ef116d,
]);

/// `GENERATOR = 7 mod r` is a generator of the `q - 1` order multiplicative
/// subgroup, or in other words a primitive root of the field.
/// It's derived with SageMath with: `GF(MODULUS).primitive_element()`.
const MULTIPLICATIVE_GENERATOR: Fq = Fq::from_raw([0x07, 0x00, 0x00, 0x00]);

/// Size of the 2-adic sub-group of the field.
const S: u32 = 6;

/// GENERATOR^t where t * 2^s + 1 = r with t odd. In other words, this is a 2^s root of unity.
/// `0xc1dc060e7a91986df9879a3fbc483a898bdeab680756045992f4b5402b052f2`
const ROOT_OF_UNITY: Fq = Fq::from_raw([
    0x992f4b5402b052f2,
    0x98bdeab680756045,
    0xdf9879a3fbc483a8,
    0xc1dc060e7a91986,
]);

/// 1 / ROOT_OF_UNITY mod q
const ROOT_OF_UNITY_INV: Fq = Fq::from_raw([
    0xb6fb30a0884f0d1c,
    0x77a275910aa413c3,
    0xefc7b0c75b8cbb72,
    0xfd3ae181f12d7096,
]);

/// 1 / 2 mod q
const TWO_INV: Fq = Fq::from_raw([
    0xdfe92f46681b20a1,
    0x5d576e7357a4501d,
    0xffffffffffffffff,
    0x7fffffffffffffff,
]);

/// Generator of the t-order multiplicative subgroup.
/// Computed by exponentiating Self::MULTIPLICATIVE_GENERATOR by 2^s, where s is Self::S.
/// `0x0000000000000000000cbc21fe4561c8d63b78e780e1341e199417c8c0bb7601`
const DELTA: Fq = Fq([
    0xd91b33d24319d9e8,
    0xb81c6596ff5d6740,
    0xa463969ca14c51c1,
    0x1900960de4b7929c,
]);

const ZETA: Fq = Fq::from_raw([
    0xdf02967c1b23bd72,
    0x122e22ea20816678,
    0xa5261c028812645a,
    0x5363ad4cc05c30e0,
]);

use crate::{
    const_montgomery_4, extend_field_legendre, field_arithmetic_4, field_bits, field_specific_4,
    impl_add_binop_specify_impl, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_field, impl_from_u64, impl_from_uniform_bytes,
    impl_prime_field, impl_serde_object, impl_sub_binop_specify_output, impl_sum_prod, pow_vartime,
};
impl_binops_additive!(Fq, Fq);
impl_binops_multiplicative!(Fq, Fq);
impl_add_binop_specify_impl!(Fq);
impl_field!(Fq, dense);
impl_serde_object!(Fq);
impl_prime_field!(Fq, [u8; 32], le);
impl_sum_prod!(Fq);
extend_field_legendre!(Fq);
impl_from_uniform_bytes!(Fq, 64);
impl_from_uniform_bytes!(Fq, 48);
impl_from_u64!(Fq);
field_bits!(Fq);

const_montgomery_4!(Fq);
field_arithmetic_4!(Fq, dense);

#[cfg(feature = "derive_serde")]
crate::serialize_deserialize_primefield!(Fq, [u8; 32]);

impl Fq {
    fn sqrt(&self) -> CtOption<Self> {
        let t = [
            0x777fa4bd19a06c82,
            0xfd755db9cd5e9140,
            0xffffffffffffffff,
            0x01ffffffffffffff,
        ];

        ff::helpers::sqrt_tonelli_shanks(self, t)
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
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 48, 64);
}
