use crate::arithmetic::{adc, mac, sbb};
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::serde::SerdeObject;
use core::convert::TryInto;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::slice::Iter;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

// Number of 64 bit limbs to represent the field element
pub(crate) const NUM_BITS: u32 = 446;

// Inverter constant
const BYIL: usize = 9;

// Jabobi constant
const JACOBI_L: usize = 8;

/// Constant representing the modulus
/// p = 0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5cda8a6c7be4a7a5fe8fadffd6a2a7e8c30006b9459ffffcd300000001
const MODULUS: Fp = Fp([
    0x9ffffcd300000001,
    0xa2a7e8c30006b945,
    0xe4a7a5fe8fadffd6,
    0x443f9a5cda8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 14] = [
    0x00000001, 0x9ffffcd3, 0x0006b945, 0xa2a7e8c3, 0x8fadffd6, 0xe4a7a5fe, 0xda8a6c7b, 0x443f9a5c,
    0xf439266f, 0xa803ca76, 0x0d7f70e4, 0x0130e000, 0x00002400, 0x24000000,
];

// pub const NEGATIVE_ONE: Fp = Fp([]);

pub(crate) const MODULUS_STR: &str = "0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5cda8a6c7be4a7a5fe8fadffd6a2a7e8c30006b9459ffffcd300000001";

/// INV = -r^{-1} mod 2^64
/// `0x9ffffcd2ffffffff`
const INV: u64 = 0x9ffffcd2ffffffff;

/// Let M be the power of `2^64` nearest to `Self::MODULUS_BITS`. Then `R = M % Self::MODULUS`.
/// `R = 2^448 mod p`
/// `0x3ffffffffff03fff7a9dfffa183e9bf67e576bf526ff2f52242c7760637089cbf6a760a123e01218d68a2aaffd0ef18a000163afffffff9`
const R: Fp = Fp([
    0xa000163afffffff9,
    0x8d68a2aaffd0ef18,
    0xbf6a760a123e0121,
    0x2242c7760637089c,
    0x67e576bf526ff2f5,
    0xf7a9dfffa183e9bf,
    0x03ffffffffff03ff,
]);

/// `R^2 = 2^896 mod p`
/// `0x1a4b16581f66e3cc8bcb0f20758aec8520b6db3d7481a84c734fd363b575c23e7a42067a8ccd154b4b20c07277ae01f1d9702c6d54dc0598`
const R2: Fp = Fp([
    0xd9702c6d54dc0598,
    0x4b20c07277ae01f1,
    0x7a42067a8ccd154b,
    0x734fd363b575c23e,
    0x20b6db3d7481a84c,
    0x8bcb0f20758aec85,
    0x1a4b16581f66e3cc,
]);

/// `R^3 = 2^1792 mod p`
/// `0x1f51e40a048ddc1789010189f4df0ae1f3bc57efac4b3280b25aa8b46a40b225e5446680e4c4ea0449937d6b40e58f05c67afa3fe916dd69`
const R3: Fp = Fp([
    0xc67afa3fe916dd69,
    0x49937d6b40e58f05,
    0xe5446680e4c4ea04,
    0xb25aa8b46a40b225,
    0xf3bc57efac4b3280,
    0x89010189f4df0ae1,
    0x1f51e40a048ddc17,
]);

/// `GENERATOR = 10 mod p` is a generator of the `p - 1` order multiplicative
/// subgroup, or in other words a primitive root of the field.
const MULTIPLICATIVE_GENERATOR: Fp = Fp::from_raw([0x0a, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);

/// Size of the 2-adic sub-group of the field.
const S: u32 = 32;

/// GENERATOR^t where t * 2^s + 1 = p
/// with t odd. In other words, this
/// is a 2^s root of unity.
/// `0x2d39f8c5f9adb3f35fe3f4222db17451ddd9602a013af5276bdbe3903ec85fc889232f5c8bc6857060c75e6f399661d6c7b82d31d563091`
const ROOT_OF_UNITY: Fp = Fp::from_raw([
    0x6c7b82d31d563091,
    0x060c75e6f399661d,
    0x889232f5c8bc6857,
    0x76bdbe3903ec85fc,
    0x1ddd9602a013af52,
    0x35fe3f4222db1745,
    0x02d39f8c5f9adb3f,
]);

/// 1 / ROOT_OF_UNITY mod p
/// `0x17725d635b00cda4153eb10c7105919d012822bd86c08691803272fbc5c9f8378055eb56ae2d55f9272bf208aad57f666deaead2c693ff66`
const ROOT_OF_UNITY_INV: Fp = Fp::from_raw([
    0x6deaead2c693ff66,
    0x272bf208aad57f66,
    0x8055eb56ae2d55f9,
    0x803272fbc5c9f837,
    0x012822bd86c08691,
    0x153eb10c7105919d,
    0x17725d635b00cda4,
]);

/// 1 / 2 mod p
/// `0x12000000000012000098700006bfb8725401e53b7a1c9337a21fcd2e6d45363df253d2ff47d6ffeb5153f46180035ca2cffffe6980000001`
pub(crate) const TWO_INV: Fp = Fp::from_raw([
    0xcffffe6980000001,
    0x5153f46180035ca2,
    0xf253d2ff47d6ffeb,
    0xa21fcd2e6d45363d,
    0x5401e53b7a1c9337,
    0x0098700006bfb872,
    0x1200000000001200,
]);
/// GENERATOR^{2^s} where t * 2^s + 1 = r with t odd. In other words, this is a t root of unity.
/// `0xeacefc6504d028d42ed23fc8766d5a5f195b456887e1e0021fb760c53233e9170c23749b459b95cc6cbb5faf3754a1e1916b2007775db04`
const DELTA: Fp = Fp::from_raw([
    0x1916b2007775db04,
    0xc6cbb5faf3754a1e,
    0x70c23749b459b95c,
    0x21fb760c53233e91,
    0xf195b456887e1e00,
    0x42ed23fc8766d5a5,
    0x0eacefc6504d028d,
]);

/// `ZETA^3 = 1 mod p` where `ZETA^2 != 1 mod p`
/// `0x480000000000360001c950000d7ee0e4a803c956d01c903d720dc8ad8b38dffaf50c100004c37ffffffe`
const ZETA: Fp = Fp::from_raw([
    0x100004c37ffffffe,
    0xc8ad8b38dffaf50c,
    0xc956d01c903d720d,
    0x50000d7ee0e4a803,
    0x00000000360001c9,
    0x0000000000004800,
    0x0000000000000000,
]);

/// NEG_ONE ; -1 mod p
pub(crate) const NEG_ONE: Fp = Fp::from_raw([
    0x9ffffcd300000000,
    0xa2a7e8c30006b945,
    0xe4a7a5fe8fadffd6,
    0x443f9a5cda8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);

use crate::{
    extend_field_legendre, field_arithmetic_7, field_bits, impl_add_binop_specify_impl,
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_field, impl_from_u64,
    impl_from_uniform_bytes, impl_prime_field, impl_serde_object, impl_sub_binop_specify_output,
    impl_sum_prod, pow_vartime,
};
impl_binops_additive!(Fp, Fp);
impl_binops_multiplicative!(Fp, Fp);
impl_add_binop_specify_impl!(Fp);
impl_field!(Fp, sparse);
impl_serde_object!(Fp);
impl_prime_field!(Fp, ReprFp, le);
impl_sum_prod!(Fp);
extend_field_legendre!(Fp);
impl_from_u64!(Fp);
impl_from_uniform_bytes!(Fp, 64);
impl_from_uniform_bytes!(Fp, 72);
impl_from_uniform_bytes!(Fp, 112);

field_arithmetic_7!(Fp);

#[cfg(target_pointer_width = "64")]
field_bits!(Fp);
#[cfg(not(target_pointer_width = "64"))]
field_bits!(Fp);

#[cfg(feature = "derive_serde")]
crate::serialize_deserialize_primefield!(Fp, ReprFp);

impl Fp {
    fn sqrt(&self) -> CtOption<Self> {
        /// `(t - 1) // 2` where t * 2^s + 1 = p with t odd.
        const T_MINUS1_OVER2: [u64; 7] = [
            0x80035ca2cffffe69,
            0x47d6ffeb5153f461,
            0x6d45363df253d2ff,
            0x7a1c9337a21fcd2e,
            0x06bfb8725401e53b,
            0x0000120000987000,
            0x0000000012000000,
        ];
        ff::helpers::sqrt_tonelli_shanks(self, T_MINUS1_OVER2)
    }
}

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

/// Canonical little-endian representation of a `Fp` element.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
pub struct ReprFp {
    #[cfg_attr(feature = "derive_serde", serde(with = "serde_arrays"))]
    pub repr: [u8; Fp::SIZE],
}

impl ReprFp {
    /// Returns an iterator over the bytes of the canonical representation of the element.
    pub fn iter(&self) -> Iter<'_, u8> {
        self.repr.iter()
    }
}

impl Default for ReprFp {
    fn default() -> Self {
        ReprFp {
            repr: [0u8; Fp::SIZE],
        }
    }
}

impl AsRef<[u8]> for ReprFp {
    fn as_ref(&self) -> &[u8] {
        self.repr.as_ref()
    }
}

impl AsMut<[u8]> for ReprFp {
    fn as_mut(&mut self) -> &mut [u8] {
        self.repr.as_mut()
    }
}
impl From<[u8; Fp::SIZE]> for ReprFp {
    fn from(repr: [u8; Fp::SIZE]) -> Self {
        Self { repr }
    }
}

impl From<ReprFp> for [u8; Fp::SIZE] {
    fn from(repr: ReprFp) -> Self {
        repr.repr
    }
}

impl core::ops::Index<usize> for ReprFp {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.repr[index]
    }
}

impl core::ops::Index<core::ops::Range<usize>> for ReprFp {
    type Output = [u8];

    fn index(&self, index: core::ops::Range<usize>) -> &Self::Output {
        &self.repr[index]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fp, "field_arithmetic");
    crate::field_testing_suite!(Fp, "conversion");
    crate::field_testing_suite!(Fp, "serialization");
    crate::field_testing_suite!(Fp, "quadratic_residue");
    crate::field_testing_suite!(Fp, "bits");
    crate::field_testing_suite!(Fp, "serialization_check");
    crate::field_testing_suite!(Fp, "constants", MODULUS_STR);
    crate::field_testing_suite!(Fp, "sqrt");
    crate::field_testing_suite!(Fp, "zeta");
    crate::field_testing_suite!(Fp, "from_uniform_bytes", 64);
}
