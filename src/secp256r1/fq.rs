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
/// q = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
const MODULUS: Fq = Fq([
    0xf3b9cac2fc632551,
    0xbce6faada7179e84,
    0xffffffffffffffff,
    0xffffffff00000000,
]);

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 8] = [
    0xfc63_2551,
    0xf3b9_cac2,
    0xa717_9e84,
    0xbce6_faad,
    0xffff_ffff,
    0xffff_ffff,
    0x0000_0000,
    0xffff_ffff,
];

///Constant representing the modulus as static str
const MODULUS_STR: &str = "0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551";

/// INV = -(q^{-1} mod 2^64) mod 2^64
const INV: u64 = 0xccd1c8aaee00bc4f;

/// R = 2^256 mod q
/// 0xffffffff00000000000000004319055258e8617b0c46353d039cdaaf
const R: Fq = Fq([
    0x0c46353d039cdaaf,
    0x4319055258e8617b,
    0x0000000000000000,
    0xffffffff,
]);

/// R^2 = 2^512 mod q
/// 0x66e12d94f3d956202845b2392b6bec594699799c49bd6fa683244c95be79eea2
const R2: Fq = Fq([
    0x83244c95be79eea2,
    0x4699799c49bd6fa6,
    0x2845b2392b6bec59,
    0x66e12d94f3d95620,
]);

/// R^3 = 2^768 mod q
/// 0x503a54e76407be652543b9246ba5e93f111f28ae0c0555c9ac8ebec90b65a624
const R3: Fq = Fq([
    0xac8ebec90b65a624,
    0x111f28ae0c0555c9,
    0x2543b9246ba5e93f,
    0x503a54e76407be65,
]);

/// `GENERATOR = 7 mod r` is a generator of the `q - 1` order multiplicative
/// subgroup, or in other words a primitive root of the field.
/// It's derived with SageMath with: `GF(MODULUS).primitive_element()`.
const MULTIPLICATIVE_GENERATOR: Fq = Fq::from_raw([0x07, 0x00, 0x00, 0x00]);

/// Size of the 2-adic sub-group of the field.
const S: u32 = 4;

/// GENERATOR^t where t * 2^s + 1 = r with t odd. In other words, this is a 2^s root of unity.
/// `ffc97f062a770992ba807ace842a3dfc1546cad004378daf0592d7fbb41e6602`
const ROOT_OF_UNITY: Fq = Fq::from_raw([
    0x0592d7fbb41e6602,
    0x1546cad004378daf,
    0xba807ace842a3dfc,
    0xffc97f062a770992,
]);

/// 1 / ROOT_OF_UNITY mod q
/// `a0a66a5562d46f2ac645fa0458131caee3ac117c794c4137379c7f0657c73764`
const ROOT_OF_UNITY_INV: Fq = Fq::from_raw([
    0x379c7f0657c73764,
    0xe3ac117c794c4137,
    0xc645fa0458131cae,
    0xa0a66a5562d46f2a,
]);

/// 1 / 2 mod q
const TWO_INV: Fq = Fq::from_raw([
    0x79dce5617e3192a9,
    0xde737d56d38bcf42,
    0x7fffffffffffffff,
    0x7fffffff80000000,
]);

const ZETA: Fq = Fq::from_raw([
    0x7cbf87ff12884e21,
    0x9405335ce9c83e1d,
    0x4e786d0777fd6aef,
    0x52891d43d946a035,
]);

/// Generator of the t-order multiplicative subgroup.
/// Computed by exponentiating Self::MULTIPLICATIVE_GENERATOR by 2^s, where s is Self::S.
const DELTA: Fq = Fq::from_raw([0x1e39a5057d81, 0, 0, 0]);

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
        // 7fffffff800000007fffffffffffffffde737d56d38bcf4279dce5617e3192a
        let t = [
            0x279dce5617e3192a,
            0xfde737d56d38bcf4,
            0x07ffffffffffffff,
            0x7fffffff8000000,
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
    crate::field_testing_suite!(Fq, "serialization_check");
    crate::field_testing_suite!(Fq, "constants", MODULUS_STR);
    crate::field_testing_suite!(Fq, "sqrt");
    crate::field_testing_suite!(Fq, "zeta");
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 64);
}
