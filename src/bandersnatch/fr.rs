use crate::arithmetic::{adc, mac, macx, sbb};
use crate::extend_field_legendre;
use crate::ff::{FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// This represents an element of $\mathbb{F}_r$ where
///
/// `r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1`
///
/// is the scalar field of the bandersnatch curve.
// The internal representation of this type is four 64-bit unsigned
// integers in little-endian order. `Fr` values are always in
// Montgomery form; i.e., Fr(a) = aR mod q, with R = 2^256.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Fr(pub(crate) [u64; 4]);

#[cfg(feature = "derive_serde")]
crate::serialize_deserialize_32_byte_primefield!(Fr);

/// Constant representing the modulus
/// r = 0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1
const MODULUS: Fr = Fr([
    0x74fd06b52876e7e1,
    0xff8f870074190471,
    0x0cce760202687600,
    0x1cfb69d4ca675f52,
]);

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 8] = [
    0x2876_e7e1,
    0x74fd_06b5,
    0x7419_0471,
    0xff8f_8700,
    0x0268_7600,
    0x0cce_7602,
    0xca67_5f52,
    0x1cfb_69d4,
];

///Constant representing the modulus as static str
const MODULUS_STR: &str = "0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1";

/// INV = -(r^{-1} mod 2^64) mod 2^64
const INV: u64 = 0xF19F22295CC063DF;

/// R = 2^256 mod r
/// 0x1824B159ACC5056F998C4FEFECBC4FF80383C7FC5F37DC745817CA56BC48C0F8
const R: Fr = Fr([
    0x5817CA56BC48C0F8,
    0x0383C7FC5F37DC74,
    0x998C4FEFECBC4FF8,
    0x1824B159ACC5056F,
]);

/// R^2 = 2^512 mod r
/// 0xAE793DDB14AEC7DAA9E6DAEC0055CEA40FA7CA27FECB938DBB4F5D658DB47CB
const R2: Fr = Fr([
    0xDBB4F5D658DB47CB,
    0x40FA7CA27FECB938,
    0xAA9E6DAEC0055CEA,
    0xAE793DDB14AEC7D,
]);

/// R^3 = 2^768 mod r
/// 0x53A3B49F57751131126587105341936B0CE06DAEDDD77691EB6E3EB79377BD1
const R3: Fr = Fr([
    0x1EB6E3EB79377BD1,
    0xB0CE06DAEDDD7769,
    0x1126587105341936,
    0x53A3B49F5775113,
]);

/// `GENERATOR = 7 mod r` is a generator of the `r - 1` order multiplicative
/// subgroup, or in other words a primitive root of the field.
/// It's derived with SageMath with: `GF(MODULUS).primitive_element()`.
const GENERATOR: Fr = Fr::from_raw([0x07, 0x00, 0x00, 0x00]);

/// GENERATOR^t where t * 2^s + 1 = r with t odd. In other words, this is a 2^s root of unity.
/// t = 409655274805673363120685472720202858103411121670017820368325103335302739775
/// `19470B7EFE802F9B36B6675F52C7008234BB3E0CB7ED22AEC65A62A1234BD960`
const ROOT_OF_UNITY: Fr = Fr::from_raw([
    0xC65A62A1234BD960,
    0x34BB3E0CB7ED22AE,
    0x36B6675F52C70082,
    0x19470B7EFE802F9B,
]);

/// 1 / ROOT_OF_UNITY mod r
/// `12BD24C544CCEE9E935FC01BF5106936E62D68F4CC287C26DE5B1F3C2F8A9200`
const ROOT_OF_UNITY_INV: Fr = Fr::from_raw([
    0xDE5B1F3C2F8A9200,
    0xE62D68F4CC287C26,
    0x935FC01BF5106936,
    0x12BD24C544CCEE9E,
]);

/// 1 / 2 mod r
/// E7DB4EA6533AFA906673B0101343B007FC7C3803A0C8238BA7E835A943B73F1
// 0x1a8321a98684c180beebb2019390cad6e0d54b593ee3286311383e6c85936f08
const TWO_INV: Fr = Fr::from_raw([
    0xba7e835a943b73f1,
    0x7fc7c3803a0c8238,
    0x06673b0101343b00,
    0xe7db4ea6533afa9,
]);

// C6A4BE1AB577A673E9D28FF76D10E05F34C2C9AC1E5639B32882FC4D84C3E8E
const ZETA: Fr = Fr::from_raw([
    0x32882FC4D84C3E8E,
    0xF34C2C9AC1E5639B,
    0x3E9D28FF76D10E05,
    0xC6A4BE1AB577A67,
]);

/// Generator of the t-order multiplicative subgroup.
/// Computed by exponentiating Self::MULTIPLICATIVE_GENERATOR by 2^s, where s is Self::S.
const DELTA: Fr = Fr::from_raw([0x1E39A5057D81, 0, 0, 0]);

use crate::{
    field_arithmetic, field_bits, field_common, field_specific, impl_add_binop_specify_output,
    impl_binops_additive, impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_from_u64, impl_sub_binop_specify_output, impl_sum_prod,
};
impl_binops_additive!(Fr, Fr);
impl_binops_multiplicative!(Fr, Fr);
field_common!(
    Fr,
    MODULUS,
    INV,
    MODULUS_STR,
    TWO_INV,
    ROOT_OF_UNITY_INV,
    DELTA,
    ZETA,
    R,
    R2,
    R3
);
impl_from_u64!(Fr, R2);
field_arithmetic!(Fr, MODULUS, INV, dense);
impl_sum_prod!(Fr);

#[cfg(target_pointer_width = "64")]
field_bits!(Fr, MODULUS);
#[cfg(not(target_pointer_width = "64"))]
field_bits!(Fr, MODULUS, MODULUS_LIMBS_32);

impl Fr {
    pub const fn size() -> usize {
        32
    }
}

impl ff::Field for Fr {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

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

    fn double(&self) -> Self {
        self.double()
    }

    #[inline(always)]
    fn square(&self) -> Self {
        self.square()
    }

    /// Returns the multiplicative inverse of the
    /// element. If it is zero, the method fails.
    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }

    fn pow_vartime<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one();
        let mut found_one = false;
        for e in exp.as_ref().iter().rev() {
            for i in (0..64).rev() {
                if found_one {
                    res = res.square();
                }

                if ((*e >> i) & 1) == 1 {
                    found_one = true;
                    res *= self;
                }
            }
        }
        res
    }

    fn sqrt(&self) -> CtOption<Self> {
        // 7fffffff800000007fffffffffffffffde737d56d38bcf4279dce5617e3192a
        let tm1d2 = [
            0x279dce5617e3192a,
            0xfde737d56d38bcf4,
            0x07ffffffffffffff,
            0x7fffffff8000000,
        ];

        ff::helpers::sqrt_tonelli_shanks(self, tm1d2)
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff::helpers::sqrt_ratio_generic(num, div)
    }
}

impl ff::PrimeField for Fr {
    type Repr = [u8; 32];

    const NUM_BITS: u32 = 253;
    const CAPACITY: u32 = 252;
    const MODULUS: &'static str = MODULUS_STR;
    const MULTIPLICATIVE_GENERATOR: Self = GENERATOR;
    const ROOT_OF_UNITY: Self = ROOT_OF_UNITY;
    const ROOT_OF_UNITY_INV: Self = ROOT_OF_UNITY_INV;
    const TWO_INV: Self = TWO_INV;
    const DELTA: Self = DELTA;
    const S: u32 = 5;

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

        // Convert to Montgomery form by computi
        // (a.R^0 * R^2) / R = a.R
        tmp *= &R2;

        CtOption::new(tmp, Choice::from(is_some))
    }

    fn to_repr(&self) -> Self::Repr {
        let tmp: [u64; 4] = (*self).into();
        let mut res = [0; 32];
        res[0..8].copy_from_slice(&tmp[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp[3].to_le_bytes());

        res
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr()[0] & 1)
    }
}

impl FromUniformBytes<64> for Fr {
    /// Converts a 512-bit little endian integer into
    /// an `Fr` by reducing by the modulus.
    fn from_uniform_bytes(bytes: &[u8; 64]) -> Self {
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
}

impl WithSmallOrderMulGroup<3> for Fr {
    const ZETA: Self = ZETA;
}

extend_field_legendre!(Fr);

#[cfg(test)]
mod test {
    use super::*;
    use ff::Field;
    use rand_core::OsRng;

    #[test]
    fn test_zeta() {
        assert_eq!(Fr::ZETA * Fr::ZETA * Fr::ZETA, Fr::ONE);
        assert_ne!(Fr::ZETA * Fr::ZETA, Fr::ONE);
    }

    #[test]
    fn test_sqrt() {
        // NB: TWO_INV is standing in as a "random" field element
        // let v = (Fr::TWO_INV).square().sqrt().unwrap();
        // assert!(v == Fr::TWO_INV || (-v) == Fr::TWO_INV);

        for _ in 0..10000 {
            let a = Fr::random(OsRng);
            let mut b = a;
            b = b.square();

            let b = b.sqrt().unwrap();
            let mut negb = b;
            negb = negb.neg();

            assert!(a == b || a == negb);
        }
    }

    #[test]
    fn test_constants() {
        assert_eq!(
            Fr::MODULUS,
            "0x1cfb69d4ca675f520cce760202687600ff8f87007419047174fd06b52876e7e1",
        );
    }

    #[test]
    fn test_delta() {
        assert_eq!(Fr::DELTA, Fr::MULTIPLICATIVE_GENERATOR.pow([1u64 << Fr::S]));
    }

    #[test]
    fn test_root_of_unity() {
        assert_eq!(Fr::ROOT_OF_UNITY.pow_vartime([1 << Fr::S]), Fr::one());
    }

    #[test]
    fn test_inv_root_of_unity() {
        assert_eq!(Fr::ROOT_OF_UNITY_INV, Fr::ROOT_OF_UNITY.invert().unwrap());
    }

    #[test]
    fn test_field() {
        crate::tests::field::random_field_tests::<Fr>("bandersnatch scalar".to_string());
    }

    #[test]
    fn test_conversion() {
        crate::tests::field::random_conversion_tests::<Fr>("bandersnatch scalar".to_string());
    }

    #[test]
    fn test_serialization() {
        crate::tests::field::random_serialization_test::<Fr>("bandersnatch scalar".to_string());
        #[cfg(feature = "derive_serde")]
        crate::tests::field::random_serde_test::<Fr>("bandersnatch scalar".to_string());
    }

    #[test]
    fn test_quadratic_residue() {
        crate::tests::field::random_quadratic_residue_test::<Fr>();
    }
}
