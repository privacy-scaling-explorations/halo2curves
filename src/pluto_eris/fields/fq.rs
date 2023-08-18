use crate::arithmetic::{adc, mac, sbb};
use crate::ff::{FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::{
    field_arithmetic_7_limbs, field_bits_7_limbs, field_common_7_limbs, impl_from_u64_7_limbs,
};
use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    impl_sum_prod,
};
use core::convert::TryInto;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::slice::Iter;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

/// This represents an element of $\mathbb{F}_r$ where
///
/// `q = 0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5c7a8a6c7be4a775fe8e177fd69ca7e85d60050af41ffffcd300000001`
///
/// is the scalar field of the Pluto curve.
/// The internal representation of this type is seven 64-bit unsigned
/// integers in little-endian order which account for the 446 bits required to be represented.
///`Fq` values are always in Montgomery form; i.e., Fq(a) = aR mod r, with R = 2^448.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
pub struct Fq(pub(crate) [u64; 7]);

/// Constant representing the modulus
/// q = 0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5c7a8a6c7be4a775fe8e177fd69ca7e85d60050af41ffffcd300000001
const MODULUS: Fq = Fq([
    0x1ffffcd300000001,
    0x9ca7e85d60050af4,
    0xe4a775fe8e177fd6,
    0x443f9a5c7a8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);

/// The modulus as u32 limbs.
#[cfg(not(target_pointer_width = "64"))]
const MODULUS_LIMBS_32: [u32; 14] = [
    0x00000001, 0x1ffffcd3, 0x60050af4, 0x9ca7e85d, 0x8e177fd6, 0xe4a775fe, 0x7a8a6c7b, 0x443f9a5c,
    0xf439266f, 0xa803ca76, 0x0d7f70e4, 0x0130e000, 0x00002400, 0x24000000,
];

const MODULUS_STR: &str = "0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5c7a8a6c7be4a775fe8e177fd69ca7e85d60050af41ffffcd300000001";

/// INV = -(q^{-1} mod 2^64) mod 2^64
/// `0x1ffffcd2ffffffff`
const INV: u64 = 0x1ffffcd2ffffffff;

/// Let M be the power of `2^64` nearest to `Self::MODULUS_BITS`. Then `R = M % Self::MODULUS`.
/// `R = 2^448 mod r`
/// `0x3ffffffffff03fff7a9dfffa183e9bf67e576bf526ff2f52242c778a637089cbf6bc60a1d5b8121b768a5725fdcb3532000163afffffff9`
const R: Fq = Fq([
    0x2000163afffffff9,
    0xb768a5725fdcb353,
    0xbf6bc60a1d5b8121,
    0x2242c778a637089c,
    0x67e576bf526ff2f5,
    0xf7a9dfffa183e9bf,
    0x3ffffffffff03ff,
]);

/// `R^2 = 2^896 mod r`
/// `0x50d7c998f46144ee436895a5a630ff544d51e923f64695651da4da1c97f716419bd905e6e4ff6c2bc64e865fe4552ad740808c831022522`
const R2: Fq = Fq([
    0x740808c831022522,
    0xbc64e865fe4552ad,
    0x19bd905e6e4ff6c2,
    0x51da4da1c97f7164,
    0x44d51e923f646956,
    0xe436895a5a630ff5,
    0x050d7c998f46144e,
]);

/// `R^3 = 2^1792 mod r`
/// `0x134bf4061fa616fe429d2480f5681584bc6d6a21ad4d314c2dc504c19a9fc4f4b2773f0ce3a202a0351e5ad3c72210cff3694c87d62c95d4`
const R3: Fq = Fq([
    0xf3694c87d62c95d4,
    0x351e5ad3c72210cf,
    0xb2773f0ce3a202a0,
    0x2dc504c19a9fc4f4,
    0xbc6d6a21ad4d314c,
    0x429d2480f5681584,
    0x134bf4061fa616fe,
]);

// TODO Why 7? 19 is proposed in the original PR
/// `GENERATOR = 7 mod q` is a generator of the `q - 1` order multiplicative
/// subgroup, or in other words a primitive root of the field.
const GENERATOR: Fq = Fq::from_raw([0x07, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]);

/// Size of the 2-adic sub-group of the field.
const S: u32 = 32;

/// GENERATOR^t where t * 2^s + 1 = q
/// with t odd. In other words, this
/// is a 2^s root of unity.
/// `0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5c7a8a6c7be4a775fe8e177fd69ca7e85d60050af41ffffcd3`
const ROOT_OF_UNITY: Fq = Fq::from_raw([
    0x60050af41ffffcd3,
    0x8e177fd69ca7e85d,
    0x7a8a6c7be4a775fe,
    0xf439266f443f9a5c,
    0x0d7f70e4a803ca76,
    0x000024000130e000,
    0x24000000,
]);

/// 1 / 2 mod q
/// `0x12000000000012000098700006bfb8725401e53b7a1c9337a21fcd2e3d45363df253baff470bbfeb4e53f42eb002857a0ffffe6980000001`
const TWO_INV: Fq = Fq::from_raw([
    0x0ffffe6980000001,
    0x4e53f42eb002857a,
    0xf253baff470bbfeb,
    0xa21fcd2e3d45363d,
    0x5401e53b7a1c9337,
    0x0098700006bfb872,
    0x1200000000001200,
]);

/// 1 / ROOT_OF_UNITY mod q
/// `0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5c7a8a6c7be4a775fe8e177fd69ca7e85d60050af41ffffcd200000001`
const ROOT_OF_UNITY_INV: Fq = Fq::from_raw([
    0x1ffffcd200000001,
    0x9ca7e85d60050af4,
    0xe4a775fe8e177fd6,
    0x443f9a5c7a8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);

/// GENERATOR^{2^s} where t * 2^s + 1 = r with t odd. In other words, this is a t root of unity.
/// 0x657946fe07116ceca983fe28713a2b257ab7a7866c95121e727f3776c3e84cb0a14f6a7f83f8cdaeadb479c657bdf2de4589640faf72e67
const DELTA: Fq = Fq::from_raw([
    0xe4589640faf72e67,
    0xeadb479c657bdf2d,
    0x0a14f6a7f83f8cda,
    0xe727f3776c3e84cb,
    0x57ab7a7866c95121,
    0xca983fe28713a2b2,
    0x657946fe07116ce,
]);

/// `ZETA^3 = 1 mod q` where `ZETA^2 != 1 mod r`
const ZETA: Fq = Fq::from_raw([3u64, 0, 0, 0, 0, 0, 0]);

impl_binops_additive!(Fq, Fq);
impl_binops_multiplicative!(Fq, Fq);
field_common_7_limbs!(
    Fq,
    FqRepr,
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
impl_sum_prod!(Fq);
impl_from_u64_7_limbs!(Fq, R2);
field_arithmetic_7_limbs!(Fq, MODULUS, INV, sparse);

#[cfg(target_pointer_width = "64")]
field_bits_7_limbs!(Fq, MODULUS);
#[cfg(not(target_pointer_width = "64"))]
field_bits_7_limbs!(Fq, MODULUS, MODULUS_LIMBS_32);

impl Fq {
    pub const fn size() -> usize {
        56
    }
}

impl ff::Field for Fq {
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

    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    fn invert(&self) -> CtOption<Self> {
        // TODO: Compute this constant.
        let tmp = self.pow([
            0x43e1f593efffffff,
            0x2833e84879b97091,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ]);

        CtOption::new(tmp, !self.ct_eq(&Self::zero()))
    }

    fn sqrt(&self) -> CtOption<Self> {
        /// `(t - 1) // 2` where t * 2^s + 1 = p with t odd.
        const T_MINUS1_OVER2: [u64; 7] = [
            0xb002857a0ffffe69,
            0x470bbfeb4e53f42e,
            0x3d45363df253baff,
            0x7a1c9337a21fcd2e,
            0x06bfb8725401e53b,
            0x0000120000987000,
            0x0000000012000000,
        ];
        ff::helpers::sqrt_tonelli_shanks(self, T_MINUS1_OVER2)
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff::helpers::sqrt_ratio_generic(num, div)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct FqRepr {
    pub repr: [u8; 56],
}

impl FqRepr {
    fn iter(&self) -> Iter<'_, u8> {
        self.repr.iter()
    }
}

impl Default for FqRepr {
    fn default() -> Self {
        FqRepr { repr: [0u8; 56] }
    }
}

impl AsRef<[u8]> for FqRepr {
    fn as_ref(&self) -> &[u8] {
        self.repr.as_ref()
    }
}

impl AsMut<[u8]> for FqRepr {
    fn as_mut(&mut self) -> &mut [u8] {
        self.repr.as_mut()
    }
}
impl From<[u8; 56]> for FqRepr {
    fn from(repr: [u8; 56]) -> Self {
        Self { repr }
    }
}

impl ff::PrimeField for Fq {
    type Repr = FqRepr;

    const NUM_BITS: u32 = 446;
    const CAPACITY: u32 = 446;
    const MODULUS: &'static str = MODULUS_STR;
    const MULTIPLICATIVE_GENERATOR: Self = GENERATOR;
    const ROOT_OF_UNITY: Self = ROOT_OF_UNITY;
    const ROOT_OF_UNITY_INV: Self = ROOT_OF_UNITY_INV;
    const TWO_INV: Self = TWO_INV;
    const DELTA: Self = DELTA;
    const S: u32 = S;

    fn from_repr(repr: Self::Repr) -> CtOption<Self> {
        let mut tmp = Fq([0, 0, 0, 0, 0, 0, 0]);
        let repr = repr.repr;

        tmp.0[0] = u64::from_le_bytes(repr[0..8].try_into().unwrap());
        tmp.0[1] = u64::from_le_bytes(repr[8..16].try_into().unwrap());
        tmp.0[2] = u64::from_le_bytes(repr[16..24].try_into().unwrap());
        tmp.0[3] = u64::from_le_bytes(repr[24..32].try_into().unwrap());
        tmp.0[4] = u64::from_le_bytes(repr[32..40].try_into().unwrap());
        tmp.0[5] = u64::from_le_bytes(repr[40..48].try_into().unwrap());
        tmp.0[6] = u64::from_le_bytes(repr[48..56].try_into().unwrap());

        // Try to subtract the modulus
        let (_, borrow) = sbb(tmp.0[0], MODULUS.0[0], 0);
        let (_, borrow) = sbb(tmp.0[1], MODULUS.0[1], borrow);
        let (_, borrow) = sbb(tmp.0[2], MODULUS.0[2], borrow);
        let (_, borrow) = sbb(tmp.0[3], MODULUS.0[3], borrow);
        let (_, borrow) = sbb(tmp.0[4], MODULUS.0[4], borrow);
        let (_, borrow) = sbb(tmp.0[5], MODULUS.0[5], borrow);
        let (_, borrow) = sbb(tmp.0[6], MODULUS.0[6], borrow);

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
        let tmp = Fq::montgomery_reduce(&[
            self.0[0], self.0[1], self.0[2], self.0[3], self.0[4], self.0[5], self.0[6], 0, 0, 0,
            0, 0, 0, 0,
        ]);

        let mut res = [0; 56];
        res[0..8].copy_from_slice(&tmp.0[0].to_le_bytes());
        res[8..16].copy_from_slice(&tmp.0[1].to_le_bytes());
        res[16..24].copy_from_slice(&tmp.0[2].to_le_bytes());
        res[24..32].copy_from_slice(&tmp.0[3].to_le_bytes());
        res[32..40].copy_from_slice(&tmp.0[3].to_le_bytes());
        res[40..48].copy_from_slice(&tmp.0[3].to_le_bytes());
        res[48..56].copy_from_slice(&tmp.0[3].to_le_bytes());
        res.into()
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr().repr[0] & 1)
    }
}

impl FromUniformBytes<64> for Fq {
    /// Converts a 512-bit little endian integer into
    /// an `Fq` by reducing by the modulus.
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

impl WithSmallOrderMulGroup<3> for Fq {
    const ZETA: Self = ZETA;
}

#[cfg(test)]
mod test {
    use crate::serde::SerdeObject;

    use super::*;
    use ark_std::{end_timer, start_timer};
    use ff::Field;
    use rand::SeedableRng;
    use rand_core::OsRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_sqrt() {
        let v = (Fq::TWO_INV).square().sqrt().unwrap();
        assert!(v == Fq::TWO_INV || (-v) == Fq::TWO_INV);

        for _ in 0..10000 {
            let a = Fq::random(OsRng);
            let mut b = a;
            b = b.square();

            let b = b.sqrt().unwrap();
            let mut negb = b;
            negb = negb.neg();

            assert!(a == b || a == negb);
        }
    }

    #[test]
    fn test_field() {
        crate::tests::field::random_field_tests::<Fq>("bn256 scalar".to_string());
    }

    #[test]
    fn test_delta() {
        assert_eq!(Fq::DELTA, GENERATOR.pow([1u64 << Fq::S]));
        assert_eq!(Fq::DELTA, Fq::MULTIPLICATIVE_GENERATOR.pow([1u64 << Fq::S]));
    }

    #[test]
    fn test_from_u512() {
        assert_eq!(
            Fq::from_raw([
                0x0000000000000000,
                0x0000000000000000,
                0x0000000000000000,
                0x7e7140b5196b9e6f,
                0x9abac9e4157b6172,
                0xf04bc41062fd7322,
                0x1185fa9c9fef6326,
            ]),
            Fq::from_u512([
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
    #[cfg(feature = "bits")]
    fn test_bits() {
        crate::tests::field::random_bits_tests::<Fq>("Fq".to_string());
    }

    #[test]
    fn test_serialization() {
        crate::tests::field::random_serialization_test::<Fq>("Fq".to_string());
        #[cfg(feature = "derive_serde")]
        crate::tests::field::random_serde_test::<Fq>("Fq".to_string());
    }

    fn is_less_than(x: &[u64; 7], y: &[u64; 7]) -> bool {
        match x[6].cmp(&y[6]) {
            core::cmp::Ordering::Less => return true,
            core::cmp::Ordering::Greater => return false,
            _ => {}
        }
        match x[5].cmp(&y[5]) {
            core::cmp::Ordering::Less => return true,
            core::cmp::Ordering::Greater => return false,
            _ => {}
        }
        match x[4].cmp(&y[4]) {
            core::cmp::Ordering::Less => return true,
            core::cmp::Ordering::Greater => return false,
            _ => {}
        }
        match x[3].cmp(&y[3]) {
            core::cmp::Ordering::Less => return true,
            core::cmp::Ordering::Greater => return false,
            _ => {}
        }
        match x[2].cmp(&y[2]) {
            core::cmp::Ordering::Less => return true,
            core::cmp::Ordering::Greater => return false,
            _ => {}
        }
        match x[1].cmp(&y[1]) {
            core::cmp::Ordering::Less => return true,
            core::cmp::Ordering::Greater => return false,
            _ => {}
        }
        x[0].lt(&y[0])
    }

    #[test]
    fn test_serialization_check() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        let start = start_timer!(|| "serialize Fq");
        // failure check
        for _ in 0..1000000 {
            let rand_word = [(); 7].map(|_| rng.next_u64());
            let a = Fq(rand_word);
            let rand_bytes = a.to_raw_bytes();
            match is_less_than(&rand_word, &MODULUS.0) {
                false => {
                    assert!(Fq::from_raw_bytes(&rand_bytes).is_none());
                }
                _ => {
                    assert_eq!(Fq::from_raw_bytes(&rand_bytes), Some(a));
                }
            }
        }
        end_timer!(start);
    }

    #[test]
    fn bench_fq_from_u16() {
        let repeat = 10000000;
        let mut rng = ark_std::test_rng();
        let base = (0..repeat).map(|_| (rng.next_u32() % (1 << 16)) as u64);

        let timer = start_timer!(|| format!("generate {} Bn256 scalar field elements", repeat));
        let _res: Vec<_> = base.map(|b| Fq::from(b)).collect();

        end_timer!(timer);
    }
}
