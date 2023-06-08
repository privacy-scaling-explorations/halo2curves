use crate::arithmetic::mul_512;
use crate::arithmetic::sbb;
use crate::arithmetic::CurveEndo;
use crate::arithmetic::EndoParameters;
use crate::bn256::Fq;
use crate::bn256::Fq2;
use crate::bn256::Fr;
use crate::endo;
use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::Curve;
use crate::group::{cofactor::CofactorGroup, prime::PrimeCurveAffine, Group, GroupEncoding};
use crate::{
    batch_add, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output, new_curve_impl,
};
use crate::{Coordinates, CurveAffine, CurveAffineExt, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::convert::TryInto;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    false,
    Fq,
    Fr,
    (G1_GENERATOR_X,G1_GENERATOR_Y),
    G1_B,
    "bn256_g1",
);

new_curve_impl!(
    (pub),
    G2,
    G2Affine,
    false,
    Fq2,
    Fr,
    (G2_GENERATOR_X, G2_GENERATOR_Y),
    G2_B,
    "bn256_g2",
);

impl CurveAffineExt for G1Affine {
    batch_add!();

    fn into_coordinates(self) -> (Self::Base, Self::Base) {
        (self.x, self.y)
    }
}

impl CurveAffineExt for G2Affine {
    batch_add!();

    fn into_coordinates(self) -> (Self::Base, Self::Base) {
        (self.x, self.y)
    }
}

const G1_GENERATOR_X: Fq = Fq::one();
const G1_GENERATOR_Y: Fq = Fq::from_raw([2, 0, 0, 0]);
const G1_B: Fq = Fq::from_raw([3, 0, 0, 0]);

const G2_B: Fq2 = Fq2 {
    c0: Fq::from_raw([
        0x3267e6dc24a138e5,
        0xb5b4c5e559dbefa3,
        0x81be18991be06ac3,
        0x2b149d40ceb8aaae,
    ]),
    c1: Fq::from_raw([
        0xe4a2bd0685c315d2,
        0xa74fa084e52d1852,
        0xcd2cafadeed8fdf4,
        0x009713b03af0fed4,
    ]),
};

const G2_GENERATOR_X: Fq2 = Fq2 {
    c0: Fq::from_raw([
        0x46debd5cd992f6ed,
        0x674322d4f75edadd,
        0x426a00665e5c4479,
        0x1800deef121f1e76,
    ]),
    c1: Fq::from_raw([
        0x97e485b7aef312c2,
        0xf1aa493335a9e712,
        0x7260bfb731fb5d25,
        0x198e9393920d483a,
    ]),
};

const G2_GENERATOR_Y: Fq2 = Fq2 {
    c0: Fq::from_raw([
        0x4ce6cc0166fa7daa,
        0xe3d1e7690c43d37b,
        0x4aab71808dcb408f,
        0x12c85ea5db8c6deb,
    ]),

    c1: Fq::from_raw([
        0x55acdadcd122975b,
        0xbc4b313370b38ef3,
        0xec9e99ad690c3395,
        0x090689d0585ff075,
    ]),
};
// Generated using https://github.com/ConsenSys/gnark-crypto/blob/master/ecc/utils.go
// with `bn256::Fr::ZETA`
// See https://github.com/demining/Endomorphism-Secp256k1/blob/main/README.md
// to have more details about the endomorphism.
const ENDO_PARAMS: EndoParameters = EndoParameters {
    // round(b2/n)
    gamma1: [
        0x7a7bd9d4391eb18du64,
        0x4ccef014a773d2cfu64,
        0x0000000000000002u64,
        0u64,
    ],
    // round(-b1/n)
    gamma2: [0xd91d232ec7e0b3d7u64, 0x0000000000000002u64, 0u64, 0u64],
    b1: [0x8211bbeb7d4f1128u64, 0x6f4d8248eeb859fcu64, 0u64, 0u64],
    b2: [0x89d3256894d213e3u64, 0u64, 0u64, 0u64],
};

endo!(G1, Fr, ENDO_PARAMS);

impl group::cofactor::CofactorGroup for G1 {
    type Subgroup = G1;

    fn clear_cofactor(&self) -> Self {
        *self
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        1.into()
    }
}

impl CofactorGroup for G2 {
    type Subgroup = G2;

    fn clear_cofactor(&self) -> Self {
        // "0x30644e72e131a029b85045b68181585e06ceecda572a2489345f2299c0f9fa8d"
        let e: [u8; 32] = [
            0x30, 0x64, 0x4e, 0x72, 0xe1, 0x31, 0xa0, 0x29, 0xb8, 0x50, 0x45, 0xb6, 0x81, 0x81,
            0x58, 0x5e, 0x06, 0xce, 0xec, 0xda, 0x57, 0x2a, 0x24, 0x89, 0x34, 0x5f, 0x22, 0x99,
            0xc0, 0xf9, 0xfa, 0x8d,
        ];

        // self * COFACTOR_G2
        let mut acc = G2::identity();
        for bit in e
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = G2::conditional_select(&acc, &(acc + self), bit);
        }
        acc
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        unimplemented!();
    }

    fn is_torsion_free(&self) -> Choice {
        // "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001"
        let e: [u8; 32] = [
            0x30, 0x64, 0x4e, 0x72, 0xe1, 0x31, 0xa0, 0x29, 0xb8, 0x50, 0x45, 0xb6, 0x81, 0x81,
            0x58, 0x5d, 0x28, 0x33, 0xe8, 0x48, 0x79, 0xb9, 0x70, 0x91, 0x43, 0xe1, 0xf5, 0x93,
            0xf0, 0x00, 0x00, 0x01,
        ];

        // self * GROUP_ORDER;

        let mut acc = G2::identity();
        for bit in e
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = G2::conditional_select(&acc, &(acc + self), bit);
        }
        acc.is_identity()
    }
}

#[cfg(test)]
mod tests {

    use crate::arithmetic::CurveEndo;
    use crate::bn256::{Fr, G1, G2};
    use crate::CurveExt;
    use ff::Field;
    use ff::PrimeField;
    use ff::WithSmallOrderMulGroup;
    use rand_core::OsRng;

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<G1>();
        crate::tests::curve::curve_tests::<G2>();
    }

    #[test]
    fn test_endo() {
        let g = G1::generator();
        assert_eq!(g * Fr::ZETA, g.endo());
        let g = G2::generator();
        assert_eq!(g * Fr::ZETA, g.endo());
        for _ in 0..100000 {
            let k = Fr::random(OsRng);
            let (k1, k1_neg, k2, k2_neg) = G1::decompose_scalar(&k);
            if k1_neg & k2_neg {
                assert_eq!(k, -Fr::from_u128(k1) + Fr::ZETA * Fr::from_u128(k2))
            } else if k1_neg {
                assert_eq!(k, -Fr::from_u128(k1) - Fr::ZETA * Fr::from_u128(k2))
            } else if k2_neg {
                assert_eq!(k, Fr::from_u128(k1) + Fr::ZETA * Fr::from_u128(k2))
            } else {
                assert_eq!(k, Fr::from_u128(k1) - Fr::ZETA * Fr::from_u128(k2))
            }
        }
    }

    #[test]
    fn test_serialization() {
        crate::tests::curve::random_serialization_test::<G1>();
        crate::tests::curve::random_serialization_test::<G2>();
        #[cfg(feature = "derive_serde")]
        {
            crate::tests::curve::random_serde_test::<G1>();
            crate::tests::curve::random_serde_test::<G2>();
        }
    }
}
