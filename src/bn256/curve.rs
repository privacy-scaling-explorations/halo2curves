use crate::arithmetic::mul_512;
use crate::bn256::Fq;
use crate::bn256::Fq2;
use crate::bn256::Fr;
use crate::{Coordinates, CurveAffine, CurveAffineExt, CurveExt, Group};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, PrimeField};
use group::Curve;
use group::{cofactor::CofactorGroup, prime::PrimeCurveAffine, Group as _, GroupEncoding};
use pasta_curves::arithmetic::FieldExt;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    batch_add, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output, new_curve_impl,
};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    G1Compressed,
    Fq::size(),
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
    G2Compressed,
    Fq2::size(),
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
const ENDO_G1: [u64; 4] = [
    0x7a7bd9d4391eb18du64,
    0x4ccef014a773d2cfu64,
    0x0000000000000002u64,
    0u64,
];
const ENDO_G2: [u64; 4] = [0xd91d232ec7e0b3d7u64, 0x0000000000000002u64, 0u64, 0u64];
const ENDO_MINUS_B1: [u64; 4] = [0x8211bbeb7d4f1128u64, 0x6f4d8248eeb859fcu64, 0u64, 0u64];
const ENDO_B2: [u64; 4] = [0x89d3256894d213e3u64, 0u64, 0u64, 0u64];
const ENDO_BETA: Fr = Fr::from_raw([
    0x8b17ea66b99c90ddu64,
    0x5bfc41088d8daaa7u64,
    0xb3c4d79d41a91758u64,
    0x0u64,
]);

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

trait CurveEndo: CurveExt {
    fn endomorphism_base(&self) -> Self;
    fn endomorphism_scalars(k: &Self::ScalarExt) -> (u128, u128);
}

impl CurveEndo for G1 {
    fn endomorphism_base(&self) -> Self {
        Self {
            x: self.x * Self::Base::ZETA,
            y: -self.y,
            z: self.z,
        }
    }

    fn endomorphism_scalars(k: &Self::ScalarExt) -> (u128, u128) {
        #[cfg(feature = "asm")]
        let input = Fr::montgomery_reduce(&[k.0[0], k.0[1], k.0[2], k.0[3], 0, 0, 0, 0]).0;

        #[cfg(not(feature = "asm"))]
        let input = Fr::montgomery_reduce(k.0[0], k.0[1], k.0[2], k.0[3], 0, 0, 0, 0).0;

        let c1_512 = mul_512(ENDO_G2, input);
        let c2_512 = mul_512(ENDO_G1, input);

        let c1_hi = [c1_512[4], c1_512[5], c1_512[6], c1_512[7]];
        let c2_hi = [c2_512[4], c2_512[5], c2_512[6], c2_512[7]];

        let q1_512 = mul_512(c1_hi, ENDO_MINUS_B1);
        let q2_512 = mul_512(c2_hi, ENDO_B2);

        let q1_lo = Self::ScalarExt::from_raw([q1_512[0], q1_512[1], q1_512[2], q1_512[3]]);
        let q2_lo = Self::ScalarExt::from_raw([q2_512[0], q2_512[1], q2_512[2], q2_512[3]]);

        let k1 = q2_lo - q1_lo;
        let k2 = (k1 * ENDO_BETA) + k;

        (k2.get_lower_128(), k1.get_lower_128())
    }
}

impl CurveEndo for G2 {
    fn endomorphism_base(&self) -> Self {
        unimplemented!();
    }

    fn endomorphism_scalars(_: &Self::ScalarExt) -> (u128, u128) {
        unimplemented!();
    }
}

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

    use crate::bn256::{
        curve::{CurveEndo, ENDO_BETA},
        Fr, G1Affine, G1, G2,
    };
    use ff::Field;
    use rand_core::OsRng;

    use crate::CurveExt;

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<G1>();
        crate::tests::curve::curve_tests::<G2>();
    }

    #[test]
    fn test_endo_consistency() {
        let g = G1::generator();
        assert_eq!(g * (-ENDO_BETA), g.endo());
    }

    #[test]
    fn test_endomorphism() {
        use crate::FieldExt;

        let scalar = Fr::random(OsRng);
        let point = G1Affine::random(OsRng);

        let expected = point * scalar;

        let (part1, part2) = G1::endomorphism_scalars(&scalar);

        let k1 = Fr::from_u128(part1);
        let k2 = Fr::from_u128(part2);

        let t1 = point * k1;
        let base = G1::endomorphism_base(&point.into());

        let t2 = base * k2;
        let result = t1 + t2;

        let res_affine: G1Affine = result.into();
        let exp_affine: G1Affine = expected.into();

        assert_eq!(res_affine, exp_affine);
    }

    #[test]
    fn test_serialization() {
        crate::tests::curve::random_serialization_test::<G1>();
        crate::tests::curve::random_serialization_test::<G2>();
    }
}

impl group::UncompressedEncoding for G1Affine {
    type Uncompressed = G1Compressed;

    fn from_uncompressed(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn from_uncompressed_unchecked(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        unimplemented!();
    }
}

impl group::UncompressedEncoding for G2Affine {
    type Uncompressed = G2Compressed;

    fn from_uncompressed(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn from_uncompressed_unchecked(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        unimplemented!();
    }
}
