//! This module provides an implementation of the $\mathbb{G}_2$ group of BLS12-381.

use core::cmp;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use group::{prime::PrimeCurveAffine, Curve, Group, GroupEncoding};
use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use pasta_curves::arithmetic::{Coordinates, CurveAffine, CurveExt};
#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

use ff::WithSmallOrderMulGroup;
use ff::{Field, PrimeField};
use group::cofactor::CofactorGroup;

use crate::bls12_381::fp::Fp;
use crate::bls12_381::fp2::Fp2;
use crate::bls12_381::Scalar;
use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

new_curve_impl!(
    (pub),
    G2Projective,
    G2Affine,
    false,
    Fp2,
    Scalar,
    (G2_GENERATOR_X, G2_GENERATOR_Y),
    G2_B,
    "bls12_381_g2",
);

const G2_GENERATOR_X: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0xf5f2_8fa2_0294_0a10,
        0xb3f5_fb26_87b4_961a,
        0xa1a8_93b5_3e2a_e580,
        0x9894_999d_1a3c_aee9,
        0x6f67_b763_1863_366b,
        0x0581_9192_4350_bcd7,
    ]),
    c1: Fp::from_raw_unchecked([
        0xa5a9_c075_9e23_f606,
        0xaaa0_c59d_bccd_60c3,
        0x3bb1_7e18_e286_7806,
        0x1b1a_b6cc_8541_b367,
        0xc2b6_ed0e_f215_8547,
        0x1192_2a09_7360_edf3,
    ]),
};

const G2_GENERATOR_Y: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0x4c73_0af8_6049_4c4a,
        0x597c_fa1f_5e36_9c5a,
        0xe7e6_856c_aa0a_635a,
        0xbbef_b5e9_6e0d_495f,
        0x07d3_a975_f0ef_25a2,
        0x0083_fd8e_7e80_dae5,
    ]),
    c1: Fp::from_raw_unchecked([
        0xadc0_fc92_df64_b05d,
        0x18aa_270a_2b14_61dc,
        0x86ad_ac6a_3be4_eba0,
        0x7949_5c4e_c93d_a33a,
        0xe717_5850_a43c_caed,
        0x0b2b_c2a1_63de_1bf2,
    ]),
};

const G2_B: Fp2 = Fp2 {
    c0: Fp::from_raw_unchecked([
        0xaa27_0000_000c_fff3,
        0x53cc_0032_fc34_000a,
        0x478f_e97a_6b0a_807f,
        0xb1d3_7ebe_e6ba_24d7,
        0x8ec9_733b_bf78_ab2f,
        0x09d6_4551_3d83_de7e,
    ]),
    c1: Fp::from_raw_unchecked([
        0xaa27_0000_000c_fff3,
        0x53cc_0032_fc34_000a,
        0x478f_e97a_6b0a_807f,
        0xb1d3_7ebe_e6ba_24d7,
        0x8ec9_733b_bf78_ab2f,
        0x09d6_4551_3d83_de7e,
    ]),
};

impl G2Projective {
    fn psi(&self) -> G2Projective {
        // 1 / ((u+1) ^ ((q-1)/3))
        let psi_coeff_x = Fp2 {
            c0: Fp::zero(),
            c1: Fp::from_raw_unchecked([
                0x890dc9e4867545c3,
                0x2af322533285a5d5,
                0x50880866309b7e2c,
                0xa20d1b8c7e881024,
                0x14e4f04fe2db9068,
                0x14e56d3f1564853a,
            ]),
        };
        // 1 / ((u+1) ^ (p-1)/2)
        let psi_coeff_y = Fp2 {
            c0: Fp::from_raw_unchecked([
                0x3e2f585da55c9ad1,
                0x4294213d86c18183,
                0x382844c88b623732,
                0x92ad2afd19103e18,
                0x1d794e4fac7cf0b9,
                0x0bd592fc7d825ec8,
            ]),
            c1: Fp::from_raw_unchecked([
                0x7bcfa7a25aa30fda,
                0xdc17dec12a927e7c,
                0x2f088dd86b4ebef1,
                0xd1ca2087da74d4a7,
                0x2da2596696cebc1d,
                0x0e2b7eedbbfd87d2,
            ]),
        };

        G2Projective {
            // x = frobenius(x)/((u+1)^((p-1)/3))
            x: self.x.frobenius_map() * psi_coeff_x,
            // y = frobenius(y)/(u+1)^((p-1)/2)
            y: self.y.frobenius_map() * psi_coeff_y,
            // z = frobenius(z)
            z: self.z.frobenius_map(),
        }
    }

    fn psi2(&self) -> G2Projective {
        // 1 / 2 ^ ((q-1)/3)
        let psi2_coeff_x = Fp2 {
            c0: Fp::from_raw_unchecked([
                0xcd03c9e48671f071,
                0x5dab22461fcda5d2,
                0x587042afd3851b95,
                0x8eb60ebe01bacb9e,
                0x03f97d6e83d050d2,
                0x18f0206554638741,
            ]),
            c1: Fp::zero(),
        };

        G2Projective {
            // x = frobenius^2(x)/2^((p-1)/3); note that q^2 is the order of the field.
            x: self.x * psi2_coeff_x,
            // y = -frobenius^2(y); note that q^2 is the order of the field.
            y: self.y.neg(),
            // z = z
            z: self.z,
        }
    }
    /// Multiply `self` by `crate::BLS_X`, using double and add.
    fn mul_by_x(&self) -> G2Projective {
        let mut xself = G2Projective::identity();
        // NOTE: in BLS12-381 we can just skip the first bit.
        let mut x = crate::bls12_381::BLS_X >> 1;
        let mut acc = *self;
        while x != 0 {
            acc = acc.double();
            if x % 2 == 1 {
                xself += acc;
            }
            x >>= 1;
        }
        // finally, flip the sign
        if crate::bls12_381::BLS_X_IS_NEGATIVE {
            xself = -xself;
        }
        xself
    }
}

impl CofactorGroup for G2Projective {
    type Subgroup = G2Projective;

    fn clear_cofactor(&self) -> Self {
        let t1 = self.mul_by_x(); // [x] P
        let t2 = self.psi(); // psi(P)

        self.double().psi2() // psi^2(2P)
            + (t1 + t2).mul_by_x() // psi^2(2P) + [x^2] P + [x] psi(P)
            - t1 // psi^2(2P) + [x^2 - x] P + [x] psi(P)
            - t2 // psi^2(2P) + [x^2 - x] P + [x - 1] psi(P)
            - self // psi^2(2P) + [x^2 - x - 1] P + [x - 1] psi(P)
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        unimplemented!();
    }

    fn is_torsion_free(&self) -> Choice {
        // Algorithm from Section 4 of https://eprint.iacr.org/2021/1130
        // Updated proof of correctness in https://eprint.iacr.org/2022/352
        //
        // Check that psi(P) == [x] P
        self.psi().ct_eq(&self.mul_by_x())
    }
}

#[cfg(test)]
mod tests {
    use crate::bls12_381::G2Projective;

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<G2Projective>();
    }

    #[test]
    fn test_serialization() {
        crate::tests::curve::random_serialization_test::<G2Projective>();
        #[cfg(feature = "derive_serde")]
        {
            crate::tests::curve::random_serde_test::<G2Projective>();
        }
    }
}
