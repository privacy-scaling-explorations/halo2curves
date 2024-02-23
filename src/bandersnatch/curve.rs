use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use crate::hash_to_curve::sswu_hash_to_curve;
use crate::bandersnatch::Fp;
use crate::bandersnatch::Fr;
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

impl group::cofactor::CofactorGroup for Bandersnatch {
    type Subgroup = Bandersnatch;

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

// Reference: https://eprint.iacr.org/2021/1152.pdf
// SW x: a76451786f95a802c0982bbd0abd68e41b92adc86c8859b4f44679b21658710
const BANDERSNATCH_GENERATOR_X: Fp = Fp::from_raw([
    0x4f44679b21658710,
    0x41b92adc86c8859b,
    0x2c0982bbd0abd68e,
    0xa76451786f95a80,
]);

// SW y: 44d150c8b4bd14f79720d021a839e7b7eb4ee43844b30243126a72ac2375490a
const BANDERSNATCH_GENERATOR_Y: Fp = Fp::from_raw([
    0x126a72ac2375490a,
    0xeb4ee43844b30243,
    0x9720d021a839e7b7,
    0x44d150c8b4bd14f7,
]);

//  − 3763200000
// 73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFE1FB22001
const BANDERSNATCH_A: Fp = Fp::from_raw([
    0xFFFFFFFE1FB22001,
    0x53BDA402FFFE5BFE,
    0x3339D80809A1D805,
    0x73EDA753299D7D48,
]);

// − 78675968000000
// 73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFB870D2E00001
const BANDERSNATCH_B: Fp = Fp::from_raw([
    0xFFFFB870D2E00001,
    0x53BDA402FFFE5BFE,
    0x3339D80809A1D805,
    0x73EDA753299D7D48,
]);

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

new_curve_impl!(
    (pub),
    Bandersnatch,
    BandersnatchAffine,
    true,
    Fp,
    Fr,
    (BANDERSNATCH_GENERATOR_X,BANDERSNATCH_GENERATOR_Y),
    BANDERSNATCH_A,
    BANDERSNATCH_B,
    "bandersnatch",
    |curve_id, domain_prefix| sswu_hash_to_curve(curve_id, domain_prefix, Bandersnatch::SSVDW_Z),
);

impl Bandersnatch {
    // TODO: switch to bandersnatch param
    // Optimal Z with: <https://datatracker.ietf.org/doc/html/rfc9380#sswu-z-code>
    // 0xffffffff00000001000000000000000000000000fffffffffffffffffffffff5
    // Z = -10 (reference: <https://www.rfc-editor.org/rfc/rfc9380.html#section-8.2>)
    const SSVDW_Z: Fp = Fp::from_raw([
        0xfffffffffffffff5,
        0x00000000ffffffff,
        0x0000000000000000,
        0xffffffff00000001,
    ]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::group::Curve;
    use crate::bandersnatch::{Fp, Fr, Bandersnatch};
    use ff::FromUniformBytes;
    use rand_core::OsRng;

    #[test]
    fn test_hash_to_curve() {
        crate::tests::curve::hash_to_curve_test::<Bandersnatch>();
    }

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<Bandersnatch>();
    }

    #[test]
    fn test_serialization() {
        crate::tests::curve::random_serialization_test::<Bandersnatch>();
        #[cfg(feature = "derive_serde")]
        crate::tests::curve::random_serde_test::<Bandersnatch>();
    }

    #[test]
    fn ecdsa_example() {
        fn mod_n(x: Fp) -> Fr {
            let mut x_repr = [0u8; 32];
            x_repr.copy_from_slice(x.to_repr().as_ref());
            let mut x_bytes = [0u8; 64];
            x_bytes[..32].copy_from_slice(&x_repr[..]);
            Fr::from_uniform_bytes(&x_bytes)
        }

        let g = Bandersnatch::generator();

        for _ in 0..1000 {
            // Generate a key pair
            let sk = Fr::random(OsRng);
            let pk = (g * sk).to_affine();

            // Generate a valid signature
            // Suppose `m_hash` is the message hash
            let msg_hash = Fr::random(OsRng);

            let (r, s) = {
                // Draw arandomness
                let k = Fr::random(OsRng);
                let k_inv = k.invert().unwrap();

                // Calculate `r`
                let r_point = (g * k).to_affine().coordinates().unwrap();
                let x = r_point.x();
                let r = mod_n(*x);

                // Calculate `s`
                let s = k_inv * (msg_hash + (r * sk));

                (r, s)
            };

            {
                // Verify
                let s_inv = s.invert().unwrap();
                let u_1 = msg_hash * s_inv;
                let u_2 = r * s_inv;

                let v_1 = g * u_1;
                let v_2 = pk * u_2;

                let r_point = (v_1 + v_2).to_affine().coordinates().unwrap();
                let x_candidate = r_point.x();
                let r_candidate = mod_n(*x_candidate);

                assert_eq!(r, r_candidate);
            }
        }
    }
}
