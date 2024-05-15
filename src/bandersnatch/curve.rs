use crate::arithmetic::mul_512;
use crate::arithmetic::sbb;
use crate::arithmetic::CurveEndo;
use crate::arithmetic::EndoParameters;
use crate::derive::curve::{IDENTITY_MASK, IDENTITY_SHIFT, SIGN_MASK, SIGN_SHIFT};
use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::Curve;
use crate::group::{prime::PrimeCurveAffine, Group, GroupEncoding};
use crate::bandersnatch::Fp;
use crate::bandersnatch::Fr;
use crate::{
    endo, impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};


// TODO: change this to bandersnatch parameters.
const ENDO_PARAMS_BANDERSNATCH: EndoParameters = EndoParameters {
    gamma1: [0xd91d232ec7e0b3d2, 0x2, 0, 0],
    gamma2: [0x5398fd0300ff655f, 0x4ccef014a773d2d2, 0x02, 0],
    b1: [0x89d3256894d213e2, 0, 0, 0],
    b2: [0x0be4e1541221250b, 0x6f4d8248eeb859fd, 0, 0],
};


// endo!(Bandersnatch, Fr, ENDO_PARAMS_BANDERSNATCH);


impl group::cofactor::CofactorGroup for Bandersnatch {
    type Subgroup = Bandersnatch;

    fn clear_cofactor(&self) -> Self {
        self * Fr::from(4)
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



new_curve_impl!(
    (pub),
    Bandersnatch,
    BandersnatchAffine,
    Fp,
    Fr,
    (BANDERSNATCH_GENERATOR_X,BANDERSNATCH_GENERATOR_Y),
    BANDERSNATCH_A,
    BANDERSNATCH_B,
    "secp256r1",
    |domain_prefix| crate::hash_to_curve::hash_to_curve(domain_prefix, Bandersnatch::default_hash_to_curve_suite()),
);



impl Bandersnatch {
    // TODO: switch to bandersnatch param
    // Optimal Z with: <https://datatracker.ietf.org/doc/html/rfc9380#sswu-z-code>
    // Script: https://github.com/privacy-scaling-explorations/halo2curves/pull/139#issuecomment-1965257028
    // Z = 7 (reference: <https://www.rfc-editor.org/rfc/rfc9380.html#section-8.2>)
    const SSVDW_Z: Fp = Fp::from_raw([
        0x0000000000000007,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    ]);

    fn default_hash_to_curve_suite() -> crate::hash_to_curve::Suite<Self, sha2::Sha256, 64> {
        crate::hash_to_curve::Suite::<Bandersnatch, sha2::Sha256, 64>::new(
            b"bandersnatch_XMD:SHA-256_SVDW_RO_", // TODO: investigate what to put here for bandersnatch
            Self::SSVDW_Z,
            crate::hash_to_curve::Method::SVDW,
        )
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use group::UncompressedEncoding;
    crate::curve_testing_suite!(Bandersnatch);
    // crate::curve_testing_suite!(Bandersnatch, "endo_consistency");
    // crate::curve_testing_suite!(Bandersnatch, "endo");
    crate::curve_testing_suite!(
        Bandersnatch,
        "constants",
        Fp::MODULUS,
        BANDERSNATCH_A,
        BANDERSNATCH_B,
        BANDERSNATCH_GENERATOR_X,
        BANDERSNATCH_GENERATOR_Y,
        Fr::MODULUS
    );
}
