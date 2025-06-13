#[cfg(not(feature = "std"))]
extern crate alloc;
#[cfg(not(feature = "std"))]
use alloc::boxed::Box;

use core::{
    cmp,
    fmt::Debug,
    iter::Sum,
    ops::{Add, Mul, Neg, Sub},
};

use rand_core::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    ff::{Field, PrimeField, WithSmallOrderMulGroup},
    group::{prime::PrimeCurveAffine, Curve, Group, GroupEncoding},
    impl_binops_additive, impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, new_curve_impl,
    secp256k1::{Fp, Fq},
    Coordinates, CurveAffine, CurveExt,
};

const SECQ_GENERATOR_X: Fq = Fq::from_raw([
    0xA24288E37702EDA6,
    0x3134E45A097781A6,
    0xB6B06C87A2CE32E2,
    0x76C39F5585CB160E,
]);

const SECQ_GENERATOR_Y: Fq = Fq::from_raw([
    0xA4120DDAD952677F,
    0xD18983D26E8DC055,
    0xDC2D265A8E82A7F7,
    0x3FFC646C7B2918B5,
]);

const SECQ_A: Fq = Fq::from_raw([0, 0, 0, 0]);
const SECQ_B: Fq = Fq::from_raw([7, 0, 0, 0]);

new_curve_impl!(
    (pub),
    Secq256k1,
    Secq256k1Affine,
    Fq,
    Fp,
    (SECQ_GENERATOR_X, SECQ_GENERATOR_Y),
    SECQ_A,
    SECQ_B,
    "secq256k1",
    |domain_prefix| crate::hash_to_curve::hash_to_curve(domain_prefix, Secq256k1::default_hash_to_curve_suite()),
    crate::serde::CompressedFlagConfig::Extra,
    standard_sign
);

impl group::cofactor::CofactorGroup for Secq256k1 {
    type Subgroup = Secq256k1;

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

impl Secq256k1 {
    const SVDW_Z: Fq = Fq::ONE;

    fn default_hash_to_curve_suite() -> crate::hash_to_curve::Suite<Self, sha2::Sha256, 48> {
        crate::hash_to_curve::Suite::<Self, sha2::Sha256, 48>::new(
            b"secq256k1_XMD:SHA-256_SVDW_RO_",
            Self::SVDW_Z,
            crate::hash_to_curve::Method::SVDW,
        )
    }
}

#[cfg(test)]
mod test {

    #[cfg(feature = "std")]
    use crate::serde::SerdeObject;

    use group::UncompressedEncoding;
    use rand_core::OsRng;

    use super::*;

    crate::curve_testing_suite!(Secq256k1);
    crate::curve_testing_suite!(Secq256k1, "endo_consistency");
    crate::curve_testing_suite!(
        Secq256k1,
        "constants",
        Fq::MODULUS,
        SECQ_A,
        SECQ_B,
        SECQ_GENERATOR_X,
        SECQ_GENERATOR_Y,
        Fp::MODULUS
    );
}
