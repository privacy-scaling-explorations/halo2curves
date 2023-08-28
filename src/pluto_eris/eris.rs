use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use crate::hash_to_curve::svdw_hash_to_curve;
use crate::pluto_eris::fields::{fp::Fp, fq::Fq};
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

const ERIS_GENERATOR_X: Fq = Fq::from_raw([
    0x1ffffcd2ffffffff,
    0x9ca7e85d60050af4,
    0xe4a775fe8e177fd6,
    0x443f9a5c7a8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);
const ERIS_GENERATOR_Y: Fq = Fq::from_raw([
    0x0000000000000007,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
]);

const ERIS_A: Fq = Fq::from_raw([0, 0, 0, 0, 0, 0, 0]);
const ERIS_B: Fq = Fq::from_raw([0x39, 0, 0, 0, 0, 0, 0]);

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

impl group::cofactor::CofactorGroup for Eris {
    type Subgroup = Eris;

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

new_curve_impl!(
    (pub),
    Eris,
    ErisAffine,
    false,
    Fq,
    Fp,
    (ERIS_GENERATOR_X,ERIS_GENERATOR_Y),
    ERIS_A,
    ERIS_B,
    "eris",
    |curve_id, domain_prefix| svdw_hash_to_curve(curve_id, domain_prefix, Eris::SVDW_Z),
);

impl Eris {
    const SVDW_Z: Fq = Fq::from_raw([0x04, 0, 0, 0, 0, 0, 0]);
}

#[test]
fn test_curve() {
    crate::tests::curve::curve_tests::<Eris>();
}

#[test]
fn test_hash_to_curve() {
    crate::tests::curve::hash_to_curve_test::<Eris>();
}

#[test]
fn test_serialization() {
    crate::tests::curve::random_serialization_test::<Eris>();
    #[cfg(feature = "derive_serde")]
    crate::tests::curve::random_serde_test::<Eris>();
}

#[test]
fn test_endo_consistency() {
    let g = Eris::generator();
    assert_eq!(g * Fp::ZETA, g.endo());
}
