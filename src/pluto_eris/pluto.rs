use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
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

const PLUTO_GENERATOR_X: Fq = Fq::from_raw([
    0x9ffffcd2ffffffff,
    0xa2a7e8c30006b945,
    0xe4a7a5fe8fadffd6,
    0x443f9a5cda8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);
const PLUTO_GENERATOR_Y: Fq = Fq::from_raw([
    0x0000000000000007,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
    0x0000000000000000,
]);

const PLUTO_B: Fq = Fq::from_raw([0x39, 0, 0, 0, 0, 0, 0]);

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

impl group::cofactor::CofactorGroup for Pluto {
    type Subgroup = Pluto;

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
    Pluto,
    PlutoAffine,
    false,
    Fq,
    Fp,
    (PLUTO_GENERATOR_X,PLUTO_GENERATOR_Y),
    PLUTO_B,
    "pluto",
);

#[test]
fn test_curve() {
    crate::tests::curve::curve_tests::<Pluto>();
}

#[test]
fn test_serialization() {
    crate::tests::curve::random_serialization_test::<Pluto>();
    #[cfg(feature = "derive_serde")]
    crate::tests::curve::random_serde_test::<Pluto>();
}

#[test]
fn test_endo_consistency() {
    let g = Pluto::generator();
    assert_eq!(g * Fp::ZETA, g.endo());
}
