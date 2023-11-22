use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::Curve;
use crate::group::{prime::PrimeCurveAffine, Group, GroupEncoding};
use crate::hash_to_curve::simple_svdw_hash_to_curve;
use crate::secp256r1::{Fp, Fq};
use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
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
    Secq256r1,
    Secq256r1Affine,
    true,
    Fq,
    Fp,
    (SECQ_GENERATOR_X, SECQ_GENERATOR_Y),
    SECQ_A,
    SECQ_B,
    "secq256r1",
    |curve_id, domain_prefix| simple_svdw_hash_to_curve(curve_id, domain_prefix, Secq256r1::SSVDW_Z),
);

impl group::cofactor::CofactorGroup for Secq256r1 {
    type Subgroup = Secq256r1;

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

impl Secq256r1 {
    const SSVDW_Z: Fq = Fq::ONE;
}

#[cfg(test)]
mod tests {
    use crate::secq256r1::Fq;
    use crate::CurveExt;
    use ff::WithSmallOrderMulGroup;

    use super::Secq256r1;

    #[test]
    fn test_hash_to_curve() {
        crate::tests::curve::hash_to_curve_test::<Secq256r1>();
    }

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<Secq256r1>();
    }

    #[test]
    fn test_endo_consistency() {
        let g = Secq256r1::generator();
        assert_eq!(g * Fq::ZETA, g.endo());
    }

    #[test]
    fn test_serialization() {
        crate::tests::curve::random_serialization_test::<Secq256r1>();
        #[cfg(feature = "derive_serde")]
        crate::tests::curve::random_serde_test::<Secq256r1>();
    }
}
