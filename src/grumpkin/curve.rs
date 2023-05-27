use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::Curve;
use crate::group::{prime::PrimeCurveAffine, Group, GroupEncoding};
use crate::grumpkin::Fq;
use crate::grumpkin::Fr;
use crate::hash_to_curve::svdw_map_to_curve;
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
    (G1_GENERATOR_X, G1_GENERATOR_Y),
    G1_B,
    "grumpkin_g1",
    |curve_id, domain_prefix| svdw_map_to_curve(curve_id, domain_prefix, Fq::ONE),
);

impl CurveAffineExt for G1Affine {
    batch_add!();

    fn into_coordinates(self) -> (Self::Base, Self::Base) {
        (self.x, self.y)
    }
}

// Parameters in montgomery form taken from
// https://github.com/AztecProtocol/barretenberg/blob/97ccf76c42db581a8b8f8bfbcffe8ca015a3dd22/cpp/src/barretenberg/ecc/curves/grumpkin/grumpkin.hpp#L14
const G1_GENERATOR_X: Fq = Fq::one();
const G1_GENERATOR_Y: Fq = Fq([
    0x11b2dff1448c41d8,
    0x23d3446f21c77dc3,
    0xaa7b8cf435dfafbb,
    0x14b34cf69dc25d68,
]);
const G1_B: Fq = Fq([
    0xdd7056026000005a,
    0x223fa97acb319311,
    0xcc388229877910c0,
    0x034394632b724eaa,
]);

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

#[cfg(test)]
mod tests {
    use crate::grumpkin::{Fr, G1};
    use crate::CurveExt;
    use ff::WithSmallOrderMulGroup;

    #[test]
    fn test_hash_to_curve() {
        crate::tests::curve::hash_to_curve_test::<G1>();
    }

    #[test]
    fn test_curve() {
        crate::tests::curve::curve_tests::<G1>();
    }

    #[test]
    fn test_endo_consistency() {
        let g = G1::generator();
        assert_eq!(g * Fr::ZETA, g.endo());
    }

    #[test]
    fn test_serialization() {
        crate::tests::curve::random_serialization_test::<G1>();
        #[cfg(feature = "derive_serde")]
        crate::tests::curve::random_serde_test::<G1>();
    }
}
