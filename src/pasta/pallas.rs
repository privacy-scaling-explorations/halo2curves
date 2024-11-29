use super::fp::Fp;
use super::fq::Fq;
use crate::group::{cofactor::CofactorGroup, prime::PrimeCurveAffine, Curve, Group, GroupEncoding};
use crate::{
    impl_binops_additive, impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, new_curve_impl,
};
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, PrimeField, WithSmallOrderMulGroup};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

new_curve_impl!(
    (pub),
    Pallas,
    PastaAffine,
    Fp,
    Fq,
    (- Fp::ONE, Fp::from_raw([2,0,0,0])),
    Fp::ZERO,
    Fp::from_raw([5,0,0,0]),
    "pasta",
    |domain_prefix| crate::hash_to_curve::hash_to_curve(domain_prefix, Pallas::default_hash_to_curve_suite()),
    crate::serde::CompressedFlagConfig::SingleSpare,
    standard_sign
);

impl CofactorGroup for Pallas {
    type Subgroup = Pallas;

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

impl Pallas {
    /// Z = -13
    pub const SVDW_Z: Fp = Fp::from_raw([
        0x992d30ecfffffff4,
        0x224698fc094cf91b,
        0x0000000000000000,
        0x4000000000000000,
    ]);

    fn default_hash_to_curve_suite() -> crate::hash_to_curve::Suite<Self, sha2::Sha256, 48> {
        crate::hash_to_curve::Suite::<Pallas, sha2::Sha256, 48>::new(
            b"pallas:SHA-256_SVDW_RO_",
            Self::SVDW_Z,
            crate::hash_to_curve::Method::SVDW,
        )
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use crate::curve_testing_suite;
    use crate::serde::SerdeObject;
    use group::UncompressedEncoding;
    use rand_core::OsRng;

    curve_testing_suite!(
        Pallas,
        "constants",
        Fp::MODULUS,
        Fp::ZERO,
        Fp::from_raw([5, 0, 0, 0]),
        -Fp::ONE,
        Fp::from_raw([2, 0, 0, 0]),
        Fq::MODULUS
    );

    curve_testing_suite!(Pallas);
    curve_testing_suite!(Pallas, "endo_consistency");
    curve_testing_suite!(Pallas, "ecdsa_example");
}
