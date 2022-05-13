use crate::secp256k1::Fp;
use crate::secp256k1::Fq;
use crate::{Coordinates, CurveAffine, CurveAffineExt, CurveExt, Group};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{Field, PrimeField};
use group::Curve;
use group::{prime::PrimeCurveAffine, Group as _, GroupEncoding};

use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

impl Secp256k1 {
    fn endomorphism_base(&self) -> Self {
        unimplemented!();
    }
}

impl group::cofactor::CofactorGroup for Secp256k1 {
    type Subgroup = Secp256k1;

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

// Reference: https://neuromancer.sk/std/secg/secp256k1
const SECP_GENERATOR_X: Fp = Fp::from_raw([
    0x59F2815B16F81798,
    0x029BFCDB2DCE28D9,
    0x55A06295CE870B07,
    0x79BE667EF9DCBBAC,
]);
const SECP_GENERATOR_Y: Fp = Fp::from_raw([
    0x9C47D08FFB10D4B8,
    0xFD17B448A6855419,
    0x5DA4FBFC0E1108A8,
    0x483ADA7726A3C465,
]);
const SECP_B: Fp = Fp::from_raw([7, 0, 0, 0]);

new_curve_impl!(
    (pub),
    Secp256k1,
    Secp256k1Affine,
    Secp256k1Compressed,
    Fp,
    Fq,
    (SECP_GENERATOR_X,SECP_GENERATOR_Y),
    SECP_B,
    "secp256k1",
);

impl CurveAffineExt for Secp256k1Affine {
    batch_add!();
}
