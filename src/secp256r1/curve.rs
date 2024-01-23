use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use crate::hash_to_curve::sswu_hash_to_curve;
use crate::secp256r1::Fp;
use crate::secp256r1::Fq;
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

impl group::cofactor::CofactorGroup for Secp256r1 {
    type Subgroup = Secp256r1;

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

// Reference: https://neuromancer.sk/std/secg/secp256r1
const SECP_GENERATOR_X: Fp = Fp::from_raw([
    0xF4A13945D898C296,
    0x77037D812DEB33A0,
    0xF8BCE6E563A440F2,
    0x6B17D1F2E12C4247,
]);

const SECP_GENERATOR_Y: Fp = Fp::from_raw([
    0xCBB6406837BF51F5,
    0x2BCE33576B315ECE,
    0x8EE7EB4A7C0F9E16,
    0x4FE342E2FE1A7F9B,
]);

const SECP_A: Fp = Fp::from_raw([
    0xFFFFFFFFFFFFFFFC,
    0x00000000FFFFFFFF,
    0x0000000000000000,
    0xFFFFFFFF00000001,
]);
const SECP_B: Fp = Fp::from_raw([
    0x3BCE3C3E27D2604B,
    0x651D06B0CC53B0F6,
    0xB3EBBD55769886BC,
    0x5AC635D8AA3A93E7,
]);

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    new_curve_impl,
};

new_curve_impl!(
    (pub),
    Secp256r1,
    Secp256r1Affine,
    true,
    Fp,
    Fq,
    (SECP_GENERATOR_X,SECP_GENERATOR_Y),
    SECP_A,
    SECP_B,
    "secp256r1",
    |curve_id, domain_prefix| sswu_hash_to_curve(curve_id, domain_prefix, Secp256r1::SSVDW_Z),
);

impl Secp256r1 {
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
mod test {
    use super::*;
    crate::curve_testing_suite!(Secp256r1);
    crate::curve_testing_suite!(Secp256r1, "ecdsa_example");
    crate::curve_testing_suite!(
        Secp256r1,
        "constants",
        Fp::MODULUS,
        SECP_A,
        SECP_B,
        SECP_GENERATOR_X,
        SECP_GENERATOR_Y,
        Fq::MODULUS
    );
}
