use super::fq::Fq;
use super::fq2::Fq2;
use super::fr::Fr;
use crate::arithmetic::{BaseExt, Coordinates, CurveAffine, CurveExt, Group};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use group::{
    cofactor::{CofactorCurve, CofactorGroup},
    prime::{PrimeCurve, PrimeCurveAffine, PrimeGroup},
    Curve as _, Group as _, GroupEncoding,
};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    Fq,
    Fr,
    32,
    "bn256_g1"
);

new_curve_impl!(
    (pub),
    G2,
    G2Affine,
    Fq2,
    Fr,
    64,
    "bn256_g2"
);

impl G1 {
    const fn curve_constant_b() -> Fq {
        Fq::from_raw([3, 0, 0, 0])
    }

    pub fn generator() -> Self {
        const TWO: Fq = Fq::from_raw([2, 0, 0, 0]);

        Self {
            x: Fq::one(),
            y: TWO,
            z: Fq::one(),
        }
    }
}

impl G1Affine {
    const fn curve_constant_b() -> Fq {
        Fq::from_raw([3, 0, 0, 0])
    }

    pub fn generator() -> Self {
        const TWO: Fq = Fq::from_raw([2, 0, 0, 0]);

        Self {
            x: Fq::one(),
            y: TWO,
            infinity: Choice::from(0u8),
        }
    }
}

impl G2 {
    const fn curve_constant_b() -> Fq2 {
        Fq2 {
            c0: Fq::from_raw([
                0x3267e6dc24a138e5,
                0xb5b4c5e559dbefa3,
                0x81be18991be06ac3,
                0x2b149d40ceb8aaae,
            ]),
            c1: Fq::from_raw([
                0xe4a2bd0685c315d2,
                0xa74fa084e52d1852,
                0xcd2cafadeed8fdf4,
                0x009713b03af0fed4,
            ]),
        }
    }

    pub fn generator() -> Self {
        let x = Fq2 {
            c0: Fq::from_raw([
                0x46debd5cd992f6ed,
                0x674322d4f75edadd,
                0x426a00665e5c4479,
                0x1800deef121f1e76,
            ]),
            c1: Fq::from_raw([
                0x97e485b7aef312c2,
                0xf1aa493335a9e712,
                0x7260bfb731fb5d25,
                0x198e9393920d483a,
            ]),
        };
        let y = Fq2 {
            c0: Fq::from_raw([
                0x4ce6cc0166fa7daa,
                0xe3d1e7690c43d37b,
                0x4aab71808dcb408f,
                0x12c85ea5db8c6deb,
            ]),

            c1: Fq::from_raw([
                0x55acdadcd122975b,
                0xbc4b313370b38ef3,
                0xec9e99ad690c3395,
                0x090689d0585ff075,
            ]),
        };

        Self {
            x,
            y,
            z: Fq2::one(),
        }
    }
}

impl G2Affine {
    const fn curve_constant_b() -> Fq2 {
        G2::curve_constant_b()
    }

    pub fn generator() -> Self {
        let x = Fq2 {
            c0: Fq::from_raw([
                0x46debd5cd992f6ed,
                0x674322d4f75edadd,
                0x426a00665e5c4479,
                0x1800deef121f1e76,
            ]),
            c1: Fq::from_raw([
                0x97e485b7aef312c2,
                0xf1aa493335a9e712,
                0x7260bfb731fb5d25,
                0x198e9393920d483a,
            ]),
        };
        let y = Fq2 {
            c0: Fq::from_raw([
                0x4ce6cc0166fa7daa,
                0xe3d1e7690c43d37b,
                0x4aab71808dcb408f,
                0x12c85ea5db8c6deb,
            ]),

            c1: Fq::from_raw([
                0x55acdadcd122975b,
                0xbc4b313370b38ef3,
                0xec9e99ad690c3395,
                0x090689d0585ff075,
            ]),
        };

        Self {
            x,
            y,
            infinity: Choice::from(0u8),
        }
    }
}
