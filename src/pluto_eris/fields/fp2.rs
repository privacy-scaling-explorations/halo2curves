use super::fp::{Fp, MODULUS_STR};
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::Legendre;
use core::convert::TryInto;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::cmp::Ordering;
use std::ops::MulAssign;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

/// -ALPHA is a quadratic non-residue in Fp. Fp2 = Fp[X]/(X^2 + ALPHA)
/// We introduce the variable u such that u^2 = -ALPHA

/// U_SQUARE = -5
/// 0x24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5cda8a6c7be4a7a5fe8fadffd6a2a7e8c30006b9459ffffcd2fffffffc
const U_SQUARE: Fp = Fp::from_raw([
    0x9ffffcd2fffffffc,
    0xa2a7e8c30006b945,
    0xe4a7a5fe8fadffd6,
    0x443f9a5cda8a6c7b,
    0xa803ca76f439266f,
    0x0130e0000d7f70e4,
    0x2400000000002400,
]);

use crate::{
    field_quadratic_ext, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output, impl_sum_prod,
};
impl_binops_additive!(Fp2, Fp2);
impl_binops_multiplicative!(Fp2, Fp2);
impl_sum_prod!(Fp2);

const EXT_ZETA: Fp = Fp::from_raw([
    0x8ffff80f80000002,
    0xd9fa5d8a200bc439,
    0x1b50d5e1ff708dc8,
    0xf43f8cddf9a5c478,
    0xa803ca76be3924a5,
    0x0130e0000d7f28e4,
    0x2400000000002400,
]);

pub(crate) const V_CUBE_0: Fp =
    // 0xcdb6db6db6dc3b6dbda9924971b3a9ace4a7f2a7bcb449573cd928ee056022c3f6072240ebe2483833bf7b35b701d98ddb6da4b5b6db6e8
    Fp::from_raw([
        0xddb6da4b5b6db6e8,
        0x833bf7b35b701d98,
        0x3f6072240ebe2483,
        0x73cd928ee056022c,
        0xce4a7f2a7bcb4495,
        0xdbda9924971b3a9a,
        0x0cdb6db6db6dc3b6,
    ]);
pub(crate) const V_CUBE_1: Fp =
    // 0x7b6db6db6db756db71cc2492776bcc3489319197d79f5f3457b57ef5366ce1a8c6d1148d5a5491bb523fb0536dcde8eeb6db62d36db6db3
    Fp::from_raw([
        0xeb6db62d36db6db3,
        0xb523fb0536dcde8e,
        0x8c6d1148d5a5491b,
        0x457b57ef5366ce1a,
        0x489319197d79f5f3,
        0xb71cc2492776bcc3,
        0x07b6db6db6db756d,
    ]);

field_quadratic_ext!(Fp2, Fp, U_SQUARE, V_CUBE_0, V_CUBE_1, 112, 56, 446, EXT_ZETA);

impl Field for Fp2 {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn random(mut rng: impl RngCore) -> Self {
        Fp2 {
            c0: Fp::random(&mut rng),
            c1: Fp::random(&mut rng),
        }
    }

    fn is_zero(&self) -> Choice {
        self.c0.is_zero() & self.c1.is_zero()
    }

    fn square(&self) -> Self {
        self.square()
    }

    fn double(&self) -> Self {
        self.double()
    }

    fn sqrt(&self) -> CtOption<Self> {
        // Algorithm 10, https://eprint.iacr.org/2012/685.pdf

        // Aux elements. Described in PRECOMPUTATION of Algorithm 10.
        // As element of Fp2: E = 0 +  U *
        // 0x13e275a1fa6a13af7a82a3d83bc9e63a667c70cf991a36e603b21f15823a404a021848271d63f0875d232408689b4c6c67153f9701e19938
        const E: Fp2 = Fp2 {
            c0: Fp::ZERO,
            c1: Fp::from_raw([
                0x67153f9701e19938,
                0x5d232408689b4c6c,
                0x021848271d63f087,
                0x03b21f15823a404a,
                0x667c70cf991a36e6,
                0x7a82a3d83bc9e63a,
                0x13e275a1fa6a13af,
            ]),
        };

        // As element of Fp2: f = 5 + 0 * U
        // 0x5
        const F: Fp2 = Fp2 {
            c0: Fp::from_raw([0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
            c1: Fp::ZERO,
        };
        let neg_one = -Fp2::ONE;

        // Algorithm (not constant time)
        let b = self.pow_vartime([
            // (p-1)/4 =
            // 0x900000000000900004c3800035fdc392a00f29dbd0e499bd10fe69736a29b1ef929e97fa3eb7ff5a8a9fa30c001ae5167ffff34c0000000
            0x67ffff34c0000000,
            0xa8a9fa30c001ae51,
            0xf929e97fa3eb7ff5,
            0xd10fe69736a29b1e,
            0x2a00f29dbd0e499b,
            0x004c3800035fdc39,
            0x0900000000000900,
        ]);

        let b_2 = b.square();
        let mut b_2_q = b_2;
        b_2_q.frobenius_map(1);

        let a0 = b_2_q * b_2;
        if a0 == neg_one {
            CtOption::new(a0, Choice::from(0))
        } else {
            let mut x = b;
            x.frobenius_map(1);
            if x * b == Fp2::ONE {
                let x0 = (b_2 * self).c0.sqrt().unwrap();
                x.c0.mul_assign(x0);
                x.c1.mul_assign(x0);
                CtOption::new(x, Choice::from(1))
            } else {
                let x0 = (self * b_2 * F).sqrt().unwrap();
                x *= x0 * E;
                CtOption::new(x, Choice::from(1))
            }
        }
    }

    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        ff::helpers::sqrt_ratio_generic(num, div)
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fp2, "field_arithmetic");
    crate::field_testing_suite!(Fp2, "conversion");
    crate::field_testing_suite!(Fp2, "serialization");
    crate::field_testing_suite!(Fp2, "quadratic_residue");
    crate::field_testing_suite!(Fp2, "sqrt");
    crate::field_testing_suite!(Fp2, "zeta", Fp);
    // extension field-specific
    crate::field_testing_suite!(Fp2, "f2_tests", Fp);
    crate::field_testing_suite!(
        Fp2,
        "frobenius",
        // Frobenius endomorphism power parameter for extension field
        //  ϕ: E → E
        //  (x, y) ↦ (x^p, y^p)
        // p: modulus of base field (Here, Fp::MODULUS)
        [
            0x9ffffcd300000001,
            0xa2a7e8c30006b945,
            0xe4a7a5fe8fadffd6,
            0x443f9a5cda8a6c7b,
            0xa803ca76f439266f,
            0x0130e0000d7f70e4,
            0x2400000000002400,
        ]
    );

    #[test]
    fn test_fp2_squaring() {
        // u + 1
        let mut a = Fp2 {
            c0: Fp::one(),
            c1: Fp::one(),
        };
        // (u + 1) ^2 = 1 + u^2 + 2u = -4 + 2u
        a.square_assign();
        let minus_4 = -Fp::from(4u64);
        assert_eq!(
            a,
            Fp2 {
                c0: minus_4,
                c1: Fp::one() + Fp::one(),
            }
        );

        // u
        let mut a = Fp2 {
            c0: Fp::zero(),
            c1: Fp::one(),
        };
        // u^2
        a.square_assign();
        assert_eq!(
            a,
            Fp2 {
                c0: U_SQUARE,
                c1: Fp::zero(),
            }
        );
    }

    #[test]
    fn test_fp2_mul_nonresidue() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        let nqr = crate::pluto_eris::fields::fp6::V_CUBE;
        for _ in 0..1000 {
            let mut a = Fp2::random(&mut rng);
            let mut b = a;
            a.mul_by_nonresidue();
            b.mul_assign(&nqr);

            assert_eq!(a, b);
        }
    }
}
