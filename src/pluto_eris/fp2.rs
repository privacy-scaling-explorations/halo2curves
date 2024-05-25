use super::fp::Fp;
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::Legendre;
use core::convert::TryInto;
use core::ops::{Add, Neg, Sub};
use rand::RngCore;
use std::cmp::Ordering;
use std::ops::MulAssign;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_calls, impl_binops_multiplicative, impl_binops_multiplicative_mixed,
    impl_sub_binop_specify_output, impl_sum_prod, impl_tower2, impl_tower2_common,
};
impl_tower2_common!(Fp, Fp2, serde);
impl_tower2!(Fp, Fp2, ReprFp2);
impl_binops_additive!(Fp2, Fp2);
impl_binops_multiplicative!(Fp2, Fp2);
impl_binops_calls!(Fp2);
impl_sum_prod!(Fp2);

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

impl Fp2 {
    pub fn double_assign(&mut self) {
        self.c0 = self.c0.double();
        self.c1 = self.c1.double();
    }

    // TODO: This is a naive method using 4 multiplications
    pub fn mul_assign(&mut self, other: &Self) {
        // r0 = s0 * s0 + U_SQUARE * s1 * o1
        // r1 = s0 * o1 - s1 * o0

        let t0 = self.c0 * other.c0;
        let t1 = self.c0 * other.c1;
        let t2 = self.c1 * other.c0;
        let t3 = self.c1 * other.c1;

        self.c0 = t0 + U_SQUARE * t3;
        self.c1 = t1 + t2
    }

    // TODO: This is a naive method using 3 multiplications
    pub fn square_assign(&mut self) {
        // r0 = s0^2 + U_SQUARE * s1^2
        // r1 = 2* s0s1

        let ab = self.c0 * self.c1;
        let a2 = self.c0 * self.c0;
        let b2 = self.c1 * self.c1;

        self.c1 = ab.double();
        self.c0 = a2 + U_SQUARE * b2;
    }

    // conjucate by negating c1
    pub fn conjugate(&mut self) {
        self.c1 = -self.c1;
    }

    pub fn frobenius_map(&mut self, power: usize) {
        //TODO Replace with constant time version if needed
        if power % 2 != 0 {
            self.conjugate()
        }
    }

    /// Multiply this element by cubic nonresidue: V_CUBE = 57/(u+3)
    pub fn mul_by_nonresidue(&mut self) {
        // (x + y * u) * 57/(u + 3)
        self.mul_assign(&super::fp6::V_CUBE)
    }

    pub fn invert(&self) -> CtOption<Self> {
        let mut t1 = self.c1;
        t1 = t1.square();
        t1 *= U_SQUARE;
        let mut t0 = self.c0;
        t0 = t0.square();
        //t0 = c0^2 - U_SQUARE c1^2
        t0 -= &t1;
        t0.invert().map(|t| {
            let mut tmp = Fp2 {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0 *= &t;
            tmp.c1 *= &t;
            tmp.c1 = -tmp.c1;

            tmp
        })
    }

    /// Norm of Fp2 as extension field in u over Fp
    fn norm(&self) -> Fp {
        // norm = self * self.conjugate()
        let t0 = self.c0.square();
        let t1 = self.c1.square() * U_SQUARE;
        t1 - t0
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
        const NEG_ONE: Fp2 = Fp2 {
            c0: Fp::ZERO.sub_const(&Fp::ONE),
            c1: Fp::ZERO,
        };
        if a0 == NEG_ONE {
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
}

const ZETA: Fp = Fp::from_raw([
    0x8ffff80f80000002,
    0xd9fa5d8a200bc439,
    0x1b50d5e1ff708dc8,
    0xf43f8cddf9a5c478,
    0xa803ca76be3924a5,
    0x0130e0000d7f28e4,
    0x2400000000002400,
]);

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
        let nqr = crate::pluto_eris::fp6::V_CUBE;
        for _ in 0..1000 {
            let mut a = Fp2::random(&mut rng);
            let mut b = a;
            a.mul_by_nonresidue();
            b.mul_assign(&nqr);

            assert_eq!(a, b);
        }
    }
}
