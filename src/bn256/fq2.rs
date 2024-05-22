use super::fq::Fq;
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::Legendre;
use core::convert::TryInto;
use core::ops::{Add, Neg, Sub};
use rand::RngCore;
use std::cmp::Ordering;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_calls, impl_binops_multiplicative, impl_binops_multiplicative_mixed,
    impl_sub_binop_specify_output, impl_sum_prod, impl_tower2, impl_tower2_common,
};
impl_tower2_common!(Fq, Fq2, serde);
impl_tower2!(Fq, Fq2, ReprFq2);
impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);
impl_binops_calls!(Fq2);
impl_sum_prod!(Fq2);

impl Fq2 {
    pub fn double_assign(&mut self) {
        self.c0 = self.c0.double();
        self.c1 = self.c1.double();
    }

    pub fn mul_assign(&mut self, other: &Self) {
        let t1 = self.c0 * other.c0;
        let t2 = self.c1 * other.c1;

        self.c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - (t1 + t2);
        self.c0 = t1 - t2;
    }

    pub fn square_assign(&mut self) {
        let a = self.c0 + self.c1;
        let b = self.c0 - self.c1;
        let c = self.c0.double();
        self.c0 = a * b;
        self.c1 = c * self.c1;
    }

    // conjugate by negating c1
    pub fn conjugate(&mut self) {
        self.c1 = -self.c1;
    }

    pub fn frobenius_map(&mut self, power: usize) {
        if power % 2 != 0 {
            self.conjugate()
        }
    }

    /// Multiply this element by quadratic nonresidue 9 + u.
    pub fn mul_by_nonresidue(&self) -> Self {
        // (xu+y)(u+9) = (9x+y)u+(9y-x)
        let t0 = self.c0;
        let t1 = self.c1;

        // 8*x*i + 8*y
        let t = self.double().double().double();

        Self {
            // 9*y
            c0: t.c0 + t0 - t1,
            // (9*x + y)
            c1: t.c1 + t0 + t1,
        }
    }

    pub fn invert(&self) -> CtOption<Self> {
        (self.c0.square() + self.c1.square())
            .invert()
            .map(|t| Self {
                c0: self.c0 * t,
                c1: self.c1 * -t,
            })
    }

    /// Norm of Fq2 as extension field in i over Fq
    #[inline]
    fn norm(&self) -> Fq {
        self.c0.square() + self.c1.square()
    }

    pub fn sqrt(&self) -> CtOption<Self> {
        // Algorithm 9, https://eprint.iacr.org/2012/685.pdf

        if self.is_zero().into() {
            CtOption::new(Self::ZERO, Choice::from(1))
        } else {
            // a1 = self^((q - 3) / 4)
            // 0xc19139cb84c680a6e14116da060561765e05aa45a1c72a34f082305b61f3f51
            let u: [u64; 4] = [
                0x4f082305b61f3f51,
                0x65e05aa45a1c72a3,
                0x6e14116da0605617,
                0x0c19139cb84c680a,
            ];
            let mut a1 = self.pow(u);
            let mut alpha = a1;

            alpha.square_assign();
            alpha.mul_assign(self);
            let mut a0 = alpha;
            a0.frobenius_map(1);
            a0.mul_assign(&alpha);

            let neg1 = Fq2 {
                c0: Fq::ZERO.sub_const(&Fq::ONE),
                c1: Fq::zero(),
            };

            if a0 == neg1 {
                CtOption::new(a0, Choice::from(0))
            } else {
                a1.mul_assign(self);

                if alpha == neg1 {
                    a1.mul_assign(&Fq2 {
                        c0: Fq::zero(),
                        c1: Fq::one(),
                    });
                } else {
                    alpha += &Fq2::ONE;
                    // alpha = alpha^((q - 1) / 2)
                    // 0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3
                    let u: [u64; 4] = [
                        0x9e10460b6c3e7ea3,
                        0xcbc0b548b438e546,
                        0xdc2822db40c0ac2e,
                        0x183227397098d014,
                    ];
                    alpha = alpha.pow(u);
                    a1.mul_assign(&alpha);
                }
                CtOption::new(a1, Choice::from(1))
            }
        }
    }
}

impl FromUniformBytes<96> for Fq2 {
    fn from_uniform_bytes(bytes: &[u8; 96]) -> Self {
        let c0: [u8; 48] = bytes[..48].try_into().unwrap();
        let c1: [u8; 48] = bytes[48..].try_into().unwrap();
        Self::new(Fq::from_uniform_bytes(&c1), Fq::from_uniform_bytes(&c0))
    }
}

const ZETA: Fq = Fq::from_raw([
    0x5763473177fffffe,
    0xd4f263f1acdb5c4f,
    0x59e26bcea0d48bac,
    0x0000000000000000,
]);

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fq2, "field_arithmetic");
    crate::field_testing_suite!(Fq2, "conversion");
    crate::field_testing_suite!(Fq2, "serialization");
    crate::field_testing_suite!(Fq2, "quadratic_residue");
    crate::field_testing_suite!(Fq2, "sqrt");
    crate::field_testing_suite!(Fq2, "zeta", Fq);
    // extension field-specific
    crate::field_testing_suite!(Fq2, "f2_tests", Fq);
    crate::field_testing_suite!(
        Fq2,
        "frobenius",
        // Frobenius endomorphism power parameter for extension field
        //  ϕ: E → E
        //  (x, y) ↦ (x^p, y^p)
        // p: modulus of base field (Here, Fq::MODULUS)
        [
            0x3c208c16d87cfd47,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ]
    );

    #[test]
    fn test_fq2_squaring() {
        let mut a = Fq2 {
            c0: Fq::one(),
            c1: Fq::one(),
        }; // u + 1
        a.square_assign();
        assert_eq!(
            a,
            Fq2 {
                c0: Fq::zero(),
                c1: Fq::one() + Fq::one(),
            }
        ); // 2u

        let mut a = Fq2 {
            c0: Fq::zero(),
            c1: Fq::one(),
        }; // u
        a.square_assign();
        assert_eq!(a, {
            let neg1 = -Fq::one();
            Fq2 {
                c0: neg1,
                c1: Fq::zero(),
            }
        }); // -1
    }

    #[test]
    fn test_fq2_mul_nonresidue() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        let nine = Fq::one().double().double().double() + Fq::one();
        let nqr = Fq2 {
            c0: nine,
            c1: Fq::one(),
        };

        for _ in 0..1000 {
            let e = Fq2::random(&mut rng);

            let a0 = e.mul_by_nonresidue();
            let a1 = e * nqr;

            assert_eq!(a0, a1);
        }
    }
}
