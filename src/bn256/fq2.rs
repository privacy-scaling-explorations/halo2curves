use super::fq::{Fq, MODULUS_STR, NEGATIVE_ONE};
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::Legendre;
use core::convert::TryInto;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use std::cmp::Ordering;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

use crate::{
    field_ext_common, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output, impl_sum_prod,
};
impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);
impl_sum_prod!(Fq2);

// The cuadratic nonresidue used to construct the extension `Fq2:Fq`;
const U_SQUARE: Fq = NEGATIVE_ONE;

const EXT_ZETA: Fq = Fq::from_raw([
    0x5763473177fffffe,
    0xd4f263f1acdb5c4f,
    0x59e26bcea0d48bac,
    0x0000000000000000,
]);

// The cubic nonresidue used to construct the extension `Fq6:Fq2`
// V_CUBE = u + 9
pub(crate) const V_CUBE_0: Fq =
    //  9
    Fq::from_raw([9, 0, 0, 0]);

pub(crate) const V_CUBE_1: Fq =
    // 1
    Fq::ONE;

field_ext_common!(Fq2, Fq, U_SQUARE, V_CUBE_0, V_CUBE_1, 64, 32, 254, EXT_ZETA, false);

impl Field for Fq2 {
    const ZERO: Self = Self::zero();
    const ONE: Self = Self::one();

    fn random(mut rng: impl RngCore) -> Self {
        Fq2 {
            c0: Fq::random(&mut rng),
            c1: Fq::random(&mut rng),
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
                c0: NEGATIVE_ONE,
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
            let mut a = Fq2::random(&mut rng);
            let mut b = a;
            a.mul_by_nonresidue();
            b.mul_assign(&nqr);

            assert_eq!(a, b);
        }
    }
}
