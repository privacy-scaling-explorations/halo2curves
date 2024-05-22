use super::fq::Fq;
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::Legendre;
use core::convert::TryInto;
use core::ops::{Add, Neg, Sub};
use rand::RngCore;
use std::cmp::Ordering;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_calls, impl_binops_multiplicative, impl_binops_multiplicative_mixed,
    impl_sub_binop_specify_output, impl_sum_prod, impl_tower2, impl_tower2_common,
    impl_tower2_from_uniform_bytes,
};
impl_tower2_common!(Fq, Fq2, serde);
impl_tower2!(Fq, Fq2, ReprFq2);
impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);
impl_binops_calls!(Fq2);
impl_sum_prod!(Fq2);
impl_tower2_from_uniform_bytes!(Fq, Fq2, 128);

impl Fq2 {
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

    #[inline]
    pub fn conjugate(&self) -> Self {
        Self {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    pub fn frobenius_map(&mut self, power: usize) {
        if power % 2 != 0 {
            *self = self.conjugate();
        }
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        Self {
            c0: self.c0 - self.c1,
            c1: self.c0 + self.c1,
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

    #[inline]
    fn norm(&self) -> Fq {
        self.c0.square() + self.c1.square()
    }

    pub fn sqrt(&self) -> CtOption<Self> {
        // Algorithm 9, https://eprint.iacr.org/2012/685.pdf
        CtOption::new(Self::zero(), self.is_zero()).or_else(|| {
            // a1 = self^((p - 3) / 4)
            let a1 = self.pow_vartime([
                0xee7f_bfff_ffff_eaaa,
                0x07aa_ffff_ac54_ffff,
                0xd9cc_34a8_3dac_3d89,
                0xd91d_d2e1_3ce1_44af,
                0x92c6_e9ed_90d2_eb35,
                0x0680_447a_8e5f_f9a6,
            ]);

            // alpha = a1^2 * self = self^((p - 3) / 2 + 1) = self^((p - 1) / 2)
            let alpha = a1.square() * self;

            // x0 = self^((p + 1) / 4)
            let x0 = a1 * self;

            // In the event that alpha = -1, the element is order p - 1 and so
            // we're just trying to get the square of an element of the subfield
            // Fq. This is given by x0 * u, since u = sqrt(-1). Since the element
            // x0 = a + bu has b = 0, the solution is therefore au.
            CtOption::new(
                Self {
                    c0: -x0.c1,
                    c1: x0.c0,
                },
                alpha.ct_eq(&Self::one().neg()),
            )
            // Otherwise, the correct solution is (1 + alpha)^((q - 1) // 2) * x0
            .or_else(|| {
                CtOption::new(
                    (alpha + Self::one()).pow_vartime([
                        0xdcff_7fff_ffff_d555,
                        0x0f55_ffff_58a9_ffff,
                        0xb398_6950_7b58_7b12,
                        0xb23b_a5c2_79c2_895f,
                        0x258d_d3db_21a5_d66b,
                        0x0d00_88f5_1cbf_f34d,
                    ]) * x0,
                    Choice::from(1),
                )
            })
            // Only return the result if it's really the square root (and so
            // self is actually quadratic nonresidue)
            .and_then(|sqrt| CtOption::new(sqrt, sqrt.square().ct_eq(self)))
        })
    }
}

const ZETA: Fq = Fq::from_raw([
    0x2e01fffffffefffe,
    0xde17d813620a0002,
    0xddb3a93be6f89688,
    0xba69c6076a0f77ea,
    0x5f19672fdf76ce51,
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
}
