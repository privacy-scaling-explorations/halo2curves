use super::fq::Fq;
use crate::ff::{Field, FromUniformBytes, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::quadratic::{QuadExtField, QuadExtFieldArith, SQRT};
use crate::ff_ext::{ExtField, Legendre};
use core::convert::TryInto;
use std::cmp::Ordering;
use subtle::{Choice, CtOption};

crate::impl_binops_additive!(Fq2, Fq2);
crate::impl_binops_multiplicative!(Fq2, Fq2);
crate::impl_binops_calls!(Fq2);
crate::impl_sum_prod!(Fq2);
crate::impl_tower2!(Fq, Fq2);
crate::impl_tower2_from_uniform_bytes!(Fq, Fq2, 128);

pub type Fq2 = QuadExtField<Fq>;
impl QuadExtFieldArith for Fq2 {
    type Base = Fq;
    const SQRT: SQRT<Fq> = SQRT::Algorithm9 {
        q_minus_3_over_4: &[
            0xee7fbfffffffeaaa,
            0x07aaffffac54ffff,
            0xd9cc34a83dac3d89,
            0xd91dd2e13ce144af,
            0x92c6e9ed90d2eb35,
            0x0680447a8e5ff9a6,
        ],
        q_minus_1_over_2: &[
            0xdcff7fffffffd555,
            0x0f55ffff58a9ffff,
            0xb39869507b587b12,
            0xb23ba5c279c2895f,
            0x258dd3db21a5d66b,
            0x0d0088f51cbff34d,
        ],
    };

    fn square_assign(el: &mut QuadExtField<Self::Base>) {
        let a = el.c0 + el.c1;
        let b = el.c0 - el.c1;
        let c = el.c0.double();
        el.c0 = a * b;
        el.c1 = c * el.c1;
    }
}

impl ExtField for Fq2 {
    const NON_RESIDUE: Self = Fq2::new(Fq::ONE, Fq::ONE);

    fn mul_by_nonresidue(&self) -> Self {
        Self {
            c0: self.c0 - self.c1,
            c1: self.c0 + self.c1,
        }
    }

    fn frobenius_map(&mut self, power: usize) {
        if power % 2 != 0 {
            self.conjugate();
        }
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
    crate::field_testing_suite!(
        Fq2,
        "frobenius",
        // Frobenius endomorphism power parameter for extension field
        //  ϕ: E → E
        //  (x, y) ↦ (x^p, y^p)
        // p: modulus of base field (Here, Fq::MODULUS)
        Fq::MODULUS_LIMBS
    );

    #[test]
    fn test_fq2_mul_nonresidue() {
        let e = Fq2::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fq2::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
