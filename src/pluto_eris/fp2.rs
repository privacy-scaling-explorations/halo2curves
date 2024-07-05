use super::fp::Fp;
use crate::ff::{Field, PrimeField, WithSmallOrderMulGroup};
use crate::ff_ext::quadratic::{QuadExtField, QuadExtFieldArith, SQRT};
use crate::ff_ext::{ExtField, Legendre};
use core::convert::TryInto;
use std::cmp::Ordering;
use subtle::{Choice, CtOption};

crate::impl_binops_additive!(Fp2, Fp2);
crate::impl_binops_multiplicative!(Fp2, Fp2);
crate::impl_binops_calls!(Fp2);
crate::impl_sum_prod!(Fp2);
crate::impl_tower2!(Fp, Fp2);

pub type Fp2 = QuadExtField<Fp>;

impl QuadExtFieldArith for Fp2 {
    type Base = Fp;
    const SQRT: SQRT<Fp> = SQRT::Algorithm10 {
        precompute_e: Fp2 {
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
        },
        precompute_f: Fp2 {
            c0: Fp::from_raw([0x05, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00]),
            c1: Fp::ZERO,
        },
        q_minus_1_over_4: &[
            0x67ffff34c0000000,
            0xa8a9fa30c001ae51,
            0xf929e97fa3eb7ff5,
            0xd10fe69736a29b1e,
            0x2a00f29dbd0e499b,
            0x004c3800035fdc39,
            0x0900000000000900,
        ],
    };
}

impl ExtField for Fp2 {
    const NON_RESIDUE: Self = Fp2 {
        c0: Fp::from_raw([
            0xddb6da4b5b6db6e8,
            0x833bf7b35b701d98,
            0x3f6072240ebe2483,
            0x73cd928ee056022c,
            0xce4a7f2a7bcb4495,
            0xdbda9924971b3a9a,
            0x0cdb6db6db6dc3b6,
        ]),

        c1: Fp::from_raw([
            0xeb6db62d36db6db3,
            0xb523fb0536dcde8e,
            0x8c6d1148d5a5491b,
            0x457b57ef5366ce1a,
            0x489319197d79f5f3,
            0xb71cc2492776bcc3,
            0x07b6db6db6db756d,
        ]),
    };
    fn frobenius_map(&mut self, power: usize) {
        if power % 2 != 0 {
            self.conjugate();
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;
    crate::field_testing_suite!(Fp2, "field_arithmetic");
    crate::field_testing_suite!(Fp2, "conversion");
    crate::field_testing_suite!(Fp2, "serialization");
    crate::field_testing_suite!(Fp2, "quadratic_residue");
    // crate::field_testing_suite!(Fp2, "sqrt");
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
        Fp::MODULUS_LIMBS
    );

    #[test]
    fn test_fq2_mul_nonresidue() {
        let e = Fp2::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fp2::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
