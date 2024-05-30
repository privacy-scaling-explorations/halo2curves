use crate::ff_ext::ExtField;
use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

impl_field!(
    bn256_base,
    Fq,
    modulus = "30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47",
    mul_gen = "3",
    zeta = "30644e72e131a0295e6dd9e7e0acccb0c28f069fbb966e3de4bd44e5607cfd48",
    from_uniform = [64, 48],
    endian = "little",
);

crate::extend_field_legendre!(Fq);
crate::impl_binops_calls!(Fq);
crate::impl_binops_additive!(Fq, Fq);
crate::impl_binops_multiplicative!(Fq, Fq);
crate::field_bits!(Fq);
crate::serialize_deserialize_primefield!(Fq);
crate::impl_from_u64!(Fq);

use ff::Field;
const NEGATIVE_ONE: Fq = Fq::ZERO.sub_const(&Fq::ONE);
impl ExtField for Fq {
    const NON_RESIDUE: Self = NEGATIVE_ONE;
    fn mul_by_nonresidue(&self) -> Self {
        self.neg()
    }
    fn frobenius_map(&mut self, _: usize) {}
}

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fq, "field_arithmetic");
    crate::field_testing_suite!(Fq, "conversion");
    crate::field_testing_suite!(Fq, "serialization");
    crate::field_testing_suite!(Fq, "quadratic_residue");
    crate::field_testing_suite!(Fq, "bits");
    crate::field_testing_suite!(Fq, "serialization_check");
    crate::field_testing_suite!(Fq, "constants");
    crate::field_testing_suite!(Fq, "sqrt");
    crate::field_testing_suite!(Fq, "zeta");
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 64, 48);
    #[test]
    fn test_fq_mul_nonresidue() {
        let e = Fq::random(rand_core::OsRng);
        let a0 = e.mul_by_nonresidue();
        let a1 = e * Fq::NON_RESIDUE;
        assert_eq!(a0, a1);
    }
}
