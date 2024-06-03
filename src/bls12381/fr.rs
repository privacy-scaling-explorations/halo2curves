use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

impl_field!(
    bls12381_scalar,
    Fr,
    modulus = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    mul_gen = "7",
    zeta = "ac45a4010001a40200000000ffffffff",
    from_uniform = [64],
    endian = "little",
);

crate::extend_field_legendre!(Fr);
crate::impl_binops_calls!(Fr);
crate::impl_binops_additive!(Fr, Fr);
crate::impl_binops_multiplicative!(Fr, Fr);
crate::field_bits!(Fr);
crate::serialize_deserialize_primefield!(Fr);
crate::impl_from_u64!(Fr);

#[cfg(test)]
mod test {
    use super::*;
    crate::field_testing_suite!(Fr, "field_arithmetic");
    crate::field_testing_suite!(Fr, "conversion");
    crate::field_testing_suite!(Fr, "serialization");
    crate::field_testing_suite!(Fr, "quadratic_residue");
    crate::field_testing_suite!(Fr, "bits");
    crate::field_testing_suite!(Fr, "serialization_check");
    crate::field_testing_suite!(Fr, "constants");
    crate::field_testing_suite!(Fr, "sqrt");
    crate::field_testing_suite!(Fr, "zeta");
    crate::field_testing_suite!(Fr, "from_uniform_bytes", 64);
}
