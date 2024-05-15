use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    extend_field_legendre, field_bits, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_calls, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_from_u64, impl_sub_binop_specify_output,
    serialize_deserialize_primefield,
};

impl_field!(
    bn256_base,
    Fr,
    modulus = "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    mul_gen = "7",
    zeta = "c6a4be1ab577a673e9d28ff76d10e05f34c2c9ac1e5639b32882fc4d84c3e8e",
    from_uniform = [64],
);

extend_field_legendre!(Fr);
impl_binops_calls!(Fr);
impl_binops_additive!(Fr, Fr);
impl_binops_multiplicative!(Fr, Fr);
field_bits!(Fr);
serialize_deserialize_primefield!(Fr);
impl_from_u64!(Fr);

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
