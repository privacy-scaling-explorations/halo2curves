use core::convert::TryInto;
use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    extend_field_legendre, field_bits, impl_add_binop_specify_output, impl_binops_additive,
    impl_binops_additive_specify_output, impl_binops_calls, impl_binops_multiplicative,
    impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
    serialize_deserialize_primefield,
};

impl_field!(
    secp256r1_scalar,
    Fq,
    modulus = "ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551",
    mul_gen = "7",
    zeta = "52891d43d946a0354e786d0777fd6aef9405335ce9c83e1d7cbf87ff12884e21",
    from_uniform = [48, 64],
);

extend_field_legendre!(Fq);
impl_binops_calls!(Fq);
impl_binops_additive!(Fq, Fq);
impl_binops_multiplicative!(Fq, Fq);
field_bits!(Fq);
serialize_deserialize_primefield!(Fq);
crate::impl_from_u64!(Fq);

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
    crate::field_testing_suite!(Fq, "from_uniform_bytes", 48, 64);
}
