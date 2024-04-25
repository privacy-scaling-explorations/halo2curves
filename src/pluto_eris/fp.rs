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
    pluto_eris_fp,
    Fp,
    modulus = "24000000000024000130e0000d7f70e4a803ca76f439266f443f9a5cda8a6c7be4a7a5fe8fadffd6a2a7e8c30006b9459ffffcd300000001",
    mul_gen = "a",
    zeta = "480000000000360001c950000d7ee0e4a803c956d01c903d720dc8ad8b38dffaf50c100004c37ffffffe",
    from_uniform = [64, 72, 112],
);

extend_field_legendre!(Fp);
impl_binops_calls!(Fp);
impl_binops_additive!(Fp, Fp);
impl_binops_multiplicative!(Fp, Fp);
field_bits!(Fp);
serialize_deserialize_primefield!(Fp);
crate::impl_from_u64!(Fp);

#[cfg(test)]
mod test {

    use super::*;
    crate::field_testing_suite!(Fp, "field_arithmetic");
    crate::field_testing_suite!(Fp, "conversion");
    crate::field_testing_suite!(Fp, "serialization");
    crate::field_testing_suite!(Fp, "quadratic_residue");
    crate::field_testing_suite!(Fp, "bits");
    crate::field_testing_suite!(Fp, "serialization_check");
    crate::field_testing_suite!(Fp, "constants");
    crate::field_testing_suite!(Fp, "sqrt");
    crate::field_testing_suite!(Fp, "zeta");
    crate::field_testing_suite!(Fp, "from_uniform_bytes", 64, 72, 112);
}
