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
    Fq,
    modulus = "30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47",
    mul_gen = "3",
    zeta = "30644e72e131a0295e6dd9e7e0acccb0c28f069fbb966e3de4bd44e5607cfd48",
    from_uniform = [64, 48],
    endian = "little",
);

extend_field_legendre!(Fq);
impl_binops_calls!(Fq);
impl_binops_additive!(Fq, Fq);
impl_binops_multiplicative!(Fq, Fq);
field_bits!(Fq);
serialize_deserialize_primefield!(Fq);
impl_from_u64!(Fq);

#[cfg(test)]
mod test {
    use super::Fq;
    use crate::{arith_test, constants_test, legendre_test, serde_test, test, test_uniform_bytes};

    constants_test!(Fq);

    arith_test!(Fq);
    legendre_test!(Fq);
    test!(arith, Fq, sqrt_test, 1000);

    serde_test!(Fq PrimeFieldBits);
    test_uniform_bytes!(Fq, 1000, L 64, L 48);
}
