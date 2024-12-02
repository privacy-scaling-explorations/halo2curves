use std::convert::TryInto;

use halo2derive::impl_field;
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    extend_field_legendre, field_bits, impl_binops_additive, impl_binops_calls,
    impl_binops_multiplicative, impl_from_bool, impl_from_u64, serialize_deserialize_primefield,
};

// Fq: Vesta base field and Pasta scalar field.
impl_field!(
    vesta_base,
    Fq,
    modulus = "40000000000000000000000000000000224698fc0994a8dd8c46eb2100000001",
    mul_gen = "5",
    zeta = "06819a58283e528e511db4d81cf70f5a0fed467d47c033af2aa9d2e050aa0e4f",
    from_uniform = [48, 64],
    endian = "little",
);

extend_field_legendre!(Fq);
impl_binops_calls!(Fq);
impl_binops_additive!(Fq, Fq);
impl_binops_multiplicative!(Fq, Fq);
impl_from_u64!(Fq);
impl_from_bool!(Fq);
field_bits!(Fq);
serialize_deserialize_primefield!(Fq);
