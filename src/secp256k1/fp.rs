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
    secp256k1_base,
    Fp,
    modulus = "fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f",
    mul_gen = "3",
    zeta = "7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee",
    from_uniform = [48, 64],
    endian = "little",
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
    use super::Fp;
    use crate::{arith_test, constants_test, legendre_test, serde_test, test, test_uniform_bytes};

    constants_test!(Fp);

    arith_test!(Fp);
    legendre_test!(Fp);
    test!(arith, Fp, sqrt_test, 1000);

    serde_test!(Fp PrimeFieldBits);
    test_uniform_bytes!(Fp, 1000, L 64, L 48);
}
