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
    bn256_scalar,
    Fr,
    modulus = "30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001",
    mul_gen = "7",
    zeta = "30644e72e131a029048b6e193fd84104cc37a73fec2bc5e9b8ca0b2d36636f23",
    from_uniform = [64, 48],
    endian = "little",
);

extend_field_legendre!(Fr);
impl_binops_calls!(Fr);
impl_binops_additive!(Fr, Fr);
impl_binops_multiplicative!(Fr, Fr);
field_bits!(Fr);
serialize_deserialize_primefield!(Fr);

#[cfg(feature = "bn256-table")]
pub use table::FR_TABLE;
#[cfg(not(feature = "bn256-table"))]
crate::impl_from_u64!(Fr);
#[cfg(feature = "bn256-table")]
impl From<u64> for Fr {
    fn from(val: u64) -> Fr {
        if val < 65536 {
            FR_TABLE[val as usize]
        } else {
            Self([val, 0, 0, 0]) * Fr::R2
        }
    }
}

#[cfg(feature = "bn256-table")]
#[rustfmt::skip]
mod table;
#[cfg(feature = "bn256-table")]
#[cfg(test)]
mod table_tests;

#[cfg(test)]
mod test {
    use super::Fr;
    use crate::{arith_test, constants_test, legendre_test, serde_test, test, test_uniform_bytes};

    constants_test!(Fr);

    arith_test!(Fr);
    legendre_test!(Fr);
    test!(arith, Fr, sqrt_test, 1000);

    serde_test!(Fr PrimeFieldBits);
    test_uniform_bytes!(Fr, 1000, L 64, L 48);
}
