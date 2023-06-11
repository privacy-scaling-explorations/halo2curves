use ff::PrimeField;
use num_bigint::BigUint;
use num_traits::Num;

pub mod curve;
pub mod field;

pub(crate) fn from_hex<F: PrimeField>(hex: impl AsRef<str>) -> F {
    let hex = hex
        .as_ref()
        .strip_prefix("0x")
        .unwrap_or_else(|| hex.as_ref());
    let oct = BigUint::from_str_radix(hex, 16).unwrap().to_string();
    F::from_str_vartime(oct.as_str()).unwrap()
}
