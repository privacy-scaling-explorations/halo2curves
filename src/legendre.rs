use ff::{Field, PrimeField};

pub trait Legendre: Field {
    type BasePrimeField: PrimeField;

    // This is (p-1)/2 where p is the modulus of the base prime field
    fn legendre_exp() -> &'static Vec<u64>;

    fn norm(&self) -> &Self::BasePrimeField;

    fn legendre(&self) -> Self::BasePrimeField {
        self.norm().pow(Self::legendre_exp())
    }
}

#[macro_export]
macro_rules! prime_field_legendre {
    ($field:ident ) => {

    lazy_static::lazy_static! {
        static ref LE_AS_BIGUINT: num_bigint::BigUint = num_bigint::BigUint::from_bytes_le((-<$field as ff::Field>::ONE).to_repr().as_ref())/2usize ;
        static ref LEGENDRE_EXP: Vec<u64> = LE_AS_BIGUINT.to_u64_digits();
    }
        impl crate::legendre::Legendre for $field {
            type BasePrimeField = Self;

            fn legendre_exp() -> &'static Vec<u64> {
                &*LEGENDRE_EXP
            }
            fn norm(&self) -> &Self::BasePrimeField {
                &self
            }
        }
    };
}
