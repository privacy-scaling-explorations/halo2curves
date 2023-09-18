use ff::{Field, PrimeField};
use subtle::{Choice, ConstantTimeEq};

pub trait Legendre: Field {
    type BasePrimeField: PrimeField;

    // This is (p-1)/2 where p is the modulus of the base prime field
    fn legendre_exp() -> &'static [u64];

    fn norm(&self) -> Self::BasePrimeField;

    #[inline]
    fn legendre(&self) -> Self::BasePrimeField {
        self.norm().pow(Self::legendre_exp())
    }

    #[inline]
    fn ct_quadratic_residue(&self) -> Choice {
        self.legendre().ct_eq(&Self::BasePrimeField::ONE)
    }

    #[inline]
    fn ct_quadratic_non_residue(&self) -> Choice {
        self.legendre().ct_eq(&-Self::BasePrimeField::ONE)
    }
}

#[macro_export]
macro_rules! prime_field_legendre {
    ($field:ident ) => {
        impl $crate::legendre::Legendre for $field {
            type BasePrimeField = Self;

            #[inline]
            fn legendre_exp() -> &'static [u64] {
                lazy_static::lazy_static! {
                    // (p-1) / 2
                    static ref LEGENDRE_EXP: Vec<u64> =
                        (num_bigint::BigUint::from_bytes_le((-<$field as ff::Field>::ONE).to_repr().as_ref())/2usize).to_u64_digits();
                }
                &*LEGENDRE_EXP
            }

            #[inline]
            fn norm(&self) -> Self::BasePrimeField {
                self.clone()
            }
        }
    };
}
