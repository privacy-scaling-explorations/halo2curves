use core::ops::Shl;
use num_bigint::BigUint;
use num_traits::{One, ToPrimitive};

fn decompose(e: &BigUint, number_of_limbs: usize, limb_size: usize) -> Vec<BigUint> {
    let mask = &(BigUint::one().shl(limb_size) - 1usize);
    (0usize..)
        .step_by(limb_size)
        .take(number_of_limbs)
        .map(|shift| ((e >> shift) & mask))
        .collect::<Vec<_>>()
}

pub(crate) fn big_to_limbs(e: &BigUint, number_of_limbs: usize) -> Vec<u64> {
    decompose(e, number_of_limbs, 64)
        .iter()
        .map(|x| x.to_u64().unwrap())
        .collect()
}

pub(crate) fn big_to_token(e: &BigUint, number_of_limbs: usize) -> proc_macro2::TokenStream {
    let limbs = big_to_limbs(e, number_of_limbs);
    quote::quote! {[#(#limbs,)*]}
}

pub(crate) fn mod_inv(e: &BigUint, modulus: &BigUint) -> BigUint {
    e.modpow(&(modulus - BigUint::from(2u64)), modulus)
}
