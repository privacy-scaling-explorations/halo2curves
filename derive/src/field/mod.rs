mod arith;
#[cfg(feature = "asm")]
mod asm;

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{Num, One};
use proc_macro::TokenStream;
use proc_macro2::Span;
use quote::quote;
use syn::Token;

struct FieldConfig {
    identifier: String,
    field: syn::Ident,
    modulus: BigUint,
    mul_gen: BigUint,
    zeta: BigUint,
    endian: String,
    from_uniform: Vec<usize>,
}

impl syn::parse::Parse for FieldConfig {
    fn parse(input: syn::parse::ParseStream<'_>) -> syn::Result<Self> {
        let identifier: syn::Ident = input.parse()?;
        let identifier = identifier.to_string();
        input.parse::<syn::Token![,]>()?;

        let field: syn::Ident = input.parse()?;
        input.parse::<syn::Token![,]>()?;

        let get_big = |is_key: &str| -> Result<BigUint, syn::Error> {
            let key: syn::Ident = input.parse()?;
            assert_eq!(key.to_string(), is_key);
            input.parse::<Token![=]>()?;
            let n: syn::LitStr = input.parse()?;
            let n = BigUint::from_str_radix(&n.value(), 16)
                .map_err(|err| syn::Error::new(Span::call_site(), err.to_string()))?;
            input.parse::<Token![,]>()?;
            Ok(n)
        };

        let get_str = |is_key: &str| -> Result<String, syn::Error> {
            let key: syn::Ident = input.parse()?;
            assert_eq!(key.to_string(), is_key);
            input.parse::<Token![=]>()?;
            let n: syn::LitStr = input.parse()?;
            let n = n.value();
            input.parse::<Token![,]>()?;
            Ok(n)
        };

        let get_usize_list = |is_key: &str| -> Result<Vec<usize>, syn::Error> {
            let key: syn::Ident = input.parse()?;
            assert_eq!(key.to_string(), is_key);
            input.parse::<Token![=]>()?;

            // chatgpt
            let content;
            syn::bracketed!(content in input);
            let punctuated: syn::punctuated::Punctuated<syn::LitInt, Token![,]> =
                content.parse_terminated(syn::LitInt::parse)?;
            let values = punctuated
                .into_iter()
                .map(|lit| lit.base10_parse::<usize>())
                .collect::<Result<Vec<_>, _>>()?;
            input.parse::<Token![,]>()?;
            Ok(values)
        };

        let modulus = get_big("modulus")?;
        let mul_gen = get_big("mul_gen")?;
        let zeta = get_big("zeta")?;
        let from_uniform = get_usize_list("from_uniform")?;
        let endian = get_str("endian")?;
        assert!(endian == "little" || endian == "big");
        assert!(input.is_empty());

        Ok(FieldConfig {
            identifier,
            field,
            modulus,
            mul_gen,
            zeta,
            from_uniform,
            endian,
        })
    }
}

pub(crate) fn impl_field(input: TokenStream) -> TokenStream {
    use crate::utils::{big_to_token, mod_inv};
    let FieldConfig {
        identifier,
        field,
        modulus,
        mul_gen,
        zeta,
        from_uniform,
        endian,
    } = syn::parse_macro_input!(input as FieldConfig);
    let _ = identifier;

    let num_bits = modulus.bits() as u32;
    let limb_size = 64;
    let num_limbs = ((num_bits - 1) / limb_size + 1) as usize;
    let size = num_limbs * 8;
    let modulus_limbs = crate::utils::big_to_limbs(&modulus, num_limbs);
    let modulus_str = format!("0x{}", modulus.to_str_radix(16));
    let modulus_limbs_ident = quote! {[#(#modulus_limbs,)*]};

    let modulus_limbs_32 = crate::utils::big_to_limbs_32(&modulus, num_limbs * 2);
    let modulus_limbs_32_ident = quote! {[#(#modulus_limbs_32,)*]};

    let to_token = |e: &BigUint| big_to_token(e, num_limbs);
    let half_modulus = (&modulus - 1usize) >> 1;
    let half_modulus = to_token(&half_modulus);

    // binary modulus
    let t = BigUint::from(1u64) << (num_limbs * limb_size as usize);
    // r1 = mont(1)
    let r1: BigUint = &t % &modulus;
    let mont = |v: &BigUint| (v * &r1) % &modulus;
    // r2 = mont(r)
    let r2: BigUint = (&r1 * &r1) % &modulus;
    // r3 = mont(r^2)
    let r3: BigUint = (&r1 * &r1 * &r1) % &modulus;

    let r1 = to_token(&r1);
    let r2 = to_token(&r2);
    let r3 = to_token(&r3);

    // inv = -(r^{-1} mod 2^64) mod 2^64
    let mut inv64 = 1u64;
    for _ in 0..63 {
        inv64 = inv64.wrapping_mul(inv64);
        inv64 = inv64.wrapping_mul(modulus_limbs[0]);
    }
    inv64 = inv64.wrapping_neg();

    let mut by_inverter_constant: usize = 2;
    loop {
        let t = BigUint::from(1u64) << (62 * by_inverter_constant - 64);
        if t > modulus {
            break;
        }
        by_inverter_constant += 1;
    }

    let mut jacobi_constant: usize = 1;
    loop {
        let t = BigUint::from(1u64) << (64 * jacobi_constant - 31);
        if t > modulus {
            break;
        }
        jacobi_constant += 1;
    }

    let mut s: u32 = 0;
    let mut t = &modulus - BigUint::one();
    while t.is_even() {
        t >>= 1;
        s += 1;
    }

    let two_inv = mod_inv(&BigUint::from(2usize), &modulus);

    let sqrt_impl = {
        if &modulus % 16u64 == BigUint::from(1u64) {
            let tm1o2 = ((&t - 1usize) * &two_inv) % &modulus;
            let tm1o2 = big_to_token(&tm1o2, num_limbs);
            quote! {
                fn sqrt(&self) -> subtle::CtOption<Self> {
                    ff::helpers::sqrt_tonelli_shanks(self, #tm1o2)
                }
            }
        } else if &modulus % 4u64 == BigUint::from(3u64) {
            let exp = (&modulus + 1usize) >> 2;
            let exp = big_to_token(&exp, num_limbs);
            quote! {
                fn sqrt(&self) -> subtle::CtOption<Self> {
                    use subtle::ConstantTimeEq;
                    let t = self.pow(#exp);
                    subtle::CtOption::new(t, t.square().ct_eq(self))
                }
            }
        } else {
            panic!("unsupported modulus")
        }
    };

    let root_of_unity = mul_gen.modpow(&t, &modulus);
    let root_of_unity_inv = mod_inv(&root_of_unity, &modulus);
    let delta = mul_gen.modpow(&(BigUint::one() << s), &modulus);

    let root_of_unity = to_token(&mont(&root_of_unity));
    let root_of_unity_inv = to_token(&mont(&root_of_unity_inv));
    let two_inv = to_token(&mont(&two_inv));
    let mul_gen = to_token(&mont(&mul_gen));
    let delta = to_token(&mont(&delta));
    let zeta = to_token(&mont(&zeta));

    let endian = match endian.as_str() {
        "little" => {
            quote! { LE }
        }
        "big" => {
            quote! { BE }
        }
        _ => {
            unreachable!()
        }
    };

    let impl_field = quote! {
        #[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
        pub struct #field(pub(crate) [u64; #num_limbs]);

        impl core::fmt::Debug for #field {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                use ff::PrimeField;
                let tmp = self.to_repr();
                write!(f, "0x")?;
                for &b in tmp.as_ref().iter().rev() {
                    write!(f, "{:02x}", b)?;
                }
                Ok(())
            }
        }

        impl ConstantTimeEq for #field {
            fn ct_eq(&self, other: &Self) -> Choice {
                Choice::from(
                    self.0
                        .iter()
                        .zip(other.0)
                        .all(|(a, b)| bool::from(a.ct_eq(&b))) as u8,
                )
            }
        }

        impl ConditionallySelectable for #field {
            fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
                let limbs = (0..#num_limbs)
                    .map(|i| u64::conditional_select(&a.0[i], &b.0[i], choice))
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                #field(limbs)
            }
        }

        impl core::cmp::PartialOrd for #field {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl core::cmp::Ord for #field {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                use ff::PrimeField;
                let left = self.to_repr();
                let right = other.to_repr();
                left.as_ref().iter()
                    .zip(right.as_ref().iter())
                    .rev()
                    .find_map(|(left_byte, right_byte)| match left_byte.cmp(right_byte) {
                        core::cmp::Ordering::Equal => None,
                        res => Some(res),
                    })
                    .unwrap_or(core::cmp::Ordering::Equal)
            }
        }

        impl<T: ::core::borrow::Borrow<#field>> ::core::iter::Sum<T> for #field {
            fn sum<I: Iterator<Item = T>>(iter: I) -> Self {
                iter.fold(Self::zero(), |acc, item| acc + item.borrow())
            }
        }

        impl<T: ::core::borrow::Borrow<#field>> ::core::iter::Product<T> for #field {
            fn product<I: Iterator<Item = T>>(iter: I) -> Self {
                iter.fold(Self::one(), |acc, item| acc * item.borrow())
            }
        }

        impl crate::encoding::endian::EndianRepr for #field {
            const ENDIAN: crate::encoding::endian::Endian = crate::encoding::endian::Endian::#endian;

            fn to_bytes(&self) -> Vec<u8> {
                self.to_bytes().to_vec()
            }

            fn from_bytes(bytes: &[u8]) -> subtle::CtOption<Self> {
                #field::from_bytes(bytes[..#field::SIZE].try_into().unwrap())
            }
        }

        impl #field {
            pub const SIZE: usize = #num_limbs * 8;
            pub const NUM_LIMBS: usize = #num_limbs;
            pub(crate) const MODULUS_LIMBS: [u64; Self::NUM_LIMBS] = #modulus_limbs_ident;
            pub(crate) const MODULUS_LIMBS_32: [u32; Self::NUM_LIMBS*2] = #modulus_limbs_32_ident;
            const R: Self = Self(#r1);
            const R2: Self = Self(#r2);
            const R3: Self = Self(#r3);

            /// Returns zero, the additive identity.
            #[inline(always)]
            pub const fn zero() -> #field {
                #field([0; Self::NUM_LIMBS])
            }

            /// Returns one, the multiplicative identity.
            #[inline(always)]
            pub const fn one() -> #field {
                Self::R
            }

            /// Converts from an integer represented in little endian
            /// into its (congruent) `$field` representation.
            pub const fn from_raw(val: [u64; Self::NUM_LIMBS]) -> Self {
                Self(val).mul_const(&Self::R2)
            }

            /// Attempts to convert a <#endian>-endian byte representation of
            /// a scalar into a `$field`, failing if the input is not canonical.
            pub fn from_bytes(bytes: &[u8; Self::SIZE]) -> subtle::CtOption<Self> {
                use crate::encoding::endian::EndianRepr;
                let mut el = #field::default();
                #field::ENDIAN.from_bytes(bytes, &mut el.0);
                subtle::CtOption::new(el * Self::R2, subtle::Choice::from(Self::is_less_than_modulus(&el.0) as u8))
            }


            /// Converts an element of `$field` into a byte representation in
            /// <#endian>-endian byte order.
            pub fn to_bytes(&self) -> [u8; Self::SIZE] {
                use crate::encoding::endian::EndianRepr;
                let el = self.from_mont();
                let mut res = [0; Self::SIZE];
                #field::ENDIAN.to_bytes(&mut res, &el);
                res.into()
            }


            // Returns the Jacobi symbol, where the numerator and denominator
            // are the element and the characteristic of the field, respectively.
            // The Jacobi symbol is applicable to odd moduli
            // while the Legendre symbol is applicable to prime moduli.
            // They are equivalent for prime moduli.
            #[inline(always)]
            fn jacobi(&self) -> i64 {
                crate::ff_ext::jacobi::jacobi::<#jacobi_constant>(&self.0, &#modulus_limbs_ident)
            }


            #[inline(always)]
            pub(crate) fn is_less_than_modulus(limbs: &[u64; Self::NUM_LIMBS]) -> bool {
                let borrow = limbs.iter().enumerate().fold(0, |borrow, (i, limb)| {
                    crate::arithmetic::sbb(*limb, Self::MODULUS_LIMBS[i], borrow).1
                });
                (borrow as u8) & 1 == 1
            }

            /// Returns whether or not this element is strictly lexicographically
            /// larger than its negation.
            pub fn lexicographically_largest(&self) -> Choice {
                const HALF_MODULUS: [u64; #num_limbs]= #half_modulus;
                let tmp = self.from_mont();
                let borrow = tmp
                    .into_iter()
                    .zip(HALF_MODULUS.into_iter())
                    .fold(0, |borrow, (t, m)| crate::arithmetic::sbb(t, m, borrow).1);
                !Choice::from((borrow as u8) & 1)
            }
        }

        impl ff::Field for #field {
            const ZERO: Self = Self::zero();
            const ONE: Self = Self::one();

            fn random(mut rng: impl RngCore) -> Self {
                let mut wide = [0u8; Self::SIZE * 2];
                rng.fill_bytes(&mut wide);
                <#field as ff::FromUniformBytes<{ #field::SIZE * 2 }>>::from_uniform_bytes(&wide)
            }

            #[inline(always)]
            #[must_use]
            fn double(&self) -> Self {
                self.double()
            }

            #[inline(always)]
            #[must_use]
            fn square(&self) -> Self {
                self.square()
            }

            // Returns the multiplicative inverse of the element. If it is zero, the method fails.
            #[inline(always)]
            fn invert(&self) -> CtOption<Self> {
                const BYINVERTOR: crate::ff_ext::inverse::BYInverter<#by_inverter_constant> =
                crate::ff_ext::inverse::BYInverter::<#by_inverter_constant>::new(&#modulus_limbs_ident, &#r2);

                if let Some(inverse) = BYINVERTOR.invert::<{ Self::NUM_LIMBS }>(&self.0) {
                    subtle::CtOption::new(Self(inverse), subtle::Choice::from(1))
                } else {
                    subtle::CtOption::new(Self::zero(), subtle::Choice::from(0))
                }
            }

            #sqrt_impl

            fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
                ff::helpers::sqrt_ratio_generic(num, div)
            }
        }
    };

    let impl_prime_field = quote! {

        // TODO use ::core::borrow::Borrow or AsRef
        impl From<#field> for crate::encoding::Repr<{ #field::SIZE }> {
            fn from(value: #field) -> crate::encoding::Repr<{ #field::SIZE }> {
                use ff::PrimeField;
                value.to_repr()
            }
        }

        impl<'a> From<&'a #field> for crate::encoding::Repr<{ #field::SIZE }> {
            fn from(value: &'a #field) -> crate::encoding::Repr<{ #field::SIZE }> {
                use ff::PrimeField;
                value.to_repr()
            }
        }

        impl ff::PrimeField for #field {
            const NUM_BITS: u32 = #num_bits;
            const CAPACITY: u32 = #num_bits-1;
            const TWO_INV :Self = Self(#two_inv);
            const MULTIPLICATIVE_GENERATOR: Self = Self(#mul_gen);
            const S: u32 = #s;
            const ROOT_OF_UNITY: Self = Self(#root_of_unity);
            const ROOT_OF_UNITY_INV: Self = Self(#root_of_unity_inv);
            const DELTA: Self = Self(#delta);
            const MODULUS: &'static str = #modulus_str;

            type Repr = crate::encoding::Repr<{ #field::SIZE }>;

            fn from_u128(v: u128) -> Self {
                Self::R2 * Self(
                    [v as u64, (v >> 64) as u64]
                        .into_iter()
                        .chain(std::iter::repeat(0))
                        .take(Self::NUM_LIMBS)
                        .collect::<Vec<_>>()
                        .try_into()
                        .unwrap(),
                )
            }

            fn from_repr(repr: Self::Repr) -> subtle::CtOption<Self> {
                let mut el = #field::default();
                crate::encoding::endian::Endian::LE.from_bytes(repr.as_ref(), &mut el.0);
                subtle::CtOption::new(el * Self::R2, subtle::Choice::from(Self::is_less_than_modulus(&el.0) as u8))
            }

            fn to_repr(&self) -> Self::Repr {
                use crate::encoding::endian::Endian;
                let el = self.from_mont();
                let mut res = [0; #size];
                crate::encoding::endian::Endian::LE.to_bytes(&mut res, &el);
                res.into()
            }

            fn is_odd(&self) -> Choice {
                Choice::from(self.to_repr()[0] & 1)
            }
        }
    };

    #[cfg(feature = "asm")]
    let impl_arith = {
        if num_limbs == 4 && num_bits < 256 {
            println!("implementing asm, {}", identifier);
            asm::limb4::impl_arith(&field, inv64)
        } else {
            arith::impl_arith(&field, num_limbs, inv64)
        }
    };
    #[cfg(not(feature = "asm"))]
    let impl_arith = arith::impl_arith(&field, num_limbs, inv64);

    let impl_arith_always_const = arith::impl_arith_always_const(&field, num_limbs, inv64);

    let impl_from_uniform_bytes = from_uniform
        .iter()
        .map(|input_size| {
            assert!(*input_size >= size);
            assert!(*input_size <= size*2);
            quote! {
                impl ff::FromUniformBytes<#input_size> for #field {
                    fn from_uniform_bytes(bytes: &[u8; #input_size]) -> Self {
                        let mut wide = [0u8; Self::SIZE * 2];
                        wide[..#input_size].copy_from_slice(bytes);
                        let (a0, a1) = wide.split_at(Self::SIZE);

                        let a0: [u64; Self::NUM_LIMBS] = (0..Self::NUM_LIMBS)
                            .map(|off| u64::from_le_bytes(a0[off * 8..(off + 1) * 8].try_into().unwrap()))
                            .collect::<Vec<_>>()
                            .try_into()
                            .unwrap();
                        let a0 = #field(a0);

                        let a1: [u64; Self::NUM_LIMBS] = (0..Self::NUM_LIMBS)
                            .map(|off| u64::from_le_bytes(a1[off * 8..(off + 1) * 8].try_into().unwrap()))
                            .collect::<Vec<_>>()
                            .try_into()
                            .unwrap();
                        let a1 = #field(a1);

                        // enforce non assembly impl since asm is likely to be optimized for sparse fields
                        a0.mul_const(&Self::R2) + a1.mul_const(&Self::R3)

                    }
                }
            }
        })
        .collect::<proc_macro2::TokenStream>();

    let impl_zeta = quote! {
        impl ff::WithSmallOrderMulGroup<3> for #field {
            const ZETA: Self = Self(#zeta);
        }
    };

    let output = quote! {
        #impl_arith
        #impl_arith_always_const
        #impl_field
        #impl_prime_field
        #impl_from_uniform_bytes
        #impl_zeta
    };

    output.into()
}
