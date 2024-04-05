#[macro_export]
macro_rules! impl_from_uniform_bytes {
    ($field:ident, $size:expr) => {
        impl FromUniformBytes<$size> for $field {
            fn from_uniform_bytes(bytes: &[u8; $size]) -> Self {
                let limbs = (0..($size - 1) / 8 + 1)
                    .map(|i| u64::from_le_bytes(bytes[i * 8..(i + 1) * 8].try_into().unwrap()))
                    .collect::<Vec<_>>();

                let lo = (limbs[..NUM_LIMBS]).try_into().unwrap();
                let hi = limbs
                    .into_iter()
                    .skip(NUM_LIMBS)
                    .map(|limb| limb)
                    .chain(std::iter::repeat(0))
                    .take(NUM_LIMBS)
                    .collect::<Vec<_>>();

                Self::montgomery_form(lo, R2) + Self::montgomery_form(hi.try_into().unwrap(), R3)
            }
        }
    };
}

#[macro_export]
macro_rules! pow_vartime {
    (dense) => {
        fn pow_vartime<S: AsRef<[u64]>>(&self, exp: S) -> Self {
            let mut res = Self::one();
            let mut found_one = false;
            for e in exp.as_ref().iter().rev() {
                for i in (0..64).rev() {
                    if found_one {
                        res = res.square();
                    }

                    if ((*e >> i) & 1) == 1 {
                        found_one = true;
                        res *= self;
                    }
                }
            }
            res
        }
    };
    (sparse) => {
        // Uses fallback implementation
    };
}

#[macro_export]
macro_rules! impl_from_u64 {
    ($field:ident) => {
        impl From<u64> for $field {
            fn from(val: u64) -> $field {
                let limbs = std::iter::once(val)
                    .chain(std::iter::repeat(0))
                    .take(NUM_LIMBS)
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();

                $field(limbs) * R2
            }
        }
    };
}

#[macro_export]
macro_rules! impl_prime_field {
    (
        $field:ident,
        $repr:ty,
        $endian:tt
    ) => {
        fn to_repr(el: &$field) -> $repr {
            let tmp: [u64; NUM_LIMBS] = (*el).into();
            let mut res = [0; $field::SIZE];
            tmp.iter().enumerate().for_each(|(i, limb)| {
                let off = i * 8;
                res[off..off + 8].copy_from_slice(&limb.to_le_bytes());
            });

            res.into()
        }

        fn from_repr(repr: $repr) -> CtOption<$field> {
            let mut tmp = $field::default();
            let borrow = tmp.0.iter_mut().enumerate().fold(0, |borrow, (i, limb)| {
                let off = i * 8;
                *limb = u64::from_le_bytes(repr[off..off + 8].try_into().unwrap());
                sbb(*limb, MODULUS.0[i], borrow).1
            });

            // If the element is smaller than MODULUS then the
            // subtraction will underflow, producing a borrow value
            // of 0xffff...ffff. Otherwise, it'll be zero.
            let is_some = (borrow as u8) & 1;

            // Convert to Montgomery form by computing
            // (a.R^0 * R^2) / R = a.R
            tmp *= &R2;

            CtOption::new(tmp, Choice::from(is_some))
        }

        macro_rules! impl_repr {
            ($field_inner:ident,$repr_inner:ty, le) => {
                fn to_repr(&self) -> $repr {
                    to_repr(self)
                }

                fn from_repr(repr: $repr) -> CtOption<Self> {
                    from_repr(repr)
                }
            };

            ($field_inner:ident,$repr_inner:ty, be) => {
                fn to_repr(&self) -> $repr {
                    let mut repr = to_repr(self);
                    repr.reverse();
                    repr
                }

                fn from_repr(mut repr: $repr) -> CtOption<Self> {
                    repr.reverse();
                    from_repr(repr)
                }
            };
        }

        impl ff::PrimeField for $field {
            type Repr = $repr;

            const MODULUS: &'static str = MODULUS_STR;
            const MULTIPLICATIVE_GENERATOR: Self = MULTIPLICATIVE_GENERATOR;
            const TWO_INV: Self = TWO_INV;
            const ROOT_OF_UNITY: Self = ROOT_OF_UNITY;
            const ROOT_OF_UNITY_INV: Self = ROOT_OF_UNITY_INV;
            const DELTA: Self = DELTA;
            const NUM_BITS: u32 = NUM_BITS;
            const CAPACITY: u32 = Self::NUM_BITS - 1;
            const S: u32 = S;

            fn from_u128(v: u128) -> Self {
                Self::from_raw(
                    [v as u64, (v >> 64) as u64]
                        .into_iter()
                        .chain(std::iter::repeat(0))
                        .take(NUM_LIMBS)
                        .collect::<Vec<_>>()
                        .try_into()
                        .unwrap(),
                )
            }

            fn is_odd(&self) -> Choice {
                Choice::from(self.to_repr()[0] & 1)
            }

            impl_repr!($field, $repr, $endian);
        }

        impl From<$field> for $repr {
            fn from(value: $field) -> $repr {
                value.to_repr()
            }
        }

        impl<'a> From<&'a $field> for $repr {
            fn from(value: &'a $field) -> $repr {
                value.to_repr()
            }
        }

        impl WithSmallOrderMulGroup<3> for $field {
            const ZETA: Self = ZETA;
        }
    };
}

#[macro_export]
macro_rules! impl_field {
    (
        $field:ident,
        $field_type:ident
    ) => {
        // Number of 64 bit limbs to represent the field element
        pub(crate) const NUM_LIMBS: usize = ((NUM_BITS - 1) / 64 + 1) as usize;

        // The internal representation of this type is four 64-bit unsigned
        // integers in little-endian order. `Fq` values are always in
        // Montgomery form; i.e., Fq(a) = aR mod q, with R = 2^256.
        #[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
        pub struct $field(pub(crate) [u64; NUM_LIMBS]);

        /// Bernstein-Yang modular multiplicative inverter created for the modulus equal to
        /// the characteristic of the field to invert positive integers in the Montgomery form.
        const BYINVERTOR: $crate::ff_ext::inverse::BYInverter<BYIL> =
            $crate::ff_ext::inverse::BYInverter::<BYIL>::new(&MODULUS.0, &R2.0);

        impl std::fmt::Debug for $field {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                let tmp = self.to_repr();
                write!(f, "0x")?;
                for &b in tmp.iter().rev() {
                    write!(f, "{:02x}", b)?;
                }
                Ok(())
            }
        }

        impl ConstantTimeEq for $field {
            fn ct_eq(&self, other: &Self) -> Choice {
                Choice::from(
                    self.0
                        .iter()
                        .zip(other.0)
                        .all(|(a, b)| bool::from(a.ct_eq(&b))) as u8,
                )
            }
        }

        impl ConditionallySelectable for $field {
            fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
                let limbs = (0..NUM_LIMBS)
                    .map(|i| u64::conditional_select(&a.0[i], &b.0[i], choice))
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap();
                $field(limbs)
            }
        }

        impl core::cmp::PartialOrd for $field {
            fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl core::cmp::Ord for $field {
            fn cmp(&self, other: &Self) -> core::cmp::Ordering {
                let left = self.to_repr();
                let right = other.to_repr();
                left.iter()
                    .zip(right.iter())
                    .rev()
                    .find_map(|(left_byte, right_byte)| match left_byte.cmp(right_byte) {
                        core::cmp::Ordering::Equal => None,
                        res => Some(res),
                    })
                    .unwrap_or(core::cmp::Ordering::Equal)
            }
        }

        impl ff::Field for $field {
            const ZERO: Self = Self::zero();
            const ONE: Self = Self::one();

            fn random(mut rng: impl RngCore) -> Self {
                let mut wide = [0u8; Self::SIZE * 2];
                rng.fill_bytes(&mut wide);
                Self::from_uniform_bytes(&wide)
            }

            #[must_use]
            fn double(&self) -> Self {
                self.double()
            }

            #[inline(always)]
            #[must_use]
            fn square(&self) -> Self {
                self.square()
            }

            fn invert(&self) -> CtOption<Self> {
                self.invert()
            }

            fn sqrt(&self) -> CtOption<Self> {
                self.sqrt()
            }

            fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
                ff::helpers::sqrt_ratio_generic(num, div)
            }

            pow_vartime!($field_type);
        }

        impl $field {
            pub const SIZE: usize = NUM_LIMBS * 8;

            /// Returns zero, the additive identity.
            #[inline]
            pub const fn zero() -> $field {
                $field([0; NUM_LIMBS])
            }

            /// Returns one, the multiplicative identity.
            #[inline]
            pub const fn one() -> $field {
                R
            }

            /// Converts from an integer represented in little endian
            /// into its (congruent) `$field` representation.
            pub const fn from_raw(val: [u64; NUM_LIMBS]) -> Self {
                Self::montgomery_form(val, R2)
            }

            /// Attempts to convert a little-endian byte representation of
            /// a scalar into a `Fr`, failing if the input is not canonical.
            pub fn from_bytes(bytes: &[u8; Self::SIZE]) -> CtOption<$field> {
                let bytes: <Self as ff::PrimeField>::Repr = (*bytes).into();
                let z = <Self as ff::PrimeField>::from_repr(bytes);
                z
            }

            /// Converts an element of `Fr` into a byte representation in
            /// little-endian byte order.
            pub fn to_bytes(&self) -> [u8; Self::SIZE] {
                <Self as ff::PrimeField>::to_repr(self).into()
            }

            /// Returns the multiplicative inverse of the
            /// element. If it is zero, the method fails.
            #[inline(always)]
            pub fn invert(&self) -> CtOption<Self> {
                if let Some(inverse) = BYINVERTOR.invert(&self.0) {
                    CtOption::new(Self(inverse), Choice::from(1))
                } else {
                    CtOption::new(Self::zero(), Choice::from(0))
                }
            }

            // Returns the Jacobi symbol, where the numerator and denominator
            // are the element and the characteristic of the field, respectively.
            // The Jacobi symbol is applicable to odd moduli
            // while the Legendre symbol is applicable to prime moduli.
            // They are equivalent for prime moduli.
            #[inline(always)]
            pub fn jacobi(&self) -> i64 {
                $crate::ff_ext::jacobi::jacobi::<JACOBI_L>(&self.0, &MODULUS.0)
            }

            // Lexicographic comparison of Montgomery forms.
            #[inline(always)]
            fn is_less_than(x: &$field, y: &$field) -> bool {
                let borrow =
                    x.0.iter()
                        .zip(y.0.iter())
                        .fold(0, |borrow, (x, y)| sbb(*x, *y, borrow).1);

                borrow >> 63 == 1
            }
        }
    };
}

#[macro_export]
macro_rules! impl_serde_object {
    ($field:ident) => {
        impl SerdeObject for $field {
            fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
                debug_assert_eq!(bytes.len(), Self::SIZE);

                let inner = (0..NUM_LIMBS)
                    .map(|off| {
                        u64::from_le_bytes(bytes[off * 8..(off + 1) * 8].try_into().unwrap())
                    })
                    .collect::<Vec<_>>();
                Self(inner.try_into().unwrap())
            }
            fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
                if bytes.len() != Self::SIZE {
                    return None;
                }
                let elt = Self::from_raw_bytes_unchecked(bytes);
                Self::is_less_than(&elt, &MODULUS).then(|| elt)
            }
            fn to_raw_bytes(&self) -> Vec<u8> {
                let mut res = Vec::with_capacity(NUM_LIMBS * 4);
                for limb in self.0.iter() {
                    res.extend_from_slice(&limb.to_le_bytes());
                }
                res
            }
            fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
                let inner = [(); NUM_LIMBS].map(|_| {
                    let mut buf = [0; 8];
                    reader.read_exact(&mut buf).unwrap();
                    u64::from_le_bytes(buf)
                });
                Self(inner)
            }
            fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
                let mut inner = [0u64; NUM_LIMBS];
                for limb in inner.iter_mut() {
                    let mut buf = [0; 8];
                    reader.read_exact(&mut buf)?;
                    *limb = u64::from_le_bytes(buf);
                }
                let elt = Self(inner);
                Self::is_less_than(&elt, &MODULUS)
                    .then(|| elt)
                    .ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            "input number is not less than field modulus",
                        )
                    })
            }
            fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
                for limb in self.0.iter() {
                    writer.write_all(&limb.to_le_bytes())?;
                }
                Ok(())
            }
        }
    };
}

#[macro_export]
macro_rules! field_bits {
    // For #[cfg(target_pointer_width = "64")]
    ($field:ident) => {
        #[cfg(feature = "bits")]
        #[cfg_attr(docsrs, doc(cfg(feature = "bits")))]
        impl ff::PrimeFieldBits for $field {
            type ReprBits = [u64; NUM_LIMBS];

            fn to_le_bits(&self) -> ff::FieldBits<Self::ReprBits> {
                let bytes: [u8; Self::SIZE] = self.to_repr().into();

                let limbs = (0..NUM_LIMBS)
                    .map(|off| {
                        u64::from_le_bytes(bytes[off * 8..(off + 1) * 8].try_into().unwrap())
                    })
                    .collect::<Vec<_>>();

                ff::FieldBits::new(limbs.try_into().unwrap())
            }

            fn char_le_bits() -> ff::FieldBits<Self::ReprBits> {
                ff::FieldBits::new(MODULUS.0)
            }
        }
    };
    // For #[cfg(not(target_pointer_width = "64"))]
    ($field:ident) => {
        #[cfg(feature = "bits")]
        #[cfg_attr(docsrs, doc(cfg(feature = "bits")))]
        impl ff::PrimeFieldBits for $field {
            type ReprBits = [u32; NUM_LIMBS * 2];

            fn to_le_bits(&self) -> ff::FieldBits<Self::ReprBits> {
                let bytes = self.to_repr();

                let limbs = (0..NUM_LIMBS * 2)
                    .map(|off| {
                        u64::from_le_bytes(bytes[off * 4..(off + 1) * 4].try_into().unwrap())
                    })
                    .collect::<Vec<_>>();

                ff::FieldBits::new(limbs.try_into().unwrap())
            }

            fn char_le_bits() -> ff::FieldBits<Self::ReprBits> {
                ff::FieldBits::new(MODULUS_LIMBS_32)
            }
        }
    };
}

/// A macro to help define serialization and deserialization for prime field implementations
/// that use 32-byte representations. This assumes the concerned type implements PrimeField
/// (for from_repr, to_repr).
#[macro_export]
macro_rules! serialize_deserialize_primefield {
    ($field:ident, $repr:ty) => {
        impl ::serde::Serialize for $field {
            fn serialize<S: ::serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                let bytes = &self.to_repr();
                if serializer.is_human_readable() {
                    hex::serde::serialize(bytes, serializer)
                } else {
                    bytes.serialize(serializer)
                }
            }
        }

        use ::serde::de::Error as _;
        impl<'de> ::serde::Deserialize<'de> for $field {
            fn deserialize<D: ::serde::Deserializer<'de>>(
                deserializer: D,
            ) -> Result<Self, D::Error> {
                let bytes = if deserializer.is_human_readable() {
                    let bytes: [u8; $field::SIZE] = hex::serde::deserialize(deserializer)?;
                    bytes.into()
                } else {
                    <$repr>::deserialize(deserializer)?
                };
                Option::from(Self::from_repr(bytes)).ok_or_else(|| {
                    D::Error::custom("deserialized bytes don't encode a valid field element")
                })
            }
        }
    };
}
