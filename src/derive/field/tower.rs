#[macro_export]
macro_rules! impl_tower2_common {
    (
        $field:ident,
        $tower:ident,
        serde
    ) => {
        #[derive(Copy, Clone, Debug, Eq, PartialEq, Default)]
        #[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
        pub struct $tower {
            pub c0: $field,
            pub c1: $field,
        }

        impl_tower2_common!($field, $tower, implementations);
    };

    (
        $field:ident,
        $tower:ident

    ) => {
        #[derive(Copy, Clone, Debug, Eq, PartialEq, Default)]
        pub struct $tower {
            pub c0: $field,
            pub c1: $field,
        }

        impl_tower2_common!($field, $tower, implementations);
    };

    (
        $field:ident,
        $tower:ident,
        implementations
    ) => {
        impl ::core::ops::Neg for $tower {
            type Output = $tower;

            #[inline]
            fn neg(self) -> $tower {
                -&self
            }
        }

        impl ConditionallySelectable for $tower {
            fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
                $tower {
                    c0: $field::conditional_select(&a.c0, &b.c0, choice),
                    c1: $field::conditional_select(&a.c1, &b.c1, choice),
                }
            }
        }

        impl ConstantTimeEq for $tower {
            fn ct_eq(&self, other: &Self) -> Choice {
                self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
            }
        }

        impl $tower {
            #[inline]
            pub const fn zero() -> $tower {
                $tower {
                    c0: $field::zero(),
                    c1: $field::zero(),
                }
            }

            #[inline]
            pub const fn one() -> $tower {
                $tower {
                    c0: $field::one(),
                    c1: $field::zero(),
                }
            }

            #[inline]
            pub const fn new(c0: $field, c1: $field) -> Self {
                $tower { c0, c1 }
            }

            pub fn double(&self) -> Self {
                Self {
                    c0: self.c0.double(),
                    c1: self.c1.double(),
                }
            }

            pub fn add(&self, other: &Self) -> Self {
                Self {
                    c0: self.c0.add(&other.c0),
                    c1: self.c1.add(&other.c1),
                }
            }

            pub fn sub(&self, other: &Self) -> Self {
                Self {
                    c0: self.c0.sub(&other.c0),
                    c1: self.c1.sub(&other.c1),
                }
            }

            pub fn neg(&self) -> Self {
                Self {
                    c0: self.c0.neg(),
                    c1: self.c1.neg(),
                }
            }

            pub fn mul(&self, other: &Self) -> Self {
                let mut t = *other;
                t.mul_assign(self);
                t
            }

            pub fn square(&self) -> Self {
                let mut t = *self;
                t.square_assign();
                t
            }
        }
    };
}

#[macro_export]
macro_rules! impl_tower2_from_uniform_bytes {
    (
        $field:ident,
        $tower:ident,
        $size:expr
    ) => {
        impl FromUniformBytes<{ $size }> for $tower {
            fn from_uniform_bytes(bytes: &[u8; $size]) -> Self {
                assert!($size % 2 == 0);
                const SIZE: usize = $size / 2;
                let c0: [u8; SIZE] = bytes[SIZE..].try_into().unwrap();
                let c1: [u8; SIZE] = bytes[..SIZE].try_into().unwrap();
                Self::new(
                    $field::from_uniform_bytes(&c0),
                    $field::from_uniform_bytes(&c1),
                )
            }
        }
    };
}

#[macro_export]
macro_rules! impl_tower2 {
    (
        $field:ident,
        $tower:ident,
        $repr:ident
    ) => {
        /// `$tower` elements are ordered lexicographically.
        impl Ord for $tower {
            #[inline(always)]
            fn cmp(&self, other: &$tower) -> Ordering {
                match self.c1.cmp(&other.c1) {
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Less => Ordering::Less,
                    Ordering::Equal => self.c0.cmp(&other.c0),
                }
            }
        }

        impl PartialOrd for $tower {
            #[inline(always)]
            fn partial_cmp(&self, other: &$tower) -> Option<Ordering> {
                Some(self.cmp(other))
            }
        }

        impl From<$tower> for [u8; $tower::SIZE] {
            fn from(value: $tower) -> [u8; $tower::SIZE] {
                value.to_bytes()
            }
        }

        impl<'a> From<&'a $tower> for [u8; $tower::SIZE] {
            fn from(value: &'a $tower) -> [u8; $tower::SIZE] {
                value.to_bytes()
            }
        }

        impl $tower {
            pub const SIZE: usize = $field::SIZE * 2;

            /// Attempts to convert a little-endian byte representation of
            /// a scalar into a `$field`, failing if the input is not canonical.
            pub fn from_bytes(bytes: &[u8; $field::SIZE * 2]) -> CtOption<$tower> {
                let c0 = $field::from_bytes(bytes[0..$field::SIZE].try_into().unwrap());
                let c1 =
                    $field::from_bytes(bytes[$field::SIZE..$field::SIZE * 2].try_into().unwrap());
                CtOption::new(
                    $tower {
                        c0: c0.unwrap(),
                        c1: c1.unwrap(),
                    },
                    c0.is_some() & c1.is_some(),
                )
            }

            /// Converts an element of `$field` into a byte representation in
            /// little-endian byte order.
            #[allow(clippy::wrong_self_convention)]
            pub fn to_bytes(&self) -> [u8; $field::SIZE * 2] {
                let mut res = [0u8; $field::SIZE * 2];
                let c0_bytes = self.c0.to_bytes();
                let c1_bytes = self.c1.to_bytes();
                res[0..$field::SIZE].copy_from_slice(&c0_bytes[..]);
                res[$field::SIZE..$field::SIZE * 2].copy_from_slice(&c1_bytes[..]);
                res
            }

            #[inline]
            pub fn lexicographically_largest(&self) -> Choice {
                // If this element's c1 coefficient is lexicographically largest
                // then it is lexicographically largest. Otherwise, in the event
                // the c1 coefficient is zero and the c0 coefficient is
                // lexicographically largest, then this element is lexicographically
                // largest.

                self.c1.lexicographically_largest()
                    | (self.c1.is_zero() & self.c0.lexicographically_largest())
            }
        }

        impl ff::Field for $tower {
            const ZERO: Self = Self::zero();
            const ONE: Self = Self::one();

            fn random(mut rng: impl RngCore) -> Self {
                $tower {
                    c0: $field::random(&mut rng),
                    c1: $field::random(&mut rng),
                }
            }

            fn is_zero(&self) -> Choice {
                self.c0.is_zero() & self.c1.is_zero()
            }

            fn square(&self) -> Self {
                self.square()
            }

            fn double(&self) -> Self {
                self.double()
            }

            fn sqrt(&self) -> CtOption<Self> {
                self.sqrt()
            }

            fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
                ff::helpers::sqrt_ratio_generic(num, div)
            }

            fn invert(&self) -> CtOption<Self> {
                self.invert()
            }
        }

        impl PrimeField for $tower {
            type Repr = $repr;

            const MODULUS: &'static str = <$field as PrimeField>::MODULUS;
            const MULTIPLICATIVE_GENERATOR: Self = $tower {
                c0: $field::MULTIPLICATIVE_GENERATOR,
                c1: $field::ZERO,
            };
            const NUM_BITS: u32 = $field::NUM_BITS;
            const CAPACITY: u32 = $field::NUM_BITS;
            const S: u32 = $field::S;

            // TODO: Check that we can just 0 this and forget.
            const ROOT_OF_UNITY: Self = $tower::ZERO;
            const ROOT_OF_UNITY_INV: Self = $tower::ZERO;
            const DELTA: Self = $tower::ZERO;

            const TWO_INV: Self = $tower {
                c0: $field::TWO_INV,
                c1: $field::ZERO,
            };

            fn from_repr(repr: Self::Repr) -> CtOption<Self> {
                let c0: [u8; $field::SIZE] = repr.0[..$field::SIZE].try_into().unwrap();
                let c0: <$field as PrimeField>::Repr = c0.into();
                let c0 = $field::from_repr(c0);

                let c1: [u8; $field::SIZE] = repr.0[$field::SIZE..].try_into().unwrap();
                let c1: <$field as PrimeField>::Repr = c1.into();
                let c1 = $field::from_repr(c1);

                CtOption::new($tower::new(c0.unwrap(), c1.unwrap()), Choice::from(1))
            }

            fn to_repr(&self) -> Self::Repr {
                let mut res = [0u8; Self::SIZE];
                let c0 = self.c0.to_repr();
                let c1 = self.c1.to_repr();
                res[0..$field::SIZE].copy_from_slice(&c0.as_ref()[..]);
                res[$field::SIZE..$field::SIZE * 2].copy_from_slice(&c1.as_ref()[..]);
                $repr(res.into())
            }

            fn is_odd(&self) -> Choice {
                self.c0.is_odd() | (self.c0.is_zero() & self.c1.is_odd())
            }
        }

        impl WithSmallOrderMulGroup<3> for $tower {
            // $field::ZETA ^2
            const ZETA: Self = $tower {
                c0: ZETA,
                c1: $field::ZERO,
            };
        }

        impl From<u64> for $tower {
            fn from(val: u64) -> Self {
                $tower {
                    c0: $field::from(val),
                    c1: $field::ZERO,
                }
            }
        }

        impl Legendre for $tower {
            fn legendre(&self) -> i64 {
                self.norm().legendre()
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct $repr([u8; $field::SIZE * 2]);

        impl Default for $repr {
            fn default() -> Self {
                Self([0u8; $field::SIZE * 2])
            }
        }

        impl AsMut<[u8]> for $repr {
            fn as_mut(&mut self) -> &mut [u8] {
                &mut self.0
            }
        }

        impl AsRef<[u8]> for $repr {
            fn as_ref(&self) -> &[u8] {
                &self.0
            }
        }

        impl FromUniformBytes<{ $field::SIZE * 2 }> for $tower {
            fn from_uniform_bytes(bytes: &[u8; $field::SIZE * 2]) -> Self {
                Self::new($field::from_uniform_bytes(bytes), $field::ZERO)
            }
        }

        impl FromUniformBytes<{ $field::SIZE * 2 * 2 }> for $tower {
            fn from_uniform_bytes(bytes: &[u8; $field::SIZE * 2 * 2]) -> Self {
                let c0: [u8; $field::SIZE * 2] = bytes[..$field::SIZE * 2].try_into().unwrap();
                let c1: [u8; $field::SIZE * 2] = bytes[$field::SIZE * 2..].try_into().unwrap();
                Self::new(
                    $field::from_uniform_bytes(&c0),
                    $field::from_uniform_bytes(&c1),
                )
            }
        }

        impl $crate::serde::SerdeObject for $tower {
            fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
                debug_assert_eq!(bytes.len(), $field::SIZE * 2);
                let [c0, c1] = [0, $field::SIZE]
                    .map(|i| $field::from_raw_bytes_unchecked(&bytes[i..i + $field::SIZE]));
                Self { c0, c1 }
            }
            fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
                if bytes.len() != $field::SIZE * 2 {
                    return None;
                }
                let [c0, c1] =
                    [0, $field::SIZE].map(|i| $field::from_raw_bytes(&bytes[i..i + $field::SIZE]));
                c0.zip(c1).map(|(c0, c1)| Self { c0, c1 })
            }
            fn to_raw_bytes(&self) -> Vec<u8> {
                let mut res = Vec::with_capacity($field::SIZE * 2);
                for limb in self.c0.0.iter().chain(self.c1.0.iter()) {
                    res.extend_from_slice(&limb.to_le_bytes());
                }
                res
            }
            fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
                let [c0, c1] = [(); 2].map(|_| $field::read_raw_unchecked(reader));
                Self { c0, c1 }
            }
            fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
                let c0 = $field::read_raw(reader)?;
                let c1 = $field::read_raw(reader)?;
                Ok(Self { c0, c1 })
            }
            fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
                self.c0.write_raw(writer)?;
                self.c1.write_raw(writer)
            }
        }
    };
}

#[macro_export]
macro_rules! impl_tower6 {
    (
        $field:ident,
        $tower2:ident,
        $tower6:ident,
        $frobenius_c1:expr,
        $frobenius_c2:expr
    ) => {
        #[derive(Copy, Clone, Debug, Eq, PartialEq, Default)]
        pub struct $tower6 {
            pub c0: $tower2,
            pub c1: $tower2,
            pub c2: $tower2,
        }

        impl ::core::ops::Neg for $tower6 {
            type Output = $tower6;

            #[inline]
            fn neg(self) -> $tower6 {
                -&self
            }
        }

        impl subtle::ConditionallySelectable for $tower6 {
            fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
                $tower6 {
                    c0: $tower2::conditional_select(&a.c0, &b.c0, choice),
                    c1: $tower2::conditional_select(&a.c1, &b.c1, choice),
                    c2: $tower2::conditional_select(&a.c2, &b.c2, choice),
                }
            }
        }

        impl subtle::ConstantTimeEq for $tower6 {
            fn ct_eq(&self, other: &Self) -> subtle::Choice {
                self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1) & self.c2.ct_eq(&other.c2)
            }
        }

        impl $tower6 {
            #[inline]
            pub const fn zero() -> Self {
                $tower6 {
                    c0: $tower2::ZERO,
                    c1: $tower2::ZERO,
                    c2: $tower2::ZERO,
                }
            }

            #[inline]
            pub const fn one() -> Self {
                $tower6 {
                    c0: $tower2::ONE,
                    c1: $tower2::ZERO,
                    c2: $tower2::ZERO,
                }
            }

            #[inline]
            pub fn double(&self) -> Self {
                Self {
                    c0: self.c0.double(),
                    c1: self.c1.double(),
                    c2: self.c2.double(),
                }
            }

            #[inline]
            pub fn add(&self, other: &Self) -> Self {
                Self {
                    c0: self.c0 + other.c0,
                    c1: self.c1 + other.c1,
                    c2: self.c2 + other.c2,
                }
            }

            #[inline]
            pub fn sub(&self, other: &Self) -> Self {
                Self {
                    c0: self.c0 - other.c0,
                    c1: self.c1 - other.c1,
                    c2: self.c2 - other.c2,
                }
            }

            #[inline]
            pub fn neg(&self) -> Self {
                Self {
                    c0: -self.c0,
                    c1: -self.c1,
                    c2: -self.c2,
                }
            }

            pub fn mul(&self, other: &Self) -> Self {
                let mut t = *other;
                t.mul_assign(self);
                t
            }

            pub fn square(&self) -> Self {
                let mut t = *self;
                t.square_assign();
                t
            }

            pub fn random(mut rng: impl RngCore) -> Self {
                $tower6 {
                    c0: $tower2::random(&mut rng),
                    c1: $tower2::random(&mut rng),
                    c2: $tower2::random(&mut rng),
                }
            }
        }

        impl $tower6 {
            pub fn mul_assign(&mut self, other: &Self) {
                let a_a = self.c0 * other.c0;
                let b_b = self.c1 * other.c1;
                let c_c = self.c2 * other.c2;

                let t1 = (other.c1 + other.c2) * (self.c1 + self.c2) - (c_c + b_b);
                let t1 = t1.mul_by_nonresidue() + a_a;

                let t3 = (other.c0 + other.c2) * (self.c0 + self.c2) - (a_a - b_b + c_c);

                let t2 = (other.c0 + other.c1) * (self.c0 + self.c1) - (a_a + b_b);
                let t2 = t2 + c_c.mul_by_nonresidue();

                self.c0 = t1;
                self.c1 = t2;
                self.c2 = t3;
            }

            pub fn square_assign(&mut self) {
                let s0 = self.c0.square();
                let s1 = (self.c0 * &self.c1).double();
                let s2 = (self.c0 - self.c1 + &self.c2).square();
                let s3 = (self.c1 * self.c2).double();
                let s4 = self.c2.square();

                self.c0 = s3.mul_by_nonresidue() + s0;
                self.c1 = s4.mul_by_nonresidue() + s1;
                self.c2 = s1 + s2 + s3 - s0 - s4;
            }

            pub fn frobenius_map(&mut self, power: usize) {
                self.c0.frobenius_map(power);
                self.c1.frobenius_map(power);
                self.c2.frobenius_map(power);
                self.c1.mul_assign(&$frobenius_c1[power % 6]);
                self.c2.mul_assign(&$frobenius_c2[power % 6]);
            }

            pub fn mul_by_nonresidue(&self) -> Self {
                let c0 = self.c2.mul_by_nonresidue();
                let c1 = self.c0;
                let c2 = self.c1;
                Self { c0, c1, c2 }
            }

            pub fn mul_by_1(&self, c1: &$tower2) -> Self {
                let b_b = self.c1 * c1;

                let t1 = c1 * (self.c1 + self.c2) - b_b;
                let t1 = t1.mul_by_nonresidue();
                let t2 = c1 * (self.c0 + self.c1) - b_b;

                Self {
                    c0: t1,
                    c1: t2,
                    c2: b_b,
                }
            }

            pub fn mul_by_01(&self, c0: &$tower2, c1: &$tower2) -> Self {
                let a_a = self.c0 * c0;
                let b_b = self.c1 * c1;

                let t1 = *c1 * (self.c1 + self.c2) - b_b;
                let t1 = a_a + t1.mul_by_nonresidue();
                let t3 = *c0 * (self.c0 + self.c2) - a_a + b_b;
                let t2 = (c0 + c1) * (self.c0 + self.c1) - a_a - b_b;

                Self {
                    c0: t1,
                    c1: t2,
                    c2: t3,
                }
            }

            pub fn invert(&self) -> CtOption<Self> {
                let c0 = self.c2.mul_by_nonresidue() * self.c1.neg() + self.c0.square();
                let c1 = self.c2.square().mul_by_nonresidue() - (self.c0 * self.c1);
                let c2 = self.c1.square() - (self.c0 * self.c2);

                let t = (self.c2 * c1) + (self.c1 * c2);
                let t = t.mul_by_nonresidue() + (self.c0 * c0);

                t.invert().map(|t| $tower6 {
                    c0: t * c0,
                    c1: t * c1,
                    c2: t * c2,
                })
            }
        }

        impl Field for $tower6 {
            const ZERO: Self = Self::zero();
            const ONE: Self = Self::one();

            fn random(rng: impl RngCore) -> Self {
                $tower6::random(rng)
            }

            fn is_zero(&self) -> subtle::Choice {
                self.c0.is_zero() & self.c1.is_zero()
            }

            fn square(&self) -> Self {
                self.square()
            }

            fn double(&self) -> Self {
                self.double()
            }

            fn sqrt(&self) -> CtOption<Self> {
                unimplemented!()
            }

            fn sqrt_ratio(_num: &Self, _div: &Self) -> (subtle::Choice, Self) {
                unimplemented!()
            }

            fn invert(&self) -> CtOption<Self> {
                self.invert()
            }
        }
    };
}

#[macro_export]
macro_rules! impl_tower12 {
    (
        $field:ident,
        $tower2:ident,
        $tower6:ident,
        $tower12:ident,
        $frobenius:expr
    ) => {
        impl ff::Field for $tower12 {
            const ZERO: Self = Self::zero();
            const ONE: Self = Self::one();

            fn random(mut rng: impl rand_core::RngCore) -> Self {
                Self {
                    c0: $tower6::random(&mut rng),
                    c1: $tower6::random(&mut rng),
                }
            }

            fn is_zero(&self) -> Choice {
                self.c0.is_zero() & self.c1.is_zero()
            }

            fn square(&self) -> Self {
                self.square()
            }

            fn double(&self) -> Self {
                self.double()
            }

            fn sqrt(&self) -> CtOption<Self> {
                unimplemented!()
            }

            fn sqrt_ratio(_num: &Self, _div: &Self) -> (Choice, Self) {
                unimplemented!()
            }

            fn invert(&self) -> CtOption<Self> {
                self.invert()
            }
        }

        impl $tower12 {
            pub fn mul_assign(&mut self, other: &Self) {
                let t1 = self.c0 * other.c0;
                let t2 = self.c1 * other.c1;

                self.c1 = (self.c0 + self.c1) * (other.c0 + other.c1) - (t1 + t2);
                self.c0 = t1 + t2.mul_by_nonresidue();
            }

            pub fn square_assign(&mut self) {
                let ab = self.c0 * self.c1;
                let c0c1 = self.c0 + self.c1;
                let c0 = (self.c1.mul_by_nonresidue() + self.c0) * c0c1 - ab;

                self.c1 = ab.double();
                self.c0 = c0 - ab.mul_by_nonresidue();
            }

            pub fn mul_by_014(&mut self, c0: &$tower2, c1: &$tower2, c4: &$tower2) {
                let aa = self.c0.mul_by_01(c0, c1);
                let bb = self.c1.mul_by_1(c4);
                self.c1 = (self.c1 + &self.c0).mul_by_01(c0, &(c1 + c4)) - (aa + bb);
                self.c0 = bb.mul_by_nonresidue() + aa;
            }

            pub fn mul_by_034(&mut self, c0: &$tower2, c3: &$tower2, c4: &$tower2) {
                let t0 = $tower6 {
                    c0: self.c0.c0 * c0,
                    c1: self.c0.c1 * c0,
                    c2: self.c0.c2 * c0,
                };
                let t1 = self.c1.mul_by_01(c3, c4);
                self.c1 = (self.c0 + self.c1).mul_by_01(&(c0 + c3), c4) - t0 - t1;
                self.c0 = t0 + t1.mul_by_nonresidue();
            }

            pub fn invert(&self) -> CtOption<Self> {
                let c0 = self.c0.square();
                let c1 = self.c1.square().mul_by_nonresidue();
                (c0 - c1).invert().map(|t| $tower12 {
                    c0: self.c0 * t,
                    c1: self.c1 * -t,
                })
            }

            pub fn conjugate(&mut self) {
                self.c1 = -self.c1;
            }

            pub fn frobenius_map(&mut self, power: usize) {
                self.c0.frobenius_map(power);
                self.c1.frobenius_map(power);

                self.c1.c0.mul_assign(&$frobenius[power % 12]);
                self.c1.c1.mul_assign(&$frobenius[power % 12]);
                self.c1.c2.mul_assign(&$frobenius[power % 12]);
            }

            pub fn cyclotomic_square(&mut self) {
                fn fp4_square(c0: &mut $tower2, c1: &mut $tower2, a0: &$tower2, a1: &$tower2) {
                    let t0 = a0.square();
                    let t1 = a1.square();
                    *c0 = t1.mul_by_nonresidue() + t0;
                    *c1 = (a0 + a1).square() - t0 - t1;
                }

                let mut t3 = $tower2::zero();
                let mut t4 = $tower2::zero();
                let mut t5 = $tower2::zero();
                let mut t6 = $tower2::zero();
                fp4_square(&mut t3, &mut t4, &self.c0.c0, &self.c1.c1);

                self.c0.c0 = (t3 - self.c0.c0).double() + t3;
                self.c1.c1 = (t4 + self.c1.c1).double() + t4;

                fp4_square(&mut t3, &mut t4, &self.c1.c0, &self.c0.c2);
                fp4_square(&mut t5, &mut t6, &self.c0.c1, &self.c1.c2);

                self.c0.c1 = (t3 - self.c0.c1).double() + t3;
                self.c1.c2 = (t4 + self.c1.c2).double() + t4;

                let t3 = t6.mul_by_nonresidue();
                self.c1.c0 = (t3 + self.c1.c0).double() + t3;
                self.c0.c2 = (t5 - self.c0.c2).double() + t5;
            }
        }
    };
}
