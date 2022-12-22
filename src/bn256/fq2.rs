use super::fq::{Fq, NEGATIVE_ONE};
use super::LegendreSymbol;
use core::convert::TryInto;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use pasta_curves::arithmetic::{FieldExt, Group, SqrtRatio};
use rand::RngCore;
use std::cmp::Ordering;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// An element of Fq2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}

/// `Fq2` elements are ordered lexicographically.
impl Ord for Fq2 {
    #[inline(always)]
    fn cmp(&self, other: &Fq2) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl PartialOrd for Fq2 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fq2) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl ConditionallySelectable for Fq2 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fq2 {
            c0: Fq::conditional_select(&a.c0, &b.c0, choice),
            c1: Fq::conditional_select(&a.c1, &b.c1, choice),
        }
    }
}

impl ConstantTimeEq for Fq2 {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.c0.ct_eq(&other.c0) & self.c1.ct_eq(&other.c1)
    }
}

impl Default for Fq2 {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl From<Fq2> for [u8; 64] {
    fn from(value: Fq2) -> [u8; 64] {
        value.to_bytes()
    }
}

impl<'a> From<&'a Fq2> for [u8; 64] {
    fn from(value: &'a Fq2) -> [u8; 64] {
        value.to_bytes()
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    #[inline]
    fn neg(self) -> Fq2 {
        -&self
    }
}

impl<'a> Neg for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn neg(self) -> Fq2 {
        self.neg()
    }
}

impl<'a, 'b> Sub<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn sub(self, rhs: &'b Fq2) -> Fq2 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn add(self, rhs: &'b Fq2) -> Fq2 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn mul(self, rhs: &'b Fq2) -> Fq2 {
        self.mul(rhs)
    }
}

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
};
impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);

impl Fq2 {
    pub const fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 { c0, c1 }
    }

    pub const fn size() -> usize {
        64
    }
    /// Attempts to convert a little-endian byte representation of
    /// a scalar into a `Fq`, failing if the input is not canonical.
    pub fn from_bytes(bytes: &[u8; 64]) -> CtOption<Fq2> {
        let c0 = Fq::from_bytes(bytes[0..32].try_into().unwrap());
        let c1 = Fq::from_bytes(bytes[32..64].try_into().unwrap());
        CtOption::new(
            Fq2 {
                c0: c0.unwrap(),
                c1: c1.unwrap(),
            },
            c0.is_some() & c1.is_some(),
        )
    }

    /// Converts an element of `Fq` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes(&self) -> [u8; 64] {
        let mut res = [0u8; 64];
        let c0_bytes = self.c0.to_bytes();
        let c1_bytes = self.c1.to_bytes();
        res[0..32].copy_from_slice(&c0_bytes[..]);
        res[32..64].copy_from_slice(&c1_bytes[..]);
        res
    }

    pub fn legendre(&self) -> LegendreSymbol {
        self.norm().legendre()
    }

    pub fn mul_assign(&mut self, other: &Self) {
        let mut t1 = self.c0 * other.c0;
        let mut t0 = self.c0 + self.c1;
        let t2 = self.c1 * other.c1;
        self.c1 = other.c0 + other.c1;
        self.c0 = t1 - t2;
        t1 += t2;
        t0 *= self.c1;
        self.c1 = t0 - t1;
    }

    pub fn square_assign(&mut self) {
        let ab = self.c0 * self.c1;
        let c0c1 = self.c0 + self.c1;
        let mut c0 = -self.c1;
        c0 += self.c0;
        c0 *= c0c1;
        c0 -= ab;
        self.c1 = ab.double();
        self.c0 = c0 + ab;
    }

    pub fn double(&self) -> Self {
        Self {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }

    pub fn double_assign(&mut self) {
        self.c0 = self.c0.double();
        self.c1 = self.c1.double();
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

    pub fn neg(&self) -> Self {
        Self {
            c0: self.c0.neg(),
            c1: self.c1.neg(),
        }
    }

    // conjucate by negating c1
    pub fn conjugate(&mut self) {
        self.c1 = -self.c1;
    }

    pub fn frobenius_map(&mut self, power: usize) {
        self.c1 *= &FROBENIUS_COEFF_FQ2_C1[power % 2];
    }

    /// Multiply this element by quadratic nonresidue 9 + u.
    pub fn mul_by_nonresidue(&mut self) {
        // (xi+y)(i+9) = (9x+y)i+(9y-x)
        let t0 = self.c0;
        let t1 = self.c1;

        // 8*x*i + 8*y
        self.double_assign();
        self.double_assign();
        self.double_assign();

        // 9*y
        self.c0 += &t0;
        // (9*y - x)
        self.c0 -= &t1;

        // (9*x)i
        self.c1 += &t1;
        // (9*x + y)
        self.c1 += &t0;
    }

    // Multiply this element by ξ where ξ=i+9
    pub fn mul_by_xi(&mut self) {
        // (xi+y)(i+9) = (9x+y)i+(9y-x)
        let t0 = self.c0;
        let t1 = self.c1;

        // 8*x*i + 8*y
        self.double_assign();
        self.double_assign();
        self.double_assign();

        // 9*y
        self.c0 += &t0;
        // (9*y - x)
        self.c0 -= &t1;

        // (9*x)i
        self.c1 += &t1;
        // (9*x + y)
        self.c1 += &t0;
    }

    /// Norm of Fq2 as extension field in i over Fq
    pub fn norm(&self) -> Fq {
        let mut t0 = self.c0;
        let mut t1 = self.c1;
        t0 = t0.square();
        t1 = t1.square();
        t1 + t0
    }

    pub fn invert(&self) -> CtOption<Self> {
        let mut t1 = self.c1;
        t1 = t1.square();
        let mut t0 = self.c0;
        t0 = t0.square();
        t0 += &t1;
        t0.invert().map(|t| {
            let mut tmp = Fq2 {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0 *= &t;
            tmp.c1 *= &t;
            tmp.c1 = -tmp.c1;

            tmp
        })
    }
}

impl Field for Fq2 {
    fn random(mut rng: impl RngCore) -> Self {
        Fq2 {
            c0: Fq::random(&mut rng),
            c1: Fq::random(&mut rng),
        }
    }

    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
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
        // Algorithm 9, https://eprint.iacr.org/2012/685.pdf

        if self.is_zero().into() {
            CtOption::new(Self::zero(), Choice::from(1))
        } else {
            // a1 = self^((q - 3) / 4)
            // 0xc19139cb84c680a6e14116da060561765e05aa45a1c72a34f082305b61f3f51
            let u: [u64; 4] = [
                0x4f082305b61f3f51,
                0x65e05aa45a1c72a3,
                0x6e14116da0605617,
                0x0c19139cb84c680a,
            ];
            let mut a1 = self.pow(&u);
            let mut alpha = a1;

            alpha.square_assign();
            alpha.mul_assign(self);
            let mut a0 = alpha;
            a0.frobenius_map(1);
            a0.mul_assign(&alpha);

            let neg1 = Fq2 {
                c0: NEGATIVE_ONE,
                c1: Fq::zero(),
            };

            if a0 == neg1 {
                CtOption::new(a0, Choice::from(0))
            } else {
                a1.mul_assign(self);

                if alpha == neg1 {
                    a1.mul_assign(&Fq2 {
                        c0: Fq::zero(),
                        c1: Fq::one(),
                    });
                } else {
                    alpha += &Fq2::one();
                    // alpha = alpha^((q - 1) / 2)
                    // 0x183227397098d014dc2822db40c0ac2ecbc0b548b438e5469e10460b6c3e7ea3
                    let u: [u64; 4] = [
                        0x9e10460b6c3e7ea3,
                        0xcbc0b548b438e546,
                        0xdc2822db40c0ac2e,
                        0x183227397098d014,
                    ];
                    alpha = alpha.pow(&u);
                    a1.mul_assign(&alpha);
                }
                CtOption::new(a1, Choice::from(1))
            }
        }
    }

    fn invert(&self) -> CtOption<Self> {
        self.invert()
    }
}

impl From<bool> for Fq2 {
    fn from(bit: bool) -> Fq2 {
        if bit {
            Fq2::one()
        } else {
            Fq2::zero()
        }
    }
}

impl From<u64> for Fq2 {
    fn from(val: u64) -> Self {
        Fq2 {
            c0: Fq::from(val),
            c1: Fq::zero(),
        }
    }
}

impl FieldExt for Fq2 {
    const MODULUS: &'static str =
        "0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47";

    const ROOT_OF_UNITY_INV: Self = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };
    const DELTA: Self = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };
    const TWO_INV: Self = Fq2 {
        c0: Fq::from_raw([
            0x9e10460b6c3e7ea4,
            0xcbc0b548b438e546,
            0xdc2822db40c0ac2e,
            0x183227397098d014,
        ]),
        c1: Fq([0, 0, 0, 0]),
    };
    const ZETA: Self = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };

    /// Converts a 512-bit little endian integer into
    /// a `Fq` by reducing by the modulus.
    fn from_bytes_wide(bytes: &[u8; 64]) -> Self {
        Self::new(Fq::from_bytes_wide(bytes), Fq::zero())
    }

    fn from_u128(v: u128) -> Self {
        Fq2 {
            c0: Fq::from_raw([v as u64, (v >> 64) as u64, 0, 0]),
            c1: Fq::zero(),
        }
    }

    fn get_lower_128(&self) -> u128 {
        self.c0.get_lower_128()
    }

    // /// Writes this element in its normalized, little endian form into a buffer.
    // fn write<W: Write>(&self, writer: &mut W) -> io::Result<()> {
    //     let compressed = self.to_bytes();
    //     writer.write_all(&compressed[..])
    // }

    // /// Reads a normalized, little endian represented field element from a
    // /// buffer.
    // fn read<R: Read>(reader: &mut R) -> io::Result<Self> {
    //     let mut compressed = [0u8; 64];
    //     reader.read_exact(&mut compressed[..])?;
    //     Option::from(Self::from_bytes(&compressed))
    //         .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "invalid point encoding in proof"))
    // }
}

impl SqrtRatio for Fq2 {
    const T_MINUS1_OVER2: [u64; 4] = [0, 0, 0, 0];

    fn pow_by_t_minus1_over2(&self) -> Self {
        unimplemented!();
    }

    fn get_lower_32(&self) -> u32 {
        unimplemented!();
    }

    #[cfg(feature = "sqrt-table")]
    fn sqrt_ratio(num: &Self, div: &Self) -> (Choice, Self) {
        unimplemented!();
    }

    #[cfg(feature = "sqrt-table")]
    fn sqrt_alt(&self) -> (Choice, Self) {
        unimplemented!();
    }
}

impl Group for Fq2 {
    type Scalar = Fq2;

    fn group_zero() -> Self {
        Self::zero()
    }
    fn group_add(&mut self, rhs: &Self) {
        *self += *rhs;
    }
    fn group_sub(&mut self, rhs: &Self) {
        *self -= *rhs;
    }
    fn group_scale(&mut self, by: &Self::Scalar) {
        *self *= *by;
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Fq2Bytes([u8; 64]);

impl Default for Fq2Bytes {
    fn default() -> Self {
        Self([0u8; 64])
    }
}

impl AsMut<[u8]> for Fq2Bytes {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl AsRef<[u8]> for Fq2Bytes {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl ff::PrimeField for Fq2 {
    type Repr = Fq2Bytes;

    const NUM_BITS: u32 = 254;
    const CAPACITY: u32 = 253;

    const S: u32 = 0;

    fn from_repr(repr: Self::Repr) -> CtOption<Self> {
        let c0 = Fq::from_bytes(&repr.0[..32].try_into().unwrap());
        let c1 = Fq::from_bytes(&repr.0[32..].try_into().unwrap());
        // Disallow overflow representation
        CtOption::new(Fq2::new(c0.unwrap(), c1.unwrap()), Choice::from(1))
    }

    fn to_repr(&self) -> Self::Repr {
        Fq2Bytes(self.to_bytes())
    }

    fn is_odd(&self) -> Choice {
        Choice::from(self.to_repr().as_ref()[0] & 1)
    }

    fn multiplicative_generator() -> Self {
        unimplemented!()
    }

    fn root_of_unity() -> Self {
        unimplemented!()
    }
}

impl crate::serde::SerdeObject for Fq2 {
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        debug_assert_eq!(bytes.len(), 64);
        let [c0, c1] = [0, 32].map(|i| Fq::from_raw_bytes_unchecked(&bytes[i..i + 32]));
        Self { c0, c1 }
    }
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 64 {
            return None;
        }
        let [c0, c1] = [0, 32].map(|i| Fq::from_raw_bytes(&bytes[i..i + 32]));
        c0.zip(c1).map(|(c0, c1)| Self { c0, c1 })
    }
    fn to_raw_bytes(&self) -> Vec<u8> {
        let mut res = Vec::with_capacity(64);
        for limb in self.c0.0.iter().chain(self.c1.0.iter()) {
            res.extend_from_slice(&limb.to_le_bytes());
        }
        res
    }
    fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
        let [c0, c1] = [(); 2].map(|_| Fq::read_raw_unchecked(reader));
        Self { c0, c1 }
    }
    fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let c0 = Fq::read_raw(reader)?;
        let c1 = Fq::read_raw(reader)?;
        Ok(Self { c0, c1 })
    }
    fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.c0.write_raw(writer)?;
        self.c1.write_raw(writer)
    }
}

pub const FROBENIUS_COEFF_FQ2_C1: [Fq; 2] = [
    // Fq(-1)**(((q^0) - 1) / 2)
    // it's 1 in Montgommery form
    Fq([
        0xd35d438dc58f0d9d,
        0x0a78eb28f5c70b3d,
        0x666ea36f7879462c,
        0x0e0a77c19a07df2f,
    ]),
    // Fq(-1)**(((q^1) - 1) / 2)
    Fq([
        0x68c3488912edefaa,
        0x8d087f6872aabf4f,
        0x51e1a24709081231,
        0x2259d6b14729c0fa,
    ]),
];

#[cfg(test)]
use rand::SeedableRng;
#[cfg(test)]
use rand_xorshift::XorShiftRng;

#[test]
fn test_ser() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let a0 = Fq2::random(&mut rng);
    let a_bytes = a0.to_bytes();
    let a1 = Fq2::from_bytes(&a_bytes).unwrap();
    assert_eq!(a0, a1);
}

#[test]
fn test_fq2_ordering() {
    let mut a = Fq2 {
        c0: Fq::zero(),
        c1: Fq::zero(),
    };

    let mut b = a;

    assert!(a.cmp(&b) == Ordering::Equal);
    b.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Equal);
    b.c1 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Less);
    a.c1 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Greater);
    b.c0 += &Fq::one();
    assert!(a.cmp(&b) == Ordering::Equal);
}

#[test]
fn test_fq2_basics() {
    assert_eq!(
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        },
        Fq2::zero()
    );
    assert_eq!(
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        },
        Fq2::one()
    );
    assert_eq!(Fq2::zero().is_zero().unwrap_u8(), 1);
    assert_eq!(Fq2::one().is_zero().unwrap_u8(), 0);
    assert_eq!(
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::one(),
        }
        .is_zero()
        .unwrap_u8(),
        0
    );
}

#[test]
fn test_fq2_squaring() {
    let mut a = Fq2 {
        c0: Fq::one(),
        c1: Fq::one(),
    }; // u + 1
    a.square_assign();
    assert_eq!(
        a,
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::one() + Fq::one(),
        }
    ); // 2u

    let mut a = Fq2 {
        c0: Fq::zero(),
        c1: Fq::one(),
    }; // u
    a.square_assign();
    assert_eq!(a, {
        let neg1 = -Fq::one();
        Fq2 {
            c0: neg1,
            c1: Fq::zero(),
        }
    }); // -1
}

#[test]
fn test_fq2_mul_nonresidue() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);
    let nine = Fq::one().double().double().double() + Fq::one();
    let nqr = Fq2 {
        c0: nine,
        c1: Fq::one(),
    };

    for _ in 0..1000 {
        let mut a = Fq2::random(&mut rng);
        let mut b = a;
        a.mul_by_nonresidue();
        b.mul_assign(&nqr);

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq2_legendre() {
    assert_eq!(LegendreSymbol::Zero, Fq2::zero().legendre());
    // i^2 = -1
    let mut m1 = Fq2::one();
    m1 = m1.neg();
    assert_eq!(LegendreSymbol::QuadraticResidue, m1.legendre());
    m1.mul_by_nonresidue();
    assert_eq!(LegendreSymbol::QuadraticNonResidue, m1.legendre());
}

#[test]
pub fn test_sqrt() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..10000 {
        let a = Fq2::random(&mut rng);
        if a.legendre() == LegendreSymbol::QuadraticNonResidue {
            assert!(bool::from(a.sqrt().is_none()));
        }
    }

    for _ in 0..10000 {
        let a = Fq2::random(&mut rng);
        let mut b = a;
        b.square_assign();
        assert_eq!(b.legendre(), LegendreSymbol::QuadraticResidue);

        let b = b.sqrt().unwrap();
        let mut negb = b;
        negb = negb.neg();

        assert!(a == b || a == negb);
    }

    let mut c = Fq2::one();
    for _ in 0..10000 {
        let mut b = c;
        b.square_assign();
        assert_eq!(b.legendre(), LegendreSymbol::QuadraticResidue);

        b = b.sqrt().unwrap();

        if b != c {
            b = b.neg();
        }

        assert_eq!(b, c);

        c += &Fq2::one();
    }
}

#[test]
fn test_frobenius() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    for _ in 0..100 {
        for i in 0..14 {
            let mut a = Fq2::random(&mut rng);
            let mut b = a;

            for _ in 0..i {
                a = a.pow(&[
                    0x3c208c16d87cfd47,
                    0x97816a916871ca8d,
                    0xb85045b68181585d,
                    0x30644e72e131a029,
                ]);
            }
            b.frobenius_map(i);

            assert_eq!(a, b);
        }
    }
}

#[test]
fn test_field() {
    crate::tests::field::random_field_tests::<Fq2>("fq2".to_string());
}

#[test]
fn test_serialization() {
    crate::tests::field::random_serialization_test::<Fq2>("fq2".to_string());
}
