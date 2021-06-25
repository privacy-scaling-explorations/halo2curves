use super::fq::Fq;
use super::fq2::Fq2;
use super::fr::Fr;
use crate::arithmetic::{BaseExt, Coordinates, CurveAffine, CurveExt, Group};
use crate::G1Prepared;
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use group::{
    cofactor::{CofactorCurve, CofactorGroup},
    prime::{PrimeCurve, PrimeCurveAffine, PrimeGroup},
    Curve as _, Group as _, GroupEncoding, UncompressedEncoding,
};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    Fq,
    Fr,
    32,
    "bn256_g1"
);

// pub trait UncompressedEncoding: Sized {
//     type Uncompressed: Default + AsRef<[u8]> + AsMut<[u8]>;

//     /// Attempts to deserialize an element from its uncompressed encoding.
//     fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self>;

//     /// Attempts to deserialize an uncompressed element, not checking if the element is in
//     /// the correct subgroup.
//     ///
//     /// **This is dangerous to call unless you trust the bytes you are reading; otherwise,
//     /// API invariants may be broken.** Please consider using
//     /// [`UncompressedEncoding::from_uncompressed`] instead.
//     fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self>;

//     /// Converts this element into its uncompressed encoding, so long as it's not
//     /// the point at infinity.
//     fn to_uncompressed(&self) -> Self::Uncompressed;
// }

impl G1 {
    const fn curve_constant_b() -> Fq {
        Fq::from_raw([3, 0, 0, 0])
    }

    pub fn generator() -> Self {
        const TWO: Fq = Fq::from_raw([2, 0, 0, 0]);

        Self {
            x: Fq::one(),
            y: TWO,
            z: Fq::one(),
        }
    }
}

impl G1Affine {
    fn random(mut rng: impl RngCore) -> Self {
        loop {
            let x = Fq::random(&mut rng);
            let ysign = (rng.next_u32() % 2) as u8;

            let x3 = x.square() * x;
            let y = (x3 + G1::curve_constant_b()).sqrt();
            if let Some(y) = Option::<Fq>::from(y) {
                let sign = y.to_bytes()[0] & 1;
                let y = if ysign ^ sign == 0 { y } else { -y };

                let p = Self {
                    x,
                    y,
                    infinity: Choice::from(0u8),
                };
                return p;
            }
        }
    }

    const fn curve_constant_b() -> Fq {
        Fq::from_raw([3, 0, 0, 0])
    }

    pub fn generator() -> Self {
        const TWO: Fq = Fq::from_raw([2, 0, 0, 0]);

        Self {
            x: Fq::one(),
            y: TWO,
            infinity: Choice::from(0u8),
        }
    }

    // pub fn prepare(&self) -> G1Prepared {
    //     G1Prepared::from_affine(*self)
    // }
}

pub struct G1Compressed([u8; 32]);

impl std::fmt::Debug for G1Compressed {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G1Compressed {
    fn default() -> Self {
        G1Compressed([0; 32])
    }
}

impl AsRef<[u8]> for G1Compressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G1Compressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

pub struct G1Uncompressed([u8; 96]);

impl std::fmt::Debug for G1Uncompressed {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G1Uncompressed {
    fn default() -> Self {
        G1Uncompressed([0; 96])
    }
}

impl AsRef<[u8]> for G1Uncompressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G1Uncompressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl UncompressedEncoding for G1Affine {
    type Uncompressed = G1Uncompressed;

    fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        unimplemented!();
    }
}

impl GroupEncoding for G1 {
    type Repr = G1Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        G1Affine::from_bytes(bytes).map(Self::from)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        G1Affine::from_bytes(bytes).map(Self::from)
    }

    fn to_bytes(&self) -> Self::Repr {
        G1Affine::from(self).to_bytes()
    }
}

impl GroupEncoding for G1Affine {
    type Repr = G1Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        let bytes = &bytes.0;
        let mut tmp = *bytes;
        let ysign = Choice::from(tmp[32 - 1] >> 7);
        tmp[32 - 1] &= 0b0111_1111;

        Fq::from_bytes(&tmp).and_then(|x| {
            CtOption::new(Self::identity(), x.ct_is_zero() & (!ysign)).or_else(|| {
                let x3 = x.square() * x;
                (x3 + G1::curve_constant_b()).sqrt().and_then(|y| {
                    let sign = Choice::from(y.to_bytes()[0] & 1);

                    let y = Fq::conditional_select(&y, &-y, ysign ^ sign);

                    CtOption::new(
                        G1Affine {
                            x,
                            y,
                            infinity: Choice::from(0u8),
                        },
                        Choice::from(1u8),
                    )
                })
            })
        })
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_bytes(bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        // TODO: not constant time
        if bool::from(self.is_identity()) {
            G1Compressed::default()
        } else {
            let (x, y) = (self.x, self.y);
            let sign = (y.to_bytes()[0] & 1) << 7;
            let mut xbytes = x.to_bytes();
            xbytes[32 - 1] |= sign;
            G1Compressed(xbytes)
        }
    }
}

new_curve_impl!(
    (pub),
    G2,
    G2Affine,
    Fq2,
    Fr,
    64,
    "bn256_g2"
);

pub struct G2Compressed([u8; 64]);

impl std::fmt::Debug for G2Compressed {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G2Compressed {
    fn default() -> Self {
        G2Compressed([0; 64])
    }
}

impl AsRef<[u8]> for G2Compressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G2Compressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

pub struct G2Uncompressed([u8; 128]);

impl std::fmt::Debug for G2Uncompressed {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0[..].fmt(f)
    }
}

impl Default for G2Uncompressed {
    fn default() -> Self {
        G2Uncompressed([0; 128])
    }
}

impl AsRef<[u8]> for G2Uncompressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for G2Uncompressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl UncompressedEncoding for G2Affine {
    type Uncompressed = G2Uncompressed;

    fn from_uncompressed(bytes: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn from_uncompressed_unchecked(bytes: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        unimplemented!();
    }
}

impl GroupEncoding for G2 {
    type Repr = G2Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        G2Affine::from_bytes(bytes).map(Self::from)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        G2Affine::from_bytes(bytes).map(Self::from)
    }

    fn to_bytes(&self) -> Self::Repr {
        G2Affine::from(self).to_bytes()
    }
}

impl GroupEncoding for G2Affine {
    type Repr = G2Compressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        let bytes = &bytes.0;
        let mut tmp = *bytes;
        let ysign = Choice::from(tmp[64 - 1] >> 7);
        tmp[64 - 1] &= 0b0111_1111;

        Fq2::from_bytes(&tmp).and_then(|x| {
            CtOption::new(Self::identity(), x.ct_is_zero() & (!ysign)).or_else(|| {
                let x3 = x.square() * x;
                (x3 + G2::curve_constant_b()).sqrt().and_then(|y| {
                    let sign = Choice::from(y.to_bytes()[0] & 1);

                    let y = Fq2::conditional_select(&y, &-y, ysign ^ sign);

                    CtOption::new(
                        G2Affine {
                            x,
                            y,
                            infinity: Choice::from(0u8),
                        },
                        Choice::from(1u8),
                    )
                })
            })
        })
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_bytes(bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        // TODO: not constant time
        if bool::from(self.is_identity()) {
            G2Compressed::default()
        } else {
            let (x, y) = (self.x, self.y);
            let sign = (y.to_bytes()[0] & 1) << 7;
            let mut xbytes = x.to_bytes();
            xbytes[64 - 1] |= sign;
            G2Compressed(xbytes)
        }
    }
}

impl G2 {
    const fn curve_constant_b() -> Fq2 {
        Fq2 {
            c0: Fq::from_raw([
                0x3267e6dc24a138e5,
                0xb5b4c5e559dbefa3,
                0x81be18991be06ac3,
                0x2b149d40ceb8aaae,
            ]),
            c1: Fq::from_raw([
                0xe4a2bd0685c315d2,
                0xa74fa084e52d1852,
                0xcd2cafadeed8fdf4,
                0x009713b03af0fed4,
            ]),
        }
    }

    pub fn generator() -> Self {
        let x = Fq2 {
            c0: Fq::from_raw([
                0x46debd5cd992f6ed,
                0x674322d4f75edadd,
                0x426a00665e5c4479,
                0x1800deef121f1e76,
            ]),
            c1: Fq::from_raw([
                0x97e485b7aef312c2,
                0xf1aa493335a9e712,
                0x7260bfb731fb5d25,
                0x198e9393920d483a,
            ]),
        };
        let y = Fq2 {
            c0: Fq::from_raw([
                0x4ce6cc0166fa7daa,
                0xe3d1e7690c43d37b,
                0x4aab71808dcb408f,
                0x12c85ea5db8c6deb,
            ]),

            c1: Fq::from_raw([
                0x55acdadcd122975b,
                0xbc4b313370b38ef3,
                0xec9e99ad690c3395,
                0x090689d0585ff075,
            ]),
        };

        Self {
            x,
            y,
            z: Fq2::one(),
        }
    }
}

impl G2Affine {
    fn random(mut rng: impl RngCore) -> Self {
        loop {
            let x = Fq2::random(&mut rng);
            let ysign = (rng.next_u32() % 2) as u8;

            let x3 = x.square() * x;
            let y = (x3 + G2::curve_constant_b()).sqrt();
            if let Some(y) = Option::<Fq2>::from(y) {
                let sign = y.to_bytes()[0] & 1;
                let y = if ysign ^ sign == 0 { y } else { -y };

                let p = Self {
                    x,
                    y,
                    infinity: Choice::from(0u8),
                };
                return p;
            }
        }
    }

    const fn curve_constant_b() -> Fq2 {
        G2::curve_constant_b()
    }

    pub fn generator() -> Self {
        let x = Fq2 {
            c0: Fq::from_raw([
                0x46debd5cd992f6ed,
                0x674322d4f75edadd,
                0x426a00665e5c4479,
                0x1800deef121f1e76,
            ]),
            c1: Fq::from_raw([
                0x97e485b7aef312c2,
                0xf1aa493335a9e712,
                0x7260bfb731fb5d25,
                0x198e9393920d483a,
            ]),
        };
        let y = Fq2 {
            c0: Fq::from_raw([
                0x4ce6cc0166fa7daa,
                0xe3d1e7690c43d37b,
                0x4aab71808dcb408f,
                0x12c85ea5db8c6deb,
            ]),

            c1: Fq::from_raw([
                0x55acdadcd122975b,
                0xbc4b313370b38ef3,
                0xec9e99ad690c3395,
                0x090689d0585ff075,
            ]),
        };

        Self {
            x,
            y,
            infinity: Choice::from(0u8),
        }
    }
}

#[cfg(test)]
use rand::SeedableRng;
#[cfg(test)]
use rand_xorshift::XorShiftRng;

#[test]
fn test_is_on_curve() {
    assert!(bool::from(G1Affine::identity().is_on_curve()));
    assert!(bool::from(G1Affine::generator().is_on_curve()));
    assert!(bool::from(G1::identity().is_on_curve()));
    assert!(bool::from(G1::generator().is_on_curve()));

    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let a = G1::random(&mut rng);
    assert!(bool::from(a.is_on_curve()));
}

#[test]
#[allow(clippy::eq_op)]
fn test_affine_point_equality() {
    let a = G1Affine::generator();
    let b = G1Affine::identity();

    assert!(a == a);
    assert!(b == b);
    assert!(a != b);
    assert!(b != a);
}

#[test]
#[allow(clippy::eq_op)]
fn test_projective_point_equality() {
    let a = G1::generator();
    let b = G1::identity();

    assert!(a == a);
    assert!(b == b);
    assert!(a != b);
    assert!(b != a);
}

#[test]
fn test_conditionally_select_affine() {
    let a = G1Affine::generator();
    let b = G1Affine::identity();

    assert_eq!(G1Affine::conditional_select(&a, &b, Choice::from(0u8)), a);
    assert_eq!(G1Affine::conditional_select(&a, &b, Choice::from(1u8)), b);
}

#[test]
fn test_conditionally_select_projective() {
    let a = G1::generator();
    let b = G1::identity();

    assert_eq!(G1::conditional_select(&a, &b, Choice::from(0u8)), a);
    assert_eq!(G1::conditional_select(&a, &b, Choice::from(1u8)), b);
}

#[test]
fn test_projective_to_affine() {
    let a = G1::generator();
    let b = G1::identity();

    assert!(bool::from(G1Affine::from(a).is_on_curve()));
    assert!(!bool::from(G1Affine::from(a).is_identity()));
    assert!(bool::from(G1Affine::from(b).is_on_curve()));
    assert!(bool::from(G1Affine::from(b).is_identity()));
}

#[test]
fn test_affine_to_projective() {
    let a = G1Affine::generator();
    let b = G1Affine::identity();

    assert!(bool::from(G1::from(a).is_on_curve()));
    assert!(!bool::from(G1::from(a).is_identity()));
    assert!(bool::from(G1::from(b).is_on_curve()));
    assert!(bool::from(G1::from(b).is_identity()));
}

#[test]
fn test_doubling() {
    {
        let tmp = G1::identity().double();
        assert!(bool::from(tmp.is_identity()));
        assert!(bool::from(tmp.is_on_curve()));
    }
    {
        let tmp = G1::generator().double();
        assert!(!bool::from(tmp.is_identity()));
        assert!(bool::from(tmp.is_on_curve()));
    }
}

#[test]
fn test_random() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let a = G1::random(&mut rng);

    assert!(bool::from(!a.is_identity()));
    assert!(bool::from(a.is_on_curve()));
}

#[test]
fn test_projective_addition() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let a = G1::identity();
    let b = G1::identity();
    let c = a + b;
    assert!(bool::from(c.is_identity()));
    assert!(bool::from(c.is_on_curve()));
    let c = a - b;
    assert!(bool::from(c.is_identity()));
    assert!(bool::from(c.is_on_curve()));

    let a = G1::identity();
    let a = -a;
    assert!(bool::from(a.is_on_curve()));
    assert!(bool::from(a.is_identity()));

    let a = G1::random(&mut rng);
    assert!(a == a + G1::identity());
    assert!(a == G1::identity() + a);
    assert!(-a == G1::identity() - a);

    let a = G1::identity();
    let a = a.double();
    assert!(bool::from(c.is_on_curve()));
    assert!(bool::from(a.is_identity()));

    let a = G1::random(&mut rng);
    assert!(a.double() - a == a);

    let a = G1::random(&mut rng);
    let b = G1::random(&mut rng);
    let c = G1::random(&mut rng);
    assert!(a + b == b + a);
    assert!(a - b == -(b - a));
    assert!(c + (a + b) == a + (c + b));
    assert!((a - b) - c == (a - c) - b);

    let a = G1::generator().double().double(); // 4P
    let b = G1::generator().double(); // 2P
    let c = a + b;

    let mut d = G1::generator();
    for _ in 0..5 {
        d += G1::generator();
    }

    assert!(c == d);
    assert!(!bool::from(c.is_identity()));
    assert!(bool::from(c.is_on_curve()));
    assert!(!bool::from(d.is_identity()));
    assert!(bool::from(d.is_on_curve()));
}

#[test]
fn test_mixed_addition() {
    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    let a = G1Affine::identity();
    let b = G1::identity();
    let c = a + b;
    assert!(bool::from(c.is_identity()));
    assert!(bool::from(c.is_on_curve()));
    let c = a - b;
    assert!(bool::from(c.is_identity()));
    assert!(bool::from(c.is_on_curve()));

    let a = G1::identity();
    let a = -a;
    assert!(bool::from(a.is_on_curve()));
    assert!(bool::from(a.is_identity()));
    let a = G1Affine::identity();
    let a = -a;
    assert!(bool::from(a.is_on_curve()));
    assert!(bool::from(a.is_identity()));

    let a = G1Affine::random(&mut rng);
    assert!(a.to_curve() == a + G1Affine::identity());

    let a = G1::random(&mut rng);
    assert!(a.double() - a == a);

    let a = G1Affine::random(&mut rng);
    let b = G1::random(&mut rng);
    let c = G1::random(&mut rng);
    assert!(a + b == b + a);
    assert!(a - b == -(b - a));
    assert!(c + (a + b) == a + (c + b));
    assert!((a - b) - c == (a - c) - b);
}

#[test]
fn test_projective_scalar_multiplication() {
    let g = G1::generator();
    let a = Fr::from_raw([
        0x2b56_8297_a56d_a71c,
        0xd8c3_9ecb_0ef3_75d1,
        0x435c_38da_67bf_bf96,
        0x8088_a050_26b6_59b2,
    ]);
    let b = Fr::from_raw([
        0x785f_dd9b_26ef_8b85,
        0xc997_f258_3769_5c18,
        0x4c8d_bc39_e7b7_56c1,
        0x70d9_b6cc_6d87_df20,
    ]);
    let c = a * b;

    assert_eq!((g * a) * b, g * c);
}

#[test]
fn test_affine_scalar_multiplication() {
    let g = G1Affine::generator();
    let a = Fr::from_raw([
        0x2b56_8297_a56d_a71c,
        0xd8c3_9ecb_0ef3_75d1,
        0x435c_38da_67bf_bf96,
        0x8088_a050_26b6_59b2,
    ]);
    let b = Fr::from_raw([
        0x785f_dd9b_26ef_8b85,
        0xc997_f258_3769_5c18,
        0x4c8d_bc39_e7b7_56c1,
        0x70d9_b6cc_6d87_df20,
    ]);
    let c = a * b;

    assert_eq!(G1Affine::from(g * a) * b, g * c);
}

#[test]
fn test_batch_normalize() {
    let a = G1::generator().double();
    let b = a.double();
    let c = b.double();

    for a_identity in (0..1).map(|n| n == 1) {
        for b_identity in (0..1).map(|n| n == 1) {
            for c_identity in (0..1).map(|n| n == 1) {
                let mut v = [a, b, c];
                if a_identity {
                    v[0] = G1::identity()
                }
                if b_identity {
                    v[1] = G1::identity()
                }
                if c_identity {
                    v[2] = G1::identity()
                }

                let mut t = [
                    G1Affine::identity(),
                    G1Affine::identity(),
                    G1Affine::identity(),
                ];
                let expected = [
                    G1Affine::from(v[0]),
                    G1Affine::from(v[1]),
                    G1Affine::from(v[2]),
                ];

                G1::batch_normalize(&v[..], &mut t[..]);

                assert_eq!(&t[..], &expected[..]);
            }
        }
    }
}

// #[test]
// fn test_is_torsion_free() {
//     let a = G1Affine {
//         x: Fq::from_raw_unchecked([
//             0x0aba_f895_b97e_43c8,
//             0xba4c_6432_eb9b_61b0,
//             0x1250_6f52_adfe_307f,
//             0x7502_8c34_3933_6b72,
//             0x8474_4f05_b8e9_bd71,
//             0x113d_554f_b095_54f7,
//         ]),
//         y: Fq::from_raw_unchecked([
//             0x73e9_0e88_f5cf_01c0,
//             0x3700_7b65_dd31_97e2,
//             0x5cf9_a199_2f0d_7c78,
//             0x4f83_c10b_9eb3_330d,
//             0xf6a6_3f6f_07f6_0961,
//             0x0c53_b5b9_7e63_4df3,
//         ]),
//         infinity: Choice::from(0u8),
//     };
//     assert!(!bool::from(a.is_torsion_free()));

//     assert!(bool::from(G1Affine::identity().is_torsion_free()));
//     assert!(bool::from(G1Affine::generator().is_torsion_free()));
// }

// #[test]
// fn test_mul_by_x() {
//     // multiplying by `x` a point in G1 is the same as multiplying by
//     // the equivalent scalar.
//     let generator = G1::generator();
//     let x = if crate::BLS_X_IS_NEGATIVE {
//         -Scalar::from(crate::BLS_X)
//     } else {
//         Scalar::from(crate::BLS_X)
//     };
//     assert_eq!(generator.mul_by_x(), generator * x);

//     let point = G1::generator() * Scalar::from(42);
//     assert_eq!(point.mul_by_x(), point * x);
// }

// #[test]
// fn test_clear_cofactor() {
//     // the generator (and the identity) are always on the curve,
//     // even after clearing the cofactor
//     let generator = G1::generator();
//     assert!(bool::from(generator.clear_cofactor().is_on_curve()));
//     let id = G1::identity();
//     assert!(bool::from(id.clear_cofactor().is_on_curve()));

//     let z = Fq::from_raw_unchecked([
//         0x3d2d1c670671394e,
//         0x0ee3a800a2f7c1ca,
//         0x270f4f21da2e5050,
//         0xe02840a53f1be768,
//         0x55debeb597512690,
//         0x08bd25353dc8f791,
//     ]);

//     let point = G1 {
//         x: Fq::from_raw_unchecked([
//             0x48af5ff540c817f0,
//             0xd73893acaf379d5a,
//             0xe6c43584e18e023c,
//             0x1eda39c30f188b3e,
//             0xf618c6d3ccc0f8d8,
//             0x0073542cd671e16c,
//         ]) * z,
//         y: Fq::from_raw_unchecked([
//             0x57bf8be79461d0ba,
//             0xfc61459cee3547c3,
//             0x0d23567df1ef147b,
//             0x0ee187bcce1d9b64,
//             0xb0c8cfbe9dc8fdc1,
//             0x1328661767ef368b,
//         ]),
//         z: z.square() * z,
//     };

//     assert!(bool::from(point.is_on_curve()));
//     assert!(!bool::from(G1Affine::from(point).is_torsion_free()));
//     let cleared_point = point.clear_cofactor();
//     assert!(bool::from(cleared_point.is_on_curve()));
//     assert!(bool::from(G1Affine::from(cleared_point).is_torsion_free()));

//     // in BLS12-381 the cofactor in G1 can be
//     // cleared multiplying by (1-x)
//     let h_eff = Scalar::from(1) + Scalar::from(crate::BLS_X);
//     assert_eq!(point.clear_cofactor(), point * h_eff);
// }
