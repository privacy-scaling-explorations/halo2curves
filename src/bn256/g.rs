use crate::arithmetic::{
    BaseExt, Coordinates, CurveAffine, CurveExt, FieldExt, Group, LinearCombinationEngine,
};
use crate::bn256::Fq;
use crate::bn256::Fq2;
use crate::bn256::Fr;
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use group::{
    cofactor::CofactorGroup, prime::PrimeCurveAffine, Curve as _, Group as _, GroupEncoding,
};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

new_curve_impl!(
    (pub),
    G1,
    G1Affine,
    G1Compressed,
    Fq,
    Fr,
    (G1_GENERATOR_X,G1_GENERATOR_Y),
    G1_B,
    "bn256_g1"
);

new_curve_impl!(
    (pub),
    G2,
    G2Affine,
    G2Compressed,
    Fq2,
    Fr,
    (G2_GENERATOR_X, G2_GENERATOR_Y),
    G2_B,
    "bn256_g2"
);

const G1_GENERATOR_X: Fq = Fq::one();
const G1_GENERATOR_Y: Fq = Fq::from_raw([2, 0, 0, 0]);
const G1_B: Fq = Fq::from_raw([3, 0, 0, 0]);

impl group::cofactor::CofactorGroup for G1 {
    type Subgroup = G1;

    fn clear_cofactor(&self) -> Self {
        *self
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, 1.into())
    }

    fn is_torsion_free(&self) -> Choice {
        1.into()
    }
}

const G2_B: Fq2 = Fq2 {
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
};

const G2_GENERATOR_X: Fq2 = Fq2 {
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

const G2_GENERATOR_Y: Fq2 = Fq2 {
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

impl CofactorGroup for G2 {
    type Subgroup = G2;

    fn clear_cofactor(&self) -> Self {
        // "0x30644e72e131a029b85045b68181585e06ceecda572a2489345f2299c0f9fa8d"
        let e: [u8; 32] = [
            0x30, 0x64, 0x4e, 0x72, 0xe1, 0x31, 0xa0, 0x29, 0xb8, 0x50, 0x45, 0xb6, 0x81, 0x81,
            0x58, 0x5e, 0x06, 0xce, 0xec, 0xda, 0x57, 0x2a, 0x24, 0x89, 0x34, 0x5f, 0x22, 0x99,
            0xc0, 0xf9, 0xfa, 0x8d,
        ];

        // self * COFACTOR_G2
        let mut acc = G2::identity();
        for bit in e
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = G2::conditional_select(&acc, &(acc + self), bit);
        }
        acc
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        unimplemented!();
    }

    fn is_torsion_free(&self) -> Choice {
        // "0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001"
        let e: [u8; 32] = [
            0x30, 0x64, 0x4e, 0x72, 0xe1, 0x31, 0xa0, 0x29, 0xb8, 0x50, 0x45, 0xb6, 0x81, 0x81,
            0x58, 0x5d, 0x28, 0x33, 0xe8, 0x48, 0x79, 0xb9, 0x70, 0x91, 0x43, 0xe1, 0xf5, 0x93,
            0xf0, 0x00, 0x00, 0x01,
        ];

        // self * GROUP_ORDER;

        let mut acc = G2::identity();
        for bit in e
            .iter()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(1)
        {
            acc = acc.double();
            acc = G2::conditional_select(&acc, &(acc + self), bit);
        }
        acc.is_identity()
    }
}

impl G2 {
    pub fn random(mut rng: impl RngCore) -> Self {
        let mut point = <Self as group::Group>::random(&mut rng);
        point = point.clear_cofactor();
        point
    }
}

pub struct PointCombination {
    sum: G1,
    base: Fr,
    lc: Fr,
}

impl LinearCombinationEngine for PointCombination {
    type Lhs = G1;
    type Rhs = Fr;

    fn new(base: Self::Rhs) -> Self {
        PointCombination {
            sum: G1::identity(),
            base,
            lc: Fr::one(),
        }
    }

    fn result(&mut self) -> Self::Lhs {
        let res = self.sum;
        self.sum = G1::identity();
        self.lc = Fr::one();
        res
    }

    fn add(&mut self, coeff: Self::Lhs) {
        self.sum += coeff * self.lc;
        self.lc *= self.base;
    }

    fn add_with_aux(&mut self, elem: Self::Lhs, aux: Self::Rhs) {
        self.sum += elem * self.lc * aux;
        self.lc *= self.base;
    }

    fn combine(base: Self::Rhs, coeffs: Vec<Self::Lhs>) -> Self::Lhs {
        let mut lc = Self::new(base);
        for coeff in coeffs.iter() {
            lc.add(*coeff);
        }
        lc.result()
    }
}

#[cfg(test)]
mod tests {

    use crate::bn256::{Fr, G1Affine, PointCombination, G1, G2};
    use ff::Field;

    use crate::arithmetic::{CurveAffine, CurveExt, FieldExt};
    use group::{cofactor::CofactorGroup, prime::PrimeCurveAffine};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    fn is_on_curve<G: CurveExt>() {
        assert!(bool::from(G::identity().is_on_curve()));
        assert!(bool::from(G::generator().is_on_curve()));
        assert!(bool::from(G::identity().is_on_curve()));
        assert!(bool::from(G::generator().is_on_curve()));

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..100 {
            let point = G::random(&mut rng);
            assert!(bool::from(point.is_on_curve()));
            let affine_point: G::AffineExt = point.into();
            assert!(bool::from(affine_point.is_on_curve()));
        }
    }

    fn equality<G: CurveExt>() {
        let a = G::generator();
        let b = G::identity();

        assert!(a == a);
        assert!(b == b);
        assert!(a != b);
        assert!(b != a);

        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        for _ in 0..100 {
            let a = G::random(&mut rng);
            let b = G::random(&mut rng);

            assert!(a == a);
            assert!(b == b);
            assert!(a != b);
            assert!(b != a);

            let a: G::AffineExt = a.into();
            let b: G::AffineExt = b.into();

            assert!(a == a);
            assert!(b == b);
            assert!(a != b);
            assert!(b != a);
        }
    }

    fn projective_to_affine_affine_to_projective<G: CurveExt>() {
        let a = G::generator();
        let b = G::identity();

        assert!(bool::from(G::AffineExt::from(a).is_on_curve()));
        assert!(!bool::from(G::AffineExt::from(a).is_identity()));
        assert!(bool::from(G::AffineExt::from(b).is_on_curve()));
        assert!(bool::from(G::AffineExt::from(b).is_identity()));

        let a = G::AffineExt::generator();
        let b = G::AffineExt::identity();

        assert!(bool::from(G::from(a).is_on_curve()));
        assert!(!bool::from(G::from(a).is_identity()));
        assert!(bool::from(G::from(b).is_on_curve()));
        assert!(bool::from(G::from(b).is_identity()));
    }

    fn projective_addition<G: CurveExt>() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let a = G::identity();
        let b = G::identity();
        let c = a + b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        let c = a - b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));

        let a = G::identity();
        let a = -a;
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(a.is_identity()));

        let a = G::random(&mut rng);
        assert!(a == a + G::identity());
        assert!(a == G::identity() + a);
        assert!(-a == G::identity() - a);

        let a = G::identity();
        let a = a.double();
        assert!(bool::from(c.is_on_curve()));
        assert!(bool::from(a.is_identity()));

        let a = G::generator();
        let a = a.double();
        assert!(bool::from(c.is_on_curve()));
        assert_eq!(a, G::generator() + G::generator());

        let a = G::random(&mut rng);
        assert!(a.double() - a == a);

        let a = G::random(&mut rng);
        let b = G::random(&mut rng);
        let c = G::random(&mut rng);
        assert!(a + b == b + a);
        assert!(a - b == -(b - a));
        assert!(c + (a + b) == a + (c + b));
        assert!((a - b) - c == (a - c) - b);

        let a = G::generator().double().double(); // 4P
        let b = G::generator().double(); // 2P
        let c = a + b;

        let mut d = G::generator();
        for _ in 0..5 {
            d += G::generator();
        }

        assert!(c == d);
        assert!(!bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        assert!(!bool::from(d.is_identity()));
        assert!(bool::from(d.is_on_curve()));
    }

    fn mixed_addition<G: CurveExt>() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let a = G::identity();
        let b = G::AffineRepr::identity();
        let c = a + b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));
        let c = a - b;
        assert!(bool::from(c.is_identity()));
        assert!(bool::from(c.is_on_curve()));

        let a = G::identity();
        let a = -a;
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(a.is_identity()));
        let a = G::AffineExt::identity();
        let a = -a;
        assert!(bool::from(a.is_on_curve()));
        assert!(bool::from(a.is_identity()));

        let a: G::AffineExt = G::random(&mut rng).into();
        assert!(a.to_curve() == a + G::AffineExt::identity());

        let a = G::random(&mut rng);
        assert!(a.double() - a == a);

        let a = G::random(&mut rng);
        let b: G::AffineExt = G::random(&mut rng).into();
        let c0 = a + b;
        let c1 = a + G::from(b);
        assert_eq!(c0, c1);
    }

    fn batch_normalize<G: CurveExt>() {
        let a = G::generator().double();
        let b = a.double();
        let c = b.double();

        for a_identity in (0..1).map(|n| n == 1) {
            for b_identity in (0..1).map(|n| n == 1) {
                for c_identity in (0..1).map(|n| n == 1) {
                    let mut v = [a, b, c];
                    if a_identity {
                        v[0] = G::identity()
                    }
                    if b_identity {
                        v[1] = G::identity()
                    }
                    if c_identity {
                        v[2] = G::identity()
                    }

                    let mut t = [
                        G::AffineExt::identity(),
                        G::AffineExt::identity(),
                        G::AffineExt::identity(),
                    ];
                    let expected = [
                        G::AffineExt::from(v[0]),
                        G::AffineExt::from(v[1]),
                        G::AffineExt::from(v[2]),
                    ];

                    G::batch_normalize(&v[..], &mut t[..]);

                    assert_eq!(&t[..], &expected[..]);
                }
            }
        }
    }

    fn multiplication<G: CurveExt>() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);

        let s1 = G::ScalarExt::random(&mut rng);
        let s2 = G::ScalarExt::random(&mut rng);

        let t0 = G::identity() * s1;
        assert!(bool::from(t0.is_identity()));

        let a = G::random(&mut rng);
        let t0 = a * G::ScalarExt::one();
        assert_eq!(a, t0);

        let t0 = a * G::ScalarExt::zero();
        assert!(bool::from(t0.is_identity()));

        let t0 = a * s1 + a * s2;

        let s3 = s1 + s2;
        let t1 = a * s3;

        assert_eq!(t0, t1);

        let mut t0 = a * s1;
        let mut t1 = a * s2;
        t0 += t1;
        let s3 = s1 + s2;
        t1 = a * s3;
        assert_eq!(t0, t1);
    }

    #[test]
    fn test_cofactor() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        let a = G2::random(&mut rng);
        assert!(bool::from(a.is_torsion_free()));
        let mut a = <G2 as group::Group>::random(&mut rng);
        assert!(!bool::from(a.is_torsion_free()));
        a = a.clear_cofactor();
        assert!(bool::from(a.is_torsion_free()));
    }

    #[test]
    fn curve_tests() {
        is_on_curve::<G1>();
        equality::<G1>();
        projective_to_affine_affine_to_projective::<G1>();
        projective_addition::<G1>();
        mixed_addition::<G1>();
        multiplication::<G1>();
        batch_normalize::<G1>();
        is_on_curve::<G2>();
        equality::<G2>();
        projective_to_affine_affine_to_projective::<G2>();
        projective_addition::<G2>();
        mixed_addition::<G2>();
        multiplication::<G2>();
        batch_normalize::<G1>();
    }

    #[test]
    fn test_lc() {
        let mut rng = XorShiftRng::from_seed([
            0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
            0xbc, 0xe5,
        ]);
        use crate::arithmetic::LinearCombinationEngine;
        let base = Fr::from_u64(10);
        let mut lc = PointCombination::new(base);
        let u0 = G1Affine::random(&mut rng).into();
        let u1 = G1Affine::random(&mut rng).into();
        let u2 = G1Affine::random(&mut rng).into();

        lc.add(u0);
        lc.add(u1);
        lc.add(u2);

        let res = u0 + u1 * base + u2 * (base * base);

        assert!(res == lc.result());
        assert!(res == PointCombination::combine(base, vec![u0, u1, u2]));

        let z0 = Fr::random(&mut rng);
        let z1 = Fr::random(&mut rng);
        let z2 = Fr::random(&mut rng);

        lc.add_with_aux(u0, z0);
        lc.add_with_aux(u1, z1);
        lc.add_with_aux(u2, z2);

        let res = u0 * z0 + u1 * (base * z1) + u2 * (base * base * z2);
        assert!(res == lc.result());
    }
}

impl group::UncompressedEncoding for G1Affine {
    type Uncompressed = G1Compressed;

    fn from_uncompressed(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn from_uncompressed_unchecked(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        unimplemented!();
    }
}

impl group::UncompressedEncoding for G2Affine {
    type Uncompressed = G2Compressed;

    fn from_uncompressed(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn from_uncompressed_unchecked(_: &Self::Uncompressed) -> CtOption<Self> {
        unimplemented!();
    }

    fn to_uncompressed(&self) -> Self::Uncompressed {
        unimplemented!();
    }
}
