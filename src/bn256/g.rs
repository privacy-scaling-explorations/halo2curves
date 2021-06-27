use super::fq::Fq;
use super::fq2::Fq2;
use super::fr::Fr;
use crate::arithmetic::{BaseExt, Coordinates, CurveAffine, CurveExt, Group};
// use crate::G1Prepared;
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use group::{prime::PrimeCurveAffine, Curve as _, Group as _, GroupEncoding};
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

// var b2 = &fe2{
// 	fe{0x3bf938e377b802a8, 0x020b1b273633535d, 0x26b7edf049755260, 0x2514c6324384a86d},
// 	fe{0x38e7ecccd1dcff67, 0x65f0b37d93ce0d3e, 0xd749d0dd22ac00aa, 0x0141b9ce4a688d4d},
// }

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

// 0x1800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed
// 0x198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c2
// 0x12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa
// 0x090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b
// 0x2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5
// 0x009713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2

#[cfg(test)]
mod tests {

    use super::{G1, G2};
    use ff::Field;

    use crate::arithmetic::{BaseExt, Coordinates, CurveAffine, CurveExt};
    use group::{prime::PrimeCurveAffine, Curve, Group, GroupEncoding};
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

        for i in 0..100 {
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

        let mut t0 = a * s1 + a * s2;

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
    fn curve_tests() {
        is_on_curve::<G1>();
        equality::<G1>();
        projective_to_affine_affine_to_projective::<G1>();
        projective_addition::<G1>();
        mixed_addition::<G1>();
        multiplication::<G1>();
        is_on_curve::<G2>();
        equality::<G2>();
        projective_to_affine_affine_to_projective::<G2>();
        projective_addition::<G2>();
        mixed_addition::<G2>();
        multiplication::<G2>();
    }
}

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
