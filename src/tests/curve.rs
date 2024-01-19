#![allow(clippy::eq_op)]

use crate::ff_ext::Legendre;
use crate::{hash_to_curve, CurveAffine, CurveExt};
use ff::PrimeField;
use num_bigint::BigUint;
use num_traits::Num;
use std::borrow::Cow;

fn fe_from_str<F: PrimeField>(string: impl AsRef<str>) -> F {
    let string = string.as_ref();
    let oct = if let Some(hex) = string.strip_prefix("0x") {
        Cow::Owned(BigUint::from_str_radix(hex, 16).unwrap().to_string())
    } else {
        Cow::Borrowed(string)
    };
    F::from_str_vartime(&oct).unwrap()
}

pub fn svdw_map_to_curve_test<G: CurveExt>(
    z: G::Base,
    precomputed_constants: [&'static str; 4],
    test_vector: impl IntoIterator<Item = (&'static str, (&'static str, &'static str))>,
) where
    <G as CurveExt>::Base: Legendre,
{
    let [c1, c2, c3, c4] = hash_to_curve::svdw_precomputed_constants::<G>(z);
    assert_eq!([c1, c2, c3, c4], precomputed_constants.map(fe_from_str));
    for (u, (x, y)) in test_vector.into_iter() {
        let u = fe_from_str(u);
        let expected = G::AffineExt::from_xy(fe_from_str(x), fe_from_str(y)).unwrap();
        let output = hash_to_curve::svdw_map_to_curve::<G>(u, c1, c2, c3, c4, z).to_affine();
        assert_eq!(output, expected);
    }
}

#[macro_export]
macro_rules! curve_testing_suite {
    ($($curve: ident),*) => {

        macro_rules! is_on_curve {
            ($c: ident) => {
                assert!(bool::from($c::identity().is_on_curve()));
                assert!(bool::from($c::generator().is_on_curve()));
                assert!(bool::from($c::identity().is_on_curve()));
                assert!(bool::from($c::generator().is_on_curve()));

                for _ in 0..100 {
                    let point = $c::random(OsRng);
                    assert!(bool::from(point.is_on_curve()));
                    let affine_point: <$c as CurveExt>::AffineExt = point.into();
                    assert!(bool::from(affine_point.is_on_curve()));
                }
            }
        }

        macro_rules! equality {
            ($c: ident) => {
                let a = $c::generator();
                let b = $c::identity();

                assert!(a == a);
                assert!(b == b);
                assert!(a != b);
                assert!(b != a);

                for _ in 0..100 {
                    let a = $c::random(OsRng);
                    let b = $c::random(OsRng);

                    assert!(a == a);
                    assert!(b == b);
                    assert!(a != b);
                    assert!(b != a);

                    let a: <$c as CurveExt>::AffineExt = a.into();
                    let b: <$c as CurveExt>::AffineExt = b.into();

                    assert!(a == a);
                    assert!(b == b);
                    assert!(a != b);
                    assert!(b != a);
                }
            }
        }

        macro_rules! projective_to_affine_affine_to_projective {
            ($c: ident) => {
                let a = $c::generator();
                let b = $c::identity();

                assert!(bool::from(<$c as CurveExt>::AffineExt::from(a).is_on_curve()));
                assert!(!bool::from(<$c as CurveExt>::AffineExt::from(a).is_identity()));
                assert!(bool::from(<$c as CurveExt>::AffineExt::from(b).is_on_curve()));
                assert!(bool::from(<$c as CurveExt>::AffineExt::from(b).is_identity()));

                let a = <$c as CurveExt>::AffineExt::generator();
                let b = <$c as CurveExt>::AffineExt::identity();

                assert!(bool::from($c::from(a).is_on_curve()));
                assert!(!bool::from($c::from(a).is_identity()));
                assert!(bool::from($c::from(b).is_on_curve()));
                assert!(bool::from($c::from(b).is_identity()));
            }
        }

        macro_rules! projective_addition {
            ($c: ident) => {
                let a = $c::identity();
                let b = $c::identity();
                let c = a + b;
                assert!(bool::from(c.is_identity()));
                assert!(bool::from(c.is_on_curve()));
                let c = a - b;
                assert!(bool::from(c.is_identity()));
                assert!(bool::from(c.is_on_curve()));

                let a = $c::identity();
                let a = -a;
                assert!(bool::from(a.is_on_curve()));
                assert!(bool::from(a.is_identity()));

                let a = $c::random(OsRng);
                assert!(a == a + $c::identity());
                assert!(a == $c::identity() + a);
                assert!(-a == $c::identity() - a);

                let a = $c::identity();
                let a = a.double();
                assert!(bool::from(c.is_on_curve()));
                assert!(bool::from(a.is_identity()));

                let a = $c::generator();
                let a = a.double();
                assert!(bool::from(c.is_on_curve()));
                assert_eq!(a, $c::generator() + $c::generator());

                let a = $c::random(OsRng);
                assert!(a.double() - a == a);

                let a = $c::random(OsRng);
                let b = $c::random(OsRng);
                let c = $c::random(OsRng);
                assert!(a + b == b + a);
                assert!(a - b == -(b - a));
                assert!(c + (a + b) == a + (c + b));
                assert!((a - b) - c == (a - c) - b);

                let a = $c::generator().double().double(); // 4P
                let b = $c::generator().double(); // 2P
                let c = a + b;

                let mut d = $c::generator();
                for _ in 0..5 {
                    d += $c::generator();
                }

                assert!(c == d);
                assert!(!bool::from(c.is_identity()));
                assert!(bool::from(c.is_on_curve()));
                assert!(!bool::from(d.is_identity()));
                assert!(bool::from(d.is_on_curve()));
            }
        }

        macro_rules! mixed_addition {
            ($c: ident) => {
                let a = $c::identity();
                let b = <$c as group::Curve>::AffineRepr::identity();
                let c = a + b;
                assert!(bool::from(c.is_identity()));
                assert!(bool::from(c.is_on_curve()));
                let c = a - b;
                assert!(bool::from(c.is_identity()));
                assert!(bool::from(c.is_on_curve()));

                let a = $c::identity();
                let a = -a;
                assert!(bool::from(a.is_on_curve()));
                assert!(bool::from(a.is_identity()));
                let a = <$c as CurveExt>::AffineExt::identity();
                let a = -a;
                assert!(bool::from(a.is_on_curve()));
                assert!(bool::from(a.is_identity()));

                let a: <$c as CurveExt>::AffineExt = $c::random(OsRng).into();
                assert!(a.to_curve() == a + <$c as CurveExt>::AffineExt::identity());

                let a = $c::random(OsRng);
                assert!(a.double() - a == a);

                let a = $c::random(OsRng);
                let b: <$c as CurveExt>::AffineExt = $c::random(OsRng).into();
                let c0 = a + b;
                let c1 = a + $c::from(b);
                assert_eq!(c0, c1);
            }
        }

        macro_rules! multiplication {
            ($c: ident) => {
                for _ in 1..1000 {
                    let s1 = <$c as CurveExt>::ScalarExt::random(OsRng);
                    let s2 = <$c as CurveExt>::ScalarExt::random(OsRng);

                    let t0 = $c::identity() * s1;
                    assert!(bool::from(t0.is_identity()));

                    let a = $c::random(OsRng);
                    let t0 = a * <$c as CurveExt>::ScalarExt::ONE;
                    assert_eq!(a, t0);

                    let t0 = a * <$c as CurveExt>::ScalarExt::ZERO;
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
            }
        }

        macro_rules! batch_normalize {
            ($c: ident) => {
                let a = $c::generator().double();
                let b = a.double();
                let c = b.double();

                for a_identity in (0..1).map(|n| n == 1) {
                    for b_identity in (0..1).map(|n| n == 1) {
                        for c_identity in (0..1).map(|n| n == 1) {
                            let mut v = [a, b, c];
                            if a_identity {
                                v[0] = $c::identity()
                            }
                            if b_identity {
                                v[1] = $c::identity()
                            }
                            if c_identity {
                                v[2] = $c::identity()
                            }

                            let mut t = [
                                <$c as CurveExt>::AffineExt::identity(),
                                <$c as CurveExt>::AffineExt::identity(),
                                <$c as CurveExt>::AffineExt::identity(),
                            ];
                            let expected = [
                                <$c as CurveExt>::AffineExt::from(v[0]),
                                <$c as CurveExt>::AffineExt::from(v[1]),
                                <$c as CurveExt>::AffineExt::from(v[2]),
                            ];

                            $c::batch_normalize(&v[..], &mut t[..]);

                            assert_eq!(&t[..], &expected[..]);
                        }
                    }
                }
            }
        }

        macro_rules! serdes {
            ($c: ident) => {
                assert!(bool::from(
                    $c::from_bytes(&$c::identity().to_bytes())
                        .unwrap()
                        .is_identity()
                ));
                assert!(bool::from(
                    <$c as CurveExt>::AffineExt::from_bytes(&<$c as CurveExt>::AffineExt::identity().to_bytes())
                        .unwrap()
                        .is_identity()
                ));
                for _ in 0..100 {
                    let projective_point = $c::random(OsRng);
                    let affine_point: <$c as CurveExt>::AffineExt = projective_point.into();
                    let projective_repr = projective_point.to_bytes();
                    let affine_repr = affine_point.to_bytes();

                    println!(
                        "{:?} \n{:?}",
                        projective_repr.as_ref(),
                        affine_repr.as_ref()
                    );

                    let projective_point_rec = $c::from_bytes(&projective_repr).unwrap();
                    let projective_point_rec_unchecked = $c::from_bytes(&projective_repr).unwrap();
                    let affine_point_rec = <$c as CurveExt>::AffineExt::from_bytes(&affine_repr).unwrap();
                    let affine_point_rec_unchecked = <$c as CurveExt>::AffineExt::from_bytes(&affine_repr).unwrap();

                    assert_eq!(projective_point, projective_point_rec);
                    assert_eq!(projective_point, projective_point_rec_unchecked);
                    assert_eq!(affine_point, affine_point_rec);
                    assert_eq!(affine_point, affine_point_rec_unchecked);
                }
            }
        }

        macro_rules! hash_to_curve_test {
            ($c: ident) => {
                let hasher = $c::hash_to_curve("test");
                let mut rng = OsRng;
                for _ in 0..1000 {
                    let message = iter::repeat_with(|| rng.next_u32().to_be_bytes())
                        .take(32)
                        .flatten()
                        .collect::<Vec<_>>();
                    assert!(bool::from(hasher(&message).is_on_curve()));
                }
            }
        }

        macro_rules! random_serialization_test {
            ($c: ident) => {
                for _ in 0..100 {
                    let projective_point = $c::random(OsRng);
                    let affine_point: <$c as CurveExt>::AffineExt = projective_point.into();

                    let projective_bytes = projective_point.to_raw_bytes();
                    let projective_point_rec = $c::from_raw_bytes(&projective_bytes).unwrap();
                    assert_eq!(projective_point, projective_point_rec);
                    let mut buf = Vec::new();
                    projective_point.write_raw(&mut buf).unwrap();
                    let projective_point_rec = $c::read_raw(&mut &buf[..]).unwrap();
                    assert_eq!(projective_point, projective_point_rec);

                    let affine_bytes = affine_point.to_raw_bytes();
                    let affine_point_rec = <$c as CurveExt>::AffineExt::from_raw_bytes(&affine_bytes).unwrap();
                    assert_eq!(affine_point, affine_point_rec);
                    let mut buf = Vec::new();
                    affine_point.write_raw(&mut buf).unwrap();
                    let affine_point_rec = <$c as CurveExt>::AffineExt::read_raw(&mut &buf[..]).unwrap();
                    assert_eq!(affine_point, affine_point_rec);
                }
            }
        }

        macro_rules! random_serde_test {
            ($c: ident) => {
                for _ in 0..100 {
                    let projective_point = $c::random(OsRng);
                    let affine_point: <$c as CurveExt>::AffineExt = projective_point.into();
                    {
                        let affine_bytes = bincode::serialize(&affine_point).unwrap();
                        let reader = std::io::Cursor::new(affine_bytes);
                        let affine_point_rec: <$c as CurveExt>::AffineExt = bincode::deserialize_from(reader).unwrap();
                        assert_eq!(projective_point.to_affine(), affine_point_rec);
                        assert_eq!(affine_point, affine_point_rec);
                    }
                    {
                        let affine_json = serde_json::to_string(&affine_point).unwrap();
                        let reader = std::io::Cursor::new(affine_json);
                        let affine_point_rec: <$c as CurveExt>::AffineExt = serde_json::from_reader(reader).unwrap();
                        assert_eq!(affine_point, affine_point_rec);
                    }
                    {
                        let projective_bytes = bincode::serialize(&projective_point).unwrap();
                        let reader = std::io::Cursor::new(projective_bytes);
                        let projective_point_rec: $c = bincode::deserialize_from(reader).unwrap();
                        assert_eq!(projective_point, projective_point_rec);
                    }
                    {
                        let projective_json = serde_json::to_string(&projective_point).unwrap();
                        let reader = std::io::Cursor::new(projective_json);
                        let projective_point_rec: $c = serde_json::from_reader(reader).unwrap();
                        assert_eq!(projective_point, projective_point_rec);
                    }
                }
            }
        }

        #[cfg(test)]
        mod tests {
            use super::*;

            use crate::ff::Field;
            use crate::group::prime::PrimeCurveAffine;
            use crate::{group::GroupEncoding, serde::SerdeObject};
            use crate::{CurveAffine, CurveExt};
            use ff::WithSmallOrderMulGroup;
            use rand_core::{OsRng, RngCore};
            use std::iter;

            #[test]
            fn test_curve() {
                $(
                    is_on_curve!($curve);
                    equality!($curve);
                    projective_to_affine_affine_to_projective!($curve);
                    projective_addition!($curve);
                    mixed_addition!($curve);
                    multiplication!($curve);
                    batch_normalize!($curve);
                    serdes!($curve);
                )*
            }

            #[test]
            fn test_hash_to_curve() {
                $(
                    hash_to_curve_test!($curve);
                )*
            }

            #[test]
            fn test_serialization() {
                $(
                    random_serialization_test!($curve);
                    #[cfg(feature = "derive_serde")]
                    random_serde_test!($curve);
                )*
            }

            #[test]
            fn test_endo_consistency() {
                $(
                    let g = $curve::generator();
                    assert_eq!(g * <$curve as CurveExt>::ScalarExt::ZETA, g.endo());

                    for _ in 0..100 {
                        let g = $curve::random(OsRng);
                        assert_eq!(g * <$curve as CurveExt>::ScalarExt::ZETA, g.endo());
                    }
                )*
            }
        }
    };
}
