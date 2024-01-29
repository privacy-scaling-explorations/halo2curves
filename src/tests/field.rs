#[macro_export]
macro_rules! field_testing_suite {
    ($field: ident, "field_arithmetic") => {
        macro_rules! random_multiplication_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    let a = $f::random(&mut $rng);
                    let b = $f::random(&mut $rng);
                    let c = $f::random(&mut $rng);

                    let mut t0 = a; // (a * b) * c
                    t0.mul_assign(&b);
                    t0.mul_assign(&c);

                    let mut t1 = a; // (a * c) * b
                    t1.mul_assign(&c);
                    t1.mul_assign(&b);

                    let mut t2 = b; // (b * c) * a
                    t2.mul_assign(&c);
                    t2.mul_assign(&a);

                    assert_eq!(t0, t1);
                    assert_eq!(t1, t2);
                }
            };
        }

        macro_rules! random_addition_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    let a = $f::random(&mut $rng);
                    let b = $f::random(&mut $rng);
                    let c = $f::random(&mut $rng);

                    let mut t0 = a; // (a + b) + c
                    t0.add_assign(&b);
                    t0.add_assign(&c);

                    let mut t1 = a; // (a + c) + b
                    t1.add_assign(&c);
                    t1.add_assign(&b);

                    let mut t2 = b; // (b + c) + a
                    t2.add_assign(&c);
                    t2.add_assign(&a);

                    assert_eq!(t0, t1);
                    assert_eq!(t1, t2);
                }
            };
        }

        macro_rules! random_subtraction_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    let a = $f::random(&mut $rng);
                    let b = $f::random(&mut $rng);

                    let mut t0 = a; // (a - b)
                    t0.sub_assign(&b);

                    let mut t1 = b; // (b - a)
                    t1.sub_assign(&a);

                    let mut t2 = t0; // (a - b) + (b - a) = 0
                    t2.add_assign(&t1);

                    assert_eq!(t2.is_zero().unwrap_u8(), 1);
                }
            };
        }

        macro_rules! random_negation_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    let a = $f::random(&mut $rng);
                    let mut b = a;
                    b = b.neg();
                    b.add_assign(&a);

                    assert_eq!(b.is_zero().unwrap_u8(), 1);
                }
            };
        }

        macro_rules! random_doubling_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    let mut a = $f::random(&mut $rng);
                    let mut b = a;
                    a.add_assign(&b);
                    b = b.double();

                    assert_eq!(a, b);
                }
            };
        }

        macro_rules! random_squaring_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    let mut a = $f::random(&mut $rng);
                    let mut b = a;
                    a.mul_assign(&b);
                    b = b.square();

                    assert_eq!(a, b);
                }
            };
        }

        macro_rules! random_inversion_tests {
            ($f: ident, $rng: expr) => {
                assert!(bool::from($f::ZERO.invert().is_none()));

                for _ in 0..1000000 {
                    let mut a = $f::random(&mut $rng);
                    let b = a.invert().unwrap(); // probablistically nonzero
                    a.mul_assign(&b);

                    assert_eq!(a, $f::ONE);
                }
            };
        }

        macro_rules! random_expansion_tests {
            ($f: ident, $rng: expr) => {
                for _ in 0..1000000 {
                    // Compare (a + b)(c + d) and (a*c + b*c + a*d + b*d)

                    let a = $f::random(&mut $rng);
                    let b = $f::random(&mut $rng);
                    let c = $f::random(&mut $rng);
                    let d = $f::random(&mut $rng);

                    let mut t0 = a;
                    t0.add_assign(&b);
                    let mut t1 = c;
                    t1.add_assign(&d);
                    t0.mul_assign(&t1);

                    let mut t2 = a;
                    t2.mul_assign(&c);
                    let mut t3 = b;
                    t3.mul_assign(&c);
                    let mut t4 = a;
                    t4.mul_assign(&d);
                    let mut t5 = b;
                    t5.mul_assign(&d);

                    t2.add_assign(&t3);
                    t2.add_assign(&t4);
                    t2.add_assign(&t5);

                    assert_eq!(t0, t2);
                }
            };
        }

        macro_rules! zero_tests {
            ($f: ident, $rng: expr) => {
                assert_eq!($f::ZERO.is_zero().unwrap_u8(), 1);
                {
                    let mut z = $f::ZERO;
                    z = z.neg();
                    assert_eq!(z.is_zero().unwrap_u8(), 1);
                }

                assert!(bool::from($f::ZERO.invert().is_none()));

                // Multiplication by zero
                {
                    let mut a = $f::random(&mut $rng);
                    a.mul_assign(&$f::ZERO);
                    assert_eq!(a.is_zero().unwrap_u8(), 1);
                }

                // Addition by zero
                {
                    let mut a = $f::random(&mut $rng);
                    let copy = a;
                    a.add_assign(&$f::ZERO);
                    assert_eq!(a, copy);
                }
            };
        }

        macro_rules! one_tests {
            ($f: ident, $rng: expr) => {
                assert!(bool::from($f::ONE.invert().is_some()));

                // Multiplication by one
                {
                    let mut a = $f::random(&mut $rng);
                    let copy = a;
                    a.mul_assign(&$f::ONE);
                    assert_eq!(a, copy);
                }

                // Addition by one
                {
                    let mut a = $f::random(&mut $rng);
                    let copy = a;
                    a.add_assign(&$f::ONE);
                    assert_eq!(a, copy + $f::ONE);
                }
            };
        }

        use ff::Field;
        use rand::SeedableRng;
        use rand_xorshift::XorShiftRng;
        use std::ops::*;

        #[test]
        fn test_field() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
                0xbc, 0xe5,
            ]);

            // normal cases
            random_multiplication_tests!($field, rng);
            random_addition_tests!($field, rng);
            random_subtraction_tests!($field, rng);
            random_negation_tests!($field, rng);
            random_doubling_tests!($field, rng);
            random_squaring_tests!($field, rng);
            random_inversion_tests!($field, rng);
            random_expansion_tests!($field, rng);

            // edge cases
            zero_tests!($field, rng);
            one_tests!($field, rng);
        }
    };

    ($field: ident, "conversion") => {
        macro_rules! random_conversion_tests {
            ($f: ident) => {
                let mut rng = XorShiftRng::from_seed([
                    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54,
                    0x06, 0xbc, 0xe5,
                ]);
                for _ in 0..1000000 {
                    let a = $f::random(&mut rng);
                    let bytes = a.to_repr();
                    let b = $f::from_repr(bytes).unwrap();
                    assert_eq!(a, b);
                }
            };
        }

        #[test]
        fn test_conversion() {
            random_conversion_tests!($field);
        }
    };

    ($field: ident, "serialization") => {
        macro_rules! random_serialization_test {
            ($f: ident) => {
                let mut rng = XorShiftRng::from_seed([
                    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54,
                    0x06, 0xbc, 0xe5,
                ]);
                for _ in 0..1000000 {
                    let a = $f::random(&mut rng);
                    let bytes = a.to_raw_bytes();
                    let b = $f::from_raw_bytes(&bytes).unwrap();
                    assert_eq!(a, b);
                    let mut buf = Vec::new();
                    a.write_raw(&mut buf).unwrap();
                    let b = $f::read_raw(&mut &buf[..]).unwrap();
                    assert_eq!(a, b);
                }
            };
        }

        #[cfg(feature = "derive_serde")]
        macro_rules! random_serde_test {
            ($f: ident) => {
                let mut rng = XorShiftRng::from_seed([
                    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54,
                    0x06, 0xbc, 0xe5,
                ]);
                for _ in 0..1000000 {
                    // byte serialization
                    let a = $f::random(&mut rng);
                    let bytes = bincode::serialize(&a).unwrap();
                    let reader = std::io::Cursor::new(bytes);
                    let b: $f = bincode::deserialize_from(reader).unwrap();
                    assert_eq!(a, b);

                    // json serialization
                    let json = serde_json::to_string(&a).unwrap();
                    let reader = std::io::Cursor::new(json);
                    let b: $f = serde_json::from_reader(reader).unwrap();
                    assert_eq!(a, b);
                }
            };
        }

        #[test]
        fn test_serialization() {
            use crate::serde::SerdeObject;
            random_serialization_test!($field);
            #[cfg(feature = "derive_serde")]
            random_serde_test!($field);
        }
    };

    ($field: ident, "quadratic_residue") => {
        macro_rules! random_quadratic_residue_test {
            ($f: ident) => {
                let mut rng = XorShiftRng::from_seed([
                    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54,
                    0x06, 0xbc, 0xe5,
                ]);
                for _ in 0..100000 {
                    let elem = $f::random(&mut rng);
                    let is_quad_res_or_zero: bool = elem.sqrt().is_some().into();
                    let is_quad_non_res: bool = elem.ct_quadratic_non_residue().into();
                    assert_eq!(!is_quad_non_res, is_quad_res_or_zero)
                }
            };
        }

        #[test]
        fn test_quadratic_residue() {
            use crate::ff_ext::Legendre;
            use ff::Field;
            use rand_core::SeedableRng;
            use rand_xorshift::XorShiftRng;
            random_quadratic_residue_test!($field);
        }
    };

    ($field: ident, "bits") => {
        #[cfg(feature = "bits")]
        macro_rules! random_bits_tests {
            ($f: ident) => {
                let mut rng = XorShiftRng::from_seed([
                    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54,
                    0x06, 0xbc, 0xe5,
                ]);
                for _ in 0..1000000 {
                    let a = $f::random(&mut rng);
                    let bytes = a.to_repr();
                    let bits = a.to_le_bits();
                    for idx in 0..bits.len() {
                        assert_eq!(bits[idx], ((bytes.as_ref()[idx / 8] >> (idx % 8)) & 1) == 1);
                    }
                }
            };
        }

        #[test]
        #[cfg(feature = "bits")]
        fn test_bits() {
            use ff::PrimeFieldBits;
            random_bits_tests!($field);
        }
    };

    ($field: ident, "serialization_check") => {
        fn is_less_than<const N: usize>(x: &[u64; N], y: &[u64; N]) -> bool {
            for i in (1..N).rev() {
                match x[i].cmp(&y[i]) {
                    core::cmp::Ordering::Less => return true,
                    core::cmp::Ordering::Greater => return false,
                    _ => {}
                }
            }
            x[0].lt(&y[0])
        }

        #[test]
        fn test_serialization_check() {
            use crate::serde::SerdeObject;
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06,
                0xbc, 0xe5,
            ]);
            const LIMBS: usize = $field::size() / 8;
            // failure check
            for _ in 0..1000000 {
                let rand_word = [(); LIMBS].map(|_| rng.next_u64());
                let a = $field(rand_word);
                let rand_bytes = a.to_raw_bytes();
                match is_less_than::<LIMBS>(&rand_word, &MODULUS.0) {
                    false => {
                        assert!($field::from_raw_bytes(&rand_bytes).is_none());
                    }
                    _ => {
                        assert_eq!($field::from_raw_bytes(&rand_bytes), Some(a));
                    }
                }
            }
        }
    };

    ($field: ident, "constants", $modulus_str: expr) => {
        #[test]
        fn test_primefield_constants() {
            assert_eq!($field::MODULUS, $modulus_str);
            assert_eq!(
                $field::ROOT_OF_UNITY_INV,
                $field::ROOT_OF_UNITY.invert().unwrap()
            );
            assert_eq!($field::from(2) * $field::TWO_INV, $field::ONE);
            if $field::S != 0 {
                assert_eq!(
                    $field::ROOT_OF_UNITY.pow_vartime([1 << $field::S]),
                    $field::one()
                );
                assert_eq!(
                    $field::DELTA,
                    $field::MULTIPLICATIVE_GENERATOR.pow([1u64 << $field::S])
                );
            }
        }
    };

    ($field: ident, "sqrt") => {
        #[test]
        fn test_sqrt() {
            use crate::ff_ext::Legendre;
            use rand_core::OsRng;

            let v = ($field::TWO_INV).square().sqrt().unwrap();
            assert!(v == $field::TWO_INV || (-v) == $field::TWO_INV);

            for _ in 0..10000 {
                let a = $field::random(OsRng);
                if a.legendre() == -1 {
                    assert!(bool::from(a.sqrt().is_none()));
                }
            }

            for _ in 0..10000 {
                let a = $field::random(OsRng);
                let mut b = a;
                b = b.square();
                assert_eq!(b.legendre(), 1);

                let b = b.sqrt().unwrap();
                let mut negb = b;
                negb = negb.neg();

                assert!(a == b || a == negb);
            }

            let mut c = $field::one();
            for _ in 0..10000 {
                let mut b = c;
                b = b.square();
                assert_eq!(b.legendre(), 1);

                b = b.sqrt().unwrap();

                if b != c {
                    b = b.neg();
                }

                assert_eq!(b, c);

                c += &$field::one();
            }
        }
    };

    ($field: ident, "zeta" $(, $base_field: ident)*) => {
        #[test]
        fn test_zeta() {
            assert_eq!($field::ZETA * $field::ZETA * $field::ZETA, $field::ONE);
            assert_ne!($field::ZETA * $field::ZETA, $field::ONE);
            $(
                let zeta = $field::new($base_field::ZETA.square(), $base_field::zero());
                assert_eq!(zeta, $field::ZETA);
            )*
        }
    };

    ($field: ident, "from_uniform_bytes", $test_vectors: expr) => {
        #[test]
        fn test_from_uniform_bytes() {
            const N_VECS: usize = 10;
            assert!($test_vectors.len() == N_VECS);

            let mut seeded_rng = XorShiftRng::seed_from_u64(0u64);
            let uniform_bytes = std::iter::from_fn(|| {
                let mut bytes = [0u8; 64];
                seeded_rng.fill_bytes(&mut bytes);
                Some(bytes)
            })
            .take(N_VECS)
            .collect::<Vec<_>>();

            for i in 0..N_VECS {
                let q = $field::from_uniform_bytes(&uniform_bytes[i]);
                assert_eq!($test_vectors[i], q);
            }
        }
    };

    ($ext_field: ident, "f2_tests", $base_field: ident) => {
        #[test]
        fn test_ser() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            let a0 = $ext_field::random(&mut rng);
            let a_bytes = a0.to_bytes();
            let a1 = $ext_field::from_bytes(&a_bytes).unwrap();
            assert_eq!(a0, a1);
        }

        #[test]
        fn test_f2_ordering() {
            let mut a = $ext_field {
                c0: $base_field::zero(),
                c1: $base_field::zero(),
            };

            let mut b = a;

            assert!(a.cmp(&b) == Ordering::Equal);
            b.c0 += &$base_field::one();
            assert!(a.cmp(&b) == Ordering::Less);
            a.c0 += &$base_field::one();
            assert!(a.cmp(&b) == Ordering::Equal);
            b.c1 += &$base_field::one();
            assert!(a.cmp(&b) == Ordering::Less);
            a.c0 += &$base_field::one();
            assert!(a.cmp(&b) == Ordering::Less);
            a.c1 += &$base_field::one();
            assert!(a.cmp(&b) == Ordering::Greater);
            b.c0 += &$base_field::one();
            assert!(a.cmp(&b) == Ordering::Equal);
        }

        #[test]
        fn test_f2_basics() {
            assert_eq!(
                $ext_field {
                    c0: $base_field::zero(),
                    c1: $base_field::zero(),
                },
                $ext_field::ZERO
            );
            assert_eq!(
                $ext_field {
                    c0: $base_field::one(),
                    c1: $base_field::zero(),
                },
                $ext_field::ONE
            );
            assert_eq!($ext_field::ZERO.is_zero().unwrap_u8(), 1);
            assert_eq!($ext_field::ONE.is_zero().unwrap_u8(), 0);
            assert_eq!(
                $ext_field {
                    c0: $base_field::zero(),
                    c1: $base_field::one(),
                }
                .is_zero()
                .unwrap_u8(),
                0
            );
        }
    };

    ($ext_field: ident, "f6_tests", $base_field: ident) => {
        #[test]
        fn test_f6_mul_nonresidue() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            let nqr = $ext_field {
                c0: $base_field::zero(),
                c1: $base_field::one(),
                c2: $base_field::zero(),
            };

            for _ in 0..1000 {
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;
                a.mul_by_nonresidue();
                b.mul_assign(&nqr);

                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_f6_mul_by_1() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..1000 {
                let c1 = $base_field::random(&mut rng);
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;

                a.mul_by_1(&c1);
                b.mul_assign(&$ext_field {
                    c0: $base_field::zero(),
                    c1,
                    c2: $base_field::zero(),
                });

                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_f6_mul_by_01() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..1000 {
                let c0 = $base_field::random(&mut rng);
                let c1 = $base_field::random(&mut rng);
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;

                a.mul_by_01(&c0, &c1);
                b.mul_assign(&$ext_field {
                    c0,
                    c1,
                    c2: $base_field::zero(),
                });

                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_squaring() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..1000 {
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;
                b.mul_assign(&a);
                a.square_assign();
                assert_eq!(a, b);
            }
        }
    };

    ($ext_field: ident, "f12_tests", $base_field_1: ident, $base_field_2: ident) => {
        #[test]
        fn test_f12_mul_by_014() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..1000 {
                let c0 = $base_field_2::random(&mut rng);
                let c1 = $base_field_2::random(&mut rng);
                let c5 = $base_field_2::random(&mut rng);
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;

                a.mul_by_014(&c0, &c1, &c5);
                b.mul_assign(&$ext_field {
                    c0: $base_field_1 {
                        c0,
                        c1,
                        c2: $base_field_2::zero(),
                    },
                    c1: $base_field_1 {
                        c0: $base_field_2::zero(),
                        c1: c5,
                        c2: $base_field_2::zero(),
                    },
                });

                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_f12_mul_by_034() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..1000 {
                let c0 = $base_field_2::random(&mut rng);
                let c3 = $base_field_2::random(&mut rng);
                let c4 = $base_field_2::random(&mut rng);
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;

                a.mul_by_034(&c0, &c3, &c4);
                b.mul_assign(&$ext_field {
                    c0: $base_field_1 {
                        c0,
                        c1: $base_field_2::zero(),
                        c2: $base_field_2::zero(),
                    },
                    c1: $base_field_1 {
                        c0: c3,
                        c1: c4,
                        c2: $base_field_2::zero(),
                    },
                });

                assert_eq!(a, b);
            }
        }

        #[test]
        fn test_squaring() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..1000 {
                let mut a = $ext_field::random(&mut rng);
                let mut b = a;
                b.mul_assign(&a);
                a.square_assign();
                assert_eq!(a, b);
            }
        }
    };

    ($ext_field: ident, "frobenius", $frobenius_param: expr) => {
        #[test]
        fn test_frobenius() {
            let mut rng = XorShiftRng::from_seed([
                0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
                0xe5,
            ]);

            for _ in 0..50 {
                for i in 0..8 {
                    let mut a = $ext_field::random(&mut rng);
                    let mut b = a;

                    for _ in 0..i {
                        a = a.pow($frobenius_param);
                    }
                    b.frobenius_map(i);

                    assert_eq!(a, b);
                }
            }
        }
    };
}
