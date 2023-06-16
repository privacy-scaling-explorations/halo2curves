use std::ops::{AddAssign, MulAssign, Neg, SubAssign};

use ark_std::test_rng;
use ff::{Field, PrimeField};
// use pasta_curves::arithmetic::FieldExt;
use rand_core::RngCore;

use crate::{bn256::Fr, maybe_u64::MaybeU64Coversion};

type MockField = Fr;

// ==========================================
// New tests added for MaybeU64 structs
// ==========================================
#[test]
fn test_conversion() {
    let mut rng = test_rng();

    for _ in 0..100 {
        let a64 = rng.next_u64();
        let a = MockField::from(a64);
        let b = a.to_full();
        let c = b.to_u64();
        assert_eq!(a, c)
    }
}

#[test]
fn multiplication() {
    let mut rng = test_rng();

    for _ in 0..100 {
        let a = MockField::random_u64(&mut rng);
        let b = MockField::random_field(&mut rng);
        let c = MockField::U64(u64::MAX);

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
}
#[test]
fn subtraction() {
    let mut rng = test_rng();

    for _ in 0..100 {
        let a = MockField::random_u64(&mut rng);
        let b = MockField::U64(u64::MAX);

        let mut t0 = a; // (a - b)
        t0.sub_assign(&b);

        let mut t1 = b; // (b - a)
        t1.sub_assign(&a);

        let mut t2 = t0; // (a - b) + (b - a) = 0
        t2.add_assign(&t1);

        assert_eq!(t2.to_u64(), MockField::from(0));
    }
}
#[test]
fn addition() {
    let mut rng = test_rng();

    for _ in 0..100 {
        let a = MockField::random_u64(&mut rng);
        let b = MockField::random_field(&mut rng);
        let c = MockField::U64(u64::MAX);

        let t1 = (a + b) + c;
        let t2 = a + (b + c);
        let t3 = (a + c) + b;
        let mut t4 = a;
        t4 += b;
        t4 += c;
        let t5 = [a, b, c].iter().sum();

        assert_eq!(t1, t2);
        assert_eq!(t1, t3);
        assert_eq!(t1, t4);
        assert_eq!(t1, t5);
    }
}

// ==========================================
// tests inherent from BN254::Fr
// ==========================================
#[test]
fn test_sqrt() {
    let mut rng = test_rng();

    let v = (MockField::TWO_INV).square().sqrt().unwrap();
    assert!(v == MockField::TWO_INV || (-v) == MockField::TWO_INV);

    for _ in 0..10000 {
        let a = MockField::random(&mut rng);
        let mut b = a;
        b = b.square();

        let b = b.sqrt().unwrap();
        let mut negb = b;
        negb = negb.neg();

        assert!(a == b || a == negb);
    }
}

#[test]
fn test_root_of_unity() {
    assert_eq!(
        MockField::ROOT_OF_UNITY.pow_vartime(&[1 << MockField::S, 0, 0, 0]),
        MockField::ONE
    );
}
#[test]
fn test_inv_root_of_unity() {
    assert_eq!(
        MockField::ROOT_OF_UNITY_INV,
        MockField::ROOT_OF_UNITY.invert().unwrap()
    );
}

#[test]
fn test_field() {
    crate::tests::field::random_field_tests::<MockField>("test field".to_string());
}

#[test]
fn test_delta() {
    assert_eq!(
        MockField::DELTA,
        MockField::MULTIPLICATIVE_GENERATOR.pow(&[1u64 << MockField::S, 0, 0, 0])
    );
}

#[test]
fn test_serialization() {
    crate::tests::field::random_serialization_test::<MockField>("fr".to_string());
}
