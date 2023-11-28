use ark_std::{end_timer, start_timer};
use ff::{Field, PrimeField};
use group::{Curve, Group};
use pasta_curves::arithmetic::CurveAffine;
use rand_core::OsRng;

use crate::{
    bn256::{Fr, G1Affine, G1},
    multicore,
};

pub fn multiexp_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let c = if bases.len() < 4 {
        1
    } else if bases.len() < 32 {
        3
    } else {
        (f64::from(bases.len() as u32)).ln().ceil() as usize
    };

    fn get_at<F: PrimeField>(segment: usize, c: usize, bytes: &F::Repr) -> usize {
        let skip_bits = segment * c;
        let skip_bytes = skip_bits / 8;

        if skip_bytes >= 32 {
            return 0;
        }

        let mut v = [0; 8];
        for (v, o) in v.iter_mut().zip(bytes.as_ref()[skip_bytes..].iter()) {
            *v = *o;
        }

        let mut tmp = u64::from_le_bytes(v);
        tmp >>= skip_bits - (skip_bytes * 8);
        tmp %= 1 << c;

        tmp as usize
    }

    let segments = (256 / c) + 1;

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            *acc = acc.double();
        }

        #[derive(Clone, Copy)]
        enum Bucket<C: CurveAffine> {
            None,
            Affine(C),
            Projective(C::Curve),
        }

        impl<C: CurveAffine> Bucket<C> {
            fn add_assign(&mut self, other: &C) {
                *self = match *self {
                    Bucket::None => Bucket::Affine(*other),
                    Bucket::Affine(a) => Bucket::Projective(a + *other),
                    Bucket::Projective(mut a) => {
                        a += *other;
                        Bucket::Projective(a)
                    }
                }
            }

            fn add(self, mut other: C::Curve) -> C::Curve {
                match self {
                    Bucket::None => other,
                    Bucket::Affine(a) => {
                        other += a;
                        other
                    }
                    Bucket::Projective(a) => other + a,
                }
            }
        }

        let mut buckets: Vec<Bucket<C>> = vec![Bucket::None; (1 << c) - 1];

        for (coeff, base) in coeffs.iter().zip(bases.iter()) {
            let coeff = get_at::<C::Scalar>(current_segment, c, coeff);
            if coeff != 0 {
                buckets[coeff - 1].add_assign(base);
            }
        }

        // Summation by parts
        // e.g. 3a + 2b + 1c = a +
        //                    (a) + b +
        //                    ((a) + b) + c
        let mut running_sum = C::Curve::identity();
        for exp in buckets.into_iter().rev() {
            running_sum = exp.add(running_sum);
            *acc += &running_sum;
        }
    }
}

pub fn multiexp_serial_2<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let c = if bases.len() < 4 {
        1
    } else if bases.len() < 32 {
        3
    } else {
        (f64::from(bases.len() as u32)).ln().ceil() as usize
    };

    let segments = (256 / c) + 1;

    for current_segment in (0..segments).rev() {
        for _ in 0..c {
            *acc = acc.double();
        }

        #[derive(Clone, Copy)]
        enum Bucket<C: CurveAffine> {
            None,
            Affine(C),
            Projective(C::Curve),
        }

        impl<C: CurveAffine> Bucket<C> {
            fn add_assign(&mut self, other: &C) {
                *self = match *self {
                    Bucket::None => Bucket::Affine(*other),
                    Bucket::Affine(a) => Bucket::Projective(a + *other),
                    Bucket::Projective(mut a) => {
                        a += *other;
                        Bucket::Projective(a)
                    }
                }
            }

            fn add(self, mut other: C::Curve) -> C::Curve {
                match self {
                    Bucket::None => other,
                    Bucket::Affine(a) => {
                        other += a;
                        other
                    }
                    Bucket::Projective(a) => other + a,
                }
            }
        }

        let mut buckets: Vec<Bucket<C>> = vec![Bucket::None; 1 << (c - 1)];

        for (coeff, base) in coeffs.iter().zip(bases.iter()) {
            let coeff = get_booth_index(current_segment, c, coeff.as_ref());
            if coeff.is_positive() {
                buckets[coeff as usize - 1].add_assign(base);
            }
            if coeff.is_negative() {
                let coeff = coeff.abs();
                buckets[coeff as usize - 1].add_assign(&base.neg());
            }
        }

        // Summation by parts
        // e.g. 3a + 2b + 1c = a +
        //                    (a) + b +
        //                    ((a) + b) + c
        let mut running_sum = C::Curve::identity();
        for exp in buckets.into_iter().rev() {
            running_sum = exp.add(running_sum);
            *acc += &running_sum;
        }
    }
}

/// Performs a small multi-exponentiation operation.
/// Uses the double-and-add algorithm with doublings shared across points.
pub fn small_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();
    let mut acc = C::Curve::identity();

    // for byte idx
    for byte_idx in (0..32).rev() {
        // for bit idx
        for bit_idx in (0..8).rev() {
            acc = acc.double();
            // for each coeff
            for coeff_idx in 0..coeffs.len() {
                let byte = coeffs[coeff_idx].as_ref()[byte_idx];
                if ((byte >> bit_idx) & 1) != 0 {
                    acc += bases[coeff_idx];
                }
            }
        }
    }

    acc
}

/// Performs a multi-exponentiation operation.
///
/// This function will panic if coeffs and bases have a different length.
///
/// This will use multithreading if beneficial.
pub fn best_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    assert_eq!(coeffs.len(), bases.len());

    let num_threads = multicore::current_num_threads();
    if coeffs.len() > num_threads {
        let chunk = coeffs.len() / num_threads;
        let num_chunks = coeffs.chunks(chunk).len();
        let mut results = vec![C::Curve::identity(); num_chunks];
        multicore::scope(|scope| {
            let chunk = coeffs.len() / num_threads;

            for ((coeffs, bases), acc) in coeffs
                .chunks(chunk)
                .zip(bases.chunks(chunk))
                .zip(results.iter_mut())
            {
                scope.spawn(move |_| {
                    multiexp_serial(coeffs, bases, acc);
                });
            }
        });
        results.iter().fold(C::Curve::identity(), |a, b| a + b)
    } else {
        let mut acc = C::Curve::identity();
        multiexp_serial(coeffs, bases, &mut acc);
        acc
    }
}

pub fn best_multiexp_2<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
    assert_eq!(coeffs.len(), bases.len());

    let num_threads = multicore::current_num_threads();
    if coeffs.len() > num_threads {
        let chunk = coeffs.len() / num_threads;
        let num_chunks = coeffs.chunks(chunk).len();
        let mut results = vec![C::Curve::identity(); num_chunks];
        multicore::scope(|scope| {
            let chunk = coeffs.len() / num_threads;

            for ((coeffs, bases), acc) in coeffs
                .chunks(chunk)
                .zip(bases.chunks(chunk))
                .zip(results.iter_mut())
            {
                scope.spawn(move |_| {
                    multiexp_serial_2(coeffs, bases, acc);
                });
            }
        });
        results.iter().fold(C::Curve::identity(), |a, b| a + b)
    } else {
        let mut acc = C::Curve::identity();
        multiexp_serial_2(coeffs, bases, &mut acc);
        acc
    }
}

fn div_ceil(a: u32, b: u32) -> u32 {
    a.checked_sub(1).map_or(0, |a| a / b + 1)
}

pub(crate) fn get_booth_index(segment: usize, window: usize, el: &[u8]) -> i32 {
    let (skip_bits, pad) = match (segment * window).checked_sub(1) {
        Some(skip_bits) => (skip_bits, false),
        None => (0, true),
    };

    let skip_bytes = skip_bits / 8;
    if skip_bytes >= 32 {
        return 0;
    }

    let mut v = [0; 4];
    for (v, o) in v.iter_mut().zip(el.iter().skip(skip_bytes)) {
        *v = *o;
    }
    let mut tmp = u32::from_le_bytes(v);
    if pad {
        tmp <<= 1; // pad left with one 0
    }
    tmp >>= skip_bits - (skip_bytes * 8);
    tmp %= 1 << (window + 1);
    // let bits = format!("T {:0>width$b}", tmp, width = window + 1);
    // println!("{}", bits);
    // tmp

    let sign = tmp & (1 << window) == 0;

    let mask = (1 << window) - 1;

    if sign {
        let idx = div_ceil(tmp, 2u32);
        idx as i32
    } else {
        let idx = !div_ceil(tmp, 2u32).saturating_sub(1) & mask;
        -(idx as i32)
    }
}

#[test]
fn get_bucket_index() {
    let window = 5;

    fn mul(a: Fr, b: Fr, window: usize) -> Fr {
        let u = b.to_repr();
        let n = div_ceil(Fr::NUM_BITS, window as u32) + 1;
        // let n = 10;

        let mut acc = Fr::ZERO;
        for i in (0..n).rev() {
            let idx = get_booth_index(i as usize, window, u.as_ref());

            let tmp = a * Fr::from(idx.abs() as u64);
            // println!("{:?}", idx);
            if idx.is_negative() {
                acc -= tmp;
            }
            if idx.is_positive() {
                acc += tmp;
            }
            if i != 0 {
                for _ in 0..window {
                    acc = acc.double();
                }
            }
        }

        acc
    }

    for _ in 0..10000 {
        let a = Fr::random(OsRng);
        let b = Fr::random(OsRng);
        let c0 = mul(a, b, window);
        let c1 = a * b;
        assert_eq!(c0, c1);
    }
}

#[test]
fn test_msm_with_booth() {
    // let n = 100;

    // let points = (0..n).map(|_| G1Affine::random(OsRng)).collect::<Vec<_>>();
    // let scalars = (0..n).map(|_| Fr::random(OsRng)).collect::<Vec<_>>();

    // let mut e0 = G1::identity();
    // multiexp_serial(&scalars[..], &points[..], &mut e0);

    // let mut e1 = G1::identity();
    // multiexp_serial_2(&scalars[..], &points[..], &mut e1);
    // assert_eq!(e0, e1);

    let n = 1 << 21;

    let points = (0..n).map(|_| G1Affine::random(OsRng)).collect::<Vec<_>>();
    let scalars = (0..n).map(|_| Fr::random(OsRng)).collect::<Vec<_>>();

    let t0 = start_timer!(|| "zcash");
    let e0 = best_multiexp(&scalars[..], &points[..]);
    end_timer!(t0);

    let t1 = start_timer!(|| "booth");
    let e1 = best_multiexp_2(&scalars[..], &points[..]);
    end_timer!(t1);
    assert_eq!(e0, e1);

    let t1 = start_timer!(|| "booth");
    let _e1 = best_multiexp_2(&scalars[..], &points[..]);
    end_timer!(t1);

    let t0 = start_timer!(|| "zcash");
    let _e0 = best_multiexp(&scalars[..], &points[..]);
    end_timer!(t0);
}
