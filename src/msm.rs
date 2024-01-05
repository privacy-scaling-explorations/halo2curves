use std::ops::Neg;

use ff::PrimeField;
use group::Group;
use pasta_curves::arithmetic::CurveAffine;

fn get_booth_index(window_index: usize, window_size: usize, el: &[u8]) -> i32 {
    // Booth encoding:
    // * step by `window` size
    // * slice by size of `window + 1``
    // * each window overlap by 1 bit
    // * append a zero bit to the least significant end
    // Indexing rule for example window size 3 where we slice by 4 bits:
    // `[0, +1, +1, +2, +2, +3, +3, +4, -4, -3, -3 -2, -2, -1, -1, 0]``
    // So we can reduce the bucket size without preprocessing scalars
    // and remembering them as in classic signed digit encoding

    let skip_bits = (window_index * window_size).saturating_sub(1);
    let skip_bytes = skip_bits / 8;

    // fill into a u32
    let mut v: [u8; 4] = [0; 4];
    for (dst, src) in v.iter_mut().zip(el.iter().skip(skip_bytes)) {
        *dst = *src
    }
    let mut tmp = u32::from_le_bytes(v);

    // pad with one 0 if slicing the least significant window
    if window_index == 0 {
        tmp <<= 1;
    }

    // remove further bits
    tmp >>= skip_bits - (skip_bytes * 8);
    // apply the booth window
    tmp &= (1 << (window_size + 1)) - 1;

    let sign = tmp & (1 << window_size) == 0;

    // div ceil by 2
    tmp = (tmp + 1) >> 1;

    // find the booth action index
    if sign {
        tmp as i32
    } else {
        ((!(tmp - 1) & ((1 << window_size) - 1)) as i32).neg()
    }
}

pub fn multiexp_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
    let coeffs: Vec<_> = coeffs.iter().map(|a| a.to_repr()).collect();

    let c = if bases.len() < 4 {
        1
    } else if bases.len() < 32 {
        3
    } else {
        (f64::from(bases.len() as u32)).ln().ceil() as usize
    };

    let number_of_windows = C::Scalar::NUM_BITS as usize / c + 1;

    for current_window in (0..number_of_windows).rev() {
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
            let coeff = get_booth_index(current_window, c, coeff.as_ref());
            if coeff.is_positive() {
                buckets[coeff as usize - 1].add_assign(base);
            }
            if coeff.is_negative() {
                buckets[coeff.unsigned_abs() as usize - 1].add_assign(&base.neg());
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

    let num_threads = rayon::current_num_threads();
    if coeffs.len() > num_threads {
        let chunk = coeffs.len() / num_threads;
        let num_chunks = coeffs.chunks(chunk).len();
        let mut results = vec![C::Curve::identity(); num_chunks];
        rayon::scope(|scope| {
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

#[cfg(test)]
mod test {

    use std::ops::Neg;

    use crate::bn256::{Fr, G1Affine, G1};
    use ark_std::{end_timer, start_timer};
    use ff::{Field, PrimeField};
    use group::{Curve, Group};
    use pasta_curves::arithmetic::CurveAffine;
    use rand_core::OsRng;

    // keeping older implementation it here for baseline comparison, debugging & benchmarking
    fn best_multiexp<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C]) -> C::Curve {
        assert_eq!(coeffs.len(), bases.len());

        let num_threads = rayon::current_num_threads();
        if coeffs.len() > num_threads {
            let chunk = coeffs.len() / num_threads;
            let num_chunks = coeffs.chunks(chunk).len();
            let mut results = vec![C::Curve::identity(); num_chunks];
            rayon::scope(|scope| {
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

    // keeping older implementation it here for baseline comparision, debugging & benchmarking
    fn multiexp_serial<C: CurveAffine>(coeffs: &[C::Scalar], bases: &[C], acc: &mut C::Curve) {
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

    #[test]
    fn test_booth_encoding() {
        fn mul(scalar: &Fr, point: &G1Affine, window: usize) -> G1Affine {
            let u = scalar.to_repr();
            let n = Fr::NUM_BITS as usize / window + 1;

            let table = (0..=1 << (window - 1))
                .map(|i| point * Fr::from(i as u64))
                .collect::<Vec<_>>();

            let mut acc = G1::identity();
            for i in (0..n).rev() {
                for _ in 0..window {
                    acc = acc.double();
                }

                let idx = super::get_booth_index(i, window, u.as_ref());

                if idx.is_negative() {
                    acc += table[idx.unsigned_abs() as usize].neg();
                }
                if idx.is_positive() {
                    acc += table[idx.unsigned_abs() as usize];
                }
            }

            acc.to_affine()
        }

        let (scalars, points): (Vec<_>, Vec<_>) = (0..10)
            .map(|_| {
                let scalar = Fr::random(OsRng);
                let point = G1Affine::random(OsRng);
                (scalar, point)
            })
            .unzip();

        for window in 1..10 {
            for (scalar, point) in scalars.iter().zip(points.iter()) {
                let c0 = mul(scalar, point, window);
                let c1 = point * scalar;
                assert_eq!(c0, c1.to_affine());
            }
        }
    }

    fn run_msm_cross<C: CurveAffine>(min_k: usize, max_k: usize) {
        let points = (0..1 << max_k)
            .map(|_| C::Curve::random(OsRng))
            .collect::<Vec<_>>();
        let mut affine_points = vec![C::identity(); 1 << max_k];
        C::Curve::batch_normalize(&points[..], &mut affine_points[..]);
        let points = affine_points;

        let scalars = (0..1 << max_k)
            .map(|_| C::Scalar::random(OsRng))
            .collect::<Vec<_>>();

        for k in min_k..=max_k {
            let points = &points[..1 << k];
            let scalars = &scalars[..1 << k];

            let t0 = start_timer!(|| format!("w/  booth k={}", k));
            let e0 = super::best_multiexp(scalars, points);
            end_timer!(t0);

            let t1 = start_timer!(|| format!("w/o booth k={}", k));
            let e1 = best_multiexp(scalars, points);
            end_timer!(t1);

            assert_eq!(e0, e1);
        }
    }

    #[test]
    fn test_msm_cross() {
        run_msm_cross::<G1Affine>(10, 18);
        // run_msm_cross::<G1Affine>(19, 23);
    }
}
