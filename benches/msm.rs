//! This benchmarks Multi Scalar Multiplication (MSM).
//! It measures `G1` from the BN256 curve.
//!
//! To run this benchmark:
//!
//!     cargo bench -- msm
//!
//! Caveat:  The multicore benchmark assumes:
//!     1. a multi-core system
//!     2. that the `multicore` feature is enabled.  It is by default.

#[macro_use]
extern crate criterion;

use criterion::{BenchmarkId, Criterion};
use ff::Field;
use ff::PrimeField;
use group::prime::PrimeCurveAffine;
use halo2curves::bn256::{Fr as Scalar, G1Affine as Point};
use halo2curves::msm::{
    best_multiexp, best_multiexp_jonathan, best_multiexp_skip_zeros, multiexp_serial,
};
use rand_core::RngCore;
use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;
use rayon::current_thread_index;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use std::time::SystemTime;

const SAMPLE_SIZE: usize = 10;
const SINGLECORE_RANGE: [u8; 6] = [3, 8, 10, 12, 14, 16];
const MULTICORE_RANGE: [u8; 9] = [3, 8, 10, 12, 14, 16, 18, 20, 22];
const SEED: [u8; 16] = [
    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc, 0xe5,
];

fn generate_coefficients_and_curvepoints(k: u8) -> (Vec<Scalar>, Vec<Point>) {
    let n: u64 = {
        assert!(k < 64);
        1 << k
    };

    println!("\n\nGenerating 2^{k} = {n} coefficients and curve points..",);
    let timer = SystemTime::now();
    let coeffs = (0..n)
        .into_par_iter()
        .map_init(
            || {
                let mut thread_seed = SEED;
                let uniq = current_thread_index().unwrap().to_ne_bytes();
                assert!(std::mem::size_of::<usize>() == 8);
                for i in 0..uniq.len() {
                    thread_seed[i] += uniq[i];
                    thread_seed[i + 8] += uniq[i];
                }
                XorShiftRng::from_seed(thread_seed)
            },
            |rng, _| Scalar::random(rng),
        )
        .collect();
    let bases = (0..n)
        .into_par_iter()
        .map_init(
            || {
                let mut thread_seed = SEED;
                let uniq = current_thread_index().unwrap().to_ne_bytes();
                assert!(std::mem::size_of::<usize>() == 8);
                for i in 0..uniq.len() {
                    thread_seed[i] += uniq[i];
                    thread_seed[i + 8] += uniq[i];
                }
                XorShiftRng::from_seed(thread_seed)
            },
            |rng, _| Point::random(rng),
        )
        .collect();
    let end = timer.elapsed().unwrap();
    println!(
        "Generating 2^{k} = {n} coefficients and curve points took: {} sec.\n\n",
        end.as_secs()
    );

    (coeffs, bases)
}

fn msm(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm");
    let max_k = *SINGLECORE_RANGE
        .iter()
        .chain(MULTICORE_RANGE.iter())
        .max()
        .unwrap_or(&16);
    let (coeffs, bases) = generate_coefficients_and_curvepoints(max_k);

    for k in SINGLECORE_RANGE {
        group
            .bench_function(BenchmarkId::new("singlecore", k), |b| {
                assert!(k < 64);
                let n: usize = 1 << k;
                let mut acc = Point::identity().into();
                b.iter(|| multiexp_serial(&coeffs[..n], &bases[..n], &mut acc));
            })
            .sample_size(10);
    }
    for k in MULTICORE_RANGE {
        group
            .bench_function(BenchmarkId::new("multicore", k), |b| {
                assert!(k < 64);
                let n: usize = 1 << k;
                b.iter(|| {
                    best_multiexp(&coeffs[..n], &bases[..n]);
                })
            })
            .sample_size(SAMPLE_SIZE);
    }
    group.finish();
}

fn gen_scalars_points(k: u8, small: bool) -> (Vec<Scalar>, Vec<Point>) {
    let points = (0..1 << k)
        .into_par_iter()
        .map_init(
            || {
                let uniq = current_thread_index().unwrap();
                assert!(std::mem::size_of::<usize>() == 8);
                XorShiftRng::seed_from_u64(uniq as u64)
            },
            |rng, _| Point::random(rng),
        )
        .collect();

    // 1 byte upper bound
    let max_val = 2u64.pow((8) as u32);

    let scalars = (0..1 << k)
        .into_par_iter()
        .map_init(
            || {
                let uniq = current_thread_index().unwrap();
                assert!(std::mem::size_of::<usize>() == 8);
                XorShiftRng::seed_from_u64(uniq as u64)
            },
            |rng, _| {
                if small {
                    let v = rng.next_u64() % max_val;
                    Scalar::from_u128(v as u128)
                } else {
                    Scalar::random(rng)
                }
            },
        )
        .collect();

    (scalars, points)
}

fn msm_cmp(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm_cmp");
    let min_k = 18;
    let max_k = 22;
    let (scalars_small, points_small) = gen_scalars_points(max_k, true);
    let (scalars_big, points_big) = gen_scalars_points(max_k, false);

    for small in [false, true] {
        let (scalars, points) = if small {
            (&scalars_small, &points_small)
        } else {
            (&scalars_big, &points_big)
        };
        for k in min_k..=max_k {
            let name = format!("msm func={}, k={}, small={}", "original", k, small);
            group
                .bench_function(BenchmarkId::new(name, k), |b| {
                    let n: usize = 1 << k;
                    b.iter(|| {
                        best_multiexp(&scalars[..n], &points[..n]);
                    })
                })
                .sample_size(10);

            let name = format!("msm func={}, k={}, small={}", "skip_zeros_edu", k, small);
            group
                .bench_function(BenchmarkId::new(name, k), |b| {
                    let n: usize = 1 << k;
                    b.iter(|| {
                        best_multiexp_skip_zeros(&scalars[..n], &points[..n]);
                    })
                })
                .sample_size(10);

            let name = format!(
                "msm func={}, k={}, small={}",
                "skip_zeros_jonathan", k, small
            );
            group
                .bench_function(BenchmarkId::new(name, k), |b| {
                    let n: usize = 1 << k;
                    b.iter(|| {
                        best_multiexp_jonathan(&scalars[..n], &points[..n]);
                    })
                })
                .sample_size(10);
        }
    }
    group.finish();
}

criterion_group!(benches, msm_cmp);
criterion_main!(benches);
