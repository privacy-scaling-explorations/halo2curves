//! This benchmarks Multi Scalar Multiplication (MSM).
//! It measures `G1` from the BN256 curve.
//!
//! To run this benchmark:
//!
//!     cargo bench -- msm
//!
//! Caveat:  `multicore` should be read as _allowing_ for multicore computation --
//! not enforcing it.
//!

#[macro_use]
extern crate criterion;

use criterion::{black_box, BenchmarkId, Criterion};
use ff::Field;
use group::prime::PrimeCurveAffine;
use halo2curves::bn256::{Fr as Scalar, G1Affine as Point};
use halo2curves::msm::{best_multiexp, multiexp_serial};
use maybe_rayon::current_thread_index;
use maybe_rayon::prelude::{IntoParallelIterator, ParallelIterator};
use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::time::{Duration, SystemTime};

const SEED: [u8; 16] = [
    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc, 0xe5,
];

const SINGLECORE_RANGE: [u8; 6] = [3, 8, 10, 12, 14, 16];

const MULTICORE_RANGE: [u8; 9] = [3, 8, 10, 12, 14, 16, 18, 20, 22];

/// This do get called twice, but the total running time entirely dominated by the larger instance.
fn get_data(k: u8) -> (Vec<Scalar>, Vec<Point>) {
    let n: u64 = {
        assert!(k < 64);
        1 << k
    };

    println!(
        "\n\nCoefficient and curve point generation starting.  {} coefficient-points pairs needed",
        n
    );
    let timer = SystemTime::now();
    let coeffs = (0..n)
        .into_par_iter()
        .map_init(
            || {
                let mut thread_seed = SEED.clone();
                let uniq = current_thread_index().unwrap().to_ne_bytes();
                assert!(std::mem::size_of::<usize>() == 8);
                for i in 0..uniq.len() {
                    thread_seed[i] += uniq[i];
                    thread_seed[i + 8] += uniq[i];
                }
                XorShiftRng::from_seed(thread_seed)
            },
            |mut rng, _| Scalar::random(&mut rng),
        )
        .collect();
    let bases = (0..n)
        .into_par_iter()
        .map_init(
            || {
                let mut thread_seed = SEED.clone();
                let uniq = current_thread_index().unwrap().to_ne_bytes();
                assert!(std::mem::size_of::<usize>() == 8);
                for i in 0..uniq.len() {
                    thread_seed[i] += uniq[i];
                    thread_seed[i + 8] += uniq[i];
                }
                XorShiftRng::from_seed(thread_seed)
            },
            |mut rng, _| Point::random(&mut rng),
        )
        .collect();
    let end = timer.elapsed().unwrap();
    println!(
        "Coefficient and curve point generation took: {} sec.\n\n",
        end.as_secs()
    );

    return (coeffs, bases);
}

fn singlecore(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm/singlecore");
    let (coeffs, bases) = get_data(*SINGLECORE_RANGE.iter().max().unwrap());
    for k in SINGLECORE_RANGE {
        group
            .bench_function(BenchmarkId::new("k", k), |b| {
                assert!(k < 64);
                let n: usize = 1 << k;

                let mut acc = Point::identity().into();

                b.iter(|| multiexp_serial(&coeffs[..n], &bases[..n], &mut black_box(acc)));
            })
            .sample_size(10);
    }
}

fn multicore(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm/multicore");
    let (coeffs, bases) = get_data(*MULTICORE_RANGE.iter().max().unwrap());
    for k in MULTICORE_RANGE {
        group
            .bench_function(BenchmarkId::new("k", k), |b| {
                assert!(k < 64);
                let n: usize = 1 << k;

                b.iter(|| {
                    best_multiexp(&coeffs[..n], &bases[..n]);
                })
            })
            .sample_size(10);
    }
}

criterion_group!(benches, singlecore, multicore);
criterion_main!(benches);
