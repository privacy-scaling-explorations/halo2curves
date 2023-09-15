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
use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;

const SEED: [u8; 16] = [
    0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc, 0xe5,
];

const SINGLECORE_RANGE: [u8; 6] = [3, 8, 10, 12, 14, 16];

const MULTICORE_RANGE: [u8; 9] = [3, 8, 10, 12, 14, 16, 18, 20, 22];

fn singlecore(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm/singlecore");
    let mut rng = XorShiftRng::from_seed(SEED);
    for k in SINGLECORE_RANGE {
        group
            .bench_function(BenchmarkId::new("k", k), |b| {
                assert!(k < 64);
                let n: u64 = 1 << k;

                let bases: Vec<_> = (0..n).map(|_| Point::random(&mut rng)).collect();
                let coeffs: Vec<_> = (0..n).map(|_| Scalar::random(&mut rng)).collect();
                let mut acc = Point::identity().into();

                b.iter(|| multiexp_serial(&coeffs, &bases, &mut black_box(acc)));
            })
            .sample_size(10);
    }
}

fn multicore(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm/multicore");
    let mut rng = XorShiftRng::from_seed(SEED);
    for k in MULTICORE_RANGE {
        group
            .bench_function(BenchmarkId::new("k", k), |b| {
                assert!(k < 64);
                let n: u64 = 1 << k;

                let bases: Vec<_> = (0..n).map(|_| Point::random(&mut rng)).collect();
                let coeffs: Vec<_> = (0..n).map(|_| Scalar::random(&mut rng)).collect();

                b.iter(|| {
                    best_multiexp(&coeffs, &bases);
                })
            })
            .sample_size(10);
    }
}

criterion_group!(benches, singlecore, multicore);
criterion_main!(benches);
