//! This benchmarks Multi Scalar Multiplication (MSM).
//! It measures `G1` from the BN256 curve.
//!
//! Benchmark with default feature `multicore` enabled:
//!
//!     cargo bench -- msm
//!
//! To run with as singlecore:
//!
//!     cargo bench --no-default-features -- msm

#[macro_use]
extern crate criterion;

use criterion::{black_box, BenchmarkId, Criterion};
use ff::Field;
use halo2curves::bn256::Fr as Scalar;
use halo2curves::bn256::G1Affine;
use halo2curves::msm::best_multiexp;
use halo2curves::CurveAffine;
use rand_core::OsRng;
use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::iter::zip;

fn random_curve_points<C: CurveAffine>(k: u8) -> Vec<G1Affine> {
    debug_assert!(k < 64);
    let n: u64 = 1 << k;

    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    (0..n).map(|_n| G1Affine::random(&mut rng)).collect()
}

#[cfg(not(feature = "multicore"))]
const RANGE: [u8; 6] = [3, 8, 10, 12 /*(Ethereum KZG / EIP 4844)*/, 14, 16];
#[cfg(feature = "multicore")]
const RANGE: [u8; 9] = [
    3, 8, 10, 12, /*(Ethereum KZG / EIP 4844)*/
    14, 16, 18, 20, 22,
];

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm");
    let rng = OsRng;
    for k in RANGE {
        group
            .bench_function(BenchmarkId::new("k", k), |b| {
                let mut g = random_curve_points::<G1Affine>(k);
                let half_len = g.len() / 2;
                let (g_lo, g_hi) = g.split_at_mut(half_len);
                let coeff_1 = Scalar::random(rng);
                let coeff_2 = Scalar::random(rng);

                b.iter(|| {
                    for (g_lo, g_hi) in zip(g_lo.iter(), g_hi.iter()) {
                        best_multiexp(&[black_box(coeff_1), black_box(coeff_2)], &[*g_lo, *g_hi]);
                    }
                })
            })
            .sample_size(10);
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
