//! This benchmark allows testing msm without depending on the `halo2_proofs`
//! crate.  This code originates in an older version of `halo2_proofs` from
//! before the `hash_to_curve` method was implemented.  It currently only uses
//! curve `Secp256k1Affine`

#[macro_use]
extern crate criterion;

use criterion::{black_box, BenchmarkId, Criterion};
use ff::Field;
use halo2_proofs::arithmetic::small_multiexp;
use halo2curves::secp256k1::Fq as Scalar;
use halo2curves::secp256k1::Secp256k1Affine;
use halo2curves::CurveAffine;
use rand_core::OsRng;
use rand_core::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::iter::zip;

fn random_curve_points<C: CurveAffine>(k: u8) -> Vec<Secp256k1Affine> {
    debug_assert!(k < 64);
    let n: u64 = 1 << k;

    let mut rng = XorShiftRng::from_seed([
        0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
        0xe5,
    ]);

    (0..n).map(|_n| Secp256k1Affine::random(&mut rng)).collect()
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("msm-alt");
    for k in 8..16 {
        group
            .bench_function(BenchmarkId::new("k", k), |b| {
                let rng = OsRng;

                let mut g = random_curve_points::<Secp256k1Affine>(k);
                let half_len = g.len() / 2;
                let (g_lo, g_hi) = g.split_at_mut(half_len);
                let coeff_1 = Scalar::random(rng);
                let coeff_2 = Scalar::random(rng);

                b.iter(|| {
                    for (g_lo, g_hi) in zip(g_lo.iter(), g_hi.iter()) {
                        small_multiexp(&[black_box(coeff_1), black_box(coeff_2)], &[*g_lo, *g_hi]);
                    }
                })
            })
            .sample_size(30);
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
