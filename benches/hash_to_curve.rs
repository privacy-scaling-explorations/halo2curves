use criterion::{black_box, criterion_group, criterion_main, Criterion};
use halo2curves::bn256::G1;
use halo2curves::secp256k1::Secp256k1;
use pasta_curves::arithmetic::CurveExt;
use rand_core::{OsRng, RngCore};
use std::iter;

fn hash_to_secp256k1(c: &mut Criterion) {
    hash_to_curve::<Secp256k1>(c, "Secp256k1");
}

fn hash_to_bn256(c: &mut Criterion) {
    hash_to_curve::<G1>(c, "Bn256");
}

fn hash_to_curve<G: CurveExt>(c: &mut Criterion, name: &'static str) {
    {
        let hasher = G::hash_to_curve("test");
        let mut rng = OsRng;
        let message = iter::repeat_with(|| rng.next_u32().to_be_bytes())
            .take(1024)
            .flatten()
            .collect::<Vec<_>>();

        c.bench_function(&format!("Hash to {}", name), move |b| {
            b.iter(|| hasher(black_box(&message)))
        });
    }
}

criterion_group!(benches, hash_to_secp256k1, hash_to_bn256);
criterion_main!(benches);
