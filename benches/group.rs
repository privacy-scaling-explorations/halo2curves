use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ff::Field;
use group::prime::PrimeCurveAffine;
use halo2curves::pluto_eris::Pluto;
use pasta_curves::arithmetic::CurveExt;
use rand_core::OsRng;

fn criterion_benchmark<G: CurveExt>(c: &mut Criterion) {
    // G1Projective
    {
        let name = "GProjective";
        let p1 = G::random(OsRng);
        let p2 = G::random(OsRng);
        let p1_affine = G::AffineExt::from(p1);
        let s = G::ScalarExt::random(OsRng);

        const N: usize = 1000;
        let v = vec![G::generator(); N];
        let mut q = vec![G::AffineExt::identity(); N];

        c.bench_function(&format!("{} check on curve", name), move |b| {
            b.iter(|| black_box(p1).is_on_curve())
        });
        c.bench_function(&format!("{} check equality", name), move |b| {
            b.iter(|| black_box(p1) == black_box(p1))
        });
        c.bench_function(&format!("{} to affine", name), move |b| {
            b.iter(|| G::AffineExt::from(black_box(p1)))
        });
        c.bench_function(&format!("{} doubling", name), move |b| {
            b.iter(|| black_box(p1).double())
        });
        c.bench_function(&format!("{} addition", name), move |b| {
            b.iter(|| black_box(p1).add(&p2))
        });
        c.bench_function(&format!("{} mixed addition", name), move |b| {
            b.iter(|| black_box(p2).add(&p1_affine))
        });
        c.bench_function(&format!("{} scalar multiplication", name), move |b| {
            b.iter(|| black_box(p1) * black_box(s))
        });
        c.bench_function(&format!("{} batch to affine n={}", name, N), move |b| {
            b.iter(|| {
                G::batch_normalize(black_box(&v), black_box(&mut q));
                black_box(&q)[0]
            })
        });
    }
}

criterion_group!(benches, criterion_benchmark<Pluto>);
criterion_main!(benches);
