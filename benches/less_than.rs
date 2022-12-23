#![allow(unused)]

use criterion::BenchmarkId;

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

#[inline(always)]
fn is_less_than(x: &[u64; 4], y: &[u64; 4]) -> bool {
    match x[3].cmp(&y[3]) {
        core::cmp::Ordering::Less => return true,
        core::cmp::Ordering::Greater => return false,
        _ => {}
    }
    match x[2].cmp(&y[2]) {
        core::cmp::Ordering::Less => return true,
        core::cmp::Ordering::Greater => return false,
        _ => {}
    }
    match x[1].cmp(&y[1]) {
        core::cmp::Ordering::Less => return true,
        core::cmp::Ordering::Greater => return false,
        _ => {}
    }
    x[0].lt(&y[0])
}

#[inline(always)]
fn check_underflow(x: &[u64; 4], y: &[u64; 4]) -> bool {
    let (_, borrow) = sbb(x[0], y[0], 0);
    let (_, borrow) = sbb(x[1], y[1], borrow);
    let (_, borrow) = sbb(x[2], y[2], borrow);
    let (_, borrow) = sbb(x[3], y[3], borrow);
    borrow >> 63 == 1
}

use criterion::{criterion_group, criterion_main, Criterion};

pub fn criterion_benchmark(c: &mut Criterion) {
    let x: [u64; 4] = [(); 4].map(|_| rand::random());
    let y: [u64; 4] = [(); 4].map(|_| rand::random());

    let mut group = c.benchmark_group("Big less than methods");

    group.bench_with_input(
        BenchmarkId::new("is_less_than", ""),
        &(x, y),
        |b, (x, y)| b.iter(|| is_less_than(x, y)),
    );

    group.bench_with_input(
        BenchmarkId::new("check_underflow", ""),
        &(x, y),
        |b, (x, y)| b.iter(|| check_underflow(x, y)),
    );
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
