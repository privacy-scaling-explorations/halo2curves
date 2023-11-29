//! This module provides "ZK Acceleration Layer" traits
//! to abstract away the execution engine for performance-critical primitives.
//!
//! The ZAL Engine is voluntarily left unconstrained
//! so that accelerator libraries are not prematurely limited.
//!
//! Terminology
//! -----------
//!
//! We use the name Backend+Engine for concrete implementations of ZalEngine.
//! For exaple H2cEngine for pure Halo2curves implementation.
//!
//! Alternative names considered were Executor or Driver however
//! - executor is already used in Rust (and the name is long)
//! - driver will be confusing as we work quite low-level with GPUs and FPGAs.
//!
//! Unfortunately Engine is used in bn256 for pairings.
//! Fortunately ZalEngine is only used in the prover
//! while "pairing engine" is only used in the verifier
//!
//! Initialization design space
//! ---------------------------
//!
//! It is recommended that ZAL backends provide:
//! - an initialization function:
//!   - either "fn new() -> ZalEngine" for simple libraries
//!   - or a builder pattern for complex initializations
//! - a shutdown function.
//!
//! The ZalEngine can be a stub type
//! and the shutdown function might be unnecessary
//! if the ZalEngine uses a global threadpool like Rayon.
//!
//! Backends might want to add as an option:
//! - The number of threads (CPU)
//! - The device(s) to run on (multi-sockets machines, multi-GPUs machines, ...)
//! - The curve (JIT-compiled backend)

use crate::msm::best_multiexp;
use pasta_curves::arithmetic::CurveAffine;

// The ZK Accel Layer API
// ---------------------------------------------------

pub(crate) trait ZalEngine{}

pub(crate) trait MsmAccel<C: CurveAffine>: ZalEngine {
    fn msm(&self, coeffs: &[C::Scalar], base: &[C]) -> C:: Curve;
}

// ZAL using Halo2curves as a backend
// ---------------------------------------------------

pub(crate) struct H2cEngine;

impl H2cEngine {
    pub(crate) fn new() -> Self {
        Self{}
    }
}

impl ZalEngine for H2cEngine{}

impl<C: CurveAffine> MsmAccel<C> for H2cEngine {
    fn msm(&self, coeffs: &[C::Scalar], bases: &[C]) -> C:: Curve {
        best_multiexp(coeffs, bases)
    }
}

// Testing
// ---------------------------------------------------

#[cfg(test)]
mod test {
    use crate::bn256::G1Affine;
    use ark_std::{end_timer, start_timer};
    use ff::Field;
    use group::{Curve, Group};
    use pasta_curves::arithmetic::CurveAffine;
    use rand_core::OsRng;
    use super::{H2cEngine, MsmAccel};

    fn run_msm_zal<C: CurveAffine>(min_k: usize, max_k: usize) {
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

            let t0 = start_timer!(|| format!("freestanding msm k={}", k));
            let e0 = super::best_multiexp(scalars, points);
            end_timer!(t0);

            let engine = H2cEngine::new();
            let t1 = start_timer!(|| format!("H2cEngine msm k={}", k));
            let e1 = engine.msm(scalars, points);
            end_timer!(t1);

            assert_eq!(e0, e1);
        }
    }

    #[test]
    fn test_msm_zal() {
        run_msm_zal::<G1Affine>(3, 14);
    }
}