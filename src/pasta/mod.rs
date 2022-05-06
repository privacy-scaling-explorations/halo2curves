pub use pasta_curves::{pallas, vesta, Ep, EpAffine, Eq, EqAffine, Fp, Fq};

impl crate::CurveAffine for EpAffine {
    fn batch_add<const COMPLETE: bool, const LOAD_POINTS: bool>(
        _: &mut [Self],
        _: &[u32],
        _: usize,
        _: usize,
        _: &[Self],
        _: &[u32],
    ) {
        unimplemented!();
    }
}

impl crate::CurveAffine for EqAffine {
    fn batch_add<const COMPLETE: bool, const LOAD_POINTS: bool>(
        _: &mut [Self],
        _: &[u32],
        _: usize,
        _: usize,
        _: &[Self],
        _: &[u32],
    ) {
        unimplemented!();
    }
}
