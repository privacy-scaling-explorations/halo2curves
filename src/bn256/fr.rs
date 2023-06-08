#[cfg(feature = "maybe_u64")]
pub use super::fr_internal::Fr as FrInternal;
#[cfg(feature = "maybe_u64")]
pub type Fr = crate::MaybeU64<FrInternal>;

#[cfg(not(feature = "maybe_u64"))]
pub use super::fr_internal::Fr;

#[cfg(feature = "maybe_u64")]
impl crate::MaybeU64<FrInternal> {
    pub const fn from_raw(data: [u64; 4]) -> Self {
        if data[1] == 0 && data[2] == 0 && data[3] == 0 {
            Self::U64(data[0])
        } else {
            Self::Full(FrInternal::from_raw(data))
        }
    }

    pub const fn size() -> usize {
        FrInternal::size()
    }
}

#[cfg(feature = "maybe_u64")]
#[cfg(test)]
mod tests {
    use ark_std::{end_timer, start_timer, test_rng};
    use ff::Field;
    use rand_core::RngCore;

    use crate::{
        bn256::{Fr, FrInternal},
        MaybeU64,
    };

    const REPEAT: usize = 10000;

    #[test]
    fn bench_maybe_u64_mul() {
        let mut rng = test_rng();

        {
            let a64 = (0..REPEAT)
                .map(|_| rng.next_u32() as u64)
                .collect::<Vec<_>>();
            let b64 = (0..REPEAT)
                .map(|_| rng.next_u32() as u64)
                .collect::<Vec<_>>();

            let start = start_timer!(|| format!(
                "Native Fr with 32 bits {} muls with initialization",
                REPEAT
            ));
            let naitve_a = a64.iter().map(|x| FrInternal::from(*x)).collect::<Vec<_>>();
            let naitve_b = b64.iter().map(|x| FrInternal::from(*x)).collect::<Vec<_>>();
            let native_c = naitve_a
                .iter()
                .zip(naitve_b.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<_>>();
            end_timer!(start);

            let start = start_timer!(|| format!(
                "maybe64 Fr with 32 bits {} muls with initialization",
                REPEAT
            ));

            let maybe64_a = a64.iter().map(|x| Fr::from(*x)).collect::<Vec<_>>();
            let maybe64_b = b64.iter().map(|x| Fr::from(*x)).collect::<Vec<_>>();
            let maybe64_c = maybe64_a
                .iter()
                .zip(maybe64_b.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<_>>();
            end_timer!(start);

            native_c
                .iter()
                .zip(maybe64_c.iter())
                .for_each(|(a, b)| assert_eq!(MaybeU64::Full(*a), *b));
        }

        {
            let a64 = (0..REPEAT).map(|_| rng.next_u64()).collect::<Vec<_>>();
            let b64 = (0..REPEAT).map(|_| rng.next_u64()).collect::<Vec<_>>();

            let start = start_timer!(|| format!(
                "Native Fr with 64 bits {} muls with initialization",
                REPEAT
            ));
            let naitve_a = a64.iter().map(|x| FrInternal::from(*x)).collect::<Vec<_>>();
            let naitve_b = b64.iter().map(|x| FrInternal::from(*x)).collect::<Vec<_>>();
            let native_c = naitve_a
                .iter()
                .zip(naitve_b.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<_>>();
            end_timer!(start);

            let start = start_timer!(|| format!(
                "maybe64 Fr with 64 bits {} muls with initialization",
                REPEAT
            ));
            let maybe64_a = a64.iter().map(|x| Fr::from(*x)).collect::<Vec<_>>();
            let maybe64_b = b64.iter().map(|x| Fr::from(*x)).collect::<Vec<_>>();
            let maybe64_c = maybe64_a
                .iter()
                .zip(maybe64_b.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<_>>();
            end_timer!(start);

            native_c
                .iter()
                .zip(maybe64_c.iter())
                .for_each(|(a, b)| assert_eq!(MaybeU64::Full(*a), *b));
        }

        {
            let naitve_a = (0..REPEAT)
                .map(|_| FrInternal::random(&mut rng))
                .collect::<Vec<_>>();
            let naitve_b = (0..REPEAT)
                .map(|_| FrInternal::random(&mut rng))
                .collect::<Vec<_>>();
            let maybe64_a = naitve_a.iter().map(|x| Fr::Full(*x)).collect::<Vec<_>>();
            let maybe64_b = naitve_b.iter().map(|x| Fr::Full(*x)).collect::<Vec<_>>();

            let start = start_timer!(|| format!(
                "Native Fr with 256 bits {} muls without initialization",
                REPEAT
            ));
            let native_c = naitve_a
                .iter()
                .zip(naitve_b.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<_>>();
            end_timer!(start);

            let start = start_timer!(|| format!(
                "maybe64 Fr with 256 bits {} muls without initialization",
                REPEAT
            ));
            let maybe64_c = maybe64_a
                .iter()
                .zip(maybe64_b.iter())
                .map(|(a, b)| a * b)
                .collect::<Vec<_>>();
            end_timer!(start);

            native_c
                .iter()
                .zip(maybe64_c.iter())
                .for_each(|(a, b)| assert_eq!(MaybeU64::Full(*a), *b));
        }
    }
}
