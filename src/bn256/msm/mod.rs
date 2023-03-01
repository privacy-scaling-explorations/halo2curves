use super::{Fr, G1Affine};
use crate::arithmetic::msm_zcash;
use crate::bn256::{msm::round::Round, G1};
use crate::group::Group;
use ff::PrimeField;
use rayon::{current_num_threads, scope};

#[macro_export]
macro_rules! div_ceil {
    ($a:expr, $b:expr) => {
        (($a - 1) / $b) + 1
    };
}

#[macro_export]
macro_rules! double_n {
    ($acc:expr, $n:expr) => {
        (0..$n).fold($acc, |acc, _| acc.double())
    };
}

#[macro_export]
macro_rules! range {
    ($index:expr, $n_items:expr) => {
        $index * $n_items..($index + 1) * $n_items
    };
}

#[macro_export]
macro_rules! index {
    ($digit:expr) => {
        ($digit & 0x7fffffff) as usize
    };
}

#[macro_export]
macro_rules! is_neg {
    ($digit:expr) => {
        sign_bit!($digit) != 0
    };
}

#[macro_export]
macro_rules! sign_bit {
    ($digit:expr) => {
        $digit & 0x80000000
    };
}

mod round;

pub struct MSM {
    signed_digits: Vec<u32>,
    sorted_positions: Vec<u32>,
    bucket_sizes: Vec<usize>,
    bucket_offsets: Vec<usize>,
    n_windows: usize,
    window: usize,
    n_buckets: usize,
    n_points: usize,
    round: Round,
}

impl MSM {
    pub fn allocate(n_points: usize, override_window: Option<usize>) -> Self {
        fn best_window(n: usize) -> usize {
            if n >= 262144 {
                15
            } else if n >= 65536 {
                12
            } else if n >= 16384 {
                11
            } else if n >= 8192 {
                10
            } else if n >= 1024 {
                9
            } else {
                7
            }
        }
        let window = match override_window {
            Some(window) => {
                let overriden = best_window(n_points);
                println!("override window from {} to {}", overriden, window);
                window
            }
            None => best_window(n_points),
        };
        let n_windows = div_ceil!(Fr::NUM_BITS as usize, window);
        let n_buckets = (1 << (window - 1)) + 1;
        let round = Round::new(n_buckets, n_points);
        MSM {
            signed_digits: vec![0u32; n_windows * n_points],
            sorted_positions: vec![0u32; n_windows * n_points],
            bucket_sizes: vec![0usize; n_windows * n_buckets],
            bucket_offsets: vec![0; n_buckets],
            n_windows,
            window,
            n_buckets,
            n_points,
            round,
        }
    }

    fn decompose(&mut self, scalars: &[Fr]) {
        pub(crate) fn get_bits(segment: usize, window: usize, bytes: &[u8]) -> u32 {
            let skip_bits = segment * window;
            let skip_bytes = skip_bits / 8;
            if skip_bytes >= 32 {
                return 0;
            }
            let mut v = [0; 4];
            for (v, o) in v.iter_mut().zip(bytes[skip_bytes..].iter()) {
                *v = *o;
            }
            let mut tmp = u32::from_le_bytes(v);
            tmp >>= skip_bits - (skip_bytes * 8);
            tmp %= 1 << window;
            tmp
        }
        let max = 1 << (self.window - 1);
        for (point_idx, scalar) in scalars.iter().enumerate() {
            let repr = scalar.to_repr();
            let mut borrow = 0u32;
            for window_idx in 0..self.n_windows {
                let windowed_digit = get_bits(window_idx, self.window, repr.as_ref()) + borrow;
                let signed_digit = if windowed_digit >= max {
                    borrow = 1;
                    ((1 << self.window) - windowed_digit) | 0x80000000
                } else {
                    borrow = 0;
                    windowed_digit
                };
                self.bucket_sizes[window_idx * self.n_buckets + index!(signed_digit)] += 1;
                self.signed_digits[window_idx * self.n_points + point_idx] = signed_digit;
            }
        }
        self.sort();
    }

    fn sort(&mut self) {
        for w_i in 0..self.n_windows {
            let sorted_positions = &mut self.sorted_positions[range!(w_i, self.n_points)];
            let bucket_sizes = &self.bucket_sizes[range!(w_i, self.n_buckets)];
            let signed_digits = &self.signed_digits[range!(w_i, self.n_points)];
            let mut offset = 0;
            for (i, size) in bucket_sizes.iter().enumerate() {
                self.bucket_offsets[i] = offset;
                offset += size;
            }
            for (sorted_idx, signed_digit) in signed_digits.iter().enumerate() {
                let bucket_idx = index!(signed_digit);
                sorted_positions[self.bucket_offsets[bucket_idx]] =
                    sign_bit!(signed_digit) | (sorted_idx as u32);
                self.bucket_offsets[bucket_idx] += 1;
            }
        }
    }

    pub fn evaluate_with(
        scalars: &[Fr],
        points: &[G1Affine],
        acc: &mut G1,
        override_window: Option<usize>,
    ) {
        let mut msm = Self::allocate(points.len(), override_window);
        msm.decompose(scalars);
        for w_i in (0..msm.n_windows).rev() {
            if w_i != msm.n_windows - 1 {
                *acc = double_n!(*acc, msm.window);
            }
            msm.round.init(
                points,
                &msm.sorted_positions[range!(w_i, msm.n_points)],
                &msm.bucket_sizes[range!(w_i, msm.n_buckets)],
            );
            let buckets = msm.round.evaluate();
            let mut running_sum = G1::identity();
            for bucket in buckets.into_iter().skip(1).rev() {
                running_sum += bucket;
                *acc += &running_sum;
            }
        }
    }

    pub fn evaluate(scalars: &[Fr], points: &[G1Affine], acc: &mut G1) {
        Self::evaluate_with(scalars, points, acc, None)
    }

    pub fn best(scalars: &[Fr], points: &[G1Affine]) -> G1 {
        assert_eq!(scalars.len(), points.len());
        let num_threads = current_num_threads();
        if scalars.len() > num_threads {
            let chunk = scalars.len() / num_threads;
            let num_chunks = scalars.chunks(chunk).len();
            let mut results = vec![G1::identity(); num_chunks];
            scope(|scope| {
                let chunk = scalars.len() / num_threads;

                for ((scalars, points), acc) in scalars
                    .chunks(chunk)
                    .zip(points.chunks(chunk))
                    .zip(results.iter_mut())
                {
                    scope.spawn(move |_| {
                        if points.len() < 1 << 8 {
                            msm_zcash(scalars, points, acc);
                        } else {
                            Self::evaluate(scalars, points, acc);
                        }
                    });
                }
            });
            results.iter().fold(G1::identity(), |a, b| a + b)
        } else {
            let mut acc = G1::identity();
            Self::evaluate(scalars, points, &mut acc);
            acc
        }
    }
}

#[cfg(test)]
mod test {
    use crate::arithmetic::msm_zcash;
    use crate::bn256::{Fr, G1Affine, G1};
    use crate::group::Group;
    use crate::serde::SerdeObject;
    use ff::Field;
    use group::Curve;
    use rand::Rng;
    use rand_core::OsRng;
    use std::fs::File;
    use std::path::Path;

    fn read_data(n: usize) -> (Vec<G1Affine>, Vec<Fr>) {
        let mut file = File::open("data.tmp").unwrap();
        (0..n)
            .map(|_| {
                let point = G1Affine::read_raw_unchecked(&mut file);
                let scalar = Fr::read_raw_unchecked(&mut file);
                (point, scalar)
            })
            .unzip()
    }

    pub(crate) fn get_data(n: usize) -> (Vec<G1Affine>, Vec<Fr>) {
        const MAX_N: usize = 1 << 22;
        assert!(n <= MAX_N);
        if Path::new("data.tmp").is_file() {
            read_data(n)
        } else {
            let mut file = File::create("data.tmp").unwrap();
            (0..MAX_N)
                .map(|_| {
                    let point = G1::random(OsRng).to_affine();
                    let scalar = Fr::random(OsRng);
                    point.write_raw(&mut file).unwrap();
                    scalar.write_raw(&mut file).unwrap();
                    (point, scalar)
                })
                .take(n)
                .unzip()
        }
    }

    #[test]
    fn test_msm() {
        let (min_k, max_k) = (4, 20);
        let (points, scalars) = get_data(1 << max_k);

        for k in min_k..=max_k {
            let mut rng = OsRng;
            let n_points = rng.gen_range(1 << (k - 1)..1 << k);
            let scalars = &scalars[..n_points];
            let points = &points[..n_points];
            let mut r0 = G1::identity();
            msm_zcash(scalars, points, &mut r0);
            let mut r1 = G1::identity();
            super::MSM::evaluate(scalars, points, &mut r1);
            assert_eq!(r0, r1);
        }
    }
}
