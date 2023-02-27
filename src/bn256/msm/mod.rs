use super::{Fr, G1Affine};
use crate::bn256::{msm::round::Round, G1};
use crate::group::Group;
use ff::PrimeField;
use rayon::{current_num_threads, scope};

#[cfg(test)]
mod pr40;
mod round;
#[cfg(test)]
mod zcash;

macro_rules! div_ceil {
    ($a:expr, $b:expr) => {
        (($a - 1) / $b) + 1
    };
}

macro_rules! double_n {
    ($acc:expr, $n:expr) => {
        (0..$n).fold($acc, |acc, _| acc.double())
    };
}

macro_rules! range {
    ($index:expr, $n_items:expr) => {
        $index * $n_items..($index + 1) * $n_items
    };
}

pub struct MSM {
    bucket_sizes: Vec<usize>,
    sorted_positions: Vec<usize>,
    bucket_indexes: Vec<usize>,
    bucket_offsets: Vec<usize>,
    n_windows: usize,
    window: usize,
    n_buckets: usize,
    n_points: usize,
    round: Round,
}

impl MSM {
    pub fn alloacate(n_points: usize) -> Self {
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
        let window = best_window(n_points);
        let n_windows = div_ceil!(Fr::NUM_BITS as usize, window);
        let n_buckets = 1 << window;
        let round = Round::new(n_buckets, n_points);
        MSM {
            bucket_indexes: vec![0usize; n_windows * n_points],
            bucket_sizes: vec![0usize; n_windows * n_buckets],
            sorted_positions: vec![0usize; n_windows * n_points],
            bucket_offsets: vec![0; n_buckets],
            n_windows,
            window,
            n_buckets,
            n_points,
            round,
        }
    }

    fn decompose(&mut self, scalars: &[Fr]) {
        pub(crate) fn get_bits(segment: usize, c: usize, bytes: &[u8]) -> u64 {
            let skip_bits = segment * c;
            let skip_bytes = skip_bits / 8;
            if skip_bytes >= 32 {
                return 0;
            }
            let mut v = [0; 8];
            for (v, o) in v.iter_mut().zip(bytes[skip_bytes..].iter()) {
                *v = *o;
            }
            let mut tmp = u64::from_le_bytes(v);
            tmp >>= skip_bits - (skip_bytes * 8);
            tmp %= 1 << c;
            tmp as u64
        }
        let scalars = scalars
            .iter()
            .map(|scalar| scalar.to_repr())
            .collect::<Vec<_>>();
        for window_idx in 0..self.n_windows {
            for (point_index, scalar) in scalars.iter().enumerate() {
                let bucket_index = get_bits(window_idx, self.window, scalar.as_ref()) as usize;
                self.bucket_sizes[window_idx * self.n_buckets + bucket_index] += 1;
                self.bucket_indexes[window_idx * self.n_points + point_index] = bucket_index;
            }
        }
        self.sort();
    }

    fn sort(&mut self) {
        for w_i in 0..self.n_windows {
            let sorted_positions = &mut self.sorted_positions[range!(w_i, self.n_points)];
            let bucket_sizes = &self.bucket_sizes[range!(w_i, self.n_buckets)];
            let bucket_indexes = &self.bucket_indexes[range!(w_i, self.n_points)];
            let mut offset = 0;
            for (i, size) in bucket_sizes.iter().enumerate() {
                self.bucket_offsets[i] = offset;
                offset += size;
            }
            for (idx, bucket_index) in bucket_indexes.iter().enumerate() {
                sorted_positions[self.bucket_offsets[*bucket_index]] = idx;
                self.bucket_offsets[*bucket_index] += 1;
            }
        }
    }

    pub fn evalulate(scalars: &[Fr], bases: &[G1Affine], acc: &mut G1) {
        let mut msm = Self::alloacate(bases.len());
        msm.decompose(scalars);
        for w_i in (0..msm.n_windows).rev() {
            if w_i != msm.n_windows - 1 {
                *acc = double_n!(*acc, msm.window);
            }
            msm.round.init(
                bases,
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

    pub fn best(coeffs: &[Fr], bases: &[G1Affine]) -> G1 {
        assert_eq!(coeffs.len(), bases.len());
        let num_threads = current_num_threads();
        if coeffs.len() > num_threads {
            let chunk = coeffs.len() / num_threads;
            let num_chunks = coeffs.chunks(chunk).len();
            let mut results = vec![G1::identity(); num_chunks];
            scope(|scope| {
                let chunk = coeffs.len() / num_threads;

                for ((coeffs, bases), acc) in coeffs
                    .chunks(chunk)
                    .zip(bases.chunks(chunk))
                    .zip(results.iter_mut())
                {
                    scope.spawn(move |_| {
                        Self::evalulate(coeffs, bases, acc);
                    });
                }
            });
            results.iter().fold(G1::identity(), |a, b| a + b)
        } else {
            let mut acc = G1::identity();
            Self::evalulate(coeffs, bases, &mut acc);
            acc
        }
    }
}

#[cfg(test)]
mod test {
    use crate::bn256::msm::pr40::{MultiExp, MultiExpContext};
    use crate::bn256::msm::zcash::{best_multiexp_zcash, msm_zcash};
    use crate::bn256::{Fr, G1Affine, G1};
    use crate::group::Group;
    use crate::serde::SerdeObject;
    use ff::Field;
    use group::Curve;
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
        let (points, scalars) = get_data(1 << 22);

        for k in 10..=22 {
            let n_points = 1 << k;
            let scalars = &scalars[..n_points];
            let points = &points[..n_points];
            println!("------ {}", k);

            let mut r0 = G1::identity();
            let time = std::time::Instant::now();
            msm_zcash(scalars, points, &mut r0);
            println!("zcash serial {:?}", time.elapsed());

            let time = std::time::Instant::now();
            let r0 = best_multiexp_zcash(scalars, points);
            println!("zcash parallel {:?}", time.elapsed());

            let time = std::time::Instant::now();
            let mut r1 = G1::identity();
            super::MSM::evalulate(scalars, points, &mut r1);
            assert_eq!(r0, r1);
            println!("this {:?}", time.elapsed());

            let time = std::time::Instant::now();
            let r1 = super::MSM::best(scalars, points);
            assert_eq!(r0, r1);
            println!("this parallel {:?}", time.elapsed());

            let time = std::time::Instant::now();
            let msm = MultiExp::new(&points);
            let mut ctx = MultiExpContext::default();
            let _ = msm.evaluate(&mut ctx, scalars, false);
            // assert_eq!(r0, r1); // fails
            println!("pr40 {:?}", time.elapsed());
        }
    }
}
