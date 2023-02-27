//! This module provides common utilities, traits and structures for group and
//! field arithmetic.
//!
//! This module is temporary, and the extension traits defined here are expected to be
//! upstreamed into the `ff` and `group` crates after some refactoring.

use pasta_curves::arithmetic::CurveAffine;

use crate::CurveExt;

pub(crate) struct EndoParameters {
    pub(crate) gamma1: [u64; 4],
    pub(crate) gamma2: [u64; 4],
    pub(crate) b1: [u64; 4],
    pub(crate) b2: [u64; 4],
}

pub trait CurveEndo: CurveExt {
    fn decompose_scalar(e: &Self::ScalarExt) -> (u128, bool, u128, bool);
}

pub trait CurveAffineExt: CurveAffine {
    fn decompose_scalar(k: &Self::ScalarExt) -> (u128, bool, u128, bool);
    fn endo(&self) -> Self;
    fn batch_add<const COMPLETE: bool, const LOAD_POINTS: bool>(
        points: &mut [Self],
        output_indices: &[u32],
        num_points: usize,
        offset: usize,
        bases: &[Self],
        base_positions: &[u32],
    );
}

/// Compute a + b + carry, returning the result and the new carry over.
#[inline(always)]
pub(crate) const fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + (b as u128) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a - (b + borrow), returning the result and the new borrow.
#[inline(always)]
pub(crate) const fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let ret = (a as u128).wrapping_sub((b as u128) + ((borrow >> 63) as u128));
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c) + carry, returning the result and the new carry over.
#[inline(always)]
pub(crate) const fn mac(a: u64, b: u64, c: u64, carry: u64) -> (u64, u64) {
    let ret = (a as u128) + ((b as u128) * (c as u128)) + (carry as u128);
    (ret as u64, (ret >> 64) as u64)
}

/// Compute a + (b * c), returning the result and the new carry over.
#[inline(always)]
pub(crate) const fn macx(a: u64, b: u64, c: u64) -> (u64, u64) {
    let res = (a as u128) + ((b as u128) * (c as u128));
    (res as u64, (res >> 64) as u64)
}

/// Compute a * b, returning the result.
#[inline(always)]
pub(crate) fn mul_512(a: [u64; 4], b: [u64; 4]) -> [u64; 8] {
    let (r0, carry) = macx(0, a[0], b[0]);
    let (r1, carry) = macx(carry, a[0], b[1]);
    let (r2, carry) = macx(carry, a[0], b[2]);
    let (r3, carry_out) = macx(carry, a[0], b[3]);

    let (r1, carry) = macx(r1, a[1], b[0]);
    let (r2, carry) = mac(r2, a[1], b[1], carry);
    let (r3, carry) = mac(r3, a[1], b[2], carry);
    let (r4, carry_out) = mac(carry_out, a[1], b[3], carry);

    let (r2, carry) = macx(r2, a[2], b[0]);
    let (r3, carry) = mac(r3, a[2], b[1], carry);
    let (r4, carry) = mac(r4, a[2], b[2], carry);
    let (r5, carry_out) = mac(carry_out, a[2], b[3], carry);

    let (r3, carry) = macx(r3, a[3], b[0]);
    let (r4, carry) = mac(r4, a[3], b[1], carry);
    let (r5, carry) = mac(r5, a[3], b[2], carry);
    let (r6, carry_out) = mac(carry_out, a[3], b[3], carry);

    [r0, r1, r2, r3, r4, r5, r6, carry_out]
}

#[cfg(test)]
mod test {
    use super::CurveEndo;
    use crate::bn256::G1;
    use ff::Field;
    use pasta_curves::Ep;
    use pasta_curves::Eq;
    use rand_core::OsRng;

    // naive glv multiplication implementation
    fn glv_mul<C: CurveEndo>(point: C, scalar: &C::ScalarExt) -> C {
        const WINDOW: usize = 3;
        // decompose scalar and convert to wnaf representation
        let (k1, k1_neg, k2, k2_neg) = C::decompose_scalar(scalar);

        let mut k1_wnaf: Vec<i64> = Vec::new();
        let mut k2_wnaf: Vec<i64> = Vec::new();
        wnaf::form(&mut k1_wnaf, k1.to_le_bytes(), WINDOW);
        wnaf::form(&mut k2_wnaf, k2.to_le_bytes(), WINDOW);

        let n = std::cmp::max(k1_wnaf.len(), k2_wnaf.len());
        k1_wnaf.resize(n, 0);
        k2_wnaf.resize(n, 0);

        // prepare tables
        let two_p = point.double();
        // T1 = {P, 3P, 5P, ...}
        let mut table_k1 = vec![point];
        // T2 = {λP, 3λP, 5λP, ...}
        let mut table_k2 = vec![point.endo()];
        for i in 1..WINDOW - 1 {
            table_k1.push(table_k1[i - 1] + two_p);
            table_k2.push(table_k1[i].endo())
        }
        if !k2_neg {
            table_k2.iter_mut().for_each(|p| *p = -*p);
        }
        if k1_neg {
            table_k1.iter_mut().for_each(|p| *p = -*p);
        }
        // TODO: batch affine tables for mixed add?

        macro_rules! add {
            ($acc:expr, $e:expr, $table:expr) => {
                let idx = ($e.abs() >> 1) as usize;
                $acc += if $e.is_positive() {
                    $table[idx]
                } else if $e.is_negative() {
                    -$table[idx]
                } else {
                    C::identity()
                };
            };
        }

        // apply simultaneus double add
        k1_wnaf
            .iter()
            .rev()
            .zip(k2_wnaf.iter().rev())
            .fold(C::identity(), |acc, (e1, e2)| {
                let mut acc = acc.double();
                add!(acc, e1, table_k1);
                add!(acc, e2, table_k2);
                acc
            })
    }

    fn run_glv_mul_test<C: CurveEndo>() {
        for _ in 0..10000 {
            let point = C::random(OsRng);
            let scalar = C::ScalarExt::random(OsRng);
            let r0 = point * scalar;
            let r1 = glv_mul(point, &scalar);
            assert_eq!(r0, r1);
        }
    }

    #[test]
    fn test_glv_mul() {
        run_glv_mul_test::<G1>();
        // run_glv_mul_test::<Eq>();
        // run_glv_mul_test::<Ep>(git);
    }

    #[test]
    fn test_wnaf_form() {
        use rand::Rng;
        fn from_wnaf(wnaf: &[i64]) -> u128 {
            wnaf.iter().rev().fold(0, |acc, next| {
                let mut acc = acc * 2;
                acc += *next as u128;
                acc
            })
        }
        for w in 2..64 {
            for e in 0..=u16::MAX {
                let mut wnaf = vec![];
                wnaf::form(&mut wnaf, e.to_le_bytes(), w);
                assert_eq!(e as u128, from_wnaf(&wnaf));
            }
        }
        let mut wnaf = vec![];
        for w in 2..64 {
            for e in u128::MAX - 10000..=u128::MAX {
                wnaf::form(&mut wnaf, e.to_le_bytes(), w);
                assert_eq!(e, from_wnaf(&wnaf));
            }
        }
        for w in 2..10 {
            for _ in 0..10000 {
                let e: u128 = OsRng.gen();
                wnaf::form(&mut wnaf, e.to_le_bytes(), w);
                assert_eq!(e as u128, from_wnaf(&wnaf));
            }
        }
    }

    // taken from zkcrypto/group
    mod wnaf {
        use std::convert::TryInto;

        #[derive(Debug, Clone)]
        struct LimbBuffer<'a> {
            buf: &'a [u8],
            cur_idx: usize,
            cur_limb: u64,
            next_limb: u64,
        }

        impl<'a> LimbBuffer<'a> {
            fn new(buf: &'a [u8]) -> Self {
                let mut ret = Self {
                    buf,
                    cur_idx: 0,
                    cur_limb: 0,
                    next_limb: 0,
                };

                // Initialise the limb buffers.
                ret.increment_limb();
                ret.increment_limb();
                ret.cur_idx = 0usize;
                ret
            }

            fn increment_limb(&mut self) {
                self.cur_idx += 1;
                self.cur_limb = self.next_limb;
                match self.buf.len() {
                    // There are no more bytes in the buffer; zero-extend.
                    0 => self.next_limb = 0,

                    // There are fewer bytes in the buffer than a u64 limb; zero-extend.
                    x @ 1..=7 => {
                        let mut next_limb = [0; 8];
                        next_limb[..x].copy_from_slice(self.buf);
                        self.next_limb = u64::from_le_bytes(next_limb);
                        self.buf = &[];
                    }

                    // There are at least eight bytes in the buffer; read the next u64 limb.
                    _ => {
                        let (next_limb, rest) = self.buf.split_at(8);
                        self.next_limb = u64::from_le_bytes(next_limb.try_into().unwrap());
                        self.buf = rest;
                    }
                }
            }

            fn get(&mut self, idx: usize) -> (u64, u64) {
                assert!([self.cur_idx, self.cur_idx + 1].contains(&idx));
                if idx > self.cur_idx {
                    self.increment_limb();
                }
                (self.cur_limb, self.next_limb)
            }
        }

        /// Replaces the contents of `wnaf` with the w-NAF representation of a little-endian
        /// scalar.
        pub(crate) fn form<S: AsRef<[u8]>>(wnaf: &mut Vec<i64>, c: S, window: usize) {
            // Required by the NAF definition
            debug_assert!(window >= 2);
            // Required so that the NAF digits fit in i64
            debug_assert!(window < 64);

            let bit_len = c.as_ref().len() * 8;

            wnaf.truncate(0);
            wnaf.reserve(bit_len + 1);

            // Initialise the current and next limb buffers.
            let mut limbs = LimbBuffer::new(c.as_ref());

            let width = 1u64 << window;
            let window_mask = width - 1;

            let mut pos = 0;
            let mut carry = 0;
            while pos <= bit_len {
                // Construct a buffer of bits of the scalar, starting at bit `pos`
                let u64_idx = pos / 64;
                let bit_idx = pos % 64;
                let (cur_u64, next_u64) = limbs.get(u64_idx);
                let bit_buf = if bit_idx + window < 64 {
                    // This window's bits are contained in a single u64
                    cur_u64 >> bit_idx
                } else {
                    // Combine the current u64's bits with the bits from the next u64
                    (cur_u64 >> bit_idx) | (next_u64 << (64 - bit_idx))
                };

                // Add the carry into the current window
                let window_val = carry + (bit_buf & window_mask);

                if window_val & 1 == 0 {
                    // If the window value is even, preserve the carry and emit 0.
                    // Why is the carry preserved?
                    // If carry == 0 and window_val & 1 == 0, then the next carry should be 0
                    // If carry == 1 and window_val & 1 == 0, then bit_buf & 1 == 1 so the next carry should be 1
                    wnaf.push(0);
                    pos += 1;
                } else {
                    wnaf.push(if window_val < width / 2 {
                        carry = 0;
                        window_val as i64
                    } else {
                        carry = 1;
                        (window_val as i64).wrapping_sub(width as i64)
                    });
                    wnaf.extend(std::iter::repeat(0).take(window - 1));
                    pos += window;
                }
            }
            wnaf.truncate(wnaf.len().saturating_sub(window - 1));
        }
    }
}
