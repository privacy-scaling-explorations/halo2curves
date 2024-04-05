#[macro_export]
macro_rules! const_montgomery_4 {
    ($field:ident) => {
        impl $field {
            const fn montgomery_form(val: [u64; 4], r: $field) -> $field {
                // Converts a 4 64-bit limb value into its congruent field representation.
                // If `val` represents a 256 bit value then `r` should be R^2,
                // if `val` represents the 256 MSB of a 512 bit value, then `r` should be R^3.

                let (r0, carry) = mac(0, val[0], r.0[0], 0);
                let (r1, carry) = mac(0, val[0], r.0[1], carry);
                let (r2, carry) = mac(0, val[0], r.0[2], carry);
                let (r3, r4) = mac(0, val[0], r.0[3], carry);

                let (r1, carry) = mac(r1, val[1], r.0[0], 0);
                let (r2, carry) = mac(r2, val[1], r.0[1], carry);
                let (r3, carry) = mac(r3, val[1], r.0[2], carry);
                let (r4, r5) = mac(r4, val[1], r.0[3], carry);

                let (r2, carry) = mac(r2, val[2], r.0[0], 0);
                let (r3, carry) = mac(r3, val[2], r.0[1], carry);
                let (r4, carry) = mac(r4, val[2], r.0[2], carry);
                let (r5, r6) = mac(r5, val[2], r.0[3], carry);

                let (r3, carry) = mac(r3, val[3], r.0[0], 0);
                let (r4, carry) = mac(r4, val[3], r.0[1], carry);
                let (r5, carry) = mac(r5, val[3], r.0[2], carry);
                let (r6, r7) = mac(r6, val[3], r.0[3], carry);

                // Montgomery reduction
                let k = r0.wrapping_mul(INV);
                let (_, carry) = mac(r0, k, MODULUS.0[0], 0);
                let (r1, carry) = mac(r1, k, MODULUS.0[1], carry);
                let (r2, carry) = mac(r2, k, MODULUS.0[2], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[3], carry);
                let (r4, carry2) = adc(r4, 0, carry);

                let k = r1.wrapping_mul(INV);
                let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
                let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
                let (r5, carry2) = adc(r5, carry2, carry);

                let k = r2.wrapping_mul(INV);
                let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
                let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
                let (r6, carry2) = adc(r6, carry2, carry);

                let k = r3.wrapping_mul(INV);
                let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
                let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
                let (r7, carry2) = adc(r7, carry2, carry);

                // Result may be within MODULUS of the correct value
                let (d0, borrow) = sbb(r4, MODULUS.0[0], 0);
                let (d1, borrow) = sbb(r5, MODULUS.0[1], borrow);
                let (d2, borrow) = sbb(r6, MODULUS.0[2], borrow);
                let (d3, borrow) = sbb(r7, MODULUS.0[3], borrow);
                let (_, borrow) = sbb(carry2, 0, borrow);
                let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
                let (d3, _) = adc(d3, MODULUS.0[3] & borrow, carry);

                $field([d0, d1, d2, d3])
            }
        }
    };
}

#[macro_export]
macro_rules! field_arithmetic_4 {
    ($field:ident, $field_type:ident) => {
        field_specific_4!($field, $field_type);

        impl $field {
            /// Doubles this field element.
            #[inline]
            pub const fn double(&self) -> $field {
                self.add(self)
            }

            /// Squares this element.
            #[inline]
            pub const fn square(&self) -> $field {
                let (r1, carry) = mac(0, self.0[0], self.0[1], 0);
                let (r2, carry) = mac(0, self.0[0], self.0[2], carry);
                let (r3, r4) = mac(0, self.0[0], self.0[3], carry);

                let (r3, carry) = mac(r3, self.0[1], self.0[2], 0);
                let (r4, r5) = mac(r4, self.0[1], self.0[3], carry);

                let (r5, r6) = mac(r5, self.0[2], self.0[3], 0);

                let r7 = r6 >> 63;
                let r6 = (r6 << 1) | (r5 >> 63);
                let r5 = (r5 << 1) | (r4 >> 63);
                let r4 = (r4 << 1) | (r3 >> 63);
                let r3 = (r3 << 1) | (r2 >> 63);
                let r2 = (r2 << 1) | (r1 >> 63);
                let r1 = r1 << 1;

                let (r0, carry) = mac(0, self.0[0], self.0[0], 0);
                let (r1, carry) = adc(0, r1, carry);
                let (r2, carry) = mac(r2, self.0[1], self.0[1], carry);
                let (r3, carry) = adc(0, r3, carry);
                let (r4, carry) = mac(r4, self.0[2], self.0[2], carry);
                let (r5, carry) = adc(0, r5, carry);
                let (r6, carry) = mac(r6, self.0[3], self.0[3], carry);
                let (r7, _) = adc(0, r7, carry);

                $field::montgomery_reduce(&[r0, r1, r2, r3, r4, r5, r6, r7])
            }

            /// Multiplies `rhs` by `self`, returning the result.
            #[inline(always)]
            #[unroll::unroll_for_loops]
            #[allow(unused_assignments)]
            pub const fn mul(&self, rhs: &Self) -> Self {
                // Fast Coarsely Integrated Operand Scanning (CIOS) as described
                // in Algorithm 2 of EdMSM: https://eprint.iacr.org/2022/1400.pdf
                //
                // Cannot use the fast version (algorithm 2) if
                // modulus_high_word >= (WORD_SIZE - 1) / 2 - 1 = (2^64 - 1)/2 - 1

                if MODULUS.0[3] < (u64::MAX / 2) {
                    const N: usize = 4;
                    let mut t: [u64; N] = [0u64; N];
                    let mut c_2: u64;
                    for i in 0..4 {
                        let mut c: u64 = 0u64;
                        for j in 0..4 {
                            (t[j], c) = mac(t[j], self.0[j], rhs.0[i], c);
                        }
                        c_2 = c;

                        let m = t[0].wrapping_mul(INV);
                        (_, c) = macx(t[0], m, MODULUS.0[0]);

                        for j in 1..4 {
                            (t[j - 1], c) = mac(t[j], m, MODULUS.0[j], c);
                        }
                        (t[N - 1], _) = adc(c_2, c, 0);
                    }

                    if bigint_geq(&t, &MODULUS.0) {
                        let mut borrow = 0;
                        (t[0], borrow) = sbb(t[0], MODULUS.0[0], borrow);
                        (t[1], borrow) = sbb(t[1], MODULUS.0[1], borrow);
                        (t[2], borrow) = sbb(t[2], MODULUS.0[2], borrow);
                        (t[3], borrow) = sbb(t[3], MODULUS.0[3], borrow);
                    }
                    $field(t)
                } else {
                    // Schoolbook multiplication

                    let (r0, carry) = mac(0, self.0[0], rhs.0[0], 0);
                    let (r1, carry) = mac(0, self.0[0], rhs.0[1], carry);
                    let (r2, carry) = mac(0, self.0[0], rhs.0[2], carry);
                    let (r3, r4) = mac(0, self.0[0], rhs.0[3], carry);

                    let (r1, carry) = mac(r1, self.0[1], rhs.0[0], 0);
                    let (r2, carry) = mac(r2, self.0[1], rhs.0[1], carry);
                    let (r3, carry) = mac(r3, self.0[1], rhs.0[2], carry);
                    let (r4, r5) = mac(r4, self.0[1], rhs.0[3], carry);

                    let (r2, carry) = mac(r2, self.0[2], rhs.0[0], 0);
                    let (r3, carry) = mac(r3, self.0[2], rhs.0[1], carry);
                    let (r4, carry) = mac(r4, self.0[2], rhs.0[2], carry);
                    let (r5, r6) = mac(r5, self.0[2], rhs.0[3], carry);

                    let (r3, carry) = mac(r3, self.0[3], rhs.0[0], 0);
                    let (r4, carry) = mac(r4, self.0[3], rhs.0[1], carry);
                    let (r5, carry) = mac(r5, self.0[3], rhs.0[2], carry);
                    let (r6, r7) = mac(r6, self.0[3], rhs.0[3], carry);

                    $field::montgomery_reduce(&[r0, r1, r2, r3, r4, r5, r6, r7])
                }
            }

            /// Subtracts `rhs` from `self`, returning the result.
            #[inline]
            pub const fn sub(&self, rhs: &Self) -> Self {
                let (d0, borrow) = sbb(self.0[0], rhs.0[0], 0);
                let (d1, borrow) = sbb(self.0[1], rhs.0[1], borrow);
                let (d2, borrow) = sbb(self.0[2], rhs.0[2], borrow);
                let (d3, borrow) = sbb(self.0[3], rhs.0[3], borrow);

                // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
                // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
                let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
                let (d3, _) = adc(d3, MODULUS.0[3] & borrow, carry);

                $field([d0, d1, d2, d3])
            }

            /// Negates `self`.
            #[inline]
            pub const fn neg(&self) -> Self {
                // Subtract `self` from `MODULUS` to negate. Ignore the final
                // borrow because it cannot underflow; self is guaranteed to
                // be in the field.
                let (d0, borrow) = sbb(MODULUS.0[0], self.0[0], 0);
                let (d1, borrow) = sbb(MODULUS.0[1], self.0[1], borrow);
                let (d2, borrow) = sbb(MODULUS.0[2], self.0[2], borrow);
                let (d3, _) = sbb(MODULUS.0[3], self.0[3], borrow);

                // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
                // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
                let mask =
                    (((self.0[0] | self.0[1] | self.0[2] | self.0[3]) == 0) as u64).wrapping_sub(1);

                $field([d0 & mask, d1 & mask, d2 & mask, d3 & mask])
            }

            /// Montgomery reduce where last 4 registers are 0
            #[inline(always)]
            pub(crate) const fn montgomery_reduce_short(r: &[u64; 4]) -> $field {
                // The Montgomery reduction here is based on Algorithm 14.32 in
                // Handbook of Applied Cryptography
                // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

                let k = r[0].wrapping_mul(INV);
                let (_, r0) = macx(r[0], k, MODULUS.0[0]);
                let (r1, r0) = mac(r[1], k, MODULUS.0[1], r0);
                let (r2, r0) = mac(r[2], k, MODULUS.0[2], r0);
                let (r3, r0) = mac(r[3], k, MODULUS.0[3], r0);

                let k = r1.wrapping_mul(INV);
                let (_, r1) = macx(r1, k, MODULUS.0[0]);
                let (r2, r1) = mac(r2, k, MODULUS.0[1], r1);
                let (r3, r1) = mac(r3, k, MODULUS.0[2], r1);
                let (r0, r1) = mac(r0, k, MODULUS.0[3], r1);

                let k = r2.wrapping_mul(INV);
                let (_, r2) = macx(r2, k, MODULUS.0[0]);
                let (r3, r2) = mac(r3, k, MODULUS.0[1], r2);
                let (r0, r2) = mac(r0, k, MODULUS.0[2], r2);
                let (r1, r2) = mac(r1, k, MODULUS.0[3], r2);

                let k = r3.wrapping_mul(INV);
                let (_, r3) = macx(r3, k, MODULUS.0[0]);
                let (r0, r3) = mac(r0, k, MODULUS.0[1], r3);
                let (r1, r3) = mac(r1, k, MODULUS.0[2], r3);
                let (r2, r3) = mac(r2, k, MODULUS.0[3], r3);

                // Result may be within MODULUS of the correct value
                (&$field([r0, r1, r2, r3])).sub(&MODULUS)
            }
        }

        impl From<$field> for [u64; 4] {
            fn from(elt: $field) -> [u64; 4] {
                // Turn into canonical form by computing
                // (a.R) / R = a
                $field::montgomery_reduce_short(&elt.0).0
            }
        }
    };
}

#[macro_export]
macro_rules! field_specific_4 {
    ($field:ident, sparse) => {
        impl $field {
            /// Adds `rhs` to `self`, returning the result.
            #[inline]
            pub const fn add(&self, rhs: &Self) -> Self {
                let (d0, carry) = adc(self.0[0], rhs.0[0], 0);
                let (d1, carry) = adc(self.0[1], rhs.0[1], carry);
                let (d2, carry) = adc(self.0[2], rhs.0[2], carry);
                let (d3, _) = adc(self.0[3], rhs.0[3], carry);

                // Attempt to subtract the modulus, to ensure the value
                // is smaller than the modulus.
                (&$field([d0, d1, d2, d3])).sub(&MODULUS)
            }

            #[inline(always)]
            pub(crate) const fn montgomery_reduce(r: &[u64; 8]) -> $field {
                // The Montgomery reduction here is based on Algorithm 14.32 in
                // Handbook of Applied Cryptography
                // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

                let k = r[0].wrapping_mul(INV);
                let (_, carry) = mac(r[0], k, MODULUS.0[0], 0);
                let (r1, carry) = mac(r[1], k, MODULUS.0[1], carry);
                let (r2, carry) = mac(r[2], k, MODULUS.0[2], carry);
                let (r3, carry) = mac(r[3], k, MODULUS.0[3], carry);
                let (r4, carry2) = adc(r[4], 0, carry);

                let k = r1.wrapping_mul(INV);
                let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
                let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
                let (r5, carry2) = adc(r[5], carry2, carry);

                let k = r2.wrapping_mul(INV);
                let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
                let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
                let (r6, carry2) = adc(r[6], carry2, carry);

                let k = r3.wrapping_mul(INV);
                let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
                let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
                let (r7, _) = adc(r[7], carry2, carry);

                // Result may be within MODULUS of the correct value
                (&$field([r4, r5, r6, r7])).sub(&MODULUS)
            }
        }
    };
    ($field:ident, dense) => {
        impl $field {
            /// Adds `rhs` to `self`, returning the result.
            #[inline]
            pub const fn add(&self, rhs: &Self) -> Self {
                let (d0, carry) = adc(self.0[0], rhs.0[0], 0);
                let (d1, carry) = adc(self.0[1], rhs.0[1], carry);
                let (d2, carry) = adc(self.0[2], rhs.0[2], carry);
                let (d3, carry) = adc(self.0[3], rhs.0[3], carry);

                // Attempt to subtract the modulus, to ensure the value
                // is smaller than the modulus.
                let (d0, borrow) = sbb(d0, MODULUS.0[0], 0);
                let (d1, borrow) = sbb(d1, MODULUS.0[1], borrow);
                let (d2, borrow) = sbb(d2, MODULUS.0[2], borrow);
                let (d3, borrow) = sbb(d3, MODULUS.0[3], borrow);
                let (_, borrow) = sbb(carry, 0, borrow);

                let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
                let (d3, _) = adc(d3, MODULUS.0[3] & borrow, carry);

                $field([d0, d1, d2, d3])
            }

            #[inline(always)]
            pub(crate) const fn montgomery_reduce(r: &[u64; 8]) -> Self {
                // The Montgomery reduction here is based on Algorithm 14.32 in
                // Handbook of Applied Cryptography
                // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

                let k = r[0].wrapping_mul(INV);
                let (_, carry) = mac(r[0], k, MODULUS.0[0], 0);
                let (r1, carry) = mac(r[1], k, MODULUS.0[1], carry);
                let (r2, carry) = mac(r[2], k, MODULUS.0[2], carry);
                let (r3, carry) = mac(r[3], k, MODULUS.0[3], carry);
                let (r4, carry2) = adc(r[4], 0, carry);

                let k = r1.wrapping_mul(INV);
                let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
                let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
                let (r5, carry2) = adc(r[5], carry2, carry);

                let k = r2.wrapping_mul(INV);
                let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
                let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
                let (r6, carry2) = adc(r[6], carry2, carry);

                let k = r3.wrapping_mul(INV);
                let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
                let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
                let (r7, carry2) = adc(r[7], carry2, carry);

                // Result may be within MODULUS of the correct value
                let (d0, borrow) = sbb(r4, MODULUS.0[0], 0);
                let (d1, borrow) = sbb(r5, MODULUS.0[1], borrow);
                let (d2, borrow) = sbb(r6, MODULUS.0[2], borrow);
                let (d3, borrow) = sbb(r7, MODULUS.0[3], borrow);
                let (_, borrow) = sbb(carry2, 0, borrow);

                let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
                let (d3, _) = adc(d3, MODULUS.0[3] & borrow, carry);

                $field([d0, d1, d2, d3])
            }
        }
    };
}
