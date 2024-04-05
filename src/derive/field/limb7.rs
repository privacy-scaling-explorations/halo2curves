#[macro_export]
macro_rules! field_arithmetic_7 {
    ($field:ident) => {
        impl From<$field> for [u64; NUM_LIMBS] {
            fn from(elt: $field) -> [u64; NUM_LIMBS] {
                $field::montgomery_reduce(&[
                    elt.0[0], elt.0[1], elt.0[2], elt.0[3], elt.0[4], elt.0[5], elt.0[6], 0, 0, 0,
                    0, 0, 0, 0,
                ])
                .0
            }
        }

        impl $field {
            const fn montgomery_form(val: [u64; NUM_LIMBS], r: $field) -> $field {
                let (r0, carry) = mac(0, val[0], r.0[0], 0);
                let (r1, carry) = mac(0, val[0], r.0[1], carry);
                let (r2, carry) = mac(0, val[0], r.0[2], carry);
                let (r3, carry) = mac(0, val[0], r.0[3], carry);
                let (r4, carry) = mac(0, val[0], r.0[4], carry);
                let (r5, carry) = mac(0, val[0], r.0[5], carry);
                let (r6, r7) = mac(0, val[0], r.0[6], carry);

                let (r1, carry) = mac(r1, val[1], r.0[0], 0);
                let (r2, carry) = mac(r2, val[1], r.0[1], carry);
                let (r3, carry) = mac(r3, val[1], r.0[2], carry);
                let (r4, carry) = mac(r4, val[1], r.0[3], carry);
                let (r5, carry) = mac(r5, val[1], r.0[4], carry);
                let (r6, carry) = mac(r6, val[1], r.0[5], carry);
                let (r7, r8) = mac(r7, val[1], r.0[6], carry);

                let (r2, carry) = mac(r2, val[2], r.0[0], 0);
                let (r3, carry) = mac(r3, val[2], r.0[1], carry);
                let (r4, carry) = mac(r4, val[2], r.0[2], carry);
                let (r5, carry) = mac(r5, val[2], r.0[3], carry);
                let (r6, carry) = mac(r6, val[2], r.0[4], carry);
                let (r7, carry) = mac(r7, val[2], r.0[5], carry);
                let (r8, r9) = mac(r8, val[2], r.0[6], carry);

                let (r3, carry) = mac(r3, val[3], r.0[0], 0);
                let (r4, carry) = mac(r4, val[3], r.0[1], carry);
                let (r5, carry) = mac(r5, val[3], r.0[2], carry);
                let (r6, carry) = mac(r6, val[3], r.0[3], carry);
                let (r7, carry) = mac(r7, val[3], r.0[4], carry);
                let (r8, carry) = mac(r8, val[3], r.0[5], carry);
                let (r9, r10) = mac(r9, val[3], r.0[6], carry);

                let (r4, carry) = mac(r4, val[4], r.0[0], 0);
                let (r5, carry) = mac(r5, val[4], r.0[1], carry);
                let (r6, carry) = mac(r6, val[4], r.0[2], carry);
                let (r7, carry) = mac(r7, val[4], r.0[3], carry);
                let (r8, carry) = mac(r8, val[4], r.0[4], carry);
                let (r9, carry) = mac(r9, val[4], r.0[5], carry);
                let (r10, r11) = mac(r10, val[4], r.0[6], carry);

                let (r5, carry) = mac(r5, val[5], r.0[0], 0);
                let (r6, carry) = mac(r6, val[5], r.0[1], carry);
                let (r7, carry) = mac(r7, val[5], r.0[2], carry);
                let (r8, carry) = mac(r8, val[5], r.0[3], carry);
                let (r9, carry) = mac(r9, val[5], r.0[4], carry);
                let (r10, carry) = mac(r10, val[5], r.0[5], carry);
                let (r11, r12) = mac(r11, val[5], r.0[6], carry);

                let (r6, carry) = mac(r6, val[6], r.0[0], 0);
                let (r7, carry) = mac(r7, val[6], r.0[1], carry);
                let (r8, carry) = mac(r8, val[6], r.0[2], carry);
                let (r9, carry) = mac(r9, val[6], r.0[3], carry);
                let (r10, carry) = mac(r10, val[6], r.0[4], carry);
                let (r11, carry) = mac(r11, val[6], r.0[5], carry);
                let (r12, r13) = mac(r12, val[6], r.0[6], carry);

                // Montgomery reduction
                let k = r0.wrapping_mul(INV);
                let (_, carry) = mac(r0, k, MODULUS.0[0], 0);
                let (r1, carry) = mac(r1, k, MODULUS.0[1], carry);
                let (r2, carry) = mac(r2, k, MODULUS.0[2], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[3], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[4], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[5], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[6], carry);
                let (r7, carry2) = adc(r7, 0, carry);

                let k = r1.wrapping_mul(INV);
                let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
                let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[4], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[5], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[6], carry);
                let (r8, carry2) = adc(r8, carry2, carry);

                let k = r2.wrapping_mul(INV);
                let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
                let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[4], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[5], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[6], carry);
                let (r9, carry2) = adc(r9, carry2, carry);

                let k = r3.wrapping_mul(INV);
                let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
                let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[4], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[5], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[6], carry);
                let (r10, carry2) = adc(r10, carry2, carry);

                let k = r4.wrapping_mul(INV);
                let (_, carry) = mac(r4, k, MODULUS.0[0], 0);
                let (r5, carry) = mac(r5, k, MODULUS.0[1], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[2], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[3], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[4], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[5], carry);
                let (r10, carry) = mac(r10, k, MODULUS.0[6], carry);
                let (r11, carry2) = adc(r11, carry2, carry);

                let k = r5.wrapping_mul(INV);
                let (_, carry) = mac(r5, k, MODULUS.0[0], 0);
                let (r6, carry) = mac(r6, k, MODULUS.0[1], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[2], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[3], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[4], carry);
                let (r10, carry) = mac(r10, k, MODULUS.0[5], carry);
                let (r11, carry) = mac(r11, k, MODULUS.0[6], carry);
                let (r12, carry2) = adc(r12, carry2, carry);

                let k = r6.wrapping_mul(INV);
                let (_, carry) = mac(r6, k, MODULUS.0[0], 0);
                let (r7, carry) = mac(r7, k, MODULUS.0[1], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[2], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[3], carry);
                let (r10, carry) = mac(r10, k, MODULUS.0[4], carry);
                let (r11, carry) = mac(r11, k, MODULUS.0[5], carry);
                let (r12, carry) = mac(r12, k, MODULUS.0[6], carry);
                let (r13, carry2) = adc(r13, carry2, carry);

                // Result may be within MODULUS of the correct value
                let (d0, borrow) = sbb(r7, MODULUS.0[0], 0);
                let (d1, borrow) = sbb(r8, MODULUS.0[1], borrow);
                let (d2, borrow) = sbb(r9, MODULUS.0[2], borrow);
                let (d3, borrow) = sbb(r10, MODULUS.0[3], borrow);
                let (d4, borrow) = sbb(r11, MODULUS.0[4], borrow);
                let (d5, borrow) = sbb(r12, MODULUS.0[5], borrow);
                let (d6, borrow) = sbb(r13, MODULUS.0[6], borrow);
                let (_, borrow) = sbb(carry2, 0, borrow);
                let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
                let (d3, carry) = adc(d3, MODULUS.0[3] & borrow, carry);
                let (d4, carry) = adc(d4, MODULUS.0[4] & borrow, carry);
                let (d5, carry) = adc(d5, MODULUS.0[5] & borrow, carry);
                let (d6, _) = adc(d6, MODULUS.0[6] & borrow, carry);

                $field([d0, d1, d2, d3, d4, d5, d6])
            }

            #[inline(always)]
            pub(crate) const fn montgomery_reduce(r: &[u64; 14]) -> $field {
                // The Montgomery reduction here is based on Algorithm 14.32 in
                // Handbook of Applied Cryptography
                // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

                let k = r[0].wrapping_mul(INV);
                let (_, carry) = mac(r[0], k, MODULUS.0[0], 0);
                let (r1, carry) = mac(r[1], k, MODULUS.0[1], carry);
                let (r2, carry) = mac(r[2], k, MODULUS.0[2], carry);
                let (r3, carry) = mac(r[3], k, MODULUS.0[3], carry);
                let (r4, carry) = mac(r[4], k, MODULUS.0[4], carry);
                let (r5, carry) = mac(r[5], k, MODULUS.0[5], carry);
                let (r6, carry) = mac(r[6], k, MODULUS.0[6], carry);
                let (r7, carry2) = adc(r[7], 0, carry);

                let k = r1.wrapping_mul(INV);
                let (_, carry) = mac(r1, k, MODULUS.0[0], 0);
                let (r2, carry) = mac(r2, k, MODULUS.0[1], carry);
                let (r3, carry) = mac(r3, k, MODULUS.0[2], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[3], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[4], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[5], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[6], carry);
                let (r8, carry2) = adc(r[8], carry2, carry);

                let k = r2.wrapping_mul(INV);
                let (_, carry) = mac(r2, k, MODULUS.0[0], 0);
                let (r3, carry) = mac(r3, k, MODULUS.0[1], carry);
                let (r4, carry) = mac(r4, k, MODULUS.0[2], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[3], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[4], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[5], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[6], carry);
                let (r9, carry2) = adc(r[9], carry2, carry);

                let k = r3.wrapping_mul(INV);
                let (_, carry) = mac(r3, k, MODULUS.0[0], 0);
                let (r4, carry) = mac(r4, k, MODULUS.0[1], carry);
                let (r5, carry) = mac(r5, k, MODULUS.0[2], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[3], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[4], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[5], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[6], carry);
                let (r10, carry2) = adc(r[10], carry2, carry);

                let k = r4.wrapping_mul(INV);
                let (_, carry) = mac(r4, k, MODULUS.0[0], 0);
                let (r5, carry) = mac(r5, k, MODULUS.0[1], carry);
                let (r6, carry) = mac(r6, k, MODULUS.0[2], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[3], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[4], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[5], carry);
                let (r10, carry) = mac(r10, k, MODULUS.0[6], carry);
                let (r11, carry2) = adc(r[11], carry2, carry);

                let k = r5.wrapping_mul(INV);
                let (_, carry) = mac(r5, k, MODULUS.0[0], 0);
                let (r6, carry) = mac(r6, k, MODULUS.0[1], carry);
                let (r7, carry) = mac(r7, k, MODULUS.0[2], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[3], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[4], carry);
                let (r10, carry) = mac(r10, k, MODULUS.0[5], carry);
                let (r11, carry) = mac(r11, k, MODULUS.0[6], carry);
                let (r12, carry2) = adc(r[12], carry2, carry);

                let k = r6.wrapping_mul(INV);
                let (_, carry) = mac(r6, k, MODULUS.0[0], 0);
                let (r7, carry) = mac(r7, k, MODULUS.0[1], carry);
                let (r8, carry) = mac(r8, k, MODULUS.0[2], carry);
                let (r9, carry) = mac(r9, k, MODULUS.0[3], carry);
                let (r10, carry) = mac(r10, k, MODULUS.0[4], carry);
                let (r11, carry) = mac(r11, k, MODULUS.0[5], carry);
                let (r12, carry) = mac(r12, k, MODULUS.0[6], carry);
                let (r13, _) = adc(r[13], carry2, carry);
                // Result may be within MODULUS of the correct value
                (&$field([r7, r8, r9, r10, r11, r12, r13])).sub(&MODULUS)
            }

            /// Adds `rhs` to `self`, returning the result.
            #[inline]
            pub const fn add(&self, rhs: &Self) -> Self {
                let (d0, carry) = adc(self.0[0], rhs.0[0], 0);
                let (d1, carry) = adc(self.0[1], rhs.0[1], carry);
                let (d2, carry) = adc(self.0[2], rhs.0[2], carry);
                let (d3, carry) = adc(self.0[3], rhs.0[3], carry);
                let (d4, carry) = adc(self.0[4], rhs.0[4], carry);
                let (d5, carry) = adc(self.0[5], rhs.0[5], carry);
                let (d6, _) = adc(self.0[6], rhs.0[6], carry);

                // Attempt to subtract the modulus, to ensure the value
                // is smaller than the modulus.
                (&$field([d0, d1, d2, d3, d4, d5, d6])).sub(&MODULUS)
            }

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
                let (r3, carry) = mac(0, self.0[0], self.0[3], carry);
                let (r4, carry) = mac(0, self.0[0], self.0[4], carry);
                let (r5, carry) = mac(0, self.0[0], self.0[5], carry);
                let (r6, r7) = mac(0, self.0[0], self.0[6], carry);

                let (r3, carry) = mac(r3, self.0[1], self.0[2], 0);
                let (r4, carry) = mac(r4, self.0[1], self.0[3], carry);
                let (r5, carry) = mac(r5, self.0[1], self.0[4], carry);
                let (r6, carry) = mac(r6, self.0[1], self.0[5], carry);
                let (r7, r8) = mac(r7, self.0[1], self.0[6], carry);

                let (r5, carry) = mac(r5, self.0[2], self.0[3], 0);
                let (r6, carry) = mac(r6, self.0[2], self.0[4], carry);
                let (r7, carry) = mac(r7, self.0[2], self.0[5], carry);
                let (r8, r9) = mac(r8, self.0[2], self.0[6], carry);

                let (r7, carry) = mac(r7, self.0[3], self.0[4], 0);
                let (r8, carry) = mac(r8, self.0[3], self.0[5], carry);
                let (r9, r10) = mac(r9, self.0[3], self.0[6], carry);

                let (r9, carry) = mac(r9, self.0[4], self.0[5], 0);
                let (r10, r11) = mac(r10, self.0[4], self.0[6], carry);

                let (r11, r12) = mac(r11, self.0[5], self.0[6], 0);

                let r13 = r12 >> 63;
                let r12 = (r12 << 1) | (r11 >> 63);
                let r11 = (r11 << 1) | (r10 >> 63);
                let r10 = (r10 << 1) | (r9 >> 63);
                let r9 = (r9 << 1) | (r8 >> 63);
                let r8 = (r8 << 1) | (r7 >> 63);
                let r7 = (r7 << 1) | (r6 >> 63);
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
                let (r7, carry) = adc(0, r7, carry);
                let (r8, carry) = mac(r8, self.0[4], self.0[4], carry);
                let (r9, carry) = adc(0, r9, carry);
                let (r10, carry) = mac(r10, self.0[5], self.0[5], carry);
                let (r11, carry) = adc(0, r11, carry);
                let (r12, carry) = mac(r12, self.0[6], self.0[6], carry);
                let (r13, _) = adc(0, r13, carry);

                $field::montgomery_reduce(&[
                    r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13,
                ])
            }

            /// Multiplies `rhs` by `self`, returning the result.
            #[inline]
            pub const fn mul(&self, rhs: &Self) -> $field {
                // Schoolbook multiplication
                let (r0, carry) = mac(0, self.0[0], rhs.0[0], 0);
                let (r1, carry) = mac(0, self.0[0], rhs.0[1], carry);
                let (r2, carry) = mac(0, self.0[0], rhs.0[2], carry);
                let (r3, carry) = mac(0, self.0[0], rhs.0[3], carry);
                let (r4, carry) = mac(0, self.0[0], rhs.0[4], carry);
                let (r5, carry) = mac(0, self.0[0], rhs.0[5], carry);
                let (r6, r7) = mac(0, self.0[0], rhs.0[6], carry);

                let (r1, carry) = mac(r1, self.0[1], rhs.0[0], 0);
                let (r2, carry) = mac(r2, self.0[1], rhs.0[1], carry);
                let (r3, carry) = mac(r3, self.0[1], rhs.0[2], carry);
                let (r4, carry) = mac(r4, self.0[1], rhs.0[3], carry);
                let (r5, carry) = mac(r5, self.0[1], rhs.0[4], carry);
                let (r6, carry) = mac(r6, self.0[1], rhs.0[5], carry);
                let (r7, r8) = mac(r7, self.0[1], rhs.0[6], carry);

                let (r2, carry) = mac(r2, self.0[2], rhs.0[0], 0);
                let (r3, carry) = mac(r3, self.0[2], rhs.0[1], carry);
                let (r4, carry) = mac(r4, self.0[2], rhs.0[2], carry);
                let (r5, carry) = mac(r5, self.0[2], rhs.0[3], carry);
                let (r6, carry) = mac(r6, self.0[2], rhs.0[4], carry);
                let (r7, carry) = mac(r7, self.0[2], rhs.0[5], carry);
                let (r8, r9) = mac(r8, self.0[2], rhs.0[6], carry);

                let (r3, carry) = mac(r3, self.0[3], rhs.0[0], 0);
                let (r4, carry) = mac(r4, self.0[3], rhs.0[1], carry);
                let (r5, carry) = mac(r5, self.0[3], rhs.0[2], carry);
                let (r6, carry) = mac(r6, self.0[3], rhs.0[3], carry);
                let (r7, carry) = mac(r7, self.0[3], rhs.0[4], carry);
                let (r8, carry) = mac(r8, self.0[3], rhs.0[5], carry);
                let (r9, r10) = mac(r9, self.0[3], rhs.0[6], carry);

                let (r4, carry) = mac(r4, self.0[4], rhs.0[0], 0);
                let (r5, carry) = mac(r5, self.0[4], rhs.0[1], carry);
                let (r6, carry) = mac(r6, self.0[4], rhs.0[2], carry);
                let (r7, carry) = mac(r7, self.0[4], rhs.0[3], carry);
                let (r8, carry) = mac(r8, self.0[4], rhs.0[4], carry);
                let (r9, carry) = mac(r9, self.0[4], rhs.0[5], carry);
                let (r10, r11) = mac(r10, self.0[4], rhs.0[6], carry);

                let (r5, carry) = mac(r5, self.0[5], rhs.0[0], 0);
                let (r6, carry) = mac(r6, self.0[5], rhs.0[1], carry);
                let (r7, carry) = mac(r7, self.0[5], rhs.0[2], carry);
                let (r8, carry) = mac(r8, self.0[5], rhs.0[3], carry);
                let (r9, carry) = mac(r9, self.0[5], rhs.0[4], carry);
                let (r10, carry) = mac(r10, self.0[5], rhs.0[5], carry);
                let (r11, r12) = mac(r11, self.0[5], rhs.0[6], carry);

                let (r6, carry) = mac(r6, self.0[6], rhs.0[0], 0);
                let (r7, carry) = mac(r7, self.0[6], rhs.0[1], carry);
                let (r8, carry) = mac(r8, self.0[6], rhs.0[2], carry);
                let (r9, carry) = mac(r9, self.0[6], rhs.0[3], carry);
                let (r10, carry) = mac(r10, self.0[6], rhs.0[4], carry);
                let (r11, carry) = mac(r11, self.0[6], rhs.0[5], carry);
                let (r12, r13) = mac(r12, self.0[6], rhs.0[6], carry);

                $field::montgomery_reduce(&[
                    r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13,
                ])
            }

            /// Subtracts `rhs` from `self`, returning the result.
            #[inline]
            pub const fn sub(&self, rhs: &Self) -> Self {
                let (d0, borrow) = sbb(self.0[0], rhs.0[0], 0);
                let (d1, borrow) = sbb(self.0[1], rhs.0[1], borrow);
                let (d2, borrow) = sbb(self.0[2], rhs.0[2], borrow);
                let (d3, borrow) = sbb(self.0[3], rhs.0[3], borrow);
                let (d4, borrow) = sbb(self.0[4], rhs.0[4], borrow);
                let (d5, borrow) = sbb(self.0[5], rhs.0[5], borrow);
                let (d6, borrow) = sbb(self.0[6], rhs.0[6], borrow);

                // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
                // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
                let (d0, carry) = adc(d0, MODULUS.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, MODULUS.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, MODULUS.0[2] & borrow, carry);
                let (d3, carry) = adc(d3, MODULUS.0[3] & borrow, carry);
                let (d4, carry) = adc(d4, MODULUS.0[4] & borrow, carry);
                let (d5, carry) = adc(d5, MODULUS.0[5] & borrow, carry);
                let (d6, _) = adc(d6, MODULUS.0[6] & borrow, carry);

                $field([d0, d1, d2, d3, d4, d5, d6])
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
                let (d3, borrow) = sbb(MODULUS.0[3], self.0[3], borrow);
                let (d4, borrow) = sbb(MODULUS.0[4], self.0[4], borrow);
                let (d5, borrow) = sbb(MODULUS.0[5], self.0[5], borrow);
                let (d6, _) = sbb(MODULUS.0[6], self.0[6], borrow);

                // `tmp` could be `MODULUS` if `self` was zero. Create a mask that is
                // zero if `self` was zero, and `u64::max_value()` if self was nonzero.
                let mask = (((self.0[0]
                    | self.0[1]
                    | self.0[2]
                    | self.0[3]
                    | self.0[4]
                    | self.0[5]
                    | self.0[6])
                    == 0) as u64)
                    .wrapping_sub(1);

                $field([
                    d0 & mask,
                    d1 & mask,
                    d2 & mask,
                    d3 & mask,
                    d4 & mask,
                    d5 & mask,
                    d6 & mask,
                ])
            }
        }
    };
}
