macro_rules! assembly_field {
    ($field:ident, $modulus:ident, $inv:ident) => {
        impl $field {
            /// Returns zero, the additive identity.
            #[inline]
            pub const fn zero() -> $field {
                $field([0, 0, 0, 0])
            }

            /// Returns one, the multiplicative identity.
            #[inline]
            pub const fn one() -> $field {
                R
            }

            /// Doubles this field element.
            #[inline]
            pub fn double(&self) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;
                unsafe {
                    asm!(
                        // load a array to former registers
                        "mov r8, qword ptr [{a_ptr} + 0]",
                        "mov r9, qword ptr [{a_ptr} + 8]",
                        "mov r10, qword ptr [{a_ptr} + 16]",
                        "mov r11, qword ptr [{a_ptr} + 24]",

                        // // add a array and b array with carry
                        "add r8, r8",
                        "adcx r9, r9",
                        "adcx r10, r10",
                        "adcx r11, r11",

                        // copy result array to latter registers
                        "mov r12, r8",
                        "mov r13, r9",
                        "mov r14, r10",
                        "mov r15, r11",

                        // mod reduction
                        "sub r12, qword ptr [{m_ptr} + 0]",
                        "sbb r13, qword ptr [{m_ptr} + 8]",
                        "sbb r14, qword ptr [{m_ptr} + 16]",
                        "sbb r15, qword ptr [{m_ptr} + 24]",

                        // if carry copy former registers to out areas
                        "cmovc r12, r8",
                        "cmovc r13, r9",
                        "cmovc r14, r10",
                        "cmovc r15, r11",

                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        a_ptr = in(reg) self.0.as_ptr(),
                        out("r8") _,
                        out("r9") _,
                        out("r10") _,
                        out("r11") _,
                        out("r12") r0,
                        out("r13") r1,
                        out("r14") r2,
                        out("r15") r3,
                        options(pure, readonly, nostack)
                    );
                }
                $field([r0, r1, r2, r3])
            }

            fn from_u512(limbs: [u64; 8]) -> $field {
                // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
                // with the higher bits multiplied by 2^256. Thus, we perform two reductions
                //
                // 1. the lower bits are multiplied by R^2, as normal
                // 2. the upper bits are multiplied by R^2 * 2^256 = R^3
                //
                // and computing their sum in the field. It remains to see that arbitrary 256-bit
                // numbers can be placed into Montgomery form safely using the reduction. The
                // reduction works so long as the product is less than R=2^256 multiplied by
                // the modulus. This holds because for any `c` smaller than the modulus, we have
                // that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
                // reduction always works so long as `c` is in the field; in this case it is either the
                // constant `R2` or `R3`.
                let d0 = $field([limbs[0], limbs[1], limbs[2], limbs[3]]);
                let d1 = $field([limbs[4], limbs[5], limbs[6], limbs[7]]);
                // Convert to Montgomery form
                d0 * R2 + d1 * R3
            }

            /// Converts from an integer represented in little endian
            /// into its (congruent) `Fr` representation.
            pub const fn from_raw(val: [u64; 4]) -> $field {
                let (r0, carry) = mac(0, val[0], R2.0[0], 0);
                let (r1, carry) = mac(0, val[0], R2.0[1], carry);
                let (r2, carry) = mac(0, val[0], R2.0[2], carry);
                let (r3, r4) = mac(0, val[0], R2.0[3], carry);

                let (r1, carry) = mac(r1, val[1], R2.0[0], 0);
                let (r2, carry) = mac(r2, val[1], R2.0[1], carry);
                let (r3, carry) = mac(r3, val[1], R2.0[2], carry);
                let (r4, r5) = mac(r4, val[1], R2.0[3], carry);

                let (r2, carry) = mac(r2, val[2], R2.0[0], 0);
                let (r3, carry) = mac(r3, val[2], R2.0[1], carry);
                let (r4, carry) = mac(r4, val[2], R2.0[2], carry);
                let (r5, r6) = mac(r5, val[2], R2.0[3], carry);

                let (r3, carry) = mac(r3, val[3], R2.0[0], 0);
                let (r4, carry) = mac(r4, val[3], R2.0[1], carry);
                let (r5, carry) = mac(r5, val[3], R2.0[2], carry);
                let (r6, r7) = mac(r6, val[3], R2.0[3], carry);

                let k = r0.wrapping_mul(INV);
                let (_, carry) = mac(r0, k, $modulus.0[0], 0);
                let (r1, carry) = mac(r1, k, $modulus.0[1], carry);
                let (r2, carry) = mac(r2, k, $modulus.0[2], carry);
                let (r3, carry) = mac(r3, k, $modulus.0[3], carry);
                let (r4, carry2) = adc(r4, 0, carry);

                let k = r1.wrapping_mul(INV);
                let (_, carry) = mac(r1, k, $modulus.0[0], 0);
                let (r2, carry) = mac(r2, k, $modulus.0[1], carry);
                let (r3, carry) = mac(r3, k, $modulus.0[2], carry);
                let (r4, carry) = mac(r4, k, $modulus.0[3], carry);
                let (r5, carry2) = adc(r5, carry2, carry);

                let k = r2.wrapping_mul(INV);
                let (_, carry) = mac(r2, k, $modulus.0[0], 0);
                let (r3, carry) = mac(r3, k, $modulus.0[1], carry);
                let (r4, carry) = mac(r4, k, $modulus.0[2], carry);
                let (r5, carry) = mac(r5, k, $modulus.0[3], carry);
                let (r6, carry2) = adc(r6, carry2, carry);

                let k = r3.wrapping_mul(INV);
                let (_, carry) = mac(r3, k, $modulus.0[0], 0);
                let (r4, carry) = mac(r4, k, $modulus.0[1], carry);
                let (r5, carry) = mac(r5, k, $modulus.0[2], carry);
                let (r6, carry) = mac(r6, k, $modulus.0[3], carry);
                let (r7, _) = adc(r7, carry2, carry);

                let (d0, borrow) = sbb(r4, $modulus.0[0], 0);
                let (d1, borrow) = sbb(r5, $modulus.0[1], borrow);
                let (d2, borrow) = sbb(r6, $modulus.0[2], borrow);
                let (d3, borrow) = sbb(r7, $modulus.0[3], borrow);

                // If underflow occurred on the final limb, borrow = 0xfff...fff, otherwise
                // borrow = 0x000...000. Thus, we use it as a mask to conditionally add the modulus.
                let (d0, carry) = adc(d0, $modulus.0[0] & borrow, 0);
                let (d1, carry) = adc(d1, $modulus.0[1] & borrow, carry);
                let (d2, carry) = adc(d2, $modulus.0[2] & borrow, carry);
                let (d3, _) = adc(d3, $modulus.0[3] & borrow, carry);

                $field([d0, d1, d2, d3])
            }

            /// Squares this element.
            #[inline]
            pub fn square(&self) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;
                unsafe {
                    asm!(
                        // schoolbook multiplication
                        //    *    |   a0    |   a1    |   a2    |   a3
                        //    b0   | b0 * a0 | b0 * a1 | b0 * a2 | b0 * a3
                        //    b1   | b1 * a0 | b1 * a1 | b1 * a2 | b1 * a3
                        //    b2   | b2 * a0 | b2 * a1 | b2 * a2 | b2 * a3
                        //    b3   | b3 * a0 | b3 * a1 | b3 * a2 | b3 * a3

                        // load value to registers
                        "mov r13, qword ptr [{a_ptr} + 0]",
                        "mov r14, qword ptr [{a_ptr} + 8]",
                        "mov r15, qword ptr [{a_ptr} + 16]",

                        // `a0`
                        "mov rdx, r13",

                        // a0 * b0
                        "mulx r9, r8, r13",

                        // a0 * b1
                        "mulx r10, rax, r14",
                        "add r9, rax",

                        // a0 * b2
                        "mulx r11, rax, r15",
                        "adcx r10, rax",

                        // a0 * b3
                        "mulx r12, rax, qword ptr [{a_ptr} + 24]",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // `a1`
                        "mov rdx, r14",

                        // a1 * b0
                        "mulx rcx, rax, r13",
                        "add r9, rax",
                        "adcx r10, rcx",
                        "adc r11, 0",

                        // a1 * b1
                        "mulx rcx, rax, r14",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",
                        "xor r13, r13",

                        // a1 * b2
                        "mulx rcx, rax, r15",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",
                        "xor r14, r14",

                        // a1 * b3
                        "mulx rcx, rax, qword ptr [{a_ptr} + 24]",
                        "add r12, rax",
                        "adcx r13, rcx",
                        "adc r14, 0",

                        // `a2`
                        "mov rdx, r15",

                        // a2 * b0
                        "mulx rcx, rax, qword ptr [{a_ptr} + 0]",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",

                        // a2 * b1
                        "mulx rcx, rax, qword ptr [{a_ptr} + 8]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // a2 * b2
                        "mulx rcx, rax, r15",
                        "add r12, rax",
                        "adcx r13, rcx",
                        "adc r14, 0",
                        "xor r15, r15",

                        // a2 * b3
                        "mulx rcx, rax, qword ptr [{a_ptr} + 24]",
                        "add r13, rax",
                        "adcx r14, rcx",
                        "adc r15, 0",

                        // `a3`
                        "mov rdx, qword ptr [{a_ptr} + 24]",

                        // a3 * b0
                        "mulx rcx, rax, qword ptr [{a_ptr} + 0]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // a3 * b1
                        "mulx rcx, rax, qword ptr [{a_ptr} + 8]",
                        "add r12, rax",
                        "adcx r13, rcx",
                        "adc r14, 0",

                        // a3 * b2
                        "mulx rcx, rax, qword ptr [{a_ptr} + 16]",
                        "add r13, rax",
                        "adcx r14, rcx",
                        "adc r15, 0",

                        // a3 * b3
                        "mulx rcx, rax, qword ptr [{a_ptr} + 24]",
                        "add r14, rax",
                        "adc r15, rcx",

                        // montgomery reduction
                        // r8 ~ r15

                        // `r8` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r8",

                        // r8' * m0
                        "mulx rcx, rax, qword ptr [{m_ptr} + 0]",
                        "add r8, rax",
                        "adcx r9, rcx",
                        "adc r10, 0",

                        // r8' * m1
                        "mulx rcx, rax, qword ptr [{m_ptr} + 8]",
                        "add r9, rax",
                        "adcx r10, rcx",
                        "adc r11, 0",

                        // r8' * m2
                        "mulx rcx, rax, qword ptr [{m_ptr} + 16]",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",

                        // r8' * m3
                        "mulx rcx, rax, qword ptr [{m_ptr} + 24]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // `r9` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r9",

                        // r9' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r9, rcx",
                        "adcx r10, rax",
                        "adc r11, 0",

                        // r9' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r10, rcx",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // r9' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r9' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // `r10` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r10",

                        // r10' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r10, rcx",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // r10' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r10' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // r10' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r13, rcx",
                        "adcx r14, rax",
                        "adc r15, 0",

                        // `r11` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r11",

                        // r11' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r11' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // r11' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r13, rcx",
                        "adcx r14, rax",
                        "adc r15, 0",

                        // r11' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r14, rcx",
                        "adcx r15, rax",

                        // reduction if limbs is greater then mod
                        "mov r8, r12",
                        "mov r9, r13",
                        "mov r10, r14",
                        "mov r11, r15",

                        "sub r8, qword ptr [{m_ptr} + 0]",
                        "sbb r9, qword ptr [{m_ptr} + 8]",
                        "sbb r10, qword ptr [{m_ptr} + 16]",
                        "sbb r11, qword ptr [{m_ptr} + 24]",

                        "cmovc r8, r12",
                        "cmovc r9, r13",
                        "cmovc r10, r14",
                        "cmovc r11, r15",

                        "mov r12, r8",
                        "mov r13, r9",
                        "mov r14, r10",
                        "mov r15, r11",

                        "sub r12, qword ptr [{m_ptr} + 0]",
                        "sbb r13, qword ptr [{m_ptr} + 8]",
                        "sbb r14, qword ptr [{m_ptr} + 16]",
                        "sbb r15, qword ptr [{m_ptr} + 24]",

                        "cmovc r12, r8",
                        "cmovc r13, r9",
                        "cmovc r14, r10",
                        "cmovc r15, r11",

                        a_ptr = in(reg) self.0.as_ptr(),
                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        inv = const $inv,
                        out("rax") _,
                        out("rcx") _,
                        out("rdx") _,
                        out("r8") _,
                        out("r9") _,
                        out("r10") _,
                        out("r11") _,
                        out("r12") r0,
                        out("r13") r1,
                        out("r14") r2,
                        out("r15") r3,
                        options(pure, readonly, nostack)
                    )
                }

                $field([r0, r1, r2, r3])
            }

            #[inline(always)]
            fn montgomery_reduce(a: &[u64; 8]) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;

                unsafe {
                    asm!(
                        // The Montgomery reduction here is based on Algorithm 14.32 in
                        // Handbook of Applied Cryptography
                        // <https://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

                        "mov r8, qword ptr [{a_ptr} + 0]",
                        "mov r9, qword ptr [{a_ptr} + 8]",
                        "mov r10, qword ptr [{a_ptr} + 16]",
                        "mov r11, qword ptr [{a_ptr} + 24]",
                        "mov r12, qword ptr [{a_ptr} + 32]",
                        "mov r13, qword ptr [{a_ptr} + 40]",
                        "mov r14, qword ptr [{a_ptr} + 48]",
                        "mov r15, qword ptr [{a_ptr} + 56]",

                        // `r8` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r8",

                        // r8' * m0
                        "mulx rcx, rax, qword ptr [{m_ptr} + 0]",
                        "add r8, rax",
                        "adcx r9, rcx",
                        "adc r10, 0",

                        // r8' * m1
                        "mulx rcx, rax, qword ptr [{m_ptr} + 8]",
                        "add r9, rax",
                        "adcx r10, rcx",
                        "adc r11, 0",

                        // // r8' * m2
                        "mulx rcx, rax, qword ptr [{m_ptr} + 16]",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",

                        // // r8' * m3
                        "mulx rcx, rax, qword ptr [{m_ptr} + 24]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // `r9` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r9",

                        // r9' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r9, rcx",
                        "adcx r10, rax",
                        "adc r11, 0",

                        // r9' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r10, rcx",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // r9' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r9' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // `r10` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r10",

                        // r10' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r10, rcx",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // r10' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r10' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // r10' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r13, rcx",
                        "adcx r14, rax",
                        "adc r15, 0",

                        // `r11` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r11",

                        // r11' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r11' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // r11' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r13, rcx",
                        "adcx r14, rax",
                        "adc r15, 0",

                        // r11' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r14, rcx",
                        "adcx r15, rax",

                        // reduction if limbs is greater then mod
                        "mov r8, r12",
                        "mov r9, r13",
                        "mov r10, r14",
                        "mov r11, r15",

                        "sub r8, qword ptr [{m_ptr} + 0]",
                        "sbb r9, qword ptr [{m_ptr} + 8]",
                        "sbb r10, qword ptr [{m_ptr} + 16]",
                        "sbb r11, qword ptr [{m_ptr} + 24]",

                        "cmovc r8, r12",
                        "cmovc r9, r13",
                        "cmovc r10, r14",
                        "cmovc r11, r15",

                        "mov r12, r8",
                        "mov r13, r9",
                        "mov r14, r10",
                        "mov r15, r11",

                        "sub r12, qword ptr [{m_ptr} + 0]",
                        "sbb r13, qword ptr [{m_ptr} + 8]",
                        "sbb r14, qword ptr [{m_ptr} + 16]",
                        "sbb r15, qword ptr [{m_ptr} + 24]",

                        "cmovc r12, r8",
                        "cmovc r13, r9",
                        "cmovc r14, r10",
                        "cmovc r15, r11",

                        a_ptr = in(reg) a.as_ptr(),
                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        inv = const $inv,
                        out("rax") _,
                        out("rcx") _,
                        out("rdx") _,
                        out("r8") _,
                        out("r9") _,
                        out("r10") _,
                        out("r11") _,
                        out("r12") r0,
                        out("r13") r1,
                        out("r14") r2,
                        out("r15") r3,
                        options(pure, readonly, nostack)
                    )
                }

                $field([r0, r1, r2, r3])
            }

            /// Multiplies `rhs` by `self`, returning the result.
            #[inline]
            pub fn mul(&self, rhs: &Self) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;
                unsafe {
                    asm!(
                        // schoolbook multiplication
                        //    *    |   a0    |   a1    |   a2    |   a3
                        //    b0   | b0 * a0 | b0 * a1 | b0 * a2 | b0 * a3
                        //    b1   | b1 * a0 | b1 * a1 | b1 * a2 | b1 * a3
                        //    b2   | b2 * a0 | b2 * a1 | b2 * a2 | b2 * a3
                        //    b3   | b3 * a0 | b3 * a1 | b3 * a2 | b3 * a3

                        // load value to registers
                        "mov r13, qword ptr [{b_ptr} + 0]",
                        "mov r14, qword ptr [{b_ptr} + 8]",
                        "mov r15, qword ptr [{b_ptr} + 16]",

                        // `a0`
                        "mov rdx, qword ptr [{a_ptr} + 0]",

                        // a0 * b0
                        "mulx r9, r8, r13",

                        // a0 * b1
                        "mulx r10, rax, r14",
                        "add r9, rax",

                        // a0 * b2
                        "mulx r11, rax, r15",
                        "adcx r10, rax",

                        // a0 * b3
                        "mulx r12, rax, qword ptr [{b_ptr} + 24]",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // `a1`
                        "mov rdx, [{a_ptr} + 8]",

                        // a1 * b0
                        "mulx rcx, rax, r13",
                        "add r9, rax",
                        "adcx r10, rcx",
                        "adc r11, 0",

                        // a1 * b1
                        "mulx rcx, rax, r14",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",
                        "xor r13, r13",

                        // a1 * b2
                        "mulx rcx, rax, r15",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",
                        "xor r14, r14",

                        // a1 * b3
                        "mulx rcx, rax, qword ptr [{b_ptr} + 24]",
                        "add r12, rax",
                        "adcx r13, rcx",
                        "adc r14, 0",

                        // `a2`
                        "mov rdx, [{a_ptr} + 16]",

                        // a2 * b0
                        "mulx rcx, rax, qword ptr [{b_ptr} + 0]",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",

                        // a2 * b1
                        "mulx rcx, rax, qword ptr [{b_ptr} + 8]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // a2 * b2
                        "mulx rcx, rax, r15",
                        "add r12, rax",
                        "adcx r13, rcx",
                        "adc r14, 0",
                        "xor r15, r15",

                        // a2 * b3
                        "mulx rcx, rax, qword ptr [{b_ptr} + 24]",
                        "add r13, rax",
                        "adcx r14, rcx",
                        "adc r15, 0",

                        // `a3`
                        "mov rdx, [{a_ptr} + 24]",

                        // a3 * b0
                        "mulx rcx, rax, qword ptr [{b_ptr} + 0]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // a3 * b1
                        "mulx rcx, rax, qword ptr [{b_ptr} + 8]",
                        "add r12, rax",
                        "adcx r13, rcx",
                        "adc r14, 0",

                        // a3 * b2
                        "mulx rcx, rax, qword ptr [{b_ptr} + 16]",
                        "add r13, rax",
                        "adcx r14, rcx",
                        "adc r15, 0",

                        // a3 * b3
                        "mulx rcx, rax, qword ptr [{b_ptr} + 24]",
                        "add r14, rax",
                        "adc r15, rcx",

                        // montgomery reduction
                        // r8 ~ r15

                        // `r8` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r8",

                        // r8' * m0
                        "mulx rcx, rax, qword ptr [{m_ptr} + 0]",
                        "add r8, rax",
                        "adcx r9, rcx",
                        "adc r10, 0",

                        // r8' * m1
                        "mulx rcx, rax, qword ptr [{m_ptr} + 8]",
                        "add r9, rax",
                        "adcx r10, rcx",
                        "adc r11, 0",

                        // // r8' * m2
                        "mulx rcx, rax, qword ptr [{m_ptr} + 16]",
                        "add r10, rax",
                        "adcx r11, rcx",
                        "adc r12, 0",

                        // // r8' * m3
                        "mulx rcx, rax, qword ptr [{m_ptr} + 24]",
                        "add r11, rax",
                        "adcx r12, rcx",
                        "adc r13, 0",

                        // `r9` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r9",

                        // r9' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r9, rcx",
                        "adcx r10, rax",
                        "adc r11, 0",

                        // r9' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r10, rcx",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // r9' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r9' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // `r10` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r10",

                        // r10' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r10, rcx",
                        "adcx r11, rax",
                        "adc r12, 0",

                        // r10' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r10' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // r10' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r13, rcx",
                        "adcx r14, rax",
                        "adc r15, 0",

                        // `r11` -> 0
                        "mov rdx, {inv}",
                        "mulx rax, rdx, r11",

                        // r11' * m0
                        "mulx rax, rcx, qword ptr [{m_ptr} + 0]",
                        "add r11, rcx",
                        "adcx r12, rax",
                        "adc r13, 0",

                        // r11' * m1
                        "mulx rax, rcx, qword ptr [{m_ptr} + 8]",
                        "add r12, rcx",
                        "adcx r13, rax",
                        "adc r14, 0",

                        // r11' * m2
                        "mulx rax, rcx, qword ptr [{m_ptr} + 16]",
                        "add r13, rcx",
                        "adcx r14, rax",
                        "adc r15, 0",

                        // r11' * m3
                        "mulx rax, rcx, qword ptr [{m_ptr} + 24]",
                        "add r14, rcx",
                        "adcx r15, rax",

                        // reduction if limbs is greater then mod
                        "mov r8, r12",
                        "mov r9, r13",
                        "mov r10, r14",
                        "mov r11, r15",

                        "sub r8, qword ptr [{m_ptr} + 0]",
                        "sbb r9, qword ptr [{m_ptr} + 8]",
                        "sbb r10, qword ptr [{m_ptr} + 16]",
                        "sbb r11, qword ptr [{m_ptr} + 24]",

                        "cmovc r8, r12",
                        "cmovc r9, r13",
                        "cmovc r10, r14",
                        "cmovc r11, r15",

                        "mov r12, r8",
                        "mov r13, r9",
                        "mov r14, r10",
                        "mov r15, r11",

                        "sub r12, qword ptr [{m_ptr} + 0]",
                        "sbb r13, qword ptr [{m_ptr} + 8]",
                        "sbb r14, qword ptr [{m_ptr} + 16]",
                        "sbb r15, qword ptr [{m_ptr} + 24]",

                        "cmovc r12, r8",
                        "cmovc r13, r9",
                        "cmovc r14, r10",
                        "cmovc r15, r11",

                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        a_ptr = in(reg) self.0.as_ptr(),
                        b_ptr = in(reg) rhs.0.as_ptr(),
                        inv = const $inv,
                        out("rax") _,
                        out("rcx") _,
                        out("rdx") _,
                        out("r8") _,
                        out("r9") _,
                        out("r10") _,
                        out("r11") _,
                        out("r12") r0,
                        out("r13") r1,
                        out("r14") r2,
                        out("r15") r3,
                        options(pure, readonly, nostack)
                    )
                }

                $field([r0, r1, r2, r3])
            }

            /// Subtracts `rhs` from `self`, returning the result.
            #[inline]
            pub fn sub(&self, rhs: &Self) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;
                unsafe {
                    asm!(
                        // init modulus area
                        "xor r12, r12",
                        "xor r13, r13",
                        "xor r14, r14",
                        "xor r15, r15",

                        // load a array to former registers
                        "mov r8, qword ptr [{a_ptr} + 0]",
                        "mov r9, qword ptr [{a_ptr} + 8]",
                        "mov r10, qword ptr [{a_ptr} + 16]",
                        "mov r11, qword ptr [{a_ptr} + 24]",

                        // sub a array and b array with borrow
                        "sub r8, qword ptr [{b_ptr} + 0]",
                        "sbb r9, qword ptr [{b_ptr} + 8]",
                        "sbb r10, qword ptr [{b_ptr} + 16]",
                        "sbb r11, qword ptr [{b_ptr} + 24]",

                        // if carry copy modulus
                        "cmovc r12, qword ptr [{m_ptr} + 0]",
                        "cmovc r13, qword ptr [{m_ptr} + 8]",
                        "cmovc r14, qword ptr [{m_ptr} + 16]",
                        "cmovc r15, qword ptr [{m_ptr} + 24]",

                        // mod addition
                        "add  r12, r8",
                        "adcx  r13, r9",
                        "adcx  r14, r10",
                        "adcx  r15, r11",

                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        a_ptr = in(reg) self.0.as_ptr(),
                        b_ptr = in(reg) rhs.0.as_ptr(),
                        out("r8") _,
                        out("r9") _,
                        out("r10") _,
                        out("r11") _,
                        out("r12") r0,
                        out("r13") r1,
                        out("r14") r2,
                        out("r15") r3,
                        options(pure, readonly, nostack)
                    );
                }
                $field([r0, r1, r2, r3])
            }

            /// Adds `rhs` to `self`, returning the result.
            #[inline]
            pub fn add(&self, rhs: &Self) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;
                unsafe {
                    asm!(
                        // load a array to former registers
                        "mov r8, qword ptr [{a_ptr} + 0]",
                        "mov r9, qword ptr [{a_ptr} + 8]",
                        "mov r10, qword ptr [{a_ptr} + 16]",
                        "mov r11, qword ptr [{a_ptr} + 24]",

                        // add a array and b array with carry
                        "add r8, qword ptr [{b_ptr} + 0]",
                        "adcx r9, qword ptr [{b_ptr} + 8]",
                        "adcx r10, qword ptr [{b_ptr} + 16]",
                        "adcx r11, qword ptr [{b_ptr} + 24]",

                        // copy result array to latter registers
                        "mov r12, r8",
                        "mov r13, r9",
                        "mov r14, r10",
                        "mov r15, r11",

                        // mod reduction
                        "sub r12, qword ptr [{m_ptr} + 0]",
                        "sbb r13, qword ptr [{m_ptr} + 8]",
                        "sbb r14, qword ptr [{m_ptr} + 16]",
                        "sbb r15, qword ptr [{m_ptr} + 24]",

                        // if carry copy former registers to out areas
                        "cmovc r12, r8",
                        "cmovc r13, r9",
                        "cmovc r14, r10",
                        "cmovc r15, r11",

                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        a_ptr = in(reg) self.0.as_ptr(),
                        b_ptr = in(reg) rhs.0.as_ptr(),
                        out("r8") _,
                        out("r9") _,
                        out("r10") _,
                        out("r11") _,
                        out("r12") r0,
                        out("r13") r1,
                        out("r14") r2,
                        out("r15") r3,
                        options(pure, readonly, nostack)
                    );
                }
                $field([r0, r1, r2, r3])
            }

            /// Negates `self`.
            #[inline]
            pub fn neg(&self) -> $field {
                let mut r0: u64;
                let mut r1: u64;
                let mut r2: u64;
                let mut r3: u64;
                unsafe {
                    asm!(
                        // load a array to former registers
                        "mov r8, qword ptr [{m_ptr} + 0]",
                        "mov r9, qword ptr [{m_ptr} + 8]",
                        "mov r10, qword ptr [{m_ptr} + 16]",
                        "mov r11, qword ptr [{m_ptr} + 24]",

                        "sub r8, qword ptr [{a_ptr} + 0]",
                        "sbb r9, qword ptr [{a_ptr} + 8]",
                        "sbb r10, qword ptr [{a_ptr} + 16]",
                        "sbb r11, qword ptr [{a_ptr} + 24]",

                        "mov r12, qword ptr [{a_ptr} + 0]",
                        "mov r13, qword ptr [{a_ptr} + 8]",
                        "mov r14, qword ptr [{a_ptr} + 16]",
                        "mov r15, qword ptr [{a_ptr} + 24]",

                        "or r12, r13",
                        "or r14, r15",
                        "or r12, r14",

                        "mov r13, 0xffffffffffffffff",
                        "cmp r12, 0x0000000000000000",
                        "cmove r13, r12",

                        "and r8, r13",
                        "and r9, r13",
                        "and r10, r13",
                        "and r11, r13",

                        a_ptr = in(reg) self.0.as_ptr(),
                        m_ptr = in(reg) $modulus.0.as_ptr(),
                        out("r8") r0,
                        out("r9") r1,
                        out("r10") r2,
                        out("r11") r3,
                        out("r12") _,
                        out("r13") _,
                        out("r14") _,
                        out("r15") _,
                        options(pure, readonly, nostack)
                    )
                }
                $field([r0, r1, r2, r3])
            }
        }
    };
}

pub(crate) use assembly_field;
