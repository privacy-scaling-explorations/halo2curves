macro_rules! assembly_field {
    (
        $field:ident,
        $modulus:ident,
        $inv:ident,
        $modulus_str:ident,
        $two_inv:ident,
        $root_of_unity_inv:ident,
        $delta:ident,
        $zeta:ident,
        $r:ident,
        $r2:ident,
        $r3:ident
    ) => {
        use std::arch::asm;

        // impl $field {
            // /// Returns zero, the additive identity.
            // #[inline]
            // pub const fn zero() -> $field {
            //     $field([0, 0, 0, 0])
            // }

            // /// Returns one, the multiplicative identity.
            // #[inline]
            // pub const fn one() -> $field {
            //     $r
            // }

            // fn from_u512(limbs: [u64; 8]) -> $field {
            //     // We reduce an arbitrary 512-bit number by decomposing it into two 256-bit digits
            //     // with the higher bits multiplied by 2^256. Thus, we perform two reductions
            //     //
            //     // 1. the lower bits are multiplied by R^2, as normal
            //     // 2. the upper bits are multiplied by R^2 * 2^256 = R^3
            //     //
            //     // and computing their sum in the field. It remains to see that arbitrary 256-bit
            //     // numbers can be placed into Montgomery form safely using the reduction. The
            //     // reduction works so long as the product is less than R=2^256 multiplied by
            //     // the modulus. This holds because for any `c` smaller than the modulus, we have
            //     // that (2^256 - 1)*c is an acceptable product for the reduction. Therefore, the
            //     // reduction always works so long as `c` is in the field; in this case it is either the
            //     // constant `R2` or `R3`.
            //     let d0 = $field([limbs[0], limbs[1], limbs[2], limbs[3]]);
            //     let d1 = $field([limbs[4], limbs[5], limbs[6], limbs[7]]);
            //     // Convert to Montgomery form
            //     d0 * $r2 + d1 * $r3
            // }

            // /// Converts from an integer represented in little endian
            // /// into its (congruent) `$field` representation.
            // pub const fn from_raw(val: [u64; 4]) -> Self {
            //     // Multiplication
            //     let (r0, carry) = mac(0, val[0], $r2.0[0], 0);
            //     let (r1, carry) = mac(0, val[0], $r2.0[1], carry);
            //     let (r2, carry) = mac(0, val[0], $r2.0[2], carry);
            //     let (r3, r4) = mac(0, val[0], $r2.0[3], carry);

            //     let (r1, carry) = mac(r1, val[1], $r2.0[0], 0);
            //     let (r2, carry) = mac(r2, val[1], $r2.0[1], carry);
            //     let (r3, carry) = mac(r3, val[1], $r2.0[2], carry);
            //     let (r4, r5) = mac(r4, val[1], $r2.0[3], carry);

            //     let (r2, carry) = mac(r2, val[2], $r2.0[0], 0);
            //     let (r3, carry) = mac(r3, val[2], $r2.0[1], carry);
            //     let (r4, carry) = mac(r4, val[2], $r2.0[2], carry);
            //     let (r5, r6) = mac(r5, val[2], $r2.0[3], carry);

            //     let (r3, carry) = mac(r3, val[3], $r2.0[0], 0);
            //     let (r4, carry) = mac(r4, val[3], $r2.0[1], carry);
            //     let (r5, carry) = mac(r5, val[3], $r2.0[2], carry);
            //     let (r6, r7) = mac(r6, val[3], $r2.0[3], carry);

            //     // Montgomery reduction (first part)
            //     let k = r0.wrapping_mul($inv);
            //     let (_, carry) = mac(r0, k, $modulus.0[0], 0);
            //     let (r1, carry) = mac(r1, k, $modulus.0[1], carry);
            //     let (r2, carry) = mac(r2, k, $modulus.0[2], carry);
            //     let (r3, carry) = mac(r3, k, $modulus.0[3], carry);
            //     let (r4, carry2) = adc(r4, 0, carry);

            //     let k = r1.wrapping_mul($inv);
            //     let (_, carry) = mac(r1, k, $modulus.0[0], 0);
            //     let (r2, carry) = mac(r2, k, $modulus.0[1], carry);
            //     let (r3, carry) = mac(r3, k, $modulus.0[2], carry);
            //     let (r4, carry) = mac(r4, k, $modulus.0[3], carry);
            //     let (r5, carry2) = adc(r5, carry2, carry);

            //     let k = r2.wrapping_mul($inv);
            //     let (_, carry) = mac(r2, k, $modulus.0[0], 0);
            //     let (r3, carry) = mac(r3, k, $modulus.0[1], carry);
            //     let (r4, carry) = mac(r4, k, $modulus.0[2], carry);
            //     let (r5, carry) = mac(r5, k, $modulus.0[3], carry);
            //     let (r6, carry2) = adc(r6, carry2, carry);

            //     let k = r3.wrapping_mul($inv);
            //     let (_, carry) = mac(r3, k, $modulus.0[0], 0);
            //     let (r4, carry) = mac(r4, k, $modulus.0[1], carry);
            //     let (r5, carry) = mac(r5, k, $modulus.0[2], carry);
            //     let (r6, carry) = mac(r6, k, $modulus.0[3], carry);
            //     let (r7, _) = adc(r7, carry2, carry);

            //     // Montgomery reduction (sub part)
            //     let (d0, borrow) = sbb(r4, $modulus.0[0], 0);
            //     let (d1, borrow) = sbb(r5, $modulus.0[1], borrow);
            //     let (d2, borrow) = sbb(r6, $modulus.0[2], borrow);
            //     let (d3, borrow) = sbb(r7, $modulus.0[3], borrow);

            //     let (d0, carry) = adc(d0, $modulus.0[0] & borrow, 0);
            //     let (d1, carry) = adc(d1, $modulus.0[1] & borrow, carry);
            //     let (d2, carry) = adc(d2, $modulus.0[2] & borrow, carry);
            //     let (d3, _) = adc(d3, $modulus.0[3] & borrow, carry);

            //     $field([d0, d1, d2, d3])
            // }

            // /// Attempts to convert a little-endian byte representation of
            // /// a scalar into a `Fr`, failing if the input is not canonical.
            // pub fn from_bytes(bytes: &[u8; 32]) -> CtOption<$field> {
            //     <Self as ff::PrimeField>::from_repr(*bytes)
            // }

            // /// Converts an element of `Fr` into a byte representation in
            // /// little-endian byte order.
            // pub fn to_bytes(&self) -> [u8; 32] {
            //     <Self as ff::PrimeField>::to_repr(self)
            // }
        // }

        // impl Group for $field {
        //     type Scalar = Self;

        //     fn group_zero() -> Self {
        //         Self::zero()
        //     }
        //     fn group_add(&mut self, rhs: &Self) {
        //         *self += *rhs;
        //     }
        //     fn group_sub(&mut self, rhs: &Self) {
        //         *self -= *rhs;
        //     }
        //     fn group_scale(&mut self, by: &Self::Scalar) {
        //         *self *= *by;
        //     }
        // }

        // impl fmt::Debug for $field {
        //     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        //         let tmp = self.to_repr();
        //         write!(f, "0x")?;
        //         for &b in tmp.iter().rev() {
        //             write!(f, "{:02x}", b)?;
        //         }
        //         Ok(())
        //     }
        // }

        // impl Default for $field {
        //     #[inline]
        //     fn default() -> Self {
        //         Self::zero()
        //     }
        // }

        // impl From<bool> for $field {
        //     fn from(bit: bool) -> $field {
        //         if bit {
        //             $field::one()
        //         } else {
        //             $field::zero()
        //         }
        //     }
        // }

        // impl From<u64> for $field {
        //     fn from(val: u64) -> $field {
        //         $field([val, 0, 0, 0]) * $r2
        //     }
        // }

        // impl ConstantTimeEq for $field {
        //     fn ct_eq(&self, other: &Self) -> Choice {
        //         self.0[0].ct_eq(&other.0[0])
        //             & self.0[1].ct_eq(&other.0[1])
        //             & self.0[2].ct_eq(&other.0[2])
        //             & self.0[3].ct_eq(&other.0[3])
        //     }
        // }

        // impl core::cmp::Ord for $field {
        //     fn cmp(&self, other: &Self) -> core::cmp::Ordering {
        //         let left = self.to_repr();
        //         let right = other.to_repr();
        //         left.iter()
        //             .zip(right.iter())
        //             .rev()
        //             .find_map(|(left_byte, right_byte)| match left_byte.cmp(right_byte) {
        //                 core::cmp::Ordering::Equal => None,
        //                 res => Some(res),
        //             })
        //             .unwrap_or(core::cmp::Ordering::Equal)
        //     }
        // }

        // impl core::cmp::PartialOrd for $field {
        //     fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        //         Some(self.cmp(other))
        //     }
        // }

        // impl ConditionallySelectable for $field {
        //     fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        //         $field([
        //             u64::conditional_select(&a.0[0], &b.0[0], choice),
        //             u64::conditional_select(&a.0[1], &b.0[1], choice),
        //             u64::conditional_select(&a.0[2], &b.0[2], choice),
        //             u64::conditional_select(&a.0[3], &b.0[3], choice),
        //         ])
        //     }
        // }

        // impl<'a> Neg for &'a $field {
        //     type Output = $field;

        //     #[inline]
        //     fn neg(self) -> $field {
        //         self.neg()
        //     }
        // }

        // impl Neg for $field {
        //     type Output = $field;

        //     #[inline]
        //     fn neg(self) -> $field {
        //         -&self
        //     }
        // }

        // impl<'a, 'b> Sub<&'b $field> for &'a $field {
        //     type Output = $field;

        //     #[inline]
        //     fn sub(self, rhs: &'b $field) -> $field {
        //         self.sub(rhs)
        //     }
        // }

        // impl<'a, 'b> Add<&'b $field> for &'a $field {
        //     type Output = $field;

        //     #[inline]
        //     fn add(self, rhs: &'b $field) -> $field {
        //         self.add(rhs)
        //     }
        // }

        // impl<'a, 'b> Mul<&'b $field> for &'a $field {
        //     type Output = $field;

        //     #[inline]
        //     fn mul(self, rhs: &'b $field) -> $field {
        //         self.mul(rhs)
        //     }
        // }

        // impl From<$field> for [u8; 32] {
        //     fn from(value: $field) -> [u8; 32] {
        //         value.to_repr()
        //     }
        // }

        // impl<'a> From<&'a $field> for [u8; 32] {
        //     fn from(value: &'a $field) -> [u8; 32] {
        //         value.to_repr()
        //     }
        // }

        // impl FieldExt for $field {
        //     const MODULUS: &'static str = $modulus_str;
        //     const TWO_INV: Self = $two_inv;
        //     const ROOT_OF_UNITY_INV: Self = $root_of_unity_inv;
        //     const DELTA: Self = $delta;
        //     const ZETA: Self = $zeta;

        //     fn from_u128(v: u128) -> Self {
        //         $field::from_raw([v as u64, (v >> 64) as u64, 0, 0])
        //     }

        //     /// Converts a 512-bit little endian integer into
        //     /// a `$field` by reducing by the modulus.
        //     fn from_bytes_wide(bytes: &[u8; 64]) -> $field {
        //         $field::from_u512([
        //             u64::from_le_bytes(bytes[0..8].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[8..16].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[16..24].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[24..32].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[32..40].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[40..48].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[48..56].try_into().unwrap()),
        //             u64::from_le_bytes(bytes[56..64].try_into().unwrap()),
        //         ])
        //     }

        //     fn get_lower_128(&self) -> u128 {
        //         let tmp = $field::montgomery_reduce(&[
        //             self.0[0], self.0[1], self.0[2], self.0[3], 0, 0, 0, 0,
        //         ]);

        //         u128::from(tmp.0[0]) | (u128::from(tmp.0[1]) << 64)
        //     }
        // }

        // impl $crate::serde::SerdeObject for $field {
        //     fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        //         debug_assert_eq!(bytes.len(), 32);
        //         let inner =
        //             [0, 8, 16, 24].map(|i| u64::from_le_bytes(bytes[i..i + 8].try_into().unwrap()));
        //         Self(inner)
        //     }
        //     fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        //         if bytes.len() != 32 {
        //             return None;
        //         }
        //         let elt = Self::from_raw_bytes_unchecked(bytes);
        //         is_less_than(&elt.0, &$modulus.0).then(|| elt)
        //     }
        //     fn to_raw_bytes(&self) -> Vec<u8> {
        //         let mut res = Vec::with_capacity(32);
        //         for limb in self.0.iter() {
        //             res.extend_from_slice(&limb.to_le_bytes());
        //         }
        //         res
        //     }
        //     fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
        //         let inner = [(); 4].map(|_| {
        //             let mut buf = [0; 8];
        //             reader.read_exact(&mut buf).unwrap();
        //             u64::from_le_bytes(buf)
        //         });
        //         Self(inner)
        //     }
        //     fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        //         let mut inner = [0u64; 4];
        //         for limb in inner.iter_mut() {
        //             let mut buf = [0; 8];
        //             reader.read_exact(&mut buf)?;
        //             *limb = u64::from_le_bytes(buf);
        //         }
        //         let elt = Self(inner);
        //         is_less_than(&elt.0, &$modulus.0)
        //             .then(|| elt)
        //             .ok_or_else(|| {
        //                 std::io::Error::new(
        //                     std::io::ErrorKind::InvalidData,
        //                     "input number is not less than field modulus",
        //                 )
        //             })
        //     }
        //     fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        //         for limb in self.0.iter() {
        //             writer.write_all(&limb.to_le_bytes())?;
        //         }
        //         Ok(())
        //     }
        // }

        // /// Lexicographic comparison of Montgomery forms.
        // #[inline(always)]
        // fn is_less_than(x: &[u64; 4], y: &[u64; 4]) -> bool {
        //     match x[3].cmp(&y[3]) {
        //         core::cmp::Ordering::Less => return true,
        //         core::cmp::Ordering::Greater => return false,
        //         _ => {}
        //     }
        //     match x[2].cmp(&y[2]) {
        //         core::cmp::Ordering::Less => return true,
        //         core::cmp::Ordering::Greater => return false,
        //         _ => {}
        //     }
        //     match x[1].cmp(&y[1]) {
        //         core::cmp::Ordering::Less => return true,
        //         core::cmp::Ordering::Greater => return false,
        //         _ => {}
        //     }
        //     x[0].lt(&y[0])
        // }

        impl $field {
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
            pub(crate) fn montgomery_reduce(a: &[u64; 8]) -> $field {
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
