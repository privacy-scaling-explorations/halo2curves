use super::fq::{Fq, FROBENIUS_COEFF_FQ2_C1, NEGATIVE_ONE};
use crate::arithmetic::BaseExt;

use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use ff::Field;
use rand::RngCore;
use std::cmp::Ordering;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

/// An element of Fq2, represented by c0 + c1 * u.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Default)]
pub struct Fq2 {
    pub c0: Fq,
    pub c1: Fq,
}

/// `Fq2` elements are ordered lexicographically.
impl Ord for Fq2 {
    #[inline(always)]
    fn cmp(&self, other: &Fq2) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl PartialOrd for Fq2 {
    #[inline(always)]
    fn partial_cmp(&self, other: &Fq2) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl ConditionallySelectable for Fq2 {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        Fq2 {
            c0: Fq::conditional_select(&a.c0, &b.c0, choice),
            c1: Fq::conditional_select(&a.c1, &b.c1, choice),
        }
    }
}

impl Fq2 {
    fn add(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.add(&other.c0),
            c1: self.c1.add(&other.c1),
        }
    }

    fn sub(&self, other: &Self) -> Self {
        Self {
            c0: self.c0.sub(&other.c0),
            c1: self.c1.sub(&other.c1),
        }
    }

    fn neg(&self) -> Self {
        Self {
            c0: self.c0.neg(),
            c1: self.c1.neg(),
        }
    }

    fn mul_assign(&mut self, other: &Self) {
        let mut t1 = self.c0 * other.c0;
        let mut t2 = self.c1 * other.c1;
        let mut t0 = self.c0 + self.c1;
        self.c1 = other.c0 + other.c1;
        self.c0 = t1 - t2;
        t1 += t2;
        t0 *= self.c1;
        self.c1 = t0 - t1;
    }

    fn mul(&self, other: &Self) -> Self {
        let mut t = other.clone();
        t.mul_assign(self);
        t
    }

    fn square_assign(&mut self) {
        let mut ab = self.c0 * self.c1;
        let mut c0c1 = self.c0 + self.c1;
        let mut c0 = -self.c1;
        c0 += self.c0;
        c0 *= c0c1;
        c0 -= ab;
        self.c1 = ab.double();
        self.c0 = c0 + ab;
    }

    fn square(&self) -> Self {
        let mut t = self.clone();
        t.square_assign();
        t
    }

    fn double(&self) -> Self {
        Fq2 {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }

    fn double_assign(&mut self) {
        self.c0 = self.c0.double();
        self.c1 = self.c1.double();
    }

    fn frobenius_map(&mut self, power: usize) {
        self.c1 *= &FROBENIUS_COEFF_FQ2_C1[power % 2];
    }

    /// Multiply this element by quadratic nonresidue 9 + u.
    pub fn mul_by_nonresidue(&mut self) {
        // (xi+y)(i+9) = (9x+y)i+(9y-x)
        let t0 = self.c0;
        let t1 = self.c1;

        // 8*x*i + 8*y
        self.double_assign();
        self.double_assign();
        self.double_assign();

        // 9*y
        self.c0 += &t0;
        // (9*y - x)
        self.c0 -= &t1;

        // (9*x)i
        self.c1 += &t1;
        // (9*x + y)
        self.c1 += &t0;
    }

    // Multiply this element by ξ where ξ=i+9
    pub fn mul_by_xi(&mut self) {
        // (xi+y)(i+9) = (9x+y)i+(9y-x)
        let t0 = self.c0;
        let t1 = self.c1;

        // 8*x*i + 8*y
        self.double_assign();
        self.double_assign();
        self.double_assign();

        // 9*y
        self.c0 += &t0;
        // (9*y - x)
        self.c0 -= &t1;

        // (9*x)i
        self.c1 += &t1;
        // (9*x + y)
        self.c1 += &t0;
    }

    /// Norm of Fq2 as extension field in i over Fq
    pub fn norm(&self) -> Fq {
        let mut t0 = self.c0;
        let mut t1 = self.c1;
        t0 = t0.square();
        t1 = t1.square();
        t1 + t0
    }

    // conjucate by negating c1
    pub fn conjugate(&mut self) {
        self.c1 -= self.c1;
    }

    fn invert(&self) -> CtOption<Self> {
        let mut t1 = self.c1;
        t1 = t1.square();
        let mut t0 = self.c0;
        t0 = t0.square();
        t0 += &t1;
        t0.invert().map(|t| {
            let mut tmp = Fq2 {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0 *= &t;
            tmp.c1 *= &t;
            tmp.c1 = -tmp.c1;

            tmp
        })
    }

    fn pow(&self, by: &[u64; 4]) -> Self {
        let mut res = Self::one();
        for e in by.iter().rev() {
            for i in (0..64).rev() {
                res = res.square();
                let mut tmp = res;
                tmp *= self;
                res.conditional_assign(&tmp, (((*e >> i) & 0x1) as u8).into());
            }
        }
        res
    }
}

impl<'a> Neg for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn neg(self) -> Fq2 {
        self.neg()
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    #[inline]
    fn neg(self) -> Fq2 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn sub(self, rhs: &'b Fq2) -> Fq2 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn add(self, rhs: &'b Fq2) -> Fq2 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fq2> for &'a Fq2 {
    type Output = Fq2;

    #[inline]
    fn mul(self, rhs: &'b Fq2) -> Fq2 {
        self.mul(rhs)
    }
}

impl Field for Fq2 {
    fn random(mut rng: impl RngCore) -> Self {
        let mut random_bytes = [0; 64];
        rng.fill_bytes(&mut random_bytes[..]);
        let c0 = Fq::from_bytes_wide(&random_bytes);
        random_bytes = [0; 64];
        rng.fill_bytes(&mut random_bytes[..]);
        let c1 = Fq::from_bytes_wide(&random_bytes);

        Fq2 { c0, c1 }
    }

    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn square(&self) -> Self {
        self.square()
    }

    fn double(&self) -> Self {
        self.double()
    }

    // fn sqrt(&self) -> Option<Self> {
    fn sqrt(&self) -> CtOption<Self> {
        // Algorithm 9, https://eprint.iacr.org/2012/685.pdf

        if self.is_zero() {
            CtOption::new(Self::zero(), Choice::from(1))
        } else {
            // a1 = self^((q - 3) / 4)
            let u: [u64; 4] = [
                0x4f082305b61f3f51,
                0x65e05aa45a1c72a3,
                0x6e14116da0605617,
                0x0c19139cb84c680a,
            ];
            let mut a1 = self.pow(&u);
            let mut alpha = a1;
            alpha.square();
            alpha.mul_assign(self);
            let mut a0 = alpha;
            a0.frobenius_map(1);
            a0.mul_assign(&alpha);

            let neg1 = Fq2 {
                c0: NEGATIVE_ONE,
                c1: Fq::zero(),
            };

            if a0 == neg1 {
                CtOption::new(a0, Choice::from(0))
            } else {
                a1.mul_assign(self);

                if alpha == neg1 {
                    a1.mul_assign(&Fq2 {
                        c0: Fq::zero(),
                        c1: Fq::one(),
                    });
                } else {
                    alpha += &Fq2::one();
                    // alpha = alpha^((q - 1) / 2)
                    let u: [u64; 4] = [
                        0x9e10460b6c3e7ea3,
                        0xcbc0b548b438e546,
                        0xdc2822db40c0ac2e,
                        0x183227397098d014,
                    ];
                    alpha = alpha.pow(&u);
                    a1.mul_assign(&alpha);
                }
                CtOption::new(a1, Choice::from(1))
            }
        }
    }

    fn invert(&self) -> CtOption<Self> {
        let mut t1 = self.c1;
        t1 = t1.square();
        let mut t0 = self.c0;
        t0 = t0.square();
        t0 += &t1;
        t0.invert().map(|t| {
            let mut tmp = Fq2 {
                c0: self.c0,
                c1: self.c1,
            };
            tmp.c0 *= &t;
            tmp.c1 *= &t;
            tmp.c1 = -tmp.c1;

            tmp
        })
    }
}

impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);
