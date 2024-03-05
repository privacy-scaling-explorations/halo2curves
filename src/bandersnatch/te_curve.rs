use crate::bandersnatch::Fp;
use crate::bandersnatch::Fr;

use crate::ff::Field;
use crate::group::{prime::PrimeCurveAffine, Group as _, GroupEncoding};

use core::borrow::Borrow;

use core::fmt;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};

use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
};
#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

// `bandersnatch` is an incomplete twisted Edwards curve. These curves have
// equations of the form: ax² + y² = 1 + dx²y².
// over some base finite field Fp.
//
// bandersnatch's curve equation: -5x² + y² = 1 + dx²y²
//
// q = 52435875175126190479447740508185965837690552500527637822603658699938581184513.
//
// a = -5.
// d = (138827208126141220649022263972958607803/
//     171449701953573178309673572579671231137) mod q
//   = 45022363124591815672509500913686876175488063829319466900776701791074614335719.
//
// Sage script to calculate these:
//
// ```text
// q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
// Fp = GF(q)
// d = (Fp(138827208126141220649022263972958607803)/Fp(171449701953573178309673572579671231137))
// ```
// These parameters and the sage script obtained from:
// <https://github.com/asanso/Bandersnatch/>

// The TE form generator is generated following Zcash's fashion:
//  "The generators of G1 and G2 are computed by finding the lexicographically
//   smallest valid x-coordinate, and its lexicographically smallest
//   y-coordinate and scaling it by the cofactor such that the result is not
//   the point at infinity."
// The SW form generator is the same TE generator converted into SW form,
// obtained from the scripts:
//   <https://github.com/zhenfeizhang/bandersnatch/blob/main/bandersnatch/script/bandersnatch.sage>

// Reference: https://eprint.iacr.org/2021/1152.pdf

// TE x: 29c132cc2c0b34c5743711777bbe42f32b79c022ad998465e1e71866a252ae18

const TE_BANDERSNATCH_GENERATOR_X: Fp = Fp::from_raw([
    0xe1e71866a252ae18,
    0x2b79c022ad998465,
    0x743711777bbe42f3,
    0x29c132cc2c0b34c5,
]);

// TE y: 2a6c669eda123e0f157d8b50badcd586358cad81eee464605e3167b6cc974166
const TE_BANDERSNATCH_GENERATOR_Y: Fp = Fp::from_raw([
    0x5e3167b6cc974166,
    0x358cad81eee46460,
    0x157d8b50badcd586,
    0x2a6c669eda123e0f,
]);

// 73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFEFFFFFFFC
const TE_A_PARAMETER: Fp = Fp::from_raw([
    0xFFFFFFFEFFFFFFFC,
    0x53BDA402FFFE5BFE,
    0x3339D80809A1D805,
    0x73EDA753299D7D48,
]);

// 6389C12633C267CBC66E3BF86BE3B6D8CB66677177E54F92B369F2F5188D58E7
const TE_D_PARAMETER: Fp = Fp::from_raw([
    0xB369F2F5188D58E7,
    0xCB66677177E54F92,
    0xC66E3BF86BE3B6D8,
    0x6389C12633C267CB,
]);

// TODO: change to bandersnatch if it's not correct, but should be
const FR_MODULUS_BYTES: [u8; 32] = [
    183, 44, 247, 214, 94, 14, 151, 208, 130, 16, 200, 204, 147, 32, 104, 166, 0, 59, 52, 1, 1, 59,
    103, 6, 169, 175, 51, 101, 234, 180, 125, 14,
];

/// This represents a Bandersnatch point in the affine `(u, v)`
/// coordinates.
#[derive(Clone, Copy, Debug, Eq)]
pub struct AffinePoint {
    u: Fp,
    v: Fp,
}

impl fmt::Display for AffinePoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Neg for AffinePoint {
    type Output = AffinePoint;

    /// This computes the negation of a point `P = (u, v)`
    /// as `-P = (-u, v)`.
    #[inline]
    fn neg(self) -> AffinePoint {
        AffinePoint {
            u: -self.u,
            v: self.v,
        }
    }
}

impl ConstantTimeEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.u.ct_eq(&other.u) & self.v.ct_eq(&other.v)
    }
}

impl PartialEq for AffinePoint {
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

#[derive(Clone, Copy, Debug, Eq)]
pub struct ExtendedPoint {
    u: Fp,
    v: Fp,
    z: Fp,
    t: Fp,
}

impl fmt::Display for ExtendedPoint {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl ConstantTimeEq for ExtendedPoint {
    fn ct_eq(&self, other: &Self) -> Choice {
        // x1/z1 == x2/z2  <==> x1 * z2 == x2 * z1
        (self.u * other.z).ct_eq(&(other.u * self.z))
            & (self.v * other.z).ct_eq(&(other.v * self.z))
    }
}

impl ConditionallySelectable for ExtendedPoint {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        ExtendedPoint {
            u: Fp::conditional_select(&a.u, &b.u, choice),
            v: Fp::conditional_select(&a.v, &b.v, choice),
            z: Fp::conditional_select(&a.z, &b.z, choice),
            t: Fp::conditional_select(&a.t, &b.t, choice),
        }
    }
}

impl PartialEq for ExtendedPoint {
    fn eq(&self, other: &Self) -> bool {
        bool::from(self.ct_eq(other))
    }
}

impl<T> Sum<T> for ExtendedPoint
where
    T: Borrow<ExtendedPoint>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

impl Neg for ExtendedPoint {
    type Output = ExtendedPoint;

    /// Computes the negation of a point `P = (U, V, Z, T)`
    /// as `-P = (-U, V, Z, -T)`.
    #[inline]
    fn neg(self) -> ExtendedPoint {
        ExtendedPoint {
            u: -self.u,
            v: self.v,
            z: self.z,
            t: -self.t,
        }
    }
}

impl ExtendedPoint {
    /// Constructs an extended point from the neutral element `(0, 1)`.
    pub const fn identity() -> Self {
        ExtendedPoint {
            u: Fp::zero(),
            v: Fp::one(),
            z: Fp::one(),
            t: Fp::zero(),
        }
    }
    /// Determines if this point is the identity.
    pub fn is_identity(&self) -> Choice {
        // If this point is the identity, then
        //     u = 0 * z = 0
        // and v = 1 * z = z
        self.u.ct_eq(&Fp::zero()) & self.v.ct_eq(&self.z)
    }

    /// Determines if this point is of small order.
    pub fn is_small_order(&self) -> Choice {
        // We only need to perform two doublings, since the 2-torsion
        // points are (0, 1) and (0, -1), and so we only need to check
        // that the u-coordinate of the result is zero to see if the
        // point is small order.
        self.double().double().u.ct_eq(&Fp::zero())
    }

    /// Determines if this point is torsion free and so is contained
    /// in the prime order subgroup.
    pub fn is_torsion_free(&self) -> Choice {
        self.multiply(&FR_MODULUS_BYTES).is_identity()
    }

    /// Determines if this point is prime order, or in other words that
    /// the smallest scalar multiplied by this point that produces the
    /// identity is `r`. This is equivalent to checking that the point
    /// is both torsion free and not the identity.
    pub fn is_prime_order(&self) -> Choice {
        self.is_torsion_free() & (!self.is_identity())
    }

    /// Multiplies this element by the cofactor `4`.
    pub fn mul_by_cofactor(&self) -> ExtendedPoint {
        self.double().double()
    }

    /// Computes the doubling of a point more efficiently than a point can
    /// be added to itself.
    pub fn double(&self) -> ExtendedPoint {
        // Doubling is more efficient (three multiplications, four squarings)
        // when we work within the projective coordinate space (U:Z, V:Z). We
        // rely on the most efficient formula, link: https://hyperelliptic.org/EFD/g1p/data/twisted/extended/doubling/dbl-2008-hwcd.

        // A = X1^2
        // B = Y1^2
        // C = 2*Z1^2
        // D = a*A
        // E = (X1+Y1)^2-A-B
        // G = D+B
        // F = G-C
        // H = D-B
        // X3 = E*F
        // Y3 = G*H
        // T3 = E*H
        // Z3 = F*G

        let a = self.u.square();
        let b = self.v.square();
        let c = self.z.square().double();
        let d = TE_A_PARAMETER * a;
        let e = (self.u + self.v).square() - a - b;
        let g = d + b;
        let f = g - c;
        let h = d - b;

        ExtendedPoint {
            u: e * f,
            v: g * h,
            z: f * g,
            t: e * h,
        }
    }

    #[inline]
    fn multiply(&self, by: &[u8; 32]) -> ExtendedPoint {
        let zero = ExtendedPoint::identity();

        let mut acc = ExtendedPoint::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the leading four bits because they're always
        // unset for Fr.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(3)
        {
            acc = acc.double();
            acc += ExtendedPoint::conditional_select(&zero, self, bit);
        }

        acc
    }
}

impl_binops_multiplicative!(ExtendedPoint, Fr);
impl_binops_additive!(ExtendedPoint, ExtendedPoint);

impl<'a, 'b> Mul<&'b Fr> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn mul(self, other: &'b Fr) -> ExtendedPoint {
        self.multiply(&other.to_bytes())
    }
}

impl<'a, 'b> Add<&'b ExtendedPoint> for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, other: &'b ExtendedPoint) -> ExtendedPoint {
        // formula: https://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html#addition-add-2008-hwcd

        // A = X1*X2
        // B = Y1*Y2
        // C = T1*d*T2
        // D = Z1*Z2
        // E = (X1+Y1)*(X2+Y2)-A-B
        // F = D-C
        // G = D+C
        // H = B-a*A
        // X3 = E*F
        // Y3 = G*H
        // T3 = E*H
        // Z3 = F*G

        let a = self.u * other.u;
        let b = self.v * other.v;
        let c = self.t * TE_D_PARAMETER * other.t;
        let d = self.z * other.z;
        let e = (self.u + self.v) * (other.u + other.v) - a - b;
        let f = d - c;
        let g = d + c;
        let h = b - TE_A_PARAMETER * a;

        let _u_r = e * f;

        ExtendedPoint {
            u: e * f,
            v: g * h,
            z: f * g,
            t: e * h,
        }
        // res.0.0 17016332093273884661
    }
}

impl<'a, 'b> Sub<&'a ExtendedPoint> for &'b ExtendedPoint {
    type Output = ExtendedPoint;

    #[allow(clippy::suspicious_arithmetic_impl)]
    fn sub(self, other: &'a ExtendedPoint) -> ExtendedPoint {
        self + (-other)
    }
}
impl<'a> Neg for &'a ExtendedPoint {
    type Output = ExtendedPoint;

    fn neg(self) -> ExtendedPoint {
        ExtendedPoint {
            u: -self.u,
            v: self.v,
            z: self.z,
            t: -self.t,
        }
    }
}

// impl<'a, 'b> Sub<&'a $name> for &'b $name {
//     type Output = $name;

//     fn sub(self, other: &'a $name) -> $name {
//         self + (-other)
//     }
// }

impl From<AffinePoint> for ExtendedPoint {
    /// Constructs an extended point (with `Z = 1`) from
    /// an affine point using the map `(u, v) => (u, v, 1, u, v)`.
    fn from(affine: AffinePoint) -> ExtendedPoint {
        ExtendedPoint {
            u: affine.u,
            v: affine.v,
            z: Fp::one(),
            t: affine.u * affine.v,
        }
    }
}

impl AffinePoint {
    /// Constructs the neutral element `(0, 1)`.
    pub const fn identity() -> Self {
        AffinePoint {
            u: Fp::zero(),
            v: Fp::one(),
        }
    }

    /// Determines if this point is the identity.
    pub fn is_identity(&self) -> Choice {
        ExtendedPoint::from(*self).is_identity()
    }

    /// Multiplies this point by the cofactor, producing an
    /// `ExtendedPoint`
    pub fn mul_by_cofactor(&self) -> ExtendedPoint {
        ExtendedPoint::from(*self).mul_by_cofactor()
    }

    /// Determines if this point is of small order.
    pub fn is_small_order(&self) -> Choice {
        ExtendedPoint::from(*self).is_small_order()
    }

    /// Determines if this point is torsion free and so is
    /// in the prime order subgroup.
    pub fn is_torsion_free(&self) -> Choice {
        ExtendedPoint::from(*self).is_torsion_free()
    }

    /// Determines if this point is prime order, or in other words that
    /// the smallest scalar multiplied by this point that produces the
    /// identity is `r`. This is equivalent to checking that the point
    /// is both torsion free and not the identity.
    pub fn is_prime_order(&self) -> Choice {
        let extended = ExtendedPoint::from(*self);
        extended.is_torsion_free() & (!extended.is_identity())
    }

    /// Converts this element into its byte representation.
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut tmp = self.v.to_bytes();
        let u = self.u.to_bytes();

        // Encode the sign of the u-coordinate in the most
        // significant bit.
        tmp[31] |= u[0] << 7;

        tmp
    }

    /// Attempts to interpret a byte representation of an
    /// affine point, failing if the element is not on
    /// the curve or non-canonical.
    pub fn from_bytes(b: [u8; 32]) -> CtOption<Self> {
        Self::from_bytes_inner(b, 1.into())
    }

    /// TODO: remove this as this mostly applies to Jubjub and not Bandersnatch
    /// Attempts to interpret a byte representation of an affine point, failing if the
    /// element is not on the curve.
    ///
    /// Most non-canonical encodings will also cause a failure. However, this API
    /// preserves (for use in consensus-critical protocols) a bug in the parsing code that
    /// caused two non-canonical encodings to be **silently** accepted:
    ///
    /// - (0, 1), which is the identity;
    /// - (0, -1), which is a point of order two.
    ///
    /// Each of these has a single non-canonical encoding in which the value of the sign
    /// bit is 1.
    ///
    /// See [ZIP 216](https://zips.z.cash/zip-0216) for a more detailed description of the
    /// bug, as well as its fix.
    pub fn from_bytes_pre_zip216_compatibility(b: [u8; 32]) -> CtOption<Self> {
        Self::from_bytes_inner(b, 0.into())
    }

    fn from_bytes_inner(mut b: [u8; 32], zip_216_enabled: Choice) -> CtOption<Self> {
        // Grab the sign bit from the representation
        let sign = b[31] >> 7;

        // Mask away the sign bit
        b[31] &= 0b0111_1111;

        // Interpret what remains as the v-coordinate
        Fp::from_bytes(&b).and_then(|v| {
            // -u^2 + v^2 = 1 + d.u^2.v^2
            // -u^2 = 1 + d.u^2.v^2 - v^2    (rearrange)
            // -u^2 - d.u^2.v^2 = 1 - v^2    (rearrange)
            // u^2 + d.u^2.v^2 = v^2 - 1     (flip signs)
            // u^2 (1 + d.v^2) = v^2 - 1     (factor)
            // u^2 = (v^2 - 1) / (1 + d.v^2) (isolate u^2)
            // We know that (1 + d.v^2) is nonzero for all v:
            //   (1 + d.v^2) = 0
            //   d.v^2 = -1
            //   v^2 = -(1 / d)   No solutions, as -(1 / d) is not a square

            let v2 = v.square();

            ((v2 - Fp::one())
                * ((Fp::one() + TE_D_PARAMETER * v2)
                    .invert()
                    .unwrap_or(Fp::zero())))
            .sqrt()
            .and_then(|u| {
                // Fix the sign of `u` if necessary
                let flip_sign = Choice::from((u.to_bytes()[0] ^ sign) & 1);
                let u_negated = -u;
                let final_u = Fp::conditional_select(&u, &u_negated, flip_sign);

                // If u == 0, flip_sign == sign_bit. We therefore want to reject the
                // encoding as non-canonical if all of the following occur:
                // - ZIP 216 is enabled
                // - u == 0
                // - flip_sign == true
                let u_is_zero = u.ct_eq(&Fp::zero());
                CtOption::new(
                    AffinePoint { u: final_u, v },
                    !(zip_216_enabled & u_is_zero & flip_sign),
                )
            })
        })
    }

    /// Attempts to interpret a batch of byte representations of affine points.
    ///
    /// Returns None for each element if it is not on the curve, or is non-canonical
    /// according to ZIP 216.
    #[cfg(feature = "alloc")]
    pub fn batch_from_bytes(items: impl Iterator<Item = [u8; 32]>) -> Vec<CtOption<Self>> {
        use ff::BatchInvert;

        #[derive(Clone, Copy, Default)]
        struct Item {
            sign: u8,
            v: Fp,
            numerator: Fp,
            denominator: Fp,
        }

        impl ConditionallySelectable for Item {
            fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
                Item {
                    sign: u8::conditional_select(&a.sign, &b.sign, choice),
                    v: Fp::conditional_select(&a.v, &b.v, choice),
                    numerator: Fp::conditional_select(&a.numerator, &b.numerator, choice),
                    denominator: Fp::conditional_select(&a.denominator, &b.denominator, choice),
                }
            }
        }

        let items: Vec<_> = items
            .map(|mut b| {
                // Grab the sign bit from the representation
                let sign = b[31] >> 7;

                // Mask away the sign bit
                b[31] &= 0b0111_1111;

                // Interpret what remains as the v-coordinate
                Fp::from_bytes(&b).map(|v| {
                    // -u^2 + v^2 = 1 + d.u^2.v^2
                    // -u^2 = 1 + d.u^2.v^2 - v^2    (rearrange)
                    // -u^2 - d.u^2.v^2 = 1 - v^2    (rearrange)
                    // u^2 + d.u^2.v^2 = v^2 - 1     (flip signs)
                    // u^2 (1 + d.v^2) = v^2 - 1     (factor)
                    // u^2 = (v^2 - 1) / (1 + d.v^2) (isolate u^2)
                    // We know that (1 + d.v^2) is nonzero for all v:
                    //   (1 + d.v^2) = 0
                    //   d.v^2 = -1
                    //   v^2 = -(1 / d)   No solutions, as -(1 / d) is not a square

                    let v2 = v.square();

                    Item {
                        v,
                        sign,
                        numerator: (v2 - Fp::one()),
                        denominator: Fp::one() + TE_D_PARAMETER * v2,
                    }
                })
            })
            .collect();

        let mut denominators: Vec<_> = items
            .iter()
            .map(|item| item.map(|item| item.denominator).unwrap_or(Fp::zero()))
            .collect();
        denominators.iter_mut().batch_invert();

        items
            .into_iter()
            .zip(denominators.into_iter())
            .map(|(item, inv_denominator)| {
                item.and_then(
                    |Item {
                         v, sign, numerator, ..
                     }| {
                        (numerator * inv_denominator).sqrt().and_then(|u| {
                            // Fix the sign of `u` if necessary
                            let flip_sign = Choice::from((u.to_bytes()[0] ^ sign) & 1);
                            let u_negated = -u;
                            let final_u = Fp::conditional_select(&u, &u_negated, flip_sign);

                            // If u == 0, flip_sign == sign_bit. We therefore want to reject the
                            // encoding as non-canonical if all of the following occur:
                            // - u == 0
                            // - flip_sign == true
                            let u_is_zero = u.ct_eq(&Fp::zero());
                            CtOption::new(AffinePoint { u: final_u, v }, !(u_is_zero & flip_sign))
                        })
                    },
                )
            })
            .collect()
    }

    /// Returns the `u`-coordinate of this point.
    pub fn get_u(&self) -> Fp {
        self.u
    }

    /// Returns the `v`-coordinate of this point.
    pub fn get_v(&self) -> Fp {
        self.v
    }

    /// Returns an `ExtendedPoint` for use in arithmetic operations.
    pub fn to_extended(&self) -> ExtendedPoint {
        ExtendedPoint {
            u: self.u,
            v: self.v,
            z: Fp::one(),
            t: self.u * self.v,
        }
    }

    /// Constructs an AffinePoint given `u` and `v` without checking
    /// that the point is on the curve.
    pub const fn from_raw_unchecked(u: Fp, v: Fp) -> AffinePoint {
        AffinePoint { u, v }
    }

    /// This is only for debugging purposes and not
    /// exposed in the public API. Checks that this
    /// point is on the curve.
    #[cfg(test)]
    fn is_on_curve_vartime(&self) -> bool {
        let u2 = self.u.square();
        let v2 = self.v.square();

        v2 - u2 == Fp::one() + TE_D_PARAMETER * u2 * v2
    }
}

#[cfg(test)]
mod tests {
    use ff::PrimeField;

    use super::{AffinePoint, TE_BANDERSNATCH_GENERATOR_X, TE_BANDERSNATCH_GENERATOR_Y};
    use crate::bandersnatch::Fr;
    use crate::bls12_381::Scalar;

    #[test]
    fn test_double() {
        let generator = AffinePoint {
            u: TE_BANDERSNATCH_GENERATOR_X,
            v: TE_BANDERSNATCH_GENERATOR_Y,
        };

        let proj_generator = generator.to_extended();

        let double_g = proj_generator.double();

        // Value is taken from arworks implementation of bandersnatch.
        assert_eq!(
            double_g.u,
            Scalar::from_str_vartime(
                "47509778783496412982820807418084268119503941123460587829794679458985081388520"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_equality_scalar_mul_double_addition() {
        let generator = AffinePoint {
            u: TE_BANDERSNATCH_GENERATOR_X,
            v: TE_BANDERSNATCH_GENERATOR_Y,
        };

        let proj_generator = generator.to_extended();

        let double_g = proj_generator + proj_generator;

        let double_g_double = proj_generator.double();

        let scalar_mul = proj_generator * Fr::from(2);

        assert!(double_g.eq(&double_g_double));
        assert!(double_g.eq(&scalar_mul));

        let minus = double_g - proj_generator;

        assert!(proj_generator.eq(&minus));

        let quadruple_g = proj_generator.double().double();

        let scalar_mul_4 = proj_generator * Fr::from(4);

        assert!(quadruple_g.eq(&scalar_mul_4));
    }
}
