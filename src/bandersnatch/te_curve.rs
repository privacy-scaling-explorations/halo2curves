use crate::bandersnatch::Fp;
use crate::bandersnatch::Fr;
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use ff::{BatchInverter, Field, PrimeField};
use group::{self, Curve};
use group::{prime::PrimeCurveAffine, GroupEncoding};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

// `bandersnatch` is an incomplete twisted Edwards curve. These curves have
// equations of the form: ax² + y² = 1 + dx²y².
// over some base finite field Fp.
//
// Bandersnatch's curve equation: -5x² + y² = 1 + dx²y²
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

const FR_MODULUS_BYTES: [u8; 32] = [
    183, 44, 247, 214, 94, 14, 151, 208, 130, 16, 200, 204, 147, 32, 104, 166, 0, 59, 52, 1, 1, 59,
    103, 6, 169, 175, 51, 101, 234, 180, 125, 14,
];

use crate::{
    impl_add_binop_specify_output, impl_binops_additive, impl_binops_additive_specify_output,
    impl_binops_multiplicative, impl_binops_multiplicative_mixed, impl_sub_binop_specify_output,
};

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
pub struct BandersnatchTE {
    pub x: Fp,
    pub y: Fp,
    pub z: Fp,
    pub t: Fp,
}

#[derive(Copy, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "derive_serde", derive(Serialize, Deserialize))]
pub struct BandersnatchTEAffine {
    pub x: Fp,
    pub y: Fp,
}

#[derive(Copy, Clone, Hash, Default)]
pub struct BandersnatchTECompressed([u8; 32]);

impl BandersnatchTE {
    /// Constructs an extended point from the neutral element `(0, 1)`.
    pub const fn identity() -> Self {
        BandersnatchTE {
            x: Fp::zero(),
            y: Fp::one(),
            z: Fp::one(),
            t: Fp::zero(),
        }
    }

    /// Determines if this point is the identity.
    pub fn is_identity(&self) -> Choice {
        // If this point is the identity, then
        //     u = 0 * z = 0
        // and v = 1 * z = z
        self.x.ct_eq(&Fp::zero()) & self.y.ct_eq(&self.z)
    }

    /// Determines if this point is torsion free and so is contained
    /// in the prime order subgroup.
    pub fn is_torsion_free(&self) -> Choice {
        self.multiply(&FR_MODULUS_BYTES).is_identity()
    }

    #[inline]
    fn multiply(&self, by: &[u8; 32]) -> BandersnatchTE {
        let zero = BandersnatchTE::identity();
        let mut acc = BandersnatchTE::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the leading three bits because they're always
        // unset for Fr.
        for bit in by
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(3)
        {
            acc = acc.double();
            acc += BandersnatchTE::conditional_select(&zero, self, bit);
        }

        acc
    }

    /// Multiplies this element by the cofactor `4`.
    pub fn mul_by_cofactor(&self) -> BandersnatchTE {
        self.double().double()
    }

    pub fn generator() -> Self {
        let generator = BandersnatchTEAffine::generator();
        Self {
            x: generator.x,
            y: generator.y,
            z: Fp::one(),
            t: generator.x * generator.y,
        }
    }

    pub fn double(&self) -> BandersnatchTE {
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

        let a = self.x.square();
        let b = self.y.square();
        let c = self.z.square().double();
        let d = TE_A_PARAMETER * a;
        let e = (self.x + self.y).square() - a - b;
        let g = d + b;
        let f = g - c;
        let h = d - b;

        BandersnatchTE {
            x: e * f,
            y: g * h,
            z: f * g,
            t: e * h,
        }
    }
}

impl BandersnatchTEAffine {
    /// Constructs the neutral element `(0, 1)`.
    pub const fn identity() -> Self {
        BandersnatchTEAffine {
            x: Fp::zero(),
            y: Fp::one(),
        }
    }

    /// Determines if this point is the identity.
    pub fn is_identity(&self) -> Choice {
        BandersnatchTE::from(*self).is_identity()
    }

    pub fn generator() -> Self {
        Self {
            x: TE_BANDERSNATCH_GENERATOR_X,
            y: TE_BANDERSNATCH_GENERATOR_Y,
        }
    }

    pub fn to_extended(&self) -> BandersnatchTE {
        BandersnatchTE {
            x: self.x,
            y: self.y,
            z: Fp::one(),
            t: self.x * self.y,
        }
    }

    // TODO: fix this to work with BandersnatchTE
    pub fn random(mut rng: impl RngCore) -> Self {
        loop {
            let y = Fp::random(&mut rng);
            let flip_sign = rng.next_u32() % 2 != 0;

            let y2 = y.square();
            let p = ((y2 - Fp::one())
                * ((Fp::one() + TE_D_PARAMETER * y2)
                    .invert()
                    .unwrap_or(Fp::zero())))
            .sqrt()
            .map(|x| BandersnatchTEAffine {
                x: if flip_sign { -x } else { x },
                y,
            });

            if p.is_some().into() {
                use crate::group::cofactor::CofactorGroup;
                let p = p.unwrap().to_curve();

                if bool::from(!p.is_identity()) {
                    return p.clear_cofactor().to_affine();
                }
            }
        }
    }

    /// Converts this element into its byte representation.
    pub fn to_bytes(&self) -> [u8; 32] {
        let mut tmp = self.y.to_bytes();
        let u = self.x.to_bytes();

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
                    BandersnatchTEAffine { x: final_u, y: v },
                    !(zip_216_enabled & u_is_zero & flip_sign),
                )
            })
        })
    }
}

// Compressed
impl std::fmt::Debug for BandersnatchTECompressed {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.0[..].fmt(f)
    }
}

impl AsRef<[u8]> for BandersnatchTECompressed {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for BandersnatchTECompressed {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

// Jacobian implementations
impl<'a> From<&'a BandersnatchTEAffine> for BandersnatchTE {
    fn from(p: &'a BandersnatchTEAffine) -> BandersnatchTE {
        p.to_curve()
    }
}

impl From<BandersnatchTEAffine> for BandersnatchTE {
    fn from(p: BandersnatchTEAffine) -> BandersnatchTE {
        p.to_curve()
    }
}

impl Default for BandersnatchTE {
    fn default() -> BandersnatchTE {
        BandersnatchTE::identity()
    }
}

impl subtle::ConstantTimeEq for BandersnatchTE {
    fn ct_eq(&self, other: &Self) -> Choice {
        (self.x * other.z).ct_eq(&(other.x * self.z))
            & (self.y * other.z).ct_eq(&(other.y * self.z))
    }
}

impl subtle::ConditionallySelectable for BandersnatchTE {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        BandersnatchTE {
            x: Fp::conditional_select(&a.x, &b.x, choice),
            y: Fp::conditional_select(&a.y, &b.y, choice),
            z: Fp::conditional_select(&a.z, &b.z, choice),
            t: Fp::conditional_select(&a.t, &b.t, choice),
        }
    }
}

impl PartialEq for BandersnatchTE {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}

impl cmp::Eq for BandersnatchTE {}

impl CurveExt for BandersnatchTE {
    type ScalarExt = Fr;
    type Base = Fp;
    type AffineExt = BandersnatchTEAffine;

    const CURVE_ID: &'static str = "BandersnatchTE";

    fn is_on_curve(&self) -> Choice {
        let affine = BandersnatchTEAffine::from(*self);
        !self.z.is_zero() & affine.is_on_curve() & (affine.x * affine.y * self.z).ct_eq(&self.t)
    }

    fn endo(&self) -> Self {
        unimplemented!();
    }

    fn jacobian_coordinates(&self) -> (Fp, Fp, Fp) {
        unimplemented!();
    }

    fn hash_to_curve<'a>(_: &'a str) -> Box<dyn Fn(&[u8]) -> Self + 'a> {
        unimplemented!();
    }

    fn a() -> Self::Base {
        unimplemented!()
    }

    fn b() -> Self::Base {
        unimplemented!()
    }

    fn new_jacobian(_x: Self::Base, _y: Self::Base, _z: Self::Base) -> CtOption<Self> {
        unimplemented!();
    }
}

impl group::Curve for BandersnatchTE {
    type AffineRepr = BandersnatchTEAffine;

    fn batch_normalize(p: &[Self], q: &mut [Self::AffineRepr]) {
        assert_eq!(p.len(), q.len());

        for (p, q) in p.iter().zip(q.iter_mut()) {
            // We use the `u` field of `AffinePoint` to store the z-coordinate being
            // inverted, and the `v` field for scratch space.
            q.x = p.z;
        }

        BatchInverter::invert_with_internal_scratch(q, |q| &mut q.x, |q| &mut q.y);

        for (p, q) in p.iter().zip(q.iter_mut()).rev() {
            let tmp = q.x;

            // Set the coordinates to the correct value
            q.x = p.x * tmp; // Multiply by 1/z
            q.y = p.y * tmp; // Multiply by 1/z
        }
    }

    fn to_affine(&self) -> Self::AffineRepr {
        // Z coordinate is always nonzero, so this is
        // its inverse.
        let zinv = self.z.invert().unwrap();

        BandersnatchTEAffine {
            x: self.x * zinv,
            y: self.y * zinv,
        }
    }
}

impl group::Group for BandersnatchTE {
    type Scalar = Fr;

    fn random(mut rng: impl RngCore) -> Self {
        BandersnatchTEAffine::random(&mut rng).to_curve()
    }

    fn generator() -> Self {
        BandersnatchTE::generator()
    }

    fn identity() -> Self {
        Self::identity()
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    #[must_use]
    fn double(&self) -> Self {
        self.double()
    }
}

impl GroupEncoding for BandersnatchTE {
    type Repr = BandersnatchTECompressed;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        BandersnatchTEAffine::from_bytes(bytes.0).map(Self::from)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        BandersnatchTEAffine::from_bytes(bytes.0).map(Self::from)
    }

    fn to_bytes(&self) -> Self::Repr {
        BandersnatchTECompressed(BandersnatchTEAffine::from(self).to_bytes())
    }
}

impl crate::serde::SerdeObject for BandersnatchTE {
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        debug_assert_eq!(bytes.len(), 4 * Fp::size());
        let [x, y, z, t] = [0, 1, 2, 3]
            .map(|i| Fp::from_raw_bytes_unchecked(&bytes[i * Fp::size()..(i + 1) * Fp::size()]));
        Self { x, y, z, t }
    }
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 4 * Fp::size() {
            return None;
        }
        let [x, y, z, t] =
            [0, 1, 2, 3].map(|i| Fp::from_raw_bytes(&bytes[i * Fp::size()..(i + 1) * Fp::size()]));
        x.zip(y).zip(z).zip(t).and_then(|(((x, y), z), t)| {
            let res = Self { x, y, z, t };
            // Check that the point is on the curve.
            bool::from(res.is_on_curve()).then_some(res)
        })
    }
    fn to_raw_bytes(&self) -> Vec<u8> {
        let mut res = Vec::with_capacity(4 * Fp::size());
        Self::write_raw(self, &mut res).unwrap();
        res
    }
    fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
        let [x, y, z, t] = [(); 4].map(|_| Fp::read_raw_unchecked(reader));
        Self { x, y, z, t }
    }
    fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let x = Fp::read_raw(reader)?;
        let y = Fp::read_raw(reader)?;
        let z = Fp::read_raw(reader)?;
        let t = Fp::read_raw(reader)?;
        Ok(Self { x, y, z, t })
    }
    fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.x.write_raw(writer)?;
        self.y.write_raw(writer)?;
        self.z.write_raw(writer)?;
        self.t.write_raw(writer)
    }
}

impl group::prime::PrimeGroup for BandersnatchTE {}

impl group::prime::PrimeCurve for BandersnatchTE {
    type Affine = BandersnatchTEAffine;
}

impl group::cofactor::CofactorCurve for BandersnatchTE {
    type Affine = BandersnatchTEAffine;
}

impl group::cofactor::CofactorGroup for BandersnatchTE {
    type Subgroup = BandersnatchTE;

    fn clear_cofactor(&self) -> Self {
        self.mul_by_cofactor()
    }

    fn into_subgroup(self) -> CtOption<Self::Subgroup> {
        CtOption::new(self, self.is_torsion_free())
    }

    fn is_torsion_free(&self) -> Choice {
        self.is_torsion_free()
    }
}

impl<'a> From<&'a BandersnatchTE> for BandersnatchTEAffine {
    fn from(p: &'a BandersnatchTE) -> BandersnatchTEAffine {
        p.to_affine()
    }
}

impl From<BandersnatchTE> for BandersnatchTEAffine {
    fn from(p: BandersnatchTE) -> BandersnatchTEAffine {
        p.to_affine()
    }
}

impl Default for BandersnatchTEAffine {
    fn default() -> BandersnatchTEAffine {
        BandersnatchTEAffine::identity()
    }
}

impl subtle::ConstantTimeEq for BandersnatchTEAffine {
    fn ct_eq(&self, other: &Self) -> Choice {
        self.x.ct_eq(&other.x) & self.y.ct_eq(&other.y)
    }
}

impl subtle::ConditionallySelectable for BandersnatchTEAffine {
    fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
        BandersnatchTEAffine {
            x: Fp::conditional_select(&a.x, &b.x, choice),
            y: Fp::conditional_select(&a.y, &b.y, choice),
        }
    }
}

impl cmp::Eq for BandersnatchTEAffine {}

impl group::GroupEncoding for BandersnatchTEAffine {
    type Repr = [u8; 32];

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_bytes(*bytes)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_bytes(*bytes)
    }

    fn to_bytes(&self) -> Self::Repr {
        self.to_bytes()
    }
}

impl crate::serde::SerdeObject for BandersnatchTEAffine {
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        debug_assert_eq!(bytes.len(), 2 * Fp::size());
        let [x, y] =
            [0, Fp::size()].map(|i| Fp::from_raw_bytes_unchecked(&bytes[i..i + Fp::size()]));
        Self { x, y }
    }
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 2 * Fp::size() {
            return None;
        }
        let [x, y] = [0, Fp::size()].map(|i| Fp::from_raw_bytes(&bytes[i..i + Fp::size()]));
        x.zip(y).and_then(|(x, y)| {
            let res = Self { x, y };
            // Check that the point is on the curve.
            bool::from(res.is_on_curve()).then_some(res)
        })
    }
    fn to_raw_bytes(&self) -> Vec<u8> {
        let mut res = Vec::with_capacity(2 * Fp::size());
        Self::write_raw(self, &mut res).unwrap();
        res
    }
    fn read_raw_unchecked<R: std::io::Read>(reader: &mut R) -> Self {
        let [x, y] = [(); 2].map(|_| Fp::read_raw_unchecked(reader));
        Self { x, y }
    }
    fn read_raw<R: std::io::Read>(reader: &mut R) -> std::io::Result<Self> {
        let x = Fp::read_raw(reader)?;
        let y = Fp::read_raw(reader)?;
        Ok(Self { x, y })
    }
    fn write_raw<W: std::io::Write>(&self, writer: &mut W) -> std::io::Result<()> {
        self.x.write_raw(writer)?;
        self.y.write_raw(writer)
    }
}

impl group::prime::PrimeCurveAffine for BandersnatchTEAffine {
    type Curve = BandersnatchTE;
    type Scalar = Fr;

    fn generator() -> Self {
        BandersnatchTEAffine::generator()
    }

    fn identity() -> Self {
        BandersnatchTEAffine::identity()
    }

    fn is_identity(&self) -> Choice {
        self.is_identity()
    }

    fn to_curve(&self) -> Self::Curve {
        BandersnatchTE {
            x: self.x,
            y: self.y,
            z: Fp::one(),
            t: self.x * self.y,
        }
    }
}

impl group::cofactor::CofactorCurveAffine for BandersnatchTEAffine {
    type Curve = BandersnatchTE;
    type Scalar = Fr;

    fn identity() -> Self {
        <Self as group::prime::PrimeCurveAffine>::identity()
    }

    fn generator() -> Self {
        <Self as group::prime::PrimeCurveAffine>::generator()
    }

    fn is_identity(&self) -> Choice {
        <Self as group::prime::PrimeCurveAffine>::is_identity(self)
    }

    fn to_curve(&self) -> Self::Curve {
        <Self as group::prime::PrimeCurveAffine>::to_curve(self)
    }
}

impl CurveAffine for BandersnatchTEAffine {
    type ScalarExt = Fr;
    type Base = Fp;
    type CurveExt = BandersnatchTE;

    fn is_on_curve(&self) -> Choice {
        let x2 = self.x.square();
        let ax2 = self.x.square() * TE_A_PARAMETER;
        let y2 = self.y.square();

        (ax2 + y2).ct_eq(&(Fp::one() + TE_D_PARAMETER * x2 * y2))
    }

    fn coordinates(&self) -> CtOption<Coordinates<Self>> {
        Coordinates::from_xy(self.x, self.y)
    }

    fn from_xy(x: Self::Base, y: Self::Base) -> CtOption<Self> {
        let p = BandersnatchTEAffine { x, y };
        CtOption::new(p, p.is_on_curve())
    }

    fn a() -> Self::Base {
        unimplemented!()
    }

    fn b() -> Self::Base {
        unimplemented!()
    }
}

impl_binops_additive!(BandersnatchTE, BandersnatchTE);
impl_binops_additive!(BandersnatchTE, BandersnatchTEAffine);
impl_binops_additive_specify_output!(BandersnatchTEAffine, BandersnatchTEAffine, BandersnatchTE);
impl_binops_additive_specify_output!(BandersnatchTEAffine, BandersnatchTE, BandersnatchTE);
impl_binops_multiplicative!(BandersnatchTE, Fr);
impl_binops_multiplicative_mixed!(BandersnatchTEAffine, Fr, BandersnatchTE);

impl<'a> Neg for &'a BandersnatchTE {
    type Output = BandersnatchTE;

    fn neg(self) -> BandersnatchTE {
        BandersnatchTE {
            x: -self.x,
            y: self.y,
            z: self.z,
            t: -self.t,
        }
    }
}

impl Neg for BandersnatchTE {
    type Output = BandersnatchTE;

    fn neg(self) -> BandersnatchTE {
        -&self
    }
}

impl<T> Sum<T> for BandersnatchTE
where
    T: core::borrow::Borrow<BandersnatchTE>,
{
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = T>,
    {
        iter.fold(Self::identity(), |acc, item| acc + item.borrow())
    }
}

impl<'a, 'b> Add<&'a BandersnatchTE> for &'b BandersnatchTE {
    type Output = BandersnatchTE;

    fn add(self, rhs: &'a BandersnatchTE) -> BandersnatchTE {
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

        let a = self.x * rhs.x;
        let b = self.y * rhs.y;
        let c = self.t * TE_D_PARAMETER * rhs.t;
        let d = self.z * rhs.z;
        let e = (self.x + self.y) * (rhs.x + rhs.y) - a - b;
        let f = d - c;
        let g = d + c;
        let h = b - TE_A_PARAMETER * a;

        let _u_r = e * f;

        BandersnatchTE {
            x: e * f,
            y: g * h,
            z: f * g,
            t: e * h,
        }
    }
}

impl<'a, 'b> Add<&'a BandersnatchTEAffine> for &'b BandersnatchTE {
    type Output = BandersnatchTE;

    fn add(self, rhs: &'a BandersnatchTEAffine) -> BandersnatchTE {
        self + rhs.to_extended()
    }
}

impl<'a, 'b> Sub<&'a BandersnatchTE> for &'b BandersnatchTE {
    type Output = BandersnatchTE;

    fn sub(self, other: &'a BandersnatchTE) -> BandersnatchTE {
        self + (-other)
    }
}

impl<'a, 'b> Sub<&'a BandersnatchTEAffine> for &'b BandersnatchTE {
    type Output = BandersnatchTE;

    fn sub(self, other: &'a BandersnatchTEAffine) -> BandersnatchTE {
        self + (-other)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<'a, 'b> Mul<&'b Fr> for &'a BandersnatchTE {
    type Output = BandersnatchTE;

    // This is a simple double-and-add implementation of point
    // multiplication, moving from most significant to least
    // significant bit of the scalar.
    //
    // We skip the leading three bits because they're always
    // unset for Fr.
    fn mul(self, other: &'b Fr) -> Self::Output {
        let mut acc = BandersnatchTE::identity();
        for bit in other
            .to_repr()
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
            .skip(3)
        {
            acc = acc.double();
            acc = BandersnatchTE::conditional_select(&acc, &(acc + self), bit);
        }

        acc
    }
}

impl<'a> Neg for &'a BandersnatchTEAffine {
    type Output = BandersnatchTEAffine;

    fn neg(self) -> BandersnatchTEAffine {
        BandersnatchTEAffine {
            x: -self.x,
            y: self.y,
        }
    }
}

impl Neg for BandersnatchTEAffine {
    type Output = BandersnatchTEAffine;

    fn neg(self) -> BandersnatchTEAffine {
        -&self
    }
}

impl<'a, 'b> Add<&'a BandersnatchTE> for &'b BandersnatchTEAffine {
    type Output = BandersnatchTE;

    fn add(self, rhs: &'a BandersnatchTE) -> BandersnatchTE {
        rhs + self
    }
}

impl<'a, 'b> Add<&'a BandersnatchTEAffine> for &'b BandersnatchTEAffine {
    type Output = BandersnatchTE;

    fn add(self, rhs: &'a BandersnatchTEAffine) -> BandersnatchTE {
        self.to_extended() + rhs.to_extended()
    }
}

impl<'a, 'b> Sub<&'a BandersnatchTEAffine> for &'b BandersnatchTEAffine {
    type Output = BandersnatchTE;

    fn sub(self, other: &'a BandersnatchTEAffine) -> BandersnatchTE {
        self + (-other)
    }
}

impl<'a, 'b> Sub<&'a BandersnatchTE> for &'b BandersnatchTEAffine {
    type Output = BandersnatchTE;

    fn sub(self, other: &'a BandersnatchTE) -> BandersnatchTE {
        self + (-other)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<'a, 'b> Mul<&'b Fr> for &'a BandersnatchTEAffine {
    type Output = BandersnatchTE;

    fn mul(self, other: &'b Fr) -> Self::Output {
        let mut acc = BandersnatchTE::identity();

        // This is a simple double-and-add implementation of point
        // multiplication, moving from most significant to least
        // significant bit of the scalar.
        //
        // We skip the leading three bits because they're always
        // unset for Fr.
        for bit in other
            .to_repr()
            .iter()
            .rev()
            .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
        {
            acc = acc.double();
            acc = BandersnatchTE::conditional_select(&acc, &(acc + self), bit);
        }

        acc
    }
}

pub trait CurveAffineExt: pasta_curves::arithmetic::CurveAffine {
    /// Unlike the `Coordinates` trait, this just returns the raw affine coordinates without checking `is_on_curve`
    fn into_coordinates(self) -> (Self::Base, Self::Base) {
        // fallback implementation
        let coordinates = self.coordinates().unwrap();
        (*coordinates.x(), *coordinates.y())
    }
}

impl CurveAffineExt for BandersnatchTEAffine {
    fn into_coordinates(self) -> (Self::Base, Self::Base) {
        (self.x, self.y)
    }
}

pub trait TwistedEdwardsCurveExt: CurveExt {
    fn a() -> <Self as CurveExt>::Base;
    fn d() -> <Self as CurveExt>::Base;
}

impl TwistedEdwardsCurveExt for BandersnatchTE {
    fn a() -> Fp {
        TE_A_PARAMETER
    }

    fn d() -> Fp {
        TE_D_PARAMETER
    }
}

pub trait TwistedEdwardsCurveAffineExt: CurveAffineExt {
    fn a() -> <Self as CurveAffine>::Base;
    fn d() -> <Self as CurveAffine>::Base;
}

impl TwistedEdwardsCurveAffineExt for BandersnatchTEAffine {
    fn a() -> Fp {
        TE_A_PARAMETER
    }

    fn d() -> Fp {
        TE_D_PARAMETER
    }
}
#[cfg(test)]
mod tests {
    use crate::bandersnatch::BandersnatchTEAffine;
    use pasta_curves::arithmetic::CurveAffine;
    use crate::bandersnatch::te_curve::TE_D_PARAMETER;
    use ff::Field;
    use crate::bandersnatch::BandersnatchTE;
    use ff::PrimeField;
    use pasta_curves::arithmetic::CurveExt;
    use crate::bandersnatch::Fr;

    #[test]
    fn test_is_on_curve() {
        assert!(bool::from(BandersnatchTEAffine::identity().is_on_curve()));
        assert!(bool::from(BandersnatchTEAffine::generator().is_on_curve()));
    }

    #[test]
    fn test_d_is_non_quadratic_residue() {
        assert!(bool::from(TE_D_PARAMETER.sqrt().is_none()));
        assert!(bool::from((-TE_D_PARAMETER).sqrt().is_none()));
        assert!(bool::from(
            (-TE_D_PARAMETER).invert().unwrap().sqrt().is_none()
        ));
    }

    #[test]
    fn test_double() {
        let p = BandersnatchTE::generator();

        assert_eq!(p.double(), p + p);

        let double_g = p.double();

        // Value is taken from arworks implementation of bandersnatch.
        assert_eq!(
            double_g.x,
            bls12_381::Scalar::from_str_vartime(
                // in this context this is Bandersnatch base field although struct name is Scalar.
                "47509778783496412982820807418084268119503941123460587829794679458985081388520"
            )
            .unwrap()
        );
    }

    #[test]
    fn test_assoc() {
        let p = BandersnatchTE::generator().mul_by_cofactor();
        assert!(bool::from(p.is_on_curve()));

        assert_eq!(
            (p * Fr::from(1000u64)) * Fr::from(3938u64),
            p * (Fr::from(1000u64) * Fr::from(3938u64)),
        );
    }
    #[test]
    fn test_equality_scalar_mul_double_addition() {
        let generator = BandersnatchTEAffine::generator();

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

    // tests failing because of random() is wrong
    // #[test]
    // fn test_curve() {
    //     crate::tests::curve::curve_tests::<BandersnatchTE>();
    // }

    // tests failing
    // #[test]
    // fn test_serialization() {
    //     crate::tests::curve::random_serialization_test::<BandersnatchTE>();
    // }
}