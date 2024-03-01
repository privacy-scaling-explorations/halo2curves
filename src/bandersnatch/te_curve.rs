/// `bandersnatch` is an incomplete twisted Edwards curve. These curves have
/// equations of the form: ax² + y² = 1 + dx²y².
/// over some base finite field Fq.
///
/// bandersnatch's curve equation: -5x² + y² = 1 + dx²y²
///
/// q = 52435875175126190479447740508185965837690552500527637822603658699938581184513.
///
/// a = -5.
/// d = (138827208126141220649022263972958607803/
///     171449701953573178309673572579671231137) mod q
///   = 45022363124591815672509500913686876175488063829319466900776701791074614335719.
///
/// Sage script to calculate these:
///
/// ```text
/// q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
/// Fq = GF(q)
/// d = (Fq(138827208126141220649022263972958607803)/Fq(171449701953573178309673572579671231137))
/// ```
/// These parameters and the sage script obtained from:
/// <https://github.com/asanso/Bandersnatch/>



use crate::ff::WithSmallOrderMulGroup;
use crate::ff::{Field, PrimeField};
use crate::group::{prime::PrimeCurveAffine, Curve, Group as _, GroupEncoding};
use crate::hash_to_curve::sswu_hash_to_curve;
use crate::bandersnatch::Fp;
use crate::bandersnatch::Fr;
use crate::{Coordinates, CurveAffine, CurveExt};
use core::cmp;
use core::fmt::Debug;
use core::iter::Sum;
use core::ops::{Add, Mul, Neg, Sub};
use rand::RngCore;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};

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
