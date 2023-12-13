#![allow(clippy::op_ref)]

use ff::{Field, FromUniformBytes, PrimeField};
use group::Group;
use pasta_curves::arithmetic::CurveExt;
use static_assertions::const_assert;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

use crate::{
    ff_ext::Legendre,
    secp256k1::{IsoSecp256k1, Secp256k1},
};

/// Hashes over a message and writes the output to all of `buf`.
/// Modified from https://github.com/zcash/pasta_curves/blob/7e3fc6a4919f6462a32b79dd226cb2587b7961eb/src/hashtocurve.rs#L11.
fn hash_to_field<F: FromUniformBytes<64>>(
    method: &str,
    curve_id: &str,
    domain_prefix: &str,
    message: &[u8],
    buf: &mut [F; 2],
) {
    assert!(domain_prefix.len() < 256);
    assert!((18 + method.len() + curve_id.len() + domain_prefix.len()) < 256);

    // Assume that the field size is 32 bytes and k is 256, where k is defined in
    // <https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html#name-security-considerations-3>.
    const CHUNKLEN: usize = 64;
    const_assert!(CHUNKLEN * 2 < 256);

    // Input block size of BLAKE2b.
    const R_IN_BYTES: usize = 128;

    let personal = [0u8; 16];
    let empty_hasher = blake2b_simd::Params::new()
        .hash_length(CHUNKLEN)
        .personal(&personal)
        .to_state();

    let b_0 = empty_hasher
        .clone()
        .update(&[0; R_IN_BYTES])
        .update(message)
        .update(&[0, (CHUNKLEN * 2) as u8, 0])
        .update(domain_prefix.as_bytes())
        .update(b"-")
        .update(curve_id.as_bytes())
        .update(b"_XMD:BLAKE2b_")
        .update(method.as_bytes())
        .update(b"_RO_")
        .update(&[(18 + method.len() + curve_id.len() + domain_prefix.len()) as u8])
        .finalize();

    let b_1 = empty_hasher
        .clone()
        .update(b_0.as_array())
        .update(&[1])
        .update(domain_prefix.as_bytes())
        .update(b"-")
        .update(curve_id.as_bytes())
        .update(b"_XMD:BLAKE2b_")
        .update(method.as_bytes())
        .update(b"_RO_")
        .update(&[(18 + method.len() + curve_id.len() + domain_prefix.len()) as u8])
        .finalize();

    let b_2 = {
        let mut empty_hasher = empty_hasher;
        for (l, r) in b_0.as_array().iter().zip(b_1.as_array().iter()) {
            empty_hasher.update(&[*l ^ *r]);
        }
        empty_hasher
            .update(&[2])
            .update(domain_prefix.as_bytes())
            .update(b"-")
            .update(curve_id.as_bytes())
            .update(b"_XMD:BLAKE2b_")
            .update(method.as_bytes())
            .update(b"_RO_")
            .update(&[(18 + method.len() + curve_id.len() + domain_prefix.len()) as u8])
            .finalize()
    };

    for (big, buf) in [b_1, b_2].iter().zip(buf.iter_mut()) {
        let mut little = [0u8; CHUNKLEN];
        little.copy_from_slice(big.as_array());
        little.reverse();
        *buf = F::from_uniform_bytes(&little);
    }
}

// Implementation of <https://datatracker.ietf.org/doc/html/rfc9380#name-simplified-swu-method>
#[allow(clippy::too_many_arguments)]
pub(crate) fn simple_svdw_map_to_curve<C>(u: C::Base, z: C::Base) -> C
where
    C: CurveExt,
{
    let zero = C::Base::ZERO;
    let one = C::Base::ONE;
    let a = C::a();
    let b = C::b();

    //1.  tv1 = u^2
    let tv1 = u.square();
    //2.  tv1 = Z * tv1
    let tv1 = z * tv1;
    //3.  tv2 = tv1^2
    let tv2 = tv1.square();
    //4.  tv2 = tv2 + tv1
    let tv2 = tv2 + tv1;
    //5.  tv3 = tv2 + 1
    let tv3 = tv2 + one;
    //6.  tv3 = B * tv3
    let tv3 = b * tv3;
    //7.  tv4 = CMOV(Z, -tv2, tv2 != 0) # tv4 = z if tv2 is 0 else tv4 = -tv2
    let tv2_is_not_zero = !tv2.ct_eq(&zero);
    let tv4 = C::Base::conditional_select(&z, &-tv2, tv2_is_not_zero);
    //8.  tv4 = A * tv4
    let tv4 = a * tv4;
    //9.  tv2 = tv3^2
    let tv2 = tv3.square();
    //10. tv6 = tv4^2
    let tv6 = tv4.square();
    //11. tv5 = A * tv6
    let tv5 = a * tv6;
    //12. tv2 = tv2 + tv5
    let tv2 = tv2 + tv5;
    //13. tv2 = tv2 * tv3
    let tv2 = tv2 * tv3;
    //14. tv6 = tv6 * tv4
    let tv6 = tv6 * tv4;
    //15. tv5 = B * tv6
    let tv5 = b * tv6;
    //16. tv2 = tv2 + tv5
    let tv2 = tv2 + tv5;
    //17.   x = tv1 * tv3
    let x = tv1 * tv3;
    //18. (is_gx1_square, y1) = sqrt_ratio(tv2, tv6)
    let (is_gx1_square, y1) = sqrt_ratio(&tv2, &tv6, &z);
    //19.   y = tv1 * u
    let y = tv1 * u;
    //20.   y = y * y1
    let y = y * y1;
    //21.   x = CMOV(x, tv3, is_gx1_square)
    let x = C::Base::conditional_select(&x, &tv3, is_gx1_square);
    //22.   y = CMOV(y, y1, is_gx1_square)
    let y = C::Base::conditional_select(&y, &y1, is_gx1_square);
    //23.  e1 = sgn0(u) == sgn0(y)
    let e1 = u.is_odd().ct_eq(&y.is_odd());
    //24.   y = CMOV(-y, y, e1) # Select correct sign of y
    let y = C::Base::conditional_select(&-y, &y, e1);
    //25.   x = x / tv4
    let x = x * tv4.invert().unwrap();
    //26. return (x, y)
    C::new_jacobian(x, y, one).unwrap()
}

// Implementation of <https://datatracker.ietf.org/doc/html/rfc9380#name-simplified-swu-method>
#[allow(clippy::type_complexity)]
pub(crate) fn simple_svdw_hash_to_curve<'a, C>(
    curve_id: &'static str,
    domain_prefix: &'a str,
    z: C::Base,
) -> Box<dyn Fn(&[u8]) -> C + 'a>
where
    C: CurveExt,
    C::Base: FromUniformBytes<64>,
{
    Box::new(move |message| {
        let mut us = [C::Base::ZERO; 2];
        hash_to_field("SSWU", curve_id, domain_prefix, message, &mut us);

        let [q0, q1]: [C; 2] = us.map(|u| simple_svdw_map_to_curve::<C>(u, z));

        let r = q0 + &q1;
        debug_assert!(bool::from(r.is_on_curve()));
        r
    })
}

// Implementation of <https://datatracker.ietf.org/doc/html/rfc9380#name-simplified-swu-for-ab-0>
#[allow(clippy::type_complexity)]
pub(crate) fn simple_svdw_hash_to_curve_secp256k1<'a>(
    _curve_id: &'static str,
    domain_prefix: &'a str,
) -> Box<dyn Fn(&[u8]) -> Secp256k1 + 'a> {
    Box::new(move |message| {
        let rp = IsoSecp256k1::hash_to_curve(domain_prefix)(message);

        let r = iso_map_secp256k1(rp);

        debug_assert!(bool::from(r.is_on_curve()));
        r
    })
}

/// 3-Isogeny Map for Secp256k1
/// Reference: <https://www.rfc-editor.org/rfc/rfc9380.html#name-3-isogeny-map-for-secp256k1>
pub fn iso_map_secp256k1(rp: IsoSecp256k1) -> Secp256k1 {
    // constants for secp256k1 iso_map computation
    const K: [[&str; 4]; 5] = [
        ["0x00", "0x00", "0x00", "0x00"],
        [
            "0x8e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38daaaaa8c7",
            "0x7d3d4c80bc321d5b9f315cea7fd44c5d595d2fc0bf63b92dfff1044f17c6581",
            "0x534c328d23f234e6e2a413deca25caece4506144037c40314ecbd0b53d9dd262",
            "0x8e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38e38daaaaa88c",
        ],
        [
            "0xd35771193d94918a9ca34ccbb7b640dd86cd409542f8487d9fe6b745781eb49b",
            "0xedadc6f64383dc1df7c4b2d51b54225406d36b641f5e41bbc52a56612a8c6d14",
            "0x00",
            "0x00",
        ],
        [
            "0x4bda12f684bda12f684bda12f684bda12f684bda12f684bda12f684b8e38e23c",
            "0xc75e0c32d5cb7c0fa9d0a54b12a0a6d5647ab046d686da6fdffc90fc201d71a3",
            "0x29a6194691f91a73715209ef6512e576722830a201be2018a765e85a9ecee931",
            "0x2f684bda12f684bda12f684bda12f684bda12f684bda12f684bda12f38e38d84",
        ],
        [
            "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffff93b",
            "0x7a06534bb8bdb49fd5e9e6632722c2989467c1bfc8e8d978dfb425d2685c2573",
            "0x6484aa716545ca2cf3a70c3fa8fe337e0a3d21162f0d6299a7bf8192bfd2a76f",
            "0x00",
        ],
    ];
    let mut k: [[<IsoSecp256k1 as CurveExt>::Base; 4]; 5] =
        [[<IsoSecp256k1 as CurveExt>::Base::from_uniform_bytes(&[0; 64]); 4]; 5];
    for i in 0..5 {
        for j in 0..4 {
            k[i][j] =
                <IsoSecp256k1 as CurveExt>::Base::from_uniform_bytes(&hex_str_to_le_bytes(K[i][j]));
        }
    }

    // convert to affine form
    let (xp, yp, zp) = rp.jacobian_coordinates();
    let (x, y) = jacobian_to_affine::<IsoSecp256k1>(xp, yp, zp);

    // iso_map logic
    let x_squared = x.square();
    let x_cubed = x * x_squared;

    let x_num = k[1][3] * x_cubed + k[1][2] * x_squared + k[1][1] * x + k[1][0];
    let x_den = x_squared + k[2][1] * x + k[2][0];

    let y_num = k[3][3] * x_cubed + k[3][2] * x_squared + k[3][1] * x + k[3][0];
    let y_den = x_cubed + k[4][2] * x_squared + k[4][1] * x + k[4][0];

    // exceptional case MUST return identity
    //   reference: <https://www.rfc-editor.org/rfc/rfc9380.html#name-simplified-swu-for-ab-0>
    if x_den.is_zero().into() || y_den.is_zero().into() {
        return Secp256k1::identity();
    }

    let x = x_num * x_den.invert().unwrap();
    let y = y * (y_num * y_den.invert().unwrap());

    Secp256k1::new_jacobian(x, y, <Secp256k1 as CurveExt>::Base::ONE).unwrap()
}

/// Converting a point from Jacobian coordinates to affine coordinates on an elliptic curve
fn jacobian_to_affine<C>(x: C::Base, y: C::Base, z: C::Base) -> (C::Base, C::Base)
where
    C: CurveExt,
{
    // identity
    if z.is_zero().into() {
        return (C::Base::ZERO, C::Base::ZERO);
    }

    let z_squared = z * z;
    let z_cubed = z_squared * z;

    let z_squared_inv = z_squared.invert().unwrap();
    let z_cubed_inv = z_cubed.invert().unwrap();

    (x * z_squared_inv, y * z_cubed_inv)
}

/// Convert hex string to little-endian bytes array of length `L`
///
/// NOTE: hex string should be prefixed with "0x"
///
/// Example:  
///
///   hex_str_to_le_bytes::<4>("0x01020304") -> [4, 3, 2, 1]
///
///   hex_str_to_le_bytes::<6>("0x01020304") -> [4, 3, 2, 1, 0, 0]
///
///   hex_str_to_le_bytes::<2>("0x01020304") -> [4, 3]
///
fn hex_str_to_le_bytes<const L: usize>(hex: &str) -> [u8; L] {
    let padded_hex_string = if hex.len() % 2 != 0 {
        format!("0{}", &hex[2..])
    } else {
        hex[2..].to_owned()
    };

    // Convert each pair of hex characters to u8 and collect into a vector
    let le_bytes: Result<Vec<u8>, _> = (0..padded_hex_string.len())
        .step_by(2)
        .rev() // Iterate in reverse order for little-endian byte order
        .map(|i| {
            u8::from_str_radix(&padded_hex_string[i..i + 2], 16)
                .map_err(|_| "Invalid hex character")
        })
        .collect();
    let le_bytes = le_bytes.expect("Invalid bytes");

    let mut result = [0; L];
    for i in 0..L.min(le_bytes.len()) {
        result[i] = le_bytes[i];
    }
    result
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn svdw_map_to_curve<C>(
    u: C::Base,
    c1: C::Base,
    c2: C::Base,
    c3: C::Base,
    c4: C::Base,
    z: C::Base,
) -> C
where
    C: CurveExt,
    C::Base: Legendre,
{
    let one = C::Base::ONE;
    let a = C::a();
    let b = C::b();

    // 1. tv1 = u^2
    let tv1 = u.square();
    // 2. tv1 = tv1 * c1
    let tv1 = tv1 * c1;
    // 3. tv2 = 1 + tv1
    let tv2 = one + tv1;
    // 4. tv1 = 1 - tv1
    let tv1 = one - tv1;
    // 5. tv3 = tv1 * tv2
    let tv3 = tv1 * tv2;
    // 6. tv3 = inv0(tv3)
    let tv3 = tv3.invert().unwrap_or(C::Base::ZERO);
    // 7. tv4 = u * tv1
    let tv4 = u * tv1;
    // 8. tv4 = tv4 * tv3
    let tv4 = tv4 * tv3;
    // 9. tv4 = tv4 * c3
    let tv4 = tv4 * c3;
    // 10. x1 = c2 - tv4
    let x1 = c2 - tv4;
    // 11. gx1 = x1^2
    let gx1 = x1.square();
    // 12. gx1 = gx1 + A
    let gx1 = gx1 + a;
    // 13. gx1 = gx1 * x1
    let gx1 = gx1 * x1;
    // 14. gx1 = gx1 + B
    let gx1 = gx1 + b;
    // 15. e1 = is_square(gx1)
    let e1 = !gx1.ct_quadratic_non_residue();
    // 16. x2 = c2 + tv4
    let x2 = c2 + tv4;
    // 17. gx2 = x2^2
    let gx2 = x2.square();
    // 18. gx2 = gx2 + A
    let gx2 = gx2 + a;
    // 19. gx2 = gx2 * x2
    let gx2 = gx2 * x2;
    // 20. gx2 = gx2 + B
    let gx2 = gx2 + b;
    // 21. e2 = is_square(gx2) AND NOT e1    # Avoid short-circuit logic ops
    let e2 = !gx2.ct_quadratic_non_residue() & (!e1);
    // 22. x3 = tv2^2
    let x3 = tv2.square();
    // 23. x3 = x3 * tv3
    let x3 = x3 * tv3;
    // 24. x3 = x3^2
    let x3 = x3.square();
    // 25. x3 = x3 * c4
    let x3 = x3 * c4;
    // 26. x3 = x3 + Z
    let x3 = x3 + z;
    // 27. x = CMOV(x3, x1, e1)    # x = x1 if gx1 is square, else x = x3
    let x = C::Base::conditional_select(&x3, &x1, e1);
    // 28. x = CMOV(x, x2, e2)    # x = x2 if gx2 is square and gx1 is not
    let x = C::Base::conditional_select(&x, &x2, e2);
    // 29. gx = x^2
    let gx = x.square();
    // 30. gx = gx + A
    let gx = gx + a;
    // 31. gx = gx * x
    let gx = gx * x;
    // 32. gx = gx + B
    let gx = gx + b;
    // 33. y = sqrt(gx)
    let y = gx.sqrt().unwrap();
    // 34. e3 = sgn0(u) == sgn0(y)
    let e3 = u.is_odd().ct_eq(&y.is_odd());
    // 35. y = CMOV(-y, y, e3)    # Select correct sign of y
    let y = C::Base::conditional_select(&-y, &y, e3);
    // 36. return (x, y)
    C::new_jacobian(x, y, one).unwrap()
}

// Implement https://datatracker.ietf.org/doc/html/rfc9380#name-sqrt_ratio-for-any-field
// Copied from ff sqrt_ratio_generic substituting F::ROOT_OF_UNITY for input Z
fn sqrt_ratio<F: PrimeField>(num: &F, div: &F, z: &F) -> (Choice, F) {
    // General implementation:
    //
    // a = num * inv0(div)
    //   = {    0    if div is zero
    //     { num/div otherwise
    //
    // b = z * a
    //   = {      0      if div is zero
    //     { z*num/div otherwise

    // Since z is non-square, a and b are either both zero (and both square), or
    // only one of them is square. We can therefore choose the square root to return
    // based on whether a is square, but for the boolean output we need to handle the
    // num != 0 && div == 0 case specifically.

    let a = div.invert().unwrap_or(F::ZERO) * num;
    let b = a * z;
    let sqrt_a = a.sqrt();
    let sqrt_b = b.sqrt();

    let num_is_zero = num.is_zero();
    let div_is_zero = div.is_zero();
    let is_square = sqrt_a.is_some();
    let is_nonsquare = sqrt_b.is_some();
    assert!(bool::from(
        num_is_zero | div_is_zero | (is_square ^ is_nonsquare)
    ));

    (
        is_square & (num_is_zero | !div_is_zero),
        CtOption::conditional_select(&sqrt_b, &sqrt_a, is_square).unwrap(),
    )
}

/// Implementation of https://www.ietf.org/archive/id/draft-irtf-cfrg-hash-to-curve-10.html#section-6.6.1
#[allow(clippy::type_complexity)]
pub(crate) fn svdw_hash_to_curve<'a, C>(
    curve_id: &'static str,
    domain_prefix: &'a str,
    z: C::Base,
) -> Box<dyn Fn(&[u8]) -> C + 'a>
where
    C: CurveExt,
    C::Base: FromUniformBytes<64> + Legendre,
{
    let [c1, c2, c3, c4] = svdw_precomputed_constants::<C>(z);

    Box::new(move |message| {
        let mut us = [C::Base::ZERO; 2];
        hash_to_field("SVDW", curve_id, domain_prefix, message, &mut us);

        let [q0, q1]: [C; 2] = us.map(|u| svdw_map_to_curve(u, c1, c2, c3, c4, z));

        let r = q0 + &q1;
        debug_assert!(bool::from(r.is_on_curve()));
        r
    })
}

pub(crate) fn svdw_precomputed_constants<C: CurveExt>(z: C::Base) -> [C::Base; 4] {
    let a = C::a();
    let b = C::b();
    let one = C::Base::ONE;
    let three = one + one + one;
    let four = three + one;
    let tmp = three * z.square() + four * a;

    // 1. c1 = g(Z)
    let c1 = (z.square() + a) * z + b;
    // 2. c2 = -Z / 2
    let c2 = -z * C::Base::TWO_INV;
    // 3. c3 = sqrt(-g(Z) * (3 * Z^2 + 4 * A))    # sgn0(c3) MUST equal 0
    let c3 = {
        let c3 = (-c1 * tmp).sqrt().unwrap();
        C::Base::conditional_select(&c3, &-c3, c3.is_odd())
    };
    // 4. c4 = -4 * g(Z) / (3 * Z^2 + 4 * A)
    let c4 = -four * c1 * tmp.invert().unwrap();

    [c1, c2, c3, c4]
}
