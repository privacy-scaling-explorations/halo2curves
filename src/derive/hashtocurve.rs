// This module implements "try and increment" hashing to short Weierstrass curves
// hash_to_field is from: https://github.com/zcash/pasta_curves

use crate::CurveExt;
use ff::{FromUniformBytes, PrimeField};
use static_assertions::const_assert;

// TODO: use simplified swu algorithm

pub fn try_and_increment<F: PrimeField, C: CurveExt<Base = F>>(t: &F) -> C {
    let mut x = *t;
    loop {
        let y2 = (x * x + C::a()) * x + C::b();
        let y = y2.sqrt();
        if bool::from(y.is_some()) {
            return C::new_jacobian(x, y.unwrap(), C::Base::ONE).unwrap();
        }
        x = x + C::Base::ONE;
    }
}

/// Hashes over a message and writes the output to all of `buf`.
pub fn hash_to_field<F: FromUniformBytes<64>>(
    curve_id: &str,
    domain_prefix: &str,
    message: &[u8],
    buf: &mut [F; 2],
) {
    assert!(domain_prefix.len() < 256);
    assert!((22 + curve_id.len() + domain_prefix.len()) < 256);

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
        .update(b"_XMD:BLAKE2b_SSWU_RO_")
        .update(&[(22 + curve_id.len() + domain_prefix.len()) as u8])
        .finalize();

    let b_1 = empty_hasher
        .clone()
        .update(b_0.as_array())
        .update(&[1])
        .update(domain_prefix.as_bytes())
        .update(b"-")
        .update(curve_id.as_bytes())
        .update(b"_XMD:BLAKE2b_SSWU_RO_")
        .update(&[(22 + curve_id.len() + domain_prefix.len()) as u8])
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
            .update(b"_XMD:BLAKE2b_SSWU_RO_")
            .update(&[(22 + curve_id.len() + domain_prefix.len()) as u8])
            .finalize()
    };

    for (big, buf) in [b_1, b_2].iter().zip(buf.iter_mut()) {
        let mut little = [0u8; CHUNKLEN];
        little.copy_from_slice(big.as_array());
        little.reverse();
        *buf = F::from_uniform_bytes(&little);
    }
}
