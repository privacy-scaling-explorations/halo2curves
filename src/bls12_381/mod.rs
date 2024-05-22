use crate::arithmetic::mul_512;
use crate::arithmetic::sbb;
use crate::{
    arithmetic::{CurveEndo, EndoParameters},
    endo,
};
pub use bls12_381::{Fp, G1Projective, Scalar};
use ff::PrimeField;
use ff::WithSmallOrderMulGroup;
use std::convert::TryInto;
use std::io::{self, Read, Write};


pub use bls12_381::*;


#[cfg(feature = "derive_serde")]
use serde::{Deserialize, Serialize};


// Obtained from https://github.com/ConsenSys/gnark-crypto/blob/master/ecc/utils.go
// See https://github.com/demining/Endomorphism-Secp256k1/blob/main/README.md
// to have more details about the endomorphism.
const ENDO_PARAMS_BLS: EndoParameters = EndoParameters {
    // round(b2/n)
    gamma2: [0x63f6e522f6cfee30u64, 0x7c6becf1e01faadd, 0x01, 0x0],
    // round(-b1/n)
    gamma1: [0x02u64, 0x0, 0x0, 0x0],
    b1: [0x01u64, 0x0, 0x0, 0x0],
    b2: [0x0000000100000000, 0xac45a4010001a402, 0x0, 0x0],
};

endo!(G1Projective, Scalar, ENDO_PARAMS_BLS);


impl crate::serde::SerdeObject for G1Affine {
    /// The purpose of unchecked functions is to read the internal memory representation
    /// of a type from bytes as quickly as possible. No sanitization checks are performed
    /// to ensure the bytes represent a valid object. As such this function should only be
    /// used internally as an extension of machine memory. It should not be used to deserialize
    /// externally provided data.
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        G1Affine::from_compressed(bytes.try_into().unwrap()).unwrap()
    }
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        Some(G1Affine::from_compressed(bytes.try_into().unwrap()).unwrap())
    }

    fn to_raw_bytes(&self) -> Vec<u8> {
        self.to_compressed().into()
    }

    /// The purpose of unchecked functions is to read the internal memory representation
    /// of a type from disk as quickly as possible. No sanitization checks are performed
    /// to ensure the bytes represent a valid object. This function should only be used
    /// internally when some machine state cannot be kept in memory (e.g., between runs)
    /// and needs to be reloaded as quickly as possible.
    fn read_raw_unchecked<R: Read>(reader: &mut R) -> Self {
        let mut buf = [0; 48];
        reader.read_exact(&mut buf).unwrap();
        G1Affine::from_compressed(&buf).unwrap()

    }
    fn read_raw<R: Read>(reader: &mut R) -> io::Result<Self> {
        let mut buf = [0; 48];
        reader.read_exact(&mut buf).unwrap();
        Ok(G1Affine::from_compressed(&buf).unwrap())
    }

    fn write_raw<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        writer.write_all(&self.to_compressed())?;
        Ok(())

    }
}


#[test]
fn test_endo() {
    use ff::Field;
    use rand_core::OsRng;

    for _ in 0..100000 {
        let k = Scalar::random(OsRng);
        let (k1, k1_neg, k2, k2_neg) = G1Projective::decompose_scalar(&k);
        if k1_neg & k2_neg {
            assert_eq!(
                k,
                -Scalar::from_u128(k1) + Scalar::ZETA * Scalar::from_u128(k2)
            )
        } else if k1_neg {
            assert_eq!(
                k,
                -Scalar::from_u128(k1) - Scalar::ZETA * Scalar::from_u128(k2)
            )
        } else if k2_neg {
            assert_eq!(
                k,
                Scalar::from_u128(k1) + Scalar::ZETA * Scalar::from_u128(k2)
            )
        } else {
            assert_eq!(
                k,
                Scalar::from_u128(k1) - Scalar::ZETA * Scalar::from_u128(k2)
            )
        }
    }

    for _ in 0..100000 {
        let k = Scalar::random(OsRng);
        let (k1, k1_neg, k2, k2_neg) = G1Projective::decompose_scalar(&k);
        if k1_neg & k2_neg {
            assert_eq!(
                k,
                -Scalar::from_u128(k1) + Scalar::ZETA * Scalar::from_u128(k2)
            )
        } else if k1_neg {
            assert_eq!(
                k,
                -Scalar::from_u128(k1) - Scalar::ZETA * Scalar::from_u128(k2)
            )
        } else if k2_neg {
            assert_eq!(
                k,
                Scalar::from_u128(k1) + Scalar::ZETA * Scalar::from_u128(k2)
            )
        } else {
            assert_eq!(
                k,
                Scalar::from_u128(k1) - Scalar::ZETA * Scalar::from_u128(k2)
            )
        }
    }
}
