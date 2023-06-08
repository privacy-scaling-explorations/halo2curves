use ff::PrimeField;
use std::{
    convert::TryInto,
    io::{Read, Result, Write},
};

use crate::{serde::SerdeObject, MaybeU64};

impl<F> SerdeObject for MaybeU64<F>
where
    F: PrimeField<Repr = [u8; 32]> + SerdeObject,
{
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self {
        debug_assert_eq!(bytes.len(), 32);
        if bytes.iter().skip(8).all(|x| *x == 0) {
            Self::U64(u64::from_le_bytes(bytes[0..8].try_into().unwrap()))
        } else {
            Self::Full(F::from_raw_bytes_unchecked(bytes))
        }
    }
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() != 32 {
            return None;
        }
        if bytes.iter().skip(8).all(|x| *x == 0) {
            Some(Self::U64(u64::from_le_bytes(
                bytes[0..8].try_into().unwrap(),
            )))
        } else {
            F::from_raw_bytes(bytes).map(Self::Full)
        }
    }

    fn to_raw_bytes(&self) -> Vec<u8> {
        match self {
            Self::Full(f) => f.to_raw_bytes(),
            Self::U64(f) => [f.to_le_bytes().as_ref(), [0u8; 24].as_ref()].concat(),
        }
    }
    fn read_raw_unchecked<R: Read>(reader: &mut R) -> Self {
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf).unwrap();
        Self::from_raw_bytes_unchecked(buf.as_ref())
    }
    fn read_raw<R: Read>(reader: &mut R) -> Result<Self> {
        let mut buf = [0u8; 32];
        reader.read_exact(&mut buf)?;
        Ok(Self::from_raw_bytes(buf.as_ref()).unwrap())
    }

    fn write_raw<W: Write>(&self, writer: &mut W) -> Result<()> {
        writer.write_all(self.to_raw_bytes().as_ref())
    }
}
