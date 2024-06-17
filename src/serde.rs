use std::io::{self, Read, Write};

#[derive(Clone, Copy, Debug)]
pub struct Repr<const T: usize>([u8; T]);

impl<const T: usize> From<[u8; T]> for Repr<T> {
    fn from(bytes: [u8; T]) -> Self {
        Self(bytes)
    }
}

impl<const T: usize> From<Repr<T>> for [u8; T] {
    fn from(repr: Repr<T>) -> Self {
        repr.0
    }
}

impl<const T: usize> Default for Repr<T> {
    fn default() -> Self {
        Self([0u8; T])
    }
}

impl<const T: usize> AsMut<[u8]> for Repr<T> {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl<const T: usize> AsRef<[u8]> for Repr<T> {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl<const T: usize> std::ops::Index<std::ops::Range<usize>> for Repr<T> {
    type Output = [u8];

    fn index(&self, range: std::ops::Range<usize>) -> &Self::Output {
        &self.0[range]
    }
}

impl<const T: usize> std::ops::Index<usize> for Repr<T> {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<const T: usize> std::ops::Index<std::ops::RangeTo<usize>> for Repr<T> {
    type Output = [u8];

    fn index(&self, range: std::ops::RangeTo<usize>) -> &Self::Output {
        &self.0[range]
    }
}

impl<const T: usize> std::ops::Index<std::ops::RangeFrom<usize>> for Repr<T> {
    type Output = [u8];

    fn index(&self, range: std::ops::RangeFrom<usize>) -> &Self::Output {
        &self.0[range]
    }
}

impl<const T: usize> std::ops::IndexMut<std::ops::Range<usize>> for Repr<T> {
    fn index_mut(&mut self, range: std::ops::Range<usize>) -> &mut Self::Output {
        &mut self.0[range]
    }
}

/// Trait for converting raw bytes to/from the internal representation of a type.
/// For example, field elements are represented in Montgomery form and serialized/deserialized without Montgomery reduction.
pub trait SerdeObject: Sized {
    /// The purpose of unchecked functions is to read the internal memory representation
    /// of a type from bytes as quickly as possible. No sanitization checks are performed
    /// to ensure the bytes represent a valid object. As such this function should only be
    /// used internally as an extension of machine memory. It should not be used to deserialize
    /// externally provided data.
    fn from_raw_bytes_unchecked(bytes: &[u8]) -> Self;
    fn from_raw_bytes(bytes: &[u8]) -> Option<Self>;

    fn to_raw_bytes(&self) -> Vec<u8>;

    /// The purpose of unchecked functions is to read the internal memory representation
    /// of a type from disk as quickly as possible. No sanitization checks are performed
    /// to ensure the bytes represent a valid object. This function should only be used
    /// internally when some machine state cannot be kept in memory (e.g., between runs)
    /// and needs to be reloaded as quickly as possible.
    fn read_raw_unchecked<R: Read>(reader: &mut R) -> Self;
    fn read_raw<R: Read>(reader: &mut R) -> io::Result<Self>;

    fn write_raw<W: Write>(&self, writer: &mut W) -> io::Result<()>;
}

pub(crate) mod endian {

    pub trait Endian {
        fn to_bytes(res: &mut [u8], el: &[u64]);
        fn from_bytes(res: &[u8], el: &mut [u64]);
    }

    pub struct LE;
    pub struct BE;

    impl Endian for LE {
        fn to_bytes(res: &mut [u8], el: &[u64]) {
            el.iter().enumerate().for_each(|(i, limb)| {
                let off = i * 8;
                res[off..off + 8].copy_from_slice(&limb.to_le_bytes());
            });
        }

        fn from_bytes(res: &[u8], el: &mut [u64]) {
            el.iter_mut().enumerate().for_each(|(i, limb)| {
                let off = i * 8;
                *limb = u64::from_le_bytes(res[off..off + 8].try_into().unwrap());
            });
        }
    }
    impl Endian for BE {
        fn to_bytes(res: &mut [u8], el: &[u64]) {
            el.iter().rev().enumerate().for_each(|(i, limb)| {
                let off = i * 8;
                res[off..off + 8].copy_from_slice(&limb.to_be_bytes());
            });
        }

        fn from_bytes(res: &[u8], el: &mut [u64]) {
            el.iter_mut().rev().enumerate().for_each(|(i, limb)| {
                let off = i * 8;
                *limb = u64::from_be_bytes(res[off..off + 8].try_into().unwrap());
            });
        }
    }
}
