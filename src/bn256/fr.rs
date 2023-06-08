#[cfg(feature = "maybe_u64")]
pub use super::fr_internal::Fr as FrInternal;
#[cfg(feature = "maybe_u64")]
pub type Fr = crate::MaybeU64<FrInternal>;

#[cfg(not(feature = "maybe_u64"))]
pub use super::fr_internal::Fr;

#[cfg(feature = "maybe_u64")]
impl crate::MaybeU64<FrInternal> {
    pub const fn from_raw(data: [u64; 4]) -> Self {
        if data[1] == 0 && data[2] == 0 && data[3] == 0 {
            Self::U64(data[0])
        } else {
            Self::Full(FrInternal::from_raw(data))
        }
    }

    pub const fn size() -> usize {
        FrInternal::size()
    }
}
