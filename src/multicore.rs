pub use maybe_rayon::{
    iter::{IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
    join, scope, Scope,
};

#[cfg(feature = "multicore")]
pub use maybe_rayon::{
    current_num_threads,
    iter::{IndexedParallelIterator, IntoParallelRefIterator},
    slice::ParallelSliceMut,
};

#[cfg(not(feature = "multicore"))]
pub fn current_num_threads() -> usize {
    1
}
