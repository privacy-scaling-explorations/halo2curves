pub use maybe_rayon::{
    iter::{IntoParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
    join, scope, Scope,
};

#[cfg(all(feature = "multicore", target_pointer_width = "64"))]
pub use maybe_rayon::{
    current_num_threads,
    iter::{IndexedParallelIterator, IntoParallelRefIterator},
    slice::ParallelSliceMut,
};

#[cfg(any(not(feature = "multicore"), target_pointer_width = "32"))]
pub fn current_num_threads() -> usize {
    1
}
