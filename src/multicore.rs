#[cfg(all(
    feature = "multicore",
    target_arch = "wasm32",
    not(target_feature = "atomics")
))]
compile_error!(
    "The multicore feature flag is not supported on wasm32 architectures without atomics"
);

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
