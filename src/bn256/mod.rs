mod curve;
mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;

#[cfg(feature = "asm")]
mod assembly;

pub use curve::*;
pub use engine::*;
pub use fq::*;
pub use fq12::*;
pub use fq2::*;
pub use fq6::*;
pub use fr::*;

#[derive(Debug, PartialEq, Eq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}
