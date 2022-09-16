mod curve;
mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;

#[cfg(all(feature = "asm", target_arch = "x86_64"))]
mod assembly;

pub use curve::*;
pub use engine::*;
pub use fq::*;
use fq2::*;
pub use fr::*;

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}
