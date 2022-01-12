mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;
mod g;

pub use engine::*;
pub use fq::*;
use fq2::*;
pub use fr::*;
pub use g::*;

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}
