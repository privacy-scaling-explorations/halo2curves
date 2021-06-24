mod ec;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;

pub use ec::*;
pub use fq::*;
pub use fr::*;

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}
