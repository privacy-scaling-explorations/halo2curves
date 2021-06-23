mod ec;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;
// mod fr;

pub use fq::*;

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}
