use ff::Field;
use std::ops::{Add, Mul};

pub trait LinearCombinationEngine {
    type Lhs: Mul<Self::Rhs, Output = Self::Lhs> + Add;
    type Rhs: Field;

    fn new(base: Self::Rhs) -> Self;
    fn result(&mut self) -> Self::Lhs;
    fn add(&mut self, elem: Self::Lhs);
    fn add_with_aux(&mut self, elem: Self::Lhs, aux: Self::Rhs);
    fn combine(base: Self::Rhs, coeffs: Vec<Self::Lhs>) -> Self::Lhs;
}
