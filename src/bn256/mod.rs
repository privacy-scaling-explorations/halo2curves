mod engine;
mod fq;
mod fq12;
mod fq2;
mod fq6;
mod fr;
mod g;

pub use fq::*;
pub use fq2::*;
pub use fr::*;
pub use g::*;
pub use engine::*;

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}

// //     let mut g1 = G1::generator();
// //     let g2 = G2::generator();
// //     g1 = g1.double();
// //     let pair21 = Bn256::pairing(&G1Affine::from(g1), &G2Affine::from(g2));

// //     assert_eq!(pair12, pair21);

// //     let g1 = G1::generator();
// //     let mut g2 = G2::generator();
// //     g2 = g2.double().double();
// //     let pair12 = Bn256::pairing(&G1Affine::from(g1), &G2Affine::from(g2));

// //     let mut g1 = G1::generator();
// //     let mut g2 = G2::generator();
// //     g1 = g1.double();
// //     g2 = g2.double();
// //     let pair21 = Bn256::pairing(&G1Affine::from(g1), &G2Affine::from(g2));

// //     assert_eq!(pair12, pair21);

// //     let mut rng = XorShiftRng::from_seed([
// //         0x59, 0x62, 0xbe, 0x5d, 0x76, 0x3d, 0x31, 0x8d, 0x17, 0xdb, 0x37, 0x32, 0x54, 0x06, 0xbc,
// //         0xe5,
// //     ]);
// //     for _ in 0..1000 {
// //         let a = Fr::random(&mut rng);
// //         let b = Fr::random(&mut rng);

// //         let mut g1 = G1::generator();
// //         g1.mul_assign(a);

// //         let mut g2 = G2::generator();
// //         g1.mul_assign(b);

// //         let pair_ab = Bn256::pairing(&G1Affine::from(g1), &G2Affine::from(g2));

// //         g1 = G1::generator();
// //         g1.mul_assign(b);

// //         g2 = G2::generator();
// //         g1.mul_assign(a);

// //         let pair_ba = Bn256::pairing(&G1Affine::from(g1), &G2Affine::from(g2));

// //         assert_eq!(pair_ab, pair_ba);
// //     }
// // }
