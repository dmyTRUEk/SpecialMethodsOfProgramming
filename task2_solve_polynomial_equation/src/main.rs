//! Solve polynomial equation.

#![feature(generic_const_exprs)]

use num::{complex::Complex64, Zero};



fn main() {
    let p = Polynomial::new([7., 2., 7., 2., 1.]); // 7 + 2x + 7x² + 2x³ + x⁴
    println!("{}", p.to_string());
    println!("{}", p.eval_at(3.));
    todo!();
}



type LinearPolynomial = Polynomial<2>;

#[derive(Debug, Clone, PartialEq)]
struct Polynomial<const N: usize> {
    k: [Complex64; N], // coefficients: a0 + a1*x + a2*x² + …
}

impl<const N: usize> Polynomial<N>
where [(); N-1]:
{
    pub fn new(k: [impl Into<Complex64> + Copy; N]) -> Self {
        Self {
            k: k
                .iter()
                .map(|&x| x.into())
                .collect::<Vec<_>>()
                .try_into().unwrap()
        }
    }

    #[allow(non_snake_case)]
    pub const fn get_N(&self) -> usize { N }

    pub fn find_all_solutions(&self) -> [Complex64; N] {
        let mut p: Polynomial<N> = self.clone();
        let mut solution: Complex64 = p.find_one_solution(0.);
        let mut solutions = vec![];
        for _ in 0..N {
            solutions.push(solution);
            let p = p.div_by(LinearPolynomial::new([solution, (1.).into()])).unwrap();
            solution = p.find_one_solution(0.);
            // solution = match N {
            //     2 => {
            //         Polynomial::<2>::find_one_solution(&p, 0.)
            //     }
            //     _ => panic!()
            // }
        }
        solutions.try_into().unwrap()
    }

    pub fn find_one_solution(&self, start_point: impl Into<Complex64>) -> Complex64 {
        todo!()
    }

    pub fn eval_at(&self, x: impl Into<Complex64>) -> Complex64 {
        let x: Complex64 = x.into();
        let mut r = Complex64::zero();
        for (i, k) in self.k.iter().enumerate() {
            r += k * x.powi(i as i32)
        }
        r
    }

    pub fn div_by(&self, rhs: LinearPolynomial) -> Option<Polynomial<{N-1}>> {
        // let mut remainder: Complex64 = 0.;
        todo!()
    }
}

impl<const N: usize> ToString for Polynomial<N> {
    fn to_string(&self) -> String {
        let mut parts = Vec::<String>::with_capacity(N);
        for (i, k) in self.k.iter().enumerate() {
            parts.push(
                match i {
                    0 => format!("{k}"),
                    1 => format!("{k} x"),
                    i if 2 <= i && i <= 9 => format!("{k} x{}", ['²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'][i-2]),
                    _ => format!("{k} x^{i}"),
                }
            );
        }
        parts.into_iter().reduce(
            |acc, el| {
                acc.to_string() + &(if el.chars().next().unwrap() != '-' {
                    format!(" + {el}")
                } else {
                    format!(" - {elm}", elm=&el[1..])
                })
            }
        ).unwrap()
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_eval_at() {
        let p = Polynomial::new([5.2, 6.1, 9.0, 0.42]);
        assert_eq!(Complex64::new(5.2, 0.), p.eval_at(0.));
        assert_eq!(Complex64::new(34.2479225, 0.), Polynomial::new([5.2, 6.1, 9.0, 0.42]).eval_at(1.45));
    }

    #[test]
    fn polynomial_div() {
        assert_eq!(
            Some(Polynomial::new([1.])),
            Polynomial::new([1., 2.]).div_by(Polynomial::new([1., 2.]))
        );

        assert_eq!(
            Some(Polynomial::new([1., 1.])),
            Polynomial::new([1., 2., 1.]).div_by(Polynomial::new([1., 1.]))
        );
        assert_eq!(
            None,
            Polynomial::new([1., 2., 1.]).div_by(Polynomial::new([1., -1.]))
        );

        assert_eq!(
            Some(Polynomial::new([1., 1.])),
            Polynomial::new([1., -2., 1.]).div_by(Polynomial::new([-1., 1.]))
        );
        assert_eq!(
            Some(Polynomial::new([-1., 1.])),
            Polynomial::new([1., -2., 1.]).div_by(Polynomial::new([1., 1.]))
        );
    }
}

