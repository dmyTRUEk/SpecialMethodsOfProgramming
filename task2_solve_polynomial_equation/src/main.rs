//! Solve polynomial equation.

use num::{
    complex::{Complex64, ComplexFloat},
    One,
    Zero,
};
use rand::{thread_rng, Rng};



fn main() {
    let p = Polynomial::new([7, 2, 7, 2, 1]); // 7 + 2x + 7x² + 2x³ + x⁴
    println!("Solving: {}", p.to_string());
    println!("find_all_solutions = {:#?}", p.find_all_solutions());
    // answers:
    // 0.035644726065908516 + 1.0816682476181705i
    // 0.035644726065908516 - 1.0816682476181705i
    // -1.035644722361672 + 2.2144580376118843i
    // -1.035644722361672 - 2.2144580376118843i
}



#[derive(Debug, Clone, PartialEq)]
struct Polynomial {
    k: Vec<Complex64>, // coefficients: a0 + a1*x + a2*x² + …
}

impl Polynomial {
    pub fn new(k: impl MyInto<Vec<Complex64>>) -> Self {
        Self { k: k.into_() }
    }

    /// returns: `C0 + C1 x`
    pub fn binomial(k: impl MyInto<[Complex64; 2]>) -> Self {
        let k: [Complex64; 2] = k.into_();
        Self { k: vec![k[0], k[1]] }
    }

    /// returns: `C + 1x`
    pub fn binomial_normalized(c: impl MyInto<Complex64>) -> Self {
        Self::binomial([c.into_(), Complex64::one()])
    }

    /// returns: `C`
    pub fn monomial(c: impl Into<Complex64>) -> Self {
        Self { k: vec![c.into()] }
    }

    pub fn one() -> Self { Polynomial::monomial(Complex64::one()) }

    pub fn get_size(&self) -> usize { self.k.len() }

    pub fn get_coefs(&self) -> Vec<Complex64> { self.k.clone() }
    pub fn get_k_at(&self, i: usize) -> Complex64 { self.k[i] }
    pub fn get_k_last(&self) -> Complex64 { self.k.last().unwrap().clone() }

    pub fn is_binomial_with_norm_k(&self) -> bool {
        self.get_size() == 2 && self.k[1] == 1.0.into()
    }

    const TOLERANCE: f64 = 1e-6;

    pub fn find_all_solutions(&self) -> Vec<Complex64> {
        let mut p: Polynomial = self.clone();
        let mut solution: Complex64;
        let mut solutions = vec![];
        while p.get_size() > 1 {
            solution = p.find_one_solution(42.);
            solutions.push(solution);
            p = p.div_by(Polynomial::binomial_normalized(-solution)).unwrap();
            if solution.im.abs() > Self::TOLERANCE {
                solutions.push(solution.conj());
                p = p.div_by(Polynomial::binomial_normalized(-solution.conj())).unwrap();
            }
            else {
                let sl = solutions.len();
                solutions[sl-1].im = 0.;
            }
        }
        solutions
    }

    pub fn find_one_solution(&self, start_point: impl Into<Complex64>) -> Complex64 {
        // by Newton's method:
        let mut x_prev = Complex64::new(f64::MAX, f64::MAX);
        let start_point: Complex64 = start_point.into();
        let mut x: Complex64 = start_point;
        for _i in 0..100 {
            x = x - self.eval_at(x) / self.eval_derivative_at(x);
            if (x - x_prev).norm_sqr() < Self::TOLERANCE { return x; }
            if x.is_nan() { return self.find_one_solution(start_point); }
            x_prev = x;
        }
        x
    }

    pub fn eval_at(&self, x: impl Into<Complex64>) -> Complex64 {
        let x: Complex64 = x.into();
        let mut r = Complex64::zero();
        for (i, k) in self.k.iter().enumerate() {
            r += k * x.powi(i as i32)
        }
        r
    }

    pub fn eval_derivative_at(&self, x: impl Into<Complex64>) -> Complex64 {
        let mut rng = thread_rng(); // TODO(optimization): don't create every time => move out of this functions => create it higher
        let dx = Complex64::new(rng.gen_range(-1. .. 1.), rng.gen_range(-1. .. 1.));
        let dx = 1e-4 * dx / dx.norm();
        let x: Complex64 = x.into();
        (self.eval_at(x+dx) - self.eval_at(x-dx)) / (2. * dx)
    }

    /// Example:
    ///
    ///  1x² + 2x + 1  |  x + 1
    ///  --------------+-------
    ///  1x² + 2x
    /// -
    ///  1x² + 1x      <- 1*(x+1) => 1x
    ///  --------                    |
    ///        1x + 1                |
    ///       -                      |
    ///        1x + 1  <- 1*(x+1) => 1
    ///        ------                |
    ///             0                V
    ///                              1x+1
    ///
    pub fn div_by(&self, rhs: Polynomial) -> Option<Polynomial> {
        if self.k == rhs.k { return Some(Polynomial::one()); }
        assert!(rhs.is_binomial_with_norm_k(), "else is unimplemented");
        let mut p = self.clone();
        let mut coefs = Vec::<Complex64>::with_capacity(self.get_size()-1);
        while p.get_size() > 1 {
            let k = p.get_k_last(); // `/ 1.0` is omitted, bc rhs is binomial with k1==1
            coefs.push(k);
            let pk = p.get_coefs();
            let pk_len = pk.len();
            let mut pk_new = pk.into_iter().take(pk_len-1).collect::<Vec<_>>();
            let pk_new_len = pk_new.len();
            pk_new[pk_new_len-1] -= k*rhs.get_k_at(0);
            p = Polynomial::new(pk_new);
        }
        let remainder: f64 = p.get_k_at(0).norm_sqr();
        if remainder > Self::TOLERANCE { return None; }
        coefs.reverse();
        assert_eq!(self.get_size()-1, coefs.len());
        Some(Polynomial::new(coefs))
    }
}

impl ToString for Polynomial {
    fn to_string(&self) -> String {
        let mut parts = Vec::<String>::with_capacity(self.get_size());
        for (i, k) in self.k.iter().enumerate() {
            let k = if k.im == 0. { k.re.to_string() } else { format!("({})", k.to_string()) };
            parts.push(
                match i {
                    0 => format!("{k}"),
                    1 => format!("{k} x"),
                    i if 2 <= i && i <= 9 => format!("{k} x{}", ['²','³','⁴','⁵','⁶','⁷','⁸','⁹'][i-2]),
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



trait MyInto<T> {
    fn into_(self) -> T;
}

impl MyInto<Complex64> for i64 {
    fn into_(self) -> Complex64 {
        (self as f64).into()
    }
}
impl MyInto<Complex64> for f64 {
    fn into_(self) -> Complex64 {
        self.into()
    }
}
impl MyInto<Complex64> for Complex64 {
    fn into_(self) -> Complex64 {
        self
    }
}

impl<const N: usize> MyInto<Vec<Complex64>> for [i64; N] {
    fn into_(self) -> Vec<Complex64> {
        self.into_iter().map(|x| (x as f64).into()).collect()
    }
}
impl<const N: usize> MyInto<Vec<Complex64>> for [f64; N] {
    fn into_(self) -> Vec<Complex64> {
        self.into_iter().map(|x| x.into()).collect()
    }
}
// impl<const N: usize> MyInto<Vec<Complex64>> for [Complex64; N] {
//     fn into_(self) -> Vec<Complex64> {
//         self.into_iter().map(|x| x.into()).collect()
//     }
// }

// impl MyInto<Vec<Complex64>> for Vec<i64> {
//     fn into_(self) -> Vec<Complex64> {
//         self.into_iter().map(|x| (x as f64).into()).collect()
//     }
// }
// impl MyInto<Vec<Complex64>> for Vec<f64> {
//     fn into_(self) -> Vec<Complex64> {
//         self.into_iter().map(|x| x.into()).collect()
//     }
// }
impl MyInto<Vec<Complex64>> for Vec<Complex64> {
    fn into_(self) -> Vec<Complex64> {
        self
    }
}

// impl<const N: usize> MyInto<[Complex64; N]> for [i64; N] {
//     fn into_(self) -> [Complex64; N] {
//         self.into_iter().map(|x| (x as f64).into()).collect::<Vec<_>>().try_into().unwrap()
//     }
// }
// impl<const N: usize> MyInto<[Complex64; N]> for [f64; N] {
//     fn into_(self) -> [Complex64; N] {
//         self.into_iter().map(|x| (x as f64).into()).collect::<Vec<_>>().try_into().unwrap()
//     }
// }
impl<const N: usize> MyInto<[Complex64; N]> for [Complex64; N] {
    fn into_(self) -> [Complex64; N] {
        self
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
            Some(Polynomial::one()),
            Polynomial::binomial_normalized(2.).div_by(Polynomial::binomial_normalized(2.))
        );

        assert_eq!(
            Some(Polynomial::binomial_normalized(3.)),
            Polynomial::new([21, 10, 1]).div_by(Polynomial::binomial_normalized(7.))
        );

        assert_eq!(
            Some(Polynomial::binomial_normalized(1.)),
            Polynomial::new([1, 2, 1]).div_by(Polynomial::binomial_normalized(1.))
        );
        assert_eq!(
            None,
            Polynomial::new([1, 2, 1]).div_by(Polynomial::binomial_normalized(-1.))
        );

        assert_eq!(
            Some(Polynomial::binomial_normalized(-1.)),
            Polynomial::new([1, -2, 1]).div_by(Polynomial::binomial_normalized(-1.))
        );
    }

    #[ignore]
    #[test]
    fn polynomial_find_all_solutions() {
        assert_eq!(
            [
                -1,
            ].into_(),
            Polynomial::binomial_normalized(1).find_all_solutions()
        );
        assert_eq!(
            [
                1,
            ].into_(),
            Polynomial::binomial_normalized(-1).find_all_solutions()
        );
        assert_eq!(
            [
                2,
                3,
                5,
                7,
            ].into_(),
            Polynomial::new([210, -247, 101, -17, 1]).find_all_solutions()
        );
    }
}

