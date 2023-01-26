//! Solution of arbitrary equation in compile time.

#![feature(
    const_cmp,
    const_eval_limit,
    const_fn_floating_point_arithmetic,
    const_trait_impl,
)]

#![const_eval_limit = "1000000"]

use std::cmp::Ordering;



const fn f(x: f64) -> f64 {
    sin(x) - 1./x
}


pub fn main() {
    const X: f64 = find_solution_by_dichotomy(0.1, 1.5);
    println!("X = {X}");
    // answer: 1.114156723022461
}


const fn find_solution_by_dichotomy(l: f64, r: f64) -> f64 {
    let m = (l + r) / 2.;
    const TOLERANCE: f64 = 1e-6;
    if abs(f(m)) < TOLERANCE { return m; }
    match f(m).partial_cmp(&0.) {
        Some(Ordering::Less)    => { find_solution_by_dichotomy(m, r) }
        Some(Ordering::Greater) => { find_solution_by_dichotomy(l, m) }
        _ => { panic!() }
    }
}





const PI: f64 = 3.1415_926_535_897;


const fn abs(x: f64) -> f64 {
    if x > 0. { x } else { -x }
}

const fn powi(x: f64, n: i8) -> f64 {
    let mut r = 1.;
    let mut i = 0;
    while i < n {
        r *= x;
        i += 1;
    }
    r
}

const fn fact(n: u8) -> u64 {
    let mut r = 1;
    let mut i = 2;
    while i <= n {
        r *= i as u64;
        i += 1;
    }
    r
}

const fn sin(x: f64) -> f64 {
    let x = x % (2.*PI); // TODO: % PI/2
    let mut r = x;
    let mut i = 3;
    let mut sign = -1;
    while i < 21 {
        r += powi(x, i) / (sign * fact(i as u8) as i64) as f64;
        sign = -sign;
        i += 2;
    }
    r
}






#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn powi_() {
        assert_eq!(1., powi(10., 0));
        assert_eq!(10., powi(10., 1));
        assert_eq!(100., powi(10., 2));
        assert_eq!(1_000., powi(10., 3));
        assert_eq!(2048., powi(2., 11))
    }

    #[test]
    fn fact_() {
        assert_eq!(120, fact(5));
        assert_eq!(720, fact(6));
    }

    #[test]
    fn sin_() {
        assert_eq!(0., sin(0.));
        assert_eq!(1., sin(PI/2.));
    }
}

