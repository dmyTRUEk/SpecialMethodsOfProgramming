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
    sin(x) - 1.0/x
}


pub fn main() {
    const X: f64 = find_solution_by_dichotomy(0.1, 1.5);
    println!("X = {X}");
    // answer: 1.114156723022461
}


const fn find_solution_by_dichotomy(l: f64, r: f64) -> f64 {
    let m: f64 = (l + r) / 2.0;
    const TOLERANCE: f64 = 1e-6;
    if abs(f(m)) < TOLERANCE { return m; }
    match f(m).partial_cmp(&0.0) {
        Some(Ordering::Less)    => { find_solution_by_dichotomy(m, r) }
        Some(Ordering::Greater) => { find_solution_by_dichotomy(l, m) }
        _ => { panic!() }
    }
}





const PI: f64 = 3.1415_926_535_897;


const fn abs(x: f64) -> f64 {
    if x > 0.0 { x } else { -x }
}

const fn powi(x: f64, n: i8) -> f64 {
    let mut r: f64 = 1.0;
    let mut i = 0;
    while i < n {
        r *= x;
        i += 1;
    }
    r
}

const fn fact(n: u8) -> u64 {
    let mut r: u64 = 1;
    let mut i = 2;
    while i <= n {
        r *= i as u64;
        i += 1;
    }
    r
}

const fn sin(x: f64) -> f64 {
    let x: f64 = x % (2.0*PI); // TODO: % PI/2
    let mut r: f64 = x;
    let mut i = 3;
    let mut sign: i64 = -1;
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
        assert_eq!(1.0, powi(10.0, 0));
        assert_eq!(10.0, powi(10.0, 1));
        assert_eq!(100.0, powi(10.0, 2));
        assert_eq!(1_000.0, powi(10.0, 3));
        assert_eq!(2048.0, powi(2.0, 11))
    }

    #[test]
    fn fact_() {
        assert_eq!(120, fact(5));
        assert_eq!(720, fact(6));
    }

    #[test]
    fn sin_() {
        assert_eq!(0.0, sin(0.0));
        assert_eq!(1.0, sin(PI/2.0));
    }
}

