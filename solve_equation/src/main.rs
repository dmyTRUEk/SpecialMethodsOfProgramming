//! Solution of arbitrary equation.

use std::cmp::Ordering;


fn f(x: f64) -> f64 {
    x.sin() - 1.0/x
}


pub fn main() {
    println!("solution_by_dichotomy: {x}", x=find_solution_by_dichotomy(0.1, 1.5));
    // answer: 1.114156723022461
}


const TOLERANCE: f64 = 1e-6;


fn find_solution_by_dichotomy(l: f64, r: f64) -> f64 {
    let m: f64 = (l + r) / 2.0;
    if f(m).abs() < TOLERANCE { return m; }
    match f(m).partial_cmp(&0.0) {
        Some(Ordering::Less)    => { find_solution_by_dichotomy(m, r) }
        Some(Ordering::Greater) => { find_solution_by_dichotomy(l, m) }
        _ => { panic!() }
    }
}

