//! Solution of arbitrary equation.

use std::cmp::Ordering;


fn f(x: f64) -> f64 {
    x.sin() - 1./x
}


pub fn main() {
    println!("solution_by_dichotomy: {x}", x=find_solution_by_dichotomy(0.1, 1.5));
    // answer: 1.114156723022461
    println!("solution_by_chords   : {x}", x=find_solution_by_chords(0.1, 1.5));
    // answer: 1.1141571408717854
}


const TOLERANCE: f64 = 1e-6;


fn find_solution_by_dichotomy(l: f64, r: f64) -> f64 {
    let m: f64 = (l + r) / 2.;
    if f(m).abs() < TOLERANCE { return m; }
    match f(m).partial_cmp(&0.) {
        Some(Ordering::Less)    => { find_solution_by_dichotomy(m, r) }
        Some(Ordering::Greater) => { find_solution_by_dichotomy(l, m) }
        _ => { panic!() }
    }
}

fn find_solution_by_chords(x0: f64, x1: f64) -> f64 {
    let mut x_n_m2 = x0;   // X_(n-2)
    let mut x_n_m1 = x1;   // X_(n-1)
    let mut x_n: f64 = 0.; // X_n
    for _ in 0..30 {
        x_n = x_n_m1 - f(x_n_m1) * (x_n_m1 - x_n_m2) / (f(x_n_m1) - f(x_n_m2));
        if (x_n - x_n_m1).abs() < TOLERANCE { return x_n; }
        x_n_m2 = x_n_m1;
        x_n_m1 = x_n;
    }
    x_n
}

