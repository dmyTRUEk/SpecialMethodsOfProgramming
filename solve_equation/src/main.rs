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
    println!("solution_by_newton   : {x}", x=find_solution_by_newton(1.));
    // answer: 1.114157140871924
    println!("solution_by_direct_iterstions: {x}", x=find_solution_by_direct_iterations(1.));
    // answer: -2.772604525254792
}


const TOLERANCE: f64 = 1e-6;


fn find_solution_by_dichotomy(l: f64, r: f64) -> f64 {
    let m = (l + r) / 2.;
    if f(m).abs() < TOLERANCE { return m; }
    match f(m).partial_cmp(&0.) {
        Some(Ordering::Less)    => { find_solution_by_dichotomy(m, r) }
        Some(Ordering::Greater) => { find_solution_by_dichotomy(l, m) }
        _ => { panic!() }
    }
}

fn find_solution_by_chords(x0: f64, x1: f64) -> f64 {
    let mut x_n_m2 = x0; // X_(n-2)
    let mut x_n_m1 = x1; // X_(n-1)
    let mut x_n    = 0.; // X_n
    for _ in 0..30 {
        x_n = x_n_m1 - f(x_n_m1) * (x_n_m1 - x_n_m2) / (f(x_n_m1) - f(x_n_m2));
        if (x_n - x_n_m1).abs() < TOLERANCE { return x_n; }
        x_n_m2 = x_n_m1;
        x_n_m1 = x_n;
    }
    x_n
}

fn find_solution_by_newton(x: f64) -> f64 {
    let mut x_n_m1 = x;  // X_(n-1)
    let mut x_n    = 0.; // X_n
    fn d(x: f64) -> f64 {
        const DELTA: f64 = 1e-3;
        (f(x+DELTA) - f(x-DELTA)) / (2.*DELTA)
    }
    for _ in 0..30 {
        x_n = x_n_m1 - f(x_n_m1) / d(x_n_m1);
        if (x_n - x_n_m1).abs() < TOLERANCE { return x_n; }
        x_n_m1 = x_n;
    }
    x_n
}

fn find_solution_by_direct_iterations(x: f64) -> f64 {
    let mut x_prev = f64::NAN;
    let mut x = x;
    for _ in 0..1000 {
        x = f(x) + x;
        if (x - x_prev).abs() < TOLERANCE { return x; }
        x_prev = x;
    }
    x
}

