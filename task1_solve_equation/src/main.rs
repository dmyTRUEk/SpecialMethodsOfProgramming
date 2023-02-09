//! Solution of arbitrary equation.

use std::cmp::Ordering;


fn f(x: f64) -> f64 {
    x.sin() - 1./x
}


pub fn main() {
    println!("solutions:");

    let (solution, f_evals) = find_solution_by_dichotomy(0.1, 1.5, &mut 0);
    println!("by dichotomy   : x = {x}\tf_evals: {fe}", x=solution, fe=f_evals);
    // answer: 1.114156723022461  , f_evals = 36

    let (solution, f_evals) = find_solution_by_chords(0.1, 1.5);
    println!("by chords      : x = {x}\tf_evals: {fe}", x=solution, fe=f_evals);
    // answer: 1.1141571408717854 , f_evals = 24

    let (solution, f_evals) = find_solution_by_newton(1.);
    println!("by newton      : x = {x}\tf_evals: {fe}", x=solution, fe=f_evals);
    // answer: 1.114157140871924  , f_evals = 12

    let (solution, f_evals) = find_solution_by_direct_iterations(1.);
    println!("by direct iters: x = {x}\tf_evals: {fe}", x=solution, fe=f_evals);
    // answer: 1.1141570186671306 , f_evals = 10
}


const TOLERANCE: f64 = 1e-6;


fn find_solution_by_dichotomy(l: f64, r: f64, f_evals: &mut u64) -> (f64, u64) {
    let m = (l + r) / 2.;
    if f(m).abs() < TOLERANCE { return (m, *f_evals); }
    *f_evals += 2;
    match f(m).partial_cmp(&0.) {
        Some(Ordering::Greater) => { find_solution_by_dichotomy(l, m, f_evals) }
        Some(Ordering::Less)    => { find_solution_by_dichotomy(m, r, f_evals) }
        _ => { panic!() }
    }
}

fn find_solution_by_chords(x0: f64, x1: f64) -> (f64, u64) {
    const MAX_ITERS: u32 = 30;
    let mut f_evals = 0;
    let mut x_n_m2 = x0; // X_(n-2)
    let mut x_n_m1 = x1; // X_(n-1)
    let mut x_n    = 0.; // X_n
    for _ in 0..MAX_ITERS {
        x_n = x_n_m1 - f(x_n_m1) * (x_n_m1 - x_n_m2) / (f(x_n_m1) - f(x_n_m2));
        f_evals += 3;
        if (x_n - x_n_m1).abs() < TOLERANCE { break }
        x_n_m2 = x_n_m1;
        x_n_m1 = x_n;
    }
    (x_n, f_evals)
}

fn find_solution_by_newton(x: f64) -> (f64, u64) {
    const MAX_ITERS: u32 = 30;
    let mut f_evals = 0;
    let mut x_n_m1 = x;  // X_(n-1)
    let mut x_n    = 0.; // X_n
    fn d(x: f64) -> f64 {
        const DELTA: f64 = 1e-3;
        (f(x+DELTA) - f(x-DELTA)) / (2.*DELTA)
    }
    for _ in 0..MAX_ITERS {
        x_n = x_n_m1 - f(x_n_m1) / d(x_n_m1);
        f_evals += 3;
        if (x_n - x_n_m1).abs() < TOLERANCE { break }
        x_n_m1 = x_n;
    }
    (x_n, f_evals)
}

fn find_solution_by_direct_iterations(x: f64) -> (f64, u64) {
    const MAX_ITERS: u32 = 1000;
    let mut f_evals = 0;
    let mut x_prev = f64::NAN;
    let mut x = x;
    for _ in 0..MAX_ITERS {
        x = -f(x) + x;
        f_evals += 1;
        if (x - x_prev).abs() < TOLERANCE { break }
        x_prev = x;
    }
    (x, f_evals)
}

