//! Find min of Rosenbrock's function with fixed precision.

use nalgebra::Vector2;


type Vec2 = Vector2<f64>;


const PRECISION: f64 = 10e-3;
const SOLUTION: (f64, f64) = (1., 1.);


fn f(p: Vec2) -> f64 {
    (1.-p.x).powi(2) + 100.*(p.y-p.x.powi(2)).powi(2)
}


fn main() {
    println!("solution by coordinate descent: {}", find_min_by_coordinate_descent(Vec2::zeros()));
    // answer: x = 0.9997922598100459 , y = 0.9995361936770234
}


fn find_min_by_coordinate_descent(point_start: Vec2) -> Vec2 {
    const MAX_ITERATION: usize = 30;
    const TOLERANCE: f64 = 10e-5;
    const DELTA_FOR_DERIVATIVE: f64 = 10e-3;

    fn find_min_along_x(p: Vec2) -> Vec2 {
        fn d(p: Vec2) -> f64 {
            let mut delta = Vec2::zeros();
            delta.x = DELTA_FOR_DERIVATIVE;
            (f(p+delta) - f(p-delta)) / (2.*DELTA_FOR_DERIVATIVE)
        }
        let mut p_n_m1 = p; // P_(n-1)
        let mut p_n    = p; // P_n
        for _ in 0..MAX_ITERATION {
            p_n.x = p_n_m1.x - f(p_n_m1) / d(p_n_m1);
            if (p_n.x - p_n_m1.x).abs() < TOLERANCE { return p_n; }
            p_n_m1 = p_n;
        }
        p_n
    }

    fn find_min_along_y(p: Vec2) -> Vec2 {
        fn d(p: Vec2) -> f64 {
            let mut delta = Vec2::zeros();
            delta.y = DELTA_FOR_DERIVATIVE;
            (f(p+delta) - f(p-delta)) / (2.*DELTA_FOR_DERIVATIVE)
        }
        let mut p_n_m1 = p; // P_(n-1)
        let mut p_n    = p; // P_n
        for _ in 0..MAX_ITERATION {
            p_n.y = p_n_m1.y - f(p_n_m1) / d(p_n_m1);
            if (p_n.y - p_n_m1.y).abs() < TOLERANCE { return p_n; }
            p_n_m1 = p_n;
        }
        p_n
    }

    let solution = Vec2::new(SOLUTION.0, SOLUTION.1);

    let mut point = point_start;
    while (point - solution).norm_squared() > PRECISION.powi(2) {
        point = find_min_along_x(point);
        point = find_min_along_y(point);
    }
    point
}

