//! Find min of Rosenbrock's function with fixed precision.

use nalgebra::Vector2;


type Vec2 = Vector2<f64>;


const PRECISION: f64 = 1e-3;
const SOLUTION: (f64, f64) = (1., 1.);


fn f(p: Vec2) -> f64 {
    (1.-p.x).powi(2) + 100.*(p.y-p.x.powi(2)).powi(2)
}


fn main() {
    println!("solution by coordinate descent:");
    println!("{}", find_min_by_coordinate_descent(Vec2::zero()));
    // answer: x = 0.9997922598100459 , y = 0.9995361936770234

    println!("solution by fastest descent:");
    println!("{}", find_min_by_fastest_descent(Vec2::zero()));
    // answer: x = 0.9995008705143117 , y = 0.9990000002451329
}


fn find_min_by_coordinate_descent(point_start: Vec2) -> Vec2 {
    const MAX_ITERATION: usize = 30;
    const DERIVATIVE_TOLERANCE: f64 = 1e-4;
    const DERIVATIVE_DELTA: f64 = 1e-2;

    fn find_min_along_x(point: Vec2) -> Vec2 {
        fn derivative_along_x(p: Vec2) -> f64 {
            let mut delta = Vec2::zero();
            delta.x = DERIVATIVE_DELTA;
            (f(p+delta) - f(p-delta)) / (2.*DERIVATIVE_DELTA)
        }
        let mut p_n_m1 = point; // P_(n-1)
        let mut p_n    = point; // P_n
        for _ in 0..MAX_ITERATION {
            p_n.x = p_n_m1.x - f(p_n_m1) / derivative_along_x(p_n_m1);
            if (p_n.x - p_n_m1.x).abs() < DERIVATIVE_TOLERANCE { return p_n; }
            p_n_m1 = p_n;
        }
        p_n
    }

    fn find_min_along_y(point: Vec2) -> Vec2 {
        fn derivative_along_y(p: Vec2) -> f64 {
            let mut delta = Vec2::zero();
            delta.y = DERIVATIVE_DELTA;
            (f(p+delta) - f(p-delta)) / (2.*DERIVATIVE_DELTA)
        }
        let mut p_n_m1 = point; // P_(n-1)
        let mut p_n    = point; // P_n
        for _ in 0..MAX_ITERATION {
            p_n.y = p_n_m1.y - f(p_n_m1) / derivative_along_y(p_n_m1);
            if (p_n.y - p_n_m1.y).abs() < DERIVATIVE_TOLERANCE { return p_n; }
            p_n_m1 = p_n;
        }
        p_n
    }

    let solution = Vec2::new(SOLUTION.0, SOLUTION.1);

    let mut point = point_start;
    while (point.x - solution.x).abs() > PRECISION || (point.y - solution.y).abs() > PRECISION {
        point = find_min_along_x(point);
        point = find_min_along_y(point);
    }
    point
}


fn find_min_by_fastest_descent(point_start: Vec2) -> Vec2 {
    fn find_min_along_gradient(point: Vec2) -> Vec2 {
        const MAX_ITERATION: usize = 10;
        const DERIVATIVE_TOLERANCE: f64 = 1e-3;
        const DERIVATIVE_DELTA: f64 = 1e-4;
        const STEP_SCALE: f64 = 3e-3;

        fn derivative_along_direction(p: Vec2, dir: Vec2) -> f64 {
            let delta = dir.normalize() * DERIVATIVE_DELTA;
            (f(p+delta) - f(p-delta)) / (2.*DERIVATIVE_DELTA)
        }
        fn grad(p: Vec2) -> Vec2 {
            Vec2::new(
                derivative_along_direction(p, Vec2::identity_along_x()),
                derivative_along_direction(p, Vec2::identity_along_y()),
            )
        }
        let grad_dir = grad(point);
        let mut p_n_m1 = point; // P_(n-1)
        let mut p_n    = point; // P_n
        for _ in 0..MAX_ITERATION {
            p_n = p_n_m1 - STEP_SCALE * grad_dir * f(p_n_m1) / derivative_along_direction(p_n_m1, grad_dir);
            if (p_n - p_n_m1).norm_squared() < DERIVATIVE_TOLERANCE.powi(2) { return p_n; }
            p_n_m1 = p_n;
        }
        p_n
    }

    let solution = Vec2::new(SOLUTION.0, SOLUTION.1);

    let mut point = point_start;
    while (point.x - solution.x).abs() > PRECISION || (point.y - solution.y).abs() > PRECISION {
        // println!("{}", point);
        // println!("x = {}\ny = {}\n", point.x, point.y);
        point = find_min_along_gradient(point);
    }
    point
}





trait Vec2Exts {
    fn zero() -> Self;
    fn identity_along_x() -> Self;
    fn identity_along_y() -> Self;
}
impl Vec2Exts for Vec2 {
    fn zero()             -> Self { Vec2::new(0., 0.) }
    fn identity_along_x() -> Self { Vec2::new(1., 0.) }
    fn identity_along_y() -> Self { Vec2::new(0., 1.) }
}

