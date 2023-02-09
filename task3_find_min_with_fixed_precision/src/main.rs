//! Find min of Rosenbrock's function with fixed precision.

use std::{cmp::Ordering, iter::Sum, ops::Div};

use nalgebra::Vector2;


// TODO(refactor): rename to `Vec2d`
type Vec2 = Vector2<f64>;


const PRECISION: f64 = 1e-3;
const SOLUTION: (f64, f64) = (1., 1.);


fn f(p: Vec2) -> f64 {
    (1.-p.x).powi(2) + 100.*(p.y-p.x.powi(2)).powi(2)
}


fn main() {
    let point_start = Vec2::new(-1.7, 1.7);

    println!("solution by coordinate descent:");
    println!("{}", find_min_by_coordinate_descent(point_start));
    // answer: x = 0.9995093629951177 , y = 0.9990616715937082

    println!("solution by fastest descent:");
    println!("{}", find_min_by_fastest_descent(point_start));
    // answer: x = 1.0004988735118054 , y = 1.0009999995954986

    println!("solution by downhill simplex:");
    println!("{}", find_min_by_downhill_simplex(point_start));
    // answer: x = 1.0004704501887318 , y = 1.000914952298626
}


fn find_min_by_coordinate_descent(point_start: Vec2) -> Vec2 {
    const MAX_ITERATION: usize = 20;
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

    let solution_exact = Vec2::new(SOLUTION.0, SOLUTION.1);

    let mut point = point_start;
    while !is_precise_enough(point, solution_exact) {
        point = find_min_along_x(point);
        if is_precise_enough(point, solution_exact) { break; }
        point = find_min_along_y(point);
    }
    point
}


fn find_min_by_fastest_descent(point_start: Vec2) -> Vec2 {
    fn find_min_along_gradient(point: Vec2) -> Vec2 {
        const MAX_ITERATION: usize = 10;
        const DERIVATIVE_TOLERANCE: f64 = 1e-3;
        const DERIVATIVE_DELTA: f64 = 1e-4;
        const STEP_SCALE: f64 = 2e-3;

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

    let solution_exact = Vec2::new(SOLUTION.0, SOLUTION.1);

    let mut point = point_start;
    while !is_precise_enough(point, solution_exact) {
        point = find_min_along_gradient(point);
    }
    point
}


fn find_min_by_downhill_simplex(point_start: Vec2) -> Vec2 {
    const INITIAL_SIMPLEX_SCALE: f64 = 2.0;
    const LERP_T: f64 = 0.5;

    let solution_exact = Vec2::new(SOLUTION.0, SOLUTION.1);

    let points: [Vec2; 3] = [
        point_start - INITIAL_SIMPLEX_SCALE*Vec2::identity_along_x(),
        point_start + INITIAL_SIMPLEX_SCALE*Vec2::identity_along_y(),
        point_start + INITIAL_SIMPLEX_SCALE*Vec2::new(1., -1.),
    ];
    let mut points: [(Vec2, f64); 3] = points.map(|p| (p, f(p)));
    trait GetPoints<T, const N: usize> {
        fn get_points(&self) -> [T; N];
    }
    impl GetPoints<Vec2, 3> for [(Vec2, f64); 3] {
        fn get_points(&self) -> [Vec2; 3] {
            self.map(|pv| pv.0)
        }
    }
    trait GetValues<T, const N: usize> {
        fn get_values(&self) -> [T; N];
    }
    impl GetValues<f64, 3> for [(Vec2, f64); 3] {
        fn get_values(&self) -> [f64; 3] {
            self.map(|pv| pv.1)
        }
    }
    while !is_precise_enough(points.get_points().avg(), solution_exact) {
        let index_of_max = points.get_values().index_of_max().unwrap();
        let (point_max, value_at_max) = points[index_of_max];

        let other_points_indices = match index_of_max {
            0 => [1, 2],
            1 => [0, 2],
            2 => [0, 1],
            _ => unreachable!()
        };
        let point_a = points[other_points_indices[0]].0;
        let point_b = points[other_points_indices[1]].0;

        let point_symmetric = point_max.mirror_relative_to(point_a, point_b);
        let value_at_symmetric = f(point_symmetric);

        points[index_of_max] = match value_at_symmetric.partial_cmp(&value_at_max).unwrap() {
            Ordering::Less => {
                (point_symmetric, value_at_symmetric)
            }
            Ordering::Greater | Ordering::Equal => {
                let point_ab = (point_a + point_b) / 2.;
                let point_lerp = point_max.lerp(&point_ab, LERP_T);
                let value_at_lerp = f(point_lerp);
                (point_lerp, value_at_lerp)
            }
        };
    }
    points.get_points().avg()
}


fn is_precise_enough(solution_found: Vec2, solution_exact: Vec2) -> bool {
    (solution_found.x - solution_exact.x).abs() < PRECISION
    &&
    (solution_found.y - solution_exact.y).abs() < PRECISION
}





trait Vec2Exts {
    fn zero() -> Self;
    fn identity_along_x() -> Self;
    fn identity_along_y() -> Self;
}
impl Vec2Exts for Vec2 {
    /// Returns `Vec2 { x: 0, y: 0 }`.
    fn zero()             -> Self { Vec2::new(0., 0.) }
    /// Returns `Vec2 { x: 1, y: 0 }`.
    fn identity_along_x() -> Self { Vec2::new(1., 0.) }
    /// Returns `Vec2 { x: 0, y: 1 }`.
    fn identity_along_y() -> Self { Vec2::new(0., 1.) }
}

trait Avg<T> {
    /// Calculates average.
    fn avg(self) -> T;
}
impl<T, const N: usize> Avg<T> for [T; N]
where T: Div<f64, Output = T> + Sum<T> + Sized
{
    /// Calculates average of elements for fixed-size array.
    fn avg(self) -> T {
        let self_len = self.len();
        self.into_iter().sum::<T>() / self_len as f64
    }
}

trait MirrorRelativeTo {
    fn mirror_relative_to(&self, a: Vec2, b: Vec2) -> Vec2;
}
impl MirrorRelativeTo for Vec2 {
    fn mirror_relative_to(&self, a: Vec2, b: Vec2) -> Vec2 {
        let ab = (a + b) / 2.;
        self + 2.*(ab - self)
    }
}

trait IndexOfMax {
    /// Returns index of element with maximum value.
    fn index_of_max(&self) -> Option<usize>;
}
impl<T: PartialOrd, const N: usize> IndexOfMax for [T; N] {
    /// Returns index of element with maximum value for fixed-size array.
    fn index_of_max(&self) -> Option<usize> {
        match self.len() {
            0 => None,
            1 => Some(0),
            _ => {
                let mut index_of_max = 0;
                for i in 1..self.len() {
                    if self[i] > self[index_of_max] {
                        index_of_max = i;
                    }
                }
                Some(index_of_max)
            }
        }
    }
}





#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mirror_relative_to() {
        assert_eq!(Vec2::new(-1., 0.), Vec2::new(1., 0.).mirror_relative_to(Vec2::new(0., 1.), Vec2::new(0., -1.)));
        assert_eq!(Vec2::new(-1., 0.), Vec2::new(1., 0.).mirror_relative_to(Vec2::new(0., -1.), Vec2::new(0., 1.)));
        assert_eq!(Vec2::new(1., 0.), Vec2::new(-1., 0.).mirror_relative_to(Vec2::new(0., 1.), Vec2::new(0., -1.)));
        assert_eq!(Vec2::new(1., 0.), Vec2::new(-1., 0.).mirror_relative_to(Vec2::new(0., -1.), Vec2::new(0., 1.)));

        assert_eq!(Vec2::new(0., -1.), Vec2::new(0., 1.).mirror_relative_to(Vec2::new(1., 0.), Vec2::new(-1., 0.)));
        assert_eq!(Vec2::new(0., -1.), Vec2::new(0., 1.).mirror_relative_to(Vec2::new(-1., 0.), Vec2::new(1., 0.)));
        assert_eq!(Vec2::new(0., 1.), Vec2::new(0., -1.).mirror_relative_to(Vec2::new(1., 0.), Vec2::new(-1., 0.)));
        assert_eq!(Vec2::new(0., 1.), Vec2::new(0., -1.).mirror_relative_to(Vec2::new(-1., 0.), Vec2::new(1., 0.)));

        assert_eq!(Vec2::new(1., 1.), Vec2::new(0., 0.).mirror_relative_to(Vec2::new(1., 0.), Vec2::new(0., 1.)));
        assert_eq!(Vec2::new(2., 2.), Vec2::new(0., 0.).mirror_relative_to(Vec2::new(2., 0.), Vec2::new(0., 2.)));
    }
}

