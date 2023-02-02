//! Find min of function.

use std::{cmp::Ordering, iter::Sum, ops::Div};

use nalgebra::Vector2;


// TODO(refactor): rename to `Vec2d`
type Vec2 = Vector2<f64>;


const PRECISION: f64 = 1e-3;


fn f(p: Vec2) -> f64 {
    let (x, y) = (p.x, p.y);
    ( 1. + (x+y+1.).powi(2) * (19.-14.*x+3.*x.powi(2)-14.*y+6.*x*y+3.*y.powi(2)) )
    *
    ( 30. + (2.*x-3.*y).powi(2) * (18.-32.*x+12.*x.powi(2)+48.*y-36.*x*y+27.*y.powi(2)) )
}


fn main() {
    println!("solution by downhill simplex:");
    println!("{}", find_min_by_downhill_simplex(Vec2::zero()));
    // answer: x = -0.0004520396141742822 , y = -1.0000038171601773
}


fn find_min_by_downhill_simplex(point_start: Vec2) -> Vec2 {
    const INITIAL_SIMPLEX_SCALE: f64 = 2.0;
    const LERP_T: f64 = 0.5;

    let mut points_prev: [Vec2; 3] = [Vec2::nan(); 3];
    let mut points: [Vec2; 3] = [
        point_start - INITIAL_SIMPLEX_SCALE*Vec2::identity_along_x(),
        point_start + INITIAL_SIMPLEX_SCALE*Vec2::identity_along_y(),
        point_start + INITIAL_SIMPLEX_SCALE*Vec2::new(1., -1.),
    ];
    while !is_close_enough(points.avg(), points_prev.avg()) {
        let (index_of_max, value_at_max) = points.index_and_value_of_max_by_key(f).unwrap();
        let point_max = points[index_of_max];
        let other_points_indices = match index_of_max {
            0 => [1, 2],
            1 => [0, 2],
            2 => [0, 1],
            _ => unreachable!()
        };

        let point_a = points[other_points_indices[0]];
        let point_b = points[other_points_indices[1]];

        let point_symmetric = point_max.mirror_relative_to(point_a, point_b);
        let value_at_point_symmetric = f(point_symmetric);

        points_prev = points;
        points[index_of_max] = match value_at_point_symmetric.partial_cmp(&value_at_max).unwrap() {
            Ordering::Less => {
                point_symmetric
            }
            Ordering::Greater | Ordering::Equal => {
                let point_ab = (point_a + point_b) / 2.;
                point_max.lerp(&point_ab, LERP_T)
            }
        };
    }
    points.avg()
}


fn is_close_enough(a: Vec2, b: Vec2) -> bool {
    (a.x - b.x).abs() < PRECISION/2.
    &&
    (a.y - b.y).abs() < PRECISION/2.
}





trait Vec2Exts {
    fn zero() -> Self;
    fn identity_along_x() -> Self;
    fn identity_along_y() -> Self;
    fn nan() -> Self;
}
impl Vec2Exts for Vec2 {
    /// Returns `Vec2 { x: 0, y: 0 }`.
    fn zero()             -> Self { Vec2::new(0., 0.) }
    /// Returns `Vec2 { x: 1, y: 0 }`.
    fn identity_along_x() -> Self { Vec2::new(1., 0.) }
    /// Returns `Vec2 { x: 0, y: 1 }`.
    fn identity_along_y() -> Self { Vec2::new(0., 1.) }
    /// Returns `Vec2 { x: NaN, y: NaN }`.
    fn nan() -> Self { Vec2::new(f64::NAN, f64::NAN) }
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

trait IndexAndValueOfMaxByKey<T, R> {
    /// Returns index of element with maximum value mapped by `f` and value of `f` at it.
    fn index_and_value_of_max_by_key(self, f: fn(T) -> R) -> Option<(usize, R)>;
}
impl<T: Clone, R: PartialOrd + Clone, const N: usize> IndexAndValueOfMaxByKey<T, R> for [T; N] {
    /// Returns index of element with maximum value mapped by `f` and value of `f` at it for fixed-size array.
    fn index_and_value_of_max_by_key(self, f: fn(T) -> R) -> Option<(usize, R)> {
        match self.len() {
            0 => None,
            1 => Some((0, f(self[0].clone()))),
            _ => {
                let array_mapped = self.map(f);
                let index_of_max = array_mapped.index_of_max().unwrap();
                Some((index_of_max, array_mapped[index_of_max].clone()))
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

