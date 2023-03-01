//! Extensions.

use std::{ops::Div, iter::Sum};

use rand::{Rng, rngs::ThreadRng};

use crate::float_type::float;


pub trait IndexOfMinWithFloor<T> {
    fn index_of_min_with_ceil(&self, ceil: T) -> Option<usize>;
}
impl IndexOfMinWithFloor<float> for Vec<float> {
    fn index_of_min_with_ceil(&self, ceil: float) -> Option<usize> {
        let mut option_index_of_min = None;
        for i in 0..self.len() {
            if self[i] >= ceil || !self[i].is_finite() { continue }
            match option_index_of_min {
                None => {
                    option_index_of_min = Some(i);
                }
                Some(index_of_min) if self[i] < self[index_of_min] => {
                    option_index_of_min = Some(i);
                }
                _ => {}
            }
        }
        option_index_of_min
    }
}

pub trait IndexOfMax<T> {
    fn index_of_max(&self) -> Option<usize>;
}
impl IndexOfMax<float> for Vec<float> {
    fn index_of_max(&self) -> Option<usize> {
        let mut option_index_of_max = None;
        for i in 0..self.len() {
            if !self[i].is_finite() { continue }
            match option_index_of_max {
                None => {
                    option_index_of_max = Some(i);
                }
                Some(index_of_max) if self[i] > self[index_of_max] => {
                    option_index_of_max = Some(i);
                }
                _ => {}
            }
        }
        option_index_of_max
    }
}


pub trait ExtGenFromArray<T, const N: usize> {
    fn gen_from_array(&mut self, array: [T; N]) -> T;
}
impl<T: Copy, const N: usize> ExtGenFromArray<T, N> for ThreadRng {
    fn gen_from_array(&mut self, array: [T; N]) -> T {
        let index = self.gen_range(0 .. array.len());
        array[index]
    }
}


pub trait Avg<T> {
    /// Calculates average.
    fn avg(self) -> T;
}
impl<T, const N: usize> Avg<T> for [T; N]
where T: Div<f64, Output = T> + Sum<T> + Sized
{
    /// Calculates average of elements for fixed-size array.
    fn avg(self) -> T {
        let len: float = self.len() as float;
        let sum: T = self.into_iter().sum::<T>();
        sum / len
    }
}
impl<T> Avg<T> for Vec<T>
where T: Div<f64, Output = T> + Sum<T> + Sized
{
    /// Calculates average of elements for dynamic-size array.
    fn avg(self) -> T {
        let len: float = self.len() as float;
        let sum: T = self.into_iter().sum::<T>();
        sum / len
    }
}





#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_of_min_with_ceil() {
        assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520., 8494893.].index_of_min_with_ceil(42.));
        assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520., float::NAN].index_of_min_with_ceil(42.));
        assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520., float::NAN].index_of_min_with_ceil(-100.));
        assert_eq!(None, vec![43., 520., 8494893.].index_of_min_with_ceil(42.));
        assert_eq!(None, vec![43., 520., 8494893., 42.].index_of_min_with_ceil(42.));
    }
}

