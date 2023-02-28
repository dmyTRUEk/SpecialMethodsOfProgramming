//! Extensions.

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


pub trait ExtGenFromArray<T, const N: usize> {
    fn gen_from_array(&mut self, array: [T; N]) -> T;
}
impl<T: Copy, const N: usize> ExtGenFromArray<T, N> for ThreadRng {
    fn gen_from_array(&mut self, array: [T; N]) -> T {
        let index = self.gen_range(0 .. array.len());
        array[index]
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

