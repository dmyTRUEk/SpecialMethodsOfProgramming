//! A lot of useful extensions.

pub trait ArrayMath {
    fn add(&self, rhs: Self) -> Self;
    fn sub(&self, rhs: Self) -> Self;
    fn scale(&self, rhs: f64) -> Self;
    fn unscale(&self, rhs: f64) -> Self;
}
impl ArrayMath for Vec<f64> {
    fn add(&self, rhs: Self) -> Self {
        assert_eq!(self.len(), rhs.len());
        self.iter().zip(rhs).map(|(l, r)| l + r).collect()
    }
    fn sub(&self, rhs: Self) -> Self {
        assert_eq!(self.len(), rhs.len());
        self.iter().zip(rhs).map(|(l, r)| l - r).collect()
    }
    fn scale(&self, rhs: f64) -> Self {
        self.iter().map(|l| l * rhs).collect()
    }
    fn unscale(&self, rhs: f64) -> Self {
        self.iter().map(|l| l / rhs).collect()
    }
}


pub trait Avg<T> {
    /// Calculates average.
    fn avg(self) -> T;
}
impl Avg<Vec<f64>> for Vec<Vec<f64>> {
    /// Calculates average of elements for dynamic-size array.
    fn avg(self) -> Vec<f64> {
        let len: f64 = self.len() as f64;
        let sum: Vec<f64> = self.into_iter().reduce(|acc, el| acc.add(el)).unwrap();
        sum.unscale(len)
    }
}


pub trait IndexOfMaxWithCeil<T> {
    fn index_of_max_with_ceil(&self, ceil : T) -> Option<usize>;
}
impl IndexOfMaxWithCeil<f64> for Vec<f64> {
    fn index_of_max_with_ceil(&self, ceil: f64) -> Option<usize> {
        let mut option_index_of_max = None;
        for i in 0..self.len() {
            if self[i] >= ceil || !self[i].is_finite() { continue }
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
pub trait IndexOfMaxWithFloor<T> {
    fn index_of_max_with_floor(&self, floor: T) -> Option<usize>;
}
impl IndexOfMaxWithFloor<f64> for Vec<f64> {
    fn index_of_max_with_floor(&self, floor: f64) -> Option<usize> {
        let mut option_index_of_max = None;
        for i in 0..self.len() {
            if self[i] <= floor || !self[i].is_finite() { continue }
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
pub trait IndexOfMinWithCeil<T> {
    fn index_of_min_with_ceil(&self, ceil : T) -> Option<usize>;
}
impl IndexOfMinWithCeil<f64> for Vec<f64> {
    fn index_of_min_with_ceil(&self, ceil: f64) -> Option<usize> {
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
pub trait IndexOfMinWithFloor<T> {
    fn index_of_min_with_floor(&self, floor: T) -> Option<usize>;
}
impl IndexOfMinWithFloor<f64> for Vec<f64> {
    fn index_of_min_with_floor(&self, floor: f64) -> Option<usize> {
        let mut option_index_of_min = None;
        for i in 0..self.len() {
            if self[i] <= floor || !self[i].is_finite() { continue }
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
impl IndexOfMax<f64> for Vec<f64> {
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
pub trait IndexOfMin<T> {
    fn index_of_min(&self) -> Option<usize>;
}
impl IndexOfMin<f64> for Vec<f64> {
    fn index_of_min(&self) -> Option<usize> {
        let mut option_index_of_min = None;
        for i in 0..self.len() {
            if !self[i].is_finite() { continue }
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


// pub trait SeparateChunksFromStart {
//     fn separate_chunks_from_start(&self, delimiter: impl ToString, chunks_size: usize) -> String;
// }
pub trait SeparateChunksFromEnd {
    fn separate_chunks_from_end(&self, delimiter: impl ToString, chunks_size: usize) -> String;
}
impl SeparateChunksFromEnd for String {
    fn separate_chunks_from_end(&self, delimiter: impl ToString, chunks_size: usize) -> String {
        let len = self.len();
        self.chars()
            .enumerate()
            .map(|(i, c)| if (len-i) % chunks_size != 0 || i == 0 { c.to_string() } else { format!("{}{}", delimiter.to_string(), c) })
            .collect()
    }
}
impl SeparateChunksFromEnd for &str {
    fn separate_chunks_from_end(&self, delimiter: impl ToString, chunks_size: usize) -> String {
        self.to_string().separate_chunks_from_end(delimiter, chunks_size)
    }
}

pub trait ToStringUnderscoreSeparated {
    fn to_string_underscore_separated(&self) -> String;
}

impl ToStringUnderscoreSeparated for u64 {
    fn to_string_underscore_separated(&self) -> String {
        self.to_string().separate_chunks_from_end("_", 3)
    }
}


pub trait ToStringWithSignificantDigits {
    fn to_string_with_significant_digits(&self, significant_digits: usize) -> String;
}
impl ToStringWithSignificantDigits for f64 {
    fn to_string_with_significant_digits(&self, precision: usize) -> String {
        let a = self.abs();
        let precision = if a > 1. {
            let n = (1. + a.log10().floor()) as usize;
            if n <= precision {
                precision - n
            } else {
                0
            }
        } else if a > 0. {
            let n = -(1. + a.log10().floor()) as usize;
            precision + n
        } else {
            0
        };
        format!("{0:.1$}", self, precision)
    }
}



#[cfg(test)]
mod index_of {
    mod max {
        mod with_ceil {
            mod without_nan {
                use super::super::super::super::IndexOfMaxWithCeil;
                #[test]
                fn ceil_between_min_and_max() {
                    assert_eq!(Some(0), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_max_with_ceil(42.));
                }
                #[test]
                fn ceil_below_min() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_max_with_ceil(-100.));
                }
                #[test]
                fn ceil_above_max() {
                    assert_eq!(Some(7), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_max_with_ceil(1000.));
                }
            }
            mod with_nan {
                use super::super::super::super::IndexOfMaxWithCeil;
                #[test]
                fn ceil_between_min_and_max() {
                    assert_eq!(Some(0), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_max_with_ceil(42.));
                }
                #[test]
                fn ceil_below_min() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_max_with_ceil(-100.));
                }
                #[test]
                fn ceil_above_max() {
                    assert_eq!(Some(7), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_max_with_ceil(1000.));
                }
            }
        }
        mod with_floor {
            mod without_nan {
                use super::super::super::super::IndexOfMaxWithFloor;
                #[test]
                fn floor_between_min_and_max() {
                    assert_eq!(Some(7), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_max_with_floor(42.));
                }
                #[test]
                fn floor_below_min() {
                    assert_eq!(Some(7), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_max_with_floor(-100.));
                }
                #[test]
                fn floor_above_max() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_max_with_floor(1000.));
                }
            }
            mod with_nan {
                use super::super::super::super::IndexOfMaxWithFloor;
                #[test]
                fn floor_between_min_and_max() {
                    assert_eq!(Some(7), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_max_with_floor(42.));
                }
                #[test]
                fn floor_below_min() {
                    assert_eq!(Some(7), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_max_with_floor(-100.));
                }
                #[test]
                fn floor_above_max() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_max_with_floor(1000.));
                }
            }
        }
    }
    mod min {
        mod with_ceil {
            mod without_nan {
                use super::super::super::super::IndexOfMinWithCeil;
                #[test]
                fn ceil_between_min_and_max() {
                    assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_min_with_ceil(42.));
                }
                #[test]
                fn ceil_below_min() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_min_with_ceil(-100.));
                }
                #[test]
                fn ceil_above_max() {
                    assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_min_with_ceil(1000.));
                }
            }
            mod with_nan {
                use super::super::super::super::IndexOfMinWithCeil;
                #[test]
                fn ceil_between_min_and_max() {
                    assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_min_with_ceil(42.));
                }
                #[test]
                fn ceil_below_min() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_min_with_ceil(-100.));
                }
                #[test]
                fn ceil_above_max() {
                    assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_min_with_ceil(1000.));
                }
            }
        }
        mod with_floor {
            mod without_nan {
                use super::super::super::super::IndexOfMinWithFloor;
                #[test]
                fn floor_between_min_and_max() {
                    assert_eq!(Some(6), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_min_with_floor(42.));
                }
                #[test]
                fn floor_below_min() {
                    assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_min_with_floor(-100.));
                }
                #[test]
                fn floor_above_max() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520.].index_of_min_with_floor(1000.));
                }
            }
            mod with_nan {
                use super::super::super::super::IndexOfMinWithFloor;
                #[test]
                fn floor_between_min_and_max() {
                    assert_eq!(Some(6), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_min_with_floor(42.));
                }
                #[test]
                fn floor_below_min() {
                    assert_eq!(Some(5), vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_min_with_floor(-100.));
                }
                #[test]
                fn floor_above_max() {
                    assert_eq!(None, vec![14., 0., 1., 4., 8., -53., 43., 520., f64::NAN].index_of_min_with_floor(1000.));
                }
            }
        }
    }
}

#[cfg(test)]
mod separate_chunks_from_end {
    use super::SeparateChunksFromEnd;
    #[test]
    fn a() {
        assert_eq!("a", "a".separate_chunks_from_end("_-", 3));
    }
    #[test]
    fn ab() {
        assert_eq!("ab", "ab".separate_chunks_from_end("_-", 3));
    }
    #[test]
    fn abc() {
        assert_eq!("abc", "abc".separate_chunks_from_end("_-", 3));
    }
    #[test]
    fn a_bcd() {
        assert_eq!("a_-bcd", "abcd".separate_chunks_from_end("_-", 3));
    }
    #[test]
    fn abcdefghijklmnopqrstuvwxyz() {
        assert_eq!("ab_-cde_-fgh_-ijk_-lmn_-opq_-rst_-uvw_-xyz", "abcdefghijklmnopqrstuvwxyz".separate_chunks_from_end("_-", 3));
    }
}

#[cfg(test)]
mod to_string_with_significant_digits {
    use super::ToStringWithSignificantDigits;
    mod _1234_5678 {
        use super::*;
        const X: f64 = 1234.5678;
        // anser isn't 1234 but 1235, bc of rounding
        // #[ignore]
        #[test]
        fn _0() {
            assert_eq!("1235", X.to_string_with_significant_digits(0));
        }
        #[test]
        fn _1() {
            assert_eq!("1235", X.to_string_with_significant_digits(1));
        }
        #[test]
        fn _2() {
            assert_eq!("1235", X.to_string_with_significant_digits(2));
        }
        #[test]
        fn _3() {
            assert_eq!("1235", X.to_string_with_significant_digits(3));
        }
        #[test]
        fn _4() {
            assert_eq!("1235", X.to_string_with_significant_digits(4));
        }
        #[test]
        fn _5() {
            assert_eq!("1234.6", X.to_string_with_significant_digits(5));
        }
        #[test]
        fn _6() {
            assert_eq!("1234.57", X.to_string_with_significant_digits(6));
        }
        #[test]
        fn _7() {
            assert_eq!("1234.568", X.to_string_with_significant_digits(7));
        }
        #[test]
        fn _8() {
            assert_eq!("1234.5678", X.to_string_with_significant_digits(8));
        }
        #[test]
        fn _9() {
            assert_eq!("1234.56780", X.to_string_with_significant_digits(9));
        }
        #[test]
        fn _10() {
            assert_eq!("1234.567800", X.to_string_with_significant_digits(10));
        }
    }
    mod _000_1234 {
        use super::*;
        const X: f64 = 0.001234;
        // #[ignore]
        #[test]
        fn _0() {
            assert_eq!("0.00", X.to_string_with_significant_digits(0));
        }
        #[test]
        fn _1() {
            assert_eq!("0.001", X.to_string_with_significant_digits(1));
        }
        #[test]
        fn _2() {
            assert_eq!("0.0012", X.to_string_with_significant_digits(2));
        }
        #[test]
        fn _3() {
            assert_eq!("0.00123", X.to_string_with_significant_digits(3));
        }
        #[test]
        fn _4() {
            assert_eq!("0.001234", X.to_string_with_significant_digits(4));
        }
        #[test]
        fn _5() {
            assert_eq!("0.0012340", X.to_string_with_significant_digits(5));
        }
        #[test]
        fn _6() {
            assert_eq!("0.00123400", X.to_string_with_significant_digits(6));
        }
        #[test]
        fn _7() {
            assert_eq!("0.001234000", X.to_string_with_significant_digits(7));
        }
    }
}

