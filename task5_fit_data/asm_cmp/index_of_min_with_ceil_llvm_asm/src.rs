// default:
pub trait IndexOfMinWithFloor<T> {
    fn index_of_min_with_ceil(&self, ceil: T) -> Option<usize>;
}
impl IndexOfMinWithFloor<f64> for Vec<f64> {
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



// combmatch:
pub trait IndexOfMinWithFloor<T> {
    fn index_of_min_with_ceil(&self, ceil: T) -> Option<usize>;
}
impl IndexOfMinWithFloor<f64> for Vec<f64> {
    fn index_of_min_with_ceil(&self, ceil: f64) -> Option<usize> {
        let mut option_index_of_min = None;
        for i in 0..self.len() {
            if self[i] >= ceil || !self[i].is_finite() { continue }
            match option_index_of_min {
                None
                | Some(..) if self[i] < self[option_index_of_min.unwrap()] => {
                    option_index_of_min = Some(i);
                }
                _ => {}
            }
        }
        option_index_of_min
    }
}

