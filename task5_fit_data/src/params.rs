//! Params type.

use std::{
    iter::Sum,
    ops::{Add, Div, Mul, Sub},
};

use crate::{
    PARAMS_DIFF_TYPE,
    fit::DiffFunctionType,
    float_type::float,
    function::Function,
    param::{Param, ParamName, ParamValue},
};


#[derive(Debug, Clone, PartialEq)]
pub struct Params {
    // TODO(optimize)?: use HashMap.
    params: Vec<Param>
}

#[allow(dead_code)]
impl Params {
    pub const fn empty() -> Self {
        Self { params: vec![] }
    }

    pub const fn new(params: Vec<Param>) -> Self {
        Self { params }
    }

    pub fn one(name: ParamName, value: ParamValue) -> Self {
        Self::new(vec![Param::new(name, value)])
    }

    pub fn from_params_array<const N: usize>(params: [Param; N]) -> Self {
        Self::new(params.to_vec())
    }

    pub const fn from_params_vec(params: Vec<Param>) -> Self {
        Self::new(params)
    }

    pub fn from_array<const N: usize>(array: [(char, float); N]) -> Self {
        Self::new(
            array
                .map(|(name, value)| Param::new(name, value))
                .to_vec()
        )
    }

    pub fn from_vec(vec: Vec<(char, float)>) -> Self {
        Self::new(
            vec.into_iter()
                .map(|(name, value)| Param::new(name, value))
                .collect()
        )
    }

    pub fn from_names_and_values_arrays<const N: usize>(names: [ParamName; N], values: [ParamValue; N]) -> Self {
        Self::from_params_vec(
            names.iter().zip(values)
                .map(|(&name, value)| Param::new(name, value))
                .collect()
        )
    }

    pub fn from_names_and_values_vecs(names: Vec<ParamName>, values: Vec<ParamValue>) -> Self {
        // assert_eq!(names.len(), values.len());
        Self::from_params_vec(
            names.iter().zip(values)
                .map(|(&name, value)| Param::new(name, value))
                .collect()
        )
    }

    pub fn gen_random_from_f(f: &Function) -> Self {
        let mut params_names = f.get_params_names();
        params_names.sort();
        params_names.dedup();
        let params = params_names.into_iter()
            .map(|name| Param::gen_random_value_with_name(name))
            .collect();
        Self::new(params)
    }

    pub fn amount(&self) -> usize {
        self.params.len()
    }

    pub fn get_all(&self) -> Vec<Param> {
        self.params.clone()
    }

    pub fn get_all_names(&self) -> Vec<ParamName> {
        self.params.iter()
            .map(|p| p.name)
            .collect()
    }

    pub fn get_all_values(&self) -> Vec<ParamValue> {
        self.params.iter()
            .map(|p| p.value)
            .collect()
    }

    pub fn get_by_index(&self, index: usize) -> Param {
        self.params[index].clone()
    }

    pub fn get_by_name_unchecked(&self, name: ParamName) -> float {
        self.params.iter()
            .find(|param| param.name == name)
            // .expect(&format!("parameter `{name}` not found in params: {params:?}"))
            .unwrap()
            .value
    }

    pub fn get_by_name_checked(&self, name: ParamName) -> Option<float> {
        self.params.iter()
           .find(|param| param.name == name)
           .map(|p| p.value)
    }

    pub fn insert(&mut self, param: Param) {
        self.params.push(param)
    }

    pub fn set_by_name(&mut self, name: ParamName, new_value: float) {
        for i in 0..self.params.len() {
            if self.params[i].name == name {
                self.params[i].value = new_value;
                return;
            }
        }
    }

    pub fn change_param_by(&mut self, name: ParamName, delta: float) {
        for i in 0..self.params.len() {
            if self.params[i].name == name {
                self.params[i].value += delta;
                return;
            }
        }
    }
    pub fn changed_param_by(&self, name: ParamName, delta: float) -> Self {
        let mut self_ = self.clone();
        for i in 0..self_.params.len() {
            if self_.params[i].name == name {
                self_.params[i].value += delta;
            }
        }
        self_
    }

    pub fn change_all_params_by(&mut self, delta: float) {
        for i in 0..self.params.len() {
            self.params[i].value += delta;
        }
    }
    pub fn changed_all_params_by(&self, delta: float) -> Self {
        let mut self_ = self.clone();
        for i in 0..self_.params.len() {
            self_.params[i].value += delta;
        }
        self_
    }

    pub fn scale_param_by(&mut self, name: ParamName, delta: float) {
        for i in 0..self.params.len() {
            if self.params[i].name == name {
                self.params[i].value *= delta;
                return;
            }
        }
    }
    pub fn scaled_param_by(&self, name: ParamName, delta: float) -> Self {
        let mut self_ = self.clone();
        for i in 0..self_.params.len() {
            if self_.params[i].name == name {
                self_.params[i].value *= delta;
            }
        }
        self_
    }

    pub fn scale_all_params_by(&mut self, delta: float) {
        for i in 0..self.params.len() {
            self.params[i].value *= delta;
        }
    }
    pub fn scaled_all_params_by(&self, delta: float) -> Self {
        let mut self_ = self.clone();
        for i in 0..self_.params.len() {
            self_.params[i].value *= delta;
        }
        self_
    }

    pub fn diff(&self, other: Self) -> float {
        // assert_eq!(self.get_all_names(), other.get_all_names());
        self.params.iter().zip(other.params)
            .map(|(param_a, param_b)| match PARAMS_DIFF_TYPE {
                DiffFunctionType::DyAbs     => (param_a.value - param_b.value).abs(),
            //  DiffFunctionType::DySquared => (param_a.value - param_b.value).powi(2),
                DiffFunctionType::DySquared => unimplemented!("bc you must also take sqrt at the end"),
                DiffFunctionType::LeastDist => unimplemented!(),
            })
            .sum()
    }
}

impl Add for Params {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        // assert_eq!(self.get_all_names(), rhs.get_all_names());
        let values: Vec<ParamValue> = self.params.iter().zip(rhs.params)
            .map(|(param_a, param_b)| param_a.value + param_b.value)
            .collect();
        Self::from_names_and_values_vecs(self.get_all_names(), values)
    }
}

impl Sub for Params {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        // assert_eq!(self.get_all_names(), rhs.get_all_names());
        let values: Vec<ParamValue> = self.params.iter().zip(rhs.params)
            .map(|(param_a, param_b)| param_a.value - param_b.value)
            .collect();
        Self::from_names_and_values_vecs(self.get_all_names(), values)
    }
}

impl Mul<float> for Params {
    type Output = Params;
    fn mul(self, rhs: float) -> Self::Output {
        self.scaled_all_params_by(rhs)
    }
}
impl Mul<Params> for float {
    type Output = Params;
    fn mul(self, rhs: Params) -> Self::Output {
        rhs * self
    }
}

impl Div<float> for Params {
    type Output = Params;
    fn div(self, rhs: float) -> Self::Output {
        self.scaled_all_params_by(1./rhs)
    }
}

impl Sum for Params {
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        let first = iter.next();
        if first.is_none() { return Params::empty() }
        let mut result: Params = first.unwrap();
        for param in iter {
            result = result + param;
        }
        result
    }
}

