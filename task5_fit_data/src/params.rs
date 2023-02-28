//! Params type.

use crate::{
    float_type::float,
    function::Function,
    param::{Param, ParamName},
};


#[derive(Debug, Clone, PartialEq)]
pub struct Params {
    // TODO(optimize)?: use HashMap.
    params: Vec<Param>
}

#[allow(dead_code)]
impl Params {
    pub const fn new() -> Self {
        Self { params: vec![] }
    }

    pub fn from_params_array<const N: usize>(params: [Param; N]) -> Self {
        Self { params: params.to_vec() }
    }

    pub const fn from_params_vec(params: Vec<Param>) -> Self {
        Self { params }
    }

    pub fn from_array<const N: usize>(array: [(char, float); N]) -> Self {
        Self {
            params: array
                .map(|(name, value)| Param::new(name, value))
                .to_vec()
        }
    }

    pub fn from_vec(vec: Vec<(char, float)>) -> Self {
        Self {
            params: vec.into_iter()
                .map(|(name, value)| Param::new(name, value))
                .collect()
        }
    }

    pub fn gen_random_from_f(f: &Function) -> Self {
        let mut params_names = f.get_params_names();
        params_names.sort();
        params_names.dedup();
        let params = params_names.into_iter()
            .map(|name| Param::gen_random_value_with_name(name))
            .collect();
        Self { params }
    }

    pub fn amount(&self) -> usize {
        self.params.len()
    }

    pub fn get_all(&self) -> Vec<Param> {
        self.params.clone()
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
}

