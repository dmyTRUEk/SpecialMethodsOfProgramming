//! Params type.

use crate::{
    float_type::float,
    function::Function,
    param::{Param, ParamName},
};


pub type Params = Vec<Param>;

// TODO(refactor): rewrite to normal struct.
// struct Params {
//     params: Vec<Param>
// }

pub trait ImplParams {
    #[allow(non_snake_case)]
    fn empty() -> Self;
    fn from_array<const N: usize>(array: [(char, float); N]) -> Self;
    fn gen_from_f(f: &Function) -> Self;
    fn get(&self, name: ParamName) -> float;
    fn set(&mut self, name: ParamName, new_value: float);
    fn change_param_by(&mut self, name: ParamName, delta: float);
}
impl ImplParams for Params {
    fn empty() -> Self { vec![] }

    fn from_array<const N: usize>(array: [(char, float); N]) -> Self {
        array
            .map(|(name, value)| Param::new(name, value))
            .to_vec()
    }

    fn gen_from_f(f: &Function) -> Self {
        let mut params_names = f.get_params_names();
        params_names.sort();
        params_names.dedup();
        let params = params_names.into_iter()
            .map(|name| Param::gen_value(name))
            .collect();
        params
    }

    fn get(&self, name: ParamName) -> float {
        self.iter()
            .find(|param| param.name == name)
            // .expect(&format!("parameter `{name}` not found in params: {params:?}"))
            .unwrap()
            .value
    }

    fn set(&mut self, name: ParamName, new_value: float) {
        for i in 0..self.len() {
            if self[i].name == name {
                self[i].value = new_value;
                return;
            }
        }
    }

    fn change_param_by(&mut self, name: ParamName, delta: float) {
        for i in 0..self.len() {
            if self[i].name == name {
                self[i].value += delta;
                return;
            }
        }
    }
}

