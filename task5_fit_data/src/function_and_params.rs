//! Function and Params struct.

use crate::{
    RESIDUAL_FUNCTION_TYPE,
    fit::ResidualFunctionType,
    float_type::float,
    function::Function,
    params::{ImplParams, Params},
    points::Points,
};


#[derive(Debug, Clone, PartialEq)]
pub struct FunctionAndParams {
    pub f: Function,
    pub params: Params,
}

impl FunctionAndParams {
    pub fn new(f: Function, params: Params) -> Self {
        Self { f, params }
    }

    pub fn gen_from_f(f: Function) -> Self {
        let params = Params::gen_from_f(&f);
        Self::new(f, params)
    }

    pub fn gen(complexity: u32) -> Self {
        let f = Function::gen(complexity);
        Self::gen_from_f(f)
    }

    // pub fn get_f(&self) -> Function {
    //     self.f.clone()
    // }
    // pub fn get_params(&self) -> Params {
    //     self.params.clone()
    // }
    // pub fn get_params_amount(&self) -> usize {
    //     self.params.len()
    // }
    // pub fn get_param(&self, name: ParamName) -> ParamValue {
    //     self.params.get(name)
    // }
    // pub fn get_param_by_index(&self, index: usize) -> Param {
    //     self.params[index].clone()
    // }

    // pub fn set_params(&mut self, params: Params) {
    //     self.params = params;
    // }
    // pub fn change_param_by(&mut self, name: ParamName, delta: ParamValue) {
    //     self.params.change_param_by(name, delta)
    // }

    pub fn f_to_string(&self) -> String {
        self.f.to_string()
    }

    #[allow(dead_code)]
    pub fn eval(&self, x: float) -> float {
        self.f.eval(x, &self.params)
    }

    pub fn eval_with_params(&self, params: &Params, x: float) -> float {
        self.f.eval(x, &params)
    }

    pub fn calc_fit_residue(&self, points: &Points) -> float {
        self.calc_fit_residue_with_params(&self.params, points)
    }

    pub fn calc_fit_residue_with_params(&self, params: &Params, points: &Points) -> float {
        const DEBUG: bool = false;
        match RESIDUAL_FUNCTION_TYPE {
            ResidualFunctionType::LeastSquares => {
                let mut res = 0.;
                for point in points {
                    let dy = self.eval_with_params(params, point.x) - point.y;
                    if DEBUG { println!("calc_fit_residue_with_params: dy = {}", dy) }
                    res += dy.powi(2);
                    if DEBUG { println!("calc_fit_residue_with_params: res = {}", res) }
                }
                res
            }
            ResidualFunctionType::LeastAbs => {
                let mut res = 0.;
                for point in points {
                    let dy = self.eval_with_params(params, point.x) - point.y;
                    res += dy.abs();
                }
                res
            }
            ResidualFunctionType::LeastDist => { todo!() }
        }
    }

    pub fn simplify(self) -> Self {
        let new_f = self.f.simplify();
        // TODO: don't randomize params' values, bc this is pretty annoying sideeffect.
        Self::new(new_f, self.params)
    }
}

impl ToString for FunctionAndParams {
    fn to_string(&self) -> String {
        let f_to_string = self.f_to_string();
        let params_str = self.params.iter()
            .map(|p| format!("{n} = {v}", n=p.name, v=p.value))
            .reduce(|acc, el| format!("{acc}, {el}"));
        match params_str {
            Some(params_str) => format!("{}, where {}", f_to_string, params_str),
            None => f_to_string,
        }
    }
}

pub trait ToStringForPlot {
    fn to_string_for_plot(&self) -> String;
}
impl ToStringForPlot for FunctionAndParams {
    fn to_string_for_plot(&self) -> String {
        let f_to_string = format!("f(x) = {}", self.f_to_string());
        let params_str = self.params.iter()
            .map(|p| format!("{n} = {v}", n=p.name, v=p.value))
            .reduce(|acc, el| format!("{acc}\n{el}"));
        match params_str {
            Some(params_str) => format!("{}\n{}", f_to_string, params_str),
            None => f_to_string,
        }
    }
}

