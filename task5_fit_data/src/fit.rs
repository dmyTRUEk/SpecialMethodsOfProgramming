//! Fits function to points.

use crate::{
    extensions::{Avg, IndexOfMax, IndexOfMinWithFloor},
    fit_params::{FIT_ALGORITHM_MIN_STEP, FIT_ALGORITHM_TYPE, FIT_RESIDUE_EVALS_MAX},
    float_type::float,
    function_and_params::FunctionAndParams,
    params::Params,
    points::Points,
    utils_io::press_enter_to_continue,
};


#[allow(dead_code)]
pub enum FitAlgorithmType {
    PatternSearch,
    DownhillSimplex,
}


pub enum DiffFunctionType {
    DySquared,
    DyAbs,
    LeastDist,
}


#[derive(Debug)]
pub struct FitResults {
    pub fit_residue: float,
    pub fit_residue_evals: u32,
}

// type FitResultsOrError = Result<FitResults, &'static str>;
type FitResultsOrNone = Option<FitResults>;


/// Returns `fit residue` and `fit_residue_evals` it finished in.
pub fn fit(f: &mut FunctionAndParams, points: &Points) -> FitResultsOrNone {
    fit_with_fit_algorith_type(f, points, FIT_ALGORITHM_TYPE)
}

pub fn fit_with_fit_algorith_type(f: &mut FunctionAndParams, points: &Points, fit_algorithm_type: FitAlgorithmType) -> FitResultsOrNone {
    match fit_algorithm_type {
        FitAlgorithmType::PatternSearch => fit_by_pattern_search_algorithm(f, points),
        FitAlgorithmType::DownhillSimplex => fit_by_downhill_simplex_algorithm(f, points),
    }
}


fn fit_by_pattern_search_algorithm(f: &mut FunctionAndParams, points: &Points) -> FitResultsOrNone {
    use crate::patter_search_params::*;
    const DEBUG: bool = false;
    let f_params_amount: usize = f.params.amount();
    let mut step: float = INITIAL_STEP;
    let mut fit_residue_evals = 0;
    if f_params_amount > 0 {
        while step > FIT_ALGORITHM_MIN_STEP && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("f.f = {}", f.f_to_string());
                println!("f.params = {:#?}", f.params);
                println!("step = {}", step);
            }

            let res_at_current_params: float = f.calc_fit_residue(points);
            fit_residue_evals += 1;
            if DEBUG { println!("res_at_current_params = {}", res_at_current_params) }
            // if !res_at_current_params.is_finite() { return Err("`res_at_current_params` isn't finite") }
            if !res_at_current_params.is_finite() { return None }

            let mut residues_at_shifted_params = Vec::<float>::with_capacity(2 * f_params_amount);

            for param in f.params.get_all() {
                // TODO(optimization)?: try mutating self.
                let mut params = f.params.clone();
                // TODO(optimization)?: unroll for loop by hands.
                for delta in [-step, step] {
                    let new_param_value = param.value + delta;
                    // if !new_param_value.is_finite() { return Err(format!("`param.value + delta` isn't finite")) }
                    // TODO(optimization)?: just remove this whole if.
                    if !new_param_value.is_finite() {
                        residues_at_shifted_params.push(float::NAN);
                        continue;
                    }
                    params.set_by_name(param.name, new_param_value);
                    let res = f.calc_fit_residue_with_params(&params, points);
                    fit_residue_evals += 1;
                    residues_at_shifted_params.push(if res.is_finite() { res } else { float::NAN });
                }
            }
            if DEBUG { println!("res_at_shifted_params = {:?}", residues_at_shifted_params) }
            // assert_eq!(2 * f_params_amount, residues_at_shifted_params.len());
            // if res_at_shifted_params.iter().any(|r| !r.is_finite()) { return Err(format!("one of `res_at_shifted_params` isn't finite")) }

            match residues_at_shifted_params.index_of_min_with_ceil(res_at_current_params) {
                Some(index_of_min) => {
                    if DEBUG { println!("INCREASE STEP") }
                    // TODO(optimization)?: use better approach?
                    let param_name = f.params.get_by_index(index_of_min as usize / 2).name;
                    let delta = if index_of_min % 2 == 0 { -step } else { step };
                    f.params.change_param_by(param_name, delta);
                    step *= ALPHA;
                }
                None => {
                    if DEBUG { println!("DECREASE STEP") }
                    step *= BETA;
                }
            }

            if DEBUG { println!("\n\n") }
        }
        if fit_residue_evals == FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("{}", "!".repeat(21));
                println!("HIT MAX_ITERS!!!");
                press_enter_to_continue();
            }
            // return Err("hit max iters");
            return None;
        }
        if DEBUG { println!("finished in {} iters", fit_residue_evals) }
    }
    fit_residue_evals += 1;
    Some(FitResults {
        fit_residue: f.calc_fit_residue(points),
        fit_residue_evals,
    })
}



pub fn fit_by_downhill_simplex_algorithm(f: &mut FunctionAndParams, points: &Points) -> FitResultsOrNone {
    use crate::downhill_simplex_params::*;
    const DEBUG: bool = false;
    const LERP_TS: [float; 15] = [0.5, 0.45, 0.55, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.01, 0.99, 0.001, 0.999];

    fn is_close_enough(params_a: Params, params_b: Params) -> bool {
        let diff = params_a.diff(params_b, PARAMS_DIFF_TYPE);
        diff < FIT_ALGORITHM_MIN_STEP
    }

    let f_params_amount: usize = f.params.amount();
    let mut fit_residue_evals = 0;
    if f_params_amount > 0 {
        let mut params_vec_prev: Vec<Params> = vec![f.params.changed_all_params_by(INITIAL_SIMPLEX_SCALE); f_params_amount+1];
        let mut params_and_ress_vec: Vec<(Params, float)> = Vec::with_capacity(f_params_amount+1);
        trait ExtParamsAndRessVec {
            fn get_params(&self) -> Vec<Params>;
            fn get_residues(&self) -> Vec<float>;
        }
        impl ExtParamsAndRessVec for Vec<(Params, float)> {
            fn get_params(&self) -> Vec<Params> {
                self.iter().map(|p| p.0.clone()).collect()
            }
            fn get_residues(&self) -> Vec<float> {
                self.iter().map(|p| p.1).collect()
            }
        }
        trait ExtParamsVec {
            fn get_all_except(&self, index: usize) -> Vec<Params>;
        }
        impl ExtParamsVec for Vec<Params> {
            fn get_all_except(&self, index: usize) -> Vec<Params> {
                let mut self_ = self.clone();
                self_.remove(index);
                self_
            }
        }
        trait ExtParams {
            fn mirror_relative_to(self, others: Vec<Params>) -> Params;
            fn lerp(self, other: Params, t: float) -> Params;
        }
        impl ExtParams for Params {
            fn mirror_relative_to(self, others: Vec<Params>) -> Params {
                let others_avg: Params = others.avg();
                self.clone() + 2. * (others_avg - self)
            }
            fn lerp(self, other: Params, t: float) -> Params {
                self*t + other*(1.-t)
            }
        }
        let mut params_and_ress_vec_push = |params: Params| {
            let fit_residue = f.calc_fit_residue_with_params(&params, points);
            fit_residue_evals += 1;
            params_and_ress_vec.push((params, fit_residue));
        };
        params_and_ress_vec_push(f.params.changed_all_params_by(-INITIAL_SIMPLEX_SCALE/(f_params_amount as float)));
        for i in 0..f_params_amount {
            let param_name = f.params.get_by_index(i).name;
            params_and_ress_vec_push(f.params.changed_param_by(param_name, INITIAL_SIMPLEX_SCALE));
        }
        // assert_eq!(f_params_amount+1, params_and_ress_vec.len());
        // assert_eq!(f_params_amount+1, params_vec_prev.len());
        // TODO(optimization): `while !is_close_enough` -> `loop { … let step = …; if step < MIN_STEP { break } }`.
        while !is_close_enough(params_and_ress_vec.get_params().avg(), params_vec_prev.avg()) && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
            let index_of_max = params_and_ress_vec.get_residues().index_of_max();
            // if index_of_max.is_none() { return Err("`fit_residue` at all `params_vec` is NaN or Inf") }
            if index_of_max.is_none() { return None }
            let index_of_max = index_of_max.unwrap();

            // same, but in other manner:
            // let index_of_max = match params_and_ress_vec.get_residues().index_of_max() {
            //     None => return Err("function residue at all `params_vec` is NaN or Inf".to_string()),
            //     Some(index_of_max) => index_of_max
            // };

            let (params_max, value_at_params_max) = params_and_ress_vec[index_of_max].clone();
            let params_other: Vec<Params> = params_and_ress_vec.get_params().get_all_except(index_of_max);
            // assert_eq!(f_params_amount, params_other.len());

            let params_symmetric = params_max.clone().mirror_relative_to(params_other.clone());
            let value_at_params_symmetric = f.calc_fit_residue_with_params(&params_symmetric, points);
            fit_residue_evals += 1;

            params_vec_prev = params_and_ress_vec.get_params();

            params_and_ress_vec[index_of_max] = if value_at_params_symmetric < value_at_params_max && value_at_params_symmetric.is_finite() {
                (params_symmetric, value_at_params_symmetric)
            } else {
                let mut option_params_value: Option<(Params, float)> = None;
                for lerp_t in LERP_TS {
                    let params_lerp = params_max.clone().lerp(params_other.clone().avg(), lerp_t);
                    let value_at_params_lerp = f.calc_fit_residue_with_params(&params_lerp, points);
                    fit_residue_evals += 1;
                    if value_at_params_lerp.is_finite() {
                        option_params_value = Some((params_lerp, value_at_params_lerp));
                        break
                    }
                }
                match option_params_value {
                    // None => return Err("`fit_residue` at all `LERP_TS` is NaN or Inf"),
                    None => return None,
                    Some(params_value) => params_value
                }
            };
        }
        if fit_residue_evals == FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("{}", "!".repeat(21));
                println!("HIT MAX_ITERS!!!");
                press_enter_to_continue();
            }
            // return Err("hit max iters");
            return None;
        }
        f.params = params_and_ress_vec.get_params().avg();
    }
    fit_residue_evals += 1;
    Some(FitResults {
        fit_residue: f.calc_fit_residue(points),
        fit_residue_evals,
    })
}





#[cfg(test)]
mod tests {
    use super::*;

    use rand::{thread_rng, Rng};

    use crate::{function::Function, params::Params, point::Point};

    const TOLERANCE: float = 1e-3;

    mod by_pattern_search {
        use super::*;

        #[allow(unused_must_use)]
        #[test]
        fn k_x() {
            let points = vec![
                Point::new(0., 0.),
                Point::new(1., 2.),
                Point::new(2., 4.),
            ];
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let k = rng.gen_range(-5. ..= 5.);
                let mut f = FunctionAndParams::new(
                    Function::from_str("k*x").unwrap(),
                    Params::from_array([('k', k)]),
                );
                fit_with_fit_algorith_type(&mut f, &points, FitAlgorithmType::PatternSearch);
                assert!((2. - f.params.get_by_name_checked('k').unwrap()).abs() < TOLERANCE);
            }
        }

        #[allow(unused_must_use)]
        #[test]
        fn x_x() {
            let points = vec![
                Point::new(0., 0.),
                Point::new(1., 1.),
                Point::new(2., 4.),
                Point::new(3., 9.),
            ];
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let a = rng.gen_range(-5. ..= 5.);
                let mut f = FunctionAndParams::new(
                    Function::from_str("a*x*x").unwrap(),
                    Params::from_array([('a', a)]),
                );
                fit_with_fit_algorith_type(&mut f, &points, FitAlgorithmType::PatternSearch);
                assert!((1. - f.params.get_by_name_checked('a').unwrap()).abs() < TOLERANCE);
            }
        }

        #[allow(unused_must_use)]
        #[test]
        fn x_squared() {
            let points = vec![
                Point::new(0., 0.),
                Point::new(1., 1.),
                Point::new(2., 4.),
                Point::new(3., 9.),
            ];
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let a = rng.gen_range(-5. ..= 5.);
                let mut f = FunctionAndParams::new(
                    Function::from_str("a*x^2").unwrap(),
                    Params::from_array([('a', a)]),
                );
                fit_with_fit_algorith_type(&mut f, &points, FitAlgorithmType::PatternSearch);
                assert!((1. - f.params.get_by_name_checked('a').unwrap()).abs() < TOLERANCE);
            }
        }
    }

    mod by_downhill_simplex {
        use super::*;

        #[allow(unused_must_use)]
        #[test]
        fn k_x() {
            let points = vec![
                Point::new(0., 0.),
                Point::new(1., 2.),
                Point::new(2., 4.),
            ];
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let k = rng.gen_range(-5. ..= 5.);
                let mut f = FunctionAndParams::new(
                    Function::from_str("k*x").unwrap(),
                    Params::from_array([('k', k)]),
                );
                fit_with_fit_algorith_type(&mut f, &points, FitAlgorithmType::DownhillSimplex).unwrap();
                assert!((2. - f.params.get_by_name_checked('k').unwrap()).abs() < TOLERANCE);
            }
        }

        #[allow(unused_must_use)]
        #[test]
        fn x_x() {
            let points = vec![
                Point::new(0., 0.),
                Point::new(1., 1.),
                Point::new(2., 4.),
                Point::new(3., 9.),
            ];
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let a = rng.gen_range(-5. ..= 5.);
                let mut f = FunctionAndParams::new(
                    Function::from_str("a*x*x").unwrap(),
                    Params::from_array([('a', a)]),
                );
                fit_with_fit_algorith_type(&mut f, &points, FitAlgorithmType::DownhillSimplex);
                assert!((1. - f.params.get_by_name_checked('a').unwrap()).abs() < TOLERANCE);
            }
        }

        #[allow(unused_must_use)]
        #[test]
        fn x_squared() {
            let points = vec![
                Point::new(0., 0.),
                Point::new(1., 1.),
                Point::new(2., 4.),
                Point::new(3., 9.),
            ];
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let a = rng.gen_range(-5. ..= 5.);
                let mut f = FunctionAndParams::new(
                    Function::from_str("a*x^2").unwrap(),
                    Params::from_array([('a', a)]),
                );
                fit_with_fit_algorith_type(&mut f, &points, FitAlgorithmType::DownhillSimplex);
                assert!((1. - f.params.get_by_name_checked('a').unwrap()).abs() < TOLERANCE);
            }
        }
    }
}

