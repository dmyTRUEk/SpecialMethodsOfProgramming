//! Fits function to points.

use crate::{
    FIT_ALGORITHM_TYPE,
    FIT_MAX_ITERS,
    MIN_STEP,
    extensions::IndexOfMinWithFloor,
    float_type::float,
    function_and_params::FunctionAndParams,
    params::ImplParams,
    points::Points,
    utils_io::press_enter_to_continue,
};


pub enum FitAlgorithmType {
    PatternSearch,
    DonwhillSimplex,
}


pub enum ResidualFunctionType {
    LeastSquares,
    LeastAbs,
    LeastDist,
}


type FitResult = Result<(float, u32), String>;


/// Returns fit residue and iters it finished in.
pub fn fit(f: &mut FunctionAndParams, points: &Points) -> FitResult {
    match FIT_ALGORITHM_TYPE {
        FitAlgorithmType::PatternSearch => fit_by_pattern_search_algorithm(f, points),
        FitAlgorithmType::DonwhillSimplex => fit_by_downhill_simplex_algorithm(f, points),
    }
}



fn fit_by_pattern_search_algorithm(f: &mut FunctionAndParams, points: &Points) -> FitResult {
    const DEBUG: bool = false;
    const ALPHA: float = 2.;         // step increase coefficient
    const BETA : float = 1. / ALPHA; // step decrease coefficient
    let mut step: float = 1.;
    let mut iters = 0;
    if f.params.len() > 0 {
        while step > MIN_STEP && iters < FIT_MAX_ITERS {
            iters += 1;
            if DEBUG {
                println!("f.f = {}", f.f_to_string());
                println!("f.params = {:?}", f.params);
                println!("step = {}", step);
            }

            let res_at_current_params: float = f.calc_fit_residue(points);
            if DEBUG { println!("res_at_current_params = {}", res_at_current_params) }
            if !res_at_current_params.is_finite() { return Err(format!("`res_at_current_params` isn't finite")) }

            let mut residues_at_shifted_params = Vec::<float>::with_capacity(2 * f.params.len());

            for i in 0..f.params.len() {
                let param = f.params[i].clone();
                // TODO(optimize)?: try mutating self.
                let mut params = f.params.clone();
                // TODO(optimize)?: unroll for loop by hands.
                for delta in [-step, step] {
                    let new_param_value = param.value + delta;
                    // if !new_param_value.is_finite() { return Err(format!("`param.value + delta` isn't finite")) }
                    if !new_param_value.is_finite() {
                        residues_at_shifted_params.push(float::NAN);
                        continue;
                    }
                    params.set(param.name, new_param_value);
                    let res = f.calc_fit_residue_with_params(&params, points);
                    residues_at_shifted_params.push(if res.is_finite() { res } else { float::NAN });
                }
            }
            if DEBUG { println!("res_at_shifted_params = {:?}", residues_at_shifted_params) }
            assert_eq!(2 * f.params.len(), residues_at_shifted_params.len());
            // if res_at_shifted_params.iter().any(|r| !r.is_finite()) { return Err(format!("one of `res_at_shifted_params` isn't finite")) }

            match residues_at_shifted_params.index_of_min_with_ceil(res_at_current_params) {
                None => {
                    if DEBUG { println!("DECREASE STEP") }
                    step *= BETA;
                }
                Some(index_of_min) => {
                    if DEBUG { println!("INCREASE STEP") }
                    let param_name = f.params[index_of_min / 2].name;
                    let delta = if index_of_min % 2 == 0 { -step } else { step };
                    f.params.change_param_by(param_name, delta);
                    step *= ALPHA;
                }
            }
            if DEBUG { println!("\n\n") }
        }
        if iters == FIT_MAX_ITERS {
            if DEBUG {
                println!("{}", "!".repeat(21));
                println!("HIT MAX_ITERS!!!");
                press_enter_to_continue();
            }
            return Err(format!("hit max iters"));
        }
        if DEBUG { println!("finished in {} iters", iters) }
    }
    Ok((f.calc_fit_residue(points), iters))
}



#[allow(unused_variables)]
pub fn fit_by_downhill_simplex_algorithm(f: &mut FunctionAndParams, points: &Points) -> FitResult {
    todo!()
}





#[cfg(test)]
mod tests {
    use super::*;

    use rand::{thread_rng, Rng};

    use crate::{function::Function, param::Param, params::Params, point::Point};

    const TOLERANCE: float = 1e-3;

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
            fit(&mut f, &points);
            assert!((2. - f.params[0].value).abs() < TOLERANCE);
        }
    }

    #[allow(unused_must_use)]
    #[test]
    fn x_x() {
        let mut f = FunctionAndParams::new(
            Function::Mul {
                lhs: box Function::Param { name: 'a' },
                rhs: box Function::Mul {
                    lhs: box Function::X,
                    rhs: box Function::X
                }
            },
            vec![ Param::new('a', 0.5) ],
        );
        let points = vec![
            Point::new(0., 0.),
            Point::new(1., 1.),
            Point::new(2., 4.),
            Point::new(3., 9.),
        ];
        fit(&mut f, &points);
        assert!((1. - f.params[0].value).abs() < TOLERANCE);
    }

    #[allow(unused_must_use)]
    #[test]
    fn x_squared() {
        let mut f = FunctionAndParams::new(
            Function::Mul {
                lhs: box Function::Param { name: 'a' },
                rhs: box Function::Sq {
                    value: box Function::X
                }
            },
            vec![ Param::new('a', 0.5) ],
        );
        let points = vec![
            Point::new(0., 0.),
            Point::new(1., 1.),
            Point::new(2., 4.),
            Point::new(3., 9.),
        ];
        fit(&mut f, &points);
        assert!((1. - f.params[0].value).abs() < TOLERANCE);
    }
}

