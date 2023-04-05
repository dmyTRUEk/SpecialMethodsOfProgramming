//! Simple deconvolution.

use std::{
    env,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use rayon::prelude::{IntoParallelIterator, ParallelIterator};


mod fit_params {
    use super::*;
    pub const INITIAL_VALUES: float = 0.;
    pub const FIT_ALGORITHM_MIN_STEP: float = 1e-4;
    pub const FIT_RESIDUE_EVALS_MAX: u64 = 100_000_000;
    pub const FIT_ALGORITHM_TYPE    : FitAlgorithmType = FitAlgorithmType::PatternSearch;
    pub const RESIDUAL_FUNCTION_TYPE: DiffFunctionType = DiffFunctionType::DySqr;
}

mod patter_search_params {
    use super::*;
    pub const INITIAL_STEP: float = 0.2;
    pub const ALPHA: float = 2.;         // step increase coefficient
    pub const BETA : float = 1. / ALPHA; // step decrease coefficient
}

mod downhill_simplex_params {
    use super::*;
    pub const INITIAL_SIMPLEX_SCALE: float = 0.815;
    pub const PARAMS_DIFF_TYPE: DiffFunctionType = DiffFunctionType::DySqr;
}

const DECONVOLUTION_SOLVER_TYPE: DeconvolutionSolverType = DeconvolutionSolverType::Simple;


fn main() {
    let args: Vec<_> = env::args().collect();
    let (filepathstr_instrument, filepathstr_spectrum): (&str, &str) = match &args[..] {
        [_, filepathstr_instrument, filepathstr_spectrum] => (filepathstr_instrument, filepathstr_spectrum),
        [_, _] => panic!("Expected two filename, provided only one."),
        [_] => panic!("Filenames not provided."),
        [] => unreachable!("Unexpected CLI args number."),
        _ => panic!("Too many CLI args.")
    };

    print!("loading instrument from `{}`...", filepathstr_instrument); flush();
    let points_instrument = load_data_y(filepathstr_instrument);
    println!(" done");

    print!("loading spectrum from `{}`...", filepathstr_spectrum); flush();
    let points_spectrum = load_data_y(filepathstr_spectrum);
    println!(" done");

    // TODO: warning if points in instr more than in spectrum.

    println!("FIT_RESIDUE_EVALS_MAX  = {}", fit_params::FIT_RESIDUE_EVALS_MAX.to_string_beautiful());
    println!("FIT_ALGORITHM_MIN_STEP = {:.2e}", fit_params::FIT_ALGORITHM_MIN_STEP);

    let deconvolve_results = DECONVOLUTION_SOLVER_TYPE.deconvolve(points_instrument, points_spectrum);
    dbg!(&deconvolve_results);
    let deconvolve_results = deconvolve_results.unwrap();
    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals.to_string_beautiful());

    let file_instrument = Path::new(filepathstr_instrument);
    let file_spectrum   = Path::new(filepathstr_spectrum);
    assert_eq!(
        file_instrument.parent().unwrap().canonicalize().unwrap().to_str().unwrap(),
        file_spectrum  .parent().unwrap().canonicalize().unwrap().to_str().unwrap()
    );
    let filepath_output = file_spectrum.with_file_name(format!(
        "results_{}_{}.dat",
        file_instrument.file_stem().unwrap().to_str().unwrap(),
        file_spectrum.file_stem().unwrap().to_str().unwrap()
    ));
    let mut file_output = File::create(filepath_output).unwrap();

    // TODO(refactor): `zip_exact`.
    // assert_eq!((1010..=1089).count(), deconvolve_results.points.len());
    // for (x, point) in (1010..=1089).zip(deconvolve_results.points) {
    //     writeln!(file_output, "{x}\t{p}", p=point).unwrap();
    // }
    use std::io::Write;
    for i in 0..deconvolve_results.points.len() {
        let point = deconvolve_results.points[i];
        writeln!(file_output, "{i}\t{p}", p=point).unwrap();
    }
}


pub fn convolve(points_instrument: &Vec<float>, points_spectrum_original: &Vec<float>) -> Vec<float> {
    assert!(points_instrument.len() % 2 == 1, "points_instrument.len() = {}", points_instrument.len());
    let mut points_convolved = Vec::<float>::with_capacity(points_spectrum_original.len());
    for i in 0..points_spectrum_original.len() {
        let mut point_convolved = 0.;
        for j in 0..points_instrument.len() {
            let d: i32 = j as i32 - points_instrument.len() as i32 / 2;
            let pii: i32 = j as i32;     // points_instrument_index
            let psi: i32 = i as i32 - d; // points_spectrum_index
            let is_pii_in_range: bool = 0 <= pii && pii < points_instrument.len() as i32;
            let is_psi_in_range: bool = 0 <= psi && psi < points_spectrum_original.len() as i32;
            if is_pii_in_range && is_psi_in_range {
                let point_instrument        = points_instrument       [pii as usize];
                let point_spectrum_original = points_spectrum_original[psi as usize];
                point_convolved += point_spectrum_original * point_instrument;
            }
        }
        points_convolved.push(point_convolved);
    }
    assert_eq!(points_spectrum_original.len(), points_convolved.len());
    points_convolved
}


#[allow(dead_code)]
pub enum DeconvolutionSolverType {
    Simple,
    Fourier,
}
impl DeconvolutionSolverType {
    // TODO: make separate types for `points_instrument` & `points_spectrum` to prevent mixing them up.
    pub fn deconvolve(&self, points_instrument: Vec<float>, points_spectrum: Vec<float>) -> FitResultsOrError {
        match DECONVOLUTION_SOLVER_TYPE {
            DeconvolutionSolverType::Simple  => Self::deconvolve_simple (points_instrument, points_spectrum),
            DeconvolutionSolverType::Fourier => Self::deconvolve_fourier(points_instrument, points_spectrum),
        }
    }

    fn deconvolve_simple(points_instrument: Vec<float>, points_spectrum: Vec<float>) -> FitResultsOrError {
        fn residue_function(params: &Vec<float>, points_instrument: &Vec<f64>, points_spectrum: &Vec<f64>) -> float {
            assert_eq!(points_spectrum.len(), params.len());
            let points_convolved = convolve(&points_instrument, params);
            fit_params::RESIDUAL_FUNCTION_TYPE.calc_diff(&points_spectrum, &points_convolved)
        }
        fit_params::FIT_ALGORITHM_TYPE.fit(residue_function, points_instrument, points_spectrum)
    }

    #[allow(unused)]
    fn deconvolve_fourier(points_instrument: Vec<float>, points_spectrum: Vec<float>) -> FitResultsOrError {
        unimplemented!()
    }
}


pub enum DiffFunctionType {
    DySqr,
    DyAbs,
    DySqrPerEl,
    DyAbsPerEl,
    LeastDist,
}
impl DiffFunctionType {
    pub fn calc_diff(&self, points_1: &Vec<float>, points_2: &Vec<float>) -> float {
        assert_eq!(points_1.len(), points_2.len());
        match self {
            DiffFunctionType::DySqr => {
                let mut res: float = 0.;
                for (point_1, point_2) in points_1.into_iter().zip(points_2) {
                    let delta = point_2 - point_1;
                    res += delta.powi(2);
                }
                res.sqrt()
            }
            DiffFunctionType::DyAbs => {
                let mut res: float = 0.;
                for (point_1, point_2) in points_1.into_iter().zip(points_2) {
                    let delta = point_2 - point_1;
                    res += delta.abs();
                }
                res
            }
            DiffFunctionType::DySqrPerEl => DiffFunctionType::DySqr.calc_diff(points_1, points_2) / points_1.len() as float,
            DiffFunctionType::DyAbsPerEl => DiffFunctionType::DyAbs.calc_diff(points_1, points_2) / points_1.len() as float,
            DiffFunctionType::LeastDist => { unimplemented!() }
        }
    }
}

#[derive(Debug)]
pub struct FitResults {
    pub points: Vec<float>,
    pub fit_residue: float,
    pub fit_residue_evals: u64,
}
// type FitResultsOrNone = Option<FitResults>;
type FitResultsOrError = Result<FitResults, &'static str>;

pub enum FitAlgorithmType {
    PatternSearch,
    DownhillSimplex,
}
impl FitAlgorithmType {
    pub fn fit(
        &self,
        residue_function: fn(&Vec<float>, &Vec<float>, &Vec<float>) -> float,
        points_instrument: Vec<float>,
        points_spectrum  : Vec<float>,
    ) -> FitResultsOrError {
        match &self {
            FitAlgorithmType::PatternSearch   => Self::fit_by_pattern_search_algorithm  (residue_function, points_instrument, points_spectrum),
            FitAlgorithmType::DownhillSimplex => Self::fit_by_downhill_simplex_algorithm(residue_function, points_instrument, points_spectrum),
        }
    }

    fn fit_by_pattern_search_algorithm(
        residue_function: fn(&Vec<float>, &Vec<float>, &Vec<float>) -> float,
        points_instrument: Vec<float>,
        points_spectrum  : Vec<float>,
    ) -> FitResultsOrError {
        use crate::{fit_params::*, patter_search_params::*};
        const DEBUG: bool = false;

        let f_params_amount: usize = points_spectrum.len();
        if f_params_amount == 0 {
            return Err("too few params");
            // return None;
        }

        type Params = Vec<float>;
        let mut params: Params = vec![INITIAL_VALUES; f_params_amount];
        let mut step: float = INITIAL_STEP;
        let mut fit_residue_evals: u64 = 0;

        while step > FIT_ALGORITHM_MIN_STEP && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("params = {:#?}", params);
                println!("step = {}", step);
            }

            let res_at_current_params: float = residue_function(&params, &points_instrument, &points_spectrum);
            fit_residue_evals += 1;
            if DEBUG { println!("res_at_current_params = {}", res_at_current_params) }
            if !res_at_current_params.is_finite() { return Err("`res_at_current_params` isn't finite") }
            // if !res_at_current_params.is_finite() { return None }

            // let mut ress_at_shifted_params = Vec::<float>::with_capacity(2 * f_params_amount);
            // for i in 0..params.len() {
            //     let mut params_new = params.clone();
            //     for delta in [-step, step] {
            //         let param_new = params_new[i] + delta;
            //         // if !param_new.is_finite() { return Err("`param.value + delta` isn't finite") }
            //         // TODO(optimization)?: remove `.is_finite()` check, bc it already will be "done" when calculating residue function.
            //         if !param_new.is_finite() || param_new < 0. {
            //             ress_at_shifted_params.push(float::NAN);
            //             continue;
            //         }
            //         let old_param_value = params_new[i];
            //         params_new[i] = param_new;
            //         let res = residue_function(&params_new);
            //         params_new[i] = old_param_value;
            //         fit_residue_evals += 1;
            //         ress_at_shifted_params.push(if res.is_finite() { res } else { float::NAN });
            //     }
            // }

            let (fit_residue_evals_extra, ress_at_shifted_params): (Vec<u64>, Vec<float>) =
                (0..2*params.len())
                    // .into_iter()
                    .into_par_iter()
                    .map(|i| -> (u64, float) {
                        let delta = if i % 2 == 0 { -step } else { step };
                        let param_new = params[i/2] + delta;
                        // if !param_new.is_finite() { return Err("`param.value + delta` isn't finite") }
                        // TODO(optimization)?: remove `.is_finite()` check, bc it already will be "done" when calculating residue function.
                        if !param_new.is_finite() || param_new < 0. {
                            (0, float::NAN)
                        } else {
                            let mut params_new = params.clone();
                            params_new[i/2] = param_new;
                            let res = residue_function(&params_new, &points_instrument, &points_spectrum);
                            (1, if res.is_finite() { res } else { float::NAN })
                        }
                        // returns tuple of `residue_function_evals` and `residue_result`.
                    })
                    .unzip();
            fit_residue_evals += fit_residue_evals_extra.iter().sum::<u64>();

            if DEBUG { println!("res_at_shifted_params = {:?}", ress_at_shifted_params) }
            assert_eq!(2 * f_params_amount, ress_at_shifted_params.len());
            // if res_at_shifted_params.iter().any(|r| !r.is_finite()) { return Err(format!("one of `res_at_shifted_params` isn't finite")) }

            match ress_at_shifted_params.index_of_min_with_ceil(res_at_current_params) {
                Some(index_of_min) => {
                    if DEBUG { println!("INCREASE STEP") }
                    let param_index = index_of_min as usize / 2;
                    let delta = if index_of_min % 2 == 0 { -step } else { step };
                    params[param_index] += delta;
                    step *= ALPHA;
                }
                None => {
                    if DEBUG { println!("DECREASE STEP") }
                    step *= BETA;
                }
            }

            if DEBUG { println!("\n\n") }
        }
        if fit_residue_evals >= FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("{}", "!".repeat(21));
                println!("HIT MAX_ITERS!!!");
                press_enter_to_continue();
            }
            return Err("hit max iters");
            // return None;
        }
        if DEBUG { println!("finished in {} iters", fit_residue_evals) }
        let fit_residue = residue_function(&params, &points_instrument, &points_spectrum);
        fit_residue_evals += 1;
        Ok(FitResults {
            points: params,
            fit_residue,
            fit_residue_evals,
        })
    }


    fn fit_by_downhill_simplex_algorithm(
        residue_function: fn(&Vec<float>, &Vec<float>, &Vec<float>) -> float,
        points_instrument: Vec<float>,
        points_spectrum  : Vec<float>,
    ) -> FitResultsOrError {
        use crate::{downhill_simplex_params::*, fit_params::*};
        const DEBUG: bool = false;
        const LERP_TS: [float; 15] = [0.5, 0.45, 0.55, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.01, 0.99, 0.001, 0.999];

        let f_params_amount: usize = points_spectrum.len();
        if f_params_amount == 0 {
            return Err("too few params");
            // return None;
        }

        type Params = Vec<float>;
        fn is_close_enough(params_a: &Params, params_b: &Params) -> bool {
            let diff = PARAMS_DIFF_TYPE.calc_diff(params_a, params_b);
            diff < FIT_ALGORITHM_MIN_STEP / (params_a.len() as float).powi(1)
        }

        let mut fit_residue_evals = 0;
        let mut params_prev_prev: Params = vec![INITIAL_VALUES+INITIAL_SIMPLEX_SCALE; f_params_amount];
        let mut params_prev_this: Params      = vec![INITIAL_VALUES-INITIAL_SIMPLEX_SCALE; f_params_amount];
        //                                 data  residue
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
            fn is_all_positive(&self) -> bool;
        }
        impl ExtParams for Params {
            fn mirror_relative_to(self, others: Vec<Params>) -> Params {
                let others_avg: Params = others.avg();
                self.add(others_avg.sub(self.clone()).scale(2.))
            }
            fn lerp(self, other: Params, t: float) -> Params {
                self.scale(t).add(other.scale(1.-t))
            }
            fn is_all_positive(&self) -> bool {
                self.iter().all(|&x| x >= 0.)
            }
        }
        let mut params_and_ress_vec_push = |params: Params| {
            let fit_residue = residue_function(&params, &points_instrument, &points_spectrum);
            fit_residue_evals += 1;
            params_and_ress_vec.push((params, fit_residue));
        };
        params_and_ress_vec_push(vec![INITIAL_VALUES-INITIAL_SIMPLEX_SCALE/(f_params_amount as float); f_params_amount]);
        for i in 0..f_params_amount {
            let mut params = vec![INITIAL_VALUES; f_params_amount];
            params[i] += INITIAL_SIMPLEX_SCALE;
            params_and_ress_vec_push(params);
        }
        assert_eq!(f_params_amount+1, params_and_ress_vec.len());
        assert_eq!(f_params_amount, params_prev_this.len());
        while !is_close_enough(&params_prev_this, &params_prev_prev) && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
            let index_of_max = params_and_ress_vec.get_residues().index_of_max();
            if index_of_max.is_none() { return Err("`fit_residue` at all `params_vec` is NaN or Inf") }
            // if index_of_max.is_none() { return None }
            let index_of_max = index_of_max.unwrap();

            let (params_max, value_at_params_max) = params_and_ress_vec[index_of_max].clone();
            let params_other: Vec<Params> = params_and_ress_vec.get_params().get_all_except(index_of_max);
            // assert_eq!(f_params_amount, params_other.len());

            let params_symmetric = params_max.clone().mirror_relative_to(params_other.clone());
            let value_at_params_symmetric = if params_symmetric.is_all_positive() {
                fit_residue_evals += 1;
                residue_function(&params_symmetric, &points_instrument, &points_spectrum)
            } else {
                float::NAN
            };

            params_prev_prev = params_prev_this;
            params_and_ress_vec[index_of_max] = if value_at_params_symmetric < value_at_params_max && value_at_params_symmetric.is_finite() {
                (params_symmetric, value_at_params_symmetric)
            } else {
                let mut option_params_value: Option<(Params, float)> = None;
                for lerp_t in LERP_TS {
                    let params_lerp = params_max.clone().lerp(params_other.clone().avg(), lerp_t);
                    let value_at_params_lerp = if params_lerp.is_all_positive() {
                        fit_residue_evals += 1;
                        residue_function(&params_lerp, &points_instrument, &points_spectrum)
                    } else {
                        float::NAN
                    };
                    if value_at_params_lerp.is_finite() {
                        option_params_value = Some((params_lerp, value_at_params_lerp));
                        break
                    }
                }
                match option_params_value {
                    None => return Err("`fit_residue` at all `LERP_TS` is NaN or Inf"),
                    // None => return None,
                    Some(params_value) => params_value
                }
            };
            params_prev_this = params_and_ress_vec[index_of_max].0.clone();
        }
        if fit_residue_evals >= FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("{}", "!".repeat(21));
                println!("HIT MAX_ITERS!!!");
                press_enter_to_continue();
            }
            return Err("hit max iters");
            // return None;
        }
        let points = params_and_ress_vec.get_params().avg();
        let fit_residue = residue_function(&points, &points_instrument, &points_spectrum);
        fit_residue_evals += 1;
        Ok(FitResults {
            points,
            fit_residue,
            fit_residue_evals,
        })
    }

}


pub fn load_data_y(filename: &str) -> Vec<float> {
    let file = File::open(filename).expect(&format!("Unable to open file: `{}`", filename));
    let lines = BufReader::new(file).lines();
    let mut points = Vec::<float>::with_capacity(20);
    for line in lines {
        let line = line.unwrap();
        let line = line.trim();
        if line == "" { continue }
        let (_, y) = line.split_once(&[' ', '\t']).unwrap();
        let y = y.trim();
        let y = y.replace(",", ".");
        let y = y.parse().unwrap();
        points.push(y);
    }
    points
}


#[allow(non_camel_case_types)]
type float = f64;



trait Math {
    fn add(&self, rhs: Self) -> Self;
    fn sub(&self, rhs: Self) -> Self;
    fn scale(&self, rhs: float) -> Self;
    fn unscale(&self, rhs: float) -> Self;
}
impl Math for Vec<float> {
    fn add(&self, rhs: Self) -> Self {
        assert_eq!(self.len(), rhs.len());
        self.iter().zip(rhs).map(|(l, r)| l + r).collect()
    }
    fn sub(&self, rhs: Self) -> Self {
        assert_eq!(self.len(), rhs.len());
        self.iter().zip(rhs).map(|(l, r)| l - r).collect()
    }
    fn scale(&self, rhs: float) -> Self {
        self.iter().map(|l| l * rhs).collect()
    }
    fn unscale(&self, rhs: float) -> Self {
        self.iter().map(|l| l / rhs).collect()
    }
}



pub trait IndexOfMinMaxWithCeilFloor<T> {
    fn index_of_min_with_ceil (&self, ceil : T) -> Option<usize>;
    fn index_of_min_with_floor(&self, floor: T) -> Option<usize>;
    fn index_of_max_with_ceil (&self, ceil : T) -> Option<usize>;
    fn index_of_max_with_floor(&self, floor: T) -> Option<usize>;
}
impl IndexOfMinMaxWithCeilFloor<float> for Vec<float> {
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
    fn index_of_min_with_floor(&self, floor: float) -> Option<usize> {
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
    fn index_of_max_with_ceil(&self, ceil: float) -> Option<usize> {
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
    fn index_of_max_with_floor(&self, floor: float) -> Option<usize> {
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

pub trait IndexOfMinMax<T> {
    fn index_of_min(&self) -> Option<usize>;
    fn index_of_max(&self) -> Option<usize>;
}
impl IndexOfMinMax<float> for Vec<float> {
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


pub trait Avg<T> {
    /// Calculates average.
    fn avg(self) -> T;
}
impl Avg<Vec<float>> for Vec<Vec<float>> {
    /// Calculates average of elements for dynamic-size array.
    fn avg(self) -> Vec<float> {
        let len: float = self.len() as float;
        let sum: Vec<f64> = self.into_iter().reduce(|acc, el| acc.add(el)).unwrap();
        sum.unscale(len)
    }
}


pub trait ExtStringSeparateChunks {
    fn separate_chunks_from_end(&self, delimiter: impl ToString, chunks_size: usize) -> String;
}
impl ExtStringSeparateChunks for String {
    fn separate_chunks_from_end(&self, delimiter: impl ToString, chunks_size: usize) -> String {
        let len = self.len();
        self.chars()
            // .rev()
            // .enumerate()
            // .map(|(i, c)| if len-i % chunks_size != 0 { c.to_string() } else { format!("{}{}", delimiter.to_string(), c) })
            // .rev()
            .enumerate()
            .map(|(i, c)| if (len-i) % chunks_size != 0 || i == 0 { c.to_string() } else { format!("{}{}", delimiter.to_string(), c) })
            .collect()
    }
}
impl ExtStringSeparateChunks for &str {
    fn separate_chunks_from_end(&self, delimiter: impl ToString, chunks_size: usize) -> String {
        self.to_string().separate_chunks_from_end(delimiter, chunks_size)
    }
}

pub trait ToStringBeautiful {
    fn to_string_beautiful(&self) -> String;
}

impl ToStringBeautiful for u64 {
    fn to_string_beautiful(&self) -> String {
        self.to_string().separate_chunks_from_end("_", 3)
    }
}



pub fn press_enter_to_continue() {
    print("PRESS ENTER TO CONTINUE.");
    wait_for_enter();
}

pub fn wait_for_enter() {
    use std::io::stdin;
    let mut line: String = String::new();
    stdin().read_line(&mut line).unwrap();
}

pub fn print(msg: impl ToString) {
    print!("{}", msg.to_string());
    flush();
}

pub fn flush() {
    use std::io::{Write, stdout};
    stdout().flush().unwrap();
}





#[cfg(test)]
mod tests {

    mod separate_chunks_from_end {
        use super::super::*;
        #[test]
        fn _a() {
            assert_eq!("a", "a".separate_chunks_from_end("_-", 3));
        }
        #[test]
        fn _ab() {
            assert_eq!("ab", "ab".separate_chunks_from_end("_-", 3));
        }
        #[test]
        fn _abc() {
            assert_eq!("abc", "abc".separate_chunks_from_end("_-", 3));
        }
        #[test]
        fn _a_bcd() {
            assert_eq!("a_-bcd", "abcd".separate_chunks_from_end("_-", 3));
        }
        #[test]
        fn _abcdefghijklmnopqrstuvwxyz() {
            assert_eq!("ab_-cde_-fgh_-ijk_-lmn_-opq_-rst_-uvw_-xyz", "abcdefghijklmnopqrstuvwxyz".separate_chunks_from_end("_-", 3));
        }
    }

    mod index_of {
        mod min {
            mod with_ceil {
                mod without_nan {
                    use super::super::super::super::super::*;
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
                    use super::super::super::super::super::*;
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
                    use super::super::super::super::super::*;
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
                    use super::super::super::super::super::*;
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
        mod max {
            mod with_ceil {
                mod without_nan {
                    use super::super::super::super::super::*;
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
                    use super::super::super::super::super::*;
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
                    use super::super::super::super::super::*;
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
                    use super::super::super::super::super::*;
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
    }

    mod convolve {
        use super::super::{DiffFunctionType, convolve, float};
        mod instrument_is_identity {
            use super::*;
            const POINTS_INSTRUMENT: [float; 1] = [1.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 20]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 20], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
        mod instrument_is_delta3 {
            use super::*;
            const POINTS_INSTRUMENT: [float; 3] = [0., 1., 0.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 20]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 20], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
        mod instrument_is_delta7 {
            use super::*;
            const POINTS_INSTRUMENT: [float; 7] = [0., 0., 0., 1., 0., 0., 0.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 20]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 20], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = points_spectrum_original;
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
        mod instrument_is_triangle5 {
            use super::*;
            const POINTS_INSTRUMENT: [float; 5] = [0., 0.5, 1., 0.5, 0.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = [vec![0.; 9], vec![0.5, 1., 0.5], vec![0.; 9]].concat();
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 20]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = [vec![1., 0.5], vec![0.; 19]].concat();
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 20], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = [vec![0.; 19], vec![0.5, 1.]].concat();
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-6;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = [vec![0.; 5], vec![0.5, 1., 0.5], vec![0.; 4], vec![0.5, 1., 0.5], vec![0.; 5]].concat();
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = [vec![1., 0.5], vec![0.; 4], vec![0.5, 1., 0.5], vec![0.; 11]].concat();
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_original = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_convolved_actual = convolve(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
                    let points_convolved_expected = [vec![0.; 11], vec![0.5, 1., 0.5], vec![0.; 4], vec![0.5, 1.]].concat();
                    println!("points_convolved_expected = {:?}", points_convolved_expected);
                    println!("points_convolved_actual = {:?}", points_convolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
                    println!("diff = {}", diff);
                    assert!(diff < 0.1);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
    }

    mod deconvolve_simple {
        use super::super::{DeconvolutionSolverType::Simple as DST, DiffFunctionType, float};
        mod instrument_is_identity {
            use super::*;
            const POINTS_INSTRUMENT: [float; 1] = [1.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1.], vec![0.; 20]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 20], vec![1.]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
        mod instrument_is_delta3 {
            use super::*;
            const POINTS_INSTRUMENT: [float; 3] = [0., 1., 0.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1.], vec![0.; 20]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 20], vec![1.]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
        mod instrument_is_delta7 {
            use super::*;
            const POINTS_INSTRUMENT: [float; 7] = [0., 0., 0., 1., 0., 0., 0.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1.], vec![0.; 20]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 20], vec![1.]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let points_deconvolved_expected = points_spectrum_convolved.clone();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
        mod instrument_is_triangle5 {
            use super::*;
            const POINTS_INSTRUMENT: [float; 5] = [0., 0.5, 1., 0.5, 0.];
            mod original_spectrum_is_delta_21 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 9], vec![0.5, 1., 0.5], vec![0.; 9]].concat();
                    let points_deconvolved_expected = [vec![0.; 10], vec![1.], vec![0.; 10]].concat();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1., 0.5], vec![0.; 19]].concat();
                    let points_deconvolved_expected = [vec![1.], vec![0.; 20]].concat();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 19], vec![0.5, 1.]].concat();
                    let points_deconvolved_expected = [vec![0.; 20], vec![1.]].concat();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
            mod original_spectrum_is_two_deltas_20 {
                use super::*;
                const EPSILON: float = 1e-4;
                #[test]
                fn at_center() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 5], vec![0.5, 1., 0.5], vec![0.; 4], vec![0.5, 1., 0.5], vec![0.; 5]].concat();
                    let points_deconvolved_expected = [vec![0.; 6], vec![1.], vec![0.; 6], vec![1.], vec![0.; 6]].concat();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_left() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![1., 0.5], vec![0.; 4], vec![0.5, 1., 0.5], vec![0.; 11]].concat();
                    let points_deconvolved_expected = [vec![1.], vec![0.; 6], vec![1.], vec![0.; 12]].concat();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
                #[test]
                fn at_right() {
                    println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
                    let points_spectrum_convolved = [vec![0.; 11], vec![0.5, 1., 0.5], vec![0.; 4], vec![0.5, 1.]].concat();
                    let points_deconvolved_expected = [vec![0.; 12], vec![1.], vec![0.; 6], vec![1.]].concat();
                    let deconvolve_results = DST.deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                    let points_deconvolved_actual = deconvolve_results.points;
                    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals);
                    println!("points_deconvolved_expected = {:?}", points_deconvolved_expected);
                    println!("points_deconvolved_actual = {:?}", points_deconvolved_actual);
                    let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_deconvolved_expected, &points_deconvolved_actual);
                    assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
                }
            }
        }
    }

}

