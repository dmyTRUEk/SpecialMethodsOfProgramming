//! Deconvolution.

use std::{
    env,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use rayon::prelude::{IntoParallelIterator, ParallelIterator};


mod deconvolution_params {
    use super::*;
    // pub const DECONVOLUTION: Deconvolution = Deconvolution::PerPoint {
    //     diff_function_type: DiffFunctionType::DySqr,
    //     // antispikes: None,
    //     antispikes: Some(Antispikes {
    //         antispikes_type: AntispikesType::DySqr,
    //         antispikes_k: 1.0,
    //     })
    // };
    // pub const DECONVOLUTION_SOLVER_TYPE: DeconvolutionType = DeconvolutionType::PerPoint;
    // pub const FUNCTION_TO_MINIMIZE: FunctionToMinimize = FunctionToMinimize {
    //     diff_function_type: DiffFunctionType::DySqr,
    //     // antispikes: None
    //     antispikes: Some(Antispikes {
    //         antispikes_type: AntispikesType::DySqr,
    //         antispikes_k: 1.0
    //     })
    // };
    pub const DECONVOLUTION: Deconvolution = Deconvolution::Exponents {
        diff_function_type: DiffFunctionType::DySqr,
        exponents_amount: 2,
        // initial_values: None,
        initial_values: Some(&[
            0.01, 2.0, 0.2,
            0.007, -1.5, 0.2,
            // 0.007, -3.0, 0.2,
        ]),
    };
}

mod fit_params {
    use super::*;
    pub const INITIAL_VALUES: float = 0.0015;
    pub const FIT_ALGORITHM_TYPE: FitAlgorithmType = FitAlgorithmType::PatternSearch;
    // pub const FIT_RESIDUE_GOAL   : float = 1e-1; // for Pattern Search
    pub const FIT_ALGORITHM_MIN_STEP: float = 1e-4; // for Pattern Search & Downhill Simplex
    pub const FIT_RESIDUE_EVALS_MAX : u64 = 100_000_000;
}

mod patter_search_params {
    use super::*;
    pub const INITIAL_STEP: float = 0.2;
    pub const ALPHA: float = 1.1;        // step increase coefficient
    pub const BETA : float = 1. / ALPHA; // step decrease coefficient
}

mod downhill_simplex_params {
    use super::*;
    pub const INITIAL_SIMPLEX_SCALE: float = 0.815;
    pub const PARAMS_DIFF_TYPE: DiffFunctionType = DiffFunctionType::DySqr;
}


fn main() {
    let args: Vec<_> = env::args().collect();
    let (filepathstr_instrument, filepathstr_spectrum): (&str, &str) = match &args[..] {
        [_, filepathstr_instrument, filepathstr_spectrum] => (filepathstr_instrument, filepathstr_spectrum),
        [_, _] => panic!("Expected two filename, provided only one."),
        [_] => panic!("Filenames not provided."),
        [] => unreachable!("Unexpected CLI args number."),
        _ => panic!("Too many CLI args.")
    };

    print!("Loading instrumental spectrum  from `{}`...", filepathstr_instrument); flush();
    let points_instrument = load_data_y(filepathstr_instrument);
    println!(" done");

    print!("Loading spectrum to deconvolve from `{}`...", filepathstr_spectrum); flush();
    let points_spectrum = load_data_y(filepathstr_spectrum);
    println!(" done");
    let points_spectrum_len: usize = points_spectrum.len();

    // TODO: warning if points in instr more than in spectrum.
    assert!(points_spectrum.len() > points_instrument.len());

    println!("FIT_ALGORITHM_TYPE    : {:#?}", fit_params::FIT_ALGORITHM_TYPE);
    println!("FIT_ALGORITHM_MIN_STEP: {:.2e}", fit_params::FIT_ALGORITHM_MIN_STEP);
    // if fit_params::FIT_ALGORITHM_TYPE == FitAlgorithmType::PatternSearch {
    //     println!("FIT_RESIDUE_GOAL     : {:.2e}", fit_params::FIT_RESIDUE_GOAL);
    // }
    println!("FIT_RESIDUE_EVALS_MAX : {}", fit_params::FIT_RESIDUE_EVALS_MAX.to_string_underscore_separated());

    let deconvolution_data: DeconvolutionData = DeconvolutionData::new(
        points_instrument,
        points_spectrum,
        deconvolution_params::DECONVOLUTION,
    );
    let deconvolve_results = deconvolution_data.deconvolve(fit_params::FIT_ALGORITHM_TYPE);
    dbg!(&deconvolve_results);
    let deconvolve_results = deconvolve_results.unwrap();
    println!("fit_residue_evals = {}", deconvolve_results.fit_residue_evals.to_string_underscore_separated());

    if let Deconvolution::Exponents { exponents_amount, .. } = deconvolution_params::DECONVOLUTION {
        let mut s: String = [r"f_{", &exponents_amount.to_string(), r"}\left(x\right)="].concat();
        s += &deconvolve_results.params
            .chunks(3).into_iter()
            .map(|parts| {
                let (amplitude, tau, shift) = (parts[0], parts[1], parts[2]);
                [
                    r"\left\{\frac{x}{",
                    &format!("{points_spectrum_len:.3}"),
                    r"}",
                    if tau < 0. { ">" } else { "<" },
                    &format!("{shift:.3}"),
                    r":",
                    &format!("{amplitude:.3}"),
                    r"e^{",
                    &format!("{tau:.3}"),
                    r"\left(\frac{x}{",
                    &format!("{points_spectrum_len:.3}"),
                    r"}-",
                    &format!("{shift:.3}"),
                    r"\right)},0\right\}",
                ].concat()
            })
            .reduce(|acc, el| format!("{acc}+{el}")).unwrap();
        println!("{s}");
    }

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

    use std::io::Write;
    // assert_eq!((1010..=1089).count(), deconvolve_results.points.len());
    // // TODO(refactor): `zip_exact`.
    // for (x, point) in (1010..=1089).zip(deconvolve_results.points) {
    //     writeln!(file_output, "{x}\t{p}", p=point).unwrap();
    // }
    for i in 0..deconvolve_results.params.len() {
        let point = deconvolve_results.params[i];
        writeln!(file_output, "{i}\t{p}", p=point).unwrap();
    }
    drop(file_output);

    let filepath_output_convolved = file_spectrum.with_file_name(format!(
        "results_{}_{}_convolved.dat",
        file_instrument.file_stem().unwrap().to_str().unwrap(),
        file_spectrum.file_stem().unwrap().to_str().unwrap()
    ));
    let mut file_output_convolved = File::create(filepath_output_convolved).unwrap();
    let convolved: Vec<float> = match deconvolution_data.deconvolution {
        Deconvolution::PerPoint { .. } => convolve_per_point(&deconvolution_data.points_instrument, &deconvolve_results.params),
        Deconvolution::Exponents { .. } => convolve_exponents(&deconvolution_data.points_instrument, &deconvolve_results.params, deconvolution_data.points_spectrum.len()),
        Deconvolution::Fourier {} => unimplemented!(),
    };
    for i in 0..convolved.len() {
        let point = convolved[i];
        writeln!(file_output_convolved, "{i}\t{p}", p=point).unwrap();
    }
    drop(file_output_convolved);
}


type DeconvolutionResultOrError = FitResultsOrError;

/// Deconvolution type and it's corresponding params.
pub enum Deconvolution {
    /// aka Simple
    PerPoint {
        diff_function_type: DiffFunctionType,
        antispikes: Option<Antispikes>,
    },
    Exponents {
        diff_function_type: DiffFunctionType,
        exponents_amount: usize,
        initial_values: Option<&'static [float]>
    },
    Fourier {
        // unimplemented
    },
}


pub struct DeconvolutionData {
    points_instrument: Vec<float>,
    points_spectrum: Vec<float>,
    // TODO(feat): points_instrument & points_spectrum step.
    deconvolution: Deconvolution,
}
impl DeconvolutionData {
    pub const fn new(
        points_instrument: Vec<float>,
        points_spectrum: Vec<float>,
        deconvolution: Deconvolution,
    ) -> Self {
        Self {
            points_instrument,
            points_spectrum,
            deconvolution,
        }
    }

    pub fn deconvolve(&self, fit_algorithm_type: FitAlgorithmType) -> DeconvolutionResultOrError {
        fit_algorithm_type.fit(&self)
    }

    /// depending on the `self.deconvolution` `params` is:
    /// - PerPoint: list of values at that point
    /// - Exponents: list of (amplitude, tau, shift)
    /// - Fourier: unimplemented
    pub fn calc_residue_function(&self, params: &Vec<float>) -> float {
        match &self.deconvolution {
            Deconvolution::PerPoint { diff_function_type, antispikes } => {
                assert_eq!(self.points_spectrum.len(), params.len());
                let points_spectrum_original = params;
                let points_convolved = convolve_per_point(&self.points_instrument, points_spectrum_original);
                diff_function_type.calc_diff_with_antispikes(&self.points_spectrum, &points_convolved, antispikes)
            }
            Deconvolution::Exponents { diff_function_type, exponents_amount, .. } => {
                assert_eq!(exponents_amount * 3, params.len());
                // let exponents: Vec<ExponentFunction> = params
                //     .chunks(3).into_iter()
                //     .map(|parts| ExponentFunction::new(parts[0], parts[1], parts[2]))
                //     .collect();
                // assert_eq!(*exponents_amount, exponents.len());
                let points_convolved = convolve_exponents(&self.points_instrument, &params, self.points_spectrum.len());
                // let points_convolved = convolve_exponents(&self.points_instrument, params, self.points_spectrum.len());
                diff_function_type.calc_diff(&self.points_spectrum, &points_convolved)
            }
            Deconvolution::Fourier {} => {
                unimplemented!()
            }
        }
    }

    pub fn get_params_amount(&self) -> usize {
        match self.deconvolution {
            Deconvolution::PerPoint { .. } => self.points_spectrum.len(),
            Deconvolution::Exponents { exponents_amount, .. } => exponents_amount * 3,
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

    pub fn get_default_params(&self) -> Vec<float> {
        match &self.deconvolution {
            Deconvolution::PerPoint { .. } => vec![fit_params::INITIAL_VALUES; self.get_params_amount()],
            Deconvolution::Exponents { initial_values, .. } => {
                let initial_values = match initial_values {
                    None => vec![fit_params::INITIAL_VALUES; self.get_params_amount()],
                    Some(initial_values) => initial_values.to_vec(),
                };
                assert_eq!(self.get_params_amount(), initial_values.len());
                initial_values
            }
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

    pub fn is_params_ok(&self, params: &Vec<float>) -> bool {
        match self.deconvolution {
            Deconvolution::PerPoint { .. } => params.into_iter().all(|&x| x >= 0.),
            Deconvolution::Exponents { .. } => params.into_iter().enumerate().all(|(i, &x)| match i % 3 {
                0 => x >= 0.,
                1 => true,
                2 => true,
                _ => unreachable!()
            }),
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

}


pub fn convolve_exponents(
    points_instrument: &Vec<float>,
    // exponents: &Vec<ExponentFunction>, // TODO(refactor,optimization): remove `Vec<ExponentFunction>`, make just `Vec<float>` and iter them by `.chunks()`.
    params: &Vec<float>,
    points_spectrum_convolved_len: usize,
) -> Vec<float> {
    let mut points_spectrum_original = vec![0.; points_spectrum_convolved_len];
    let exponents: Vec<ExponentFunction> = params
        .chunks(3).into_iter()
        .map(|parts| ExponentFunction::new(parts[0], parts[1], parts[2]))
        .collect();
    for i in 0..points_spectrum_convolved_len {
        let sum: float = exponents.iter()
            .map(|exponent| exponent.eval_at(i as float / points_spectrum_convolved_len as float))
            .sum();
        // let sum: float = params
        //     .chunks(3).into_iter()
        //     .map(|parts| {
        //         let x: float = i as float / points_spectrum_convolved_len as float;
        //         // let exp_func = ExponentFunction::new(parts[0], parts[1], parts[2]);
        //         // exp_func.eval_at(x)
        //         ExponentFunction::eval_at_static(parts[0], parts[1], parts[2], x)
        //     })
        //     .sum();
        points_spectrum_original[i] = sum;
    }
    let points_convolved = convolve_per_point(points_instrument, &points_spectrum_original);
    points_convolved
}


pub fn convolve_per_point(points_instrument: &Vec<float>, points_spectrum_original: &Vec<float>) -> Vec<float> {
    assert!(points_instrument.len() % 2 == 1, "points_instrument.len() = {}", points_instrument.len()); // why?
    let mut points_convolved = vec![0.; points_spectrum_original.len()];
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
                point_convolved += point_instrument * point_spectrum_original;
            }
        }
        points_convolved[i] = point_convolved;
    }
    points_convolved
}


#[derive(Clone)]
pub struct Antispikes {
    antispikes_type: AntispikesType,
    antispikes_k: float,
}
impl Antispikes {
    fn calc(&self, points_1: &Vec<float>, points_2: &Vec<float>) -> float {
        self.antispikes_k * self.antispikes_type.calc(points_1, points_2)
    }
}


#[derive(Clone)]
pub enum AntispikesType {
    DySqr,
    DyAbs,
}
impl AntispikesType {
    fn calc(&self, points_1: &Vec<float>, points_2: &Vec<float>) -> float {
        assert_eq!(points_1.len(), points_2.len());
        match self {
            Self::DySqr => {
                let mut res: float = 0.;
                for points in [points_1, points_2] {
                    for point_prev_next in points.windows(2).into_iter() {
                        let (point_prev, point_next) = (point_prev_next[0], point_prev_next[1]);
                        let delta = point_next - point_prev;
                        let delta = delta.powi(2); // TODO(refactor): rename var
                        res += delta;
                    }
                }
                res.sqrt()
            }
            Self::DyAbs => {
                let mut res: float = 0.;
                for points in [points_1, points_2] {
                    for point_prev_next in points.windows(2).into_iter() {
                        let (point_prev, point_next) = (point_prev_next[0], point_prev_next[1]);
                        let delta = point_next - point_prev;
                        let delta = delta.abs(); // TODO(refactor): rename var
                        res += delta;
                    }
                }
                res
            }
        }
    }
}


pub struct ExponentFunction {
    pub amplitude: float,
    pub tau: float,
    pub shift: float,
}
impl ExponentFunction {
    pub const fn new(amplitude: float, tau: float, shift: float) -> Self {
        Self { amplitude, tau, shift }
    }
    pub fn eval_at(&self, x: float) -> float {
        let in_exp: float = self.tau * (x - self.shift);
        // if (if self.tau <= 0. { 1. } else { -1. }) * (x - self.shift) > 0. {
        if in_exp <= 0. {
            // self.amplitude * exp(self.tau * (x - self.shift))
            self.amplitude * exp(in_exp)
        } else {
            0.
        }
    }
    // pub fn eval_at_static(amplitude: float, tau: float, shift: float, x: float) -> float {
    //     let in_exp: float = tau * (x - shift);
    //     // if (if self.tau <= 0. { 1. } else { -1. }) * (x - self.shift) > 0. {
    //     if in_exp <= 0. {
    //         // self.amplitude * exp(self.tau * (x - self.shift))
    //         amplitude * exp(in_exp)
    //     } else {
    //         0.
    //     }
    // }
}


#[derive(Clone, Copy)]
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
            Self::DySqr => {
                let mut res: float = 0.;
                for (point_1, point_2) in points_1.into_iter().zip(points_2) {
                    let delta = point_2 - point_1;
                    let delta = delta.powi(2); // TODO(refactor): rename var
                    res += delta;
                }
                res.sqrt()
            }
            Self::DyAbs => {
                let mut res: float = 0.;
                for (point_1, point_2) in points_1.into_iter().zip(points_2) {
                    let delta = point_2 - point_1;
                    let delta = delta.abs(); // TODO(refactor): rename var
                    res += delta;
                }
                res
            }
            Self::DySqrPerEl => Self::DySqr.calc_diff(points_1, points_2) / points_1.len() as float,
            Self::DyAbsPerEl => Self::DyAbs.calc_diff(points_1, points_2) / points_1.len() as float,
            Self::LeastDist => { unimplemented!() }
        }
    }

    pub fn calc_diff_with_antispikes(&self, points_1: &Vec<float>, points_2: &Vec<float>, antispikes: &Option<Antispikes>) -> float {
        let diff_main: float = self.calc_diff(points_1, points_2);
        let diff_antispikes: float = antispikes.as_ref().map_or(
            0.,
            |antispikes| antispikes.calc(points_1, points_2)
        );
        diff_main + diff_antispikes
    }
}


#[derive(Debug)]
pub struct FitResults {
    pub params: Vec<float>,
    pub fit_residue: float,
    pub fit_residue_evals: u64,
}
// type FitResultsOrNone = Option<FitResults>;
type FitResultsOrError = Result<FitResults, &'static str>;

#[derive(Debug)]
pub enum FitAlgorithmType {
    PatternSearch,
    DownhillSimplex,
}
impl FitAlgorithmType {
    pub fn fit(&self, deconvolution_data: &DeconvolutionData) -> FitResultsOrError {
        match &self {
            FitAlgorithmType::PatternSearch   => Self::fit_by_pattern_search_algorithm  (deconvolution_data),
            FitAlgorithmType::DownhillSimplex => Self::fit_by_downhill_simplex_algorithm(deconvolution_data),
        }
    }

    fn fit_by_pattern_search_algorithm(deconvolution_data: &DeconvolutionData) -> FitResultsOrError {
        use crate::{fit_params::*, patter_search_params::*};
        const DEBUG: bool = false;

        let f_params_amount: usize = deconvolution_data.get_params_amount();
        if f_params_amount == 0 {
            return Err("too few params");
            // return None;
        }

        type Params = Vec<float>;
        let mut params: Params = deconvolution_data.get_default_params();
        let mut step: float = INITIAL_STEP;
        let mut fit_residue_evals: u64 = 0;

        while step > FIT_ALGORITHM_MIN_STEP && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
        // while residue_function(&params, &points_instrument, &points_spectrum) > FIT_RESIDUE_GOAL && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
            if DEBUG {
                println!("params = {:#?}", params);
                println!("step = {}", step);
            }

            let res_at_current_params: float = deconvolution_data.calc_residue_function(&params);
            fit_residue_evals += 1;
            if DEBUG { println!("res_at_current_params = {}", res_at_current_params) }
            if !res_at_current_params.is_finite() { return Err("`res_at_current_params` isn't finite") }
            // if !res_at_current_params.is_finite() { return None }

            let (fit_residue_evals_extra, ress_at_shifted_params): (Vec<u64>, Vec<float>) =
                (0..2*params.len())
                    // .into_iter()
                    .into_par_iter()
                    .map(|i| -> (u64, float) {
                        let delta = if i % 2 == 0 { -step } else { step };
                        let param_new = params[i/2] + delta;
                        // if !param_new.is_finite() { return Err("`param.value + delta` isn't finite") }
                        // TODO(optimization)?: remove `.is_finite()` check, bc it already will be "done" when calculating residue function.
                        let mut params_new = params.clone();
                        params_new[i/2] = param_new;
                        if !param_new.is_finite() || !deconvolution_data.is_params_ok(&params_new) {
                            (0, float::NAN)
                        } else {
                            let res = deconvolution_data.calc_residue_function(&params_new);
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
            return Err("hit max evals");
            // return None;
        }
        if DEBUG { println!("finished in {} iters", fit_residue_evals) }
        let fit_residue = deconvolution_data.calc_residue_function(&params);
        fit_residue_evals += 1;
        Ok(FitResults {
            params,
            fit_residue,
            fit_residue_evals,
        })
    }


    fn fit_by_downhill_simplex_algorithm(deconvolution_data: &DeconvolutionData) -> FitResultsOrError {
        use crate::{downhill_simplex_params::*, fit_params::*};
        const DEBUG: bool = false;
        const LERP_TS: [float; 15] = [0.5, 0.45, 0.55, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9, 0.01, 0.99, 0.001, 0.999];

        let f_params_amount: usize = deconvolution_data.get_params_amount();
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
        todo!("rewrite using new architecture");
        let mut params_prev_prev: Params = vec![INITIAL_VALUES+INITIAL_SIMPLEX_SCALE; f_params_amount];
        let mut params_prev_this: Params = vec![INITIAL_VALUES-INITIAL_SIMPLEX_SCALE; f_params_amount];
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
        }
        impl ExtParams for Params {
            fn mirror_relative_to(self, others: Vec<Params>) -> Params {
                let others_avg: Params = others.avg();
                self.add(others_avg.sub(self.clone()).scale(2.))
            }
            fn lerp(self, other: Params, t: float) -> Params {
                self.scale(t).add(other.scale(1.-t))
            }
        }
        let mut params_and_ress_vec_push = |params: Params| {
            let fit_residue = deconvolution_data.calc_residue_function(&params);
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
            let value_at_params_symmetric = if deconvolution_data.is_params_ok(&params_symmetric) {
                fit_residue_evals += 1;
                deconvolution_data.calc_residue_function(&params_symmetric)
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
                    let value_at_params_lerp = if deconvolution_data.is_params_ok(&params_lerp) {
                        fit_residue_evals += 1;
                        deconvolution_data.calc_residue_function(&params_lerp)
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
            return Err("hit max evals");
            // return None;
        }
        let points = params_and_ress_vec.get_params().avg();
        let fit_residue = deconvolution_data.calc_residue_function(&points);
        fit_residue_evals += 1;
        Ok(FitResults {
            params: points,
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



pub fn exp(x: float) -> float {
    x.exp()
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
    // fn separate_chunks_from_start(&self, delimiter: impl ToString, chunks_size: usize) -> String;
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

pub trait ToStringUnderscoreSeparated {
    fn to_string_underscore_separated(&self) -> String;
}

impl ToStringUnderscoreSeparated for u64 {
    fn to_string_underscore_separated(&self) -> String {
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
        use crate::ExtStringSeparateChunks;
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

    mod index_of {
        mod min {
            mod with_ceil {
                mod without_nan {
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
                    use crate::IndexOfMinMaxWithCeilFloor;
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
        mod per_point {
            use crate::{DiffFunctionType, convolve_per_point, float};
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = points_spectrum_original.clone();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = [vec![0.; 9], vec![0.5, 1., 0.5], vec![0.; 9]].concat();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = [vec![1., 0.5], vec![0.; 19]].concat();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = [vec![0.; 19], vec![0.5, 1.]].concat();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = [vec![0.; 5], vec![0.5, 1., 0.5], vec![0.; 4], vec![0.5, 1., 0.5], vec![0.; 5]].concat();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = [vec![1., 0.5], vec![0.; 4], vec![0.5, 1., 0.5], vec![0.; 11]].concat();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                        let points_convolved_expected = [vec![0.; 11], vec![0.5, 1., 0.5], vec![0.; 4], vec![0.5, 1.]].concat();
                        let points_convolved_actual = convolve_per_point(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
        mod exponents {
            // TODO:
            // use crate::{DiffFunctionType, ExponentFunction, convolve_exponents, float};
            // mod instrument_is_identity {
            //     use super::*;
            //     const POINTS_INSTRUMENT: [float; 1] = [1.];
            //     #[test]
            //     fn original_spectrum_is_one_exp_1_1_1() {
            //         const EPSILON: float = 1e-6;
            //         const RES_LEN: usize = 10;
            //         println!("POINTS_INSTRUMENT = {:?}", POINTS_INSTRUMENT);
            //         let exponents = vec![ExponentFunction::new(1., 1., 1.)];
            //         let points_convolved_expected = [vec![0.; RES_LEN]].concat();
            //         let points_convolved_actual = convolve_exponents(&POINTS_INSTRUMENT.to_vec(), &exponents, RES_LEN);
            //         println!("points_convolved_expected = {:?}", points_convolved_expected);
            //         println!("points_convolved_actual = {:?}", points_convolved_actual);
            //         let diff = DiffFunctionType::DySqrPerEl.calc_diff(&points_convolved_expected, &points_convolved_actual);
            //         println!("diff = {}", diff);
            //         assert!(diff < 0.1);
            //         assert!(diff < EPSILON, "expected `diff < EPSILON, but diff={} and EPSILON={}`", diff, EPSILON);
            //     }
            // }
        }
    }

    mod deconvolve {
        mod per_point {
            use crate::{Deconvolution, DeconvolutionData, DeconvolutionResultOrError, DiffFunctionType, FitAlgorithmType, float};
            const FIT_ALGORITHM_TYPE: FitAlgorithmType = FitAlgorithmType::PatternSearch;
            const DECONVOLUTION: Deconvolution = Deconvolution::PerPoint {
                diff_function_type: DiffFunctionType::DySqr,
                antispikes: None,
            };
            fn deconvolve(points_instrument: Vec<float>, points_spectrum: Vec<float>) -> DeconvolutionResultOrError {
                let deconvolution_data: DeconvolutionData = DeconvolutionData::new(
                    points_instrument,
                    points_spectrum,
                    DECONVOLUTION
                );
                deconvolution_data.deconvolve(FIT_ALGORITHM_TYPE)
            }
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
                        let deconvolve_results = deconvolve(POINTS_INSTRUMENT.to_vec(), points_spectrum_convolved).unwrap();
                        let points_deconvolved_actual = deconvolve_results.params;
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
}

