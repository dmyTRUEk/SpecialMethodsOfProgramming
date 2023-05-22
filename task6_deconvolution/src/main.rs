//! Deconvolution.

use std::{
    cmp::Ordering,
    env,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
};

use rand::{Rng, rngs::ThreadRng, thread_rng};
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

mod aliases_method_to_function;
mod extensions;
mod float_type;
mod utils_io;

use aliases_method_to_function::exp;
use extensions::{Avg, ArrayMath, IndexOfMax, IndexOfMinWithCeil, ToStringUnderscoreSeparated, ToStringWithSignificantDigits};
use float_type::float;
use utils_io::{flush, press_enter_to_continue};


mod deconvolution_params {
    use super::{Deconvolution, DiffFunctionType, float};

    pub const EXPONENTS_AMOUNT: usize = 2; // only for `Deconvolution::Exponents`.

    pub const DECONVOLUTION: Deconvolution = {

        // Deconvolution::PerPoint {
        //     diff_function_type: DiffFunctionType::DySqr,
        //     // antispikes: None,
        //     antispikes: Some(Antispikes {
        //         antispikes_type: AntispikesType::DySqr,
        //         antispikes_k: 1.,
        //     }),
        //     initial_value: 0.,
        // }

        // Deconvolution::Exponents {
        //     diff_function_type: DiffFunctionType::DySqr,
        //     // exponents_amount: 2,
        //     initial_values: [
        //         // 30., -10., 30.,
        //         // 30., 1., -2.,
        //         1., 1., 1.,
        //         1., 1., 1.,
        //     ],
        // }

        // Deconvolution::SatExp_DecExp {
        //     // WARNING?: set `FIT_ALGORITHM_MIN_STEP` to `1e-3`.
        //     diff_function_type: DiffFunctionType::DySqr,
        //     //               a   s   t1  t2
        //     // initial_values: [1., 0., 1., 1.],
        //     // initial_values: [0.04, -12., 1., 30.], // fr: 1.870
        //     // initial_values: [1., 1., 1., 1.],
        //     initial_values: [0.1, -10., 1., 10.],
        // }

        // Deconvolution::Two_SatExp_DecExp {
        //     // WARNING?: set `FIT_ALGORITHM_MIN_STEP` to `1e-3`.
        //     diff_function_type: DiffFunctionType::DySqr,
        //     initial_values: [
        //         100., -10., 100., 10.,
        //         100., -10., 100., 10.,
        //     ],
        // }

        // Deconvolution::SatExp_DecExpPlusConst {
        //     diff_function_type: DiffFunctionType::DySqr,
        //     initial_values: [0.1, -1., 1e-2, 0.1, 10.],
        //     allow_tb_less_than_ta: false,
        // }

        // Deconvolution::SatExp_TwoDecExp {
        //     diff_function_type: DiffFunctionType::DySqr,
        //     // initial_values: [0.1, -10., 0.01, 1., 1.],
        //     initial_values: [0.02, -9., 6e-6, 35., 8.],
        // }

        // Deconvolution::SatExp_TwoDecExpPlusConst {
        //     diff_function_type: DiffFunctionType::DySqr,
        //     // initial_values: [0.1, -5., 1e-2, 0.1, 10., 20.],    // ../data3/AIS3col_e650.dat
        //     initial_values: [0.01, -10., 0.1, 5e-3, 5., 20.],      // ../data3/AIS3fil_e650.dat
        //     // initial_values: [0.03, -8., 0.085, 5e-3, 6., 18.3], // ../data3/AIS3fil_e650.dat
        // }

        Deconvolution::SatExp_TwoDecExp_SeparateConsts {
            diff_function_type: DiffFunctionType::DySqr,
            // initial_values: [0.1, -10., 0.01, 1., 1.],
            initial_values: [0.02, 0.02, -9., 0.1, 35., 100.],
        }

    };

    pub const TRY_RANDOMIZED_INITIAL_VALUES: bool = false;
    pub const INITIAL_VALUES_RANDOM_SCALE: float = 10.;
    pub const CHANGE_SING_PROBABILITY: float = 0.05;
}

mod output_params {
    pub const SIGNIFICANT_DIGITS: usize = 4;
}

mod fit_params {
    use super::{FitAlgorithmType, float};
    // pub const INITIAL_VALUES: float = 0.0015;
    pub const FIT_ALGORITHM_TYPE: FitAlgorithmType = FitAlgorithmType::PatternSearch;
    // pub const FIT_RESIDUE_GOAL   : float = 1e-1; // for Pattern Search
    pub const FIT_ALGORITHM_MIN_STEP: float = 1e-4; // for Pattern Search & Downhill Simplex
    pub const FIT_RESIDUE_EVALS_MAX : u64 = 1_000_000;
    pub const FIT_RESIDUE_MAX_VALUE : float = 1e6;
}

mod pattern_search_params {
    use super::float;
    pub const INITIAL_STEP: float = 1.;
    pub const ALPHA: float = 1.1;        // step increase coefficient
    pub const BETA : float = 1. / ALPHA; // step decrease coefficient
}

mod downhill_simplex_params {
    use super::{DiffFunctionType, float};
    pub const INITIAL_SIMPLEX_SCALE: float = 0.815;
    pub const PARAMS_DIFF_TYPE: DiffFunctionType = DiffFunctionType::DySqr;
}


fn main() {
    let args: Vec<_> = env::args().collect();
    let (filepathstr_instrument, filepathstr_measured): (&str, &str) = match &args[..] {
        [_, filepathstr_instrument, filepathstr_measured] => (filepathstr_instrument, filepathstr_measured),
        [_, _] => panic!("Expected two filename, provided only one."),
        [_] => panic!("Filenames not provided."),
        [] => unreachable!("Unexpected CLI args number."),
        _ => panic!("Too many CLI args.")
    };

    print!("Loading instrumental spectrum  from `{}`...", filepathstr_instrument); flush();
    let instrument = Spectrum::load_from_file(filepathstr_instrument);
    println!(" done");

    print!("Loading spectrum to deconvolve from `{}`...", filepathstr_measured); flush();
    let measured = Spectrum::load_from_file(filepathstr_measured);
    println!(" done");

    // TODO: warning if points in instr more than in spectrum.
    // assert!(measured.points.len() > instrument.points.len());

    println!("FIT_ALGORITHM_TYPE    : {:#?}", fit_params::FIT_ALGORITHM_TYPE);
    println!("FIT_ALGORITHM_MIN_STEP: {:.2e}", fit_params::FIT_ALGORITHM_MIN_STEP);
    // if fit_params::FIT_ALGORITHM_TYPE == FitAlgorithmType::PatternSearch {
    //     println!("FIT_RESIDUE_GOAL     : {:.2e}", fit_params::FIT_RESIDUE_GOAL);
    // }
    println!("FIT_RESIDUE_EVALS_MAX : {}", fit_params::FIT_RESIDUE_EVALS_MAX.to_string_underscore_separated());

    let file_instrument = Path::new(filepathstr_instrument);
    let file_spectrum   = Path::new(filepathstr_measured);
    assert_eq!(
        file_instrument.parent().unwrap().canonicalize().unwrap().to_str().unwrap(),
        file_spectrum  .parent().unwrap().canonicalize().unwrap().to_str().unwrap()
    );
    let filepath_output = file_spectrum.with_file_name(format!(
        "results_{}_{}.dat",
        file_instrument.file_stem().unwrap().to_str().unwrap(),
        file_spectrum.file_stem().unwrap().to_str().unwrap()
    ));
    let filepathstr_output: &str = filepath_output.to_str().unwrap();
    let filepath_output_convolved = file_spectrum.with_file_name(format!(
        "results_{}_{}_convolved.dat",
        file_instrument.file_stem().unwrap().to_str().unwrap(),
        file_spectrum.file_stem().unwrap().to_str().unwrap()
    ));
    let filepathstr_output_convolved: &str = filepath_output_convolved.to_str().unwrap();

    let deconvolution = deconvolution_params::DECONVOLUTION;

    let deconvolution_data: DeconvolutionData = DeconvolutionData {
        instrument,
        measured,
        deconvolution,
    }.aligned_steps_to_smaller();


    println!();
    let fit_residue_with_initial_values = deconvolution_data.calc_residue_function(&deconvolution_data.get_initial_params());
    println!("fit_residue @ initial_values: {}", fit_residue_with_initial_values);
    println!();

    let deconvolve_results = deconvolution_data.deconvolve(fit_params::FIT_ALGORITHM_TYPE);
    match deconvolve_results {
        Err(err) => println!("ERROR: {}", err),
        Ok(ref deconvolution_results_unwrapped) => {
            output_results(
                &deconvolution_data,
                deconvolution_results_unwrapped,
                filepathstr_output,
                filepathstr_output_convolved,
            );
        }
    }
    if !deconvolution_params::TRY_RANDOMIZED_INITIAL_VALUES { return }

    println!();
    println!("------- NOW TRYING RANDOM INITIAL VALUES -------");
    println!();

    let mut rng = thread_rng();
    let mut best_fit_residue: float = if deconvolve_results.is_ok() { deconvolve_results.unwrap().fit_residue } else { float::MAX };
    let mut initial_values_tried: u64 = 0;
    loop {
        initial_values_tried += 1;
        let mut deconvolution_data = deconvolution_data.clone();
        fn randomize_array<const N: usize>(array: &mut [float; N], rng: &mut ThreadRng) {
            for i in 0..N {
                let is_change_sign: bool = rng.gen_bool(deconvolution_params::CHANGE_SING_PROBABILITY);
                let random_scale: float = rng.gen_range(
                    1./deconvolution_params::INITIAL_VALUES_RANDOM_SCALE ..= deconvolution_params::INITIAL_VALUES_RANDOM_SCALE
                );
                array[i] *= if is_change_sign { -1. } else { 1. } * random_scale;
            }
        }
        match deconvolution_data.deconvolution {
            Deconvolution::PerPoint { .. } => panic!("there is no need to try different initial params"),
            Deconvolution::Exponents { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::SatExp_DecExp { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::Two_SatExp_DecExp { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::SatExp_DecExpPlusConst { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::SatExp_TwoDecExp { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::SatExp_TwoDecExpPlusConst { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::SatExp_TwoDecExp_SeparateConsts { ref mut initial_values, .. } => randomize_array(initial_values, &mut rng),
            Deconvolution::Fourier {} => unimplemented!(),
        }
        let deconvolution_results = deconvolution_data.deconvolve(fit_params::FIT_ALGORITHM_TYPE);
        match deconvolution_results {
            Ok(deconvolution_results_unwrapped) if deconvolution_results_unwrapped.fit_residue < best_fit_residue => {
                best_fit_residue = deconvolution_results_unwrapped.fit_residue;
                println!("{}", "-".repeat(42));
                println!("initial_values_tried: {}", initial_values_tried);
                // let Deconvolution::Exponents { initial_values, .. } = deconvolution_data.deconvolution else { unreachable!() };
                // dbg!(initial_values);
                output_results(
                    &deconvolution_data,
                    &deconvolution_results_unwrapped,
                    filepathstr_output,
                    filepathstr_output_convolved,
                );
                println!("{}", "-".repeat(42));
            }
            _ => {
                println!(
                    "fit_residue: {}",
                    deconvolution_results.as_ref()
                        .map(|dr| format!("{:.4}", dr.fit_residue))
                        .unwrap_or_else(|err| format!("Error: {err}"))
                );
            },
        }
    }
}


fn output_results(
    deconvolution_data: &DeconvolutionData,
    deconvolution_results: &FitResults,
    filepathstr_output: &str,
    filepathstr_output_convolved: &str,
) {
    println!("deconvolution_results = {deconvolution_results:#?}");
    // println!("fit_residue_evals = {}", deconvolution_results.fit_residue_evals.to_string_underscore_separated());

    let params = &deconvolution_results.params;

    let desmos_function_str = deconvolution_data.deconvolution.to_desmos_function(&params);
    if let Ok(ref desmos_function_str) = desmos_function_str {
        println!("{}", desmos_function_str);
    }

    // let mut file_output = File::create(filepath_output).unwrap();
    // assert_eq!((1010..=1089).count(), deconvolve_results.points.len());
    // // TODO(refactor): `zip_exact`.
    // for (x, point) in (1010..=1089).zip(deconvolve_results.points) {
    //     writeln!(file_output, "{x}\t{p}", p=point).unwrap();
    // }
    let fit_residue_and_evals_msg = format!(
        "fit residue {fr:.3} achieved in {fre} fit residue function evals",
        fr=deconvolution_results.fit_residue,
        fre=deconvolution_results.fit_residue_evals,
    );
    // TODO(refactor):
    // - extract function name into separate method?
    // - extract common logic
    match &deconvolution_data.deconvolution {
        Deconvolution::PerPoint { .. } => {
            let sd_deconvolved = Spectrum {
                points: deconvolution_results.params.clone(),
                step: deconvolution_data.get_step(),
                x_start: deconvolution_data.measured.x_start,
            };
            sd_deconvolved.write_to_file(filepathstr_output);
            todo!("also print `fit_residue` and `fit_residue_evals`");
        }
        self_ @ Deconvolution::Exponents { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            for parts in deconvolution_results.params.chunks(3) {
                let (amplitude, shift, tau) = (parts[0], parts[1], parts[2]);
                writeln!(file_output, "amplitude={amplitude}").unwrap();
                writeln!(file_output, "shift={shift}").unwrap();
                writeln!(file_output, "tau={tau}").unwrap();
            }
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        self_ @ Deconvolution::SatExp_DecExp { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            let (amplitude, shift, tau_a, tau_b) = (params[0], params[1], params[2], params[3]);
            writeln!(file_output, "amplitude={amplitude}").unwrap();
            writeln!(file_output, "shift={shift}").unwrap();
            writeln!(file_output, "tau_a={tau_a}").unwrap();
            writeln!(file_output, "tau_b={tau_b}").unwrap();
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        self_ @ Deconvolution::Two_SatExp_DecExp { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            let (amplitude_1, shift_1, tau_a1, tau_b1) = (params[0], params[1], params[2], params[3]);
            let (amplitude_2, shift_2, tau_a2, tau_b2) = (params[4], params[5], params[6], params[7]);
            writeln!(file_output, "amplitude_1={amplitude_1}").unwrap();
            writeln!(file_output, "shift_1={shift_1}").unwrap();
            writeln!(file_output, "tau_a1={tau_a1}").unwrap();
            writeln!(file_output, "tau_b1={tau_b1}").unwrap();
            writeln!(file_output, "amplitude_2={amplitude_2}").unwrap();
            writeln!(file_output, "shift_2={shift_2}").unwrap();
            writeln!(file_output, "tau_a2={tau_a2}").unwrap();
            writeln!(file_output, "tau_b2={tau_b2}").unwrap();
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        self_ @ Deconvolution::SatExp_DecExpPlusConst { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            let (amplitude, shift, height, tau_a, tau_b) = (params[0], params[1], params[2], params[3], params[4]);
            writeln!(file_output, "amplitude={amplitude}").unwrap();
            writeln!(file_output, "shift={shift}").unwrap();
            writeln!(file_output, "height={height}").unwrap();
            writeln!(file_output, "tau_a={tau_a}").unwrap();
            writeln!(file_output, "tau_b={tau_b}").unwrap();
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        self_ @ Deconvolution::SatExp_TwoDecExp { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            let (amplitude, shift, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4]);
            writeln!(file_output, "amplitude={amplitude}").unwrap();
            writeln!(file_output, "shift={shift}").unwrap();
            writeln!(file_output, "tau_a={tau_a}").unwrap();
            writeln!(file_output, "tau_b={tau_b}").unwrap();
            writeln!(file_output, "tau_c={tau_c}").unwrap();
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        self_ @ Deconvolution::SatExp_TwoDecExpPlusConst { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            let (amplitude, shift, height, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
            writeln!(file_output, "amplitude={amplitude}").unwrap();
            writeln!(file_output, "shift={shift}").unwrap();
            writeln!(file_output, "height={height}").unwrap();
            writeln!(file_output, "tau_a={tau_a}").unwrap();
            writeln!(file_output, "tau_b={tau_b}").unwrap();
            writeln!(file_output, "tau_c={tau_c}").unwrap();
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        self_ @ Deconvolution::SatExp_TwoDecExp_SeparateConsts { .. } => {
            let mut file_output = File::create(filepathstr_output).unwrap();
            writeln!(file_output, "{name} params ({fit_residue_and_evals_msg}):", name=self_.get_name()).unwrap();
            let (amplitude_b, amplitude_c, shift, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
            writeln!(file_output, "amplitude_b={amplitude_b}").unwrap();
            writeln!(file_output, "amplitude_c={amplitude_c}").unwrap();
            writeln!(file_output, "shift={shift}").unwrap();
            writeln!(file_output, "tau_a={tau_a}").unwrap();
            writeln!(file_output, "tau_b={tau_b}").unwrap();
            writeln!(file_output, "tau_c={tau_c}").unwrap();
            if let Ok(desmos_function_str) = desmos_function_str {
                writeln!(file_output, "{desmos_function_str}").unwrap();
            }
        }
        Deconvolution::Fourier {} => unimplemented!(),
    }

    let convolved_points: Vec<float> = deconvolution_data.convolve_from_params(&deconvolution_results.params);
    let convolved = Spectrum {
        points: convolved_points,
        x_start: deconvolution_data.measured.x_start,
        step: deconvolution_data.measured.step,
    };
    convolved.write_to_file(filepathstr_output_convolved);
}


type DeconvolutionResultOrError = FitResultsOrError;

/// Deconvolution type and it's corresponding params.
#[allow(non_camel_case_types)]
#[derive(Debug, Clone, PartialEq)]
pub enum Deconvolution {
    /// aka Simple
    /// [y0, y1, y2, ...]
    PerPoint {
        diff_function_type: DiffFunctionType,
        antispikes: Option<Antispikes>,
        initial_value: float, // [y0, y1, y2, ...]
    },
    /// a1*exp(-(x-s1)/t1) + ...
    Exponents {
        diff_function_type: DiffFunctionType,
        // exponents_amount: usize,
        // initial_values: &'a [float],
        // initial_values: Vec<float>,
        // TODO(refactor): save initial params as struct with named fields, not stupid array.
        initial_values: [float; 3*deconvolution_params::EXPONENTS_AMOUNT], // ai, si, ti, ...
    },
    /// a * (1-exp(-(x-s)/ta)) * exp(-(x-s)/tb)
    SatExp_DecExp {
        diff_function_type: DiffFunctionType,
        initial_values: [float; 4], // a, s, ta, tb
        // initial_values: { a: float, s: float, tau_a: float, tau_b: float },
    },
    /// a1 * (1-exp(-(x-s1)/ta1)) * exp(-(x-s1)/tb1) + a2 * (1-exp(-(x-s2)/ta2)) * exp(-(x-s2)/tb2)
    Two_SatExp_DecExp {
        diff_function_type: DiffFunctionType,
        initial_values: [float; 8], // a1, s1, ta1, tb1, a2, s2, ta2, tb2
    },
    /// a * (1-exp(-(x-s)/ta)) * (exp(-(x-s)/tb) + h)
    SatExp_DecExpPlusConst {
        diff_function_type: DiffFunctionType,
        // initial_values: [float; 4], // s, h, ta, tb
        initial_values: [float; 5], // a, s, h, ta, tb
        allow_tb_less_than_ta: bool,
    },
    /// a * (1-exp(-(x-s)/ta)) * (exp(-(x-s)/tb) + exp(-(x-s)/tc))
    SatExp_TwoDecExp {
        diff_function_type: DiffFunctionType,
        initial_values: [float; 5], // a, s, ta, tb, tc
    },
    /// a * (1-exp(-(x-s)/ta)) * (exp(-(x-s)/tb) + exp(-(x-s)/tc) + h)
    SatExp_TwoDecExpPlusConst {
        diff_function_type: DiffFunctionType,
        initial_values: [float; 6], // a, s, h, ta, tb, tc
    },
    /// (1-exp(-(x-s)/ta)) * (b*exp(-(x-s)/tb) + c*exp(-(x-s)/tc))
    SatExp_TwoDecExp_SeparateConsts {
        diff_function_type: DiffFunctionType,
        initial_values: [float; 6], // b, c, s, ta, tb, tc
    },
    Fourier {
        // unimplemented
    },
}
impl<'a> Deconvolution {
    pub const SAT_DEC_EXP_AMPLITUDE_SCALE: float = 1. / 1.;

    pub fn get_name(&self) -> &'static str {
        match self {
            Deconvolution::PerPoint { .. } => todo!(),
            Deconvolution::Exponents { .. } => "exponents",
            Deconvolution::SatExp_DecExp { .. } => "saturated decaying exponential",
            Deconvolution::Two_SatExp_DecExp { .. } => "two saturated decaying exponentials",
            Deconvolution::SatExp_DecExpPlusConst { .. } => "saturated decaying exponential plus const",
            Deconvolution::SatExp_TwoDecExp { .. } => "saturated exponential and two decaying exponentials",
            Deconvolution::SatExp_TwoDecExpPlusConst { .. } => "saturated exponential and two decaying exponentials plus const",
            Deconvolution::SatExp_TwoDecExp_SeparateConsts { .. } => "saturated exponential and two decaying exponentials with individual amplitudes",
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

    pub fn params_to_points(
        &self,
        params: &Vec<float>,
        points_len: usize,
        (x_start, x_end): (float, float),
    ) -> Vec<float> {
        assert!(points_len > 1);
        assert!(x_start < x_end);

        // TODO(optimization): measure perf, if this is slowing the program
        fn i_to_x(
            i: usize,
            points_len: usize,
            (x_start, x_end): (float, float),
        ) -> float {
            let x_range: float = x_end - x_start;
            let t: float = (i as float) / ((points_len - 1) as float);
            let x: float = t * x_range + x_start;
            x
        }

        match self {
            Self::PerPoint { .. } => params.to_vec(),
            Self::Exponents { .. } => {
                assert_eq!(deconvolution_params::EXPONENTS_AMOUNT * 3, params.len());
                let exponents: Vec<ExponentFunction> = params
                    .chunks(3).into_iter()
                    .map(|parts| ExponentFunction { amplitude: parts[0], shift: parts[1], tau: parts[2] } )
                    .collect();
                let mut points = vec![0.; points_len];
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let sum: float = exponents.iter()
                        .map(|exponent| exponent.eval_at(x))
                        .sum();
                    points[i] = sum;
                }
                points
            }
            Self::SatExp_DecExp { .. } => {
                let mut points = vec![0.; points_len];
                let (amplitude, shift, tau_a, tau_b) = (params[0], params[1], params[2], params[3]);
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let x_m_shift: float = x - shift;
                    // let y = if x >= shift
                    let y = if x_m_shift >= 0. {
                        // "optimization" (don't work): precalc `1/tau_a` & `1/tau_b`.
                        Self::SAT_DEC_EXP_AMPLITUDE_SCALE * amplitude * (1. - exp(-x_m_shift/tau_a)) * exp(-x_m_shift/tau_b)
                    } else {
                        0.
                    };
                    points[i] = y;
                }
                points
            }
            Self::Two_SatExp_DecExp { .. } => {
                let mut points = vec![0.; points_len];
                let (amplitude_1, shift_1, tau_a1, tau_b1) = (params[0], params[1], params[2], params[3]);
                let (amplitude_2, shift_2, tau_a2, tau_b2) = (params[4], params[5], params[6], params[7]);
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let x_m_shift_1: float = x - shift_1;
                    let x_m_shift_2: float = x - shift_2;
                    let y1 = if x_m_shift_1 >= 0. {
                        Self::SAT_DEC_EXP_AMPLITUDE_SCALE * amplitude_1 * (1. - exp(-(x_m_shift_1)/tau_a1)) * exp(-(x_m_shift_1)/tau_b1)
                    } else {
                        0.
                    };
                    let y2 = if x_m_shift_2 >= 0. {
                        Self::SAT_DEC_EXP_AMPLITUDE_SCALE * amplitude_2 * (1. - exp(-(x_m_shift_2)/tau_a2)) * exp(-(x_m_shift_2)/tau_b2)
                    } else {
                        0.
                    };
                    points[i] = y1 + y2;
                }
                points
            }
            Self::SatExp_DecExpPlusConst { .. } => {
                let mut points = vec![0.; points_len];
                let (amplitude, shift, height, tau_a, tau_b) = (params[0], params[1], params[2], params[3], params[4]);
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let x_m_shift: float = x - shift;
                    let y = if x_m_shift >= 0. {
                        amplitude * (1. - exp(-x_m_shift/tau_a)) * (exp(-x_m_shift/tau_b) + height)
                    } else {
                        0.
                    };
                    points[i] = y;
                }
                points
            }
            Self::SatExp_TwoDecExp { .. } => {
                let mut points = vec![0.; points_len];
                let (amplitude, shift, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4]);
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let x_m_shift: float = x - shift;
                    let y = if x_m_shift >= 0. {
                        amplitude * (1. - exp(-x_m_shift/tau_a)) * (exp(-x_m_shift/tau_b) + exp(-x_m_shift/tau_c))
                    } else {
                        0.
                    };
                    points[i] = y;
                }
                points
            }
            Self::SatExp_TwoDecExpPlusConst { .. } => {
                let mut points = vec![0.; points_len];
                let (amplitude, shift, height, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let x_m_shift: float = x - shift;
                    let y = if x_m_shift >= 0. {
                        amplitude * (1. - exp(-x_m_shift/tau_a)) * (exp(-x_m_shift/tau_b) + exp(-x_m_shift/tau_c) + height)
                    } else {
                        0.
                    };
                    points[i] = y;
                }
                points
            }
            Self::SatExp_TwoDecExp_SeparateConsts { .. } => {
                let mut points = vec![0.; points_len];
                let (b, c, shift, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
                for i in 0..points_len {
                    let x: float = i_to_x(i, points_len, (x_start, x_end));
                    let x_m_shift: float = x - shift;
                    let y = if x_m_shift >= 0. {
                        (1. - exp(-x_m_shift/tau_a)) * (b*exp(-x_m_shift/tau_b) + c*exp(-x_m_shift/tau_c))
                    } else {
                        0.
                    };
                    points[i] = y;
                }
                points
            }
            Self::Fourier {} => unimplemented!(),
        }
    }

    pub fn to_desmos_function(&self, params: &Vec<float>) -> Result<String, &'static str> {
        use output_params::SIGNIFICANT_DIGITS as SD;
        match self {
            Deconvolution::PerPoint { .. } => Err("impossible to build this function"),
            Deconvolution::Exponents { .. } => {
                Ok([
                    r"f_{",
                    &deconvolution_params::EXPONENTS_AMOUNT.to_string(),
                    r"}\left(x\right)=",
                    &params
                        .chunks(3).into_iter()
                        .map(|parts| {
                            let (amplitude, shift, tau) = (parts[0], parts[1], parts[2]);
                            let neg_shift = -shift;
                            assert_ne!(0., tau);
                            [
                                r"\left\{x",
                                if tau > 0. { ">" } else { "<" },
                                &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                                r":\frac{",
                                &format!("{}", amplitude.to_string_with_significant_digits(SD)),
                                r"}{",
                                &format!("{}", (1./ExponentFunction::AMPLITUDE_SCALE).to_string_with_significant_digits(SD)),
                                r"}e^{-\frac{x",
                                if neg_shift.is_sign_positive() { "+" } else { "" },
                                &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                                r"}{",
                                &format!("{}", tau.to_string_with_significant_digits(SD)),
                                r"}},0\right\}",
                            ].concat()
                        })
                        .reduce(|acc, el| format!("{acc}+{el}")).unwrap(),
                ].concat())
            }
            Deconvolution::SatExp_DecExp { .. } => {
                let (amplitude, shift, tau_a, tau_b) = (params[0], params[1], params[2], params[3]);
                let neg_shift = -shift;
                Ok([
                    // y=a\left(1-e^{-\frac{x-s}{t_{1}}}\right)\left(e^{-\frac{x-s-d}{t_{2}}}\right)\left\{x>s\right\}
                    r"y=\frac{",
                    &format!("{}", amplitude.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", (1./Deconvolution::SAT_DEC_EXP_AMPLITUDE_SCALE).to_string_with_significant_digits(SD)),
                    r"}\left(1-e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b.to_string_with_significant_digits(SD)),
                    r"}}\right)\left\{x\ge",
                    &format!("{}", shift.to_string_with_significant_digits(SD)),
                    r"\right\}",
                ].concat())
            }
            Deconvolution::Two_SatExp_DecExp { .. } => {
                let (amplitude_1, shift_1, tau_a1, tau_b1) = (params[0], params[1], params[2], params[3]);
                let (amplitude_2, shift_2, tau_a2, tau_b2) = (params[4], params[5], params[6], params[7]);
                let neg_shift_1 = -shift_1;
                let neg_shift_2 = -shift_2;
                Ok([
                    r"y=\frac{",
                    &format!("{}", amplitude_1.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", (1./Deconvolution::SAT_DEC_EXP_AMPLITUDE_SCALE).to_string_with_significant_digits(SD)),
                    r"}\left(1-e^{-\frac{x",
                    if neg_shift_1.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift_1.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a1.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(e^{-\frac{x",
                    if neg_shift_1.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift_1.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b1.to_string_with_significant_digits(SD)),
                    r"}}\right)\left\{x\ge",
                    &format!("{}", shift_1.to_string_with_significant_digits(SD)),
                    r"\right\}+\frac{",
                    &format!("{}", amplitude_2.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", (1./Deconvolution::SAT_DEC_EXP_AMPLITUDE_SCALE).to_string_with_significant_digits(SD)),
                    r"}\left(1-e^{-\frac{x",
                    if neg_shift_2.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift_2.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a2.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(e^{-\frac{x",
                    if neg_shift_2.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift_2.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b2.to_string_with_significant_digits(SD)),
                    r"}}\right)\left\{x\ge",
                    &format!("{}", shift_2.to_string_with_significant_digits(SD)),
                    r"\right\}",
                ].concat())
            }
            Deconvolution::SatExp_DecExpPlusConst { .. } => {
                let (amplitude, shift, height, tau_a, tau_b) = (params[0], params[1], params[2], params[3], params[4]);
                let neg_shift = -shift;
                Ok([
                    r"y=",
                    &format!("{}", amplitude.to_string_with_significant_digits(SD)),
                    r"\left(1-e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b.to_string_with_significant_digits(SD)),
                    r"}}+",
                    &format!("{}", height.to_string_with_significant_digits(SD)),
                    r"\right)\left\{x\ge",
                    &format!("{}", shift.to_string_with_significant_digits(SD)),
                    r"\right\}",
                ].concat())
            }
            Deconvolution::SatExp_TwoDecExp { .. } => {
                let (amplitude, shift, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4]);
                let neg_shift = -shift;
                Ok([
                    r"y=",
                    &format!("{}", amplitude.to_string_with_significant_digits(SD)),
                    r"\left(1-e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b.to_string_with_significant_digits(SD)),
                    r"}}+e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_c.to_string_with_significant_digits(SD)),
                    r"}}\right)\left\{x\ge",
                    &format!("{}", shift.to_string_with_significant_digits(SD)),
                    r"\right\}",
                ].concat())
            }
            Deconvolution::SatExp_TwoDecExpPlusConst { .. } => {
                let (amplitude, shift, height, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
                let neg_shift = -shift;
                Ok([
                    r"y=",
                    &format!("{}", amplitude.to_string_with_significant_digits(SD)),
                    r"\left(1-e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b.to_string_with_significant_digits(SD)),
                    r"}}+e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_c.to_string_with_significant_digits(SD)),
                    r"}}",
                    if height.is_sign_positive() { "+" } else { "" },
                    &format!("{}", height.to_string_with_significant_digits(SD)),
                    r"\right)\left\{x\ge",
                    &format!("{}", shift.to_string_with_significant_digits(SD)),
                    r"\right\}",
                ].concat())
            }
            Deconvolution::SatExp_TwoDecExp_SeparateConsts { .. } => {
                let (b, c, shift, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
                let neg_shift = -shift;
                Ok([
                    r"y=\left(1-",
                    // if a.is_sign_positive() { "-" } else { "+" },
                    // &format!("{}", a.abs().to_string_with_significant_digits(SD)),
                    r"e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_a.to_string_with_significant_digits(SD)),
                    r"}}\right)\left(",
                    &format!("{}", b.to_string_with_significant_digits(SD)),
                    r"e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_b.to_string_with_significant_digits(SD)),
                    r"}}",
                    if c.is_sign_positive() { "+" } else { "" },
                    &format!("{}", c.to_string_with_significant_digits(SD)),
                    r"e^{-\frac{x",
                    if neg_shift.is_sign_positive() { "+" } else { "" },
                    &format!("{}", neg_shift.to_string_with_significant_digits(SD)),
                    r"}{",
                    &format!("{}", tau_c.to_string_with_significant_digits(SD)),
                    r"}}\right)\left\{x\ge",
                    &format!("{}", shift.to_string_with_significant_digits(SD)),
                    r"\right\}",
                ].concat())
            }
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }
}


#[derive(Debug, Clone, PartialEq)]
pub struct DeconvolutionData {
    pub instrument: Spectrum,
    pub measured: Spectrum,
    pub deconvolution: Deconvolution,
}
impl DeconvolutionData {
    pub fn assert_steps_is_aligned(&self) {
        assert_eq!(self.instrument.step, self.measured.step);
    }
    pub fn assert_x_starts_is_aligned(&self) {
        assert_eq!(self.instrument.x_start, self.measured.x_start);
    }

    pub fn get_step(&self) -> float {
        self.assert_steps_is_aligned();
        self.instrument.step
    }

    pub fn get_x_start(&self) -> float {
        self.assert_x_starts_is_aligned();
        self.instrument.x_start
    }

    /// Make [`step`] in [`instrument`] and [`measured`] same,
    /// towards smaller step (more points in total).
    ///
    /// [`step`]: SpectrumData::step
    /// [`instrument`]: DeconvolutionData::instrument
    /// [`measured`]: DeconvolutionData::measured
    pub fn aligned_steps_to_smaller(mut self) -> Self {
        match self.instrument.step.partial_cmp(&self.measured.step) {
            Some(Ordering::Equal) => return self,
            Some(Ordering::Less) => {
                self.measured = self.measured.recalculated_with_step(self.instrument.step);
            }
            Some(Ordering::Greater) => {
                self.instrument = self.instrument.recalculated_with_step(self.measured.step);
            }
            None => panic!("One of the steps is `NaN`")
        };
        self.assert_steps_is_aligned();
        self
    }

    /// Make [`step`] in [`instrument`] and [`measured`] same,
    /// towards bigger step (less points in total).
    ///
    /// [`step`]: SpectrumData::step
    /// [`instrument`]: DeconvolutionData::instrument
    /// [`measured`]: DeconvolutionData::measured
    pub fn aligned_steps_to_bigger(mut self) -> Self {
        match self.instrument.step.partial_cmp(&self.measured.step) {
            Some(Ordering::Equal) => return self,
            Some(Ordering::Less) => {
                self.instrument = self.instrument.recalculated_with_step(self.measured.step);
            }
            Some(Ordering::Greater) => {
                self.measured = self.measured.recalculated_with_step(self.instrument.step);
            }
            None => panic!("One of the steps is `NaN`")
        };
        self.assert_steps_is_aligned();
        self
    }

    // TODO: align_steps_to_bigger, align_steps_to_smaller

    pub fn deconvolve(&self, fit_algorithm_type: FitAlgorithmType) -> DeconvolutionResultOrError {
        self.assert_steps_is_aligned();
        fit_algorithm_type.fit(&self)
    }

    /// depending on the `self.deconvolution` `params` is:
    /// - PerPoint: list of values at that point
    /// - Exponents: list of (amplitude, shift, tau)
    /// - Fourier: unimplemented
    pub fn calc_residue_function(&self, params: &Vec<float>) -> float {
        let points_convolved: Vec<float> = self.convolve_from_params(params);
        assert_eq!(self.get_params_amount(), params.len());
        match &self.deconvolution {
            Deconvolution::PerPoint { diff_function_type, antispikes, .. } => {
                diff_function_type.calc_diff_with_antispikes(&self.measured.points, &points_convolved, antispikes)
            }
            Deconvolution::Exponents { diff_function_type, .. }
            | Deconvolution::SatExp_DecExp { diff_function_type, .. }
            | Deconvolution::Two_SatExp_DecExp { diff_function_type, .. }
            | Deconvolution::SatExp_DecExpPlusConst { diff_function_type, .. }
            | Deconvolution::SatExp_TwoDecExp { diff_function_type, .. }
            | Deconvolution::SatExp_TwoDecExpPlusConst { diff_function_type, .. }
            | Deconvolution::SatExp_TwoDecExp_SeparateConsts { diff_function_type, .. }
            => {
                diff_function_type.calc_diff(&self.measured.points, &points_convolved)
            }
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

    pub fn get_params_amount(&self) -> usize {
        match self.deconvolution {
            Deconvolution::PerPoint { .. } => self.measured.points.len(),
            Deconvolution::Exponents { initial_values, .. } => initial_values.len(),
            Deconvolution::SatExp_DecExp { initial_values, .. } => initial_values.len(),
            Deconvolution::Two_SatExp_DecExp { initial_values, .. } => initial_values.len(),
            Deconvolution::SatExp_DecExpPlusConst { initial_values, .. } => initial_values.len(),
            Deconvolution::SatExp_TwoDecExp { initial_values, .. } => initial_values.len(),
            Deconvolution::SatExp_TwoDecExpPlusConst { initial_values, .. } => initial_values.len(),
            Deconvolution::SatExp_TwoDecExp_SeparateConsts { initial_values, .. } => initial_values.len(),
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

    pub fn get_initial_params(&self) -> Vec<float> {
        let initial_params: Vec<float> = match &self.deconvolution {
            Deconvolution::PerPoint { initial_value, .. } => vec![*initial_value; self.get_params_amount()],
            Deconvolution::Exponents { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::SatExp_DecExp { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::Two_SatExp_DecExp { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::SatExp_DecExpPlusConst { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::SatExp_TwoDecExp { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::SatExp_TwoDecExpPlusConst { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::SatExp_TwoDecExp_SeparateConsts { initial_values, .. } => initial_values.to_vec(),
            Deconvolution::Fourier {} => unimplemented!(),
        };
        assert_eq!(self.get_params_amount(), initial_params.len());
        initial_params
    }

    pub fn is_params_ok(&self, params: &Vec<float>) -> bool {
        match self.deconvolution {
            Deconvolution::PerPoint { .. } => params.into_iter().all(|&x| x >= 0.),
            // Deconvolution::Exponents { .. } => params.into_iter().enumerate().all(|(i, &x)| match i % 3 {
            //     0 => x >= 0.,
            //     1 => true,
            //     2 => true,
            //     _ => unreachable!()
            // }),
            Deconvolution::Exponents { .. } => params.chunks(3).into_iter().all(|parts| {
                let (amplitude, _tau, _shift) = (parts[0], parts[1], parts[2]);
                amplitude >= 0.
            }),
            Deconvolution::SatExp_DecExp { .. } => {
                let (amplitude, _, tau_a, tau_b) = (params[0], params[1], params[2], params[3]);
                amplitude >= 0. && tau_a >= 0. && tau_b >= 0.
            }
            Deconvolution::Two_SatExp_DecExp { .. } => {
                let (amplitude_1, _, tau_a1, tau_b1) = (params[0], params[1], params[2], params[3]);
                let (amplitude_2, _, tau_a2, tau_b2) = (params[4], params[5], params[6], params[7]);
                amplitude_1 >= 0. && tau_a1 >= 0. && tau_b1 >= 0. &&
                amplitude_2 >= 0. && tau_a2 >= 0. && tau_b2 >= 0.
            }
            Deconvolution::SatExp_DecExpPlusConst { allow_tb_less_than_ta, .. } => {
                let (amplitude, _, height, tau_a, tau_b) = (params[0], params[1], params[2], params[3], params[4]);
                amplitude >= 0. && height >= 0. && tau_a >= 0. && tau_b >= 0. && if allow_tb_less_than_ta { true } else { tau_a < tau_b }
            }
            Deconvolution::SatExp_TwoDecExp { .. } => {
                let (amplitude, _, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4]);
                amplitude >= 0. && tau_a >= 0. && tau_b >= 0. && tau_c >= 0.
            }
            Deconvolution::SatExp_TwoDecExpPlusConst { .. } => {
                let (amplitude, _, height, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
                amplitude >= 0. && height >= 0. && tau_a >= 0. && tau_b >= 0. && tau_c >= 0.
            }
            Deconvolution::SatExp_TwoDecExp_SeparateConsts { .. } => {
                let (b, c, _, tau_a, tau_b, tau_c) = (params[0], params[1], params[2], params[3], params[4], params[5]);
                b >= 0. && c >= 0. && tau_a >= 0. && tau_b >= 0. && tau_c >= 0.
            }
            Deconvolution::Fourier {} => unimplemented!(),
        }
    }

    pub fn convolve_from_params(&self, params: &Vec<float>) -> Vec<float> {
        // convert `params` into `points` ("deconvolved"):
        let points_deconvolved: Vec<float> = self.deconvolution.params_to_points(
            &params,
            self.measured.points.len(),
            (self.measured.x_start, self.measured.get_x_end())
        );
        self.convolve_from_points(&points_deconvolved)
    }

    pub fn convolve_from_points(&self, points_deconvolved: &Vec<float>) -> Vec<float> {
        let points_convolved: Vec<float> = convolve_by_points(&self.instrument.points, &points_deconvolved);
        assert_eq!(self.measured.points.len(), points_convolved.len());
        points_convolved
    }
}


/// Must be used only in `tests` & `DeconvolutionData::convolve()`.
pub fn convolve_by_points(points_instrument: &Vec<float>, points_deconvolved: &Vec<float>) -> Vec<float> {
    let points_instrument_len: usize = points_instrument.len();
    let points_deconvolved_len: usize = points_deconvolved.len();
    assert!(points_instrument_len % 2 == 1, "points_instrument_len = {}", points_instrument_len); // why?
    let mut points_convolved = vec![0.; points_deconvolved_len];
    for i in 0..points_deconvolved_len {
        let mut point_convolved = 0.;
        for j in 0..points_instrument_len {
            let d: i32 = j as i32 - points_instrument_len as i32 / 2;
            let pii: i32 = j as i32;     // points_instrument_index
            let psi: i32 = i as i32 - d; // points_spectrum_index
            let is_pii_in_range: bool = 0 <= pii && pii < points_instrument_len as i32;
            let is_psi_in_range: bool = 0 <= psi && psi < points_deconvolved_len as i32;
            if is_pii_in_range && is_psi_in_range {
                let point_instrument  = points_instrument [pii as usize];
                let point_deconvolved = points_deconvolved[psi as usize];
                point_convolved += point_instrument * point_deconvolved;
            }
        }
        points_convolved[i] = point_convolved;
    }
    points_convolved
}


#[derive(Debug, Clone, PartialEq)]
pub struct Spectrum {
    pub points: Vec<float>,
    pub step: float,
    pub x_start: float,
}
impl Spectrum {
    pub fn get_x_range(&self) -> float {
        // self.get_x_end() - self.x_start
        self.step * (self.points.len().saturating_sub(1) as float)
    }

    pub fn get_x_end(&self) -> float {
        // self.start_x + self.step * (self.points.len() as float)
        self.get_x_from_index(self.points.len())
    }

    pub fn get_xy_from_index(&self, i: usize) -> (float, float) {
        (self.get_x_from_index(i), self.get_y_from_index(i))
    }

    pub fn get_x_from_index(&self, i: usize) -> float {
        self.x_start + self.step * (i as float)
    }

    pub fn get_y_from_index(&self, i: usize) -> float {
        self.points[i]
    }

    pub fn get_indices_of_closest_to_lhs_rhs(&self, x: float) -> (usize, usize) {
        assert!(self.x_start <= x && x <= self.get_x_end());
        let x_from_start: float = x - self.x_start;
        let index_as_float: float = x_from_start / self.step;
        (index_as_float.floor() as usize, index_as_float.ceil() as usize)
    }

    fn get_points_len_after_recalc_with_step(&self, step_new: float) -> usize {
        assert!(self.points.len() > 1);
        ((self.get_x_range()) / step_new).floor() as usize + 1
    }

    pub fn recalculated_with_step(mut self, step_new: float) -> Self {
        if !step_new.is_finite() { panic!() }
        let self_old = self.clone();
        self.points = vec![];
        self.step = step_new;
        let points_len_after_recalc: usize = self_old.get_points_len_after_recalc_with_step(step_new);
        for i in 0..points_len_after_recalc {
            let x: float = self.get_x_from_index(i);
            let (index_of_closest_lhs, index_of_closest_rhs): (usize, usize) = self_old.get_indices_of_closest_to_lhs_rhs(x);
            let y: float = if index_of_closest_lhs == index_of_closest_rhs {
                self_old.points[index_of_closest_lhs]
            } else {
                assert_eq!(1, index_of_closest_rhs - index_of_closest_lhs);
                let t: float = (x - self_old.get_x_from_index(index_of_closest_lhs)) / self_old.step;
                (1.-t) * self_old.points[index_of_closest_lhs] + t * self_old.points[index_of_closest_rhs]
            };
            self.points.push(y);
        }
        self
    }

    pub fn write_to_file(&self, filename: &str) {
        self.write_to_file_with_separators(filename, ".", "\t");
    }
    pub fn write_to_file_with_separators(&self, filepath: &str, decimal_point_symbol: &str, numbers_separator: &str) {
        let mut file_output = File::create(filepath).unwrap();
        for i in 0..self.points.len() {
            let (x, y) = self.get_xy_from_index(i);
            writeln!(
                file_output,
                "{x}{numbers_separator}{y}",
                x=format!("{x}").replace('.', decimal_point_symbol),
                y=format!("{y}").replace('.', decimal_point_symbol),
            ).unwrap();
        }
    }

    pub fn load_from_file(filename: &str) -> Self {
        Self::try_load_from_file(filename).unwrap()
    }
    pub fn try_load_from_file(filename: &str) -> Result<Self, &'static str> {
        let Ok(file) = File::open(filename) else { return Err("Unable to open file") };
        let lines = BufReader::new(file).lines();
        let mut x_start: Option<float> = None;
        let mut x_prev: Option<float> = None;
        let mut step: Option<float> = None;
        let mut ys = Vec::<float>::with_capacity(20);
        for line in lines.into_iter() {
            let Ok(line) = line else { return Err("Unable to unwrap line") };
            let line = line.trim();
            if line == "" { continue }
            let Some((x, y)) = line.split_once([' ', '\t']) else { return Err("Unable to split line once at space or tab.") };
            let Ok(x) = x.trim().replace(',', ".").parse::<float>() else { return Err("Unable to parse `x`") };
            match x_start {
                None => {
                    x_start = Some(x);
                }
                Some(x_start) => {
                    match step {
                        None => {
                            step = Some(x - x_start);
                        }
                        Some(step) => {
                            assert!(((x - x_prev.unwrap()) - step).abs() < 1e-6);
                        }
                    }
                    x_prev = Some(x);
                }
            }
            let Ok(y) = y.trim().replace(',', ".").parse() else { return Err("Unable to parse `y`") };
            ys.push(y);
        }
        let Some(x_start) = x_start else { return Err("`start_x` is None") };
        let Some(step) = step else { return Err("`step` is None") };
        Ok(Spectrum {
            points: ys,
            x_start,
            step,
        })
    }
}


#[derive(Debug, Clone, PartialEq)]
pub struct Antispikes {
    antispikes_type: AntispikesType,
    antispikes_k: float,
}
impl Antispikes {
    pub fn calc(&self, points_1: &Vec<float>, points_2: &Vec<float>) -> float {
        self.antispikes_k * self.antispikes_type.calc(points_1, points_2)
    }
}


#[derive(Debug, Clone, PartialEq)]
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
    pub const AMPLITUDE_SCALE: float = 0.001;

    // pub const fn new(amplitude: float, shift: float, tau: float) -> Self {
    //     Self { amplitude, shift, tau }
    // }

    pub fn eval_at(&self, x: float) -> float {
        Self::eval_at_(self.amplitude, self.tau, self.shift, x)
    }

    /// This is to prevent memory "segmentation":
    /// [`ExponentFunction`] have 3 floats, but whole struct will be aligned to 4 floats (i guess?)
    /// + 1 float as arg => 5 floats in memory,
    /// whereas this method uses only 4 floats, as expected.
    ///
    /// Also this maybe improves cache locality a tiny bit (no extra ghost float in memory).
    ///
    /// Unfortunately, no performance gain was measured.
    pub fn eval_at_(amplitude: float, tau: float, shift: float, x: float) -> float {
        // "optimization" (i think it won't work): somehow precalc `1/tau`.
        let in_exp = -(x - shift) / tau;
        if in_exp <= 0. {
            Self::AMPLITUDE_SCALE * amplitude * exp(in_exp)
        } else {
            0.
        }
    }
}


#[derive(Debug, Clone, PartialEq)]
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
        use crate::{fit_params::*, pattern_search_params::*};
        const DEBUG: bool = false;

        let f_params_amount: usize = deconvolution_data.get_params_amount();
        if f_params_amount == 0 {
            return Err("too few params");
            // return None;
        }

        type Params = Vec<float>;
        let mut params: Params = deconvolution_data.get_initial_params();
        let mut step: float = INITIAL_STEP;
        let mut fit_residue_evals: u64 = 0;

        let mut res_at_current_params: float = deconvolution_data.calc_residue_function(&params);
        fit_residue_evals += 1;
        if DEBUG { println!("res_at_current_params = {}", res_at_current_params) }
        if !res_at_current_params.is_finite() { return Err("`res_at_current_params` isn't finite") }
        // if !res_at_current_params.is_finite() { return None }
        if res_at_current_params >= FIT_RESIDUE_MAX_VALUE { return Err("`res_at_current_params` is too big") }

        while step > FIT_ALGORITHM_MIN_STEP && fit_residue_evals < FIT_RESIDUE_EVALS_MAX {
        // while residue_function(&params, &points_instrument, &points_spectrum) > FIT_RESIDUE_GOAL && fit_residue_evals < FIT_RESIDUE_EVALS_MAX
            if DEBUG {
                println!("params = {:#?}", params);
                println!("step = {}", step);
            }

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

                    res_at_current_params = ress_at_shifted_params[index_of_min];
                    if DEBUG { println!("res_at_current_params = {}", res_at_current_params) }
                    if !res_at_current_params.is_finite() { return Err("`res_at_current_params` isn't finite") }
                    // if !res_at_current_params.is_finite() { return None }
                    if res_at_current_params >= FIT_RESIDUE_MAX_VALUE { return Err("`res_at_current_params` is too big") }

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
        let fit_residue = res_at_current_params;
        Ok(FitResults {
            params,
            fit_residue,
            fit_residue_evals,
        })
    }


    #[allow(unused)]
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
        // unimplemented!("must rewrite using all new methods and fields");
        // let mut params_prev_prev: Params = vec![INITIAL_VALUES+INITIAL_SIMPLEX_SCALE; f_params_amount];
        // let mut params_prev_this: Params = vec![INITIAL_VALUES-INITIAL_SIMPLEX_SCALE; f_params_amount];
        let mut params_prev_prev: Params = deconvolution_data.get_initial_params().into_iter().map(|p| p + INITIAL_SIMPLEX_SCALE).collect::<Vec<float>>();
        let mut params_prev_this: Params = deconvolution_data.get_initial_params().into_iter().map(|p| p - INITIAL_SIMPLEX_SCALE).collect::<Vec<float>>();
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
        // params_and_ress_vec_push(vec![INITIAL_VALUES-INITIAL_SIMPLEX_SCALE/(f_params_amount as float); f_params_amount]);
        params_and_ress_vec_push(deconvolution_data.get_initial_params().into_iter().map(|p| p - INITIAL_SIMPLEX_SCALE/(f_params_amount as float)).collect::<Vec<float>>());
        for i in 0..f_params_amount {
            // let mut params = vec![INITIAL_VALUES; f_params_amount];
            let mut params = deconvolution_data.get_initial_params();
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





#[cfg(test)]
mod convolve {
    mod per_point {
        use crate::{DiffFunctionType, convolve_by_points, float};
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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
                    let points_convolved_actual = convolve_by_points(&POINTS_INSTRUMENT.to_vec(), &points_spectrum_original);
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

#[cfg(test)]
mod deconvolution_data {
    use crate::{Deconvolution, DeconvolutionData, DiffFunctionType, Spectrum};
    mod align_steps_to_smaller {
        use super::*;
        #[test]
        fn align_i0_4_to_m0_2() {
            assert_eq!(
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0., 0., 0.4999999999999999, 1., 0.5000000000000001, 0., 0., 0.],
                        step: 0.2,
                        x_start: 0.7,
                    },
                    measured: Spectrum {
                        points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                        step: 0.2,
                        x_start: 0.3,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                },
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0., 1., 0., 0.,],
                        step: 0.4,
                        x_start: 0.7,
                    },
                    measured: Spectrum {
                        points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                        step: 0.2,
                        x_start: 0.3,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                }.aligned_steps_to_smaller()
            );
        }
        #[test]
        fn align_m0_4_to_i0_2() {
            assert_eq!(
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                        step: 0.2,
                        x_start: 0.5,
                    },
                    measured: Spectrum {
                        points: vec![0., 0., 0., 0.4999999999999999, 1., 0.5000000000000007, 0., 0., 0.],
                        step: 0.2,
                        x_start: 0.9,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                },
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                        step: 0.2,
                        x_start: 0.5,
                    },
                    measured: Spectrum {
                        points: vec![0., 0., 1., 0., 0.,],
                        step: 0.4,
                        x_start: 0.9,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                }.aligned_steps_to_smaller()
            );
        }
    }
    mod align_steps_to_bigger {
        use super::*;
        #[test]
        fn align_m0_2_to_i0_4() {
            assert_eq!(
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0., 1., 0., 0.,],
                        step: 0.4,
                        x_start: 0.1,
                    },
                    measured: Spectrum {
                        points: vec![0., 0.2, 0.4, 0.2, 0.],
                        step: 0.4,
                        x_start: 0.5,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                },
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0., 1., 0., 0.,],
                        step: 0.4,
                        x_start: 0.1,
                    },
                    measured: Spectrum {
                        points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                        step: 0.2,
                        x_start: 0.5,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                }.aligned_steps_to_bigger()
            );
        }
        #[test]
        fn align_i0_2_to_m0_4() {
            assert_eq!(
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0.2, 0.4, 0.2, 0.],
                        step: 0.4,
                        x_start: 0.5,
                    },
                    measured: Spectrum {
                        points: vec![0., 0., 1., 0., 0.,],
                        step: 0.4,
                        x_start: 0.9,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                },
                DeconvolutionData {
                    instrument: Spectrum {
                        points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                        step: 0.2,
                        x_start: 0.5,
                    },
                    measured: Spectrum {
                        points: vec![0., 0., 1., 0., 0.,],
                        step: 0.4,
                        x_start: 0.9,
                    },
                    deconvolution: Deconvolution::PerPoint {
                        diff_function_type: DiffFunctionType::DySqr,
                        antispikes: None,
                        initial_value: 0.,
                    },
                }.aligned_steps_to_bigger()
            );
        }
    }
}

#[cfg(test)]
mod deconvolve {
    mod per_point {
        use crate::{Deconvolution, DeconvolutionData, DeconvolutionResultOrError, DiffFunctionType, FitAlgorithmType, Spectrum, float};
        const FIT_ALGORITHM_TYPE: FitAlgorithmType = FitAlgorithmType::PatternSearch;
        const DECONVOLUTION: Deconvolution = Deconvolution::PerPoint {
            diff_function_type: DiffFunctionType::DySqr,
            antispikes: None,
            initial_value: 0.,
        };
        fn deconvolve(points_instrument: Vec<float>, points_spectrum: Vec<float>) -> DeconvolutionResultOrError {
            let instrument: Spectrum = Spectrum {
                points: points_instrument,
                step: 1.,
                x_start: 0.,
            };
            let measured: Spectrum = Spectrum {
                points: points_spectrum,
                step: 1.,
                x_start: 0.,
            };
            let deconvolution_data: DeconvolutionData = DeconvolutionData {
                instrument,
                measured,
                deconvolution: DECONVOLUTION,
            };
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

#[cfg(test)]
mod spectrum {
    use crate::Spectrum;
    mod get_x_from_index {
        use super::*;
        #[test]
        fn _5_0() {
            assert_eq!(
                0.7,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_from_index(0)
            );
        }
        #[test]
        fn _5_1() {
            assert_eq!(
                1.1,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_from_index(1)
            );
        }
        #[test]
        fn _5_2() {
            assert_eq!(
                1.5,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_from_index(2)
            );
        }
        #[test]
        fn _5_3() {
            assert_eq!(
                1.9000000000000001,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_from_index(3)
            );
        }
        #[test]
        fn _5_4() {
            assert_eq!(
                2.3,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_from_index(4)
            );
        }
        #[test]
        fn _5_5() {
            assert_eq!(
                2.7,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_from_index(5)
            );
        }
    }
    mod get_x_end {
        use super::*;
        #[test]
        fn _0() {
            assert_eq!(
                0.7,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_end()
            );
        }
        #[test]
        fn _1() {
            assert_eq!(
                1.1,
                Spectrum {
                    points: vec![0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_end()
            );
        }
        #[test]
        fn _2() {
            assert_eq!(
                1.5,
                Spectrum {
                    points: vec![0., 0.1],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_end()
            );
        }
        #[test]
        fn _3() {
            assert_eq!(
                1.9000000000000001,
                Spectrum {
                    points: vec![0., 0.1, 0.2],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_end()
            );
        }
        #[test]
        fn _4() {
            assert_eq!(
                2.3,
                Spectrum {
                    points: vec![0., 0.1, 0.2, 0.1],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_end()
            );
        }
        #[test]
        fn _5() {
            assert_eq!(
                2.7,
                Spectrum {
                    points: vec![0., 0.1, 0.2, 0.1, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_end()
            );
        }
    }
    mod get_x_range {
        use super::*;
        #[test]
        fn _0() {
            assert_eq!(
                0.,
                Spectrum {
                    points: vec![],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_range()
            );
        }
        #[test]
        fn _1() {
            assert_eq!(
                0.,
                Spectrum {
                    points: vec![0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_range()
            );
        }
        #[test]
        fn _2() {
            assert_eq!(
                0.4,
                Spectrum {
                    points: vec![0., 0.1],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_range()
            );
        }
        #[test]
        fn _3() {
            assert_eq!(
                0.8,
                Spectrum {
                    points: vec![0., 0.1, 0.2],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_range()
            );
        }
        #[test]
        fn _4() {
            assert_eq!(
                1.2000000000000002,
                Spectrum {
                    points: vec![0., 0.1, 0.2, 0.1],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_range()
            );
        }
        #[test]
        fn _5() {
            assert_eq!(
                1.6,
                Spectrum {
                    points: vec![0., 0.1, 0.2, 0.1, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_x_range()
            );
        }
    }
    #[allow(non_snake_case)]
    mod get_indices_of_closest_to_lhs_rhs {
        use super::*;
        #[test]
        fn _5__0_1() {
            for x in [0.8, 0.9, 1.] {
                dbg!(x);
                assert_eq!(
                    (0, 1),
                    Spectrum {
                        points: vec![0., 0.1, 0.2, 0.1, 0.],
                        step: 0.4,
                        x_start: 0.7,
                    }.get_indices_of_closest_to_lhs_rhs(x)
                );
            }
        }
        #[test]
        fn _5__1_2() {
            for x in [1.2, 1.3, 1.4] {
                dbg!(x);
                assert_eq!(
                    (1, 2),
                    Spectrum {
                        points: vec![0., 0.1, 0.2, 0.1, 0.],
                        step: 0.4,
                        x_start: 0.7,
                    }.get_indices_of_closest_to_lhs_rhs(x)
                );
            }
        }
        #[test]
        fn _5__2_3() {
            for x in [1.6, 1.7, 1.8] {
                dbg!(x);
                assert_eq!(
                    (2, 3),
                    Spectrum {
                        points: vec![0., 0.1, 0.2, 0.1, 0.],
                        step: 0.4,
                        x_start: 0.7,
                    }.get_indices_of_closest_to_lhs_rhs(x)
                );
            }
        }
        #[test]
        fn _5__3_4() {
            for x in [2., 2.1, 2.2] {
                dbg!(x);
                assert_eq!(
                    (3, 4),
                    Spectrum {
                        points: vec![0., 0.1, 0.2, 0.1, 0.],
                        step: 0.4,
                        x_start: 0.7,
                    }.get_indices_of_closest_to_lhs_rhs(x)
                );
            }
        }
    }
    #[allow(non_snake_case)]
    mod get_points_len_after_recalc_with_step {
        #[test]
        fn _2__0_2() {
            assert_eq!(
                6, // dx: 0. 0.2 0.4 0.6 0.8 1.0
                Spectrum {
                    // dx:       0.  1.
                    points: vec![0., 10.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.2)
            );
        }
        #[test]
        fn _2__0_199() {
            assert_eq!(
                6, // dx: 0. 0.199 0.398 0.597 0.796 0.995
                Spectrum {
                    // dx:       0.  1.
                    points: vec![0., 10.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.199)
            );
        }
        #[test]
        fn _2__0_201() {
            assert_eq!(
                5, // dx: 0. 0.201 0.402 0.603 0.804
                Spectrum {
                    // dx:       0.  1.
                    points: vec![0., 10.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.201)
            );
        }
        #[test]
        fn _3__0_2() {
            assert_eq!(
                11, // dx: 0. 0.2 0.4 0.6 0.8 1. 1.2 1.4 1.6 1.8 2.0
                Spectrum {
                    // dx:       0.  1.   2.
                    points: vec![0., 10., 20.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.2)
            );
        }
        #[test]
        fn _3__0_199() {
            assert_eq!(
                11,
                Spectrum {
                    // dx:       0.  1.   2.
                    points: vec![0., 10., 20.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.199)
            );
        }
        #[test]
        fn _3__0_201() {
            assert_eq!(
                10,
                Spectrum {
                    // dx:       0.  1.   2.
                    points: vec![0., 10., 20.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.201)
            );
        }
        #[test]
        fn _4__0_2() {
            assert_eq!(
                16, // dx: 0. 0.2 0.4 0.6 0.8 1. 1.2 1.4 1.6 1.8 2. 2.2 2.4 2.6 2.8 3.0
                Spectrum {
                    // dx:       0.  1.   2.   3.
                    points: vec![0., 10., 20., 30.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.2)
            );
        }
        #[test]
        fn _4__0_199() {
            assert_eq!(
                16,
                Spectrum {
                    // dx:       0.  1.   2.   3.
                    points: vec![0., 10., 20., 30.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.199)
            );
        }
        #[test]
        fn _4__0_201() {
            assert_eq!(
                15,
                Spectrum {
                    // dx:       0.  1.   2.   3.
                    points: vec![0., 10., 20., 30.],
                    step: 1.,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.201)
            );
        }
        use super::*;
        #[test]
        fn _5__0_2() {
            assert_eq!(
                9,
                Spectrum {
                    points: vec![0., 0.2, 0.4, 0.2, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.2)
            );
        }
        #[test]
        fn _5__0_199() {
            assert_eq!(
                9,
                Spectrum {
                    points: vec![0., 0.2, 0.4, 0.2, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.199)
            );
        }
        #[test]
        fn _5__0_201() {
            assert_eq!(
                8,
                Spectrum {
                    points: vec![0., 0.2, 0.4, 0.2, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.get_points_len_after_recalc_with_step(0.201)
            );
        }
    }
    mod recalculate_with_step {
        use super::*;
        #[test]
        fn _5_into_9() {
            assert_eq!(
                Spectrum {
                    // points: vec![0., 0.1, 0.2, 0.3, 0.4, 0.3, 0.2, 0.1, 0.],
                    points: vec![0., 0.09999999999999998, 0.2, 0.3, 0.4, 0.30000000000000004, 0.2, 0.10000000000000003, 1.5543122344752193e-16],
                    step: 0.2,
                    x_start: 0.7,
                },
                Spectrum {
                    points: vec![0., 0.2, 0.4, 0.2, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.recalculated_with_step(0.2)
            );
        }
        #[test]
        fn _2_into_6() {
            assert_eq!(
                Spectrum {
                    // dx:       0. 0.2 0.4 0.6 0.8 1.0
                    // points: vec![0., 2., 4., 6., 8., 10.],
                    points: vec![0., 1.9999999999999996, 4.000000000000002, 6.000000000000001, 8., 10.],
                    step: 0.2,
                    x_start: 0.7,
                },
                Spectrum {
                    // dx:       0.  1.
                    points: vec![0., 10.],
                    step: 1.,
                    x_start: 0.7,
                }.recalculated_with_step(0.2)
            );
        }
        #[test]
        fn _9_into_4() {
            assert_eq!(
                Spectrum {
                    points: vec![0., 0.3, 0.8999999999999997, 0.2, 0.],
                    step: 0.8,
                    x_start: 0.7,
                },
                Spectrum {
                    points: vec![0., 0.1, 0.3, 0.5, 0.9, 0.6, 0.2, 0.1, 0.],
                    step: 0.4,
                    x_start: 0.7,
                }.recalculated_with_step(0.8)
            );
        }
    }
}

