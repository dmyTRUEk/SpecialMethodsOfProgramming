//! Param struct.

use rand::{Rng, thread_rng};

use crate::{extensions::ExtGenFromArray, float_type::float, FUNCTION_PARAM_VALUE_MAX, FUNCTION_PARAM_VALUE_MIN};


pub type ParamName = char;
pub type ParamValue = float; // use `ParamValue` instead of `float` where needed.


pub const PARAMETER_NAMES: [ParamName; 20] = [
    'a', 'b', 'c', 'd',
    // 'e', // removed to not confuse it with e = 2.71….
    // 'f',
    'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q',
    // 'r', // removed to not confuse it with polar coordinates radius.
    's', 't', 'u', 'v', 'w',
    // 'x', 'y', 'z' // removed to not confuse them with variables.
];


#[derive(Debug, Clone, PartialEq)]
pub struct Param {
    pub name: ParamName,
    pub value: ParamValue,
}

#[allow(dead_code)]
impl Param {
    pub const fn new(name: ParamName, value: ParamValue) -> Self {
        Self { name, value }
    }

    pub fn gen_random_name_with_value(value: ParamValue) -> Self {
        Self::new(Self::gen_random_name(), value)
    }

    pub fn gen_random_value_with_name(name: ParamName) -> Self {
        Self::new(name, Self::gen_random_value())
    }

    pub fn gen_random_name_and_value() -> Self {
        Self::new(Self::gen_random_name(), Self::gen_random_value())
    }

    pub fn gen_random_name() -> ParamName {
        thread_rng().gen_from_array(PARAMETER_NAMES)
    }

    pub fn gen_random_value() -> ParamValue {
        thread_rng().gen_range(FUNCTION_PARAM_VALUE_MIN ..= FUNCTION_PARAM_VALUE_MAX)
    }
}

