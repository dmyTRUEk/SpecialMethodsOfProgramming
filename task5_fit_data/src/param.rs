//! Param struct.

use rand::{Rng, thread_rng};

use crate::{extensions::ExtGenFromArray, float_type::float};


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

impl Param {
    pub const fn new(name: ParamName, value: ParamValue) -> Self {
        Self { name, value }
    }

    #[allow(dead_code)]
    pub fn gen_name(value: ParamValue) -> Self {
        let mut rng = thread_rng();
        Self {
            name: rng.gen_from_array(PARAMETER_NAMES),
            value,
        }
    }

    pub fn gen_value(name: ParamName) -> Self {
        let mut rng = thread_rng();
        Self {
            name,
            value: rng.gen_range(-3. ..= 5.),
        }
    }

    #[allow(dead_code)]
    pub fn gen_name_and_value() -> Self {
        let mut rng = thread_rng();
        Self {
            name: rng.gen_from_array(PARAMETER_NAMES),
            value: rng.gen_range(-3. ..= 5.),
        }
    }
}

