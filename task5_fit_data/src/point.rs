//! Points struct.

use crate::float_type::float;


#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    pub x: float,
    pub y: float,
}

impl Point {
    pub const fn new(x: float, y: float) -> Self {
        Self { x, y }
    }
}

