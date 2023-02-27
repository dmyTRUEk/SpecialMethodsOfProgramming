//! Points type.

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use crate::{float_type::float, point::Point};


pub type Points = Vec<Point>;

pub trait ImplPoints {
    fn from_array<const N: usize>(array: [(float, float); N]) -> Self;
    fn load_from_file(filename: &str) -> Self;
}
impl ImplPoints for Points {
    fn from_array<const N: usize>(array: [(float, float); N]) -> Self {
        array.into_iter()
            .map(|(x, y)| Point::new(x, y))
            .collect()
    }

    fn load_from_file(filename: &str) -> Self {
        let file = BufReader::new(File::open(filename).expect(&format!("can't open file `{}`", filename)));
        let mut points: Points = vec![];
        for line in file.lines() {
            let line = line.unwrap();
            let parts: Vec<_> = line.split(['\t', ' ']).collect();
            assert_eq!(2, parts.len());
            let x = parts[0];
            let y = parts[1];
            let x = x.replace(',', ".");
            let y = y.replace(',', ".");
            let x = x.parse().expect(&format!("unable to parse into float: `{}`", x));
            let y = y.parse().expect(&format!("unable to parse into float: `{}`", y));
            let point = Point::new(x, y);
            points.push(point);
        }
        points
    }
}

