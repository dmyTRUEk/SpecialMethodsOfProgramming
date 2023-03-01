const X_MAX: f64 = 1.37552;

fn fr(x: f64) -> Result<f64, &'static str> {
    if x > X_MAX {
        Ok(x)
    } else {
        Err("this is my error msg")
    }
}

fn fo(x: f64) -> Option<f64> {
    if x > X_MAX {
        Some(x)
    } else {
        None
    }
}

pub fn o(x: f64) -> Option<f64> {
    fr(x).ok()
    // fo(x)
}

