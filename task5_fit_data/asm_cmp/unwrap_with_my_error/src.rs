//! Decompilation in collab with @hedghg.

pub fn unwrap_with_my_error(n: Option<usize>) -> Result<usize, String> {
    // if n.is_none() { return Err("my custom error".to_string()) }
    // Ok(n.unwrap() / 10)

    match n {
        None => return Err("my custom error".to_string()),
        Some(index_of_max) => Ok(index_of_max / 10),
    }
}

