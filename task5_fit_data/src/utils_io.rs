//! IO utils.

use std::io::{Write, stdin, stdout};


pub fn press_enter_to_continue() {
    print("PRESS ENTER TO CONTINUE.");
    wait_for_enter();
}

pub fn wait_for_enter() {
    let mut line: String = String::new();
    stdin().read_line(&mut line).unwrap();
}

pub fn print(msg: impl ToString) {
    print!("{}", msg.to_string());
    stdout().flush().unwrap();
}

