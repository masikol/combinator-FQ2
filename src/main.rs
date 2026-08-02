
use std::process::ExitCode;

use combinator_fq2::run_combinator_fq2;


fn main() -> ExitCode {
    match run_combinator_fq2() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err_msg) => {
            eprintln!("{}", err_msg);
            ExitCode::FAILURE
        }
    }
}

