
mod args;
mod find_overlap;


use std::fs;
use std::process::ExitCode;

use args::Args;
use find_overlap::{find_overlap_s2s, find_overlap_e2s, find_overlap_e2e};


fn main() -> ExitCode {

    let args = match Args::parse() {
        Ok(args) => args,
        Err(err_msg) => {
            eprintln!("{}", err_msg);
            return ExitCode::FAILURE;
        }
    };
    println!("{:?}", args);

    if create_outdir(&args).is_err() {
        return ExitCode::FAILURE;
    }

    let overlap: usize = find_overlap_s2s(
        &String::from("AGTCaaaaaaaaa"),
        &String::from("AGTCttttttttttttttt"),
        3,
        5
    );
    println!("overlap = {overlap}");

    let overlap: usize = find_overlap_e2s(
        &String::from("aaaaaaaaaaAGTC"),
        &String::from("AGTCttttttttttttt"),
        3,
        5
    );
    println!("overlap = {overlap}");

    let overlap: usize = find_overlap_e2e(
        &String::from("aaaaaaaAGTC"),
        &String::from("tttttttttttttAGTC"),
        3,
        5
    );
    println!("overlap = {overlap}");

    ExitCode::SUCCESS
}

fn create_outdir(args: &Args) -> Result<(), ()> {
    if ! args.outdir_path.is_dir() {
        if let Err(error) = fs::create_dir(&args.outdir_path) {
            eprintln!(
                "Error: cannot create directory `{}`",
                args.outdir_path.display()
            );
            eprintln!("Reason: {}", error);
            return Err(());
        }
    }
    Ok(())
}
