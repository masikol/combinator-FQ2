
use std::path::PathBuf;

use clap::Parser;


#[derive(Parser)]
#[command(
    name = "combinator_fq2",
    version = "0.1.0",
    about = "A program to find adjacent contigs by matching their ends of length k.",
    long_about = None,
    author = "Maksim Sikolenko",
    help_template = "{name} {version}\nBy {author}\n\n{about}\n\nUSAGE:\n    {usage}\n\n{all-args}"
)]
struct RawArgs {

    /// Input fasta file
    #[arg(value_parser)]
    input_fpath: String,

    /// Minimum k to test.
    /// Integer > 0; Default is 21 bp
    #[arg(
        short = 'i',
        long,
        default_value_t = 21,
        value_parser 
    )]
    mink: usize,

    /// Maximum k to test.
    /// Integer > 0; Default is 127 bp
    #[arg(
        short = 'a',
        long,
        default_value_t = 127,
        value_parser 
    )]
    maxk:usize,

    /// Single k to test.
    /// If speified, `-i` and `-a` options are ignored.
    /// Integer > 0. Disabled by default
    #[arg(
        short = 'k',
        long = "k-mer",
        value_parser 
    )]
    k: Option<usize>,

    /// Output directory.
    /// Default value: `combinator-result`
    #[arg(
        short = 'o',
        long = "outdir",
        default_value_t = String::from("combinator-result"),
        value_parser 
    )]
    outdir_path: String,
}

impl RawArgs {
    fn validate(&self) -> Result<(), String> {

        if self.mink == 0 {
            return Err("Error: mink must not be zero".to_string());
        }
        if self.maxk == 0 {
            return Err("Error: mink must not be zero".to_string());
        }
        if let Some(k) = self.k {
            if k == 0 {
                return Err("Error: k must not be zero".to_string());
            }
        }

        if self.mink > self.maxk {
            return Err(
                format!(
                    "Error: mink ({}) is greater than maxk ({}).",
                    self.mink,
                    self.maxk
                )
            );
        }

        if !PathBuf::from(&self.input_fpath).is_file() {
            return Err(
                format!(
                    "Error: input file does not exist: `{}`",
                    self.input_fpath
                )
            );
        }

        Ok(())
    }
}


#[derive(Debug)]
pub struct Args {
    // TODO: add -f/force option
    pub mink: usize,
    pub maxk: usize,
    pub input_fpath: PathBuf,
    pub outdir_path: PathBuf,
}

impl Args {

    pub fn parse() -> Result<Args, String> {
        let raw_args = RawArgs::parse();
        if let Err(err_msg) = raw_args.validate() {
            eprintln!("{}", err_msg);
            return Err(err_msg);
        }

        let mut mink = raw_args.mink;
        let mut maxk = raw_args.maxk;
        if let Some(k) = raw_args.k {
            mink = k;
            maxk = k;
        }

        Ok(Args {
            mink: mink,
            maxk: maxk,
            input_fpath: PathBuf::from(&raw_args.input_fpath),
            outdir_path: PathBuf::from(&raw_args.outdir_path),
        })
    }
}
