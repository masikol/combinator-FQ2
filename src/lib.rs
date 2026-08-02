//! # combinator-FQ2
//!
//! `combinator_fq2` is a program thet detects adjacent contigs
//! by comparing their k-mer termini: starts and ends.
//!
//! It might be useful when finishing genome assemblies after, e.g. SPAdes and A5 assemblers.
//!
//! Please see the details on GitHub:
//! [https://github.com/masikol/combinator-FQ2](https://github.com/masikol/combinator-FQ2).

pub mod args;
pub mod output;
pub mod overlaps;
pub mod fasta_reader;
pub mod contig_record;
pub mod assign_multiplicity;

mod iupac;
mod spades;
mod revcompl;
mod find_overlap;
mod cov_summarizer;


pub use args::Args;
pub use fasta_reader::FastaReader;
pub use contig_record::ContigRecord;
pub use overlaps::detect_adjacent_contigs;
pub use assign_multiplicity::assign_multiplicity;
pub use output::{
    write_adjacency_table,
    write_full_log,
    write_summary,
    FULL_LOG_FILENAME,
    ADJ_TABLE_FILENAME,
    SUMMARY_FILENAME,
};


use std::fs;
use std::path::PathBuf;


/// Runs the combinator_fq2 pipeline end-to-end.
///
/// Parses the CLI arguments, creates the output directory, reads the
/// input contigs from the FASTA file, detects adjacent contigs,
/// assigns multiplicities, and writes the full matching log,
/// adjacency table, and summary into the output directory.
///
/// # Errors
///
/// Returns `Err` holding a `String` on any pipeline failure:
/// invalid or unreadable input, an output-directory problem,
/// or an output-write error.
pub fn run_combinator_fq2() -> Result<(), String> {
    let args = Args::parse()?;
    println!("{:?}", args);

    create_outdir(&args)?;

    let mut contig_records = read_contig_records(&args)?;
    let overlaps = detect_adjacent_contigs(&contig_records, &args);
    assign_multiplicity(&mut contig_records, &overlaps);

    write_adjacency_table(&contig_records, &overlaps, &args)?;
    write_full_log(&contig_records, &overlaps, &args)?;
    write_summary(&contig_records, &overlaps, &args)?;

    Ok(())
}

// Create outdir if it doesn’t exist
fn create_outdir(args: &Args) -> Result<(), String> {
    if ! args.outdir_path.is_dir() {
        if let Err(error) = fs::create_dir(&args.outdir_path) {
            return Err(format!(
                "Error: cannot create directory `{}`\nReason: {}",
                args.outdir_path.display(),
                error
            ));
        }
    }
    // Return Err if any output file exists and --force is false
    if ! args.force {
        err_if_output_exists(args)?;
    }
    Ok(())
}

// Return Err if any output file exists and --force is false
fn err_if_output_exists(args: &Args) -> Result<(), String> {
    let all_out_file_names = [
        FULL_LOG_FILENAME,
        ADJ_TABLE_FILENAME,
        SUMMARY_FILENAME,
    ];

    for filename in all_out_file_names {
        let fpath: PathBuf = args.outdir_path.join(filename);
        if fpath.is_file() {
            return Err(format!(
                "Output file {:?} already exists.\n\
                Cowardly refusing to overwrite.\n\
                Use -f / --force to overwrite.",
                fpath.display()
            ));
        }
    }

    Ok(())
}

/// Reads contig sequences from the input FASTA file, keeping only
/// records whose length is at least `args.maxk`.
///
/// # Errors
///
/// Returns `Err` if:
/// - the FASTA file cannot be opened;
/// - a sequence record fails IUPAC validation (e.g. when computing its
///   reverse-complement end during conversion to a `ContigRecord`);
/// - no input sequence reaches the minimum length `args.maxk`.
pub fn read_contig_records(args: &Args) -> Result<Vec<ContigRecord>, String> {

    let reader = FastaReader::open(&args.input_fpath);
    if let Err(error) = reader {
        return Err(
            format!(
                "Error. Cannot open fasta file `{}`: {}",
                args.input_fpath.display(),
                error
            )
        );
    }

    let mut contig_records: Vec<ContigRecord> = Vec::new();
    let reader = reader.unwrap();

    for seq_record in reader {
        match seq_record {
            Ok(seq_record) => {
                let seq_len = seq_record.seq.len();
                if seq_len >= args.maxk {
                    contig_records.push(
                        ContigRecord::from(seq_record, args.maxk)?
                    );
                }
            },
            Err(err_str) => {
                return Err(err_str);
            }
        }
    }

    if contig_records.is_empty() {
        return Err(format!(
            "Error: found no input sequences that have length >= {}",
            args.maxk
        ));
    }

    Ok(contig_records)
}
