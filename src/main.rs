
mod args;
mod iupac;
mod spades;
mod output;
mod revcompl;
mod overlaps;
mod fasta_reader;
mod find_overlap;
mod contig_record;
mod cov_summarizer;
mod assign_multiplicity;


use std::fs;
use std::path::PathBuf;
use std::process::ExitCode;

use args::Args;
use output as out;
use fasta_reader::FastaReader;
use assign_multiplicity as amu;
use contig_record::ContigRecord;
use overlaps::detect_adjacent_contigs;


fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err_msg) => {
            eprintln!("{}", err_msg);
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let args = Args::parse()?;
    println!("{:?}", args);

    create_outdir(&args)?;

    let mut contig_records = read_contig_records(&args)?;
    let overlaps = detect_adjacent_contigs(&contig_records, &args);
    amu::assign_multiplty(&mut contig_records, &overlaps);

    out::write_full_log(&contig_records, &overlaps, &args)?;
    out::write_adjacency_table(&contig_records, &overlaps, &args)?;
    out::write_summary(&contig_records, &overlaps, &args)?;

    Ok(())
}

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
    if ! args.force {
        err_if_output_exists(args)?;
    }
    Ok(())
}

fn err_if_output_exists(args: &Args) -> Result<(), String> {
    let all_out_file_names = [
        out::FULL_LOG_FILENAME,
        out::ADJ_TABLE_FILENAME,
        out::SUMMARY_FILENAME,
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

fn read_contig_records(args: &Args) -> Result<Vec<ContigRecord>, String> {

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
