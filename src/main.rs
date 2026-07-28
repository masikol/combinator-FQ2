
mod args;
mod iupac;
mod output;
mod revcompl;
mod overlaps;
mod fasta_reader;
mod find_overlap;
mod contig_record;
mod cov_summarizer;
mod assign_multiplicity;


use std::fs;
use std::process::ExitCode;

use args::Args;
use output as out;
use fasta_reader::FastaReader;
use assign_multiplicity as amu;
use contig_record::ContigRecord;
use overlaps::{OverlapCollection, detect_adjacent_contigs};


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

    let contig_records = read_contig_records(&args);
    if let Err(err_str) = contig_records {
        eprintln!("{}", err_str);
        return ExitCode::FAILURE;
    }
    let mut contig_records: Vec<ContigRecord> = contig_records.unwrap();

    println!("{:?}", contig_records);

    for r in &contig_records {
        println!("{:?}", r);
    }

    let overlaps: OverlapCollection = detect_adjacent_contigs(&contig_records, &args);
    for (i, ovl_vec) in &overlaps.collection {
        println!("Overlaps for contig #{:?}", i);
        for ovl in ovl_vec {
            println!("{:?}", ovl);
        }
    }

    amu::assign_multiplty(&mut contig_records, &overlaps);

    // Write full matching log
    let out_result = out::write_full_log(
        &contig_records,
        &overlaps,
        &args
    );
    if let Err(err_str) = out_result {
        eprintln!("{}", err_str);
        return ExitCode::FAILURE;
    }

    // Write adjacency table
    let out_result = out::write_adjacency_table(
        &contig_records,
        &overlaps,
        &args
    );
    if let Err(err_str) = out_result {
        eprintln!("{}", err_str);
        return ExitCode::FAILURE;
    }

    // Write summary
    let out_result = out::write_summary(
        &contig_records,
        &overlaps,
        &args
    );
    if let Err(err_str) = out_result {
        eprintln!("{}", err_str);
        return ExitCode::FAILURE;
    }

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

    // TODO: use collect() instead
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
