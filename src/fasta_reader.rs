
use std::fs::File;
use std::path::PathBuf;
use std::io::{BufReader, BufRead, Lines, Error as IOError};

use crate::seq_record::SeqRecord;
use crate::iupac::IUPACValidator;

type FastaLines = Lines<BufReader<File>>;


// TODO: remove test func
pub fn read_test_fasta(file_path: &PathBuf) -> Result<(), String> {
    let reader = FastaReader::open(file_path);
    if let Err(error) = reader {
        return Err(
            format!("Error. Cannot open fasta file `{}`: {}", file_path.display(), error)
        );
    }
    let reader = reader.unwrap();
    for data_chunk in reader {
        match data_chunk {
            Ok(seq_record) => {
                println!("Name: `{}`, Seq: `{}`", seq_record.name, seq_record.seq);
            },
            Err(err_str) => {return Err(err_str)}
        }
    }
    Ok(())
}


struct FastaReader {
    file_lines: FastaLines,
    next_header_line: String,
    end_of_file_reached: bool,
    seq_validator: IUPACValidator,
}

impl FastaReader {
    fn open(file_path: &PathBuf) -> Result<FastaReader, String> {
        let open_result = File::open(file_path);
        if let Err(error) = open_result {
            return Err(
                format!(
                    "Problem opening the file `{}`: {:?}",
                    file_path.display(),
                    error
                )
            );
        }

        let mut file_lines = BufReader::new(
            open_result.unwrap()
        ).lines();

        let next_header_line = FastaReader::read_first_line(&mut file_lines)?;

        let reader = FastaReader {
            file_lines,
            next_header_line: next_header_line,
            end_of_file_reached: false,
            seq_validator: IUPACValidator::new()
        };

        Ok(reader)
    }

    fn read_first_line(file_lines: &mut FastaLines) -> Result<String, String> {
        let next_line: Option<Result<String, IOError>> = file_lines.next();
        if next_line.is_none() {
            return Err(
                "Error: empty input file or invalid fasta data".to_string()
            );
        }

        let next_line: Result<String, IOError> = next_line.unwrap();
        if let Err(error) = next_line {
            return Err(
                format!("Error: cannot read fasta line: {}", error)
            );
        }

        let next_line: String = next_line.unwrap();
        let first_char: Option<char> = next_line.chars().nth(0);
        if first_char.is_none() {
            return Err(
                "Error: invalid input fasta data: the first line is empty".to_string()
            );
        }

        let first_char: char = first_char.unwrap();
        if first_char != '>' {
            return Err(
                "Error: invalid input fasta data: \
                the first line does not start with '>'".to_string()
            );
        }

        Ok(next_line)
    }
}

impl Iterator for FastaReader {
    type Item = Result<SeqRecord, String>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.end_of_file_reached {
            return None;
        }

        let mut next_record_reached = false;
        let mut curr_seq = String::from("");
        let mut curr_header_line = String::from("");

        while !next_record_reached {

            let next_line: Option<Result<String, IOError>> = self.file_lines.next();
            if next_line.is_none() {
                curr_header_line = self.next_header_line.clone();
                self.end_of_file_reached = true;
                break;
            }

            let next_line: Result<String, IOError> = next_line.unwrap();
            if let Err(error) = next_line {
                return Some(Err(
                    format!("Error: cannot read fasta line: {}", error)
                ))
            }

            let mut next_line: String = next_line.unwrap().trim().to_string();

            let first_char: Option<char> = next_line.chars().nth(0);
            if first_char.is_none() {
                return Some(Err(
                    "Error: unexpected empty line in the input fasta file".to_string()
                ));
            }
            let first_char = first_char.unwrap();

            if first_char != '>' {
                // Append sequence line to curr_seq
                if let Err(err_str) = self.seq_validator.validate(&next_line) {
                    return Some(Err(err_str));
                }
                next_line.make_ascii_uppercase();
                curr_seq.push_str(&next_line);
            } else {
                // We have reached the next record.
                // Save it’s header line to self.next_header_line
                //   and end the loop
                curr_header_line = self.next_header_line.clone();
                self.next_header_line = String::from(next_line);
                next_record_reached = true;
            }
        }

        let seq_name = curr_header_line[1..].to_string();

        Some(Ok(SeqRecord {
            name: seq_name,
            seq: curr_seq,
        }))
    }
}
