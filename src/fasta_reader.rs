
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

        if next_line == ">" {
            return Err(
                "Error: invalid fasta file: empty header line".to_string()
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

        if curr_seq.is_empty() {
            return Some(Err(
                "Error: empty fasta sequence encoutered".to_string()
            ));
        }

        let seq_name = curr_header_line[1..].to_string();
        if seq_name.is_empty() {
            return Some(Err(
                "Error: invalid fasta file: empty header line".to_string()
            ));
        }

        Some(Ok(SeqRecord {
            name: seq_name,
            seq: curr_seq,
        }))
    }
}


#[cfg(test)]
mod tests_fasta_reader {
    use super::*;

    fn test_path(name: &str) -> PathBuf {
        PathBuf::from("test_data").join("fasta_reader").join(name)
    }

    fn collect_records(path: &PathBuf) -> Result<Vec<SeqRecord>, String> {
        let reader = FastaReader::open(path)?;
        reader.collect()
    }

    #[test]
    fn empty_file() {
        assert!(FastaReader::open(&test_path("empty.fasta")).is_err());
    }

    #[test]
    fn newline_only() {
        assert!(FastaReader::open(&test_path("newline_only.fasta")).is_err());
    }

    #[test]
    fn valid_two_seqs() {
        let records = collect_records(&test_path("two_seqs.fasta"));
        assert!(records.is_ok());
    }

    #[test]
    fn valid_single_line_seqs() {
        let records = collect_records(&test_path("single_line_seqs.fasta"));
        assert!(records.is_ok());
    }

    #[test]
    fn valid_diff_line_lengths() {
        let records = collect_records(&test_path("diff_line_lengths.fasta"));
        assert!(records.is_ok());
    }

    #[test]
    fn valid_single_seq() {
        let records = collect_records(&test_path("single_seq.fasta"));
        assert!(records.is_ok());
    }

    #[test]
    fn invalid_no_leading_gt() {
        assert!(FastaReader::open(&test_path("no_leading_gt.fasta")).is_err());
    }

    #[test]
    fn invalid_second_header_no_gt() {
        let result: Result<Vec<SeqRecord>, String> =
            FastaReader::open(&test_path("second_header_no_gt.fasta"))
                .and_then(|r| r.collect());
        assert!(result.is_err());
    }

    #[test]
    fn invalid_blank_line() {
        let mut reader = FastaReader::open(&test_path("blank_line.fasta")).unwrap();
        assert!(reader.next().unwrap().is_err());
    }

    #[test]
    fn invalid_non_iupac() {
        let mut reader = FastaReader::open(&test_path("non_iupac.fasta")).unwrap();
        assert!(reader.next().unwrap().is_err());
    }

    #[test]
    fn invalid_empty_header() {
        assert!(FastaReader::open(&test_path("empty_header.fasta")).is_err());
    }

    #[test]
    fn invalid_empty_seq() {
        let records = collect_records(&test_path("empty_seq.fasta"));
        assert!(records.is_err());
    }

    #[test]
    fn valid_single_seq_multi_line() {
        let records = collect_records(&test_path("single_seq_multi_line.fasta"));
        assert!(records.is_ok());
    }

    #[test]
    fn invalid_second_header_empty() {
        let mut reader = FastaReader::open(&test_path("second_header_empty.fasta")).unwrap();
        assert!(reader.next().unwrap().is_ok());
        assert!(reader.next().unwrap().is_err());
    }
}
