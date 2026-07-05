
// TODO: use
// mod revcompl;


use std::fs::File;
use std::path::PathBuf;
use std::io::{BufReader, BufRead, Lines};


// TODO: remove test func
pub fn read_test_fasta(file_path: &PathBuf) {
    let reader = FastaReader::open(file_path).unwrap();
    for seq_record in reader {
        println!(">{}\n{}", seq_record.name, seq_record.seq);
    }
}


struct SeqRecord {
    name: String,
    seq: String,
}

// impl SeqRecord {
//     // TODO: GC content
//     // TODO: coverage
//     // TODO: length
//     // start
//     // rcstart
//     // end
//     // rcend
// }


struct FastaReader {
    file_lines: Lines<BufReader<File>>,
    next_header_line: Option<String>,
    end_of_file_reached: bool,
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
            )
        }

        Ok(FastaReader {
            file_lines: BufReader::new(
                open_result.unwrap()
            ).lines(),
            next_header_line: None,
            end_of_file_reached: false,
        })
    }
}

impl Iterator for FastaReader {
    type Item = SeqRecord;

    fn next(&mut self) -> Option<Self::Item> {

        if self.end_of_file_reached {
            return None;
        }

        if self.next_header_line.is_none() {
            // TODO: handle None and Err
            // self.file_lines.next() returns Option<Result<std::string::String, std::io::Error>>
            let line = self.file_lines.next().unwrap().unwrap();
            // TODO: unwrap unwrap unwrap...
            if line.chars().nth(0).unwrap() != '>' {
                // TODO: handle gracely
                panic!("Error: invalid input fasta data: the first line does not start with '>'");
            }
            self.next_header_line = Some(line);
        }

        let mut curr_seq = String::from("");
        let mut curr_header_line = String::from("");

        let mut next_record_reached = false;
        while !next_record_reached {

            let next_line_option = self.file_lines.next();
            if next_line_option.is_none() {
                curr_header_line = self.next_header_line.clone().unwrap();
                self.next_header_line = None;
                self.end_of_file_reached = true;
                break;
            }

            // TODO: handle Err
            // TODO: to_string() overhead
            let mut next_line = next_line_option
                .unwrap()
                .expect("Error: cannot read fasta line")
                .trim()
                .to_string();

            // TODO: if line is ""
            // TODO: if line contains non-IUPAC chars

            if next_line.chars().nth(0).unwrap() != '>' {
                next_line.make_ascii_uppercase();
                curr_seq.push_str(&next_line);
            } else {
                // TODO: lame unwrap and clone
                // TODO: does the execution even reach it??
                curr_header_line = self.next_header_line.clone().unwrap();
                self.next_header_line = Some(String::from(next_line));
                next_record_reached = true;
            }
        }

        // TODO: lame error handling: clone(), clone()...
        let name = String::from(
            curr_header_line.clone().strip_prefix(">").unwrap_or_else(
                || panic!(
                    "Error: invalid fasta header: `{}`",
                    curr_header_line.clone()
                )
            )
        );
        return Some(SeqRecord {
            name: name,
            seq: curr_seq,
        })
    }
}
