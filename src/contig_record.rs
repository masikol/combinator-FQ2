
use std::sync::OnceLock;
use std::collections::HashSet;

use regex::Regex;

use crate::revcompl::revcompl;
use crate::iupac::nucl_bases::{
    GUANINE,
    CYTOSINE,
    STRONG,
};
use crate::fasta_reader::SeqRecord;


#[derive(Debug)]
pub struct ContigRecord {
    pub name: String,
    pub length: usize,
    pub gc_content: f64,
    pub coverage: Option<f64>,
    pub start: String,
    pub rcstart: String,
    pub end: String,
    pub rcend: String,
}

impl ContigRecord {

    pub fn from(seq_record: SeqRecord,
                maxk: usize) -> Result<ContigRecord, String> {

        let seq_len = seq_record.seq.len() as usize;

        let gc_count = calculate_gc_count(&seq_record.seq);
        let gc_content = gc_count / (seq_len as f64) * 100.0;

        let coverage = parse_coverage(&seq_record.name);

        let start_terminus = &seq_record.seq[..maxk];
        let end_terminus = &seq_record.seq[(seq_len-maxk)..];

        let start = start_terminus.to_string();
        let end = end_terminus.to_string();

        let rcstart = revcompl(start_terminus);
        if let Err(err_char) = rcstart {
            return Err(format!(
                "Reverse-complement sequence creating for {} failed on character `{}`",
                seq_record.name,
                err_char
            ));
        }
        let rcend = revcompl(end_terminus);
        if let Err(err_char) = rcend {
            return Err(format!(
                "Reverse-complement sequence creating for {} failed on character `{}`",
                seq_record.name,
                err_char
            ));
        }

        Ok(ContigRecord {
            name: seq_record.name,
            length: seq_len,
            gc_content: gc_content,
            coverage: coverage,
            start: start,
            rcstart: rcstart.unwrap(),
            end: end,
            rcend: rcend.unwrap(),
        })
    }
}


fn calculate_gc_count(seq: &String) -> f64 {
    let gc_base_set: HashSet<char> = HashSet::from_iter(
        vec![GUANINE, CYTOSINE, STRONG].into_iter()
    );

    let gc_count = seq.chars().filter(
        |c| gc_base_set.contains(c)
    ).count() as f64;

    gc_count
}


static SPADES_NAME_REGEX: OnceLock<Regex> = OnceLock::new();

fn get_spades_name_regex() -> &'static Regex {
    SPADES_NAME_REGEX.get_or_init(
        || Regex::new(
            r"^NODE_\d+_length_\d+_cov_(\d+\.(\d+)?)"
        ).unwrap()
    )
}

fn parse_coverage(seq_name: &String) -> Option<f64> {
    let re = get_spades_name_regex();
    let capture = re.captures(seq_name)?.get(1)?;
    let coverage: Result<f64, _> = capture.as_str().parse();
    if coverage.is_err() {
        panic!("Error: failed to parse coverage from header `{}`", seq_name);
    }
    let coverage = coverage.unwrap();
    Some(coverage)
}
