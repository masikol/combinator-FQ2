
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

pub const MULTIPLTY_EPSILON: f64 = 1e-6;


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
    pub multiplty: Option<f64>,
}

impl ContigRecord {

    pub fn from(seq_record: SeqRecord,
                maxk: usize) -> Result<ContigRecord, String> {

        let seq_len = seq_record.seq.len() as usize;

        let gc_count = calculate_gc_count(&seq_record.seq);
        let gc_content = calculate_gc_content(gc_count, seq_len);

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
            multiplty: None,
        })
    }

    pub fn set_multiplty_by_cov(&mut self, first_contig_cov: f64) {
        if first_contig_cov < MULTIPLTY_EPSILON {
            panic!("Error: first_contig_cov ({}) is zero or too small", first_contig_cov);
        }
        let cov_ratio: f64 = self.coverage.unwrap() / first_contig_cov;
        self.multiplty = Some(cov_ratio);
    }
}


fn calculate_gc_count(seq: &String) -> usize {
    let gc_base_set: HashSet<char> = HashSet::from_iter(
        vec![GUANINE, CYTOSINE, STRONG].into_iter()
    );

    let gc_count = seq.chars().filter(
        |c| gc_base_set.contains(c)
    ).count();

    gc_count
}

fn calculate_gc_content(gc_count: usize, seq_len: usize) -> f64 {
    (gc_count as f64) / (seq_len as f64) * 100.0
}


static SPADES_NAME_REGEX: OnceLock<Regex> = OnceLock::new();

fn get_spades_name_regex() -> &'static Regex {
    SPADES_NAME_REGEX.get_or_init(
        || Regex::new(
            r"^NODE_\d+_length_\d+_cov_(\d+(\.\d+)?)"
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


#[cfg(test)]
mod tests_parse_coverage {
    use super::parse_coverage;

    #[test]
    fn standard_spades_header() {
        assert_eq!(
            parse_coverage(&"NODE_1_length_100_cov_50.0".into()),
            Some(50.0)
        );
    }

    #[test]
    fn decimal_coverage() {
        assert_eq!(
            parse_coverage(&"NODE_5_length_250_cov_30.75".into()),
            Some(30.75)
        );
    }

    #[test]
    fn integer_coverage_with_trailing_dot() {
        assert_eq!(
            parse_coverage(&"NODE_3_length_50_cov_100.".into()),
            Some(100.0)
        );
    }

    #[test]
    fn zero_coverage() {
        assert_eq!(
            parse_coverage(&"NODE_1_length_100_cov_0.0".into()),
            Some(0.0)
        );
    }

    #[test]
    fn non_spades_header() {
        assert_eq!(parse_coverage(&">seq1".into()), None);
    }

    #[test]
    fn empty_string() {
        assert_eq!(parse_coverage(&String::new()), None);
    }

    #[test]
    fn missing_cov_part() {
        assert_eq!(parse_coverage(&"NODE_1_length_100".into()), None);
    }

    #[test]
    fn coverage_without_decimal() {
        assert_eq!(
            parse_coverage(&"NODE_1_length_100_cov_50".into()),
            Some(50.0)
        );
    }
}


#[cfg(test)]
mod tests_contig_record_from {
    use super::ContigRecord;
    use crate::fasta_reader::SeqRecord;

    fn make_seq(name: &str, seq: &str) -> SeqRecord {
        SeqRecord {
            name: name.to_string(),
            seq: seq.to_string(),
        }
    }

    #[test]
    fn standard_spades_header() {
        let record = ContigRecord::from(
            make_seq("NODE_1_length_100_cov_50.0", "ACGTACGTAC"),
            3,
        ).unwrap();
        assert_eq!(record.name, "NODE_1_length_100_cov_50.0");
        assert_eq!(record.length, 10);
        assert!((record.gc_content - 50.0).abs() < f64::EPSILON);
        assert_eq!(record.coverage, Some(50.0));
        assert_eq!(record.start, "ACG");
        assert_eq!(record.end, "TAC");
        assert_eq!(record.rcstart, "CGT");
        assert_eq!(record.rcend, "GTA");
    }

    #[test]
    fn non_spades_header() {
        let record = ContigRecord::from(
            make_seq("contig1", "GGGGCCCC"),
            2,
        ).unwrap();
        assert_eq!(record.name, "contig1");
        assert_eq!(record.coverage, None);
    }

    #[test]
    fn maxk_equals_seq_len() {
        let record = ContigRecord::from(
            make_seq("NODE_1_length_5_cov_10.0", "GATCA"),
            5,
        ).unwrap();
        assert_eq!(record.start, "GATCA");
        assert_eq!(record.end, "GATCA");
        assert_eq!(record.rcstart, "TGATC");
        assert_eq!(record.rcend, "TGATC");
    }

    #[test]
    fn non_iupac_character() {
        let result = ContigRecord::from(
            make_seq("NODE_1_length_5_cov_10.0", "GATCZ"),
            3,
        );
        assert!(result.is_err());
    }
}


#[cfg(test)]
mod tests_calculate_gc_count {
    use super::calculate_gc_count;

    #[test]
    fn empty_string() {
        assert_eq!(calculate_gc_count(&String::new()), 0);
    }

    #[test]
    fn no_gc_bases() {
        assert_eq!(calculate_gc_count(&"ATAT".into()), 0);
    }

    #[test]
    fn only_gc_bases() {
        assert_eq!(calculate_gc_count(&"GGCC".into()), 4);
    }

    #[test]
    fn mixed_bases() {
        assert_eq!(calculate_gc_count(&"ACGTS".into()), 3);
    }

    #[test]
    fn iupac_s_is_gc() {
        assert_eq!(calculate_gc_count(&"S".into()), 1);
    }

    #[test]
    fn lowercase_not_counted() {
        assert_eq!(calculate_gc_count(&"gc".into()), 0);
    }

    #[test]
    fn no_gc_in_iupac_non_gc() {
        assert_eq!(calculate_gc_count(&"RYWKM".into()), 0);
    }
}


#[cfg(test)]
mod tests_calculate_gc_content {
    use super::calculate_gc_content;

    #[test]
    fn zero_gc_count() {
        let content = calculate_gc_content(0, 100);
        assert!((content - 0.0).abs() < f64::EPSILON);
    }

    #[test]
    fn half_gc() {
        let content = calculate_gc_content(50, 100);
        assert!((content - 50.0).abs() < f64::EPSILON);
    }

    #[test]
    fn all_gc() {
        let content = calculate_gc_content(100, 100);
        assert!((content - 100.0).abs() < f64::EPSILON);
    }

    #[test]
    fn quarter_gc() {
        let content = calculate_gc_content(25, 100);
        assert!((content - 25.0).abs() < f64::EPSILON);
    }

    #[test]
    fn odd_length() {
        let content = calculate_gc_content(3, 9);
        let one_third: f64 = 1.0 / 3.0 * 100.0;
        assert!((content - one_third).abs() < f64::EPSILON);
    }

    #[test]
    fn single_base_gc() {
        let content = calculate_gc_content(1, 1);
        assert!((content - 100.0).abs() < f64::EPSILON);
    }

    #[test]
    fn single_base_no_gc() {
        let content = calculate_gc_content(0, 1);
        assert!((content - 0.0).abs() < f64::EPSILON);
    }
}


#[cfg(test)]
mod tests_set_multiplty_by_cov {
    use super::*;

    fn make_contig_record(coverage: Option<f64>) -> ContigRecord {
        ContigRecord {
            name: "test".to_string(),
            length: 100,
            gc_content: 50.0,
            coverage,
            multiplty: None,
            start: String::new(),
            rcstart: String::new(),
            end: String::new(),
            rcend: String::new(),
        }
    }

    #[test]
    #[should_panic]
    fn panics_when_coverage_is_none() {
        let mut record = make_contig_record(None);
        record.set_multiplty_by_cov(42.0);
    }

    #[test]
    #[should_panic]
    fn panics_when_first_cov_is_zero() {
        let mut record = make_contig_record(Some(50.0));
        record.set_multiplty_by_cov(0.0);
    }

    #[test]
    #[should_panic]
    fn panics_when_first_cov_is_very_small() {
        let mut record = make_contig_record(Some(50.0));
        let very_small_cov = MULTIPLTY_EPSILON / 10.0;
        record.set_multiplty_by_cov(very_small_cov);
    }

    #[test]
    fn sets_multiplicity_correctly() {
        let mut record = make_contig_record(Some(50.0));
        record.set_multiplty_by_cov(25.0);
        let result = record.multiplty.unwrap();
        let expected: f64 = 2.0;
        assert!((result - expected).abs() < f64::EPSILON);
        assert_eq!(record.multiplty, Some(expected));
    }

    #[test]
    fn fractional_multiplicity() {
        let mut record = make_contig_record(Some(10.0));
        record.set_multiplty_by_cov(100.0);
        let result = record.multiplty.unwrap();
        let expected: f64 = 0.1;
        assert!((result - expected).abs() < f64::EPSILON);
        assert_eq!(record.multiplty, Some(expected));
    }

    #[test]
    fn equal_coverage() {
        let mut record = make_contig_record(Some(42.0));
        record.set_multiplty_by_cov(42.0);
        let result = record.multiplty.unwrap();
        let expected: f64 = 1.0;
        assert!((result - expected).abs() < f64::EPSILON);
        assert_eq!(record.multiplty, Some(expected));
    }
}
