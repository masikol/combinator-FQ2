
use crate::overlaps::{Overlap, OverlapCollection};
use crate::contig_record::{ContigRecord, MULTIPLTY_EPSILON};


pub fn assign_multiplty(contig_collection: &mut Vec<ContigRecord>,
                        overlap_collection: &OverlapCollection) {
    // The function assigns multiplicity (copies of this contig in the genome) to contigs.

    // Coverage of 1-st contig can be zero.
    // In this case we cannot calculate multiplicity of contigs based on coverage.
    // Set `first_cov_is_valid` to false if the first contig has no coverage.
    let first_cov_is_valid: bool = check_if_first_cov_is_valid(&contig_collection);
    let first_contig_cov: Option<f64> = contig_collection[0].coverage;

    for i in 0..contig_collection.len() {

        // Validate coverage of the current contig
        let calc_multiplty_by_cov = first_cov_is_valid
                                    && contig_collection[i].coverage.is_some();

        if calc_multiplty_by_cov {
            // Calculate multiplicity based on coverage
            contig_collection[i].set_multiplty_by_cov(first_contig_cov);
        } else {
            // Calculate multiplicity based on overlaps
            let overlaps: &[Overlap] = overlap_collection.get(&i);
            let coverage: f64 = calc_multiplty_by_overlaps(overlaps);
            contig_collection[i].multiplty = Some(coverage);
        }
    }
}


fn check_if_first_cov_is_valid(contig_collection: &Vec<ContigRecord>) -> bool {
    let first_cov: Option<f64> = contig_collection[0].coverage;
    if first_cov.is_none() {
        return false;
    }
    let first_cov = first_cov.unwrap();

    // Set `first_cov_is_valid` to false if the first contig has zero coverage.
    // And report it.
    if first_cov < MULTIPLTY_EPSILON {
        // Coverage of 1-st contig is zero
        eprintln!(
            "\n`{}` has coverage is less than {} (the actual value is {:.7}).",
            contig_collection[0].name,
            MULTIPLTY_EPSILON,
            first_cov
        );
        eprintln!("Multiplicity of contigs will be calculated based on overlaps instead of coverage.\n");
        return false;
    }

    true
}

fn calc_multiplty_by_overlaps(overlaps: &[Overlap]) -> f64 {
    // The function for calculating multiplicity of a given contig
    //   based on number of overlaps of this contig.

    // Count overlaps associated with start
    let num_start_matches: usize = overlaps.iter()
        .filter(|ovl: &&Overlap| ovl.is_start_match())
        .count();

    // Count overlaps associated with end
    let num_end_matches: usize = overlaps.iter()
        .filter(|ovl: &&Overlap| ovl.is_end_match())
        .count();

    // Obtain multiplicity based on number of overlaps
    let multiplicity: f64 = num_start_matches.min(
        num_end_matches
    ).max(1) as f64;

    multiplicity
}


// >>> Tests >>>

#[cfg(test)]
mod tests_check_if_first_cov_is_valid {
    use crate::contig_record::ContigRecord;

    use super::check_if_first_cov_is_valid;

    fn make_contig(cov: Option<f64>) -> ContigRecord {
        ContigRecord {
            name: "test".to_string(),
            length: 100,
            gc_content: 50.0,
            coverage: cov,
            multiplty: None,
            start: String::new(),
            rcstart: String::new(),
            end: String::new(),
            rcend: String::new(),
        }
    }

    #[test]
    fn cov_is_none_returns_false() {
        let contigs = vec![make_contig(None)];
        assert!(!check_if_first_cov_is_valid(&contigs));
    }

    #[test]
    fn cov_is_zero_returns_false() {
        let contigs = vec![make_contig(Some(0.0))];
        assert!(!check_if_first_cov_is_valid(&contigs));
    }

    #[test]
    fn cov_is_very_small_returns_false() {
        let contigs = vec![make_contig(Some(1e-7))];
        assert!(!check_if_first_cov_is_valid(&contigs));
    }

    #[test]
    fn cov_equals_epsilon_returns_true() {
        let contigs = vec![make_contig(Some(1e-6))];
        assert!(check_if_first_cov_is_valid(&contigs));
    }

    #[test]
    fn cov_normal_returns_true() {
        let contigs = vec![make_contig(Some(50.0))];
        assert!(check_if_first_cov_is_valid(&contigs));
    }

    #[test]
    #[should_panic]
    fn panics_on_empty_vec() {
        let contigs: Vec<ContigRecord> = vec![];
        check_if_first_cov_is_valid(&contigs);
    }
}


#[cfg(test)]
mod tests_calc_multiplty_by_overlaps {
    use std::path::PathBuf;

    use crate::args::Args;
    use crate::fasta_reader::FastaReader;
    use crate::contig_record::ContigRecord;
    use crate::overlaps::detect_adjacent_contigs;

    use super::calc_multiplty_by_overlaps;

    #[test]
    fn test_multiplty_1_through_2() {
        let file = "test_data/assign_multiplicity/test_multiplty_1-2.fasta";
        let maxk = 18;
        let mink = 16;

        let path = PathBuf::from(file);
        let args = Args {
            mink,
            maxk,
            input_fpath: path.clone(),
            outdir_path: PathBuf::from(""),
            force: false,
        };

        let reader = FastaReader::open(&path).unwrap();
        let contigs: Vec<ContigRecord> = reader
            .filter_map(|r| r.ok())
            .map(|s| ContigRecord::from(s, maxk).unwrap())
            .collect();

        let overlap_collection = detect_adjacent_contigs(&contigs, &args);

        let expected = 1.0;
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&0));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&1));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&3));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&4));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&5));
        assert!((multiplty - expected).abs() < f64::EPSILON);

        let expected = 2.0;
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&2));
        assert!((multiplty - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn test_multiplty_1_through_3() {
        let file = "test_data/assign_multiplicity/test_multiplty_1-3.fasta";
        let maxk = 18;
        let mink = 16;

        let path = PathBuf::from(file);
        let args = Args {
            mink,
            maxk,
            input_fpath: path.clone(),
            outdir_path: PathBuf::from(""),
            force: false,
        };

        let reader = FastaReader::open(&path).unwrap();
        let contigs: Vec<ContigRecord> = reader
            .filter_map(|r| r.ok())
            .map(|s| ContigRecord::from(s, maxk).unwrap())
            .collect();

        let overlap_collection = detect_adjacent_contigs(&contigs, &args);

        let expected = 1.0;
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&0));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&1));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&3));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&4));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&5));
        assert!((multiplty - expected).abs() < f64::EPSILON);
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&6));
        assert!((multiplty - expected).abs() < f64::EPSILON);

        let expected = 2.0;
        let multiplty = calc_multiplty_by_overlaps(overlap_collection.get(&2));
        assert!((multiplty - expected).abs() < f64::EPSILON);
    }
}


#[cfg(test)]
mod tests_assign_multiplty {
    use std::path::PathBuf;

    use crate::args::Args;
    use crate::fasta_reader::FastaReader;
    use crate::contig_record::ContigRecord;
    use crate::overlaps::detect_adjacent_contigs;

    use super::assign_multiplty;

    #[test]
    fn test_multiplty_1st_valid() {
        let file = "test_data/assign_multiplicity/test_multiplty_first_valid.fasta";
        let maxk = 18;
        let mink = 16;

        let path = PathBuf::from(file);
        let args = Args {
            mink,
            maxk,
            input_fpath: path.clone(),
            outdir_path: PathBuf::from(""),
            force: false,
        };

        let reader = FastaReader::open(&path).unwrap();
        let mut contigs: Vec<ContigRecord> = reader
            .filter_map(|r| r.ok())
            .map(|s| ContigRecord::from(s, maxk).unwrap())
            .collect();

        let overlaps = detect_adjacent_contigs(&contigs, &args);
        assign_multiplty(&mut contigs, &overlaps);

        let expected = 1.0;
        assert!((contigs[0].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[2].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[3].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[5].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[6].multiplty.unwrap() - expected).abs() < f64::EPSILON);

        let expected = 2.0;
        assert!((contigs[1].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[4].multiplty.unwrap() - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn test_multiplty_1st_invalid() {
        let file = "test_data/assign_multiplicity/test_multiplty_first_invalid.fasta";
        let maxk = 18;
        let mink = 16;

        let path = PathBuf::from(file);
        let args = Args {
            mink,
            maxk,
            input_fpath: path.clone(),
            outdir_path: PathBuf::from(""),
            force: false,
        };

        let reader = FastaReader::open(&path).unwrap();
        let mut contigs: Vec<ContigRecord> = reader
            .filter_map(|r| r.ok())
            .map(|s| ContigRecord::from(s, maxk).unwrap())
            .collect();

        let overlaps = detect_adjacent_contigs(&contigs, &args);
        assign_multiplty(&mut contigs, &overlaps);

        let expected = 1.0;
        assert!((contigs[0].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[1].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[2].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[3].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[5].multiplty.unwrap() - expected).abs() < f64::EPSILON);
        assert!((contigs[6].multiplty.unwrap() - expected).abs() < f64::EPSILON);

        let expected = 2.0;
        assert!((contigs[4].multiplty.unwrap() - expected).abs() < f64::EPSILON);
    }
}
