
use crate::contig_record::ContigRecord;


/// Computes coverage statistics (min, max, mean, median) across the
/// contigs of an assembly.
///
/// Contigs without coverage are ignored. When no contig has coverage,
/// every statistic is rendered as `"NA"`.
pub struct CovSummarizer {
    cov_vector: Vec<f64>,
}

impl CovSummarizer {
    /// Collects the coverages of all `contig_collection` entries,
    /// skipping contigs without coverage.
    pub fn from(contig_collection: &Vec<ContigRecord>) -> CovSummarizer {
        CovSummarizer {
            cov_vector: contig_collection.iter()
                .map(|contig| contig.coverage)
                .filter_map(|coverage| coverage)
                .collect(),
        }
    }

    /// Minimum coverage, formatted to two decimals,
    /// or `"NA"` if there are no coverages.
    pub fn min_str(&self) -> String {
        if self.cov_vector.is_empty() {
            return String::from("NA");
        }
        format!(
            "{:.2}",
            self.cov_vector.iter()
                .min_by(|a, b| a.total_cmp(b))
                .unwrap()
        )
    }

    /// Maximum coverage, formatted to two decimals,
    /// or `"NA"` if there are no coverages.
    pub fn max_str(&self) -> String {
        if self.cov_vector.is_empty() {
            return String::from("NA");
        }
        format!(
            "{:.2}",
            self.cov_vector.iter()
                .max_by(|a, b| a.total_cmp(b))
                .unwrap()
        )
    }

    /// Mean coverage, formatted to two decimals,
    /// or `"NA"` if there are no coverages.
    ///
    /// Computed by dividing each value first, then summing, to avoid
    /// float overflow.
    pub fn mean_str(&self) -> String {
        if self.cov_vector.is_empty() {
            return String::from("NA");
        }
        // Avoiding float overflow: divide first, add then
        let contig_count = self.cov_vector.len() as f64;
        format!(
            "{:.2}",
            self.cov_vector.iter()
                .map(|cov| cov / contig_count)
                .sum::<f64>()
        )
    }

    /// Median coverage, formatted to two decimals,
    /// or `"NA"` if there are no coverages.
    pub fn median_str(&self) -> String {
        if self.cov_vector.is_empty() {
            return String::from("NA");
        }
        format!(
            "{:.2}",
            median(&self.cov_vector)
        )
    }

}


/// Median of a non-empty slice of values: the middle value for an odd
/// count, or the mean of the two middle values for an even count.
fn median(vec: &Vec<f64>) -> f64 {
    let mut vec_sorted = vec.clone();
    vec_sorted.sort_by(
        |a, b| a.total_cmp(b)
    );

    let vec_length = vec_sorted.len();

    if vec_length % 2 == 0 {
        let right_val_idx: usize = vec_length / 2;
        let left_val_idx:  usize = right_val_idx - 1;
        return (
            vec_sorted[right_val_idx] + vec_sorted[left_val_idx]
        ) / 2.0;
    } else {
        let median_val_idx: usize = vec_length / 2;
        return vec_sorted[median_val_idx];
    }
}


// >>> Tests >>>

#[cfg(test)]
mod tests_min_str {
    use super::*;

    fn make_contig(cov: Option<f64>) -> ContigRecord {
        ContigRecord {
            name: "t".into(),
            length: 10,
            gc_content: 50.0,
            coverage: cov,
            start: "AAAAAAAAAA".into(),
            rcstart: "TTTTTTTTTT".into(),
            end: "CCCCCCCCCC".into(),
            rcend: "GGGGGGGGGG".into(),
            multiplty: None,
        }
    }

    #[test]
    fn empty_returns_na() {
        let records = vec![make_contig(None)];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.min_str(), "NA");
    }

    #[test]
    fn single_element() {
        let records = vec![make_contig(Some(42.5))];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.min_str(), "42.50");
    }

    #[test]
    fn picks_minimum() {
        let records = vec![
            make_contig(Some(10.0)),
            make_contig(Some(5.0)),
            make_contig(Some(20.0)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.min_str(), "5.00");
    }
}


#[cfg(test)]
mod tests_max_str {
    use super::*;

    fn make_contig(cov: Option<f64>) -> ContigRecord {
        ContigRecord {
            name: "t".into(),
            length: 10,
            gc_content: 50.0,
            coverage: cov,
            start: "AAAAAAAAAA".into(),
            rcstart: "TTTTTTTTTT".into(),
            end: "CCCCCCCCCC".into(),
            rcend: "GGGGGGGGGG".into(),
            multiplty: None,
        }
    }

    #[test]
    fn empty_returns_na() {
        let records = vec![make_contig(None)];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.max_str(), "NA");
    }

    #[test]
    fn single_element() {
        let records = vec![make_contig(Some(42.5))];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.max_str(), "42.50");
    }

    #[test]
    fn picks_maximum() {
        let records = vec![
            make_contig(Some(10.0)),
            make_contig(Some(5.0)),
            make_contig(Some(20.0)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.max_str(), "20.00");
    }
}


#[cfg(test)]
mod tests_mean_str {
    use super::*;

    fn make_contig(cov: Option<f64>) -> ContigRecord {
        ContigRecord {
            name: "t".into(),
            length: 10,
            gc_content: 50.0,
            coverage: cov,
            start: "AAAAAAAAAA".into(),
            rcstart: "TTTTTTTTTT".into(),
            end: "CCCCCCCCCC".into(),
            rcend: "GGGGGGGGGG".into(),
            multiplty: None,
        }
    }

    #[test]
    fn empty_returns_na() {
        let records = vec![make_contig(None)];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.mean_str(), "NA");
    }

    #[test]
    fn single_element() {
        let records = vec![make_contig(Some(42.5))];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.mean_str(), "42.50");
    }

    #[test]
    fn two_elements() {
        let records = vec![
            make_contig(Some(10.0)),
            make_contig(Some(20.0)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.mean_str(), "15.00");
    }

    #[test]
    fn four_elements() {
        let records = vec![
            make_contig(Some(1.0)),
            make_contig(Some(2.0)),
            make_contig(Some(3.0)),
            make_contig(Some(4.0)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.mean_str(), "2.50");
    }

    #[test]
    fn filters_none() {
        let records = vec![
            make_contig(None),
            make_contig(Some(42.5)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.mean_str(), "42.50");
    }
}


#[cfg(test)]
mod tests_median_str {
    use super::*;

    fn make_contig(cov: Option<f64>) -> ContigRecord {
        ContigRecord {
            name: "t".into(),
            length: 10,
            gc_content: 50.0,
            coverage: cov,
            start: "AAAAAAAAAA".into(),
            rcstart: "TTTTTTTTTT".into(),
            end: "CCCCCCCCCC".into(),
            rcend: "GGGGGGGGGG".into(),
            multiplty: None,
        }
    }

    #[test]
    fn empty_returns_na() {
        let records = vec![make_contig(None)];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.median_str(), "NA");
    }

    #[test]
    fn single_element() {
        let records = vec![make_contig(Some(42.5))];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.median_str(), "42.50");
    }

    #[test]
    fn odd_count() {
        let records = vec![
            make_contig(Some(3.0)),
            make_contig(Some(1.0)),
            make_contig(Some(2.0)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.median_str(), "2.00");
    }

    #[test]
    fn even_count() {
        let records = vec![
            make_contig(Some(4.0)),
            make_contig(Some(1.0)),
            make_contig(Some(3.0)),
            make_contig(Some(2.0)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.median_str(), "2.50");
    }

    #[test]
    fn filters_none() {
        let records = vec![
            make_contig(None),
            make_contig(Some(42.5)),
        ];
        let summarizer = CovSummarizer::from(&records);
        assert_eq!(summarizer.median_str(), "42.50");
    }
}
