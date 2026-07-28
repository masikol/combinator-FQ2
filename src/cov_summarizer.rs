
use crate::contig_record::ContigRecord;


pub struct CovSummarizer {
    cov_vector: Vec<f64>,
}

impl CovSummarizer {
    pub fn from(contig_collection: &Vec<ContigRecord>) -> CovSummarizer {
        CovSummarizer {
            cov_vector: contig_collection.iter()
                .map(|contig| contig.coverage)
                .filter_map(|coverage| coverage)
                .collect(),
        }
    }

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
