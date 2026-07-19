
use std::io::{self, Write};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use crate::args::Args;
use crate::contig_record::ContigRecord;
use crate::find_overlap::{
    find_overlap_s2s,
    find_overlap_e2s,
    find_overlap_e2e,
};


#[derive(Debug, PartialEq, Eq, Hash)]
enum Terminus {
    Start,
    End,
    RcStart,
    RcEnd,
}

type ContigIdx = usize;

#[derive(Debug, Eq)]
pub struct Overlap {
    contig_i: ContigIdx,
    terminus_i: Terminus,
    contig_j: ContigIdx,
    terminus_j: Terminus,
    ovl_len: usize,
}

impl Overlap {
    /// A constructor for concise Overlap instatiation
    fn new(contig_i: ContigIdx,
           terminus_i: Terminus,
           contig_j: ContigIdx,
           terminus_j: Terminus,
           ovl_len: usize) -> Self {
        Overlap {
            contig_i: contig_i,
            terminus_i: terminus_i,
            contig_j: contig_j,
            terminus_j: terminus_j,
            ovl_len: ovl_len,
        }
    }

    pub fn is_start_match(&self) -> bool {
        (
            self.terminus_i == Terminus::Start
         && self.terminus_j == Terminus::End
        ) || (
            self.terminus_i == Terminus::Start
         && self.terminus_j == Terminus::RcStart
        )
    }

    pub fn is_end_match(&self) -> bool {
        (
            self.terminus_i == Terminus::End
         && self.terminus_j == Terminus::Start
        ) || (
            self.terminus_i == Terminus::End
         && self.terminus_j == Terminus::RcEnd
        )
    }
}

impl PartialEq for Overlap {
    fn eq(&self, other: &Self) -> bool {
        if self.contig_i != other.contig_i {
            return false;
        }
        if self.contig_j != other.contig_j {
            return false;
        }

        if self.terminus_i != other.terminus_i {
            return false;
        }
        if self.terminus_j != other.terminus_j {
            return false;
        }

        if self.ovl_len != other.ovl_len {
            return false;
        }

        true
    }
}

impl Hash for Overlap {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.contig_i.hash(state);
        self.terminus_i.hash(state);
        self.contig_j.hash(state);
        self.terminus_j.hash(state);
        self.ovl_len.hash(state);
    }
}


pub struct OverlapCollection {
    // TODO: production: remove pub
    pub collection: HashMap<ContigIdx, Vec<Overlap>>,
}


impl OverlapCollection {
    /// Creates a new empty OverlapCollection
    pub fn new() -> Self {
        OverlapCollection {
            collection: HashMap::new(),
        }
    }

    // TODO: why &ContigIdx?
    /// Returns a reference to the list of overlaps for a given contig key.
    /// Returns an empty slice if the key doesn't exist.
    pub fn get(&self, key: &ContigIdx) -> &[Overlap] {
        self.collection.get(key).map_or(&[], |v: &Vec<Overlap>| v.as_slice())
    }

    // TODO: use
    // /// Returns the number of contigs in the collection.
    // pub fn len(&self) -> usize {
    //     self.collection.len()
    // }

    // TODO: use
    // /// Returns true if the collection is empty.
    // pub fn is_empty(&self) -> bool {
    //     self.collection.is_empty()
    // }

    /// Adds an overlap to the collection for the given key.
    /// Creates a new list if the key doesn't exist yet.
    pub fn add(&mut self, key: ContigIdx, overlap: Overlap) {
        self.collection.entry(key).or_insert_with(Vec::new).push(overlap);
    }
}

pub fn detect_adjacent_contigs(contigs: &Vec<ContigRecord>,
                               args: &Args) -> OverlapCollection {
    let mink = args.mink;
    let maxk = args.maxk;
    let num_contigs: ContigIdx = contigs.len();

    let mut overlaps = OverlapCollection::new();

    // Init makeshift status bar
    print!("\n{}/{}", 0, num_contigs);
    io::stdout().flush().unwrap();

    for i in 0..num_contigs {

        let cont_i = &contigs[i];

        // === Compare start of the current contig to end of the current contig ===
        // Match variant 1
        let ovl_len = find_overlap_e2s(&cont_i.end, &cont_i.start, mink, maxk);
        if ovl_len != 0 && ovl_len != cont_i.length {
            overlaps.add(i, Overlap::new(i, Terminus::End,   i, Terminus::Start, ovl_len));
            overlaps.add(i, Overlap::new(i, Terminus::Start, i, Terminus::End,   ovl_len));
        }

        // === Compare start of the current conitg to rc-end of the current contig ===
        // Match variant 2
        let ovl_len = find_overlap_s2s(&cont_i.start, &cont_i.rcend, mink, maxk);
        if ovl_len != 0 {
            overlaps.add(i, Overlap::new(i, Terminus::Start, i, Terminus::RcEnd, ovl_len));
            overlaps.add(i, Overlap::new(i, Terminus::RcEnd, i, Terminus::Start, ovl_len));
        }

        // |=== Compare i-th contig to contigs from i+1 to N ===|
        // We do it in order not to compare pairs of contigs more than one time
        for j in i+1..num_contigs {

            let cont_j = &contigs[j];

            // === Compare i-th start to j-th end ===
            // Match variant 3
            let ovl_len = find_overlap_e2s(&cont_j.end, &cont_i.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::End,   ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End,   i, Terminus::Start, ovl_len));
            }

            // === Compare i-th end to j-th start ===
            // Match variant 4
            let ovl_len = find_overlap_e2s(&cont_i.end, &cont_j.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End,   j, Terminus::Start, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::End,   ovl_len));
            }

            // === Compare i-th start to reverse-complement j-th start ===
            // Match variant 5
            let ovl_len = find_overlap_e2s(&cont_j.rcstart, &cont_i.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::RcStart, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::RcStart, ovl_len));
            }

            // === Compare i-th end to reverse-complement j-th end ===
            // Match variant 6
            let ovl_len = find_overlap_e2s(&cont_i.end, &cont_j.rcend, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End, j, Terminus::RcEnd, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End, i, Terminus::RcEnd, ovl_len));
            }

            // === Compare i-th start to j-th start ===
            // Match variant 7
            let ovl_len = find_overlap_s2s(&cont_i.start, &cont_j.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::Start, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::Start, ovl_len));
            }

            // === Compare i-th end to j-th end ===
            // Match variant 8
            let ovl_len = find_overlap_e2e(&cont_i.end, &cont_j.end, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End, j, Terminus::End, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End, i, Terminus::End, ovl_len));
            }

            // === Compare i-th start to reverse-complement j-th end ===
            // Match variant 9
            let ovl_len = find_overlap_s2s(&cont_i.start, &cont_j.rcend, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::RcEnd,   ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End,   i, Terminus::RcStart, ovl_len));
            }

            // === Compare i-th end to reverse-complement j-th start ===
            // Match variant 10
            let ovl_len = find_overlap_e2e(&cont_i.end, &cont_j.rcstart, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End,   j, Terminus::RcStart, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::RcEnd,   ovl_len));
            }
        }

        // Update makeshift status bar
        print!("\r{}/{}", i+1, num_contigs);
        io::stdout().flush().unwrap();
    }

    println!("");

    overlaps
}


#[cfg(test)]
mod tests_detect_adjacent_contigs {
    use std::path::PathBuf;
    use std::collections::HashSet;

    use crate::args::Args;
    use crate::fasta_reader::{FastaReader, SeqRecord};
    use crate::contig_record::ContigRecord;
    use super::*;

    fn test_path_variants_3_4() -> PathBuf {
        PathBuf::from("test_data/overlaps/test_contigs_variants_3-4.fasta")
    }

    fn test_path_variant_1() -> PathBuf {
        PathBuf::from("test_data/overlaps/test_contigs_variant_1.fasta")
    }

    fn test_path_variant_2() -> PathBuf {
        PathBuf::from("test_data/overlaps/test_contigs_variant_2.fasta")
    }

    fn test_path_variants_5_6() -> PathBuf {
        PathBuf::from("test_data/overlaps/test_contigs_variants_5-6.fasta")
    }

    fn test_path_variants_7_8() -> PathBuf {
        PathBuf::from("test_data/overlaps/test_contigs_variants_7-8.fasta")
    }

    fn test_path_rc_variants_9_10() -> PathBuf {
        PathBuf::from("test_data/overlaps/test_contigs_variants_9-10.fasta")
    }

    fn read_test_contigs(file_path: &PathBuf, maxk: usize) -> Vec<ContigRecord> {
        let reader = FastaReader::open(file_path).unwrap();
        reader
            .filter_map(|result| result.ok())
            .map(|seq| ContigRecord::from(seq, maxk).unwrap())
            .collect()
    }

    fn make_args(file_path: &PathBuf, mink: usize, maxk: usize) -> Args {
        Args {
            mink,
            maxk,
            input_fpath: file_path.clone(),
            outdir_path: PathBuf::from(""),
        }
    }

    fn assert_overlap_set_eq(got: &[Overlap], expected: &[Overlap]) {
        let got_set: HashSet<&Overlap> = got.iter().collect();
        let expected_set: HashSet<&Overlap> = expected.iter().collect();
        assert_eq!(got_set, expected_set);
    }


    // ===== Tests =====

    #[test]
    fn test_detect_adjacent_contigs_variant_1() {
        // Tests match variant 1
        let file_path = test_path_variant_1();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 16, 25);

        let idx: ContigIdx = 0;

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let expected: Vec<Overlap> = vec![
            Overlap::new(idx, Terminus::Start, idx, Terminus::End,   21),
            Overlap::new(idx, Terminus::End,   idx, Terminus::Start, 21),
        ];
        assert_overlap_set_eq(&overlaps.collection[&idx], &expected);
    }

    #[test]
    fn test_detect_adjacent_contigs_variant_2() {
        // Tests match variant 2
        let file_path = test_path_variant_2();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 16, 25);

        let idx: ContigIdx = 0;

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let expected: Vec<Overlap> = vec![
            Overlap::new(idx, Terminus::Start, idx, Terminus::RcEnd, 21),
            Overlap::new(idx, Terminus::RcEnd, idx, Terminus::Start, 21),
        ];
        assert_overlap_set_eq(&overlaps.collection[&idx], &expected);
    }

    #[test]
    fn test_detect_adjacent_contigs_variants_3_4() {
        // Tests match variant 3 and variant 4
        let file_path = test_path_variants_3_4();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 16, 25);

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let n1: ContigIdx = 0;
        let n2: ContigIdx = 1;

        // Node 1
        let expected: Vec<Overlap> = vec![
            Overlap::new(n1, Terminus::Start, n2, Terminus::End,   18),
            Overlap::new(n1, Terminus::End,   n2, Terminus::Start, 16),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n1], &expected);

        // Node 2
        let expected: Vec<Overlap> = vec![
            Overlap::new(n2, Terminus::Start, n1, Terminus::End,   16),
            Overlap::new(n2, Terminus::End,   n1, Terminus::Start, 18),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n2], &expected);
    }

    #[test]
    fn test_detect_adjacent_contigs_variants_5_6() {
        // Tests match variant 5 and variant 6
        let file_path = test_path_variants_5_6();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 16, 25);

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let n1: ContigIdx = 0;
        let n2: ContigIdx = 1;

        // Node 1
        let expected: Vec<Overlap> = vec![
            Overlap::new(n1, Terminus::Start, n2, Terminus::RcStart, 18),
            Overlap::new(n1, Terminus::End,   n2, Terminus::RcEnd,   16),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n1], &expected);

        // Node 2
        let expected: Vec<Overlap> = vec![
            Overlap::new(n2, Terminus::End,   n1, Terminus::RcEnd,   16),
            Overlap::new(n2, Terminus::Start, n1, Terminus::RcStart, 18),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n2], &expected);
    }

    #[test]
    fn test_detect_adjacent_contigs_variants_7_8() {
        // Tests match variant 7 and variant 8
        let file_path = test_path_variants_7_8();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 16, 25);

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let n1: ContigIdx = 0;
        let n2: ContigIdx = 1;

        // Node 1
        let expected: Vec<Overlap> = vec![
            Overlap::new(n1, Terminus::Start, n2, Terminus::Start, 18),
            Overlap::new(n1, Terminus::End,   n2, Terminus::End,   16),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n1], &expected);

        // Node 2
        let expected: Vec<Overlap> = vec![
            Overlap::new(n2, Terminus::Start, n1, Terminus::Start, 18),
            Overlap::new(n2, Terminus::End,   n1, Terminus::End,   16),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n2], &expected);
    }

    #[test]
    fn test_detect_adjacent_contigs_variants_9_10() {
        // Tests match variant 9 and variant 10
        let file_path = test_path_rc_variants_9_10();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 16, 25);

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let n1: ContigIdx = 0;
        let n2: ContigIdx = 1;

        // Node 1
        let expected: Vec<Overlap> = vec![
            Overlap::new(n1, Terminus::Start, n2, Terminus::RcEnd,   18),
            Overlap::new(n1, Terminus::End,   n2, Terminus::RcStart, 16),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n1], &expected);

        // Node 2
        let expected: Vec<Overlap> = vec![
            Overlap::new(n2, Terminus::Start, n1, Terminus::RcEnd,   16),
            Overlap::new(n2, Terminus::End,   n1, Terminus::RcStart, 18),
        ];
        assert_overlap_set_eq(&overlaps.collection[&n2], &expected);
    }

    #[test]
    fn test_no_overlap_contigs() {
        // Single contig with four cleanly distinct terminii
        let contigs = vec![
            ContigRecord {
                name: "a".to_string(),
                length: 25,
                gc_content: 50.0,
                coverage: None,
                start:     "AAAAAAAAAAAAAAAAAAAAAAAAA".to_string(),
                rcstart:   "TTTTTTTTTTTTTTTTTTTTTTTTT".to_string(),
                end:       "CCCCCCCCCCCCCCCCCCCCCCCCC".to_string(),
                rcend:     "GGGGGGGGGGGGGGGGGGGGGGGGG".to_string(),
                multiplty: None,
            },
            ContigRecord {
                name: "b".to_string(),
                length: 25,
                gc_content: 50.0,
                coverage: None,
                start:     "HMMMMMMMMMMMMMMMMMMMMMMMM".to_string(),
                rcstart:   "HMMMMMMMMMMMMMMMMMMMMMMMM".to_string(),
                end:       "DWWWWWWWWWWWWWWWWWWWWWWWWW".to_string(),
                rcend:     "DWWWWWWWWWWWWWWWWWWWWWWWWW".to_string(),
                multiplty: None,
            },
        ];
        let args = Args {
            mink: 3,
            maxk: 5,
            input_fpath: PathBuf::from(""),
            outdir_path: PathBuf::from(""),
        };
        let overlaps = detect_adjacent_contigs(&contigs, &args);

        assert!(overlaps.collection.is_empty());
    }

    #[test]
    fn test_self_overlap_end_to_start() {
        // Build ContigRecord manually with identical start and end
        let contigs = vec![
            ContigRecord {
                name: "circular_contig".to_string(),
                length: 50,
                gc_content: 50.0,
                coverage: None,
                start:     "ACGATCGATCGATCGTAGCTAGCA".to_string(),
                rcstart:   "TGCTAGCTACGATCGATCGATCGT".to_string(),
                end:       "ACGATCGATCGATCGTAGCTAGCA".to_string(),
                rcend:     "TGCTAGCTACGATCGATCGATCGT".to_string(),
                multiplty: None,
            },
        ];
        let args = Args {
            mink: 3,
            maxk: 25,
            input_fpath: PathBuf::from(""),
            outdir_path: PathBuf::from(""),
        };

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        let expected: Vec<Overlap> = vec![
            Overlap::new(0, Terminus::Start, 0, Terminus::End,   24),
            Overlap::new(0, Terminus::End,   0, Terminus::Start, 24),
        ];
        assert_overlap_set_eq(&overlaps.collection[&0], &expected);
    }

    #[test]
    fn test_mink_equals_maxk() {
        let file_path = test_path_variants_3_4();
        let contigs = read_test_contigs(&file_path, 25);
        let args = make_args(&file_path, 25, 25);

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        // All known overlaps are <25 (16, 18, 21), so no overlap should be detected
        assert!(overlaps.collection.is_empty());
    }

    #[test]
    fn test_contig_shorter_than_maxk() {
        // Build ContigRecords manually with terminii shorter than args.maxk
        let contigs = vec![
            ContigRecord {
                name: "a".to_string(),
                length: 8,
                gc_content: 50.0,
                coverage: None,
                start: "AAAA".to_string(),
                rcstart: "TTTT".to_string(),
                end: "CCCC".to_string(),
                rcend: "GGGG".to_string(),
                multiplty: None,
            },
        ];
        let args = Args {
            mink: 100,
            maxk: 100,
            input_fpath: PathBuf::from(""),
            outdir_path: PathBuf::from(""),
        };

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        // find_overlap functions cap maxk to min(seq_len1, seq_len2) and
        // return 0 if any seq is shorter than mink (100 > 4). Empty result.
        assert!(overlaps.collection.is_empty());
    }

    #[test]
    fn test_empty_contigs() {
        let contigs: Vec<ContigRecord> = vec![];
        let args = Args {
            mink: 16,
            maxk: 25,
            input_fpath: PathBuf::from(""),
            outdir_path: PathBuf::from(""),
        };

        let overlaps = detect_adjacent_contigs(&contigs, &args);

        assert!(overlaps.collection.is_empty());
    }
}


#[cfg(test)]
mod tests_overlap_collection_get {
    use super::*;

    fn make_overlap(contig_i: ContigIdx, contig_j: ContigIdx, ovl_len: usize) -> Overlap {
        Overlap::new(contig_i, Terminus::End, contig_j, Terminus::Start, ovl_len)
    }

    #[test]
    fn get_returns_empty_slice_for_nonexistent_key() {
        let col = OverlapCollection::new();
        let result = col.get(&0);
        assert!(result.is_empty());
    }

    #[test]
    fn get_returns_overlaps_for_existing_key() {
        let mut col = OverlapCollection::new();
        let ovl1 = make_overlap(0, 1, 10);
        let ovl2 = make_overlap(0, 2, 20);
        col.add(0, ovl1);
        col.add(0, ovl2);

        let result = col.get(&0);
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], make_overlap(0, 1, 10));
        assert_eq!(result[1], make_overlap(0, 2, 20));
    }

    #[test]
    fn get_returns_all_overlaps_multiple_keys() {
        let mut col = OverlapCollection::new();
        let ovl_a = make_overlap(0, 1, 10);
        let ovl_b = make_overlap(1, 0, 15);
        col.add(0, ovl_a);
        col.add(1, ovl_b);

        let result_0 = col.get(&0);
        assert_eq!(result_0.len(), 1);
        assert_eq!(result_0[0], make_overlap(0, 1, 10));

        let result_1 = col.get(&1);
        assert_eq!(result_1.len(), 1);
        assert_eq!(result_1[0], make_overlap(1, 0, 15));
    }

    #[test]
    fn get_does_not_affect_collection() {
        let mut col = OverlapCollection::new();
        let ovl = make_overlap(0, 1, 10);
        col.add(0, ovl);

        let _result = col.get(&0);
        assert_eq!(col.get(&0).len(), 1);
    }
}
