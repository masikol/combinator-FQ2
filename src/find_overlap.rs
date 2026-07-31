
pub fn find_overlap_s2s(seq1: &String,
                        seq2: &String,
                        mink: usize,
                        maxk: usize) -> usize {
    // The function searches for identity between starts of seq1 and seq2.
    // Function regards overlap of length [mink, maxk].
    //
    // Returns 0 if overlap is less than 'mink' and
    //   length of the overlap (which is <= maxk) otherwise.

    let seq1_len: usize = seq1.len();
    let seq2_len: usize = seq2.len();

    if seq1_len < mink || seq2_len < mink {
        return 0;
    }

    let maxk = maxk.min(seq1_len).min(seq2_len);

    let mut overlap: usize = 0;

    let mut i: usize = mink;
    while i <= maxk && &seq1[..i] == &seq2[..i] {
        overlap += 1;
        i += 1;
    }

    if overlap == 0 {
        return 0;
    } else {
        return mink + overlap - 1;
    }
}


pub fn find_overlap_e2s(seq1: &String,
                        seq2: &String,
                        mink: usize,
                        maxk: usize) -> usize {
    // The function searches for identity between end of seq1 and start of seq2.
    // Function regards overlap of length [mink, maxk].
    //
    // Returns 0 if overlap is less than 'mink' and
    //   length of the overlap (which is <= maxk) otherwise.

    let seq1_len: usize = seq1.len();
    let seq2_len: usize = seq2.len();

    if seq1_len < mink || seq2_len < mink {
        return 0;
    }

    let maxk = maxk.min(seq1_len).min(seq2_len);

    let mut overlap: usize = 0;

    for i in mink..(maxk+1) {
        if &seq1[seq1_len-i..] == &seq2[..i] {
            overlap = i;
        }
    }

    return overlap
}


pub fn find_overlap_e2e(seq1: &String,
                        seq2: &String,
                        mink: usize,
                        maxk: usize) -> usize {
    // The function searches for identity between ends of seq1 and seq2.
    // Function regards overlap of length [mink, maxk].
    //
    // Returns 0 if overlap is less than 'mink' and
    //   length of the overlap (which is <= maxk) otherwise.

    let seq1_len: usize = seq1.len();
    let seq2_len: usize = seq2.len();

    if seq1_len < mink || seq2_len < mink {
        return 0;
    }

    let maxk = maxk.min(seq1_len).min(seq2_len);

    let mut overlap: usize = 0;

    let mut i: usize = mink;
    while i <= maxk && seq1[seq1_len-i..] == seq2[seq2_len-i..] {
        overlap += 1;
        i += 1;
    }

    if overlap == 0 {
        return 0;
    } else {
        return mink + overlap - 1;
    }
}


// >>> Tests >>>

#[cfg(test)]
mod tests_s2s {
    use super::*;

    #[test]
    fn no_overlap() {
        let mink = 2;
        let maxk = 4;
        let overlap = find_overlap_s2s(&"aaa".into(), &"gggg".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn overlap_below_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"Attt".into(), &"Aaa".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn overlap_at_mink() {
        let mink = 3;
        let maxk = 3;
        let overlap = find_overlap_s2s(&"AGCt".into(), &"AGCa".into(), mink, maxk);
        assert_eq!(overlap, mink);
    }

    #[test]
    fn overlap_partial() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AGTCa".into(), &"AGTCgg".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn overlap_at_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AGCTAtt".into(), &"AGCTAa".into(), mink, maxk);
        assert_eq!(overlap, maxk);
    }

    #[test]
    fn overlap_beyond_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AGCTAGCt".into(), &"AGCTAGCTaa".into(), mink, maxk);
        assert_eq!(overlap, maxk);
    }

    #[test]
    fn seq1_shorter_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AG".into(), &"AGaaa".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq2_shorter_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AGaaa".into(), &"AG".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq1_shorter_than_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AGAA".into(), &"AGAAaa".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn seq2_shorter_than_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_s2s(&"AGAAaa".into(), &"AGAA".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }
}

#[cfg(test)]
mod tests_e2s {
    use super::*;

    #[test]
    fn no_overlap() {
        let mink = 2;
        let maxk = 4;
        let overlap = find_overlap_e2s(&"aaaaa".into(), &"ggggg".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn overlap_at_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"aaGCT".into(), &"GCTtt".into(), mink, maxk);
        assert_eq!(overlap, mink);
    }

    #[test]
    fn overlap_larger_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"aAGCT".into(), &"AGCTtt".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn mink_fails_larger_succeeds() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"cAGCT".into(), &"AGCTtt".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn overlap_at_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"ttAGCTA".into(), &"AGCTAa".into(), mink, maxk);
        assert_eq!(overlap, maxk);
    }

    #[test]
    fn overlap_beyond_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"tTAGCTA".into(), &"TAGCTAaa".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq1_shorter_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"AG".into(), &"aaaAG".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq2_shorter_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"aaaAG".into(), &"AG".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq1_shorter_than_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"AGAA".into(), &"AGAAaa".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn seq2_shorter_than_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2s(&"aaAGAA".into(), &"AGAA".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }
}

#[cfg(test)]
mod tests_e2e {
    use super::*;

    #[test]
    fn no_overlap() {
        let mink = 2;
        let maxk = 4;
        let overlap = find_overlap_e2e(&"aaaa".into(), &"ggg".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn overlap_at_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"ttAGC".into(), &"aAGC".into(), mink, maxk);
        assert_eq!(overlap, mink);
    }

    #[test]
    fn overlap_partial() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"aAGCT".into(), &"ttAGCT".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn overlap_at_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"aTAGCT".into(), &"ttTAGCT".into(), mink, maxk);
        assert_eq!(overlap, maxk);
    }

    #[test]
    fn overlap_beyond_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"tCTAGCT".into(), &"aaCTAGCT".into(), mink, maxk);
        assert_eq!(overlap, maxk);
    }

    #[test]
    fn seq1_shorter_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"AG".into(), &"aaaAG".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq2_shorter_than_mink() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"aaaAG".into(), &"AG".into(), mink, maxk);
        assert_eq!(overlap, 0);
    }

    #[test]
    fn seq1_shorter_than_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"AGAA".into(), &"aaAGAA".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn seq2_shorter_than_maxk() {
        let mink = 3;
        let maxk = 5;
        let overlap = find_overlap_e2e(&"aaAGAA".into(), &"AGAA".into(), mink, maxk);
        assert_eq!(overlap, 4);
    }
}
