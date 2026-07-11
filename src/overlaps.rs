
use std::io::{self, Write};
use std::collections::HashMap;

use crate::args::Args;
use crate::contig_record::ContigRecord;
use crate::find_overlap::{
    find_overlap_s2s,
    find_overlap_e2s,
    find_overlap_e2e,
};


#[derive(Debug)]
#[derive(PartialEq)]
enum Terminus {
    Start,
    End,
    RcStart,
    RcEnd,
}

type ContigIdx = usize;

#[derive(Debug)]
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


pub struct OverlapCollection {
    // TODO: remove pub
    pub collection: HashMap<ContigIdx, Vec<Overlap>>,
}


impl OverlapCollection {
    /// Creates a new empty OverlapCollection
    pub fn new() -> Self {
        OverlapCollection {
            collection: HashMap::new(),
        }
    }

    // TODO: use
    // /// Returns a reference to the list of overlaps for a given contig key.
    // /// Returns an empty slice if the key doesn't exist.
    // pub fn get(&self, key: &ContigIdx) -> &[Overlap] {
    //     self.collection.get(key).map_or(&[], |v: &Vec<Overlap>| v.as_slice())
    // }

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

// TODO: remove
// impl fmt::Display for OverlapCollection {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         write!(f, "{:?}", self.collection)
//     }
// }

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
        let ovl_len = find_overlap_e2s(&cont_i.end, &cont_i.start, mink, maxk);
        if ovl_len != 0 && ovl_len != cont_i.length {
            overlaps.add(i, Overlap::new(i, Terminus::End,   i, Terminus::Start, ovl_len));
            overlaps.add(i, Overlap::new(i, Terminus::Start, i, Terminus::End,   ovl_len));
        }

        // === Compare start of the current conitg to rc-end of the current contig ===
        let ovl_len = find_overlap_s2s(&cont_i.end, &cont_i.rcstart, mink, maxk);
        if ovl_len != 0 {
            overlaps.add(i, Overlap::new(i, Terminus::Start, i, Terminus::RcEnd, ovl_len));
            overlaps.add(i, Overlap::new(i, Terminus::RcEnd, i, Terminus::Start, ovl_len));
        }

        // |=== Compare i-th contig to contigs from i+1 to N ===|
        // We do it in order not to compare pairs of contigs more than one time
        for j in i+1..num_contigs {

            let cont_j = &contigs[j];

            // === Compare i-th start to j-th end ===
            let ovl_len = find_overlap_e2s(&cont_j.end, &cont_i.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::End,   ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End,   i, Terminus::Start, ovl_len));
            }

            // === Compare i-th end to j-th start ===
            let ovl_len = find_overlap_e2s(&cont_i.end, &cont_j.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End,   j, Terminus::Start, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::End,   ovl_len));
            }

            // === Compare i-th start to reverse-complement j-th start ===
            let ovl_len = find_overlap_e2s(&cont_j.rcstart, &cont_i.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::RcStart, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::RcStart, ovl_len));
            }

            // === Compare i-th end to reverse-complement j-th end ===
            let ovl_len = find_overlap_e2s(&cont_i.end, &cont_j.rcend, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End, j, Terminus::RcEnd, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End, i, Terminus::RcEnd, ovl_len));
            }

            // === Compare i-th start to j-th start ===
            let ovl_len = find_overlap_s2s(&cont_i.start, &cont_j.start, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::Start, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::Start, i, Terminus::Start, ovl_len));
            }

            // === Compare i-th end to j-th end ===
            let ovl_len = find_overlap_e2e(&cont_i.end, &cont_j.end, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End, j, Terminus::End, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::End, i, Terminus::End, ovl_len));
            }

            // === Compare i-th start to reverse-complement j-th end ===
            let ovl_len = find_overlap_s2s(&cont_i.start, &cont_j.rcend, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::Start, j, Terminus::RcEnd, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::RcEnd, i, Terminus::Start, ovl_len));
            }

            // === Compare i-th end to reverse-complement j-th start ===
            let ovl_len = find_overlap_e2e(&cont_i.end, &cont_j.rcstart, mink, maxk);
            if ovl_len != 0 {
                overlaps.add(i, Overlap::new(i, Terminus::End,     j, Terminus::RcStart, ovl_len));
                overlaps.add(j, Overlap::new(j, Terminus::RcStart, i, Terminus::End,     ovl_len));
            }
        }

        // Update makeshift status bar
        print!("\r{}/{}", i+1, num_contigs);
        io::stdout().flush().unwrap();
    }

    overlaps
}
