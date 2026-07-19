
use std::fs::File;
use std::path::PathBuf;
use std::collections::HashMap;
use std::io::{Write, BufWriter};

use crate::args::Args;
use crate::contig_record::ContigRecord;
use crate::overlaps::{Terminus, Overlap, OverlapCollection};


// TODO: remove
// pub fn make_out_prefix(args: &Args) -> PathBuf {
//     let in_basename_prefix = args.input_fpath.file_prefix().unwrap();
//     args.outdir_path.join(in_basename_prefix)
// }



pub fn write_full_log(contig_collection: &Vec<ContigRecord>,
                      overlap_collection: &OverlapCollection,
                      args: &Args) -> Result<(), String>{

    // The function writes full matching log (not only adjacency-associated matches)
    //   to "full-log" file.

    // Make path to full log file
    let log_fpath: PathBuf = args.outdir_path.join(
        "combinator_full_matching_log.txt"
    );

    println!("Writing full matching log to {:?}", log_fpath);

    let file = File::create(&log_fpath);
    if file.is_err() {
        return Err(format!(
            "Error: cannot open file `{:?}` for writing", &log_fpath
        ));
    }

    let mut writer = BufWriter::new(file.unwrap());

    // Write information about discovered adjacency
    for i in 0..contig_collection.len() {
        // Write what matches start of current contig
        let log_strings = get_overlap_strings_for_log(
            contig_collection,
            overlap_collection,
            i
        );
        for s in log_strings {
            let write_result = writeln!(writer, "{}", s);
            if let Err(e) = write_result {
                return Err(format!(
                    "Error: failed to write to file `{:?}`: {}",
                    &log_fpath,
                    e
                ));
            }
        }
    }

    Ok(())
}


fn get_overlap_strings_for_log(contig_collection: &Vec<ContigRecord>,
                               overlap_collection: &OverlapCollection,
                               key: usize) -> Vec<String> {
    // TODO: update comment
    // The function extracts overlaps of `key` contigs associated with `term` terminus,
    //   converts this vector of `Overlap` instances to their string representations
    //   for full log.

    // Extract overlaps for current contig
    let overlaps: &[Overlap] = overlap_collection.get(&key);
    let mut match_strings = Vec::new(); // a list for formatted strings

    if overlaps.is_empty() {
        return match_strings; // no  overlaps found
    }
    let key2word_map = make_key2word_map();

    for ovl in overlaps {
        // If contig does not match itself
        if ovl.contig_i != ovl.contig_j {

            // Word for the first contig of the overlap
            let word1: &String = key2word_map.get(&ovl.terminus_i).unwrap();
            // Word for the second contig of the overlap
            let word2: &String = key2word_map.get(&ovl.terminus_j).unwrap();

            // Convert and append
            let fmt_str = format!(
                "{}: {} matches {} of {} with overlap of {} bp",
                contig_collection[key].name,
                word1,
                word2,
                contig_collection[ovl.contig_j].name,
                ovl.ovl_len
            );
            match_strings.push(fmt_str);
        } else {
            if ovl.terminus_i == Terminus::End && ovl.terminus_j == Terminus::Start {
                // Contig is circular
                let fmt_str = format!(
                    "{}: contig is circular with overlap of {} bp",
                    contig_collection[key].name,
                    ovl.ovl_len
                );
                match_strings.push(fmt_str);
            } else if ovl.terminus_i == Terminus::Start && ovl.terminus_j == Terminus::RcEnd {
                // Start of contig matches it's own reverse-complement end
                let fmt_str = format!(
                    "{}: start is identical to it's own rc-end with overlap of {} bp",
                    contig_collection[key].name,
                    ovl.ovl_len
                );
                match_strings.push(fmt_str);
            }
        }
    }

    match_strings
}


fn make_key2word_map() -> HashMap<Terminus, String> {
    // This dictionary maps `Terminus` to it's "word" representation
    //   for full log.
    let mut key2word_map: HashMap<Terminus, String> = HashMap::new();

    key2word_map.insert(Terminus::Start,   String::from("start"));
    key2word_map.insert(Terminus::RcStart, String::from("rc-start"));
    key2word_map.insert(Terminus::End,     String::from("end"));
    key2word_map.insert(Terminus::RcEnd,   String::from("rc-end"));

    key2word_map
}
