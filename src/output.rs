
use std::fs::File;
use std::path::PathBuf;
use std::collections::HashMap;
use std::io::{Write, BufWriter};

use crate::args::Args;
use crate::contig_record::ContigRecord;
use crate::overlaps::{Terminus, Overlap, OverlapCollection};


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
            if let Err(err) = write_result {
                return Err(format!(
                    "Error: failed to write to file `{:?}`: {}",
                    &log_fpath,
                    err
                ));
            }
        }
    }

    Ok(())
}


fn get_overlap_strings_for_log(contig_collection: &Vec<ContigRecord>,
                               overlap_collection: &OverlapCollection,
                               key: usize) -> Vec<String> {
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


pub fn write_adjacency_table(contig_collection: &Vec<ContigRecord>,
                             overlap_collection: &OverlapCollection,
                             args: &Args) -> Result<(), String> {
    // The function writes adjacency table to output TSV file.
    let sep = "\t";

    // Make path to adjacency table file
    let adj_table_fpath: PathBuf = args.outdir_path.join(
        "combinator_adjacent_contigs.tsv"
    );

    println!("Writing adjacency table to {:?}", adj_table_fpath);

    let file = File::create(&adj_table_fpath);
    if file.is_err() {
        return Err(format!(
            "Error: cannot open file `{:?}` for writing", &adj_table_fpath
        ));
    }
    let mut writer = BufWriter::new(file.unwrap());

    write_adj_table_header(&mut writer, sep, &adj_table_fpath)?;

    for (i, contig) in contig_collection.iter().enumerate() {
        let cov_str: String = contig.coverage.map_or(
            String::from("-"),
            |cov| format!("{:.2}", cov)
        );

        let start_ovl_str = get_overlap_str_for_table(
            overlap_collection,
            contig_collection,
            i,
            TermForTable::Start
        );
        let end_ovl_str = get_overlap_str_for_table(
            overlap_collection,
            contig_collection,
            i,
            TermForTable::End
        );

        let out_values: Vec<String> = vec![
            format!("{}", i + 1),
            format!("{}", contig.name),
            format!("{}", contig.length),
            cov_str,
            format!("{:.2}", contig.gc_content),
            format!("{:.2}", contig.multiplty.unwrap()), // Multiplisity must be Some now
            String::from(""), // Empty column for annotation
            start_ovl_str,
            end_ovl_str,
        ];

        let row_str = out_values.join(sep);
        let write_result = writeln!(writer, "{}", row_str);
        if let Err(err) = write_result {
            return Err(format!(
                "Error: failed to write to file `{:?}`: {}",
                adj_table_fpath,
                err
            ));
        }
    }

    Ok(())
}

enum TermForTable {
    Start,
    End,
}

fn write_adj_table_header(writer: &mut BufWriter<File>,
                          sep: &str,
                          out_fpath: &PathBuf) -> Result<(), String>{
    let col_names: Vec<&str> = vec![
        "#",
        "Contig name",
        "Length",
        "Coverage",
        "GC(%)",
        "Multiplicity",
        "Annotation",
        "Start",
        "End",
    ];
    let header_str = col_names.join(sep);
    let write_result = writeln!(writer, "{}", header_str);
    if let Err(err) = write_result {
        return Err(format!(
            "Error: failed to write to file `{:?}`: {}",
            out_fpath,
            err
        ));
    }

    Ok(())
}

fn get_overlap_str_for_table(overlap_collection: &OverlapCollection,
                             contig_collection:  &Vec<ContigRecord>,
                             key: usize,
                             term: TermForTable) -> String {
    // Function extracts overlaps of `key` contigs associated with `term` terminus
    //   and converts this collection of `Overlap`
    // to string representation  for adjacency table.

    let overlaps = match term {
        TermForTable::Start => {
            get_start_matches(overlap_collection.get(&key))
        },
        TermForTable::End => {
            get_end_matches(overlap_collection.get(&key))
        },
    };

    if overlaps.is_empty() {
        return String::from("-");
    }

    let mut ovl_strings: Vec<String> = Vec::new();
    for ovl in overlaps.iter() {
        if ovl.contig_i != ovl.contig_j {
            // Letter for the first contig of the overlap
            let letter_i = get_match_letter(&ovl.terminus_i);
            // Letter for the second contig of the overlap
            let letter_j = get_match_letter(&ovl.terminus_j);
            // Convert and append
            ovl_strings.push(format!("[{}={}({}); ovl={}]",
                letter_i,
                letter_j,
                contig_collection[ovl.contig_j].name,
                ovl.ovl_len
            ));
        } else {
            // TODO: match it it’s own RC-end?
            ovl_strings.push(format!(
                "[Circle; ovl={}]", ovl.ovl_len
            ));
        }
    }

    ovl_strings.join(" ")
}

fn get_start_matches(overlaps: &[Overlap]) -> Vec<&Overlap> {
    // The function selects "start-associated" overlaps from a collection of overlaps.
    return overlaps.iter().filter(
        |ovl| is_start_match(ovl)
    ).collect();
}

fn get_end_matches(overlaps: &[Overlap]) -> Vec<&Overlap> {
    // The function selects "end-associated" overlaps from a collection of overlaps.
    return overlaps.iter().filter(
        |ovl| is_end_match(ovl)
    ).collect();
}

fn is_start_match(ovl: &Overlap) -> bool {
    // The function returns True if overlap `ovl` is associated with start.
    (
        ovl.terminus_i == Terminus::Start && ovl.terminus_j == Terminus::End
    ) || (
        ovl.terminus_i == Terminus::Start && ovl.terminus_j == Terminus::RcStart
    )
}

fn is_end_match(ovl: &Overlap) -> bool {
    // The function returns True if overlap `ovl` is associated with end.
    (
        ovl.terminus_i == Terminus::End && ovl.terminus_j == Terminus::Start
    ) || (
        ovl.terminus_i == Terminus::End && ovl.terminus_j == Terminus::RcEnd
    )
}

fn get_match_letter(terminus: &Terminus) -> String {
    match terminus {
        Terminus::Start   => String::from("S"),
        Terminus::End     => String::from("E"),
        Terminus::RcStart => String::from("rc_S"),
        Terminus::RcEnd   => String::from("rc_S"),
    }
}
