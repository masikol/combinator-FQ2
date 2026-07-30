
use std::fs::File;
use std::path::PathBuf;
use std::collections::HashMap;
use std::io::{Write, BufWriter};

use crate::args::Args;
use crate::contig_record::ContigRecord;
use crate::cov_summarizer::CovSummarizer;
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

pub fn write_summary(contig_collection: &Vec<ContigRecord>,
                     overlap_collection: &OverlapCollection,
                     args: &Args) -> Result<(), String> {
    // Make path to adjacency table file
    let summary_fpath: PathBuf = args.outdir_path.join(
        "combinator_summary_FQ.txt"
    );

    println!("Writing summary to {:?}", summary_fpath);

    let file = File::create(&summary_fpath);
    if file.is_err() {
        return Err(format!(
            "Error: cannot open file `{:?}` for writing", &summary_fpath
        ));
    }
    let mut writer = BufWriter::new(file.unwrap());

    let out_lines = make_summary_lines(
        contig_collection,
        overlap_collection,
        args
    );

    println!("");
    for line in out_lines {
        println!("{}", line);
        let write_result = writeln!(writer, "{}", line);
        if let Err(err) = write_result {
            return Err(format!(
                "Error: failed to write to file `{:?}`: {}",
                summary_fpath,
                err
            ));
        }
    }

    Ok(())
}

fn make_summary_lines(contig_collection: &Vec<ContigRecord>,
                      overlap_collection: &OverlapCollection,
                      args: &Args) -> Vec<String> {
    let mut out_lines: Vec<String> = Vec::with_capacity(10);

    // Summary with some statistics
    out_lines.push(
        String::from("=== Summary ===")
    );

    // Path to input file
    out_lines.push(format!(
        "Input file: {:?}", args.input_fpath
    ));

    // Number of contigs processed:
    out_lines.push(format!(
        "{} contigs.", contig_collection.len()
    ));

    // Sum of contigs' lengths
    out_lines.push(format!(
        "Sum of contig lengths: {} bp",
        calc_sum_contig_lengths(contig_collection)
    ));

    // Expected length of the genome
    out_lines.push(format!(
        "Expected genome size: {} bp",
        calc_exp_genome_size(contig_collection, overlap_collection)
    ));

    // Create a coverage summarizer
    let cov_summarizer = CovSummarizer::from(contig_collection);

    // Min coverage
    out_lines.push(format!(
        "Min coverage: {}",
        cov_summarizer.min_str()
    ));

    // Max coverage
    out_lines.push(format!(
        "Max coverage: {}",
        cov_summarizer.max_str()
    ));

    // Mean coverage
    out_lines.push(format!(
        "Mean coverage: {}",
        cov_summarizer.mean_str()
    ));

    // Median coverage
    out_lines.push(format!(
        "Median coverage: {}",
        cov_summarizer.median_str()
    ));

    // LQ coefficient
    out_lines.push(format!(
        "LQ coefficient: {:.2}",
        calc_lq_coef(contig_collection, overlap_collection)
    ));
    out_lines.push(
        "-".repeat(20)
    );

    out_lines
}

fn calc_sum_contig_lengths(contig_collection: &Vec<ContigRecord>) -> usize {
    contig_collection.iter().map(
        |contig| contig.length
    ).sum()
}

fn calc_exp_genome_size(contig_collection: &Vec<ContigRecord>,
                        overlap_collection: &OverlapCollection) -> usize {
    let sum_overlap_len = calc_sum_overlap_len(
        contig_collection,
        overlap_collection
    );

    contig_collection.iter().map(
        |contig| contig.length * contig.multiplty.unwrap().round() as usize
    ).sum::<usize>() - sum_overlap_len
}

fn calc_sum_overlap_len(contig_collection: &Vec<ContigRecord>,
                        overlap_collection: &OverlapCollection) -> usize {
    let mut total_overlap_len: usize = 0;

    for (i, contig) in contig_collection.iter().enumerate() {
        let start_ovls = get_start_matches(overlap_collection.get(&i));
        let end_ovls   = get_end_matches(overlap_collection.get(&i));
        let multiplty = contig.multiplty.unwrap().round() as usize;
        let not_already_counted = |ovl: &&Overlap| ovl.contig_j >= i;

        let mut ovls_to_add: Vec<&Overlap>;

        for ovl_vector in [start_ovls, end_ovls] {
            if ovl_vector.len() <= multiplty {
                // No extra overlaps.
                // We will just add lengths of overlaps to `total_overlap_len`.
                ovls_to_add = ovl_vector.into_iter()
                    .filter(not_already_counted)
                    .collect();
            } else {
                // Some extra overlaps discovered.
                // We will consider only M longest overlaps,
                //   where M is contig's multiplicity.
                ovls_to_add = ovl_vector.into_iter()
                    .filter(not_already_counted)
                    .collect();
                ovls_to_add.sort_by_key(|overlap| overlap.ovl_len);
                ovls_to_add.reverse();
                ovls_to_add.truncate(multiplty);
            }
            total_overlap_len += ovls_to_add.iter().map(
                |ovl| ovl.ovl_len
            ).sum::<usize>();
        }
    }

    total_overlap_len
}

fn calc_lq_coef(contig_collection: &Vec<ContigRecord>,
                overlap_collection: &OverlapCollection) -> f64 {
    // The function calculates LQ-coefficient for given contigs.

    // Number of termini of a contig
    let num_contig_termini: usize = 2;
    // Total number of dead ends taking account of multiplicity
    let mut total_dead_ends: usize = 0;

    let num_cotigs: usize = contig_collection.len();

    for i in 0..num_cotigs {
        // Count overlaps associated with start
        let num_start_overlaps: usize = overlap_collection.get(&i)
            .iter()
            .filter(|ovl| is_start_match(ovl)) // TODO: or simply is_start_match?
            .count();
        // Count overlaps associated with end
        let num_end_overlaps: usize = overlap_collection.get(&i)
            .iter()
            .filter(|ovl| is_end_match(ovl)) // TODO: or simply is_start_match?
            .count();

        let num_dead_ends: usize = 2
            - num_start_overlaps.min(1)
            - num_end_overlaps.min(1);
        total_dead_ends += num_dead_ends;
    }

    // Total number of termini taking account of multiplicity
    let num_total_termini: usize = num_contig_termini * num_cotigs;

    // Calculate and return LQ coefficient
    (
        1.0
        - total_dead_ends as f64 / num_total_termini as f64
    ) * 100.0
}


#[cfg(test)]
mod tests_get_start_matches {
    use super::*;

    fn make_ovl(ti: Terminus, tj: Terminus) -> Overlap {
        Overlap {
            contig_i: 0,
            terminus_i: ti,
            contig_j: 1,
            terminus_j: tj,
            ovl_len: 10,
        }
    }

    #[test]
    fn returns_only_start_matches() {
        let overlaps: Vec<Overlap> = vec![
            make_ovl(Terminus::Start, Terminus::End),
            make_ovl(Terminus::Start, Terminus::RcStart),
            make_ovl(Terminus::End, Terminus::Start),
            make_ovl(Terminus::End, Terminus::RcEnd),
        ];
        let result = get_start_matches(&overlaps);

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].terminus_i, Terminus::Start);
        assert_eq!(result[0].terminus_j, Terminus::End);
        assert_eq!(result[1].terminus_i, Terminus::Start);
        assert_eq!(result[1].terminus_j, Terminus::RcStart);
    }

    #[test]
    fn returns_empty_for_empty_input() {
        let overlaps: Vec<Overlap> = vec![];
        let result = get_start_matches(&overlaps);

        assert!(result.is_empty());
    }

    #[test]
    fn returns_empty_when_no_start_matches() {
        let overlaps: Vec<Overlap> = vec![
            make_ovl(Terminus::End, Terminus::Start),
            make_ovl(Terminus::End, Terminus::RcEnd),
        ];
        let result = get_start_matches(&overlaps);

        assert!(result.is_empty());
    }

    #[test]
    fn returns_all_when_all_are_start_matches() {
        let overlaps: Vec<Overlap> = vec![
            make_ovl(Terminus::Start, Terminus::End),
            make_ovl(Terminus::Start, Terminus::RcStart),
        ];
        let result = get_start_matches(&overlaps);

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].terminus_i, Terminus::Start);
        assert_eq!(result[1].terminus_i, Terminus::Start);
    }
}


#[cfg(test)]
mod tests_get_end_matches {
    use super::*;

    fn make_ovl(ti: Terminus, tj: Terminus) -> Overlap {
        Overlap {
            contig_i: 0,
            terminus_i: ti,
            contig_j: 1,
            terminus_j: tj,
            ovl_len: 10,
        }
    }

    #[test]
    fn returns_only_end_matches() {
        let overlaps: Vec<Overlap> = vec![
            make_ovl(Terminus::End, Terminus::Start),
            make_ovl(Terminus::End, Terminus::RcEnd),
            make_ovl(Terminus::Start, Terminus::End),
            make_ovl(Terminus::Start, Terminus::RcStart),
        ];
        let result = get_end_matches(&overlaps);

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].terminus_i, Terminus::End);
        assert_eq!(result[0].terminus_j, Terminus::Start);
        assert_eq!(result[1].terminus_i, Terminus::End);
        assert_eq!(result[1].terminus_j, Terminus::RcEnd);
    }

    #[test]
    fn returns_empty_for_empty_input() {
        let overlaps: Vec<Overlap> = vec![];
        let result = get_end_matches(&overlaps);

        assert!(result.is_empty());
    }

    #[test]
    fn returns_empty_when_no_end_matches() {
        let overlaps: Vec<Overlap> = vec![
            make_ovl(Terminus::Start, Terminus::End),
            make_ovl(Terminus::Start, Terminus::RcStart),
        ];
        let result = get_end_matches(&overlaps);

        assert!(result.is_empty());
    }

    #[test]
    fn returns_all_when_all_are_end_matches() {
        let overlaps: Vec<Overlap> = vec![
            make_ovl(Terminus::End, Terminus::Start),
            make_ovl(Terminus::End, Terminus::RcEnd),
        ];
        let result = get_end_matches(&overlaps);

        assert_eq!(result.len(), 2);
        assert_eq!(result[0].terminus_i, Terminus::End);
        assert_eq!(result[1].terminus_i, Terminus::End);
    }
}


#[cfg(test)]
mod tests_calc_lq_coef {
    use super::*;

    fn make_contig() -> ContigRecord {
        ContigRecord {
            name: "t".into(),
            length: 10,
            gc_content: 50.0,
            coverage: None,
            start: "AAAAAAAAAA".into(),
            rcstart: "TTTTTTTTTT".into(),
            end: "CCCCCCCCCC".into(),
            rcend: "GGGGGGGGGG".into(),
            multiplty: None,
        }
    }

    fn start_match(contig_i: usize, contig_j: usize) -> Overlap {
        Overlap {
            contig_i: contig_i,
            terminus_i: Terminus::Start,
            contig_j: contig_j,
            terminus_j: Terminus::End,
            ovl_len: 10,
        }
    }

    fn end_match(contig_i: usize, contig_j: usize) -> Overlap {
        Overlap {
            contig_i: contig_i,
            terminus_i: Terminus::End,
            contig_j: contig_j,
            terminus_j: Terminus::Start,
            ovl_len: 10,
        }
    }

    #[test]
    fn no_overlaps_single_contig() {
        let contigs = vec![make_contig()];
        let overlaps = OverlapCollection::new();

        let result = calc_lq_coef(&contigs, &overlaps);
        let expected = 0.00;

        assert!((result - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn only_start_matched() {
        let contigs = vec![make_contig()];
        let mut overlaps = OverlapCollection::new();
        overlaps.add(0, start_match(0, 1));

        let result = calc_lq_coef(&contigs, &overlaps);
        let expected = 50.00;

        assert!((result - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn only_end_matched() {
        let contigs = vec![make_contig()];
        let mut overlaps = OverlapCollection::new();
        overlaps.add(0, end_match(0, 1));

        let result = calc_lq_coef(&contigs, &overlaps);
        let expected = 50.00;

        assert!((result - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn multiple_contigs_no_overlaps() {
        let contigs = vec![make_contig(), make_contig()];
        let overlaps = OverlapCollection::new();

        let result = calc_lq_coef(&contigs, &overlaps);
        let expected = 0.00;

        assert!((result - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn multiple_contigs_quarter_overlap() {
        let contigs = vec![
            make_contig(),
            make_contig(),
            make_contig(),
            make_contig(),
        ];
        let mut overlaps = OverlapCollection::new();
        overlaps.add(0, end_match(0, 1));
        overlaps.add(1, start_match(1, 0));

        let result = calc_lq_coef(&contigs, &overlaps);
        let expected = 25.00;

        assert!((result - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn both_contigs_fully_matched() {
        let contigs = vec![make_contig(), make_contig()];
        let mut overlaps = OverlapCollection::new();
        overlaps.add(0, start_match(0, 1));
        overlaps.add(0, end_match(0, 1));
        overlaps.add(1, start_match(1, 0));
        overlaps.add(1, end_match(1, 0));

        let result = calc_lq_coef(&contigs, &overlaps);
        let expected = 100.00;

        assert!((result - expected).abs() < f64::EPSILON);
    }

    #[test]
    fn empty_collection() {
        let contigs: Vec<ContigRecord> = vec![];
        let overlaps = OverlapCollection::new();

        let result = calc_lq_coef(&contigs, &overlaps);

        assert!(result.is_nan());
    }
}
