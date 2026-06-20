
mod find_overlap;

use find_overlap::{find_overlap_s2s, find_overlap_e2s, find_overlap_e2e};


fn main() {
    let overlap: usize = find_overlap_s2s(
        &String::from("AGTCaaaaaaaaa"),
        &String::from("AGTCttttttttttttttt"),
        3,
        5
    );
    println!("overlap = {overlap}");

    let overlap: usize = find_overlap_e2s(
        &String::from("aaaaaaaaaaAGTC"),
        &String::from("AGTCttttttttttttt"),
        3,
        5
    );
    println!("overlap = {overlap}");

    let overlap: usize = find_overlap_e2e(
        &String::from("aaaaaaaAGTC"),
        &String::from("tttttttttttttAGTC"),
        3,
        5
    );
    println!("overlap = {overlap}");
}
