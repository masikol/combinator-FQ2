
use std::sync::OnceLock;

use regex::Regex;


static SPADES_NAME_REGEX: OnceLock<Regex> = OnceLock::new();

pub fn get_spades_name_regex() -> &'static Regex {
    SPADES_NAME_REGEX.get_or_init(
        || Regex::new(
            r"^(NODE_\d+)_length_\d+_cov_(\d+(\.\d+)?)"
        ).unwrap()
    )
}
