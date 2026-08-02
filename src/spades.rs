
use std::sync::OnceLock;

use regex::Regex;


/// Lazily-initialized regex matching SPAdes contig headers.
///
/// Captures the node id (`NODE_1` in group 1)
/// and the coverage value (group 2),
/// e.g. `NODE_1_length_100_cov_50.0`.
static SPADES_NAME_REGEX: OnceLock<Regex> = OnceLock::new();

/// Returns the regex for parsing SPAdes-style contig headers.
///
/// The regex is compiled once and cached for the process lifetime.
pub fn get_spades_name_regex() -> &'static Regex {
    SPADES_NAME_REGEX.get_or_init(
        || Regex::new(
            r"^(NODE_\d+)_length_\d+_cov_(\d+(\.\d+)?)"
        ).unwrap()
    )
}
