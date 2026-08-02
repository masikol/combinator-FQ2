
use std::collections::HashSet;


/// IUPAC nucleotide-code validation.
///
/// Provides the IUPAC nucleic acid alphabet.
pub mod nucl_bases {
    /// Adenine (`A`).
    pub const ADENINE   : char = 'A';
    /// Thymine (`T`).
    pub const THYMINE   : char = 'T';
    /// Guanine (`G`).
    pub const GUANINE   : char = 'G';
    /// Cytosine (`C`).
    pub const CYTOSINE  : char = 'C';
    /// Purine (`R`, A or G).
    pub const PURINE    : char = 'R';
    /// Pyrimidine (`Y`, C or T).
    pub const PYRIMIDINE: char = 'Y';
    /// Weak (`W`, A or T).
    pub const WEAK      : char = 'W';
    /// Strong (`S`, C or G).
    pub const STRONG    : char = 'S';
    /// Amino (`M`, A or C).
    pub const AMINO     : char = 'M';
    /// Keto (`K`, G or T).
    pub const KETO      : char = 'K';
    /// Not A (`B`, C/G/T).
    pub const NOT_A     : char = 'B';
    /// Not C (`D`, A/G/T).
    pub const NOT_C     : char = 'D';
    /// Not G (`H`, A/C/T).
    pub const NOT_G     : char = 'H';
    /// Not T (`V`, A/C/G).
    pub const NOT_T     : char = 'V';
    /// Any base (`N`).
    pub const ANY       : char = 'N';
    /// Uracil (`U`).
    pub const URACIL    : char = 'U';
}


/// Validates sequences against the IUPAC nucleotide alphabet.
///
/// The accepted alphabet is the full IUPAC set in both upper and lower case.
pub struct IUPACValidator {
    valid_base_set: HashSet<char>,
}

impl IUPACValidator {
    /// Creates a validator accepting every IUPAC nucleotide code
    /// in both cases.
    pub fn new() -> IUPACValidator {
        let mut valid_base_set = HashSet::new();
        let valid_bases = vec![
            nucl_bases::ADENINE,
            nucl_bases::THYMINE,
            nucl_bases::GUANINE,
            nucl_bases::CYTOSINE,
            nucl_bases::PURINE,
            nucl_bases::PYRIMIDINE,
            nucl_bases::WEAK,
            nucl_bases::STRONG,
            nucl_bases::AMINO,
            nucl_bases::KETO,
            nucl_bases::NOT_A,
            nucl_bases::NOT_C,
            nucl_bases::NOT_G,
            nucl_bases::NOT_T,
            nucl_bases::ANY,
            nucl_bases::URACIL,
        ];
        for base in valid_bases {
            // to_ascii_uppercase call adds some overhead here
            //   but we tolerate it for the sake of explicity
            valid_base_set.insert(base.to_ascii_uppercase());
            valid_base_set.insert(base.to_ascii_lowercase());
        }
        IUPACValidator {
            valid_base_set: valid_base_set
        }
    }

    /// Checks that every character in `seq` is a valid IUPAC nucleotide
    /// code (either case).
    ///
    /// # Errors
    ///
    /// Returns `Err` naming the first character that is not in the
    /// alphabet.
    pub fn validate(&self, seq: &String) -> Result<(), String> {
        for c in seq.chars() {
            if !self.valid_base_set.contains(&c) {
                return Err(
                    format!("Error: non-IUPAC character encountered: {}", c)
                );
            }
        }
        Ok(())
    }
}


// >>> Tests >>>

#[cfg(test)]
mod tests_iupac_validator {
    use super::*;

    #[test]
    fn valid_dna_bases_uppercase() {
        let validator = IUPACValidator::new();
        assert_eq!(validator.validate(&"ACGT".into()), Ok(()));
    }

    #[test]
    fn valid_dna_bases_lowercase() {
        let validator = IUPACValidator::new();
        assert_eq!(validator.validate(&"acgt".into()), Ok(()));
    }

    #[test]
    fn valid_all_bases_uppercase() {
        let validator = IUPACValidator::new();
        assert_eq!(
            validator.validate(&"ATGCRYSWKMBDHVNU".into()),
            Ok(())
        );
    }

    #[test]
    fn valid_mixed_case() {
        let validator = IUPACValidator::new();
        assert_eq!(validator.validate(&"AtGcRyU".into()), Ok(()));
    }

    #[test]
    fn empty_string() {
        let validator = IUPACValidator::new();
        assert_eq!(validator.validate(&String::new()), Ok(()));
    }

    #[test]
    fn invalid_character_fails() {
        let validator = IUPACValidator::new();
        let result = validator.validate(&"ACGF".into());
        assert_eq!(
            result,
            Err("Error: non-IUPAC character encountered: F".to_string())
        );
    }

    #[test]
    fn non_letter_character_fails() {
        let validator = IUPACValidator::new();
        let result = validator.validate(&"AC GT".into());
        assert_eq!(
            result,
            Err("Error: non-IUPAC character encountered:  ".to_string())
        );
    }

    #[test]
    fn digit_character_fails() {
        let validator = IUPACValidator::new();
        let result = validator.validate(&"AC3T".into());
        assert_eq!(
            result,
            Err("Error: non-IUPAC character encountered: 3".to_string())
        );
    }

    #[test]
    fn nonascii_character_fails() {
        let validator = IUPACValidator::new();
        let result = validator.validate(&"ACўT".into());
        assert_eq!(
            result,
            Err("Error: non-IUPAC character encountered: ў".to_string())
        );
    }
}
