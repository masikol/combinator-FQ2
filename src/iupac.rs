
use std::collections::HashSet;


pub mod nucl_bases {
    pub const ADENINE   : char = 'A';
    pub const THYMINE   : char = 'T';
    pub const GUANINE   : char = 'G';
    pub const CYTOSINE  : char = 'C';
    pub const PURINE    : char = 'R';
    pub const PYRIMIDINE: char = 'Y';
    pub const WEAK      : char = 'W';
    pub const STRONG    : char = 'S';
    pub const AMINO     : char = 'M';
    pub const KETO      : char = 'K';
    pub const NOT_A     : char = 'B';
    pub const NOT_C     : char = 'D';
    pub const NOT_G     : char = 'H';
    pub const NOT_T     : char = 'V';
    pub const ANY       : char = 'N';
    pub const URACIL    : char = 'U';
}


pub struct IUPACValidator {
    valid_base_set: HashSet<char>,
}

impl IUPACValidator {
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
