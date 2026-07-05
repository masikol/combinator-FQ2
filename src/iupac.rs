
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

    pub fn validate(&self, string: &String) -> Result<(), String> {
        for c in string.chars() {
            if !self.valid_base_set.contains(&c) {
                return Err(
                    format!("Error: non-IUPAC character encountered: {}", c)
                );
            }
        }
        Ok(())
    }
}
