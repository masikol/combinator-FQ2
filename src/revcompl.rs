
use crate::iupac::nucl_bases::*;


pub fn revcompl(seq: &String) -> Result<String, char> {
    seq.chars()
        .map(make_compl_base)
        .rev()
        .collect()
}


fn make_compl_base(base: char) -> Result<char, char> {
    match base {
        ADENINE    => Ok(THYMINE),
        THYMINE    => Ok(ADENINE),
        GUANINE    => Ok(CYTOSINE),
        CYTOSINE   => Ok(GUANINE),
        PURINE     => Ok(PYRIMIDINE),
        PYRIMIDINE => Ok(PURINE),
        WEAK       => Ok(WEAK),
        STRONG     => Ok(STRONG),
        AMINO      => Ok(KETO),
        KETO       => Ok(AMINO),
        NOT_A      => Ok(NOT_T),
        NOT_T      => Ok(NOT_A),
        NOT_G      => Ok(NOT_C),
        NOT_C      => Ok(NOT_G),
        ANY        => Ok(ANY),
        URACIL     => Ok(ADENINE),
        _          => Err(base)
    }
}


#[cfg(test)]
mod test_revcompl {

    use super::revcompl;


    #[test]
    fn ok_revcompl_seq() {
        let seq      = String::from("WRTYUASDGKCBM");
        let expected = String::from("KVGMCHSTARAYW");

        let revcompl_result = revcompl(&seq);
        assert_eq!(
            revcompl_result,
            Ok(expected),
            "Reverse-complement sequence is formed incorrectly"
        );
    }


    #[test]
    fn invalid_revcompl_seq() {
        let invalid_char = 'F';
        let valid_seq    = String::from("WRTYUASDGKCBM");
        let invalid_seq  = format!("{}{}", valid_seq, invalid_char);

        let revcompl_result = revcompl(&invalid_seq);
        assert_eq!(
            revcompl_result,
            Err(invalid_char),
            "Input sequence is invalid, but `revcompl` function didn't fail"
        );
    }
}


#[cfg(test)]
mod test_make_compl_base {

    use super::make_compl_base;
    use super::{ADENINE, THYMINE};


    #[test]
    fn ok_compl_base() {
        let compl_base_result = make_compl_base(ADENINE);
        assert_eq!(
            compl_base_result,
            Ok(THYMINE),
            "Base `{}` is not complement to `{}`",
                compl_base_result.unwrap(),
                ADENINE
        );
    }


    #[test]
    fn invalid_compl_base() {
        let invalid_char: char = 'F';
        let compl_base_result = make_compl_base(invalid_char);
        assert_eq!(
            compl_base_result,
            Err(invalid_char),
            "No error on invalid char `{}`", invalid_char
        );
    }
}
