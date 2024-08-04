use crate::field;
use std::collections::HashMap;

pub struct MPolynomial {
    dictionary: HashMap<Vec<u128>, field::Felt>,
}

impl MPolynomial {
    pub fn new(dictionary: HashMap<Vec<u128>, field::Felt>) -> Self {
        Self { dictionary }
    }

    pub fn zero() -> Self {
        Self {
            dictionary: HashMap::new(),
        }
    }

    pub fn constant(element: field::Felt) -> Self {
        Self {
            dictionary: HashMap::from([(vec![0], element)]),
        }
    }

    pub fn is_zero(&self) -> bool {
        for v in self.dictionary.values() {
            if !v.is_zero() {
                return false;
            }
        }
        true
    }

    pub fn variables(num_variables: u128, fld: field::Field) -> Vec<Self> {
        todo!()
    }
}
