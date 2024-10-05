use crate::field;
use crate::polynomial::univariate;
use std::collections::HashMap;
use std::ops::{Add, Mul, Neg, Sub};

#[derive(Debug, PartialEq, Eq, Clone)]
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

    pub fn pow(&self, mut exponent: u128) -> Self {
        if self.is_zero() {
            return Self::new(HashMap::new());
        }
        let (element, value) = self.dictionary.iter().next().unwrap();
        let fld = value.field();
        let num_variables = element.len();
        let mut acc = Self::new(HashMap::from([(vec![0; num_variables], fld.one())]));
        if exponent == 0 {
            return acc;
        }
        let mut base = self.clone();

        while exponent > 1 {
            if (exponent & 1) == 1 {
                acc = &acc * &base;
            }
            exponent /= 2;
            base = &base * &base;
        }
        &acc * &base
    }

    pub fn constant(element: field::Felt) -> Self {
        Self::new(HashMap::from([(vec![0], element)]))
    }

    pub fn is_zero(&self) -> bool {
        for v in self.dictionary.values() {
            if !v.is_zero() {
                return false;
            }
        }
        true
    }

    pub fn variables(num_variables: usize, fld: field::Field) -> Vec<Self> {
        (0..num_variables)
            .into_iter()
            .map(|i| {
                let mut exponent = vec![0; num_variables];
                exponent[i] = 1;
                Self::new(HashMap::from([(exponent, fld.one())]))
            })
            .collect()
    }

    pub fn lift(polynomial: &univariate::Polynomial, variable_index: usize) -> Self {
        if polynomial.is_zero() {
            return Self::new(HashMap::new());
        }
        let fld = polynomial.field().unwrap();
        let variables = Self::variables(variable_index + 1, fld);
        let x = variables.last().unwrap();
        let degree: usize = polynomial.degree().unwrap().try_into().unwrap();
        let mut acc = Self::new(HashMap::with_capacity(degree));
        for i in 0..degree {
            acc = &acc + &(&Self::constant(polynomial.coefficients()[i]) * &x.pow(i as u128))
        }
        acc
    }

    pub fn evaluate(&self, point: &[field::Felt]) -> field::Felt {
        debug_assert!(!point.is_empty(), "cannot evaluate an empty point");
        let mut acc = point[0].field().zero();
        for (k, v) in &self.dictionary {
            let mut prod = *v;
            for i in 0..k.len() {
                prod *= point[i].pow(k[i]);
            }
            acc += prod;
        }
        acc
    }

    pub fn evaluate_symbolic(&self, point: &[univariate::Polynomial]) -> univariate::Polynomial {
        let mut acc = univariate::Polynomial::new([]);
        for (k, v) in &self.dictionary {
            let mut prod = univariate::Polynomial::new([*v]);
            for i in 0..k.len() {
                prod = &prod * &point[i].pow(k[i]);
            }
            acc += prod;
        }
        acc
    }
}

impl<'a, 'b> Add<&'b MPolynomial> for &'a MPolynomial {
    type Output = MPolynomial;

    fn add(self, other: &'b MPolynomial) -> MPolynomial {
        let num_variables = self
            .dictionary
            .keys()
            .chain(other.dictionary.keys())
            .map(|exponent| exponent.len())
            .max()
            .unwrap_or(0);
        let mut dictionary = HashMap::with_capacity(self.dictionary.len() + other.dictionary.len());
        for (exponent, coefficient) in self.dictionary.iter().chain(other.dictionary.iter()) {
            let pad = num_variables - exponent.len();
            let new_exponent: Vec<u128> = exponent
                .into_iter()
                .cloned()
                .chain(0..pad as u128)
                .collect();
            let zero = coefficient.field().zero();
            *dictionary.entry(new_exponent).or_insert(zero) += *coefficient;
        }
        MPolynomial::new(dictionary)
    }
}

impl Sub<MPolynomial> for &MPolynomial {
    type Output = MPolynomial;

    fn sub(self, other: MPolynomial) -> MPolynomial {
        self.add(&(-other))
    }
}

impl<'a, 'b> Mul<&'b MPolynomial> for &'a MPolynomial {
    type Output = MPolynomial;

    fn mul(self, other: &'b MPolynomial) -> MPolynomial {
        let num_variables = self
            .dictionary
            .keys()
            .chain(other.dictionary.keys())
            .map(|exponent| exponent.len())
            .max()
            .unwrap_or(0);
        let mut dictionary = HashMap::with_capacity(self.dictionary.len() + other.dictionary.len());
        for (exponent0, coefficient0) in &self.dictionary {
            for (exponent1, coefficient1) in &other.dictionary {
                let mut exponent = vec![0; num_variables];
                (&mut exponent[0..exponent0.len()]).copy_from_slice(&exponent0);
                exponent1
                    .into_iter()
                    .enumerate()
                    .for_each(|(i, k)| exponent[i] += k);
                let zero = coefficient0.field().zero();
                *dictionary.entry(exponent).or_insert(zero) += *coefficient0 * *coefficient1;
            }
        }
        MPolynomial::new(dictionary)
    }
}

impl Neg for MPolynomial {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.dictionary.values_mut().for_each(|c| *c = -*c);
        self
    }
}

// TODO tests https://github.com/aszepieniec/stark-anatomy/blob/master/code/test_multivariate.py
