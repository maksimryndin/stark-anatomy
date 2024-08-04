use crate::field;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Rem, Sub, SubAssign};

#[derive(Debug, Clone)]
pub struct Polynomial {
    coefficients: Vec<field::Felt>,
}

impl Polynomial {
    pub fn new(coefficients: impl Into<Vec<field::Felt>>) -> Self {
        Self {
            coefficients: coefficients.into(),
        }
    }

    // https://en.wikipedia.org/wiki/Lagrange_polynomial
    pub fn interpolate_domain(domain: &[field::Felt], values: &[field::Felt]) -> Self {
        debug_assert_eq!(
            domain.len(),
            values.len(),
            "number of elements in domain does not match number of values -- cannot interpolate"
        );
        debug_assert!(!domain.is_empty(), "cannot interpolate between zero points");
        let fld = domain[0].field();
        let x_poly = Polynomial::new([fld.zero(), fld.one()]);
        let mut acc = Polynomial::new([]);
        for (x, y) in domain.into_iter().zip(values.into_iter()) {
            let mut prod = Polynomial::new([*y]);
            for d in domain.into_iter().filter(|d| *d != x) {
                prod = &prod
                    * &(&(x_poly.clone() - Polynomial::new([*d]))
                        * &Polynomial::new([(*x - *d).inverse()]));
            }
            acc += prod;
        }
        acc
    }

    pub fn zeroifier_domain(domain: &[field::Felt]) -> Self {
        debug_assert!(!domain.is_empty(), "cannot zeroifier zero points");
        let fld = domain[0].field();
        let x_poly = Polynomial::new([fld.zero(), fld.one()]);
        let mut acc = Polynomial::new([fld.one()]);
        for d in domain.into_iter() {
            acc = &acc * &(x_poly.clone() - Polynomial::new([*d]));
        }
        acc
    }

    pub fn degree(&self) -> Option<u64> {
        self.coefficients
            .iter()
            .rposition(|felt| !felt.is_zero())
            .map(|pos| pos.try_into().expect("polynomial degree fits u64"))
    }

    pub fn is_zero(&self) -> bool {
        self.degree().is_none()
    }

    pub fn leading_coefficient(&self) -> Option<field::Felt> {
        self.degree()
            .map(|degree| self.coefficients[degree as usize])
    }

    pub fn pow(&self, mut exponent: u128) -> Self {
        if self.is_zero() {
            return Polynomial::new([]);
        }
        let fld = self.coefficients[0].field();
        let mut acc = Polynomial::new([fld.one()]);
        if exponent == 0 {
            return acc;
        }
        let mut base = self.clone();

        while exponent > 1 {
            if (exponent & 1) == 1 {
                // TODO an external buffer to reuse
                acc = &acc * &base;
            }
            exponent /= 2;
            // TODO an external buffer to reuse
            base = &base * &base;
        }
        &acc * &base
    }

    pub fn scale(&mut self, factor: field::Felt) {
        self.coefficients
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| *c = factor.pow(i as u128) * *c);
    }

    pub fn test_collinearity(points: &[(field::Felt, field::Felt)]) -> bool {
        let domain: Vec<field::Felt> = points.iter().map(|(x, _)| *x).collect();
        let values: Vec<field::Felt> = points.iter().map(|(_, y)| *y).collect();
        let poly = Polynomial::interpolate_domain(&domain, &values);
        poly.degree() <= Some(1)
    }

    fn divide(self, denominator: &Self) -> Option<(Self, Self)> {
        if denominator.is_zero() {
            return None;
        }
        let numerator_degree = self.degree();
        let denominator_degree = denominator.degree();
        if numerator_degree < denominator_degree {
            return Some((Polynomial::new([]), self));
        }
        // safe to unwrap as checks were hold earlier
        let numerator_degree = numerator_degree.unwrap();
        let denominator_degree = denominator_degree.unwrap();
        let quotient_length = numerator_degree - denominator_degree + 1;
        let fld = self.coefficients[0].field();
        let mut quotient_buf = vec![fld.zero(); quotient_length as usize];
        let mut remainder = self;
        for _ in 0..quotient_length {
            let remainder_degree = remainder.degree();
            if remainder_degree.is_none() {
                break;
            }
            let remainder_degree = remainder_degree.unwrap();
            if remainder_degree < denominator_degree {
                break;
            }
            // safe to unwrap as checks were hold earlier
            let coefficient = remainder.leading_coefficient().unwrap()
                / denominator.leading_coefficient().unwrap();
            let shift = remainder_degree - denominator_degree;
            // TODO reuse allocation subtractee_coefficients for every iteration
            let mut subtractee_coefficients = vec![fld.zero(); shift as usize + 1];
            subtractee_coefficients[shift as usize] = coefficient;
            // TODO provide an external buffer Multiplication
            let subtractee = &Polynomial::new(subtractee_coefficients) * &denominator;
            quotient_buf[shift as usize] = coefficient;
            remainder -= subtractee;
        }
        Some((Polynomial::new(quotient_buf), remainder))
    }

    fn evaluate(&self, point: field::Felt) -> field::Felt {
        let mut xi = point.field().one();
        let mut value = point.field().zero();
        for c in &self.coefficients {
            value += *c * xi;
            xi *= point;
        }
        value
    }

    fn evaluate_domain_mut(&self, domain: &mut [field::Felt]) {
        domain.iter_mut().for_each(|p| *p = self.evaluate(*p));
    }
}

impl PartialEq for Polynomial {
    fn eq(&self, other: &Self) -> bool {
        let (long, short) = if self.coefficients.len() > other.coefficients.len() {
            (self, other)
        } else {
            (other, self)
        };
        for i in 0..short.coefficients.len() {
            if long.coefficients[i] != short.coefficients[i] {
                return false;
            }
        }
        for i in short.coefficients.len()..long.coefficients.len() {
            if !long.coefficients[i].is_zero() {
                return false;
            }
        }
        true
    }
}
impl Eq for Polynomial {}

impl Add for Polynomial {
    type Output = Self;

    fn add(mut self, mut other: Self) -> Self {
        if self.degree().is_none() {
            other
        } else if other.degree().is_none() {
            self
        } else if self.coefficients.len() > other.coefficients.len() {
            self.coefficients
                .iter_mut()
                .zip(other.coefficients.into_iter())
                .for_each(|(c, o)| *c += o);
            self
        } else {
            other
                .coefficients
                .iter_mut()
                .zip(self.coefficients.into_iter())
                .for_each(|(c, o)| *c += o);
            other
        }
    }
}

impl AddAssign for Polynomial {
    fn add_assign(&mut self, other: Self) {
        let original = std::mem::replace(self, Polynomial::new([]));
        *self = original.add(other);
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.add(-other)
    }
}

impl SubAssign for Polynomial {
    fn sub_assign(&mut self, other: Self) {
        let original = std::mem::replace(self, Polynomial::new([]));
        *self = original.sub(other);
    }
}

impl<'a, 'b> Mul<&'b Polynomial> for &'a Polynomial {
    type Output = Polynomial;

    // TODO make a method which takes a buffer for result
    fn mul(self, other: &'b Polynomial) -> Polynomial {
        if self.coefficients.is_empty() || other.coefficients.is_empty() {
            Polynomial::new([])
        } else {
            let fld = self.coefficients[0].field();
            let mut buf = vec![fld.zero(); self.coefficients.len() + other.coefficients.len() - 1];
            for i in 0..self.coefficients.len() {
                if self.coefficients[i].is_zero() {
                    continue;
                }
                for j in 0..other.coefficients.len() {
                    buf[i + j] += self.coefficients[i] * other.coefficients[j];
                }
            }
            Polynomial::new(buf)
        }
    }
}

impl Div<&Polynomial> for Polynomial {
    type Output = Self;

    fn div(self, other: &Self) -> Self::Output {
        let (quotient, remainder) = self
            .divide(other)
            .expect("cannot divide by zero polynomial");
        debug_assert!(
            remainder.is_zero(),
            "cannot perform polynomial division because remainder is not zero"
        );
        quotient
    }
}

impl Rem<&Polynomial> for Polynomial {
    type Output = Self;

    fn rem(self, modulus: &Self) -> Self::Output {
        let (_, remainder) = self
            .divide(modulus)
            .expect("cannot divide by zero polynomial");
        remainder
    }
}

impl Neg for Polynomial {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        self.coefficients.iter_mut().for_each(|c| *c = -*c);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polynomial_degree() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        // f(x) = 0*x^4 + 2*x^3 + 1*x^2 + 0*x^1 + 0
        let pol = Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]);
        assert_eq!(pol.degree(), Some(3));
        let pol = Polynomial::new([fld.zero(), fld.zero(), fld.zero(), fld.zero()]);
        assert_eq!(pol.degree(), None);
    }

    #[test]
    fn polynomial_negation() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        let pol = Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]);
        assert_eq!(
            -pol,
            Polynomial::new([fld.zero(), fld.zero(), -felt1, -felt2, fld.zero()])
        );
    }

    #[test]
    fn polynomial_addition() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        let pol = Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]);
        assert_eq!(
            -pol,
            Polynomial::new([fld.zero(), fld.zero(), -felt1, -felt2, fld.zero()])
        );
    }

    #[test]
    fn polynomial_multiplication() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        // p1(x) = 0*x^4 + 2*x^3 + 1*x^2 + 0*x^1 + 0
        let pol = Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]);
        let zero = Polynomial::new([fld.zero(), fld.zero(), fld.zero(), fld.zero()]);
        assert_eq!(&pol * &zero, zero);
        // p2(x) = 5*x^2 + 2*x^1 + 3
        let pol2 = Polynomial::new([
            field::Felt::new(3, fld),
            felt2,
            field::Felt::new(5, fld),
            fld.zero(),
        ]);
        // p(x) = p1(x) * p2(x) = 10*x^5 + 9*x^4 + 8*x^3 + 3*x^2
        let prod = Polynomial::new([
            fld.zero(),
            fld.zero(),
            field::Felt::new(3, fld),
            field::Felt::new(8, fld),
            field::Felt::new(9, fld),
            field::Felt::new(10, fld),
            fld.zero(),
        ]);
        assert_eq!(&pol * &pol2, prod);
    }

    #[test]
    fn polynomial_division() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        let pol = Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]);
        let zero = Polynomial::new([fld.zero(), fld.zero(), fld.zero(), fld.zero()]);
        assert_eq!(zero.clone() / &pol, zero);
        // p(x) = p1(x) * p2(x) = 10*x^5 + 9*x^4 + 8*x^3 + 3*x^2
        let prod = Polynomial::new([
            fld.zero(),
            fld.zero(),
            field::Felt::new(3, fld),
            field::Felt::new(8, fld),
            field::Felt::new(9, fld),
            field::Felt::new(10, fld),
            fld.zero(),
        ]);
        assert_eq!(
            prod / &pol,
            Polynomial::new([
                field::Felt::new(3, fld),
                felt2,
                field::Felt::new(5, fld),
                fld.zero(),
            ])
        );
    }

    #[test]
    fn polynomial_exponentiation() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        // p1(x) = 0*x^4 + 2*x^3 + 1*x^2 + 0*x^1 + 0
        let pol = Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]);
        assert_eq!(&(&pol * &pol) * &pol, pol.pow(3));
    }

    #[test]
    fn polynomial_equality() {
        let fld = field::Field::new(field::PRIME);
        let felt1 = field::Felt::new(1, fld);
        let felt2 = field::Felt::new(2, fld);
        assert_eq!(
            Polynomial::new([fld.zero(), fld.zero(), felt1, felt2, fld.zero()]),
            Polynomial::new([fld.zero(), fld.zero(), felt1, felt2])
        );
        assert_eq!(
            Polynomial::new([fld.zero(), fld.zero(), fld.zero(), fld.zero()]),
            Polynomial::new([])
        );
    }

    #[test]
    fn interpolation() {
        let fld = field::Field::new(field::PRIME);
        let domain = [
            fld.one(),
            field::Felt::new(2, fld),
            field::Felt::new(3, fld),
        ];
        let values = [
            fld.one(),
            field::Felt::new(4, fld),
            field::Felt::new(9, fld),
        ];
        assert_eq!(
            Polynomial::new([fld.zero(), fld.zero(), fld.one()]),
            Polynomial::interpolate_domain(&domain, &values)
        );
    }
}
