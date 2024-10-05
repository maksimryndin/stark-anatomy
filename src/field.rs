use crypto_bigint::{
    modular::runtime_mod::{DynResidue, DynResidueParams},
    CtChoice, NonZero, Zero, U128, U256,
};
use std::fmt;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub};

pub const PRIME: u128 = 1 + 407 * (1 << 119);

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Felt {
    value: U128,
    field: Field,
}

impl fmt::Debug for Felt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_primitive())
    }
}

impl Felt {
    pub fn new(value: u128, field: Field) -> Self {
        let value = U128::from_u128(value);
        Self::from_U128(value, field)
    }

    #[allow(non_snake_case)]
    pub fn from_U128(value: U128, field: Field) -> Self {
        debug_assert!(value < field.p, "felt value cannot exceed field modulo");
        Self { value, field }
    }

    pub fn inverse(self) -> Felt {
        self.field.inverse(self)
    }

    pub fn is_zero(self) -> bool {
        Into::<bool>::into(self.value.is_zero())
    }

    pub fn pow(&self, exponent: u128) -> Felt {
        let base_mod = DynResidue::new(&self.value, DynResidueParams::new(&self.field.p));
        let res = base_mod.pow(&U128::from_u128(exponent));
        let value: U128 = res.retrieve();
        Felt::from_U128(value, self.field)
    }

    pub fn as_primitive(&self) -> u128 {
        let bytes = self.to_be_bytes();
        u128::from_be_bytes(bytes)
    }

    pub fn to_be_bytes(&self) -> [u8; 16] {
        let mut bytes = [0u8; 16];
        for (i, b) in self
            .value
            .as_words()
            .into_iter()
            .rev()
            .flat_map(|word| word.to_be_bytes())
            .enumerate()
        {
            bytes[i] = b;
        }
        bytes
    }

    pub fn field(self) -> Field {
        self.field
    }
}

impl Add for Felt {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        self.field.add(self, other)
    }
}

impl AddAssign for Felt {
    fn add_assign(&mut self, other: Self) {
        self.value = self.field.add(*self, other).value;
    }
}

impl Sub for Felt {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        self.field.subtract(self, other)
    }
}

impl Mul for Felt {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        self.field.multiply(self, rhs)
    }
}

impl MulAssign for Felt {
    fn mul_assign(&mut self, other: Self) {
        self.value = self.field.multiply(*self, other).value;
    }
}

impl Div for Felt {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        self.field.divide(self, rhs)
    }
}

impl Neg for Felt {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.field.negate(self)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Field {
    p: U128,
}

impl Field {
    pub fn new(p: u128) -> Self {
        debug_assert!(p != 0, "field modulo cannot be zero");
        Self {
            p: U128::from_u128(p),
        }
    }

    pub fn zero(&self) -> Felt {
        Felt::from_U128(U128::ZERO, *self)
    }

    pub fn one(&self) -> Felt {
        Felt::from_U128(U128::ONE, *self)
    }

    pub fn multiply(&self, left: Felt, right: Felt) -> Felt {
        debug_assert!(
            left.field == right.field,
            "two felts should belong to the same field"
        );
        let product: U256 = left.value.mul(&right.value);
        let modulo = NonZero::new((&self.p).into()).expect("field modulo is not zero");
        let product_words = product.rem(&modulo).to_words();
        debug_assert!(
            product_words[2] == 0,
            "mod product of two felts should fit a result felt"
        );
        let mod_prod: U128 = U128::from_words((&product_words[..2]).try_into().unwrap());
        Felt::from_U128(mod_prod, *self)
    }

    pub fn divide(&self, left: Felt, right: Felt) -> Felt {
        debug_assert!(
            left.field == right.field,
            "two felts should belong to the same field"
        );
        debug_assert!(!Into::<bool>::into(right.value.is_zero()), "divide by zero");
        let inverse = self.inverse(right);
        self.multiply(left, inverse)
    }

    pub fn add(&self, left: Felt, right: Felt) -> Felt {
        debug_assert!(
            left.field == right.field,
            "two felts should belong to the same field"
        );
        Felt::from_U128(left.value.add_mod(&right.value, &self.p), *self)
    }

    pub fn subtract(&self, left: Felt, right: Felt) -> Felt {
        debug_assert!(
            left.field == right.field,
            "two felts should belong to the same field"
        );
        Felt::from_U128(left.value.sub_mod(&right.value, &self.p), *self)
    }

    pub fn negate(&self, operand: Felt) -> Felt {
        Felt::from_U128(operand.value.neg_mod(&self.p), *self)
    }

    pub fn inverse(&self, operand: Felt) -> Felt {
        let (res, success): (U128, CtChoice) = operand.value.inv_mod(&self.p);
        debug_assert!(Into::<bool>::into(success));
        Felt::from_U128(res, *self)
    }

    pub fn generator(&self) -> Felt {
        debug_assert!(
            self.p == U128::from_u128(PRIME),
            "Do not know generator for other fields beyond 1+407*2^119"
        );
        Felt::new(85408008396924667383611388730472331217, *self)
    }

    pub fn primitive_nth_root(&self, n: u128) -> Felt {
        debug_assert!(
            self.p == U128::from_u128(PRIME),
            "Unknown field, can't return root of unity."
        );
        debug_assert!(
            n <= 1 << 119 && (n & (n - 1)) == 0,
            "Field does not have nth root of unity where n > 2^119 or not power of two."
        );
        let mut root = Felt::new(85408008396924667383611388730472331217, *self);
        let mut order: u128 = 1 << 119;
        while order != n {
            root = root.pow(2);
            order = order / 2;
        }
        root
    }

    pub fn sample(&self, byte_array: &[u8]) -> Felt {
        let mut acc = 0u128;
        for b in byte_array.into_iter() {
            let b: u128 = (*b).into();
            acc = (acc << 8) ^ b;
        }
        let p = Felt::from_U128(self.p, *self);
        Felt::new(acc % p.as_primitive(), *self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mod_add() {
        let field = Field::new(PRIME);
        let felt1 = Felt::new(PRIME - 1, field);
        let felt2 = Felt::new(PRIME - 2, field);
        assert_eq!(field.add(felt1, felt2), Felt::new(PRIME - 3, field));
    }

    #[test]
    fn mod_sub() {
        let field = Field::new(PRIME);
        let felt1 = Felt::new(PRIME - 2, field);
        let felt2 = Felt::new(PRIME - 1, field);
        assert_eq!(field.subtract(felt1, felt2), Felt::new(PRIME - 1, field));
    }

    #[test]
    fn mod_neg() {
        let field = Field::new(PRIME);
        let felt = Felt::new(PRIME - 2, field);
        assert_eq!(field.negate(felt), Felt::new(2, field));
        let felt = Felt::new(2, field);
        assert_eq!(field.negate(felt), Felt::new(PRIME - 2, field));
    }

    #[test]
    fn mod_inv() {
        let field = Field::new(PRIME);
        let felt = Felt::new(PRIME - 2, field);
        assert_eq!(
            field.inverse(felt),
            Felt::new(135248948571115190067962368383525060608, field)
        );
    }

    #[test]
    fn mod_mul() {
        let field = Field::new(PRIME);
        let felt1 = Felt::new(PRIME - 2, field);
        let felt2 = Felt::new(PRIME - 1, field);
        assert_eq!(field.multiply(felt1, felt2), Felt::new(2, field));
    }

    #[test]
    fn mod_div() {
        let field = Field::new(PRIME);
        let felt1 = Felt::new(PRIME - 2, field);
        let felt2 = Felt::new(PRIME - 1, field);
        assert_eq!(field.divide(felt1, felt2), Felt::new(2, field));
    }

    #[test]
    fn mod_pow() {
        let field = Field::new(PRIME);
        let felt1 = Felt::new(PRIME - 2, field);
        assert_eq!(felt1.pow(PRIME - 1), field.one());
    }

    #[test]
    fn bytes_representation() {
        let field = Field::new(PRIME);
        let felt = Felt::new(129, field);
        assert_eq!(
            felt.to_be_bytes(),
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 129u8]
        );
        let felt = Felt::new(PRIME - 1, field);
        assert_eq!(
            felt.to_be_bytes(),
            [203u8, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        );
    }
}
