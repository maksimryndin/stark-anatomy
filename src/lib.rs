pub mod field;
pub mod merkle;
pub mod polynomial;
pub mod fri;
pub mod fiat_shamir;

pub fn main() -> field::Field {
    field::Field::new(field::PRIME)
}
