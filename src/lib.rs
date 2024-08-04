pub mod field;
pub mod polynomial;

pub fn main() -> field::Field {
    field::Field::new(field::PRIME)
}
