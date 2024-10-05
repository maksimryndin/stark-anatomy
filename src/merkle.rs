use blake2::{Blake2b512, Digest};
use std::ops::Add;
use crate::field;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MerkleDigest(pub [u8; 64]);

impl From<[u8; 64]> for MerkleDigest {
    fn from(arr: [u8; 64]) -> Self {
        MerkleDigest(arr)
    }
}

impl Add for MerkleDigest {
    type Output = [u8; 128];

    fn add(self, other: Self) -> Self::Output {
        let mut sum = [0u8; 128];
        for i in 0..self.0.len() {
            sum[i] = self.0[i];
            sum[i + self.0.len()] = other.0[i];
        }
        sum
    }
}

pub fn hash(bytes: &[u8]) -> MerkleDigest {
    let mut hasher = Blake2b512::new();
    hasher.update(bytes);
    let arr: [u8; 64] = hasher.finalize().into();
    arr.into()
}

struct Merkle {}

impl Merkle {
    // computes the Merkle root of a given array.
    fn commit(&self, leafs: &[MerkleDigest]) -> MerkleDigest {
        debug_assert!(
            leafs.len() & (leafs.len() - 1) == 0,
            "length must be power of two"
        );
        if leafs.len() == 1 {
            leafs[0]
        } else {
            let hash0 = self.commit(&leafs[..leafs.len() / 2]);
            let hash1 = self.commit(&leafs[leafs.len() / 2..]);

            hash(&(hash0 + hash1))
        }
    }

    // computes the authentication path of an indicated leaf in the Merkle tree.
    fn open(&self, index: usize, leafs: &[MerkleDigest], path: &mut Vec<MerkleDigest>) {
        debug_assert!(
            leafs.len() & (leafs.len() - 1) == 0,
            "length must be power of two"
        );
        debug_assert!(
            index < leafs.len(),
            "cannot open invalid index"
        );
        debug_assert!(
            path.len() >= leafs.len().ilog2() as usize,
            "insufficient path capacity"
        );

        if leafs.len() == 2 {
            path.push(leafs[1 - index]);
        } else if index < leafs.len() / 2 {
            self.open(index, &leafs[..leafs.len() / 2], path);
            path.push(self.commit(&leafs[leafs.len() / 2..]));
        } else {
            self.open(leafs.len() - index, &leafs[leafs.len() / 2..], path);
            path.push(self.commit(&leafs[..leafs.len() / 2]));
        }
    }

    // verifies that a given leaf is an element of the committed vector at the given index.
    fn verify(
        &self,
        root: MerkleDigest,
        index: usize,
        path: &[MerkleDigest],
        leaf: MerkleDigest,
    ) -> bool {
        debug_assert!(
            index < (1 << path.len()),
            "cannot verify invalid index"
        );
        if path.len() == 1 {
            if index == 0 {
                root == hash(&(leaf + path[0]))
            } else {
                root == hash(&(path[0] + leaf))
            }
        } else {
            if index & 1 == 0 {
                self.verify(root, index >> 1, &path[1..], hash(&(leaf + path[0])))
            } else {
                self.verify(root, index >> 1, &path[1..], hash(&(path[0] + leaf)))
            }
        }
    }
}

pub trait ToBeBytes {
    fn to_be_bytes(&self) -> [u8; 32];
}

impl ToBeBytes for field::Felt {
    fn to_be_bytes(&self) -> [u8; 32] {
        let mut bytes = [0u8; 32];
        for (i, b) in self.to_be_bytes().into_iter().enumerate() {
            bytes[i+16] = b;
        }
        bytes
    }
}

pub fn commit<T: ToBeBytes>(data_array: &[T]) -> MerkleDigest {
    let digests: Vec<MerkleDigest> = data_array
        .into_iter()
        .map(|d| hash(&d.to_be_bytes()))
        .collect();
    Merkle {}.commit(&digests)
}

pub fn open<T: ToBeBytes>(index: usize, data_array: &[T]) -> Vec<MerkleDigest> {
    let mut path = Vec::with_capacity(data_array.len().ilog2() as usize);
    let digests: Vec<MerkleDigest> = data_array
        .into_iter()
        .map(|d| hash(&d.to_be_bytes()))
        .collect();
    Merkle {}.open(index, &digests, &mut path);
    path
}

pub fn verify<T: ToBeBytes>(
    root: MerkleDigest,
    index: usize,
    path: &[MerkleDigest],
    data_element: T,
) -> bool {
    Merkle {}.verify(root, index, path, hash(&data_element.to_be_bytes()))
}

// TODO tests