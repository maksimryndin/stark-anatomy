use crate::field;
use crate::merkle;
use crate::fiat_shamir;

pub struct Fri {
    pub offset: field::Felt,
    pub omega: field::Felt,
    pub domain_length: usize,
    pub field: field::Field,
    pub expansion_factor: usize,
    pub num_colinearity_tests: usize,
}

impl Fri {
    fn num_rounds(&self) -> usize {
        let mut codeword_length = self.domain_length;
        let mut num_rounds = 0;
        while codeword_length > self.expansion_factor && 4 * self.num_colinearity_tests < codeword_length {
            codeword_length /= 2;
            num_rounds += 1;
        }
        num_rounds
    }

    fn eval_domain(&self) -> Vec<field::Felt> {
        (0..self.domain_length).map(|i| self.offset * (self.omega.pow(i as u128))).collect()
    }

    fn prove(&self, codeword: Vec<field::Felt>, proof_stream: fiat_shamir::ProofStream) -> Vec<usize> {
        debug_assert_eq!(self.domain_length, codeword.len(), "initial codeword length does not match length of initial codeword");

        // commit phase
        let codewords = self.commit(codeword, &mut proof_stream);

        // get indices
        let top_level_indices = self.sample_indices(proof_stream.prover_fiat_shamir(), codewords[1].len(), codewords.last().unwrap().len(), self.num_colinearity_tests);
        let mut indices = top_level_indices.clone();

        // query phase
        for i in 0..codewords.len()-1 {
            indices = indices.into_iter().map(|index| index % (codewords[i].len() / 2)).collect();
            self.query(&codewords[i], &codewords[i+1], &indices, &mut proof_stream);
        }
        top_level_indices
    }

    fn commit(&self, codeword: Vec<field::Felt>, proof_stream: &mut fiat_shamir::ProofStream) -> Vec<Vec<field::Felt>> {
        let one = self.field.one();
        let two = field::Felt::new(2, self.field);
        let mut omega = self.omega;
        let mut offset = self.offset;
        let mut codewords = vec![];
        let num_rounds = self.num_rounds();

        // for each round
        for r in 0..num_rounds {

            // compute and send Merkle root
            let root = merkle::commit(&codeword);
            proof_stream.push(root);

            // prepare next round if necessary
            if r == num_rounds - 1 {
                break;
            }

            // get challenge
            let alpha = self.field.sample(&proof_stream.prover_fiat_shamir());

            // split and fold
            let next_codeword: Vec<field::Felt> = (0..codeword.len()/2).into_iter().map(|i| 
                two.inverse() * ( (one + alpha / (offset * (omega.pow(i as u128))) ) * codeword[i] + (one - alpha / (offset * (omega.pow(i as u128))) ) * codeword[codeword.len()/2 + i] )
            ).collect();

            // collect codeword
            codewords.push(codeword);

            omega = omega.pow(2);
            offset = offset.pow(2);
            
        }

        // send last codeword
        proof_stream.push(codeword.clone());

        // collect last codeword too
        codewords.push(codeword);
        codewords
    }

    fn query(&self, current_codeword: & Vec<field::Felt>, next_codeword: &Vec<field::Felt>, c_indices: Vec<usize>, proof_stream: &mut fiat_shamir::ProofStream) -> Vec<usize> {
        // infer a and b indices
        let a_indices = c_indices.clone();
        let b_indices: Vec<usize> = c_indices.iter().map(|index| index + current_codeword.len()/2).collect();

        // reveal leafs
        for s in 0..self.num_colinearity_tests {
            proof_stream.push((current_codeword[a_indices[s]], current_codeword[b_indices[s]], next_codeword[c_indices[s]]));
        }

        // reveal authentication path
        for s in 0..self.num_colinearity_tests {
            proof_stream.push(merkle::open(a_indices[s], current_codeword));
            proof_stream.push(merkle::open(b_indices[s], current_codeword));
            proof_stream.push(merkle::open(c_indices[s], next_codeword));
        }

        a_indices.into_iter().chain(b_indices.into_iter()).collect()
    }

    fn sample_index(byte_array: merkle::MerkleDigest, size: usize) -> usize {
        let mut acc = 0u64;
        for b in byte_array.0.into_iter() {
            let b: u64 = b.into();
            acc = (acc << 8) ^ b;
        }
        let acc: usize = acc.try_into().unwrap();
        acc % size
    }

    fn sample_indices(&self, seed: [u8; 16], size: usize, reduced_size: usize, number: usize) -> Vec<usize> {
        debug_assert!(number <= 2 * reduced_size, "not enough entropy in indices wrt last codeword");
        debug_assert!(number <= reduced_size, "cannot sample more indices than available in last codeword; requested: {number}, available: {reduced_size}");

        let mut indices = vec![];
        let mut reduced_indices = vec![];
        let mut counter = 0u128;
        while indices.len() < number {
            let bytes: Vec<u8> = seed.into_iter().chain(counter.to_be_bytes().into_iter()).collect();
            let index = Fri::sample_index(merkle::hash(&bytes), size);
            let reduced_index = index % reduced_size;
            counter += 1;
            if reduced_indices.binary_search(&reduced_index).is_ok() {
                indices.push(index);
                reduced_indices.push(reduced_index);
            }
        }
        indices
    }
}