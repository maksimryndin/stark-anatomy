use crate::field;

pub struct ProofStream {
    objects: Vec<Vec<field::Felt>>,
}

impl ProofStream {
    pub fn prover_fiat_shamir(&self) -> [u8; 32] {
        todo!()
    }

    pub fn push(&mut self, obj: Vec<field::Felt>) {
        self.objects.push(obj);
    }
}

// TODO use json serde

// from hashlib import shake_256
// import pickle as pickle # serialization

// class ProofStream:
//     def __init__( self ):
//         self.objects = []
//         self.read_index = 0

//     def push( self, obj ):
//         self.objects += [obj]

//     def pull( self ):
//         assert(self.read_index < len(self.objects)), "ProofStream: cannot pull object; queue empty."
//         obj = self.objects[self.read_index]
//         self.read_index += 1
//         return obj

//     def serialize( self ):
//         return pickle.dumps(self.objects)

//     def deserialize( bb ):
//         ps = ProofStream()
//         ps.objects = pickle.loads(bb)
//         return ps

//     def prover_fiat_shamir( self, num_bytes=32 ):
//         return shake_256(self.serialize()).digest(num_bytes)

//     def verifier_fiat_shamir( self, num_bytes=32 ):
//         return shake_256(pickle.dumps(self.objects[:self.read_index])).digest(num_bytes)
