use crate::matrix::NussinovMatrix;
use std::ops::Index;

pub struct Nussinov {
    matrix: NussinovMatrix,
}

#[derive(Debug, Clone)]
pub struct RNASequence {
    inner: Vec<char>,
}

impl Nussinov {
    pub fn new(sequence: &str, minimal_loop_length: usize) -> Nussinov {
        let deserialized_sequence = RNASequence::new(sequence);
        match deserialized_sequence {
            Ok(s) => Nussinov {
                matrix: NussinovMatrix::new(minimal_loop_length, s.len(), s),
            },
            Err(_) => panic!("The given RNA-sequence is invalid"),
        }
    }

    pub fn run(&mut self) {
        self.matrix.solve();
        println!("{}", self.matrix);
        let tree = self.matrix.traceback();
        dbg!(&tree);
        let paths = tree.traverse_all_paths();
        tree.to_symbols();
    }
}

impl RNASequence {
    const VALID_CHARS: [char; 4] = ['A', 'U', 'G', 'C'];

    pub fn new(rna_sequence: &str) -> Result<RNASequence, ()> {
        let sequence = rna_sequence.to_uppercase();
        if !RNASequence::is_valid(&sequence) {
            return Err(());
        }
        Ok(RNASequence {
            inner: sequence.chars().collect(),
        })
    }

    fn is_valid(sequence: &str) -> bool {
        sequence
            .chars()
            .all(|c| RNASequence::VALID_CHARS.contains(&c))
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

impl Index<usize> for &RNASequence {
    type Output = char;

    fn index(&self, index: usize) -> &Self::Output {
        &self.inner[index]
    }
}
