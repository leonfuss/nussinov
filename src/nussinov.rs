use std::ops::Index;

use crate::{
    matrix::{Matrix, Position},
    matrix_builder::NussinovMatrixBuilder,
    path_converter::NussinovPathConverter,
    path_converter::PathConverter,
    traceback_paths::{NussinovTracebackPathsBuilder, TracebackPathsBuilder},
};

pub struct Nussinov {
    minimal_loop_length: usize,
    sequence: RNASequence,
    matrix: Matrix,
}

// TODO remove Clone
#[derive(Debug, Clone)]
pub struct RNASequence(Vec<char>);

impl Nussinov {
    pub fn new(sequence: &str, minimal_loop_length: usize) -> Nussinov {
        let deserialized_sequence = RNASequence::new(sequence);
        match deserialized_sequence {
            Ok(s) => Nussinov {
                minimal_loop_length,
                matrix: Matrix::new(s.len()),
                sequence: s,
            },
            Err(_) => panic!("The given RNA-sequence is invalid"),
        }
    }

    pub fn run(&mut self) {
        println!("Analysing sequence: {:?}", self.sequence);
        println!();
        let matrix_builder = NussinovMatrixBuilder::new(&self.sequence, self.minimal_loop_length);
        self.matrix.fill(&matrix_builder);

        print!("{}", self.matrix);

        let traceback_builder = NussinovTracebackPathsBuilder::new();
        let paths = traceback_builder.build(&self.matrix);

        let path_converter = NussinovPathConverter::new(&self.sequence);
        let symbolic_paths = path_converter.convert(&paths);

        println!("{:#?}", symbolic_paths);
    }
}

impl RNASequence {
    const VALID_CHARS: [char; 4] = ['A', 'U', 'G', 'C'];
    const COMPLEMENTS: [[char; 2]; 3] = [['A', 'U'], ['G', 'C'], ['G', 'U']];

    pub fn new(rna_sequence: &str) -> Result<RNASequence, ()> {
        let sequence = rna_sequence.to_uppercase();
        if !RNASequence::is_valid(&sequence) {
            return Err(());
        }

        Ok(RNASequence(sequence.chars().collect()))
    }

    fn is_valid(sequence: &str) -> bool {
        sequence
            .chars()
            .all(|c| RNASequence::VALID_CHARS.contains(&c))
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_complement(&self, position: &Position) -> Result<bool, ()> {
        let size = self.len();
        let i = position.i;
        let j = position.j - 1;

        if i >= size || j >= size {
            return Err(());
        }
        let first_nucleotid = self.0[j];
        let second_nucleotid = self.0[i];

        let complement = RNASequence::COMPLEMENTS
            .iter()
            .map(|c| {
                (c.contains(&first_nucleotid) && c.contains(&second_nucleotid))
                    && first_nucleotid != second_nucleotid
            })
            .fold(false, |c, n| c || n);

        Ok(complement)
    }
}

impl Index<usize> for &RNASequence {
    type Output = char;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
