use crate::matrix_builder::MatrixBuilder;
use crate::matrix_builder::UnpairedType;
use array2d::Array2D;
use std::fmt::Debug;
use std::fmt::Display;
use std::ops::Index;
use std::ops::IndexMut;

pub struct Matrix(Array2D<MatrixNode>);

#[derive(Clone, Debug)]
pub struct MatrixNode {
    pub value: usize,
    pub position: Position,
    pub trace: Trace,
}

pub type Trace = Vec<TraceType>;

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub struct Position {
    pub i: usize,
    pub j: usize,
}

#[derive(Clone, Debug)]
pub enum TraceType {
    Complementary(Position),
    Unpaired(Position),
    Decomposition(Position, Position),
}

impl Display for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();

        for i in 0..self.rows() {
            for j in 0..self.columns() {
                let node = &self[Position::from(i, j)];
                let number_string = format!("{} ", node.value);
                s.push_str(&number_string);
            }
            s.push('\n');
        }
        writeln!(f, "{}", s).unwrap();

        Ok(())
    }
}

impl Matrix {
    pub fn new(size: usize) -> Matrix {
        let mut counter = 0;
        let generator = || {
            let j = counter % (size + 1);
            let i = (counter - j) / (size + 1);
            counter += 1;

            MatrixNode {
                value: 0,
                position: Position::from(i, j),
                trace: Trace::new(),
            }
        };
        let array = Array2D::filled_by_row_major(generator, size + 1, size + 1);

        Matrix(array)
    }

    pub fn columns(&self) -> usize {
        self.0.row_len()
    }

    pub fn rows(&self) -> usize {
        self.0.column_len()
    }

    pub fn fill(&mut self, builder: &dyn MatrixBuilder) {
        builder.fill(self);
    }

    pub fn root(&self) -> &MatrixNode {
        &self[Position::from(0, self.columns() - 1)]
    }
}

impl Index<Position> for Matrix {
    type Output = MatrixNode;

    fn index(&self, index: Position) -> &Self::Output {
        &self.0[index.into()]
    }
}

impl IndexMut<Position> for Matrix {
    fn index_mut(&mut self, index: Position) -> &mut Self::Output {
        &mut self.0[index.into()]
    }
}

impl Position {
    pub fn from(i: usize, j: usize) -> Position {
        Position { i, j }
    }

    pub fn is_diagonal_relation(&self, pos: &Position) -> bool {
        (self.i + 1) == pos.i && self.j == (pos.j + 1)
    }

    pub fn is_diagonal(&self) -> bool {
        self.i == self.j
    }

    pub fn get_unpaired(&self, unpaired_type: UnpairedType) -> Position {
        match unpaired_type {
            UnpairedType::Left => Position {
                i: self.i,
                j: self.j - 1,
            },
            UnpairedType::Bottom => Position {
                i: self.i + 1,
                j: self.j,
            },
        }
    }

    pub fn get_decomposition(&self, k: usize) -> (Position, Position) {
        let pos1 = Position { i: self.i, j: k };
        let pos2 = Position { i: k, j: self.j };
        (pos1, pos2)
    }

    pub fn get_complementary(&self) -> Position {
        Position {
            i: self.i + 1,
            j: self.j - 1,
        }
    }
}

impl From<&Position> for (usize, usize) {
    fn from(pos: &Position) -> Self {
        (pos.i, pos.j)
    }
}

impl From<Position> for (usize, usize) {
    fn from(pos: Position) -> Self {
        (pos.i, pos.j)
    }
}

impl Display for Position {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "({},{})", self.i, self.j)
    }
}
