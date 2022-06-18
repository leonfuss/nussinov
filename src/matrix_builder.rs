use crate::{
    matrix::{Matrix, MatrixNode, Position, Trace, TraceType},
    nussinov::RNASequence,
};

pub trait MatrixBuilder {
    fn fill<'m>(&self, matrix: &'m mut Matrix);
}

#[derive(Debug)]
pub struct NussinovMatrixBuilder<'s> {
    sequence: &'s RNASequence,
    minimal_loop_length: usize,
}

struct DiagonalMatrixIterator {
    columns: usize,
    rows: usize,
    position: usize,
    gap: usize,
}

pub enum UnpairedType {
    Left,
    Bottom,
}

impl NussinovMatrixBuilder<'_> {
    pub fn new(sequence: &RNASequence, minimal_loop_length: usize) -> NussinovMatrixBuilder {
        NussinovMatrixBuilder {
            sequence,
            minimal_loop_length,
        }
    }

    fn determine_max(&self, matrix: &Matrix, pos: &Position) -> (Trace, usize) {
        let mut possible_traces = vec![
            self.get_complementary(matrix, pos),
            NussinovMatrixBuilder::get_unpaired(matrix, pos, UnpairedType::Left),
            NussinovMatrixBuilder::get_unpaired(matrix, pos, UnpairedType::Bottom),
        ];
        if let Some(mut t) = self.get_decomposition(matrix, pos) {
            possible_traces.append(&mut t);
        }

        let mut max = 0;

        for pt in &possible_traces {
            if let Some(x) = pt {
                let (_, value) = x;
                if max < *value {
                    max = *value;
                }
            }
        }

        let (trace, _): (Vec<TraceType>, Vec<usize>) = possible_traces
            .into_iter()
            .filter_map(|x| x)
            .filter(|(_, v)| *v == max)
            .unzip();

        (trace, max)
    }

    fn get_complementary(&self, matrix: &Matrix, pos: &Position) -> Option<(TraceType, usize)> {
        let is_complement = self.sequence.is_complement(pos).unwrap_or_default();
        if is_complement && (self.minimal_loop_length + pos.i) < (pos.j - 1) {
            let node = &matrix[pos.get_complementary()];
            let trace_type = TraceType::Complementary(node.position);
            return Some((trace_type, node.value + 1));
        }
        None
    }

    fn get_decomposition(
        &self,
        matrix: &Matrix,
        pos: &Position,
    ) -> Option<Vec<Option<(TraceType, usize)>>> {
        let mut value_max = None;
        let mut nodes_max: Option<Vec<(Position, Position)>> = None;
        for k in (pos.i + 2)..(pos.j - 1) {
            let (pos1, pos2) = pos.get_decomposition(k);
            let (node1, node2) = (&matrix[pos1], &matrix[pos2]);
            let value: usize = node1.value + node2.value;

            match value_max {
                Some(v) => {
                    if value == v {
                        nodes_max
                            .as_mut()
                            .unwrap()
                            .push((node1.position, node2.position));
                    } else if value > v {
                        value_max = Some(value);
                        nodes_max = Some(vec![(node1.position, node2.position)]);
                    }
                }
                None => {
                    value_max = Some(value);
                    nodes_max = Some(vec![(node1.position, node2.position)]);
                }
            }
        }

        if let Some(v) = value_max {
            let nodes = nodes_max.unwrap();
            let traces: Vec<Option<(TraceType, usize)>> = nodes
                .into_iter()
                .map(|(pos1, pos2)| Some((TraceType::Decomposition(pos1, pos2), v)))
                .collect();
            return Some(traces);
        }

        None
    }

    fn get_unpaired(
        matrix: &Matrix,
        pos: &Position,
        unpaired_type: UnpairedType,
    ) -> Option<(TraceType, usize)> {
        let node = &matrix[pos.get_unpaired(unpaired_type)];
        let trace_type = TraceType::Unpaired(node.position);
        return Some((trace_type, node.value));
    }
}

impl MatrixBuilder for NussinovMatrixBuilder<'_> {
    fn fill<'m>(&self, matrix: &'m mut Matrix) {
        let j = matrix.columns();
        let i = matrix.rows();
        let diagonal_iter = DiagonalMatrixIterator::new(j, i);

        for position in diagonal_iter {
            let (trace, value) = self.determine_max(matrix, &position);
            let node = MatrixNode {
                trace,
                position,
                value,
            };
            matrix[position] = node;
        }
    }
}

impl DiagonalMatrixIterator {
    pub fn new(rows: usize, columns: usize) -> DiagonalMatrixIterator {
        DiagonalMatrixIterator {
            columns,
            rows,
            position: 0,
            gap: 1,
        }
    }
}

impl Iterator for DiagonalMatrixIterator {
    type Item = Position;

    fn next(&mut self) -> Option<Self::Item> {
        if (self.position + self.gap) >= self.columns {
            self.position = 0;
            self.gap += 1;
        }
        if self.gap >= self.columns {
            return None;
        }
        let i = self.position;
        let j = i + self.gap;

        if j >= self.rows {
            return None;
        }

        self.position += 1;
        Some(Position::from(i, j))
    }
}
