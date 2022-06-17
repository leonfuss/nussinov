use crate::nussinov::RNASequence;
use array2d::Array2D;
use std::collections::HashSet;

pub struct NussinovMatrix {
    sequence: RNASequence,
    minimal_loop_length: usize,
    size: usize,
    inner: Array2D<MatrixNode>,
}

#[derive(Default, Debug, Clone)]
struct MatrixNode {
    value: u8,
    position: (usize, usize),
    max: Vec<MaxType>,
}

#[derive(Debug)]
enum TracebackType {
    Complementary(TracebackNode),
    Decomposition(TracebackNode, TracebackNode),
    Unpaired(TracebackNode),
}

#[derive(Clone, Debug)]
enum MaxType {
    Complementary(usize, usize),
    Decomposition {
        first: (usize, usize),
        second: (usize, usize),
    },
    Unpaired(usize, usize),
}

#[derive(Debug)]
pub struct TracebackTree {
    root: TracebackNode,
    sequence: RNASequence,
}

#[derive(Debug)]
struct TracebackNode {
    value: u8,
    position: (usize, usize),
    children: Vec<TracebackType>,
}

impl std::fmt::Display for NussinovMatrix {
    fn fmt(&self, _f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..(self.size + 1) {
            for j in 0..(self.size + 1) {
                let value = self.inner[(i, j)].value;
                print!("{}", value);
            }
            println!();
        }
        Ok(())
    }
}

impl std::fmt::Display for TracebackTree {
    fn fmt(&self, _f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}

lazy_static! {
    static ref COMPLEMENTARY: Vec<HashSet<char>> = {
        vec![
            HashSet::from(['A', 'U']),
            HashSet::from(['G', 'C']),
            HashSet::from(['G', 'U']),
        ]
    };
}

impl NussinovMatrix {
    pub fn new(minimal_loop_length: usize, size: usize, sequence: RNASequence) -> NussinovMatrix {
        let mut counter = 0;
        let array_iter = || {
            let j = counter % (size + 1);
            let i = (counter - j) / (size + 1);
            counter += 1;

            MatrixNode {
                value: 0,
                position: (i, j),
                max: vec![],
            }
        };

        NussinovMatrix {
            sequence,
            minimal_loop_length,
            size,
            inner: Array2D::filled_by_row_major(array_iter, size + 1, size + 1),
        }
    }

    pub fn solve(&mut self) {
        for j in 1..(self.size + 1) {
            let mut i = 0;
            while j + i < (self.size + 1) {
                self.inner[(i, j + i)] = self.get_node(i, j + i);
                i += 1;
            }
        }
    }

    fn get_node(&self, i: usize, j: usize) -> MatrixNode {
        let mut max_nodes = vec![];

        if i >= j || j == 0 {
            panic!("create_node should only be called by NussinovMatrix::solve()");
        }

        let complementary_node = &self.inner[(i + 1, j - 1)];
        let is_complement = &self.is_complementary(i, j, &self.sequence);
        if *is_complement {
            max_nodes.push((
                complementary_node.value + 1,
                MaxType::Complementary(i + 1, j - 1),
            ));
        }

        if let Some((v, t)) = self.get_decomposition(i, j) {
            max_nodes.push((v, t))
        }

        max_nodes.push((self.inner[(i + 1, j)].value, MaxType::Unpaired(i + 1, j)));

        max_nodes.push((self.inner[(i, j - 1)].value, MaxType::Unpaired(i, j - 1)));

        let (value, _) = *max_nodes
            .iter()
            .max_by_key(|&&(v, _)| v)
            .expect("Max_nodes can not be empty");
        let max = max_nodes
            .iter()
            .filter(|&&(v, _)| v == value)
            .map(|(_, t)| t.clone())
            .collect();

        MatrixNode {
            value,
            position: (i, j),
            max,
        }
    }

    fn get_decomposition(&self, i: usize, j: usize) -> Option<(u8, MaxType)> {
        let mut max = None;
        let mut first = (0, 0);
        let mut second = (0, 0);

        for k in (i + 1)..(j - 1) {
            let n1 = &self.inner[(i, k)];
            let n2 = &self.inner[(k + 1, j)];
            let value = n1.value + n2.value;

            // always takes the last max
            if max.is_none() || value >= max.unwrap() {
                max = Some(value);
                first = n1.position;
                second = n2.position
            }
        }

        if max == None {
            return None;
        }

        Some((max.unwrap(), MaxType::Decomposition { first, second }))
    }

    fn is_complementary(&self, i: usize, j: usize, sequence: &RNASequence) -> bool {
        if j == 0 || (i + self.minimal_loop_length) >= (j - 1) {
            return false;
        }

        let base: HashSet<char> = HashSet::from([sequence[i], sequence[j - 1]]);

        for complement in COMPLEMENTARY.iter() {
            if complement.symmetric_difference(&base).count() == 0 {
                return true;
            }
        }

        false
    }

    pub fn traceback(&self) -> TracebackTree {
        TracebackTree::new(&self)
    }
}

impl TracebackTree {
    fn new(matrix: &NussinovMatrix) -> TracebackTree {
        let matrix_root = &matrix.inner[(0, matrix.size)];
        let root = TracebackTree::get_children(&matrix_root, &matrix.inner);

        TracebackTree {
            root,
            sequence: matrix.sequence.clone(),
        }
    }

    fn get_children(node: &MatrixNode, matrix: &Array2D<MatrixNode>) -> TracebackNode {
        let mut children = vec![];
        for child in node.max.iter() {
            match child {
                MaxType::Complementary(i, j) => {
                    let node = TracebackTree::get_children(&matrix[(*i, *j)], &matrix);
                    let t = TracebackType::Complementary(node);
                    children.push(t);
                }
                MaxType::Decomposition { first, second } => {
                    let node1 = TracebackTree::get_children(&matrix[*first], &matrix);
                    let node2 = TracebackTree::get_children(&matrix[*second], &matrix);
                    let t = TracebackType::Decomposition(node1, node2);
                    children.push(t);
                }
                MaxType::Unpaired(i, j) => {
                    let node = TracebackTree::get_children(&matrix[(*i, *j)], &matrix);
                    let t = TracebackType::Unpaired(node);
                    children.push(t);
                }
            }
        }

        TracebackNode {
            value: node.value,
            position: node.position,
            children,
        }
    }

    pub fn to_symbols(&self) {
        let paths = self.traverse_all_paths();

        dbg!(&paths);

        for mut path in paths {
            let mut chars = vec!['.'; self.sequence.len()];
            path.reverse();

            for i in 0..(path.len() - 1) {
                let first_pos = path[i];
                let second_pos = path[i + 1];
                if TracebackTree::is_diagonal(first_pos, second_pos) {
                    let _ = std::mem::replace(&mut chars[first_pos.1 - 1], ')');
                    let _ = std::mem::replace(&mut chars[first_pos.0], '(');
                }
            }
            let string = chars.iter().clone().collect::<String>();
            dbg!(&string);
            println!("{}", string);
        }
    }

    fn is_diagonal(first: (usize, usize), second: (usize, usize)) -> bool {
        let (a, b) = first;
        let (c, d) = second;
        let result = (a == (c - 1)) && (b == (d + 1));
        dbg!(first, second, result);
        result
    }

    pub fn traverse_all_paths(&self) -> Vec<Vec<(usize, usize)>> {
        let mut visited: Vec<(usize, usize)> = vec![];
        let mut paths = vec![];

        while !visited.contains(&self.root.position) {
            match TracebackTree::traverse(&self.root, &visited, &paths) {
                Ok(p) => {
                    let fist = p.first();
                    match fist {
                        Some(&f) => visited.push(f),
                        None => {}
                    }
                    paths.push(p);
                }
                Err(v) => visited.push(v),
            }
        }

        paths
    }

    fn traverse(
        root: &TracebackNode,
        visited: &Vec<(usize, usize)>,
        paths: &Vec<Vec<(usize, usize)>>,
    ) -> Result<Vec<(usize, usize)>, (usize, usize)> {
        println!("-- root: {:?}, visited: {:?}", &root.position, &visited);

        let (x, y) = root.position;
        if x == y {
            return Ok(vec![]);
        }

        for child in root.children.iter() {
            let mut node = None;
            match child {
                TracebackType::Complementary(e) if !visited.contains(&e.position) => {
                    println!("TT:Complementary: {:?}", &e.position);
                    node = Some(e);
                }
                TracebackType::Decomposition(a, b) => {
                    if !visited.contains(&a.position) {
                        println!("TT:Decomp-first: {:?}", &a.position);
                        node = Some(a);
                    } else if !visited.contains(&b.position) {
                        println!("TT:Decomp-second: {:?}", &b.position);
                        node = Some(b);
                    }
                }
                TracebackType::Unpaired(e) if (!visited.contains(&e.position)) => {
                    println!("TT:Unpaired: {:?}", &e.position);
                    node = Some(e);
                }
                _ => {}
            };

            if let Some(e) = node {
                return match TracebackTree::traverse(&e, visited, paths) {
                    Ok(mut trace) => {
                        trace.push(root.position);
                        println!("clearing stack -- trace: {:?}", trace);
                        Ok(trace)
                    }
                    Err(e) => Err(e),
                };
            }
        }

        println!("End Reached -- returning");
        return Err(root.position);
    }

    fn travel_condition(
        current_position: (usize, usize),
        next_position: (usize, usize),
        visited: &Vec<(usize, usize)>,
        paths: &Vec<Vec<(usize, usize)>>,
    ) -> bool {
        let both_visited = visited
            .windows(2)
            .any(|x| x[0] == current_position && x[1] == next_position);
        let both_traveled = paths
            .iter()
            .map(|path| {
                path.windows(2)
                    .any(|x| x[0] == current_position && x[1] == next_position)
            })
            .fold(false, |x, y| x || y);
        !(both_visited && both_traveled)
    }
}
