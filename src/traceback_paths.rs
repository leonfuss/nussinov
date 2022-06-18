use std::collections::HashMap;

use crate::matrix::{Matrix, MatrixNode, Position, TraceType};

pub type TracebackPaths = Vec<TracebackPath>;

pub type TracebackPath = Vec<TracebackPathElement>;

#[derive(Debug, Clone)]
pub enum TracebackPathElement {
    Single(Position),
    Decomposition(Vec<TracebackPathElement>, Vec<TracebackPathElement>),
}

struct TracebackTree(TracebackNode);

struct TracebackNode {
    visited: bool,
    position: Position,
    children: Trace,
    decomposition: HashMap<Position, Vec<TracebackPathElement>>,
}

type Trace = Vec<TracebackType>;

enum TracebackType {
    Complementary(TracebackNode),
    Unpaired(TracebackNode),
    Decomposition(TracebackNode, TracebackNode),
}

pub trait TracebackPathsBuilder {
    fn build(&self, matrix: &Matrix) -> TracebackPaths;
}

pub struct NussinovTracebackPathsBuilder();

impl NussinovTracebackPathsBuilder {
    pub fn new() -> NussinovTracebackPathsBuilder {
        NussinovTracebackPathsBuilder()
    }
}

impl TracebackPathsBuilder for NussinovTracebackPathsBuilder {
    fn build(&self, matrix: &Matrix) -> TracebackPaths {
        let mut tree = TracebackTree::new(matrix);

        tree.traverse_paths()
    }
}

impl TracebackTree {
    fn new(matrix: &Matrix) -> TracebackTree {
        let matrix_root = matrix.root();
        let root = TracebackTree::get_children(&matrix_root, matrix);

        TracebackTree(root)
    }

    fn get_children(node: &MatrixNode, matrix: &Matrix) -> TracebackNode {
        let children = node
            .trace
            .iter()
            .map(|trace_type| match trace_type {
                TraceType::Complementary(pos) => {
                    let child_node = &matrix[*pos];
                    let node = TracebackTree::get_children(&child_node, matrix);
                    TracebackType::Complementary(node)
                }
                TraceType::Decomposition(pos1, pos2) => {
                    let child_node1 = &matrix[*pos1];
                    let child_node2 = &matrix[*pos2];
                    let node1 = TracebackTree::get_children(child_node1, matrix);
                    let node2 = TracebackTree::get_children(child_node2, matrix);
                    TracebackType::Decomposition(node1, node2)
                }
                TraceType::Unpaired(pos) => {
                    let child_node = &matrix[*pos];
                    let node = TracebackTree::get_children(&child_node, matrix);
                    TracebackType::Unpaired(node)
                }
            })
            .collect();

        TracebackNode {
            visited: false,
            position: node.position,
            children,
            decomposition: HashMap::new(),
        }
    }

    fn traverse_paths(&mut self) -> TracebackPaths {
        TracebackTree::traverse_node_paths(&mut self.0)
    }

    fn traverse_node_paths(node: &mut TracebackNode) -> TracebackPaths {
        let mut paths = vec![];
        while !node.visited {
            if let Some(path) = TracebackTree::traverse(node) {
                let path = path.into_iter().rev().collect();
                paths.push(path);
            }
        }
        paths
    }

    fn traverse(node: &mut TracebackNode) -> Option<TracebackPath> {
        if node.children.is_empty() {
            return Some(vec![]);
        }
        let position = &node.position;
        let optional_child = node
            .children
            .iter_mut()
            .find(|trace_type| match trace_type {
                TracebackType::Complementary(n) if !n.visited => true,
                TracebackType::Unpaired(n) if !n.visited => true,
                TracebackType::Decomposition(f, _) => {
                    node.decomposition.is_empty()
                        || (node.decomposition.contains_key(&f.position)
                            && node.decomposition.get(&f.position).unwrap().len() > 0)
                        || !node.decomposition.contains_key(&f.position)
                }
                _ => false,
            });
        if let None = optional_child {
            node.visited = true;
            return None;
        }

        let trace;
        match &mut optional_child.unwrap() {
            TracebackType::Complementary(n) => {
                trace = TracebackTree::traverse(n);
            }
            TracebackType::Unpaired(n) => {
                trace = TracebackTree::traverse(n);
            }
            TracebackType::Decomposition(f, s) => {
                if !node.decomposition.contains_key(&f.position) {
                    let first = TracebackTree::traverse_node_paths(f);
                    let second = TracebackTree::traverse_node_paths(s);

                    let mut decomposition = vec![];

                    for i in 0..first.len() {
                        for j in 0..second.len() {
                            let trace = TracebackPathElement::Decomposition(
                                first[i].clone(),
                                second[j].clone(),
                            );
                            decomposition.push(trace);
                        }
                    }
                    node.decomposition.insert(f.position, decomposition);
                }

                trace = match node.decomposition.get_mut(&f.position).unwrap().pop() {
                    Some(t) => Some(vec![t]),
                    None => None,
                };
            }
        }

        match trace {
            Some(mut trace) => {
                if trace.is_empty() {
                    node.visited = true;
                }
                trace.push(TracebackPathElement::Single(*position));
                Some(trace)
            }
            None => None,
        }
    }
}
