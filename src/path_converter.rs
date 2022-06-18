use crate::{
    matrix::Position,
    nussinov::RNASequence,
    traceback_paths::TracebackPaths,
    traceback_paths::{TracebackPath, TracebackPathElement},
};

pub type SymbolicPaths = Vec<&'static str>;

pub trait PathConverter {
    fn convert(&self, paths: &TracebackPaths) -> SymbolicPaths;
}

pub struct NussinovPathConverter<'s>(&'s RNASequence);

impl NussinovPathConverter<'_> {
    pub fn new<'a>(sequence: &'a RNASequence) -> NussinovPathConverter {
        NussinovPathConverter::<'a>(sequence)
    }
}

impl PathConverter for NussinovPathConverter<'_> {
    fn convert(&self, paths: &TracebackPaths) -> SymbolicPaths {
        let mut symbolic_paths = vec![];
        for path in paths {
            let mut char_string = vec!['.'; self.0.len()];

            let flatten_path: Vec<Position> = NussinovPathConverter::flatten_path(path);

            for pos in flatten_path.windows(2).into_iter() {
                let first_pos = pos[0];
                let second_pos = pos[1];
                let diagonal = first_pos.is_diagonal_relation(&second_pos);

                if diagonal {
                    let opening_index = first_pos.i;
                    let closing_index = second_pos.j;

                    let _ = std::mem::replace(&mut char_string[opening_index], '(');
                    let _ = std::mem::replace(&mut char_string[closing_index], ')');
                }
            }

            let string: String = char_string.into_iter().collect();
            symbolic_paths.push(string)
        }

        symbolic_paths
            .into_iter()
            .map(|s| NussinovPathConverter::leak_string(s))
            .collect()
    }
}

impl NussinovPathConverter<'_> {
    fn leak_string(s: String) -> &'static str {
        Box::leak(s.into_boxed_str())
    }

    fn flatten_path(path: &TracebackPath) -> Vec<Position> {
        let mut position_path = vec![];

        for path_element in path {
            match path_element {
                TracebackPathElement::Single(p) => position_path.push(*p),
                TracebackPathElement::Decomposition(f, s) => {
                    let mut decomp_path = NussinovPathConverter::flatten_path(f);
                    decomp_path.append(&mut NussinovPathConverter::flatten_path(s));
                    position_path.append(&mut decomp_path)
                }
            }
        }

        position_path
    }
}
