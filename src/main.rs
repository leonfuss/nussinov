use crate::nussinov::Nussinov;

#[macro_use]
extern crate lazy_static;

mod matrix;
mod nussinov;

fn main() {
    let mut nussinov = Nussinov::new("UAGA", 1);
    nussinov.run();
}
