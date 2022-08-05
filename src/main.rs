use std::process::exit;

use clap::Parser;
use nussinov_cli::{nussinov::Nussinov, settings};

fn main() {
    let settings = settings::Settings::parse();

    let sequence;
    if settings.sequence.is_some() {
        sequence = settings.sequence.unwrap();
    } else if settings.file.is_some() {
        let path = settings.file.unwrap();
        sequence = std::fs::read_to_string(path).expect("File could not be read");
    } else {
        println!("Please call nussinov either with a sequence or a file path");
        exit(0);
    }

    let mut nussinov = Nussinov::new(&sequence, 1);
    nussinov.run();
}
