use nussinov_cli::nussinov::Nussinov;

fn main() {
    let mut nussinov = Nussinov::new("GGUCCGACGUAAUA", 1);
    nussinov.run();
}
