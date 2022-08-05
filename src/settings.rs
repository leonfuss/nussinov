use clap::Parser;

#[derive(Parser)]
#[clap(name = "Nussinov RNA Algorithm")]
#[clap(author = "Leon Fuss <hello@leonfuss.me")]
pub struct Settings {
    #[clap(short, long, value_parser = file_exists)]
    pub file: Option<String>,

    #[clap(short, long, value_parser)]
    pub sequence: Option<String>,
}

fn file_exists(s: &str) -> Result<String, String> {
    if std::path::Path::new(s).exists() {
        return Ok(s.into());
    }

    Err("The provided file path doesn't exist".into())
}
