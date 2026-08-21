use std::env;
use std::path::PathBuf;

fn main() -> spring_rs::Result<()> {
    let mut args = env::args_os().skip(1);
    let genome = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing genome FASTA".into()))?,
    );
    let prefix = PathBuf::from(
        args.next()
            .ok_or_else(|| spring_rs::Error::InvalidInput("missing index prefix".into()))?,
    );
    if args.next().is_some() {
        return Err(spring_rs::Error::InvalidInput(
            "expected genome FASTA and index prefix".into(),
        ));
    }
    let _index = spring_rs::build_genome_index(genome, prefix)?;
    Ok(())
}
