use std::env;
use std::path::PathBuf;

fn main() -> seqarc::Result<()> {
    let mut args = env::args_os().skip(1);
    let genome = PathBuf::from(
        args.next()
            .ok_or_else(|| seqarc::Error::InvalidInput("missing genome FASTA".into()))?,
    );
    let index_file = PathBuf::from(
        args.next()
            .ok_or_else(|| seqarc::Error::InvalidInput("missing index file".into()))?,
    );
    if args.next().is_some() {
        return Err(seqarc::Error::InvalidInput(
            "expected genome FASTA and index file".into(),
        ));
    }
    let _index = seqarc::build_genome_index(genome, index_file)?;
    Ok(())
}
