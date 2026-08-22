use assert_cmd::Command;
use predicates::prelude::*;
use seqarc::gzip_writer;
use std::fs::File;
use std::io::Write;

fn write_gzip(path: &std::path::Path, contents: &[u8]) {
    let mut writer = gzip_writer(File::create(path).unwrap(), 1);
    writer.write_all(contents).unwrap();
    writer.finish().unwrap();
}

#[test]
fn help_lists_aligned_workflow_commands() {
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("--help")
        .assert()
        .success()
        .stdout(predicate::str::contains("index"))
        .stdout(predicate::str::contains("compress"))
        .stdout(predicate::str::contains("decompress"))
        .stdout(predicate::str::contains("validate"));
}

#[test]
fn aligned_cli_round_trip_and_validation() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq.gz");
    let second_path = directory.path().join("reads.2.fastq.gz");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_first = directory.path().join("decoded.1.fastq");
    let decoded_second = directory.path().join("decoded.2.fastq");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();

    let first_sequence = &reference[20..60];
    let second_sequence = reference[20..60]
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect::<Vec<_>>();
    let first_fastq = format!(
        "@pair0/1\n{}\n+\n{}\n@pair1/1\n{}\n+\n{}\n",
        String::from_utf8_lossy(first_sequence),
        "I".repeat(first_sequence.len()),
        "N".repeat(40),
        "#".repeat(40),
    );
    let second_fastq = format!(
        "@pair0/2\n{}\n+\n{}\n@pair1/2\n{}\n+\n{}\n",
        String::from_utf8_lossy(&second_sequence),
        "J".repeat(second_sequence.len()),
        "N".repeat(40),
        "$".repeat(40),
    );
    write_gzip(&first_path, first_fastq.as_bytes());
    write_gzip(&second_path, second_fastq.as_bytes());

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--output")
        .arg(&archive_path)
        .arg("--chunk-pairs")
        .arg("1")
        .assert()
        .success()
        .stderr(predicate::str::contains("stage=index_load done"))
        .stderr(predicate::str::contains("stage=map_and_spool done"))
        .stderr(predicate::str::contains("stage=reference done"))
        .stderr(predicate::str::contains("stage=chunks done"))
        .stderr(predicate::str::contains("stage=archive done"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .assert()
        .success()
        .stdout(predicate::str::contains("validated 2 pairs byte-for-byte"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .assert()
        .success();

    assert_eq!(std::fs::read_to_string(decoded_first).unwrap(), first_fastq);
    assert_eq!(
        std::fs::read_to_string(decoded_second).unwrap(),
        second_fastq
    );
}
