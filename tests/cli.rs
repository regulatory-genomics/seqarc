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

fn write_zstd(path: &std::path::Path, contents: &[u8]) {
    let mut writer = zstd::stream::write::Encoder::new(File::create(path).unwrap(), 1).unwrap();
    writer.write_all(contents).unwrap();
    writer.finish().unwrap();
}

fn assert_directory_empty(path: &std::path::Path) {
    assert!(path.is_dir());
    assert!(std::fs::read_dir(path).unwrap().next().is_none());
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
    let reference_path = directory.path().join("reference.fa.zst");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq.gz");
    let second_path = directory.path().join("reads.2.fastq.gz");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_first = directory.path().join("decoded.1.fastq");
    let decoded_second = directory.path().join("decoded.2.fastq");
    let temp_root = directory.path().join("temporary");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    write_zstd(
        &reference_path,
        &[b">chr1\n".as_slice(), reference, b"\n"].concat(),
    );

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
    write_zstd(&second_path, second_fastq.as_bytes());

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .arg("--temp-dir")
        .arg(&temp_root)
        .assert()
        .success();
    assert_directory_empty(&temp_root);

    let failed_archive = directory.path().join("failed.seqarc");
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(directory.path().join("missing.fastq.gz"))
        .arg("--output")
        .arg(&failed_archive)
        .arg("--temp-dir")
        .arg(&temp_root)
        .assert()
        .failure();
    assert_directory_empty(&temp_root);
    assert!(!failed_archive.exists());

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
        .arg("--chunk-records")
        .arg("1")
        .arg("--id-prefix-chars")
        .arg("5")
        .arg("--temp-dir")
        .arg(&temp_root)
        .assert()
        .success()
        .stderr(predicate::str::contains("stage=index_load done"))
        .stderr(predicate::str::contains("stage=map_and_spool done"))
        .stderr(predicate::str::contains("stage=reference done"))
        .stderr(predicate::str::contains("stage=chunks done"))
        .stderr(predicate::str::contains("stage=archive done"))
        .stderr(predicate::str::contains("% aligned"))
        .stderr(predicate::str::contains("% unaligned"))
        .stderr(predicate::str::contains("% no hit"))
        .stderr(predicate::str::contains("% not compressible"));
    assert_directory_empty(&temp_root);

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("-1")
        .arg(&first_path)
        .arg("-2")
        .arg(&second_path)
        .assert()
        .failure();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("-1")
        .arg(&first_path)
        .arg("-2")
        .arg(&second_path)
        .arg("--id-prefix-chars")
        .arg("5")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 2 records against originals (PairedEnd)",
        ));

    let incomplete_output = directory.path().join("incomplete.fastq");
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&incomplete_output)
        .assert()
        .failure();
    assert!(!incomplete_output.exists());

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .arg("--threads")
        .arg("2")
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "stage=decompress start total_chunks=2",
        ))
        .stderr(predicate::str::contains(
            "stage=decompress done decoded_chunks=2 total_chunks=2",
        ));

    let expected_first = first_fastq
        .replace("@pair0/1", "@pair0")
        .replace("@pair1/1", "@pair1");
    let expected_second = second_fastq
        .replace("@pair0/2", "@pair0")
        .replace("@pair1/2", "@pair1");
    assert_eq!(
        std::fs::read_to_string(decoded_first).unwrap(),
        expected_first
    );
    assert_eq!(
        std::fs::read_to_string(decoded_second).unwrap(),
        expected_second
    );
}

#[test]
fn index_cleans_custom_temp_directory_after_failure() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("empty.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let temp_root = directory.path().join("temporary");
    std::fs::write(&reference_path, b">empty\n").unwrap();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("index")
        .arg(&reference_path)
        .arg("--output")
        .arg(&index_path)
        .arg("--temp-dir")
        .arg(&temp_root)
        .assert()
        .failure();

    assert_directory_empty(&temp_root);
    assert!(!index_path.exists());
}

#[test]
fn single_end_cli_round_trip_and_validation() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_path = directory.path().join("decoded.fastq");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let sequence = &reference[20..60];
    let fastq = format!(
        "@read0\n{}\n+\n{}\n@read1\n{}\n+\n{}\n",
        String::from_utf8_lossy(sequence),
        "I".repeat(sequence.len()),
        "N".repeat(40),
        "#".repeat(40),
    );
    std::fs::write(&input_path, fastq.as_bytes()).unwrap();

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
        .arg(&input_path)
        .arg("--output")
        .arg(&archive_path)
        .arg("--chunk-records")
        .arg("1")
        .assert()
        .success()
        .stderr(predicate::str::contains("layout=SingleEnd"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&input_path)
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 2 records against originals (SingleEnd)",
        ));

    let invalid_first = directory.path().join("invalid.1.fastq");
    let invalid_second = directory.path().join("invalid.2.fastq");
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&invalid_first)
        .arg("--mate-2")
        .arg(&invalid_second)
        .assert()
        .failure();
    assert!(!invalid_first.exists());
    assert!(!invalid_second.exists());

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_path)
        .assert()
        .success();

    assert_eq!(std::fs::read_to_string(decoded_path).unwrap(), fastq);
}

#[test]
fn embedded_index_ranges_round_trip_across_mates() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let second_path = directory.path().join("reads.2.fastq");
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
    let first_fastq = b"@pair0/1\nAAAACCCCGGGGTTTTACGT\n+\nABCDEFGHIJKLMNOPQRST\n@pair1/1\nAAAACCCCGGAGTTTTACGT\n+\nTSRQPONMLKJIHGFEDCBA\n";
    let second_fastq = b"@pair0/2\nTGCATTTTAAAACCCCGGTA\n+\nIIIIIIIIIIIIIIIIIIII\n@pair1/2\nTGCATTTGAAAACCCCGGTA\n+\nJJJJJJJJJJJJJJJJJJJJ\n";
    std::fs::write(&first_path, first_fastq).unwrap();
    std::fs::write(&second_path, second_fastq).unwrap();

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
        .arg("--index-range")
        .arg("R1:9:12")
        .arg("--index-range")
        .arg("R2:5:8")
        .arg("--index-range")
        .arg("R1:1:4")
        .arg("--chunk-records")
        .arg("1")
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .success();

    let archive = std::fs::read(&archive_path).unwrap();
    assert_eq!(u32::from_le_bytes(archive[9..13].try_into().unwrap()), 16);

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .arg("--threads")
        .arg("2")
        .assert()
        .success();

    assert_eq!(std::fs::read(&decoded_first).unwrap(), first_fastq);
    assert_eq!(std::fs::read(&decoded_second).unwrap(), second_fastq);
}
