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

fn pseudo_random_dna(length: usize) -> Vec<u8> {
    let mut state = 0x9e37_79b9u32;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
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
fn rustar_cli_builds_single_file_index_and_compresses_spliced_read() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let gtf_path = directory.path().join("genes.gtf");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_path = directory.path().join("decoded.fastq");
    let exon_1 = pseudo_random_dna(60);
    let exon_2 = pseudo_random_dna(60);
    let reference = [
        exon_1.as_slice(),
        b"GT",
        pseudo_random_dna(100).as_slice(),
        b"AG",
        exon_2.as_slice(),
    ]
    .concat();
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference.as_slice(), b"\n"].concat(),
    )
    .unwrap();
    std::fs::write(
        &gtf_path,
        b"chr1\ttest\texon\t1\t60\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr1\ttest\texon\t165\t224\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n",
    )
    .unwrap();
    let read = [&exon_1[30..], &exon_2[..30]].concat();
    let fastq = format!(
        "@spliced\n{}\n+\n{}\n",
        String::from_utf8_lossy(&read),
        "I".repeat(read.len())
    );
    std::fs::write(&input_path, &fastq).unwrap();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("index")
        .arg(&reference_path)
        .args(["--aligner", "rustar", "--gtf"])
        .arg(&gtf_path)
        .args([
            "--sjdb-overhang",
            "30",
            "--sa-index-nbases",
            "3",
            "--chr-bin-nbits",
            "8",
            "--threads",
            "1",
            "--output",
        ])
        .arg(&index_path)
        .assert()
        .success();
    assert_eq!(&std::fs::read(&index_path).unwrap()[..8], b"SRSTAR01");

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("compress")
        .args(["--index"])
        .arg(&index_path)
        .args(["--mate-1"])
        .arg(&input_path)
        .args(["--threads", "1", "--chunk-records", "1", "--output"])
        .arg(&archive_path)
        .assert()
        .success()
        .stderr(predicate::str::contains("100.00% aligned"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .args(["--mate-1"])
        .arg(&decoded_path)
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(decoded_path).unwrap(), fastq);
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
        .failure()
        .stderr(predicate::str::contains("missing.fastq.gz"));
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
fn validate_respects_record_limit() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let archive_path = directory.path().join("reads.seqarc");

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
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("-1")
        .arg(&input_path)
        .arg("--records")
        .arg("1")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 1 records (record limit 1) against originals (SingleEnd)",
        ));

    // A limit past the end of the archive behaves like a full run.
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("-1")
        .arg(&input_path)
        .arg("--records")
        .arg("5")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 2 records (record limit 5) against originals (SingleEnd)",
        ));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--records")
        .arg("0")
        .assert()
        .failure()
        .stderr(predicate::str::contains("0 is not in 1.."));
}

#[test]
fn compress_post_validates_first_records() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");

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

    let compress = |extra_args: &[&str], archive: &std::path::Path| {
        Command::cargo_bin("seqarc")
            .unwrap()
            .args(["compress", "--index"])
            .arg(&index_path)
            .arg("--mate-1")
            .arg(&input_path)
            .args(extra_args)
            .arg("--output")
            .arg(archive)
            .assert()
            .success()
    };

    // Default: the limit exceeds the archive, so every record is checked.
    let archive_path = directory.path().join("default.seqarc");
    compress(&[], &archive_path).stderr(
        predicate::str::contains("stage=post_validate done")
            .and(predicate::str::contains("records=2")),
    );

    // An explicit limit below the archive size checks only the first records.
    let limited_path = directory.path().join("limited.seqarc");
    compress(&["--post-validate-records", "1"], &limited_path).stderr(
        predicate::str::contains("stage=post_validate done")
            .and(predicate::str::contains("records=1")),
    );

    // Zero disables the check entirely.
    let disabled_path = directory.path().join("disabled.seqarc");
    compress(&["--post-validate-records", "0"], &disabled_path)
        .stderr(predicate::str::contains("post_validate").not());

    // Archives without stored identifiers still validate on the remaining
    // fields.
    let no_ids_path = directory.path().join("no-ids.seqarc");
    compress(&["--no-ids"], &no_ids_path).stderr(
        predicate::str::contains("stage=post_validate done")
            .and(predicate::str::contains("records=2")),
    );
}

#[test]
fn compress_from_stdin_skips_post_validation() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let archive_path = directory.path().join("stdin.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let sequence = &reference[20..60];
    let fastq = format!(
        "@read0\n{}\n+\n{}\n",
        String::from_utf8_lossy(sequence),
        "I".repeat(sequence.len()),
    );

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .assert()
        .success();

    // Pipes and /dev/stdin cannot be re-opened for the post-compression
    // check, so the check is skipped instead of failing the compression.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .args(["-1", "/dev/stdin"])
        .arg("--output")
        .arg(&archive_path)
        .write_stdin(fastq.as_bytes())
        .assert()
        .success()
        .stderr(predicate::str::contains(
            "stage=post_validate skipped input_not_reopenable",
        ));
}

#[test]
fn decompress_validates_output_layout_for_empty_archives() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let archive_path = directory.path().join("empty.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    std::fs::write(&input_path, b"").unwrap();

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
        .assert()
        .success();

    // A single-end archive rejects a mate-2 output even with zero records.
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .args(["-1", directory.path().join("e1.fastq").to_str().unwrap()])
        .args(["-2", directory.path().join("e2.fastq").to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "single-end archive received two output paths",
        ));
    assert!(!directory.path().join("e2.fastq").exists());
}

#[test]
fn info_reports_archive_metadata() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let archive_path = directory.path().join("reads.seqarc");

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
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("info")
        .arg(&archive_path)
        .assert()
        .success()
        .stdout(
            predicate::str::contains("format version: 0.28")
                .and(predicate::str::contains("layout: SingleEnd"))
                .and(predicate::str::contains("records: 2"))
                .and(predicate::str::contains("chunks: 1/1 scanned"))
                .and(predicate::str::contains(
                    "ids: preserved (mode: independent)",
                ))
                .and(predicate::str::contains("quality: preserved"))
                .and(predicate::str::contains("file size: "))
                .and(predicate::str::contains("use `seqarc verify`"))
                // The storage breakdown only appears with --breakdown.
                .and(predicate::str::contains("storage breakdown:").not()),
        );

    // The breakdown reconciles every on-disk section with the file size.
    let output = Command::cargo_bin("seqarc")
        .unwrap()
        .arg("info")
        .arg(&archive_path)
        .arg("--breakdown")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    for expected in [
        "storage breakdown:",
        "preamble",
        "chunk headers",
        "payload directories",
        "identifiers:",
        "qualities:",
        "records:",
        "    flags:",
        "    lengths:",
        "    positions:",
        "    noise positions:",
        "    noise:",
        "    reference skips:",
        "    emits:",
        "    reverses:",
        "    unaligned sequences:",
        "reference block:",
        "index metadata:",
        "  accounted:",
    ] {
        assert!(
            stdout.contains(expected),
            "missing {expected:?} in:\n{stdout}"
        );
    }
    // Single-end subsections carry no mate prefix.
    assert!(!stdout.contains("R1 lengths"));
    let file_size = std::fs::metadata(&archive_path).unwrap().len();
    let (accounted, total) = parse_accounted(&stdout);
    assert_eq!(accounted, file_size);
    assert_eq!(total, file_size);

    // A truncated chunk header is reported without failing the command.
    // The preamble spans 46 + reference(27) + metadata(8) = 81 bytes; keep
    // it plus only part of the 48-byte chunk header.
    let archive = std::fs::read(&archive_path).unwrap();
    let truncated_path = directory.path().join("truncated.seqarc");
    std::fs::write(&truncated_path, &archive[..100]).unwrap();
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("info")
        .arg(&truncated_path)
        .assert()
        .success()
        .stdout(predicate::str::contains("chunks: 0/1 scanned"))
        .stderr(predicate::str::contains("archive is truncated"));

    // The breakdown tolerates truncation the same way, aggregating only the
    // chunks that were scanned.
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("info")
        .arg(&truncated_path)
        .arg("--breakdown")
        .assert()
        .success()
        .stdout(predicate::str::contains("chunks: 0/1 scanned"))
        .stderr(predicate::str::contains("archive is truncated"));

    // A corrupt payload section directory is rejected: the 24-byte directory
    // starts right after the preamble (81 bytes) and chunk header (48 bytes).
    let mut corrupt = archive.clone();
    corrupt[81 + 48] ^= 0xff;
    let corrupt_path = directory.path().join("corrupt.seqarc");
    std::fs::write(&corrupt_path, &corrupt).unwrap();
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("info")
        .arg(&corrupt_path)
        .arg("--breakdown")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "chunk payload section sizes do not match the chunk header",
        ));
}

/// Parse the `accounted  X of Y bytes` line from an `info --breakdown` run.
fn parse_accounted(stdout: &str) -> (u64, u64) {
    let line = stdout
        .lines()
        .find(|line| line.trim_start().starts_with("accounted:"))
        .expect("accounted line missing");
    let body = line.trim().trim_start_matches("accounted:").trim();
    let (accounted, total) = body.split_once(" of ").expect("accounted format");
    let accounted: u64 = accounted.replace(',', "").parse().unwrap();
    let total: u64 = total
        .trim_end_matches(" bytes")
        .replace(',', "")
        .parse()
        .unwrap();
    (accounted, total)
}

#[test]
fn info_breakdown_reports_paired_subsections() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let second_path = directory.path().join("reads.2.fastq");
    let archive_path = directory.path().join("reads.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let first_sequence = &reference[20..60];
    let second_sequence: Vec<u8> = first_sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect();
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
    std::fs::write(&first_path, first_fastq.as_bytes()).unwrap();
    std::fs::write(&second_path, second_fastq.as_bytes()).unwrap();

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
        .assert()
        .success();

    let output = Command::cargo_bin("seqarc")
        .unwrap()
        .arg("info")
        .arg(&archive_path)
        .arg("--breakdown")
        .output()
        .unwrap();
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout).unwrap();
    // Paired archives group the eight subsections per mate.
    for expected in [
        "    R1 lengths:",
        "    R1 positions:",
        "    R1 noise:",
        "    R1 unaligned sequences:",
        "    R2 lengths:",
        "    R2 positions:",
        "    R2 noise:",
        "    R2 unaligned sequences:",
        // The unaligned second read contributes raw sequence bytes.
        "raw",
    ] {
        assert!(
            stdout.contains(expected),
            "missing {expected:?} in:\n{stdout}"
        );
    }
    let file_size = std::fs::metadata(&archive_path).unwrap().len();
    let (accounted, total) = parse_accounted(&stdout);
    assert_eq!(accounted, file_size);
    assert_eq!(total, file_size);
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
fn id_field_selection_round_trips_and_combines_with_prefix() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let second_path = directory.path().join("reads.2.fastq");
    let decoded_first = directory.path().join("decoded.1.fastq");
    let decoded_second = directory.path().join("decoded.2.fastq");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let sequence = &reference[20..60];
    let reverse_complement = sequence
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
        "@read0 1:N:0:ATCG+CCAG\n{}\n+\n{}\n@read1 1:N:0:GGTT+TTAA\n{}\n+\n{}\n",
        String::from_utf8_lossy(sequence),
        "I".repeat(sequence.len()),
        "N".repeat(40),
        "#".repeat(40),
    );
    let second_fastq = format!(
        "@read0 2:N:0:ATCG+CCAG\n{}\n+\n{}\n@read1 2:N:0:GGTT+TTAA\n{}\n+\n{}\n",
        String::from_utf8_lossy(&reverse_complement),
        "J".repeat(sequence.len()),
        "N".repeat(40),
        "$".repeat(40),
    );
    std::fs::write(&first_path, first_fastq.as_bytes()).unwrap();
    std::fs::write(&second_path, second_fastq.as_bytes()).unwrap();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .assert()
        .success();

    let field_only_archive = directory.path().join("field.seqarc");
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--output")
        .arg(&field_only_archive)
        .arg("--chunk-records")
        .arg("1")
        .arg("--id-field")
        .arg("1")
        .assert()
        .success()
        .stderr(predicate::str::contains("id_field=Some(1)"))
        .stderr(predicate::str::contains(
            "id-transform preview: \"@read0 1:N:0:ATCG+CCAG\" -> \"@1:N:0:ATCG+CCAG\"",
        ));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&field_only_archive)
        .arg("-1")
        .arg(&first_path)
        .arg("-2")
        .arg(&second_path)
        .assert()
        .failure();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&field_only_archive)
        .arg("-1")
        .arg(&first_path)
        .arg("-2")
        .arg(&second_path)
        .arg("--id-field")
        .arg("1")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 2 records against originals (PairedEnd)",
        ));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&field_only_archive)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&decoded_first).unwrap(),
        first_fastq
            .replace("@read0 1:N:0:ATCG+CCAG", "@1:N:0:ATCG+CCAG")
            .replace("@read1 1:N:0:GGTT+TTAA", "@1:N:0:GGTT+TTAA")
    );
    assert_eq!(
        std::fs::read_to_string(&decoded_second).unwrap(),
        second_fastq
            .replace("@read0 2:N:0:ATCG+CCAG", "@2:N:0:ATCG+CCAG")
            .replace("@read1 2:N:0:GGTT+TTAA", "@2:N:0:GGTT+TTAA")
    );

    let combined_archive = directory.path().join("combined.seqarc");
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--output")
        .arg(&combined_archive)
        .arg("--id-field")
        .arg("1")
        .arg("--id-prefix-chars")
        .arg("3")
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&combined_archive)
        .arg("-1")
        .arg(&first_path)
        .arg("-2")
        .arg(&second_path)
        .arg("--id-field")
        .arg("1")
        .arg("--id-prefix-chars")
        .arg("3")
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&combined_archive)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .arg("--force")
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&decoded_first).unwrap(),
        first_fastq
            .replace("@read0 1:N:0:ATCG+CCAG", "@1:N")
            .replace("@read1 1:N:0:GGTT+TTAA", "@1:N")
    );
    assert_eq!(
        std::fs::read_to_string(&decoded_second).unwrap(),
        second_fastq
            .replace("@read0 2:N:0:ATCG+CCAG", "@2:N")
            .replace("@read1 2:N:0:GGTT+TTAA", "@2:N")
    );

    let missing_field_archive = directory.path().join("missing-field.seqarc");
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--output")
        .arg(&missing_field_archive)
        .arg("--id-field")
        .arg("2")
        .assert()
        .failure();
    assert!(!missing_field_archive.exists());

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--output")
        .arg(&missing_field_archive)
        .arg("--no-ids")
        .arg("--id-field")
        .arg("0")
        .assert()
        .failure();
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
    assert_eq!(u16::from_le_bytes(archive[9..11].try_into().unwrap()), 0);
    assert_eq!(u16::from_le_bytes(archive[11..13].try_into().unwrap()), 28);

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

#[test]
fn embedded_index_range_may_consume_the_entire_read() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_first = directory.path().join("decoded.1.fastq");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let first_fastq = b"@read0\nAAAACCCCGGGGTTTTACGT\n+\nABCDEFGHIJKLMNOPQRST\n";
    std::fs::write(&first_path, first_fastq).unwrap();

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
        .arg("--index-range")
        .arg("R1:1:20")
        .arg("--chunk-records")
        .arg("1")
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&first_path)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--threads")
        .arg("2")
        .assert()
        .success();

    assert_eq!(std::fs::read(&decoded_first).unwrap(), first_fastq);
}

#[test]
fn paired_with_index_read_files_round_trips() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let second_path = directory.path().join("reads.2.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let i2_path = directory.path().join("index.2.fastq");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_first = directory.path().join("decoded.1.fastq");
    let decoded_second = directory.path().join("decoded.2.fastq");
    let decoded_i1 = directory.path().join("decoded.i1.fastq");
    let decoded_i2 = directory.path().join("decoded.i2.fastq");

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
    let i1_fastq = "@pair0/1\nAACCGGTT\n+\nJJJJJJJJ\n@pair1/1\nTTTTGGCC\n+\nJJJJJJJJ\n";
    let i2_fastq = "@pair0/1\nGGCCTTAA\n+\nKKKKKKKK\n@pair1/1\nACGTACGT\n+\nKKKKKKKK\n";
    std::fs::write(&first_path, first_fastq.as_bytes()).unwrap();
    std::fs::write(&second_path, second_fastq.as_bytes()).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();
    std::fs::write(&i2_path, i2_fastq).unwrap();

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
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--chunk-records")
        .arg("1")
        .arg("--profile")
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .success()
        .stderr(predicate::str::contains("index_streams=2"))
        .stderr(predicate::str::contains("index reads: 0 raw").not())
        .stderr(predicate::str::contains("index reads:"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .arg("--index-output-1")
        .arg(&decoded_i1)
        .arg("--index-output-2")
        .arg(&decoded_i2)
        .arg("--threads")
        .arg("2")
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&decoded_first).unwrap(),
        first_fastq
    );
    assert_eq!(
        std::fs::read_to_string(&decoded_second).unwrap(),
        second_fastq
    );
    assert_eq!(std::fs::read(&decoded_i1).unwrap(), i1_fastq.as_bytes());
    assert_eq!(std::fs::read(&decoded_i2).unwrap(), i2_fastq.as_bytes());

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--index-read-2")
        .arg(&i2_path)
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 2 records against originals (PairedEnd)",
        ));

    let missing_output = directory.path().join("missing-i2.fastq");
    let partial_first = directory.path().join("partial.1.fastq");
    let partial_second = directory.path().join("partial.2.fastq");
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&partial_first)
        .arg("--mate-2")
        .arg(&partial_second)
        .arg("--index-output-1")
        .arg(&missing_output)
        .assert()
        .failure();
    assert!(!partial_first.exists());
    assert!(!partial_second.exists());
    assert!(!missing_output.exists());
}

#[test]
fn index_read_descriptions_may_differ_across_files() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let second_path = directory.path().join("reads.2.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let i2_path = directory.path().join("index.2.fastq");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let sequence = &reference[20..60];
    let reverse_complement = sequence
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
        "@pair0 1:N:0:ATCG+CCAG\n{}\n+\n{}\n@pair1 1:N:0:GGTT+TTAA\n{}\n+\n{}\n",
        String::from_utf8_lossy(sequence),
        "I".repeat(sequence.len()),
        "N".repeat(40),
        "#".repeat(40),
    );
    let second_fastq = format!(
        "@pair0 2:N:0:ATCG+CCAG\n{}\n+\n{}\n@pair1 2:N:0:GGTT+TTAA\n{}\n+\n{}\n",
        String::from_utf8_lossy(&reverse_complement),
        "J".repeat(sequence.len()),
        "N".repeat(40),
        "$".repeat(40),
    );
    let i1_fastq = "@pair0 some other comment\nAACCGGTT\n+\nJJJJJJJJ\n@pair1 another comment\nTTTTGGCC\n+\nJJJJJJJJ\n";
    let i2_fastq = "@pair0 a third comment\nGGCCTTAA\n+\nKKKKKKKK\n@pair1 a fourth comment\nACGTACGT\n+\nKKKKKKKK\n";
    std::fs::write(&first_path, first_fastq.as_bytes()).unwrap();
    std::fs::write(&second_path, second_fastq.as_bytes()).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();
    std::fs::write(&i2_path, i2_fastq).unwrap();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .assert()
        .success();

    // Without a transform the differing descriptions cannot round-trip, since
    // index reads inherit mate-1's identifier during reconstruction.
    let raw_archive = directory.path().join("raw.seqarc");
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--output")
        .arg(&raw_archive)
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "transformed index read identifier does not match transformed mate 1 identifier",
        ));
    assert!(!raw_archive.exists());

    // Field selection homogenizes the identifiers, so the archive round-trips
    // and validate accepts the transformed originals.
    let field_archive = directory.path().join("field.seqarc");
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--id-field")
        .arg("0")
        .arg("--output")
        .arg(&field_archive)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&field_archive)
        .arg("-1")
        .arg(&first_path)
        .arg("-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--id-field")
        .arg("0")
        .assert()
        .success()
        .stdout(predicate::str::contains(
            "validated 2 records against originals (PairedEnd)",
        ));

    let decoded_i1 = directory.path().join("decoded.i1.fastq");
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&field_archive)
        .arg("--mate-1")
        .arg(directory.path().join("field.1.fastq"))
        .arg("--mate-2")
        .arg(directory.path().join("field.2.fastq"))
        .arg("--index-output-1")
        .arg(&decoded_i1)
        .arg("--index-output-2")
        .arg(directory.path().join("decoded.i2.fastq"))
        .arg("--force")
        .assert()
        .success();
    assert_eq!(
        std::fs::read_to_string(&decoded_i1).unwrap(),
        "@pair0\nAACCGGTT\n+\nJJJJJJJJ\n@pair1\nTTTTGGCC\n+\nJJJJJJJJ\n"
    );

    // A differing identifier proper is always rejected, even with field selection.
    let mismatched_i1 = directory.path().join("mismatched.1.fastq");
    let mismatched_i1_fastq = i1_fastq
        .replace("@pair0", "@other0")
        .replace("@pair1", "@other1");
    std::fs::write(&mismatched_i1, mismatched_i1_fastq).unwrap();
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&mismatched_i1)
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--id-field")
        .arg("0")
        .arg("--output")
        .arg(directory.path().join("mismatched.seqarc"))
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "index read identifier does not match mate 1 identifier",
        ));

    // Field selection that diverges between mate 1 and the index reads is
    // rejected after transformation.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--id-field")
        .arg("1")
        .arg("--output")
        .arg(directory.path().join("divergent.seqarc"))
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "transformed index read identifier does not match transformed mate 1 identifier",
        ));
}

#[test]
fn single_end_with_index_read_file_round_trips() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let i2_path = directory.path().join("index.2.fastq");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_path = directory.path().join("decoded.fastq");
    let decoded_i1 = directory.path().join("decoded.i1.fastq");

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
    let i1_fastq = "@read0\nAACCGGTT\n+\nJJJJJJJJ\n@read1\nTTTTGGCC\n+\nJJJJJJJJ\n";
    std::fs::write(&input_path, fastq.as_bytes()).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();

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
        .arg("--index-read-2")
        .arg(&i2_path)
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .failure();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&input_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--chunk-records")
        .arg("1")
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .success()
        .stderr(predicate::str::contains("index_streams=1"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_path)
        .arg("--index-output-1")
        .arg(&decoded_i1)
        .assert()
        .success();
    assert_eq!(std::fs::read_to_string(&decoded_path).unwrap(), fastq);
    assert_eq!(std::fs::read(&decoded_i1).unwrap(), i1_fastq.as_bytes());

    let no_index_output = directory.path().join("no-index.fastq");
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&no_index_output)
        .assert()
        .failure();
    assert!(!no_index_output.exists());

    let excess_1 = directory.path().join("excess.1.fastq");
    let excess_i2 = directory.path().join("excess.i2.fastq");
    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&excess_1)
        .arg("--index-output-1")
        .arg(&excess_i2)
        .arg("--index-output-2")
        .arg(&excess_1)
        .assert()
        .failure();
    assert!(!excess_1.exists());
    assert!(!excess_i2.exists());
}

#[test]
fn index_reads_coexist_with_ranges_and_respect_global_options() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let second_path = directory.path().join("reads.2.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let archive_path = directory.path().join("reads.seqarc");
    let decoded_first = directory.path().join("decoded.1.fastq");
    let decoded_second = directory.path().join("decoded.2.fastq");
    let decoded_i1 = directory.path().join("decoded.i1.fastq");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let first_fastq = b"@pair0/1\nAAAACCCCGGGGTTTTACGT\n+\nABCDEFGHIJKLMNOPQRST\n\
@pair1/1\nAAAACCCCGGAGTTTTACGT\n+\nTSRQPONMLKJIHGFEDCBA\n";
    let second_fastq = b"@pair0/2\nTGCATTTTAAAACCCCGGTA\n+\nIIIIIIIIIIIIIIIIIIII\n\
@pair1/2\nTGCATTTGAAAACCCCGGTA\n+\nJJJJJJJJJJJJJJJJJJJJ\n";
    let i1_fastq = "@pair0/1\nAACCGGT\n+\nJJJJJJJ\n@pair1/1\nTTTGGCC\n+\nJJJJJJJ\n";
    std::fs::write(&first_path, first_fastq).unwrap();
    std::fs::write(&second_path, second_fastq).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();

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
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--chunk-records")
        .arg("1")
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .success()
        .stderr(predicate::str::contains("index_streams=1"));

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .arg("--index-output-1")
        .arg(&decoded_i1)
        .assert()
        .success();
    assert_eq!(std::fs::read(&decoded_first).unwrap(), first_fastq);
    assert_eq!(std::fs::read(&decoded_second).unwrap(), second_fastq);
    assert_eq!(std::fs::read(&decoded_i1).unwrap(), i1_fastq.as_bytes());

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--mate-2")
        .arg(&second_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .assert()
        .success();

    let reduced_archive = directory.path().join("reduced.seqarc");
    let fasta_first = directory.path().join("decoded.1.fa");
    let fasta_second = directory.path().join("decoded.2.fa");
    let fasta_i1 = directory.path().join("decoded.i1.fa");
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
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--chunk-records")
        .arg("1")
        .arg("--no-ids")
        .arg("--quality-mode")
        .arg("illumina")
        .arg("--output")
        .arg(&reduced_archive)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("decompress")
        .arg(&reduced_archive)
        .arg("--mate-1")
        .arg(&fasta_first)
        .arg("--mate-2")
        .arg(&fasta_second)
        .arg("--index-output-1")
        .arg(&fasta_i1)
        .arg("--format")
        .arg("fasta")
        .assert()
        .success();

    let expected_first_fasta = String::from_utf8_lossy(first_fastq)
        .lines()
        .enumerate()
        .filter(|(index, _)| index % 4 == 1)
        .map(|(_, line)| format!(">\n{line}\n"))
        .collect::<String>();
    let expected_second_fasta = String::from_utf8_lossy(second_fastq)
        .lines()
        .enumerate()
        .filter(|(index, _)| index % 4 == 1)
        .map(|(_, line)| format!(">\n{line}\n"))
        .collect::<String>();
    let expected_i1_fasta = i1_fastq
        .lines()
        .enumerate()
        .filter(|(index, _)| index % 4 == 1)
        .map(|(_, line)| format!(">\n{line}\n"))
        .collect::<String>();
    assert_eq!(
        std::fs::read_to_string(&fasta_first).unwrap(),
        expected_first_fasta
    );
    assert_eq!(
        std::fs::read_to_string(&fasta_second).unwrap(),
        expected_second_fasta
    );
    assert_eq!(
        std::fs::read_to_string(&fasta_i1).unwrap(),
        expected_i1_fasta
    );
}

#[test]
fn index_read_files_reject_unequal_record_counts() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let archive_path = directory.path().join("reads.seqarc");

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
    let i1_fastq = "@read0\nAACCGGTT\n+\nJJJJJJJJ\n";
    std::fs::write(&input_path, fastq.as_bytes()).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();

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
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .failure();
    assert!(!archive_path.exists());
}

#[test]
fn index_read_files_reject_mismatched_identifiers() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let archive_path = directory.path().join("reads.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let sequence = &reference[20..60];
    let fastq = format!(
        "@read0\n{}\n+\n{}\n",
        String::from_utf8_lossy(sequence),
        "I".repeat(sequence.len()),
    );
    let i1_fastq = "@other\nAACCGGTT\n+\nJJJJJJJJ\n";
    std::fs::write(&input_path, fastq.as_bytes()).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();

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
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "index read identifier does not match mate 1 identifier",
        ));
    assert!(!archive_path.exists());
}

#[test]
fn validate_rejects_mismatched_index_originals() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let input_path = directory.path().join("reads.fastq");
    let i1_path = directory.path().join("index.1.fastq");
    let wrong_i1_path = directory.path().join("wrong-index.1.fastq");
    let archive_path = directory.path().join("reads.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    let sequence = &reference[20..60];
    let fastq = format!(
        "@read0\n{}\n+\n{}\n",
        String::from_utf8_lossy(sequence),
        "I".repeat(sequence.len()),
    );
    let i1_fastq = "@read0\nAACCGGTT\n+\nJJJJJJJJ\n";
    let wrong_i1_fastq = "@wrong\nTTTTGGCC\n+\nJJJJJJJJ\n";
    std::fs::write(&input_path, fastq.as_bytes()).unwrap();
    std::fs::write(&i1_path, i1_fastq).unwrap();
    std::fs::write(&wrong_i1_path, wrong_i1_fastq).unwrap();

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
        .arg("--index-read-1")
        .arg(&i1_path)
        .arg("--chunk-records")
        .arg("1")
        .arg("--output")
        .arg(&archive_path)
        .assert()
        .success();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&input_path)
        .arg("--index-read-1")
        .arg(&wrong_i1_path)
        .assert()
        .failure();

    Command::cargo_bin("seqarc")
        .unwrap()
        .arg("validate")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&input_path)
        .arg("--index-read-1")
        .arg(&i1_path)
        .assert()
        .success();
}

fn qvz_test_quality(length: usize, seed: u64) -> String {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (35 + (state >> 33) % 38) as u8 as char
        })
        .collect()
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => unreachable!(),
        })
        .collect()
}

#[test]
fn qvz_quality_mode_round_trips_with_lossy_quality() {
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

    // Six pairs aligned to the reference with varied quality strings so the
    // QVZ training statistics are non-degenerate.
    let read_length = 40usize;
    let mut first_fastq = String::new();
    let mut second_fastq = String::new();
    for pair in 0..6usize {
        let start = 8 + pair * 5;
        let first_sequence = &reference[start..start + read_length];
        let second_sequence = reverse_complement(&reference[start + 2..start + 2 + read_length]);
        first_fastq.push_str(&format!(
            "@pair{pair}/1\n{}\n+\n{}\n",
            String::from_utf8_lossy(first_sequence),
            qvz_test_quality(read_length, 1000 + pair as u64)
        ));
        second_fastq.push_str(&format!(
            "@pair{pair}/2\n{}\n+\n{}\n",
            String::from_utf8_lossy(&second_sequence),
            qvz_test_quality(read_length, 2000 + pair as u64)
        ));
    }
    std::fs::write(&first_path, &first_fastq).unwrap();
    std::fs::write(&second_path, &second_fastq).unwrap();

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
        .arg("--quality-mode")
        .arg("qvz")
        .arg("--qvz-ratio")
        .arg("1.0")
        .assert()
        .success();

    // Quality is lossy, so full validation against the originals must fail.
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
        .arg("decompress")
        .arg(&archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .assert()
        .success();

    let parse_fastq = |path: &std::path::Path| -> Vec<(String, String, String)> {
        let text = std::fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        (0..lines.len())
            .step_by(4)
            .map(|index| {
                (
                    lines[index].to_owned(),
                    lines[index + 1].to_owned(),
                    lines[index + 3].to_owned(),
                )
            })
            .collect()
    };

    let original_first = parse_fastq(&first_path);
    let original_second = parse_fastq(&second_path);
    let decoded_mate1 = parse_fastq(&decoded_first);
    let decoded_mate2 = parse_fastq(&decoded_second);

    assert_eq!(decoded_mate1.len(), original_first.len());
    assert_eq!(decoded_mate2.len(), original_second.len());
    let mut quality_changed = false;
    for (((original, decoded), original_second), decoded_second) in original_first
        .iter()
        .zip(&decoded_mate1)
        .zip(&original_second)
        .zip(&decoded_mate2)
    {
        // Identifiers and sequences survive exactly.
        assert_eq!(decoded.0, original.0);
        assert_eq!(decoded.1, original.1);
        assert_eq!(decoded_second.0, original_second.0);
        assert_eq!(decoded_second.1, original_second.1);
        // Qualities keep their length but are quantized.
        assert_eq!(decoded.2.len(), original.2.len());
        assert_eq!(decoded_second.2.len(), original_second.2.len());
        assert!(decoded.2.chars().all(|q| ('!'..='h').contains(&q)));
        assert!(decoded_second.2.chars().all(|q| ('!'..='h').contains(&q)));
        if decoded.2 != original.2 || decoded_second.2 != original_second.2 {
            quality_changed = true;
        }
    }
    assert!(quality_changed, "qvz mode must alter quality values");
}

#[test]
fn qvz_quality_mode_rejects_invalid_ratio_configuration() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq");
    let archive_path = directory.path().join("reads.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();
    std::fs::write(
        &first_path,
        "@r1\nACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGG\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n",
    )
    .unwrap();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .assert()
        .success();

    // Missing ratio.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--output")
        .arg(&archive_path)
        .arg("--quality-mode")
        .arg("qvz")
        .assert()
        .failure()
        .stderr(predicate::str::contains("requires --qvz-ratio"));

    // Zero ratio.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--output")
        .arg(&archive_path)
        .arg("--quality-mode")
        .arg("qvz")
        .arg("--qvz-ratio")
        .arg("0")
        .assert()
        .failure()
        .stderr(predicate::str::contains("positive"));

    // Ratio supplied without the qvz mode.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--output")
        .arg(&archive_path)
        .arg("--qvz-ratio")
        .arg("1.0")
        .assert()
        .failure()
        .stderr(predicate::str::contains(
            "--qvz-ratio requires --quality-mode qvz",
        ));

    // Incompatible with --no-quality.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["compress", "--index"])
        .arg(&index_path)
        .arg("--mate-1")
        .arg(&first_path)
        .arg("--output")
        .arg(&archive_path)
        .arg("--no-quality")
        .arg("--quality-mode")
        .arg("qvz")
        .arg("--qvz-ratio")
        .arg("1.0")
        .assert()
        .failure()
        .stderr(predicate::str::contains("--no-quality"));

    assert!(!archive_path.exists());
}

#[test]
fn compress_rejects_quality_bytes_outside_printable_ascii() {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let archive_path = directory.path().join("reads.seqarc");

    let reference =
        b"ACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGGTACAGTCTAGGCTAACGTGCTTAGCCTGATACGTTCAAGGTC";
    std::fs::write(
        &reference_path,
        [b">chr1\n".as_slice(), reference, b"\n"].concat(),
    )
    .unwrap();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["index", reference_path.to_str().unwrap(), "--output"])
        .arg(&index_path)
        .assert()
        .success();

    for bad_byte in [0x01u8, 0x7f] {
        let reads_path = directory.path().join("reads.fastq");
        let mut quality = vec![b'I'; 40];
        quality[5] = bad_byte;
        let mut record = b"@r1\nACGTTGCACTGATCGTACCTGAGTTCGATGCCATGACTGG\n+\n".to_vec();
        record.extend_from_slice(&quality);
        record.push(b'\n');
        std::fs::write(&reads_path, record).unwrap();

        Command::cargo_bin("seqarc")
            .unwrap()
            .args(["compress", "--index"])
            .arg(&index_path)
            .arg("--mate-1")
            .arg(&reads_path)
            .arg("--output")
            .arg(&archive_path)
            .assert()
            .failure()
            .stderr(predicate::str::contains(
                "quality scores must be printable Phred+33 ASCII (bytes 33..=126)",
            ));
    }

    assert!(!archive_path.exists());
}

struct PairedArchiveFixture {
    _directory: tempfile::TempDir,
    archive_path: std::path::PathBuf,
    first_fastq: String,
    second_fastq: String,
}

/// Build a two-pair, two-chunk paired archive: pair0 maps to the reference and
/// pair1 is an unmapped N-read fallback, with one record per chunk.
fn create_paired_archive_fixture() -> PairedArchiveFixture {
    let directory = tempfile::tempdir().unwrap();
    let reference_path = directory.path().join("reference.fa");
    let index_path = directory.path().join("reference.seqarc-index");
    let first_path = directory.path().join("reads.1.fastq.gz");
    let second_path = directory.path().join("reads.2.fastq.gz");
    let archive_path = directory.path().join("reads.seqarc");

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
        .args(["--threads", "1", "--chunk-records", "1"])
        .assert()
        .success();

    PairedArchiveFixture {
        _directory: directory,
        archive_path,
        first_fastq,
        second_fastq,
    }
}

fn flip_final_archive_byte(archive_path: &std::path::Path) {
    let mut bytes = std::fs::read(archive_path).unwrap();
    let last = bytes.len() - 1;
    bytes[last] ^= 0xff;
    std::fs::write(archive_path, bytes).unwrap();
}

#[test]
fn verify_detects_corrupted_chunks() {
    let fixture = create_paired_archive_fixture();

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["verify"])
        .arg(&fixture.archive_path)
        .assert()
        .success()
        .stdout(predicate::str::contains("all checksums match"));

    flip_final_archive_byte(&fixture.archive_path);

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["verify"])
        .arg(&fixture.archive_path)
        .assert()
        .failure()
        .stderr(predicate::str::contains("chunks failed their checksums"));
}

#[test]
fn decompress_range_selects_logical_records() {
    let fixture = create_paired_archive_fixture();
    let decoded_first = fixture._directory.path().join("ranged.1.fastq");
    let decoded_second = fixture._directory.path().join("ranged.2.fastq");

    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["decompress"])
        .arg(&fixture.archive_path)
        .arg("--mate-1")
        .arg(&decoded_first)
        .arg("--mate-2")
        .arg(&decoded_second)
        .args(["--range", "1:1"])
        .assert()
        .success();

    let expected_first: Vec<&str> = fixture.first_fastq.lines().skip(4).take(4).collect();
    let expected_second: Vec<&str> = fixture.second_fastq.lines().skip(4).take(4).collect();
    assert_eq!(
        std::fs::read_to_string(&decoded_first).unwrap(),
        format!("{}\n", expected_first.join("\n"))
    );
    assert_eq!(
        std::fs::read_to_string(&decoded_second).unwrap(),
        format!("{}\n", expected_second.join("\n"))
    );
}

#[test]
fn best_effort_decompress_recovers_undamaged_chunks() {
    let fixture = create_paired_archive_fixture();
    let strict_first = fixture._directory.path().join("strict.1.fastq");
    let best_effort_first = fixture._directory.path().join("recovered.1.fastq");
    let best_effort_second = fixture._directory.path().join("recovered.2.fastq");

    flip_final_archive_byte(&fixture.archive_path);

    // Strict decoding refuses the damaged archive.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["decompress"])
        .arg(&fixture.archive_path)
        .arg("--mate-1")
        .arg(&strict_first)
        .assert()
        .failure();

    // Best-effort decoding recovers pair0 from the undamaged chunk. The
    // command reports the skipped chunk and exits with a checksum failure.
    Command::cargo_bin("seqarc")
        .unwrap()
        .args(["decompress"])
        .arg(&fixture.archive_path)
        .arg("--mate-1")
        .arg(&best_effort_first)
        .arg("--mate-2")
        .arg(&best_effort_second)
        .arg("--best-effort")
        .assert()
        .failure()
        .stderr(predicate::str::contains("best-effort recovery"));

    let recovered_first: Vec<&str> = fixture.first_fastq.lines().take(4).collect();
    let recovered_second: Vec<&str> = fixture.second_fastq.lines().take(4).collect();
    assert_eq!(
        std::fs::read_to_string(&best_effort_first).unwrap(),
        format!("{}\n", recovered_first.join("\n"))
    );
    assert_eq!(
        std::fs::read_to_string(&best_effort_second).unwrap(),
        format!("{}\n", recovered_second.join("\n"))
    );
}
