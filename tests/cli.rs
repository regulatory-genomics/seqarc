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

#[cfg(feature = "rustar")]
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

#[cfg(feature = "rustar")]
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
