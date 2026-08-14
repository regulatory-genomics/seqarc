use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use spring_rs::{canonical_minimizers, PairStore, ReadRecord, ReorderOptions};
use std::hint::black_box;
use std::time::Duration;

const PAIRS: usize = 10_000;
const READ_LEN: usize = 150;

fn record(index: usize, mate: usize) -> ReadRecord {
    let mut sequence = vec![b'A'; READ_LEN];
    for (position, base) in sequence.iter_mut().enumerate() {
        *base = b"ACGT"[(position + index + mate) % 4];
    }
    sequence[index % READ_LEN] = b"ACGT"[(index / 7 + mate) % 4];
    ReadRecord::new(
        format!("read-{index}/{mate}"),
        sequence,
        Some(vec![b'I' + (index % 3) as u8; READ_LEN]),
    )
}

fn build_store(count: usize) -> (tempfile::TempDir, PairStore) {
    let temp = tempfile::tempdir().unwrap();
    let mut store = PairStore::create(temp.path()).unwrap();
    for index in 0..count {
        store.append(&record(index, 1), &record(index, 2)).unwrap();
    }
    (temp, store)
}

fn record_store(c: &mut Criterion) {
    let (_temp, mut store) = build_store(PAIRS);
    let mut group = c.benchmark_group("record_store");
    group.throughput(Throughput::Elements(PAIRS as u64));
    group.bench_function("read_first_full_pair_overread", |b| {
        b.iter(|| {
            store.reset_stats();
            for id in 0..PAIRS as u64 {
                black_box(store.read_first(id).unwrap());
            }
            black_box(store.stats())
        })
    });
    group.bench_function("read_mate1_sequence_only", |b| {
        b.iter(|| {
            store.reset_stats();
            for id in 0..PAIRS as u64 {
                black_box(store.read_sequence(id, 0).unwrap());
            }
            black_box(store.stats())
        })
    });
    group.finish();
}

fn minimizers(c: &mut Criterion) {
    let sequence = record(17, 1).sequence;
    let mut group = c.benchmark_group("minimizers");
    group.throughput(Throughput::Bytes(sequence.len() as u64));
    for (k, window) in [(15, 20), (21, 31), (31, 41)] {
        group.bench_with_input(
            BenchmarkId::new("canonical", format!("k{k}_w{window}")),
            &(k, window),
            |b, &(k, window)| {
                b.iter(|| black_box(canonical_minimizers(black_box(&sequence), k, window)))
            },
        );
    }
    group.finish();
}

fn external_sort(c: &mut Criterion) {
    let mut group = c.benchmark_group("external_sort");
    group.sample_size(10);
    for chunk_size in [10_000usize, 250_000, 2_000_000] {
        group.bench_with_input(
            BenchmarkId::new("primary_assignments", chunk_size),
            &chunk_size,
            |b, &chunk_size| {
                b.iter_batched(
                    || build_store(5_000),
                    |(temp, mut store)| {
                        let options = ReorderOptions {
                            kmer_size: 21,
                            window_size: 31,
                            sort_chunk_size: chunk_size,
                            sort_threads: 4,
                            ..ReorderOptions::default()
                        };
                        black_box(
                            spring_rs::reorder::sort_store_minimizers(
                                &mut store,
                                &options,
                                temp.path(),
                            )
                            .unwrap(),
                        );
                    },
                    criterion::BatchSize::LargeInput,
                )
            },
        );
    }
    group.finish();
}

fn zstd_levels(c: &mut Criterion) {
    let mut payload = Vec::with_capacity(32 * 1024 * 1024);
    for index in 0..200_000 {
        payload.extend(record(index, 1).sequence);
        payload.extend(record(index, 1).quality.unwrap());
    }
    let mut group = c.benchmark_group("zstd");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(20));
    group.throughput(Throughput::Bytes(payload.len() as u64));
    for level in [3, 9, 19] {
        group.bench_with_input(BenchmarkId::new("bulk", level), &level, |b, &level| {
            b.iter(|| black_box(zstd::bulk::compress(black_box(&payload), level).unwrap()))
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    record_store,
    minimizers,
    external_sort,
    zstd_levels
);
criterion_main!(benches);
