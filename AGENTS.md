# AGENTS.md

## Build & Run

```bash
cargo build
cargo run -- <fasta_file>
```

## Test

```bash
cargo test
```

- Rust edition 2024 — requires Rust 1.85+ (stable since 2025-02-20)

## Dependencies

Only `clap` (4.6.1, derive feature) for CLI parsing.

## Architecture

Single-crate binary (`src/main.rs`). Entrypoint parses CLI args (via `src/args.rs`), reads FASTA (`src/fasta_reader.rs`), and finds sequence overlaps (`src/find_overlap.rs`).

Modules:
- `args` — clap-based CLI (`-i`/`-a` for k-mer range, `-k` for single k, `-o` for outdir)
- `iupac` — IUPAC nucleotide code validation
- `seq_record` — `SeqRecord` struct (name + seq)
- `fasta_reader` — FASTA iterator with IUPAC validation
- `find_overlap` — three overlap functions: `find_overlap_s2s`, `find_overlap_e2s`, `find_overlap_e2e`
- `revcompl` — reverse complement (module currently commented out in `main.rs`, tests exist)
- `lib.rs` — empty placeholder

`main.rs` still has hardcoded example calls; CLI arg processing is wired up but the overlap functions are not yet called on real parsed data.

## Tests

- `find_overlap.rs` has 3 test modules (`tests_s2s`, `tests_e2s`, `tests_e2e`)
- `iupac.rs` has 1 test module (`tests_iupac_validator`) — 8 tests
- `fasta_reader.rs` has 1 test module (`tests_fasta_reader`) — 14 tests
  - Test fasta files live in `test_data/fasta_reader/*.fasta`
- `revcompl.rs` has 2 test modules (`test_revcompl`, `test_make_compl_base`) — not compiled because `mod revcompl` is commented out

## CLI

```
combinator_fq2 0.1.0
By Maksim Sikolenko

A program to find adjacent contigs by matching their ends of length k.

USAGE:
    combinator_fq2 <INPUT> [-i <mink>] [-a <maxk>] [-k <k-mer>] [-o <outdir>]
```
