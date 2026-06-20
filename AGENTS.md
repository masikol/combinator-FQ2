# AGENTS.md

## Build & Run

```bash
cargo build
cargo run
```

## Test

```bash
cargo test
```

- Rust edition 2024 — requires Rust 1.85+ (stable since 2025-02-20)

## Architecture

Single-crate binary (`src/main.rs`). One module: `src/find_overlap.rs` with
three public functions: `find_overlap_s2s`, `find_overlap_e2s`, `find_overlap_e2e`.

`src/lib.rs` is empty (placeholder for future library API).

## Notes

- No tests, no CI, no dependencies beyond std
- The binary prints debug info (k-mer length, slices) to stdout during overlap search
