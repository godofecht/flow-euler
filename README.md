# Flow x Project Euler

[![Verify (C backend)](https://github.com/godofecht/flow-euler/actions/workflows/verify.yml/badge.svg)](https://github.com/godofecht/flow-euler/actions/workflows/verify.yml)
[![MLIR backend](https://github.com/godofecht/flow-euler/actions/workflows/mlir.yml/badge.svg)](https://github.com/godofecht/flow-euler/actions/workflows/mlir.yml)

Solving [Project Euler](https://projecteuler.net) in [Flow](https://github.com/flooooooooooow/flow), a statically typed language that reads close to high-level math and compiles to native C and MLIR.

1007 problems solved. All 1007 known answers in `answers.txt` have Flow solutions. 262 use native C/C++ helpers for algorithms that need `__int128`, GMP, MPFR, or NTT.

## Quick start

```bash
git clone https://github.com/flooooooooooow/flow.git ~/flow
export FLOW_REPO=$HOME/flow

./scripts/run.sh 1              # run problem 001
./scripts/run.sh all            # run every solved problem
./scripts/verify.sh             # check all against answers.txt
./scripts/bench.sh              # Flow-generated C vs hand-written C
```

Each `problems/pNNN.flow` prints a single answer line.

> Point `FLOW_BIN` at `$FLOW_REPO/flow` (not a symlink). Symlinked `~/.local/bin/flow` breaks Flow's `SCRIPT_DIR` resolution.

## MLIR backend

Flow also compiles through MLIR. This repo tests both paths in CI.

```bash
brew install llvm               # provides mlir-opt, mlir-translate, llc
export LLVM_BIN=/opt/homebrew/opt/llvm/bin

./scripts/verify-mlir.sh        # verify through MLIR backend
./scripts/verify-mlir.sh 1 50   # problems 1-50 only
```

See [docs/backends.md](docs/backends.md) for details on the MLIR pipeline, current limitations, and manual compilation steps.

## Progress

| Range | Solved | Status |
|-------|--------|--------|
| 001-200 | 200 | done |
| 201-400 | 200 | done |
| 401-600 | 200 | done |
| 601-1007 | 407 | in progress |

1007 known answers in `answers.txt`. All 1007 have Flow solutions.

262 problems use native C/C++ helpers for NTT, `__int128`, GMP, or MPFR.

## Why Flow here

| Concern | Flow | Hand-written C | Python |
|---------|------|----------------|--------|
| Expressiveness | ranges, typed helpers, clear control flow | verbose loops | excellent, but slow |
| Performance | compiles to `-O3` C / MLIR | baseline | 10-100x slower on tight loops |
| Types | `i32` / `i64` where it matters | same | duck typing |
| Deployment | one native binary | same | interpreter + deps |

### Benchmark (PE 010 sieve, 2e6)

Same algorithm, both compiled with `clang -O3 -march=native`, 200 in-process rounds on Apple Silicon.

| Impl | us / run |
|------|----------|
| Flow to C | **1675** |
| Hand-written C | 1767 |

Reproduce with `./scripts/bench.sh`.

## Layout

```
problems/           one .flow file per Euler problem
problems/native/    C/C++ helpers (NTT, i128, GMP, MPFR)
data/               input grids, digit strings, tables
bench/              hand-written C twins + timing harness
scripts/            run, verify, verify-mlir, bench
answers.txt         expected answers for automated checks
docs/               architecture, backends, adding solutions
```

## Example: PE 001

```flow
function sum_multiples(limit: i64, step: i64) -> i64 {
    let mut total: i64 = 0
    for n in 0..limit step step {
        total = total + n
    }
    return total
}

function solve(limit: i64) -> i64 {
    return sum_multiples(limit, 3) + sum_multiples(limit, 5) - sum_multiples(limit, 15)
}
```

Readable arithmetic-series code. Flow lowers the stepped `for` to a tight C loop.

## Example: native helper (PE 775)

Some algorithms need `__int128` or modular arithmetic that is cleaner in C. The Flow file declares an extern and calls it:

```flow
extern {
    function p775_native() -> i64
}

function main() -> i32 {
    printf("%lld\n", p775_native())
    return 0
}
```

The C file under `problems/native/p775.c` implements the solver and exports `long long p775_native(void)`. The build system auto-links GMP and MPFR if Homebrew headers are present.

See [docs/adding-solutions.md](docs/adding-solutions.md) for the full guide.

## CI

Three GitHub Actions workflows run on every push and pull request:

- **verify.yml**: compiles every solution through the C backend and checks answers against `answers.txt`. Problems listed in `.ci-skip.txt` are skipped (too slow, crash, or broken). Runs on `macos-14` with a 120-minute timeout.
- **mlir.yml**: compiles every pure-Flow solution through the MLIR backend (`Flow -> MLIR -> LLVM IR -> llc -> clang`). Reports MLIR lowering failures separately from answer mismatches. Runs on `macos-14` with a 60-minute timeout.
- **site.yml**: generates a static GitHub Pages site showing the Flow source, generated C, generated MLIR, and output for every problem. Deploys to GitHub Pages on push to `main`.

## Documentation

- [docs/architecture.md](docs/architecture.md): repository structure, solution patterns, compilation paths
- [docs/backends.md](docs/backends.md): C and MLIR backend details, prerequisites, limitations
- [docs/adding-solutions.md](docs/adding-solutions.md): how to add a new solution, native helper conventions

## License

Solutions are educational. Project Euler problems remain (c) Project Euler.
