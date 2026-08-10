# Flow × Project Euler

[![Verify](https://github.com/godofecht/flow-euler/actions/workflows/verify.yml/badge.svg)](https://github.com/godofecht/flow-euler/actions/workflows/verify.yml)

Solving [Project Euler](https://projecteuler.net) in [Flow](https://github.com/flooooooooooow/flow) — a statically typed language that reads close to high-level math code and compiles to native C (and MLIR).

**Thesis:** Flow’s syntax is expressive enough to write number-theory algorithms clearly; its compilation path makes those programs as fast as (or faster than) hand-written C.

Reference solutions inspired by [leonlan/projecteuler](https://github.com/leonlan/projecteuler).

## Quick start

```bash
# Needs a Flow checkout (https://github.com/flooooooooooow/flow)
export FLOW_REPO=$HOME/flow          # path to the Flow repo
./scripts/run.sh 1                   # run problem 001
./scripts/run.sh all                 # run every solved problem
./scripts/verify.sh                  # check against answers.txt
./scripts/bench.sh                   # Flow-generated C vs hand-written C
```

Each `problems/pNNN.flow` prints a single answer line.

> Tip: point `FLOW_BIN` at `$FLOW_REPO/flow` (not a symlink). Symlinked `~/.local/bin/flow` breaks Flow’s `SCRIPT_DIR` resolution.

## Why Flow here

| Concern | Flow | Typical C | Typical Python |
|---------|------|-----------|----------------|
| Expressiveness | ranges, typed helpers, clear control flow | verbose loops | excellent, but slow |
| Performance | compiles to `-O3` C / MLIR | baseline | 10–100× slower on tight loops |
| Types | `i32` / `i64` where it matters | same | duck typing |
| Deployment | one native binary | same | interpreter + deps |

### Benchmark snapshot (PE 010 sieve, 2e6)

Same algorithm, both compiled with `clang -O3 -march=native`, 200 in-process rounds on Apple Silicon:

| Impl | µs / run |
|------|----------|
| Flow → C | **1675** |
| Hand-written C | 1767 |

Reproduce with `./scripts/bench.sh`.

## Progress

766 problems solved. Of the 750 problems with known answers in `answers.txt`, 728 have Flow solutions (22 remain unsolved, mostly requiring floating-point, bignum, or file input).

| Range | Solved | Status |
|-------|--------|--------|
| 001-200 | 200 | done |
| 201-400 | 200 | done |
| 401-600 | 200 | done |
| 601-1000 | 166 | in progress |

61 problems use native C/C++ helpers (`problems/native/`) for NTT, i128, or multiprecision arithmetic.

## Layout

```
problems/        # one Flow file per Euler problem
problems/native/ # C/C++ helpers for problems needing NTT, i128, or multiprecision
data/            # input grids / digit strings when needed
bench/           # C twins + timing harness
scripts/         # run / verify / bench
answers.txt      # expected answers for automated checks
```

## Example — PE 001

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

Readable arithmetic-series code; Flow lowers the stepped `for` to a tight C loop.

## License

Solutions are educational. Project Euler problems remain © Project Euler.
