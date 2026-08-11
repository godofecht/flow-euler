# Architecture

## Repository layout

```
problems/           One .flow file per Project Euler problem
problems/native/    C/C++ helpers for problems needing NTT, i128, GMP, etc.
data/               Input files (grids, digit strings, tables)
bench/              Hand-written C twins + timing harness for benchmarks
scripts/            run.sh, verify.sh, verify-mlir.sh, bench.sh
answers.txt         Expected answers, one per line: "NNN answer"
.build/             Generated C, MLIR, LLVM IR, and binaries (gitignored)
.github/workflows/  CI: verify.yml (C backend), mlir.yml (MLIR backend)
```

## Solution structure

Each problem has a Flow file `problems/pNNN.flow` that prints a single
answer line. The file is self-contained: it defines `main()` and any
helper functions it needs.

### Pure Flow solutions

Most problems (700+) are pure Flow. The transpiler lowers them to C,
clang compiles, and the binary runs.

```
pNNN.flow -> Flow transpiler (--c) -> pNNN.c -> clang -O3 -> pNNN -> answer
```

### Native-backed solutions

Problems that need algorithms difficult to express in Flow (NTT, i128,
GMP, MPFR) have a companion C or C++ file under `problems/native/`.
The Flow file declares the native function as `extern` and calls it
from `main`.

```
pNNN.flow -> Flow transpiler (--c) -> pNNN_flow.c
pNNN.c                              (hand-written native solver)
                                    -> clang -O3 pNNN_flow.c pNNN.c -> pNNN -> answer
```

The Flow file typically looks like:

```flow
extern {
    function pNNN_native() -> i64
}

function main() -> i32 {
    printf("%lld\n", pNNN_native())
    return 0
}
```

For floating-point answers, the native function returns `f64` and
`main` prints with an appropriate format string.

For answers that exceed `i64` range (e.g. p974), the native code prints
the answer directly to stdout and returns a sentinel.

## Compilation backends

Flow supports two native compilation paths:

### C backend (default)

`python3 -m flow.transpiler input.flow --c --lenient -o output.c`

Generates C99 code. Compiled with `clang -O3 -march=native`. This is
the primary path used by `scripts/run.sh` and `scripts/verify.sh`.

### MLIR backend

`python3 -m flow.transpiler input.flow --mlir --lenient -o output.mlir`

Generates MLIR (LLVM dialect). Lowered to LLVM IR via `mlir-opt` and
`mlir-translate`, then compiled with `llc` and `clang`. This path is
exercised by `scripts/verify-mlir.sh` and the `mlir.yml` workflow.

```
pNNN.flow -> Flow transpiler (--mlir) -> pNNN.mlir
          -> mlir-opt (lower to LLVM dialect) -> pNNN_lowered.mlir
          -> mlir-translate (--mlir-to-llvmir) -> pNNN.ll
          -> llc -O3 -> pNNN.s
          -> clang -> pNNN -> answer
```

The MLIR backend currently has a known limitation: `i64` loop variables
trigger a type mismatch in `arith.index_cast` (the generator hardcodes
`i32`). Problems that use `i64` for loop bounds or range iteration fail
MLIR lowering. Problems with only `i32` arithmetic or closed-form
solutions (no loops) pass. See `scripts/verify-mlir.sh` for the test
harness, which reports `mlir_fail` separately from answer mismatches.

## Scripts

### `scripts/run.sh N`

Compile and run a single problem. Auto-detects native helpers and
links them. Prints `pNNN: answer`.

### `scripts/verify.sh [first] [last]`

Run every problem with a known answer in `answers.txt` through the C
backend. Compares output against expected. Exits non-zero on any
mismatch.

### `scripts/verify-mlir.sh [first] [last]`

Same as `verify.sh` but routes through the MLIR backend. Skips
problems with native helpers (the MLIR path does not link external C
TUs). Reports `mlir_fail` for problems that fail MLIR lowering or
compilation, separately from answer mismatches.

### `scripts/bench.sh`

Benchmarks Flow-generated C against hand-written C for PE 010 (sieve
of Eratosthenes). Same algorithm, same `-O3 -march=native` flags,
200 in-process rounds.
