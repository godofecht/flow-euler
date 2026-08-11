# Backends

Flow compiles to two native backends. This repo tests both.

## C backend

The default path. Flow transpiles to C99, then `clang -O3 -march=native`
compiles to a native binary.

```
pNNN.flow --c--> pNNN.c --clang -O3--> pNNN --> answer
```

Used by `scripts/run.sh` and `scripts/verify.sh`. Tested in the
`verify.yml` GitHub Actions workflow.

All 855 solutions pass through this path.

## MLIR backend

Flow generates MLIR in the LLVM dialect. The MLIR is lowered to LLVM IR
by `mlir-opt`, translated by `mlir-translate`, then compiled by `llc`
and linked by `clang`.

```
pNNN.flow --mlir--> pNNN.mlir
  --mlir-opt--> pNNN_lowered.mlir
  --mlir-translate--> pNNN.ll
  --llc -O3--> pNNN.s
  --clang--> pNNN --> answer
```

Used by `scripts/verify-mlir.sh`. Tested in the `mlir.yml` GitHub
Actions workflow.

### Prerequisites

The MLIR backend requires LLVM tools not shipped with macOS:

```bash
brew install llvm
export LLVM_BIN=/opt/homebrew/opt/llvm/bin
```

This provides `mlir-opt`, `mlir-translate`, and `llc`.

### Current status

Problems without loops (closed-form solutions, recursive functions)
generally pass through MLIR. Problems using `i64` loop variables or
stepped `for` ranges fail MLIR lowering due to a type mismatch in the
Flow MLIR generator: `arith.index_cast` hardcodes `i32` even when the
operand is `i64`.

The `verify-mlir.sh` script reports three categories:

- `passed`: MLIR compilation succeeded and the answer matches
- `mlir_fail`: MLIR lowering or compilation failed (compiler limitation)
- `failed`: MLIR compilation succeeded but the answer was wrong

Only `failed` (answer mismatch) causes a non-zero exit. `mlir_fail` is
expected for a significant fraction of problems until the Flow compiler
fixes the `index_cast` type bug.

### Running locally

```bash
export FLOW_REPO=$HOME/flow
export LLVM_BIN=/opt/homebrew/opt/llvm/bin
./scripts/verify-mlir.sh          # all problems
./scripts/verify-mlir.sh 1 50     # problems 1-50
```

### Manual MLIR compilation

```bash
export PATH="/opt/homebrew/opt/llvm/bin:$PATH"
export PYTHONPATH=$HOME/flow/src

# Flow -> MLIR
python3 -m flow.transpiler problems/p006.flow --mlir --lenient -o /tmp/p006.mlir

# MLIR -> LLVM IR
mlir-opt \
  --convert-arith-to-llvm \
  --convert-func-to-llvm \
  --convert-cf-to-llvm \
  --convert-index-to-llvm \
  --reconcile-unrealized-casts \
  /tmp/p006.mlir > /tmp/p006_lowered.mlir

mlir-translate --mlir-to-llvmir /tmp/p006_lowered.mlir -o /tmp/p006.ll

# LLVM IR -> binary
llc -O3 /tmp/p006.ll -o /tmp/p006.s
clang /tmp/p006.s -o /tmp/p006 -lm
/tmp/p006
```
