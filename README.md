# flow-euler

Project Euler solutions in Flow, used as a compiler test suite. Expressive syntax, C-speed binaries.

**See it in your browser:** [Open the solutions site](https://godofecht.github.io/flow-euler/)

## What this is

[Project Euler](https://projecteuler.net) is a set of maths puzzles that need a program to solve. This repo solves them in [Flow](https://github.com/flooooooooooow/flow), a statically typed language (every variable has a fixed type, checked before the program runs) that reads close to high-level maths. Flow is a compiled language: a compiler turns each `.flow` file into C code, and a C compiler turns that into a fast native program.

There is one solution file for each of 1007 problems, and each one prints a single answer. A checking script runs every solution and compares its output with the known answer in `answers.txt`. That makes the repo a large test suite for the Flow compiler: if a change to Flow breaks something, a wrong answer shows up here.

## What you can learn from it

- How to turn a maths problem into an algorithm, from sums of multiples (problem 1) to number theory (Computer Science, Mathematics).
- Prime sieves, greatest common divisors and modular arithmetic, in `lib/nt.flow` and the problem files.
- How a compiler can translate one language into another (here Flow into C) and why that can still run fast.
- How automated testing works: a script runs every program and checks each answer.
- How to measure speed fairly, by timing the same algorithm written two ways.

## Run it

### In your browser

Open the [solutions site](https://godofecht.github.io/flow-euler/). The front page is a table with one row per problem: its status (pass or skip), the program's output, the expected answer, the run time, the memory used and an estimate of how the run time grows with the input size. Select a problem number to see its page, with the Flow source, the C code Flow generated, the MLIR code, and the output.

### On your own computer

The scripts were written for and tested on macOS (the automatic checks ran on macOS 14 on Apple Silicon). Linux may work but has not been checked.

**You will need:**

- [Git](https://git-scm.com/downloads), to download the code.
- [Python 3](https://www.python.org/downloads/) with `numpy` (`python3 -m pip install numpy`). The Flow compiler used here is written in Python.
- [Clang](https://clang.llvm.org/), a C compiler. On a Mac, run `xcode-select --install` to get it.
- The Flow compiler, at the version described in step 2.

The Flow compiler needs care. These scripts call Flow's Python transpiler (a compiler that outputs source code in another language) with `python3 -m flow.transpiler`. On 30 September 2026 Flow removed that Python code from its main branch, and Flow 2.0.0 onwards is under a proprietary licence. The last version these solutions were checked against is Flow commit `37681b9`, from 10 September 2026, which is under the MIT licence.

1. Download this repo and go into it:

   ```bash
   git clone https://github.com/godofecht/flow-euler.git
   cd flow-euler
   ```

2. Download Flow next to it and switch to the commit the checks last passed with:

   ```bash
   git clone https://github.com/flooooooooooow/flow.git ~/flow
   git -C ~/flow checkout 37681b9
   export FLOW_REPO=$HOME/flow
   ```

3. Run problem 1:

   ```bash
   ./scripts/run.sh 1
   ```

   When it works you see `p001: 233168`.

4. Run every problem, or check every answer against `answers.txt`:

   ```bash
   ./scripts/run.sh all
   ./scripts/verify.sh
   ./scripts/verify.sh 1 50   # problems 1 to 50 only
   ```

   `verify.sh` prints `OK` for each correct answer and ends with a line such as `passed=1007 skipped=0 failed=0`. Some of the later problems take minutes each, so the full run is long.

Point `FLOW_BIN` at `$FLOW_REPO/flow` if you set it. A symlinked `~/.local/bin/flow` breaks Flow's `SCRIPT_DIR` resolution.

## How it works

Each problem lives in its own file, `problems/pNNN.flow`, where `NNN` is the problem number. The file defines `main()` and prints one answer. `scripts/run.sh` sends the file through the Flow transpiler to get a C file, compiles that with `clang -O3 -march=native` (the highest optimisation level, tuned for your processor), runs the program and keeps the last answer line. `scripts/verify.sh` does this for every line in `answers.txt` and reports any mismatch.

Problem 1 asks for the sum of all multiples of 3 or 5 below 1000. The solution in `problems/p001.flow` uses inclusion-exclusion (add the multiples of 3 and of 5, then subtract the multiples of 15, which were counted twice) and the formula for an arithmetic series:

```flow
function sum_ap(n: i64, step: i64) -> i64 {
    let k: i64 = (n - 1) / step
    return step * k * (k + 1) / 2
}

function main() -> i32 {
    println(sum_ap(1000, 3) + sum_ap(1000, 5) - sum_ap(1000, 15))
    return 0
}
```

`i64` is a 64-bit whole number. `i32` is a 32-bit one.

| File or folder | What is in it |
|---|---|
| `problems/` | One `.flow` file per Project Euler problem |
| `lib/` | Shared Flow code: `nt.flow` (number theory) and `digits.flow` (large numbers as digit arrays) |
| `data/` | Input files some problems read: grids, digit strings, tables |
| `answers.txt` | The expected answer for each problem, one per line |
| `scripts/` | `run.sh`, `verify.sh`, `verify-mlir.sh`, `bench.sh` and the site generator |
| `bench/` | A hand-written C version of problem 10 and its Flow twin, for timing |
| `docs/` | Notes on the layout, the two compiler back ends and adding a solution |
| `FLOW_IDIOMS.md` | Writing conventions learned from this set of solutions |

## Ideas to try

1. Open `problems/p001.flow` and change each `1000` to `10`. The puzzle says the answer for 10 is 23. Run `./scripts/run.sh 1` to check.
2. Read `problems/p014.flow`, which finds the longest Collatz chain under one million. Add a `printf` that also prints the length of that chain.
3. Add a new function to `lib/nt.flow`, such as Euler's totient (how many of the numbers from 1 to n share no factor with n other than 1). Import it in a problem file the same way `problems/p005.flow` imports `gcd` and `lcm`.
4. Pick a problem that still carries its own copy of `gcd` or `isqrt` and switch it to `import euler.nt`. `FLOW_IDIOMS.md` explains when this is worth doing. Run `./scripts/verify.sh` on that problem number afterwards.
5. Write problem 10 a different way in `problems/p010.flow` and compare your run time with the version in `bench/flow_p010.flow`.
6. Follow `docs/adding-solutions.md` to write a second solution to a problem you have solved yourself, and check it against `answers.txt`.

## Technical details

### Progress

| Range | Solved |
|-------|--------|
| 001-200 | 200 |
| 201-400 | 200 |
| 401-600 | 200 |
| 601-1007 | 407 |

All 1007 known answers in `answers.txt` have Flow solutions. The live site reports 1002 passed and 5 skipped, because the site generator stops each problem after 30 seconds.

For the conventions learned from this corpus, see [FLOW_IDIOMS.md](FLOW_IDIOMS.md). Its rule: adopt an idiom when repeated real solutions show that it improves intent or reuse without changing the low-level execution contract.

### Why Flow here

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
| Flow to C | 1675 |
| Hand-written C | 1767 |

Reproduce with `./scripts/bench.sh`. The script calls `$FLOW_REPO/flow compile` and reads the generated C from `$FLOW_REPO/build/flow_p010.c`.

### Shared libraries

`lib/nt.flow` is the module `euler.nt`, with common `i64` number-theory helpers: `gcd`, `lcm`, `isqrt`, `mulmod`, `mod_pow`, `is_prime`. `lib/digits.flow` is `euler.digits`, a decimal big integer stored least significant digit first. Import from a problem file:

```flow
import euler.nt { gcd, lcm, is_prime }
```

Resolution comes from the repo-root `flow.toml`:

```toml
[paths]
euler = "lib"   # import euler.nt -> lib/nt.flow
```

`flow run` and `scripts/run.sh` walk upward from the source file to find `flow.toml`, so no working-directory tricks are needed when launching from the repo root:

```bash
FLOW_BIN=/path/to/flow FLOW_HOST=python ./scripts/run.sh 5
# or
FLOW_HOST=python /path/to/flow run problems/p005.flow
```

### C backend

The default path. Flow transpiles to C99, then `clang -O3 -march=native` compiles to a native binary.

```text
pNNN.flow -> Flow transpiler (--c) -> pNNN.c -> clang -O3 -> pNNN -> answer
```

The command is `python3 -m flow.transpiler input.flow --c --lenient -o output.c`. `scripts/run.sh` runs the transpiler directly because the Flow shell launcher uses Bash namerefs, which the Bash 3.2 shipped on macOS lacks.

### MLIR backend

Flow also compiles through MLIR (Multi-Level Intermediate Representation, part of the LLVM compiler project). This repo tests both paths in CI.

```bash
brew install llvm               # provides mlir-opt, mlir-translate, llc
export LLVM_BIN=/opt/homebrew/opt/llvm/bin

./scripts/verify-mlir.sh        # verify through MLIR backend
./scripts/verify-mlir.sh 1 50   # problems 1-50 only
```

```text
pNNN.flow -> Flow transpiler (--mlir) -> pNNN.mlir
          -> mlir-opt (lower to LLVM dialect) -> pNNN_lowered.mlir
          -> mlir-translate (--mlir-to-llvmir) -> pNNN.ll
          -> llc -O3 -> pNNN.s
          -> clang -> pNNN -> answer
```

`verify-mlir.sh` reports three categories:

- `passed`: MLIR compilation succeeded and the answer matches.
- `mlir_fail`: MLIR lowering or compilation failed (a compiler limitation).
- `failed`: MLIR compilation succeeded but the answer was wrong.

Only `failed` causes a non-zero exit. `docs/backends.md` records a known limitation: the MLIR generator hardcodes `i32` in `arith.index_cast`, so problems with `i64` loop variables or stepped `for` ranges fail lowering. See [docs/backends.md](docs/backends.md) for the manual compilation steps.

### Native helpers

`scripts/run.sh` supports a companion C or C++ file for algorithms that need `__int128`, GMP, MPFR or NTT. If `problems/native/pNNN.c` or `problems/native/pNNN.cpp` exists, the script compiles the Flow file to C and links it with that file, adding GMP and MPFR when their Homebrew headers are present. The Flow file declares the native function as `extern`:

```flow
extern {
    function p775_native() -> i64
}

function main() -> i32 {
    printf("%lld\n", p775_native())
    return 0
}
```

The current tree has no `problems/native/` folder. Every solution, including problem 775, is pure Flow. See [docs/adding-solutions.md](docs/adding-solutions.md) for the full guide.

### CI

Four GitHub Actions workflows run on `macos-14`, each cloning the Flow compiler first:

- `verify.yml`: compiles every solution through the C backend and checks answers against `answers.txt`. Problems listed in `.ci-skip.txt` are skipped (too slow, crash, or broken); the list is empty today. 120-minute timeout.
- `mlir.yml`: compiles every pure-Flow solution through the MLIR backend (`Flow -> MLIR -> LLVM IR -> llc -> clang`). Reports MLIR lowering failures separately from answer mismatches. 60-minute timeout.
- `site.yml`: runs `scripts/measure.py` and `scripts/gen-site.py` to build the GitHub Pages site with the Flow source, generated C, generated MLIR, and output for every problem, then deploys it.
- `complexity.yml`: measures time, peak memory and estimated complexity for every problem and commits `complexity.json`.

### Documentation

- [FLOW_IDIOMS.md](FLOW_IDIOMS.md): corpus-derived Flow idioms and the rules for adopting them
- [AUDIT_FLOW_IDIOM.md](AUDIT_FLOW_IDIOM.md): an audit from 9 August 2026 of how idiomatic the solutions were at that time
- [docs/architecture.md](docs/architecture.md): repository structure, solution patterns, compilation paths
- [docs/backends.md](docs/backends.md): C and MLIR backend details, prerequisites, limitations
- [docs/adding-solutions.md](docs/adding-solutions.md): how to add a new solution, native helper conventions

## Project status

[![Verify (C backend)](https://github.com/godofecht/flow-euler/actions/workflows/verify.yml/badge.svg)](https://github.com/godofecht/flow-euler/actions/workflows/verify.yml)
[![MLIR backend](https://github.com/godofecht/flow-euler/actions/workflows/mlir.yml/badge.svg)](https://github.com/godofecht/flow-euler/actions/workflows/mlir.yml)

The last CI run, on 10 September 2026, passed all 1007 problems through the C backend. Since then Flow has removed the Python transpiler these scripts call, and the workflows clone Flow's main branch. The next CI run will fail until the workflows pin a Flow commit or the scripts move to Flow's `flowc` compiler.

## Licence

You can use, copy and change this code for any purpose, as long as you keep the copyright notice. See [LICENSE](LICENSE) (MIT). The Project Euler problems themselves remain (c) Project Euler.

The Flow compiler has its own licence. Versions up to 1.0.2, and the commit `37681b9` used above, are under the MIT licence. Flow 2.0.0 and later are proprietary.

## Credits

- Problems from [Project Euler](https://projecteuler.net).
- Built with the [Flow](https://github.com/flooooooooooow/flow) language, Clang, and LLVM/MLIR.
- GMP and MPFR are used by the native helper path when installed.
