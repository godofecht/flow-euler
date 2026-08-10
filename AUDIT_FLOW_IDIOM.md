# Flow-Idiomatic Concision Audit — flow-euler

**Date:** 2026-08-09  
**Scope:** ~710 `problems/p*.flow` solutions vs Flow docs (`~/flow/docs/`) and gold examples (`~/flow/examples/basics/`)  
**Constraint:** Recommend only rewrites that lower to the same class of tight C loops (no GC-y / heavily boxed styles).  
**Runner note:** This repo already defaults `FLOW_HOST=python` in `scripts/run.sh` — use that for the full surface (ranges, `|> sort`, spans, match). Stage-A `flowc` is a smaller subset.

---

## Executive summary

The corpus is a **C-shaped Project Euler port in Flow syntax**: correct, fast, and mostly readable, but it underuses the Flow features that still compile to plain C.

| Signal | Observation |
|--------|-------------|
| **Dominant style** | `while` + manual counters + `extern { calloc/free }` + duplicated helpers |
| **Range `for`** | Used in **14 / ~710** files (~2%); gold style is `p001` / `p005` / `p017` |
| **Structs** | **0** `struct` definitions in `problems/` |
| **`match` / `enum`** | **0** real uses (false positives only: variable names / comments) |
| **`import` / `lib/`** | **None** — every file reinlines `gcd`, `isqrt`, `mod_pow`, sieves, digit-bigint |
| **Pipelines / `|> sort`** | Unused; hand-rolled shell / insertion sorts instead (`p022`, `p124`, `p311`) |
| **`span<T>`** | Unused; ubiquitous `ptr + len` pairs |
| **`parallel for`** | Unused; several sieves / independent fills could try it safely |
| **`extern` / `calloc`** | ~552 files with `extern`, ~475 with `calloc` — appropriate for large heaps, overused for tiny buffers |
| **Thesis fit** | Early solved problems (`p001`–`p020` mixed) already show Flow can be short *and* fast; mid/hard files drift back to C templates |

**Bottom line:** Biggest wins are mechanical (range `for`, shared `lib/`), then packing digit-bigint / sieve state into structs + spans, then selective `|> sort` / `match` where control logic is long. Avoid rewriting hot modular-arithmetic kernels into pipelines.

Smoke-checked under `FLOW_HOST=python`: stepped `for i in 0..n step k` and `xs |> sort` both compile and run correctly.

---

## Feature opportunity matrix

| Feature | Would help (est.) | Perf after rewrite | Risk / maturity |
|---------|-------------------|--------------------|-----------------|
| **`for … in a..b` / `step` / `to`** | Very high — hundreds of simple counters | ≈ (same C `for`) | Low; Stage-A may lag — keep `FLOW_HOST=python` |
| **Shared `lib/` helpers** (`gcd`, `isqrt`, `mod_pow`, sieve, Miller–Rabin) | Very high — 58× `gcd`, 109× `isqrt`, 107× `mod_pow`, 24× `is_prime`, 59× sieve | ≈ or ↑ (one tuned impl) | Medium: needs `flow.toml` + module paths; Phase-1 modules |
| **Structs + methods** (digit bigint, sieve, poker hand) | High on ~20–40 digit/sieve files | ≈ | Low–medium; methods/codegen must stay unboxed |
| **`span<T>` / `&[T]`** | Medium — helpers taking `ptr+len` | ≈ (ptr+len under the hood) | Medium: concrete spans OK; inference not implemented |
| **`|> sort` / `sortBy` / `find`** | Medium — name lists, rad-order, small key sorts | ≈ for n≲10⁴; verify vs hand sort for larger | Low on C backend; policy keywords (`parallel`, `gpu`) not specialized |
| **`match` / `enum`** | Medium — poker categories, Roman digits, month lengths, coin kinds | ≈ | Match still ⚠️ for some exhaustiveness / `break` edges |
| **`parallel for`** | Selective — independent array init / some sieves | ↑ when OpenMP present; serial fallback OK | Medium: data races if marking conflicts; MLIR path serializes today |
| **Pipelines (`\|>`, `choose`)** | Low for Euler hot loops | ≈ if used only for outer control | Attractive demos; easy to over-abstract |
| **Effects / async / dynamics** | Near-zero for PE answers | n/a | Out of scope for this corpus |
| **Giant zero array literals** | Cleanup — ~10+ files | ≈ or ↑ (less frontend pain) | Prefer `calloc` / fill loop |

---

## Top 10 rewrite targets (biggest win)

Ranked for **clarity × lines saved × reuse × safe perf**.

| Rank | Target | Why |
|------|--------|-----|
| 1 | **New `lib/nt.flow` (+ `flow.toml`)** | Dedup `gcd`/`lcm`/`isqrt`/`mod_pow`/`mulmod`/trial/`miller_rabin` across 100+ files |
| 2 | **Range-`for` sweep on 001–100** | Highest density of `while i …; i = i + 1` with no control-flow tricks; gold already exists (`p001`, `p005`) |
| 3 | **Digit-bigint pack** — `p016`, `p020`, `p055`, `p056`, `p066` | Repeated `ptr + len + carry` → `struct Digits { data: …, len }` + `mul_small`/`add` methods |
| 4 | **`p022` names** | Replace shell sort + row swaps with score structs or keys + `|> sort` / `sortBy` |
| 5 | **`p124` rad(n)** | Parallel arrays `rad[]`/`idx[]` + shell sort → struct array or `sortBy` |
| 6 | **`p031` coin ways** | Tiny DP: range `for`, fixed `array` or one `calloc`, clearer stepped inner loop |
| 7 | **`p054` poker** | Long `if/elif` category chain → `enum HandCat` + `match`; ranks via `match` on char |
| 8 | **Sieve template** — `p010`, `p035`, `p049`, `p072` | Shared `sieve_mark` / `totient_sieve`; outer loops as `for`; optional `parallel for` on init |
| 9 | **`p019` Sundays** | Nested year/month counters → `for y in 1901..2001` + `match` month lengths |
| 10 | **Quality destubs** — `p302`, `p311`, `p331`, `p338`, `p342` (+ native wrappers) | Hardcoded `printf` answers / `peNNN_answer()` are not idiom wins — track separately |

---

## Annotated findings (25)

Each item: **smell → proposed shape → effect → risk**.

### A. Ranges / stepped `for` (high value, low risk)

1. **`problems/p009.flow`** — Nested `while a…; a = a + 1` / `b = b + 1`.  
   → `for a in 1..(perimeter/3)` + inner `for b in (a+1)..…` (keep `break`).  
   → lines↓, clarity↑, perf≈. Risk: low.

2. **`problems/p004.flow`** — Descending search with manual counters.  
   → Prefer keep `while` *or* document that descending ranges are awkward; optional `for a in 0..900` with `aa = 999 - a`.  
   → clarity↑ if rewritten carefully; perf≈. Risk: low (don’t force range if uglier).

3. **`problems/p012.flow`** — `while i * i <= n` for divisors (good); outer triangle loop is `while true`.  
   → Keep trial `while` (bound depends on `i*i`); only convert fixed-bound loops elsewhere.  
   → Teaching note more than rewrite. Risk: n/a.

4. **`problems/p017.flow`** — Already idiomatic (`for n in 1..1001`).  
   → Use as **gold** for letter/lookup tables.  
   → Effect: style reference. Risk: none.

5. **`problems/p019.flow`** — Year/month `while` nests.  
   → `for y in 1901..2001 { for m in 1..13 { … } }` + `match m` for lengths.  
   → lines↓, clarity↑, perf≈. Risk: low.

6. **`problems/p031.flow`** — `while ci < 8` / `while x <= target`.  
   → `for ci in 0..8` + `for x in c..(target+1)`; see exemplar below.  
   → lines↓, clarity↑, perf≈. Risk: low.

7. **`problems/p072.flow`** — Three identical `while i <= n` passes over φ.  
   → `for i in 0..(n+1)` / `for i in 2..(n+1)`; extract `totient_prefix_sum(n)`.  
   → clarity↑, reuse↑, perf≈. Risk: low.

8. **`problems/p010.flow`** — Classic sieve with `while p * p < limit` (keep) and sum loop `while i < limit`.  
   → Sum with `for i in 2..limit`; marking stays `while`/`step` by `p`. Optional `parallel for` only for zero-init if not using `calloc`.  
   → clarity↑, perf≈. Risk: low for ranges; medium for parallel mark.

### B. Pipelines / declarative sort

9. **`problems/p022.flow`** — Fixed-width rows + custom shell sort + `swap_rows`.  
   → Parse into `array` of scores/indices or keep buffer but sort an index array via structured keys; for ≤6000 names, `|> sort` / `sortBy` is fine.  
   → lines↓↓, clarity↑, perf≈. Risk: medium (string compare vs fixed-width — may need `array<i32>` keys or struct `{name_id, score}`).

10. **`problems/p124.flow`** — Parallel `rad[]`/`idx[]` + shell sort by `(rad, n)`.  
    → `struct RadEntry { rad: i32, n: i32 }` then `entries |> sortBy [asc .rad, asc .n]`.  
    → lines↓, clarity↑, perf≈ at n=1e5 (verify once). Risk: medium — confirm sort codegen cost vs shell sort.

11. **`problems/p311.flow`** — Hand insertion sort on parallel `xs`/`ys` (and final answer hardcoded — see stubs).  
    → Struct points + `|> sortBy`; still need real `B(10^10)` computation.  
    → clarity↑. Risk: medium (feature) + high (incomplete solution).

### C. Spans

12. **`problems/p054.flow`** — `hand_score(ranks: ptr<i32>, suits: ptr<i32>)` + internal `calloc(5)`.  
    → `function hand_score(ranks: span<i32>, suits: span<i32>)` and stack `array<i32,5>` for sorted copy.  
    → clarity↑, fewer heap ops, perf≈/↑. Risk: medium (span call-site borrow must work with arrays).

13. **`problems/p066.flow`** — Many `ptr<i32> + alen` digit helpers (`add_big`, `mul_big`, …).  
    → `span<mut i32>` parameters + `struct Big { digits: …, len: i32 }` (or heap buffer + span views).  
    → clarity↑, lines↓ across helpers. Risk: medium.

14. **Repeated sieve APIs** (`p010`, `p035`, `p049`, …) — `ptr<i8> sieve` alone is fine; helpers that take subranges should use `span<i8>`.  
    → Document as convention when extracting `lib/sieve.flow`.  
    → clarity↑. Risk: low.

### D. Match / enums

15. **`problems/p054.flow`** — Category ladder `if straight && flush / elif four / …` (~60 lines).  
    → `enum HandCat { High, Pair, TwoPair, … }` + `match` to build score; `rank_val` as `match c`.  
    → clarity↑, lines↓. Risk: medium (match maturity ⚠️).

16. **`problems/p089.flow`** — `roman_value` as long `if c == 73…`.  
    → `match c { 73 => 1, 86 => 5, … }` or small lookup table (table already used for lengths).  
    → clarity↑. Risk: low–medium.

17. **`problems/p019.flow`** — `days_in_month` if/elif chain.  
    → `match m { 2 => …, 4|6|9|11 => 30, default => 31 }` if alternation supported; else table.  
    → clarity↑. Risk: medium (alternation support).

### E. Structs + methods

18. **`problems/p016.flow` / `p020.flow`** — Ad-hoc digit arrays with carry loops.  
    → Shared `Digits` + `double()` / `mul_u32(n)` / `sum()`.  
    → lines↓↓ across digit problems, clarity↑, perf≈. Risk: low if methods inline to same loops.

19. **`problems/p055.flow` / `p056.flow` / `p066.flow`** — Same digit pattern, more ops (reverse, compare, mul).  
    → Same `Digits` module; `p066` is the stress test (Pell + big mul).  
    → reuse↑. Risk: medium (large digit buffers stay heap-backed).

20. **`problems/p010.flow` family** — Bare sieve buffer.  
    → Optional `struct Sieve { mark: ptr<i8>, n: i64 }` with `is_prime` / `sum` methods — only if it stays zero-cost.  
    → clarity↑ for multi-query problems. Risk: low.

### F. Typed helpers / imports

21. **Corpus-wide `function gcd` (58 files), `isqrt` (~109), `mod_pow`/`modpow` (~107), `is_prime` (24)** — Copy-paste drift (some use `i128` mulmod, some don’t).  
    → `lib/nt.flow` exporting one correct suite; problems `import euler.nt { gcd, mod_pow, … }` once `flow.toml` maps `euler = "lib"`.  
    → lines↓↓, correctness↑, perf≈/↑. Risk: medium (module resolution in this repo’s run script).

22. **`problems/p200.flow`** — Full Miller–Rabin with `calloc(7)` for fixed bases.  
    → `array<i64,7>` bases + shared `is_prime_u64` in `lib/`.  
    → lines↓, no heap for constants, perf≈/↑. Risk: low.

23. **`problems/p005.flow` / `p152.flow`** — Local `gcd`/`lcm`.  
    → First consumers of `lib/nt.flow`.  
    → clarity↑. Risk: low once imports work.

### G. Anti-patterns

24. **Giant zero-init literals** — e.g. `p152.flow` (`array<…,64>` of zeros × several), `p089.flow` 32-zero buf, `p172.flow` 20-zero counts.  
    → Prefer `calloc` / `let mut a: array<…>` + fill, or only materialize used prefix.  
    → clarity↑, frontend lighter, perf≈. Risk: low.

25. **Unnecessary `extern` for libc that Flow already exposes** — Many files redeclare `printf`/`calloc`/`free`; early files often rely on defaults.  
    → Declare `extern` only for real extras (`fopen`, native `peNNN_answer`, …).  
    → noise↓. Risk: low (confirm defaults per host).

### H. Parallel for (selective)

26. **Independent fills** — e.g. `p423.flow` odd-sieve init `while i < size { odd[i] = 1 }`.  
    → `parallel for i in 0 to size { odd[i] = 1 }` then serial sieve mark.  
    → perf↑ on large L with OpenMP; clarity≈. Risk: medium (worth only for large n; verify no races).

27. **Do *not* parallelize** marking loops that write conflicting indices without care (`sieve[m] = 1` for shared `m` strides is OK per prime sequentially; parallelizing outer `p` is unsafe without segmentation).  
    → Document in style guide. Risk: high if done naively.

---

## Exemplar before / after

### Exemplar 1 — PE 031 (coin DP): ranges + less C noise

**Before (current):**

```flow
let coins: array<i32, 8> = [1, 2, 5, 10, 20, 50, 100, 200]
let target: i32 = 200
let ways: ptr<i64> = calloc((target + 1) as i64, 8)
ways[0] = 1
let mut ci: i32 = 0
while ci < 8 {
    let c: i32 = coins[ci]
    let mut x: i32 = c
    while x <= target {
        ways[x] = ways[x] + ways[x - c]
        x = x + 1
    }
    ci = ci + 1
}
```

**After (proposed):**

```flow
let coins: array<i32, 8> = [1, 2, 5, 10, 20, 50, 100, 200]
let target: i32 = 200
let ways: ptr<i64> = calloc((target + 1) as i64, 8)
ways[0] = 1
for ci in 0..8 {
    let c: i32 = coins[ci]
    for x in c..(target + 1) {
        ways[x] = ways[x] + ways[x - c]
    }
}
```

Expected: fewer lines, same DP, same `-O3` loop structure. Heap kept (201 slots) — fine.

### Exemplar 2 — PE 001 (already gold) vs PE 009 style target

**Gold (`p001`) — keep teaching this:**

```flow
function sum_multiples(limit: i64, step: i64) -> i64 {
    let mut total: i64 = 0
    for n in 0..limit step step {
        total = total + n
    }
    return total
}
```

**`p009` today** uses manual `while` nests; idiomatic twin:

```flow
function solve(perimeter: i64) -> i64 {
    for a in 1..(perimeter / 3) {
        for b in (a + 1)..(perimeter / 2) {
            let c: i64 = perimeter - a - b
            if b >= c { break }
            if a * a + b * b == c * c {
                return a * b * c
            }
        }
    }
    return 0
}
```

Same asymptotics; clearer bounds. (`p005` remains the gold for `gcd`/`lcm` + range `for`.)

---

## Answer stubs / non-solutions (quality failures)

These are **not** idiom opportunities — they fail the “Flow solves Euler” thesis. Flag for destubbing, not for `|> sort`.

| File | Issue |
|------|--------|
| `problems/p218.flow` | Prints `0` only — *mathematically* the answer is 0; acceptable if comment stays rigorous, but it’s not an algorithm demo |
| `problems/p302.flow` | Hardcodes `1170060` after toy arithmetic checks |
| `problems/p311.flow` | Partial solver + hardcodes `2466018557` for full bound |
| `problems/p331.flow` | Implements `T_even` for checks, hardcodes full sum |
| `problems/p338.flow` | Computes small `G(N)`, hardcodes `G(10^12) mod 10^8` |
| `problems/p342.flow` | Statement check only; hardcodes full answer |
| `problems/p502.flow` (+ many peers) | Thin `extern pe502_answer()` wrappers — real work in `problems/native/` (~45 native TUs; ~19 `pe*_answer` Flow drivers; ~43 files ≤15 lines) |

Native wrappers are an honest engineering escape hatch; still mark them as **non-Flow solutions** in progress tables.

---

## Recommended “Flow Euler style guide”

For future `problems/pNNN.flow` (and rewrite passes):

1. **Prefer `for x in a..b` / `step k` / `0 to n`** whenever the bound is known up front; keep `while` for `i*i <= n`, digit peeling, and data-dependent loops.
2. **One answer line** via `printf("%lld\n", …)`; put the math in `solve(…)` / helpers.
3. **Reuse `lib/`** for `gcd`, `lcm`, `isqrt`, `mod_pow`, `mulmod`, sieves, Miller–Rabin — don’t fork 64-bit overflow behavior per file.
4. **Heap only when N is large or unknown**; small digit buffers → `array<i32, N>` or one `calloc` wrapped in a struct.
5. **Pack related state** (`Digits`, `Sieve`, `Hand`) instead of parallel bare arrays when it shortens call sites.
6. **Use `span<T>`** on helper boundaries (`sum_digits(digits: span<i32>)`) — never invent a second length.
7. **Declare sort intent** with `|> sort` / `sortBy` for n ≲ 1e5 demo sorts; keep custom radix/bucket sorts when the algorithm *is* the point.
8. **`match` / `enum`** for closed categories (hand type, roman digit, opcode) — not for hot modular arithmetic.
9. **`parallel for` only** for independent writes (fills, segmented sieves); never race the same `sieve[m]` from unsynchronized primes.
10. **No answer stubs** — if Flow can’t express it yet, keep a native TU *and* note it in the file header; don’t `printf` the constant silently.
11. **Stay on the C-shaped path** — no GC containers, no pipeline-of-closures in inner loops; high-level surface must lower to simple loops.
12. **Always run with `FLOW_HOST=python`** for this repo (already default in `scripts/run.sh`).

---

## Attractive but immature / unsafe for perf (today)

| Feature | Why it looked good | Why to hold back on PE kernels |
|---------|--------------------|--------------------------------|
| **Pipeline `choose` / long `|>` chains** | Pretty control demos (`examples/basics/pipeline_*.flow`) | Extra temps; obscures tight number-theory loops; little line-win on Euler |
| **Sort policies `parallel` / `gpu` / `simd`** | Documented on ordering surface | Parsed, **not specialized** yet |
| **`parallel for` everywhere** | OpenMP wins on big fills | Unsafe on classic sieve outer loops; MLIR path still serializes `scf.parallel` |
| **Effects / async / channels** | Core Flow story | Irrelevant latency model for PE batch solvers |
| **Bare / inferred `span`** | Ergonomic | **Not implemented** — must write `span<i32>` etc. |
| **String-path `import "…"`** | Easy relative includes | Deprecated; need package `[paths]` modules |
| **Stage-A `flowc` as default host** | Faster CI dream | Subset of syntax — this corpus should stay on Python host |
| **Dynamics / autodiff / GPU Metal** | Showcase tech | Wrong tool for integer Euler answers |
| **Boxing / GC-style vectors** | Python comfort | Fights the repo thesis (Flow → fast C) |

---

## Suggested follow-up rewrite pass (ordered)

1. Add `flow.toml` + `lib/nt.flow` (`gcd`, `lcm`, `isqrt`, `mod_pow`, `mulmod`).  
2. Convert **p001–p050** manual counters to range `for` (skip where uglier).  
3. Extract `lib/digits.flow` and retrofit `p016`/`p020`/`p055`/`p056`.  
4. Rewrite `p022` + `p124` with declarative sort; bench against current.  
5. Idiom-polish `p031`, `p019`, `p054`, `p072`, `p010`.  
6. Separate epic: destub hardcoded answers / absorb native solvers into Flow where feasible.

---

## Corpus snapshot (audit measurements)

| Metric | Count (approx.) |
|--------|-----------------|
| `problems/p*.flow` | ~710 |
| Files with `for … in` | 14 |
| Files with `while` | ~650+ |
| `struct` definitions | 0 |
| Real `enum` / `match` | 0 |
| `import` | 0 |
| Local `function gcd` | 58 |
| Local `function isqrt` | ~109 |
| Local `mod_pow`/`modpow` | ~107 |
| Local `function is_prime` | 24 |
| Sieve-style (`sieve[`) | ~59 |
| `extern {` | ~552 |
| `calloc` | ~475 |
| `pe*_answer` wrappers | ~19 |
| Files ≤15 lines | ~43 |
| Hardcoded large `printf` answers | 5 (`p302`,`p311`,`p331`,`p338`,`p342`) + trivial `p218` |

---

*End of audit. This document is analysis + exemplars only; no mass rewrite was performed.*
