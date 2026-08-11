# Adding a New Solution

## 1. Find or compute the answer

Add the answer to `answers.txt`:

```
NNN answer
```

Integer answers go on one line. Floating-point answers should match the
precision Project Euler expects (usually 10-12 significant digits).

## 2. Write the Flow solution

Create `problems/pNNN.flow`. The file must define `main()` and print
the answer as the last non-empty line.

```flow
# Project Euler NNN
# Brief description.

function solve() -> i64 {
    # ... your algorithm ...
    return 42
}

function main() -> i32 {
    printf("%lld\n", solve())
    return 0
}
```

### Type conventions

- Integer answers: `i64`, print with `%lld`
- Floating-point answers: `f64`, print with `%.10f` or `%.12e`
- Large integers (beyond i64): print from native code directly

## 3. Test

```bash
./scripts/run.sh NNN
```

Output should be `pNNN: <answer>`.

## 4. Verify

```bash
./scripts/verify.sh NNN NNN
```

## 5. Native helper (if needed)

If the algorithm needs `__int128`, GMP, MPFR, NTT, or other facilities
not available in Flow, create a native C helper.

Create `problems/native/pNNN.c`:

```c
#include <stdint.h>

long long pNNN_native(void) {
    /* ... algorithm ... */
    return 42;
}
```

Update `problems/pNNN.flow` to call it:

```flow
extern {
    function pNNN_native() -> i64
}

function main() -> i32 {
    printf("%lld\n", pNNN_native())
    return 0
}
```

### Floating-point native helpers

```c
double pNNN_native(void) {
    return 3.14159265358979;
}
```

```flow
extern {
    function pNNN_native() -> f64
}

function main() -> i32 {
    printf("%.10f\n", pNNN_native())
    return 0
}
```

### C++ helpers

Name the file `problems/native/pNNN.cpp`. The build system compiles
the Flow-generated C as C, then links with the C++ TU. Use `extern "C"`
for the exported function.

### GMP / MPFR

The build system auto-detects Homebrew `gmp` and `mpfr` and links them
if the headers are present at `/opt/homebrew/include/`. Just include
the headers in your C file.

### Answers exceeding i64

For answers that don't fit in `i64` (e.g. p974's 29-digit integer),
have the native code print the answer directly and return 0:

```c
void pNNN_print(void) {
    /* compute and print the answer */
    printf("13313751171933973557517973175\n");
}
```

```flow
extern {
    function pNNN_print()
}

function main() -> i32 {
    pNNN_print()
    return 0
}
```

## 6. Data files

If the problem needs input data, put it in `data/` with a descriptive
name (e.g. `data/p673_beds.txt`). Reference it from the native helper
using a relative path.

## 7. Update README

Update the progress table in `README.md` if the solution count changed.
