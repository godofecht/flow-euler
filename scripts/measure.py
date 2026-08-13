#!/usr/bin/env python3
"""Measure time, peak memory, and estimated complexity for each problem.

Compiles each solution, runs it under /usr/bin/time -l, and statically
analyses the generated C to estimate Big-O time and space complexity.

Outputs a JSON file (default: complexity.json) with per-problem metrics.

Usage:
    python3 scripts/measure.py [--out FILE] [--limit N] [--timeout SECS]
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FLOW_REPO = os.environ.get("FLOW_REPO", str(Path.home() / "flow"))
FLOW_PY = str(Path(FLOW_REPO) / "src")

ANSWERS = {}
with open(ROOT / "answers.txt") as f:
    for line in f:
        line = line.strip()
        if line and line[0].isdigit():
            parts = line.split(None, 1)
            if len(parts) == 2:
                ANSWERS[int(parts[0])] = parts[1]

SKIP = set()
skip_file = ROOT / ".ci-skip.txt"
if skip_file.exists():
    with open(skip_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if parts:
                SKIP.add(int(parts[0]))


def has_native(n: int) -> bool:
    label = f"{n:03d}"
    return (ROOT / "problems/native" / f"p{label}.c").exists() or \
           (ROOT / "problems/native" / f"p{label}.cpp").exists()


def compile_problem(n: int) -> Path | None:
    """Compile a problem and return the executable path, or None on failure."""
    label = f"{n:03d}"
    src = ROOT / f"problems/p{label}.flow"
    if not src.exists():
        return None

    build = ROOT / ".build"
    build.mkdir(exist_ok=True)
    c_out = build / f"p{label}_measure.c"
    exe = build / f"p{label}_measure"

    native_c = ROOT / f"problems/native/p{label}.c"
    native_cpp = ROOT / f"problems/native/p{label}.cpp"
    brew_prefix = os.environ.get("HOMEBREW_PREFIX", "/opt/homebrew")

    env = os.environ.copy()
    env["PYTHONPATH"] = FLOW_PY + (":" + env.get("PYTHONPATH", ""))
    env["FLOW_HOST"] = "python"

    result = subprocess.run(
        [sys.executable, "-m", "flow.transpiler", str(src), "--c", "--lenient", "-o", str(c_out)],
        capture_output=True, text=True, env=env, timeout=60,
    )
    if not c_out.exists():
        return None

    libs = ["-lm"]
    if Path(f"{brew_prefix}/include/gmp.h").exists():
        libs += [f"-L{brew_prefix}/lib", f"-I{brew_prefix}/include", "-lgmp"]
    if Path(f"{brew_prefix}/include/mpfr.h").exists():
        libs += [f"-L{brew_prefix}/lib", f"-I{brew_prefix}/include", "-lmpfr"]

    if native_c.exists():
        cmd = ["clang", "-O3", "-march=native", f"-I{brew_prefix}/include",
               f"-L{brew_prefix}/lib", str(c_out), str(native_c), "-o", str(exe)] + libs
    elif native_cpp.exists():
        flow_o = build / f"p{label}_measure.o"
        subprocess.run(["clang", "-O3", "-march=native", "-c", str(c_out), "-o", str(flow_o)],
                       capture_output=True, timeout=60)
        cmd = ["clang++", "-O3", "-march=native", "-std=c++17", "-pthread",
               f"-I{brew_prefix}/include", f"-L{brew_prefix}/lib",
               str(flow_o), str(native_cpp), "-o", str(exe)] + libs
    else:
        cmd = ["clang", "-O3", "-march=native", str(c_out), "-o", str(exe), "-lm"]

    try:
        subprocess.run(cmd, capture_output=True, timeout=120)
    except subprocess.TimeoutExpired:
        return None
    if not exe.exists():
        return None
    return exe


def measure_runtime(exe: Path, timeout: int) -> dict:
    """Run the executable under /usr/bin/time -l and parse results."""
    try:
        proc = subprocess.run(
            ["/usr/bin/time", "-l", str(exe)],
            capture_output=True, text=True, timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {"time_ms": -1, "peak_rss_kb": -1, "output": "(timeout)"}

    # /usr/bin/time writes stats to stderr on macOS
    stderr = proc.stderr
    time_ms = -1
    peak_rss_kb = -1

    # Parse "X.XX real" for wall-clock time
    m = re.search(r"(\d+\.\d+)\s+real", stderr)
    if m:
        time_ms = int(float(m.group(1)) * 1000)

    # Parse "N  maximum resident set size" (in bytes on macOS)
    m = re.search(r"(\d+)\s+maximum resident set size", stderr)
    if m:
        peak_rss_kb = int(m.group(1)) // 1024

    # Extract last output line
    out = proc.stdout.strip()
    lines = [l for l in out.splitlines() if l.strip()]
    output = lines[-1] if lines else "(no output)"

    return {"time_ms": time_ms, "peak_rss_kb": peak_rss_kb, "output": output}


def estimate_complexity(c_path: Path) -> dict:
    """Statically analyse generated C to estimate Big-O time and space complexity.

    Looks at loop nesting depth per function, recursion, sorting calls,
    and allocation patterns. Returns approximate complexity classes.
    """
    if not c_path.exists():
        return {"time": "?", "space": "?", "loop_depth": 0, "has_recursion": False,
                "has_sort": False, "has_hashmap": False, "allocs": 0}

    code = c_path.read_text(errors="replace")

    # Strip comments and string literals for cleaner analysis.
    cleaned = []
    i = 0
    in_string = False
    in_char = False
    in_line_comment = False
    in_block_comment = False
    while i < len(code):
        c = code[i]
        if in_line_comment:
            if c == "\n":
                in_line_comment = False
                cleaned.append(c)
            i += 1
            continue
        if in_block_comment:
            if c == "*" and i + 1 < len(code) and code[i + 1] == "/":
                in_block_comment = False
                i += 2
                continue
            i += 1
            continue
        if in_string:
            if c == "\\":
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
            continue
        if in_char:
            if c == "\\":
                i += 2
                continue
            if c == "'":
                in_char = False
            i += 1
            continue
        if c == "/" and i + 1 < len(code):
            if code[i + 1] == "/":
                in_line_comment = True
                i += 2
                continue
            if code[i + 1] == "*":
                in_block_comment = True
                i += 2
                continue
        if c == '"':
            in_string = True
            i += 1
            continue
        if c == "'":
            in_char = True
            i += 1
            continue
        cleaned.append(c)
        i += 1
    clean_code = "".join(cleaned)

    # Extract function definitions with their bodies.
    # Pattern: return_type function_name(params) { body }
    # We use brace matching to find the full body.
    func_pattern = re.compile(
        r"(?:static\s+)?(?:inline\s+)?[\w\s\*]+?\s+(\w+)\s*\([^;]*\)\s*\{"
    )
    functions = {}  # name -> body text
    for m in func_pattern.finditer(clean_code):
        fname = m.group(1)
        # Skip common non-function matches
        if fname in ("if", "while", "for", "switch", "sizeof", "return"):
            continue
        brace_start = m.end() - 1
        depth = 0
        j = brace_start
        while j < len(clean_code):
            if clean_code[j] == "{":
                depth += 1
            elif clean_code[j] == "}":
                depth -= 1
                if depth == 0:
                    body = clean_code[brace_start + 1:j]
                    functions[fname] = body
                    break
            j += 1

    # Per-function loop depth: scan each function body for max for/while nesting.
    def max_loop_depth(body: str) -> int:
        max_d = 0
        d = 0
        k = 0
        while k < len(body):
            if re.match(r"\b(for|while)\s*\(", body[k:k + 20]):
                d += 1
                if d > max_d:
                    max_d = d
            if body[k] == "}":
                if d > 0:
                    d -= 1
            k += 1
        return max_d

    # Find the max loop depth across all functions.
    global_max_depth = 0
    for fname, body in functions.items():
        d = max_loop_depth(body)
        if d > global_max_depth:
            global_max_depth = d

    # Detect recursion: a function body contains a call to itself.
    has_recursion = False
    for fname, body in functions.items():
        if re.search(rf"\b{re.escape(fname)}\s*\(", body):
            has_recursion = True
            break

    # Detect sorting
    has_sort = bool(re.search(r"\bqsort\s*\(|\bsort\s*\(", clean_code))

    # Detect hashmap / hash table patterns
    has_hashmap = bool(re.search(r"\bhash\b|\bslot\b|\bcap\b|\bCAP\b", clean_code, re.IGNORECASE))

    # Count calloc/malloc calls
    allocs = len(re.findall(r"\bcalloc\s*\(|\bmalloc\s*\(", clean_code))

    # Look at calloc size arguments for space estimation.
    # If any calloc uses a variable (not constant), it's O(n) space.
    has_dynamic_alloc = False
    for m in re.finditer(r"\bcalloc\s*\(([^,]+),", clean_code):
        arg = m.group(1).strip()
        # If the first arg is not a plain integer literal, it's dynamic.
        if not re.match(r"^\d+$", arg):
            has_dynamic_alloc = True

    # Classify time complexity.
    # Recursion with branching (multiple self-calls in body) suggests exponential.
    # Recursion with single self-call suggests linear or logarithmic.
    if has_recursion:
        # Check if any recursive function has 2+ self-calls (branching recursion).
        branching = False
        for fname, body in functions.items():
            calls = re.findall(rf"\b{re.escape(fname)}\s*\(", body)
            if len(calls) >= 2:
                branching = True
                break
        if branching:
            time_complexity = "O(2^n)"
        elif global_max_depth >= 1:
            time_complexity = "O(n log n)"
        else:
            time_complexity = "O(log n)"
    elif has_sort:
        if global_max_depth >= 2:
            time_complexity = "O(n^2 log n)"
        else:
            time_complexity = "O(n log n)"
    elif global_max_depth == 0:
        time_complexity = "O(1)"
    elif global_max_depth == 1:
        time_complexity = "O(n)"
    elif global_max_depth == 2:
        time_complexity = "O(n^2)"
    elif global_max_depth == 3:
        time_complexity = "O(n^3)"
    elif global_max_depth == 4:
        time_complexity = "O(n^4)"
    else:
        time_complexity = f"O(n^{global_max_depth})"

    # Classify space complexity.
    if allocs == 0:
        space_complexity = "O(1)"
    elif has_dynamic_alloc:
        if allocs >= 5:
            space_complexity = "O(n^2)"
        else:
            space_complexity = "O(n)"
    else:
        space_complexity = "O(1)"

    return {
        "time": time_complexity,
        "space": space_complexity,
        "loop_depth": global_max_depth,
        "has_recursion": has_recursion,
        "has_sort": has_sort,
        "has_hashmap": has_hashmap,
        "allocs": allocs,
    }


def main():
    parser = argparse.ArgumentParser(description="Measure complexity of Flow Euler solutions")
    parser.add_argument("--out", default="complexity.json", help="output JSON file")
    parser.add_argument("--limit", type=int, default=0, help="limit number of problems (0 = all)")
    parser.add_argument("--timeout", type=int, default=60, help="run timeout per problem (seconds)")
    args = parser.parse_args()

    results = []
    count = 0

    for n in sorted(ANSWERS.keys()):
        if args.limit and count >= args.limit:
            break
        label = f"{n:03d}"
        count += 1

        if n in SKIP:
            results.append({
                "problem": n,
                "label": label,
                "status": "skipped",
                "time_ms": -1,
                "peak_rss_kb": -1,
                "output": "(skipped)",
                "expected": ANSWERS[n],
                "native": has_native(n),
                "complexity": {"time": "?", "space": "?"},
            })
            print(f"  p{label}: skipped", flush=True)
            continue

        print(f"  p{label}: compiling...", flush=True)
        exe = compile_problem(n)
        if exe is None:
            results.append({
                "problem": n,
                "label": label,
                "status": "compile-failed",
                "time_ms": -1,
                "peak_rss_kb": -1,
                "output": "(compile failed)",
                "expected": ANSWERS[n],
                "native": has_native(n),
                "complexity": {"time": "?", "space": "?"},
            })
            print(f"  p{label}: compile failed", flush=True)
            continue

        # Estimate complexity from generated C
        c_path = ROOT / ".build" / f"p{label}_measure.c"
        complexity = estimate_complexity(c_path)

        # Measure runtime
        print(f"  p{label}: measuring...", flush=True)
        metrics = measure_runtime(exe, args.timeout)

        expected = ANSWERS[n]
        if metrics["output"] == expected:
            status = "pass"
        elif metrics["output"].startswith("(timeout"):
            status = "timeout"
        elif metrics["output"].startswith("(no output"):
            status = "no-output"
        else:
            status = "mismatch"

        results.append({
            "problem": n,
            "label": label,
            "status": status,
            "time_ms": metrics["time_ms"],
            "peak_rss_kb": metrics["peak_rss_kb"],
            "output": metrics["output"],
            "expected": expected,
            "native": has_native(n),
            "complexity": complexity,
        })
        print(f"  p{label}: {status} time={metrics['time_ms']}ms rss={metrics['peak_rss_kb']}KB "
              f"time_O={complexity['time']} space_O={complexity['space']}", flush=True)

    out_path = ROOT / args.out
    out_path.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {len(results)} results to {out_path}")


if __name__ == "__main__":
    main()
