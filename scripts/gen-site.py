#!/usr/bin/env python3
"""Generate a static GitHub Pages site showing Flow source, generated C,
generated MLIR, and the output for every Project Euler solution.

Usage:
    python3 scripts/gen-site.py [--out DIR] [--timeout SECS] [--limit N]

Outputs HTML files under the given directory (default: site/).
Skips problems in .ci-skip.txt to avoid timeouts.
"""
import argparse
import html
import json
import os
import re
import subprocess
import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

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

FLOW_REPO = os.environ.get("FLOW_REPO", str(Path.home() / "flow"))
FLOW_PY = str(Path(FLOW_REPO) / "src")

# Load complexity data if available.
COMPLEXITY = {}
complexity_file = ROOT / "complexity.json"
if complexity_file.exists():
    with open(complexity_file) as f:
        for entry in json.load(f):
            COMPLEXITY[entry["problem"]] = entry

# Load best-known reference data.
BEST_KNOWN = {}
bk_file = ROOT / "scripts" / "best-known.json"
if bk_file.exists():
    with open(bk_file) as f:
        raw = json.load(f)
        for k, v in raw.items():
            BEST_KNOWN[int(k)] = v


def extract_description(n: int) -> str:
    """Extract the problem description from the comment header of a .flow file."""
    label = f"{n:03d}"
    src = ROOT / f"problems/p{label}.flow"
    if not src.exists():
        return ""
    text = src.read_text(errors="replace")
    lines = text.split("\n")
    comments = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            comments.append(stripped[1:].strip())
        elif stripped == "":
            continue
        else:
            break
    # Skip "Project Euler NNN" header line
    desc_parts = []
    for c in comments[1:]:
        if c.startswith("Project Euler"):
            continue
        if c == "":
            continue
        desc_parts.append(c)
    return " ".join(desc_parts)


def complexity_match(ours: str, best: str) -> str:
    """Compare our complexity string against the best known."""
    if not ours or not best or ours == "?" or best == "?":
        return "unknown"
    # Normalise for comparison
    ours_n = ours.replace(" ", "").lower()
    best_n = best.replace(" ", "").lower()
    if ours_n == best_n:
        return "match"
    # Rough ordering of common complexity classes (best to worst)
    order = ["o(1)", "o(loglogn)", "o(logn)", "o(sqrt(n))",
             "o(n)", "o(nloglogn)", "o(nlogn)", "o(n*sqrt(n))",
             "o(n^2)", "o(n^2logn)", "o(n^3)", "o(n^4)",
             "o(2^n)", "o(n!)"]
    try:
        oi = order.index(ours_n)
        bi = order.index(best_n)
    except ValueError:
        return "unknown"
    if oi <= bi:
        return "match"
    return "worse"


def transpile(src_path: Path, mode: str, out_path: Path) -> str:
    """Transpile a .flow file to C or MLIR. Returns the output text."""
    cmd = [
        sys.executable, "-m", "flow.transpiler",
        str(src_path),
        f"--{mode}",
        "--lenient",
        "-o", str(out_path),
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = FLOW_PY + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["FLOW_HOST"] = "python"
    result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=60)
    if out_path.exists():
        return out_path.read_text(errors="replace")
    err = result.stderr.strip()
    if err:
        return f"(transpilation failed: {err[:500]})"
    return "(transpilation failed)"


def run_problem(n: int, timeout: int) -> str:
    """Run a problem through scripts/run.sh and return the output."""
    env = os.environ.copy()
    env["FLOW_HOST"] = "python"
    if "FLOW_REPO" not in env:
        env["FLOW_REPO"] = str(Path.home() / "flow")
    try:
        result = subprocess.run(
            ["bash", str(ROOT / "scripts/run.sh"), str(n)],
            capture_output=True, text=True, timeout=timeout,
            cwd=str(ROOT), env=env,
        )
        out = result.stdout.strip()
        for line in out.splitlines():
            if line.startswith("p"):
                parts = line.split(None, 1)
                if len(parts) == 2:
                    return parts[1]
        err = result.stderr.strip()
        if err:
            return f"(no output: {err[:500]})"
        return out.splitlines()[-1] if out else "(no output)"
    except subprocess.TimeoutExpired:
        return "(timeout)"
    except Exception as e:
        return f"(error: {e})"


def has_native(n: int) -> bool:
    label = f"{n:03d}"
    return (ROOT / "problems/native" / f"p{label}.c").exists() or \
           (ROOT / "problems/native" / f"p{label}.cpp").exists()


def esc(text: str) -> str:
    return html.escape(text)


def page_template(n: int, flow_src: str, c_src: str, mlir_src: str,
                  output: str, expected: str, native: bool) -> str:
    label = f"{n:03d}"
    title = f"PE {label}"
    status = "PASS" if output == expected else "DIFF"
    if output.startswith("(timeout") or output.startswith("(error"):
        status = "SKIP"

    native_note = '<p class="native-note">This problem uses a native C/C++ helper. MLIR is generated from the Flow wrapper only and may not run standalone.</p>' if native else ''

    # Problem description and PE link
    desc = extract_description(n)
    desc_html = f'<p class="problem-desc">{esc(desc)}</p>' if desc else ''
    pe_link = f'<p class="pe-link"><a href="https://projecteuler.net/problem={n}">View problem on Project Euler</a></p>'

    # Complexity metrics from complexity.json (if available).
    cx = COMPLEXITY.get(n, {})
    cx_rows = ""
    our_time_o = "?"
    our_space_o = "?"
    if cx:
        time_ms = cx.get("time_ms", -1)
        rss_kb = cx.get("peak_rss_kb", -1)
        comp = cx.get("complexity", {})
        our_time_o = comp.get("time", "?")
        our_space_o = comp.get("space", "?")
        time_str = f"{time_ms} ms" if time_ms >= 0 else "n/a"
        rss_str = f"{rss_kb} KB" if rss_kb >= 0 else "n/a"
        cx_rows = f"""<tr><th>Runtime</th><td>{time_str}</td></tr>
<tr><th>Peak memory</th><td>{rss_str}</td></tr>
<tr><th>Time complexity</th><td><code>{esc(our_time_o)}</code> <span class="est-note">(estimated)</span></td></tr>
<tr><th>Space complexity</th><td><code>{esc(our_space_o)}</code> <span class="est-note">(estimated)</span></td></tr>"""

    # Best-known comparison
    bk = BEST_KNOWN.get(n, {})
    bk_time = bk.get("best_time", "?")
    bk_space = bk.get("best_space", "?")
    bk_approach = bk.get("approach", "Not curated")
    match = complexity_match(our_time_o, bk_time)
    match_label = {"match": "Optimal", "worse": "Suboptimal", "unknown": "Unknown"}[match]
    match_cls = f"cmp-{match}"

    comparison_section = f"""<h2>Performance comparison</h2>
<table class="comparison">
<thead><tr><th>Metric</th><th>Our solution</th><th>Best known</th></tr></thead>
<tbody>
<tr><th>Time complexity</th><td><code>{esc(our_time_o)}</code></td><td><code>{esc(bk_time)}</code></td></tr>
<tr><th>Space complexity</th><td><code>{esc(our_space_o)}</code></td><td><code>{esc(bk_space)}</code></td></tr>
<tr><th>Approach</th><td>Flow solution</td><td>{esc(bk_approach)}</td></tr>
<tr><th>Verdict</th><td colspan="2" class="{match_cls}">{match_label}</td></tr>
</tbody>
</table>"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title} - Flow x Project Euler</title>
<link rel="stylesheet" href="../style.css">
</head>
<body>
<nav><a href="../index.html">&larr; All problems</a></nav>
<h1>Problem {label}</h1>
{desc_html}
{pe_link}
<table class="meta">
<tr><th>Answer</th><td><code>{esc(expected)}</code></td></tr>
<tr><th>Output</th><td><code>{esc(output)}</code></td></tr>
<tr><th>Status</th><td class="status-{status.lower()}">{status}</td></tr>
<tr><th>Native helper</th><td>{"yes" if native else "no"}</td></tr>
{cx_rows}
</table>
{comparison_section}
{native_note}
<h2>Flow source</h2>
<pre><code>{esc(flow_src)}</code></pre>
<h2>Generated C</h2>
<pre><code>{esc(c_src)}</code></pre>
<h2>Generated MLIR</h2>
<pre><code>{esc(mlir_src)}</code></pre>
</body>
</html>
"""


def index_template(entries: list) -> str:
    rows = []
    for n, status, output, expected, native in entries:
        label = f"{n:03d}"
        cls = f"status-{status.lower()}"
        cx = COMPLEXITY.get(n, {})
        time_str = ""
        rss_str = ""
        time_o = ""
        space_o = ""
        if cx:
            t = cx.get("time_ms", -1)
            r = cx.get("peak_rss_kb", -1)
            comp = cx.get("complexity", {})
            time_o = comp.get("time", "")
            space_o = comp.get("space", "")
            time_str = f"{t} ms" if t >= 0 else ""
            rss_str = f"{r} KB" if r >= 0 else ""
        bk = BEST_KNOWN.get(n, {})
        bk_time = bk.get("best_time", "?")
        match = complexity_match(time_o, bk_time)
        match_label = {"match": "=", "worse": ">", "unknown": "?"}[match]
        match_cls = f"cmp-{match}"
        rows.append(
            f'<tr><td><a href="problems/p{label}.html">{label}</a></td>'
            f'<td class="{cls}">{status}</td>'
            f'<td><code>{esc(output)}</code></td>'
            f'<td><code>{esc(expected)}</code></td>'
            f'<td>{"native" if native else "flow"}</td>'
            f'<td>{time_str}</td>'
            f'<td>{rss_str}</td>'
            f'<td><code>{esc(time_o)}</code></td>'
            f'<td><code>{esc(space_o)}</code></td>'
            f'<td><code>{esc(bk_time)}</code></td>'
            f'<td class="{match_cls}">{match_label}</td></tr>'
        )
    body = "\n".join(rows)
    total = len(entries)
    passed = sum(1 for e in entries if e[1] == "PASS")
    skipped = sum(1 for e in entries if e[1] == "SKIP")
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Flow x Project Euler</title>
<link rel="stylesheet" href="style.css">
</head>
<body>
<h1>Flow x Project Euler</h1>
<p>Solving <a href="https://projecteuler.net">Project Euler</a> in <a href="https://github.com/flooooooooooow/flow">Flow</a>.</p>
<p>Each problem page shows the problem statement, Flow source, generated C, generated MLIR, the program output, and a performance comparison against the best known approach.</p>
<p>{passed} passed, {skipped} skipped, {total} total.</p>
<table>
<thead><tr><th>Problem</th><th>Status</th><th>Output</th><th>Expected</th><th>Type</th><th>Time</th><th>Memory</th><th>Time O</th><th>Space O</th><th>Best O</th><th>Cmp</th></tr></thead>
<tbody>
{body}
</tbody>
</table>
</body>
</html>
"""


CSS = """\
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; max-width: 960px; margin: 0 auto; padding: 1.5rem; color: #1a1a1a; }
nav { margin-bottom: 1rem; }
h1 { font-size: 1.5rem; }
h2 { font-size: 1.2rem; margin-top: 2rem; border-bottom: 1px solid #ddd; padding-bottom: 0.3rem; }
table { border-collapse: collapse; width: 100%; font-size: 0.9rem; }
th, td { text-align: left; padding: 0.3rem 0.5rem; border-bottom: 1px solid #eee; }
th { background: #f5f5f5; }
code { font-family: "SF Mono", "Fira Code", monospace; font-size: 0.85rem; }
pre { background: #f8f8f8; border: 1px solid #e0e0e0; border-radius: 4px; padding: 0.8rem; overflow-x: auto; font-size: 0.8rem; line-height: 1.4; }
pre code { font-size: 0.8rem; }
.meta th { width: 8rem; }
.status-pass { color: #2a7d2a; }
.status-diff { color: #c62828; font-weight: bold; }
.status-skip { color: #666; }
.native-note { color: #666; font-size: 0.85rem; font-style: italic; }
.est-note { color: #999; font-size: 0.8rem; }
a { color: #1565c0; text-decoration: none; }
a:hover { text-decoration: underline; }
.problem-desc { font-size: 1rem; line-height: 1.5; margin: 0.5rem 0 1rem; color: #333; }
.pe-link { font-size: 0.9rem; margin: 0 0 1.5rem; }
.comparison { font-size: 0.9rem; }
.comparison th { width: 8rem; }
.cmp-match { color: #2a7d2a; font-weight: bold; }
.cmp-worse { color: #c62828; font-weight: bold; }
.cmp-unknown { color: #666; }
"""


def main():
    parser = argparse.ArgumentParser(description="Generate GitHub Pages site")
    parser.add_argument("--out", default="site", help="output directory")
    parser.add_argument("--timeout", type=int, default=30, help="run timeout per problem")
    parser.add_argument("--limit", type=int, default=0, help="limit number of problems (0 = all)")
    parser.add_argument("--skip-run", action="store_true", help="skip running problems (use expected as output)")
    args = parser.parse_args()

    out_dir = Path(args.out)
    problems_dir = out_dir / "problems"
    problems_dir.mkdir(parents=True, exist_ok=True)
    build_dir = ROOT / ".build" / "site"
    build_dir.mkdir(parents=True, exist_ok=True)

    (out_dir / "style.css").write_text(CSS)
    (out_dir / ".nojekyll").write_text("")

    entries = []
    count = 0
    for n in sorted(ANSWERS.keys()):
        if args.limit and count >= args.limit:
            break
        if n in SKIP:
            entries.append((n, "SKIP", "(skipped)", ANSWERS[n], has_native(n)))
            count += 1
            continue

        label = f"{n:03d}"
        src = ROOT / f"problems/p{label}.flow"
        if not src.exists():
            entries.append((n, "SKIP", "(no source)", ANSWERS[n], False))
            count += 1
            continue

        print(f"  generating p{label}...", flush=True)
        flow_src = src.read_text(errors="replace")
        native = has_native(n)

        c_path = build_dir / f"p{label}.c"
        mlir_path = build_dir / f"p{label}.mlir"
        c_src = transpile(src, "c", c_path)
        mlir_src = transpile(src, "mlir", mlir_path)

        if args.skip_run:
            output = ANSWERS[n]
        else:
            output = run_problem(n, args.timeout)

        expected = ANSWERS[n]
        if output == expected:
            status = "PASS"
        elif output.startswith("(timeout") or output.startswith("(error") or output.startswith("(no output"):
            status = "SKIP"
        else:
            status = "DIFF"

        page = page_template(n, flow_src, c_src, mlir_src, output, expected, native)
        (problems_dir / f"p{label}.html").write_text(page)

        entries.append((n, status, output, expected, native))
        count += 1

    index = index_template(entries)
    (out_dir / "index.html").write_text(index)
    print(f"Site generated in {out_dir} ({count} problems)")


if __name__ == "__main__":
    main()
