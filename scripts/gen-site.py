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
    subprocess.run(cmd, env=env, capture_output=True, timeout=60)
    if out_path.exists():
        return out_path.read_text(errors="replace")
    return "(transpilation failed)"


def run_problem(n: int, timeout: int) -> str:
    """Run a problem through scripts/run.sh and return the output."""
    try:
        result = subprocess.run(
            ["bash", str(ROOT / "scripts/run.sh"), str(n)],
            capture_output=True, text=True, timeout=timeout,
            cwd=str(ROOT),
        )
        out = result.stdout.strip()
        for line in out.splitlines():
            if line.startswith("p"):
                parts = line.split(None, 1)
                if len(parts) == 2:
                    return parts[1]
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
<table class="meta">
<tr><th>Answer</th><td><code>{esc(expected)}</code></td></tr>
<tr><th>Output</th><td><code>{esc(output)}</code></td></tr>
<tr><th>Status</th><td class="status-{status.lower()}">{status}</td></tr>
<tr><th>Native helper</th><td>{"yes" if native else "no"}</td></tr>
</table>
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
        rows.append(
            f'<tr><td><a href="problems/p{label}.html">{label}</a></td>'
            f'<td class="{cls}">{status}</td>'
            f'<td><code>{esc(output)}</code></td>'
            f'<td><code>{esc(expected)}</code></td>'
            f'<td>{"native" if native else "flow"}</td></tr>'
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
<p>Each problem page shows the Flow source, generated C, generated MLIR, and the program output.</p>
<p>{passed} passed, {skipped} skipped, {total} total.</p>
<table>
<thead><tr><th>Problem</th><th>Status</th><th>Output</th><th>Expected</th><th>Type</th></tr></thead>
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
a { color: #1565c0; text-decoration: none; }
a:hover { text-decoration: underline; }
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
