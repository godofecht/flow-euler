#!/usr/bin/env bash
# Run one or all Project Euler Flow solutions.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# Prefer the Flow repo binary (symlink in ~/.local/bin breaks SCRIPT_DIR).
if [[ -z "${FLOW_BIN:-}" ]]; then
  if [[ -x "${FLOW_REPO:-$HOME/flow}/flow" ]]; then
    FLOW_BIN="${FLOW_REPO:-$HOME/flow}/flow"
  else
    FLOW_BIN="$(command -v flow)"
  fi
fi
export FLOW_HOST="${FLOW_HOST:-python}"

extract_answer() {
  # Print the last answer line between the runner's separators.
  # Accepts integers, decimals, fractions, scientific notation, and alphanumeric strings.
  awk '
    /^----------------------------------------$/ { inblock=!inblock; next }
    NF {
      if ($0 ~ /^[-]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ ||
          $0 ~ /^[0-9]+\/[0-9]+$/ || $0 ~ /^[A-Za-z0-9.,+-]+$/) {
        fallback=$0
        if (inblock) ans=$0
      }
    }
    END { if (ans != "") print ans; else if (fallback != "") print fallback }
  '
}

# Problems with a companion under problems/native/pNNN.{c,cpp} are compiled
# as Flow→C plus that native TU (real solvers that need NTT / i128 / etc.).
run_with_native() {
  local n="$1" src="$2" native="$3"
  local build="$ROOT/.build"
  mkdir -p "$build"
  local c_out="$build/p${n}_flow.c"
  local exe="$build/p${n}"
  local py="${FLOW_REPO:-$HOME/flow}/src"
  PYTHONPATH="$py${PYTHONPATH:+:$PYTHONPATH}" \
    python3 -m flow.transpiler "$src" --c --lenient -o "$c_out" >/dev/null
  local cc=clang
  local brew_prefix="${HOMEBREW_PREFIX:-/opt/homebrew}"
  local libs=(-lm)
  # Optional high-precision / multiprecision for heavy solvers.
  if [[ -f "$brew_prefix/include/gmp.h" ]]; then
    libs+=(-L"$brew_prefix/lib" -I"$brew_prefix/include" -lgmp)
  fi
  if [[ -f "$brew_prefix/include/mpfr.h" ]]; then
    libs+=(-L"$brew_prefix/lib" -I"$brew_prefix/include" -lmpfr)
  fi
  if [[ "$native" == *.cpp ]]; then
    # Compile Flow C as C, then link with the C++ TU (preserves extern "C").
    local flow_o="$build/p${n}_flow.o"
    clang -O3 -march=native -c "$c_out" -o "$flow_o"
    clang++ -O3 -march=native -std=c++17 -pthread -I"$brew_prefix/include" -L"$brew_prefix/lib" \
      "$flow_o" "$native" -o "$exe" "${libs[@]}"
  else
    clang -O3 -march=native -I"$brew_prefix/include" -L"$brew_prefix/lib" \
      "$c_out" "$native" -o "$exe" "${libs[@]}"
  fi
  local out
  out=$("$exe")
  # Drop everything but the last numeric-looking line.
  printf '%s\n' "$out" | awk '
    NF {
      if ($0 ~ /^[-]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/) ans=$0
      else if ($0 ~ /^[0-9]+\/[0-9]+$/) ans=$0
      else if ($0 ~ /^[A-Za-z0-9.,+-]+$/) ans=$0
    }
    END { if (ans != "") print ans }
  '
}

# The Flow shell launcher uses Bash namerefs, which are unavailable in the
# Bash 3.2 shipped on macOS runners.  Compile through the Python transpiler
# directly so local runs and GitHub Actions use the same portable route.
run_with_flow_c() {
  local n="$1" src="$2"
  local build="$ROOT/.build"
  mkdir -p "$build"
  local c_out="$build/p${n}.c"
  local exe="$build/p${n}"
  local py="${FLOW_REPO:-$HOME/flow}/src"
  PYTHONPATH="$py${PYTHONPATH:+:$PYTHONPATH}" \
    python3 -m flow.transpiler "$src" --c --lenient -o "$c_out" >/dev/null
  clang -O3 -march=native "$c_out" -o "$exe" -lm
  "$exe"
}

run_one() {
  local n src
  # Force decimal: bash treats leading-zero args as octal.
  n=$(printf '%03d' "$((10#$1))")
  src="problems/p${n}.flow"
  if [[ ! -f "$src" ]]; then
    echo "missing $src" >&2
    return 1
  fi
  local out
  local native=""
  if [[ -f "problems/native/p${n}.c" ]]; then
    native="problems/native/p${n}.c"
  elif [[ -f "problems/native/p${n}.cpp" ]]; then
    native="problems/native/p${n}.cpp"
  fi
  if [[ -n "$native" ]]; then
    out=$(run_with_native "$n" "$src" "$native")
  else
    out=$(run_with_flow_c "$n" "$src" | extract_answer)
  fi
  echo "p${n}: ${out}"
}

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <n|all>" >&2
  exit 2
fi

if [[ "$1" == "all" ]]; then
  shopt -s nullglob
  for f in problems/p*.flow; do
    n=${f#problems/p}
    n=${n%.flow}
    # strip leading zeros for display arg, but keep numeric
    run_one "$((10#$n))"
  done
else
  run_one "$1"
fi
