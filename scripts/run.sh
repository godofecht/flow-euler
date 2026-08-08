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
    inblock && NF {
      if ($0 ~ /^[-]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/) ans=$0
      else if ($0 ~ /^[0-9]+\/[0-9]+$/) ans=$0
      else if ($0 ~ /^[A-Za-z0-9.,+-]+$/) ans=$0
    }
    END { if (ans != "") print ans }
  '
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
  out=$(FLOW_HOST="$FLOW_HOST" "$FLOW_BIN" run "$src" 2>/dev/null | extract_answer)
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
