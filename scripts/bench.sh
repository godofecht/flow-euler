#!/usr/bin/env bash
# Fair in-process bench: Flow-generated C vs hand-written C, same -O3 flags.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

FLOW_REPO="${FLOW_REPO:-$HOME/flow}"
FLOW_BIN="${FLOW_BIN:-$FLOW_REPO/flow}"
export FLOW_HOST="${FLOW_HOST:-python}"
ROUNDS="${ROUNDS:-200}"

mkdir -p bench/bin

echo "== transpile Flow bench/flow_p010.flow =="
FLOW_HOST="$FLOW_HOST" "$FLOW_BIN" compile bench/flow_p010.flow >/dev/null || true
GEN="$FLOW_REPO/build/flow_p010.c"
[[ -f "$GEN" ]] || { echo "missing $GEN" >&2; exit 1; }

# Ensure time.h is present for clock()
if ! grep -q 'time.h' "$GEN"; then
  tmp=$(mktemp)
  awk 'NR==1{print; print "#include <time.h>"; next}1' "$GEN" > "$tmp"
  mv "$tmp" "$GEN"
fi

echo "== clang -O3 -march=native (Flow-generated C) =="
clang -O3 -march=native -o bench/bin/flow_p010 "$GEN"

echo "== clang -O3 -march=native (hand-written C) =="
clang -O3 -march=native -o bench/bin/c_p010 bench/c_p010.c

echo "== run ($ROUNDS in-process rounds each) =="
echo "--- Flow ---"
./bench/bin/flow_p010
echo "--- C ---"
./bench/bin/c_p010 "$ROUNDS"
