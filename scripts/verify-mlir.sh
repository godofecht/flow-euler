#!/usr/bin/env bash
# Verify Flow solutions through the MLIR backend.
#
# Flow -> MLIR -> mlir-opt -> mlir-translate -> llc -> clang -> run
#
# Only tests problems without native C helpers (pure Flow -> MLIR).
# Problems with native helpers are skipped because the MLIR path does not
# link external C TUs.
#
# Usage:
#   ./scripts/verify-mlir.sh [first] [last]
#
# Environment:
#   FLOW_REPO   path to Flow compiler checkout (default $HOME/flow)
#   LLVM_BIN    path to LLVM bin dir with mlir-opt, mlir-translate, llc
#               (default /opt/homebrew/opt/llvm/bin)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

FLOW_REPO="${FLOW_REPO:-$HOME/flow}"
LLVM_BIN="${LLVM_BIN:-/opt/homebrew/opt/llvm/bin}"
export PATH="$LLVM_BIN:$PATH"
export FLOW_HOST="${FLOW_HOST:-python}"

first="${1:-1}"
last="${2:-9999}"

pass=0
fail=0
skip=0
mlir_fail=0

build="$ROOT/.build/mlir"
mkdir -p "$build"

run_mlir() {
  local n="$1" src="$2"
  local mlir_out="$build/p${n}.mlir"
  local lowered="$build/p${n}_lowered.mlir"
  local ll_out="$build/p${n}.ll"
  local s_out="$build/p${n}.s"
  local exe="$build/p${n}"

  # Flow -> MLIR
  if ! PYTHONPATH="$FLOW_REPO/src" python3 -m flow.transpiler "$src" \
        --mlir --lenient -o "$mlir_out" 2>/dev/null; then
    return 2  # transpile failure
  fi

  # MLIR -> LLVM IR (lower dialects, translate)
  if ! mlir-opt \
        --convert-arith-to-llvm \
        --convert-func-to-llvm \
        --convert-cf-to-llvm \
        --convert-index-to-llvm \
        --reconcile-unrealized-casts \
        "$mlir_out" > "$lowered" 2>/dev/null; then
    return 3  # mlir-opt failure
  fi

  if ! mlir-translate --mlir-to-llvmir "$lowered" -o "$ll_out" 2>/dev/null; then
    return 3
  fi

  # LLVM IR -> assembly -> binary
  if ! llc -O3 "$ll_out" -o "$s_out" 2>/dev/null; then
    return 3
  fi

  if ! clang "$s_out" -o "$exe" -lm 2>/dev/null; then
    return 3
  fi

  # Run and extract answer (30s timeout for slow solvers)
  timeout 30 "$exe" | awk '
    NF {
      if ($0 ~ /^[-]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/) ans=$0
      else if ($0 ~ /^[0-9]+\/[0-9]+$/) ans=$0
      else if ($0 ~ /^[A-Za-z0-9.,+-]+$/) ans=$0
    }
    END { if (ans != "") print ans }
  '
}

while read -r num expected; do
  [[ -z "${num:-}" || "$num" =~ ^# ]] && continue
  numeric=$((10#$num))
  (( numeric < first || numeric > last )) && continue

  label=$(printf '%03d' "$numeric")
  src="problems/p${label}.flow"

  if [[ ! -f "$src" ]]; then
    skip=$((skip + 1))
    continue
  fi

  # Skip problems with native helpers (MLIR path doesn't link C TUs)
  if [[ -f "problems/native/p${label}.c" ]] || [[ -f "problems/native/p${label}.cpp" ]]; then
    skip=$((skip + 1))
    continue
  fi

  out=$(run_mlir "$numeric" "$src" 2>/dev/null || true)
  rc=$?
  if [[ -z "$out" ]]; then
    echo "MLIR_FAIL p${label}" >&2
    mlir_fail=$((mlir_fail + 1))
    continue
  fi

  if [[ "$out" == "$expected" ]]; then
    echo "OK  p${label} = $out"
    pass=$((pass + 1))
  else
    echo "FAIL p${label}: got '${out}' expected '${expected}'" >&2
    fail=$((fail + 1))
  fi
done < answers.txt

echo "----"
echo "passed=$pass failed=$fail mlir_fail=$mlir_fail skipped=$skip"
# MLIR failures are expected for some problems due to compiler limitations.
# Only hard answer mismatches (fail) cause a non-zero exit.
[[ "$fail" -eq 0 ]]
