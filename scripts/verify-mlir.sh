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

failed_labels=()
failed_actual=()
failed_expected=()
failed_reasons=()
mlir_failed_labels=()
mlir_failed_reasons=()

build="$ROOT/.build/mlir"
mkdir -p "$build"

gha_escape() {
  local value="$1"
  value="${value//'%'/'%25'}"
  value="${value//$'\r'/'%0D'}"
  value="${value//$'\n'/'%0A'}"
  printf '%s' "$value"
}

record_failure() {
  local label="$1"
  local actual="$2"
  local expected="$3"
  local reason="$4"

  failed_labels+=("$label")
  failed_actual+=("$actual")
  failed_expected+=("$expected")
  failed_reasons+=("$reason")
  fail=$((fail + 1))

  echo "FAIL p${label}: ${reason}; got '${actual}' expected '${expected}'" >&2
  if [[ "${GITHUB_ACTIONS:-}" == "true" ]]; then
    echo "::error title=Euler p${label} ${reason}::$(gha_escape "got '${actual}', expected '${expected}'")"
  fi
}

record_mlir_failure() {
  local label="$1"
  local reason="$2"

  mlir_failed_labels+=("$label")
  mlir_failed_reasons+=("$reason")
  mlir_fail=$((mlir_fail + 1))
  echo "MLIR_FAIL p${label}: ${reason}" >&2
}

print_failure_summary() {
  if [[ "$fail" -gt 0 ]]; then
    echo
    echo "FAILED PROBLEMS ($fail)"
    echo "--------------------"
    for i in "${!failed_labels[@]}"; do
      echo "p${failed_labels[$i]}: ${failed_reasons[$i]}; got '${failed_actual[$i]}' expected '${failed_expected[$i]}'"
    done
  fi

  if [[ "$mlir_fail" -gt 0 ]]; then
    echo
    echo "MLIR BACKEND LIMITATIONS ($mlir_fail)"
    echo "--------------------------------"
    for i in "${!mlir_failed_labels[@]}"; do
      echo "p${mlir_failed_labels[$i]}: ${mlir_failed_reasons[$i]}"
    done
  fi

  if [[ -n "${GITHUB_STEP_SUMMARY:-}" ]]; then
    {
      if [[ "$fail" -gt 0 ]]; then
        echo "## MLIR answer mismatches"
        echo
        echo "| Problem | Failure | Expected | Actual |"
        echo "| --- | --- | --- | --- |"
        for i in "${!failed_labels[@]}"; do
          echo "| \`p${failed_labels[$i]}\` | ${failed_reasons[$i]} | \`${failed_expected[$i]}\` | \`${failed_actual[$i]}\` |"
        done
        echo
      fi

      if [[ "$mlir_fail" -gt 0 ]]; then
        echo "## MLIR backend limitations"
        echo
        echo "| Problem | Stage |"
        echo "| --- | --- |"
        for i in "${!mlir_failed_labels[@]}"; do
          echo "| \`p${mlir_failed_labels[$i]}\` | ${mlir_failed_reasons[$i]} |"
        done
      fi
    } >> "$GITHUB_STEP_SUMMARY"
  fi
}

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

  if out=$(run_mlir "$numeric" "$src" 2>/dev/null); then
    rc=0
  else
    rc=$?
  fi

  if [[ "$rc" -ne 0 || -z "$out" ]]; then
    case "$rc" in
      2) reason="Flow to MLIR transpile failure" ;;
      3) reason="MLIR/LLVM toolchain failure" ;;
      124) reason="execution timeout" ;;
      0) reason="no answer produced" ;;
      *) reason="backend exited with status $rc" ;;
    esac
    record_mlir_failure "$label" "$reason"
    continue
  fi

  if [[ "$out" == "$expected" ]]; then
    echo "OK  p${label} = $out"
    pass=$((pass + 1))
  else
    record_failure "$label" "$out" "$expected" "answer mismatch"
  fi
done < answers.txt

print_failure_summary

echo "----"
echo "passed=$pass failed=$fail mlir_fail=$mlir_fail skipped=$skip"
# MLIR failures are expected for some problems due to compiler limitations.
# Only hard answer mismatches (fail) cause a non-zero exit.
[[ "$fail" -eq 0 ]]
