#!/usr/bin/env bash
# Verify all solved problems against answers.txt.
# Problems listed in .ci-skip.txt are skipped (too slow, crash, or broken).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

pass=0
fail=0
skip=0
first="${1:-1}"
last="${2:-9999}"

failed_labels=()
failed_actual=()
failed_expected=()
failed_reasons=()

# Build a space-separated skip list for grep-style matching.
skip_nums=""
if [[ -f "$ROOT/.ci-skip.txt" ]]; then
  while read -r snum sreason; do
    [[ -z "${snum:-}" || "$snum" =~ ^# ]] && continue
    skip_nums="$skip_nums $snum"
  done < "$ROOT/.ci-skip.txt"
fi

is_skipped() {
  local n="$1"
  for s in $skip_nums; do
    if [[ "$((10#$s))" == "$((10#$n))" ]]; then return 0; fi
  done
  return 1
}

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

print_failure_summary() {
  if [[ "$fail" -eq 0 ]]; then
    return
  fi

  echo
  echo "FAILED PROBLEMS ($fail)"
  echo "--------------------"
  for i in "${!failed_labels[@]}"; do
    echo "p${failed_labels[$i]}: ${failed_reasons[$i]}; got '${failed_actual[$i]}' expected '${failed_expected[$i]}'"
  done

  if [[ -n "${GITHUB_STEP_SUMMARY:-}" ]]; then
    {
      echo "## Euler verification failures"
      echo
      echo "| Problem | Failure | Expected | Actual |"
      echo "| --- | --- | --- | --- |"
      for i in "${!failed_labels[@]}"; do
        echo "| \`p${failed_labels[$i]}\` | ${failed_reasons[$i]} | \`${failed_expected[$i]}\` | \`${failed_actual[$i]}\` |"
      done
    } >> "$GITHUB_STEP_SUMMARY"
  fi
}

while read -r num expected; do
  [[ -z "${num:-}" || "$num" =~ ^# ]] && continue
  numeric=$((10#$num))
  (( numeric < first || numeric > last )) && continue
  if is_skipped "$numeric"; then
    label=$(printf '%03d' "$numeric")
    echo "SKIP p${label}"
    skip=$((skip + 1))
    continue
  fi

  label=$(printf '%03d' "$numeric")
  if run_output=$("./scripts/run.sh" "$numeric"); then
    out=$(printf '%s\n' "$run_output" | awk '{print $2}')
    if [[ "$out" == "$expected" ]]; then
      echo "OK  p${label} = $out"
      pass=$((pass + 1))
    else
      record_failure "$label" "$out" "$expected" "answer mismatch"
    fi
  else
    rc=$?
    record_failure "$label" "run.sh exited with status $rc" "$expected" "execution failure"
  fi
done < answers.txt

print_failure_summary

echo "----"
echo "passed=$pass skipped=$skip failed=$fail"
[[ "$fail" -eq 0 ]]
