#!/usr/bin/env bash
# Verify all solved problems against answers.txt
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

pass=0
fail=0
first="${1:-1}"
last="${2:-9999}"
while read -r num expected; do
  [[ -z "${num:-}" || "$num" =~ ^# ]] && continue
  numeric=$((10#$num))
  (( numeric < first || numeric > last )) && continue
  out=$("./scripts/run.sh" "$numeric" | awk '{print $2}')
  label=$(printf '%03d' "$numeric")
  if [[ "$out" == "$expected" ]]; then
    echo "OK  p${label} = $out"
    pass=$((pass + 1))
  else
    echo "FAIL p${label}: got '${out}' expected '${expected}'" >&2
    fail=$((fail + 1))
  fi
done < answers.txt

echo "----"
echo "passed=$pass failed=$fail"
[[ "$fail" -eq 0 ]]
