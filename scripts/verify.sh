#!/usr/bin/env bash
# Verify all solved problems against answers.txt
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

pass=0
fail=0
while read -r num expected; do
  [[ -z "${num:-}" || "$num" =~ ^# ]] && continue
  out=$("./scripts/run.sh" "$((10#$num))" | awk '{print $2}')
  label=$(printf '%03d' "$((10#$num))")
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
