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
echo "passed=$pass skipped=$skip failed=$fail"
[[ "$fail" -eq 0 ]]
