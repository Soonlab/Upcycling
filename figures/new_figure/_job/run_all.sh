#!/usr/bin/env bash
# Rebuild every figure, then run the saved-file verifier and the hard-coding scan.
# Detached use:  setsid nohup bash _job/run_all.sh > _job/run_all.log 2>&1 < /dev/null &
set -u
cd "$(dirname "$0")/.."
PY=/home/soon/miniconda3/envs/dram_env/bin/python
fail=0
for f in build_fig*.py build_supS*.py; do
  echo "=== $f"
  if ! $PY "$f"; then echo "!!! $f FAILED"; fail=$((fail+1)); fi
done
echo "=== verify_outputs"; $PY _job/verify_outputs.py || fail=$((fail+1))
echo "=== check_no_hardcoded"; $PY _job/check_no_hardcoded.py
echo "=== done, $fail failure(s)"
exit $fail
