#!/usr/bin/env bash
# Rebuild the eleven figure pages into figures_v2/ next to this script.
set -u
cd "$(dirname "$0")"
PY=/home/soon/miniconda3/envs/dram_env/bin/python
fail=0
for f in build_v2_fig*.py build_v2_supS*.py build_v2_graphical_abstract.py; do
  echo "=== $f"
  if ! $PY "$f"; then echo "!!! $f FAILED"; fail=$((fail+1)); fi
done
echo "=== done, $fail failure(s)"
exit $fail
