#!/usr/bin/env bash
#
# verify_omp_change.sh -- three-way correctness check for an uncommitted
# OpenMP change, run from the APOST-3D repo root:
#
#   bash tests/verify_omp_change.sh [NTHREADS]
#
# NTHREADS defaults to 4 -- pass the core count you want to stress-test
# (e.g. `bash tests/verify_omp_change.sh 8`).
#
# What it does:
#   1. BASELINE   = git HEAD (your uncommitted changes stashed away),
#                   built and run at NTHREADS=1.
#   2. FIXED N=1  = your current working tree (uncommitted changes back
#                   in place), built and run at NTHREADS=1.
#   3. FIXED N=K  = the same fixed build, run at NTHREADS=K.
#
# BASELINE vs FIXED N=1 should be BIT-IDENTICAL. A REDUCTION/COLLAPSE
# clause running on a single thread collapses to the exact same
# accumulation order as the original serial code -- any difference here
# means the change actually altered something, not just reordered it.
#
# FIXED N=1 vs FIXED N=K may show tiny floating-point-reordering
# differences (last printed digit or two, from summing the same terms in
# a different order across threads) -- that's expected and fine. A
# difference bigger than that, or a different atom/energy value entirely,
# is not.
#
# Your working tree is stashed only for step 1 and restored immediately
# after -- a trap guarantees the stash is popped back even if the script
# fails partway through.

set -euo pipefail

NTHREADS="${1:-4}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

WORKDIR="$(mktemp -d)"
BASE_OUT="$WORKDIR/baseline"
FIXED1_OUT="$WORKDIR/fixed_n1"
FIXEDN_OUT="$WORKDIR/fixed_n${NTHREADS}"

STASHED=0
restore_stash() {
  if [ "$STASHED" -eq 1 ]; then
    echo "(restoring your working-tree changes...)"
    git stash pop -q
    STASHED=0
  fi
}
trap restore_stash EXIT

# Lines that legitimately vary run-to-run / thread-count-to-thread-count
# (timing prints, thread-count announcements) -- not numerical results,
# so ignore them when diffing. TIMING CPU/WALL is the current format
# (timing_mod, added 2026-08-13); the old "Elapsed time" format it
# replaced is kept here too in case any output still uses it.
NOISE_PATTERN='TIMING (CPU|WALL)|Elapsed time|distributed over|threads out of'

filter() { grep -Ev "$NOISE_PATTERN" "$1"; }

echo "=== 1/3: BASELINE (git HEAD, NTHREADS=1) ==="
if ! git diff --quiet -- sources/ || ! git diff --cached --quiet -- sources/; then
  git stash push -q -m "verify_omp_change: temporary stash"
  STASHED=1
fi
bash make_compile.sh clean >/dev/null
bash make_compile.sh >/dev/null
python3 tests/run_tests.py --nthreads 1 --keep-output="$BASE_OUT"

restore_stash

echo
echo "=== 2/3: FIXED (working tree, NTHREADS=1) ==="
bash make_compile.sh clean >/dev/null
bash make_compile.sh >/dev/null
python3 tests/run_tests.py --nthreads 1 --keep-output="$FIXED1_OUT"

echo
echo "=== 3/3: FIXED at NTHREADS=$NTHREADS ==="
python3 tests/run_tests.py --nthreads "$NTHREADS" --keep-output="$FIXEDN_OUT"

echo
echo "════════════════════════════════════════════════════════════"
echo "  BASELINE vs FIXED, both NTHREADS=1 -- expect IDENTICAL"
echo "════════════════════════════════════════════════════════════"
FAIL1=0
for f in "$BASE_OUT"/*.apost; do
  name="$(basename "$f")"
  if diff -q <(filter "$f") <(filter "$FIXED1_OUT/$name") >/dev/null 2>&1; then
    echo "  identical   $name"
  else
    echo "  DIFFERS     $name"
    diff <(filter "$f") <(filter "$FIXED1_OUT/$name") || true
    FAIL1=1
  fi
done

echo
echo "════════════════════════════════════════════════════════════"
echo "  FIXED NTHREADS=1 vs FIXED NTHREADS=$NTHREADS -- small drift OK"
echo "════════════════════════════════════════════════════════════"
for f in "$FIXED1_OUT"/*.apost; do
  name="$(basename "$f")"
  if diff -q <(filter "$f") <(filter "$FIXEDN_OUT/$name") >/dev/null 2>&1; then
    echo "  identical   $name"
  else
    echo "  differs     $name  (inspect below -- should be last printed digit only)"
    diff <(filter "$f") <(filter "$FIXEDN_OUT/$name") || true
  fi
done

echo
if [ "$FAIL1" -eq 0 ]; then
  echo "PASS: baseline and fixed-at-1-thread outputs are bit-identical."
else
  echo "FAIL: baseline and fixed-at-1-thread outputs differ -- see above."
  echo "      This means the change altered results, not just their order."
fi
echo
echo "Raw .apost outputs kept in: $WORKDIR"
