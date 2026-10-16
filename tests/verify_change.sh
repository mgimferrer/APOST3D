#!/usr/bin/env bash
#
# verify_change.sh -- check that a change to the sources did not change
# the results, by running every test with a build of a base commit and a
# build of the working tree, both on one thread, and comparing outputs.
#
#   bash tests/verify_change.sh [--base REV] [--format] [--threads N]
#                               [--filter STR] [--keep]
#
#   --base REV    commit to compare against (default: HEAD). It is built
#                 in a separate git worktree; the working tree is not
#                 touched apart from its usual incremental build.
#   --format      output-format change: only the numbers must be
#                 identical (compare_outputs.py with --ulps 0, and
#                 compare_numbers.py: every number of each output, in
#                 order); changed layout or labels are listed as notes,
#                 and a change in the numbers fails. Without it, the
#                 outputs must be identical line by line (refactors,
#                 parallel regions), timings and thread counts aside.
#                 The results files (<test>.apost.json) are compared as
#                 well: identical values, or within compare_json.py's
#                 tolerance with --format.
#   --threads N   also run the working-tree build on N threads and show
#                 its differences from the one-thread run (summing in
#                 another order may change the last printed digit).
#   --filter STR  only the tests whose name contains STR.
#   --keep        keep the base worktree (removed by default).
#
# Exit code 0 when the comparison passes, 1 otherwise. The outputs stay
# in the directory printed at the end.

set -euo pipefail

BASE="HEAD"
FORMAT=0
THREADS=""
FILTER=""
KEEP=0
while [ $# -gt 0 ]; do
  case "$1" in
    --base)    BASE="$2"; shift 2 ;;
    --format)  FORMAT=1; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --filter)  FILTER="$2"; shift 2 ;;
    --keep)    KEEP=1; shift ;;
    -h|--help) sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "verify_change.sh: unknown option $1" >&2; exit 2 ;;
  esac
done

REPO="$(git -C "$(dirname "$0")" rev-parse --show-toplevel)"
BASE_SHA="$(git -C "$REPO" rev-parse --verify "$BASE^{commit}")"
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/verify_change.XXXXXX")"
WT="$WORKDIR/base-src"

cleanup() {
  if [ "$KEEP" -eq 0 ] && [ -d "$WT" ]; then
    git -C "$REPO" worktree remove --force "$WT" >/dev/null 2>&1 || true
  fi
}
trap cleanup EXIT

## the base build uses the libxc of the current build (the worktree has
## none of its own, and building one would take minutes)
if [ -z "${LIBXC_DIR:-}" ]; then
  if [ -f "$REPO/objects/libxc.mk" ]; then
    LIBXC_DIR="$(sed -n 's/^LIBXC_INC := -I\(.*\)\/include$/\1/p' "$REPO/objects/libxc.mk")"
  fi
  if [ -z "${LIBXC_DIR:-}" ]; then
    LIBXC_DIR="$REPO/libxc-$(sed -n 's/^LIBXC_VERSION="\(.*\)"/\1/p' "$REPO/compile_libxc.sh")"
  fi
fi
export LIBXC_DIR

run_suite() {   # binary, nthreads, output dir
  local args=(--no-full --no-color --nthreads "$2" --binary "$1" --output-dir "$3")
  [ -n "$FILTER" ] && args+=(--filter "$FILTER")
  python3 "$REPO/tests/run_tests.py" "${args[@]}" > "$3.log" 2>&1 || true
}

echo "Base commit : $(git -C "$REPO" log -1 --format='%h %s' "$BASE_SHA")"
echo "Compared to : working tree of $REPO"
echo "libxc       : $LIBXC_DIR"
echo "Outputs in  : $WORKDIR"
echo

echo "1. Building the base commit (worktree) ..."
git -C "$REPO" worktree add --detach "$WT" "$BASE_SHA" >/dev/null 2>&1
if ! (cd "$WT" && APOST3D_PATH="$WT" bash make_compile.sh > "$WORKDIR/build-base.log" 2>&1) \
   || [ ! -x "$WT/apost3d" ]; then
  echo "   base build failed, see $WORKDIR/build-base.log"; exit 1
fi

echo "2. Building the working tree ..."
if ! (cd "$REPO" && APOST3D_PATH="$REPO" bash make_compile.sh > "$WORKDIR/build-new.log" 2>&1); then
  echo "   working-tree build failed, see $WORKDIR/build-new.log"; exit 1
fi

echo "3. Running the tests on one thread with both builds ..."
run_suite "$WT/apost3d"   1 "$WORKDIR/base"
run_suite "$REPO/apost3d" 1 "$WORKDIR/new"

## lines that change from run to run without being results
NOISE='TIMING (CPU|WALL)|Elapsed time|distributed over|threads out of|OMP_NUM_THREADS|OpenMP threads'
filtered() { grep -Ev "$NOISE" "$1" || true; }

echo
FAIL=0
if [ "$FORMAT" -eq 0 ]; then
  echo "Base vs working tree, one thread: outputs must be identical"
  for f in "$WORKDIR"/base/*.apost; do
    name="$(basename "$f")"
    if [ ! -f "$WORKDIR/new/$name" ]; then
      echo "  MISSING     $name"; FAIL=1
    elif diff -q <(filtered "$f") <(filtered "$WORKDIR/new/$name") >/dev/null; then
      echo "  identical   $name"
    else
      echo "  DIFFERS     $name"
      diff <(filtered "$f") <(filtered "$WORKDIR/new/$name") | head -40 || true
      FAIL=1
    fi
  done
else
  echo "Base vs working tree, one thread: every number must be identical"
  python3 "$REPO/tests/compare_outputs.py" --ref-dir "$WORKDIR/base" \
    --out-dir "$WORKDIR/new" --ulps 0 --abs-tol 0 || FAIL=1
  echo
  echo "Every number of each output, in order (any layout):"
  python3 "$REPO/tests/compare_numbers.py" --old-dir "$WORKDIR/base" \
    --new-dir "$WORKDIR/new" | grep -v '^  same numbers' || true
  python3 "$REPO/tests/compare_numbers.py" --old-dir "$WORKDIR/base" \
    --new-dir "$WORKDIR/new" > /dev/null || FAIL=1
fi

## results files (<test>.apost.json), when the base wrote them too
if ls "$WORKDIR"/base/*.apost.json >/dev/null 2>&1; then
  echo
  if [ "$FORMAT" -eq 0 ]; then
    echo "Results files (.apost.json): values must be identical"
    python3 "$REPO/tests/compare_json.py" --ref-dir "$WORKDIR/base" \
      --out-dir "$WORKDIR/new" --rtol 0 | grep -v '^  same' || true
    python3 "$REPO/tests/compare_json.py" --ref-dir "$WORKDIR/base" \
      --out-dir "$WORKDIR/new" --rtol 0 > /dev/null || FAIL=1
  else
    echo "Results files (.apost.json): values within compare_json's tolerance"
    python3 "$REPO/tests/compare_json.py" --ref-dir "$WORKDIR/base" \
      --out-dir "$WORKDIR/new" | grep -v '^  same' || true
    python3 "$REPO/tests/compare_json.py" --ref-dir "$WORKDIR/base" \
      --out-dir "$WORKDIR/new" > /dev/null || FAIL=1
  fi
fi

if [ -n "$THREADS" ]; then
  echo
  echo "4. Working tree on $THREADS threads vs one thread (last-digit drift is expected)"
  run_suite "$REPO/apost3d" "$THREADS" "$WORKDIR/new-n$THREADS"
  for f in "$WORKDIR"/new/*.apost; do
    name="$(basename "$f")"
    if diff -q <(filtered "$f") <(filtered "$WORKDIR/new-n$THREADS/$name") >/dev/null 2>&1; then
      echo "  identical   $name"
    else
      echo "  differs     $name"
      diff <(filtered "$f") <(filtered "$WORKDIR/new-n$THREADS/$name") | head -20 || true
    fi
  done
fi

echo
if [ "$FAIL" -eq 0 ]; then
  echo "PASS"
else
  echo "FAIL: the change altered the results (see above)"
fi
exit "$FAIL"
