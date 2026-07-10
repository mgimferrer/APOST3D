#!/usr/bin/env bash
# ==============================================================================
# run_tests.sh — APOST-3D gfortran regression test suite
#
# Runs all compiler-testset cases, saves outputs, and (optionally) compares
# them against reference outputs to flag numerical differences.
#
# Usage:
#   export APOST3D_PATH=/path/to/APOST3D
#   bash run_tests.sh [OPTIONS]
#
# Options:
#   --ref  <dir>   Directory containing reference .apost1 files for comparison
#                  (default: $APOST3D_PATH/compiler-testset if no separate ref)
#   --out  <dir>   Where to write output files (default: ./testrun_YYYYMMDD_HHMM)
#   --nthreads <N> OMP_NUM_THREADS (default: 1)
#   --skip <name>  Skip a test by name (can be repeated)
#   --help         Show this message
#
# Example:
#   bash run_tests.sh --ref /path/to/REFERENCE-OUTPUTS/compiler-testset \
#                     --out ./testresults --nthreads 4
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOST3D_PATH="${APOST3D_PATH:-$SCRIPT_DIR}"
APOST_BIN="$APOST3D_PATH/apost3d"
INPDIR="$APOST3D_PATH/compiler-testset"
OUTDIR=""
REFDIR=""
NTHREADS=1
SKIP=()

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref)    REFDIR="$2";    shift 2 ;;
    --out)    OUTDIR="$2";    shift 2 ;;
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --skip)   SKIP+=("$2");  shift 2 ;;
    --help)
      sed -n '2,30p' "$0" | sed 's/^# \{0,2\}//'
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Default output directory
if [[ -z "$OUTDIR" ]]; then
  OUTDIR="$APOST3D_PATH/testrun_$(date +%Y%m%d_%H%M)"
fi

# If no explicit ref dir, see if we have reference files in the input dir
if [[ -z "$REFDIR" ]]; then
  if ls "$INPDIR"/*.apost1 &>/dev/null 2>&1; then
    REFDIR="$INPDIR"
  fi
fi

# ------------------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------------------
if [[ ! -x "$APOST_BIN" ]]; then
  echo "ERROR: apost3d binary not found or not executable: $APOST_BIN"
  echo "       Build it first with: make all"
  exit 1
fi
if [[ ! -d "$INPDIR" ]]; then
  echo "ERROR: Input directory not found: $INPDIR"
  exit 1
fi

mkdir -p "$OUTDIR"

# ------------------------------------------------------------------------------
# Test definitions
# get_extra_files NAME  — prints space-separated list of extra files to copy
# (bash 3.2 compatible; no associative arrays)
# ------------------------------------------------------------------------------
get_extra_files() {
  case "$1" in
    CH3F)       echo "CH3F-OSLOs.fchk CH3F-OSLOs-preortho.fchk" ;;
    FeO4-2)     echo "FeO4-2-OSLOs.fchk" ;;
    O2-CASSCF)  echo "O2-CASSCF.dm1 O2-CASSCF.dm2" ;;
    *)          echo "" ;;
  esac
}

TESTS=("H2O-T-B3LYP" "CH3F" "FeCO2-PBEPBE" "FeO4-2")

# ------------------------------------------------------------------------------
# Helper: compare output vs reference
# Extracts lines with 6+ decimal places, ignores timing, reports differences
# ------------------------------------------------------------------------------
compare_output() {
  local name="$1"
  local new="$2"
  local ref="$3"

  if [[ ! -f "$ref" ]]; then
    echo "    [SKIP compare] reference not found: $ref"
    return
  fi

  local ndiff
  ndiff=$(diff \
    <(grep -E "[0-9]\.[0-9]{6,}" "$ref" | grep -iv "elapsed\| s)") \
    <(grep -E "[0-9]\.[0-9]{6,}" "$new" | grep -iv "elapsed\| s)") \
    | grep "^[<>]" | wc -l)

  if [[ "$ndiff" -eq 0 ]]; then
    echo "    [COMPARE] ✓  Numerical values identical to reference"
  else
    # Count purely formatting differences (same number, different width)
    local fmt_diff format_only=0
    while IFS= read -r line; do
      local val_ref val_new
      val_ref=$(echo "$line" | grep "^<" | grep -oE "[-]?[0-9]+\.[0-9]+" | head -1)
      val_new=$(echo "$line"  | grep "^>" | grep -oE "[-]?[0-9]+\.[0-9]+" | head -1)
      # If the difference is only in trailing digits, it's formatting/rounding
    done < <(diff \
      <(grep -E "[0-9]\.[0-9]{6,}" "$ref" | grep -iv "elapsed\| s)") \
      <(grep -E "[0-9]\.[0-9]{6,}" "$new" | grep -iv "elapsed\| s)"))

    echo "    [COMPARE] ⚠  $ndiff differing lines vs reference — see below:"
    diff \
      <(grep -E "[0-9]\.[0-9]{6,}" "$ref" | grep -iv "elapsed\| s)") \
      <(grep -E "[0-9]\.[0-9]{6,}" "$new" | grep -iv "elapsed\| s)") \
      | grep "^[<>]" | head -20
  fi
}

# ------------------------------------------------------------------------------
# Environment setup
# ------------------------------------------------------------------------------
ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS="$NTHREADS"

echo "============================================================"
echo "  APOST-3D regression test suite"
echo "============================================================"
echo "  Binary   : $APOST_BIN"
echo "  Inputs   : $INPDIR"
echo "  Outputs  : $OUTDIR"
echo "  Reference: ${REFDIR:-'(none, comparison disabled)'}"
echo "  Threads  : $NTHREADS"
echo "  Started  : $(date)"
echo "============================================================"
echo ""

# ------------------------------------------------------------------------------
# Run tests
# ------------------------------------------------------------------------------
PASS=0; FAIL=0; SKIPPED=0

for name in "${TESTS[@]}"; do
  # Check skip list
  skip=0
  for s in "${SKIP[@]:-}"; do
    [[ "$s" == "$name" ]] && skip=1 && break
  done
  if [[ "$skip" -eq 1 ]]; then
    echo "[ SKIP ] $name"
    SKIPPED=$((SKIPPED+1))
    continue
  fi

  # Check required input files exist
  if [[ ! -f "$INPDIR/$name.fchk" || ! -f "$INPDIR/$name.inp" ]]; then
    echo "[ SKIP ] $name  (input files not found in $INPDIR)"
    SKIPPED=$((SKIPPED+1))
    continue
  fi

  echo "------------------------------------------------------------"
  echo "[ RUN  ] $name"
  echo "         $(cat "$INPDIR/$name.inp" | grep -v "^#\|^$" | tr '\n' ' ' | xargs)"

  # Set up isolated run directory
  rundir="$OUTDIR/$name"
  mkdir -p "$rundir"

  # Copy required input files
  cp "$INPDIR/$name.fchk" "$rundir/"
  cp "$INPDIR/$name.inp"  "$rundir/"
  for extra in $(get_extra_files "$name"); do
    [[ -f "$INPDIR/$extra" ]] && cp "$INPDIR/$extra" "$rundir/"
  done

  outfile="$OUTDIR/$name.apost"
  t_start=$(date +%s)

  # Run
  set +e
  (
    cd "$rundir"
    ulimit -s unlimited 2>/dev/null || true
    "$APOST_BIN" "$name" > "$outfile" 2>&1
  )
  rc=$?
  set -e

  t_end=$(date +%s)
  elapsed=$((t_end - t_start))

  if [[ "$rc" -ne 0 ]]; then
    echo "    [FAIL ] exit code $rc  (${elapsed}s)"
    echo "    Last lines of output:"
    tail -5 "$outfile" | sed 's/^/            /'
    FAIL=$((FAIL+1))
  else
    termline=$(grep -c "Normal Termination" "$outfile" || true)
    if [[ "$termline" -gt 0 ]]; then
      echo "    [PASS ] Normal Termination  (${elapsed}s)  -> $outfile"
      PASS=$((PASS+1))
      # Compare vs reference if available
      if [[ -n "$REFDIR" ]]; then
        ref_file=""
        [[ -f "$REFDIR/$name.apost1" ]] && ref_file="$REFDIR/$name.apost1"
        [[ -n "$ref_file" ]] && compare_output "$name" "$outfile" "$ref_file"
      fi
    else
      echo "    [FAIL ] No 'Normal Termination' in output  (${elapsed}s)"
      tail -5 "$outfile" | sed 's/^/            /'
      FAIL=$((FAIL+1))
    fi
  fi
  echo ""
done

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------
echo "============================================================"
echo "  SUMMARY"
echo "============================================================"
echo "  PASS   : $PASS"
echo "  FAIL   : $FAIL"
echo "  SKIPPED: $SKIPPED"
echo "  Outputs: $OUTDIR/"
echo "  Finished: $(date)"
echo "============================================================"

[[ "$FAIL" -eq 0 ]]
