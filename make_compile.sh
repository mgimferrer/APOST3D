#!/usr/bin/env bash
# ==============================================================================
# make_compile.sh — compile APOST-3D with GCC/gfortran
#
# Replaces the old Intel ifort/PGO two-step build with a single-
# step portable gfortran build. No Profile-Guided Optimisation required.
#
# Usage:
#   export APOST3D_PATH=/path/to/APOST3D   # or edit DEFAULT below
#   bash make_compile.sh [--nthreads N] [--clean] [--help]
#
# Options:
#   --nthreads <N>  OMP_NUM_THREADS for test runs (default: all logical CPUs)
#   --clean         Run 'make clean' before building
#   --help          Show this message
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOST3D_PATH="${APOST3D_PATH:-$SCRIPT_DIR}"
CLEAN=0
NTHREADS=""

# ------------------------------------------------------------------------------
# Argument parsing
# ------------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --nthreads) NTHREADS="$2"; shift 2 ;;
    --clean)    CLEAN=1; shift ;;
    --help)
      sed -n '2,15p' "$0" | sed 's/^# \{0,2\}//'
      exit 0
      ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# ------------------------------------------------------------------------------
# Auto-detect logical CPU count
# ------------------------------------------------------------------------------
if [[ -z "$NTHREADS" ]]; then
  NTHREADS=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 1)
fi

# ------------------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------------------
MAKEFILE="$APOST3D_PATH/Makefile"
LIBXC_A="$APOST3D_PATH/libxc-4.2.3/lib/libxc.a"

if [[ ! -f "$MAKEFILE" ]]; then
  echo "ERROR: Makefile not found at $APOST3D_PATH"
  echo "       Make sure APOST3D_PATH is set correctly."
  exit 1
fi

if [[ ! -f "$LIBXC_A" ]]; then
  echo "ERROR: libxc not built yet. Run first:"
  echo "       bash $APOST3D_PATH/compile_libxc.sh"
  exit 1
fi

# ------------------------------------------------------------------------------
# gfortran preflight: must exist and be >= 10 (required for
# -fallow-argument-mismatch, used throughout the Makefile).
# ------------------------------------------------------------------------------
if ! command -v gfortran &>/dev/null; then
  echo "ERROR: gfortran not found on PATH."
  echo "       Install with:"
  echo "         macOS:  brew install gcc"
  echo "         Ubuntu: sudo apt install gfortran"
  exit 1
fi

GFORTRAN_VERSION_FULL="$(gfortran --version | head -1)"
GFORTRAN_VERSION_MAJOR="$(gfortran -dumpversion | cut -d. -f1)"

if [[ ! "$GFORTRAN_VERSION_MAJOR" =~ ^[0-9]+$ ]] || (( GFORTRAN_VERSION_MAJOR < 10 )); then
  echo "ERROR: gfortran >= 10 is required (found: $GFORTRAN_VERSION_FULL)."
  echo "       -fallow-argument-mismatch was introduced in GCC 10."
  echo "       Install a newer gfortran, e.g.:"
  echo "         macOS:  brew install gcc"
  echo "         Ubuntu: sudo apt install gfortran-12"
  exit 1
fi

# ------------------------------------------------------------------------------
# Environment
# ------------------------------------------------------------------------------
export APOST3D_PATH
export OMP_NUM_THREADS="$NTHREADS"
ulimit -s unlimited 2>/dev/null || true

mkdir -p "$APOST3D_PATH/objects"

echo "============================================================"
echo "  APOST-3D gfortran build"
echo "============================================================"
echo "  APOST3D_PATH : $APOST3D_PATH"
echo "  Makefile     : $MAKEFILE"
echo "  Compiler     : $GFORTRAN_VERSION_FULL"
echo "  OMP threads  : $OMP_NUM_THREADS"
echo "  Started      : $(date)"
echo "============================================================"
echo ""

# ------------------------------------------------------------------------------
# Compiler-identity drift check.
#
# gfortran .mod files are NOT portable across compiler versions ("Cannot read
# module file ... created by a different version of GNU Fortran"). The
# Makefile's file-timestamp dependencies can't detect "same source, different
# compiler" — only a change in the *compiler itself* triggers this. So we
# stamp the compiler identity used for the last build and compare it here;
# if it changed (new gfortran version, switched machines, etc.) we force a
# clean before rebuilding instead of letting a cryptic module-version error
# surface mid-build.
# ------------------------------------------------------------------------------
COMPILER_STAMP="$APOST3D_PATH/objects/.gfortran_version"

if [[ -f "$COMPILER_STAMP" ]]; then
  PREV_VERSION="$(cat "$COMPILER_STAMP")"
  if [[ "$PREV_VERSION" != "$GFORTRAN_VERSION_FULL" ]] && [[ "$CLEAN" -ne 1 ]]; then
    echo "--- Compiler change detected ---"
    echo "  Previous build : $PREV_VERSION"
    echo "  Current        : $GFORTRAN_VERSION_FULL"
    echo "  Forcing 'make clean' to avoid stale/incompatible .mod files."
    echo ""
    CLEAN=1
  fi
fi

# ------------------------------------------------------------------------------
# Optional clean
# ------------------------------------------------------------------------------
if [[ "$CLEAN" -eq 1 ]]; then
  echo "--- make clean ---"
  make -f "$MAKEFILE" -C "$APOST3D_PATH" clean
  echo ""
fi

# Record the compiler identity used for this build (written after clean so a
# failed/interrupted build doesn't falsely mark the stamp as up to date).
echo "$GFORTRAN_VERSION_FULL" > "$COMPILER_STAMP"

# ------------------------------------------------------------------------------
# Build main binary + standalone EOS + utility
# ------------------------------------------------------------------------------
echo "--- Building apost3d, apost3d-eos, eos_aom ---"
make -f "$MAKEFILE" -C "$APOST3D_PATH" all
echo ""

# ------------------------------------------------------------------------------
# Ad-hoc code signing (required on macOS 26+ Tahoe beta for GCC-compiled bins)
# codesign is a no-op on Linux, so this is safe cross-platform.
# ------------------------------------------------------------------------------
if command -v codesign &>/dev/null; then
  echo "--- Ad-hoc code signing (macOS) ---"
  for bin in apost3d apost3d-eos eos_aom; do
    if [[ -x "$APOST3D_PATH/$bin" ]]; then
      codesign --force --sign - "$APOST3D_PATH/$bin" && \
        echo "  signed: $bin" || \
        echo "  warning: codesign failed for $bin (non-fatal)"
    fi
  done
  echo ""
fi

# ------------------------------------------------------------------------------
# Verify binaries
# ------------------------------------------------------------------------------
echo "--- Binaries produced ---"
for bin in apost3d apost3d-eos eos_aom; do
  if [[ -x "$APOST3D_PATH/$bin" ]]; then
    echo "  ✓  $APOST3D_PATH/$bin  ($(du -sh "$APOST3D_PATH/$bin" | cut -f1))"
  else
    echo "  ✗  $APOST3D_PATH/$bin  NOT FOUND — build may have failed"
  fi
done
echo ""

echo "============================================================"
echo "  Build complete: $(date)"
echo "============================================================"
echo ""
echo "To run a calculation:"
echo "  ulimit -s unlimited"
echo "  export OMP_NUM_THREADS=$NTHREADS"
echo "  cd /path/to/input/files"
echo "  $APOST3D_PATH/apost3d jobname > jobname.apost 2>&1"
echo ""
echo "To run the full regression test suite:"
echo "  bash $APOST3D_PATH/run_tests.sh \\"
echo "    --ref $APOST3D_PATH/compiler-testset \\"
echo "    --out $APOST3D_PATH/testrun_\$(date +%Y%m%d) \\"
echo "    --nthreads $NTHREADS"
