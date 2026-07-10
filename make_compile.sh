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
echo "  OMP threads  : $OMP_NUM_THREADS"
echo "  Started      : $(date)"
echo "============================================================"
echo ""

# ------------------------------------------------------------------------------
# Optional clean
# ------------------------------------------------------------------------------
if [[ "$CLEAN" -eq 1 ]]; then
  echo "--- make clean ---"
  make -f "$MAKEFILE" -C "$APOST3D_PATH" clean
  echo ""
fi

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
