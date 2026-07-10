#!/usr/bin/env bash
# ==============================================================================
# compile_libxc_gfortran.sh — build libxc-4.2.3 with GCC/gfortran
#
# Replaces compile_libxc.sh (which requires Intel icx/ifort).
# Works on macOS (Homebrew gcc) and Linux (system gcc/gfortran).
#
# Usage:
#   export APOST3D_PATH=/path/to/APOST3D
#   bash compile_libxc_gfortran.sh
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Require APOST3D_PATH
# ------------------------------------------------------------------------------
if [[ -z "${APOST3D_PATH:-}" ]]; then
  echo "ERROR: APOST3D_PATH is not set."
  echo "       Run:  export APOST3D_PATH=/path/to/APOST3D"
  exit 1
fi

LIBXCDIR="${APOST3D_PATH}/libxc-4.2.3"

# ------------------------------------------------------------------------------
# Auto-detect gcc / gfortran
# Homebrew on macOS installs versioned binaries (gcc-14, gfortran-14, …).
# Prefer the versioned name so we avoid picking up Apple's cc wrapper.
# ------------------------------------------------------------------------------
find_compiler() {
  local name="$1"
  # Try versioned names 15 down to 10
  for v in 15 14 13 12 11 10; do
    if command -v "${name}-${v}" &>/dev/null; then
      echo "${name}-${v}"
      return
    fi
  done
  # Fall back to unversioned
  if command -v "${name}" &>/dev/null; then
    echo "${name}"
    return
  fi
  echo ""
}

CC_CMD=$(find_compiler gcc)
FC_CMD=$(find_compiler gfortran)

if [[ -z "$CC_CMD" ]]; then
  echo "ERROR: gcc not found. Install with:"
  echo "  macOS:  brew install gcc"
  echo "  Ubuntu: sudo apt install gcc"
  exit 1
fi
if [[ -z "$FC_CMD" ]]; then
  echo "ERROR: gfortran not found. Install with:"
  echo "  macOS:  brew install gcc"
  echo "  Ubuntu: sudo apt install gfortran"
  exit 1
fi

echo "Using C compiler  : $CC_CMD  ($(${CC_CMD} --version | head -1))"
echo "Using FC compiler : $FC_CMD  ($(${FC_CMD} --version | head -1))"
echo ""

# ------------------------------------------------------------------------------
# Extract source (always start from a clean tree)
# ------------------------------------------------------------------------------
cd "$APOST3D_PATH"

if [[ ! -f libxc-4.2.3.tar.gz ]]; then
  echo "ERROR: libxc-4.2.3.tar.gz not found in $APOST3D_PATH"
  exit 1
fi

echo "Extracting libxc-4.2.3.tar.gz ..."
rm -rf libxc-4.2.3
tar -xzf libxc-4.2.3.tar.gz
echo ""

# ------------------------------------------------------------------------------
# Configure
# ------------------------------------------------------------------------------
cd "$LIBXCDIR"

echo "Running configure ..."
CC="$CC_CMD" \
FC="$FC_CMD" \
CFLAGS="-O3" \
FCFLAGS="-O3 -ffixed-line-length-132" \
  ./configure --prefix="$LIBXCDIR" --enable-shared=no
echo ""

# ------------------------------------------------------------------------------
# Build and install
# ------------------------------------------------------------------------------
NCPU=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4)
echo "Building with ${NCPU} parallel jobs ..."
make -j"$NCPU"
echo ""

echo "Installing into $LIBXCDIR ..."
make install
echo ""

# ------------------------------------------------------------------------------
# Copy F90 interfaces needed by Makefile_gfortran
# ------------------------------------------------------------------------------
cp src/libxc_funcs.f90 "$LIBXCDIR/"
cp src/libxc.f90       "$LIBXCDIR/"

echo "============================================================"
echo "  libxc-4.2.3 built successfully."
echo "  Headers : $LIBXCDIR/include/"
echo "  Library : $LIBXCDIR/lib/libxc.a"
echo "  F90 src : $LIBXCDIR/libxc_funcs.f90"
echo "            $LIBXCDIR/libxc.f90"
echo "============================================================"
echo ""
echo "Next step:"
echo "  mkdir -p \$APOST3D_PATH/objects"
echo "  make -f \$APOST3D_PATH/Makefile_gfortran -C \$APOST3D_PATH all"
