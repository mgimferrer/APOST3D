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

# ------------------------------------------------------------------------------
# Force the platform-native ar/ranlib.
#
# libxc's own build (autotools) picks whatever `ar` it finds first on PATH.
# On macOS this can end up being a GNU `ar` (e.g. from Homebrew binutils, or
# bundled alongside a Homebrew GCC toolchain), which writes archives in GNU
# format (a member literally named "/" holding the symbol table). Apple's
# system `ld` only understands the BSD archive format (a `__.SYMDEF` member)
# and fails with "archive member '/' not a mach-o file" when linking against
# a GNU-format .a — this is a toolchain mismatch, not a code problem, and it
# only shows up at link time (compiling and archiving both "succeed").
#
# Fix: explicitly force Apple's own /usr/bin/ar + /usr/bin/ranlib on macOS
# (always present via Xcode Command Line Tools, always BSD-format, always
# compatible with Apple's ld). On Linux, GNU ar/ranlib is correct and
# expected, so just use whatever's on PATH there.
# ------------------------------------------------------------------------------
if [[ "$(uname -s)" == "Darwin" ]]; then
  AR_CMD=/usr/bin/ar
  RANLIB_CMD=/usr/bin/ranlib
  if [[ ! -x "$AR_CMD" ]]; then
    echo "ERROR: $AR_CMD not found — install Xcode Command Line Tools:"
    echo "  xcode-select --install"
    exit 1
  fi
else
  AR_CMD=$(command -v ar || echo ar)
  RANLIB_CMD=$(command -v ranlib || echo ranlib)
fi
echo "Using AR          : $AR_CMD"
echo "Using RANLIB      : $RANLIB_CMD"
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
AR="$AR_CMD" \
RANLIB="$RANLIB_CMD" \
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
# Safety net: regenerate the archive symbol table with the platform-native
# ranlib, regardless of what AR/RANLIB the internal build actually used for
# each .a (belt-and-braces against the GNU-vs-BSD archive mismatch above).
# ------------------------------------------------------------------------------
echo "Re-indexing archives with $RANLIB_CMD ..."
for a in "$LIBXCDIR"/lib/*.a; do
  [[ -f "$a" ]] && "$RANLIB_CMD" "$a"
done
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
