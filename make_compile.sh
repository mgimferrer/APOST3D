#!/usr/bin/env bash
# ==============================================================================
# make_compile.sh — compile APOST-3D with GCC/gfortran
#
# Builds apost3d and the utilities (utils/), fetching and
# building libxc first if no usable copy is found.
#
# Usage:
#   bash make_compile.sh [clean] [NTHREADS=N] [ARCH=cpu] [help]
#   (APOST3D_PATH defaults to this script's directory)
#
# Arguments (same grammar as `make`: bare words are actions, KEY=value sets
# a parameter — and NTHREADS is spelled identically to `make test NTHREADS=N`):
#   clean         Run 'make clean' before building
#   NTHREADS=<N>  Parallel compile jobs, and the OMP_NUM_THREADS suggested
#                 for runs and tests (default: 8, or fewer if the machine
#                 has fewer CPUs; use a small value on a shared login node)
#   ARCH=<cpu>    Target CPU. Default: generic code that runs on any CPU of
#                 this architecture (the safe choice for clusters with nodes
#                 of different ages). ARCH=native: this machine's CPU only,
#                 may be faster but can stop with "Illegal instruction"
#                 elsewhere. Any gcc -march value works (e.g. x86-64-v3).
#   help          Show this message (also: --help, -h)
#
# Examples:
#   bash make_compile.sh
#   bash make_compile.sh NTHREADS=4
#   bash make_compile.sh clean ARCH=native
#   bash make_compile.sh help
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Defaults
# ------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APOST3D_PATH="${APOST3D_PATH:-$SCRIPT_DIR}"
CLEAN=0
NTHREADS=""
ARCH=""

# ------------------------------------------------------------------------------
# Argument parsing — bare words are actions, KEY=value sets a parameter,
# mirroring `make <target> VAR=value`.
# ------------------------------------------------------------------------------
for arg in "$@"; do
  case "$arg" in
    help|--help|-h)
      sed -n '2,29p' "$0" | sed 's/^# \{0,2\}//'
      exit 0
      ;;
    clean)
      CLEAN=1
      ;;
    NTHREADS=*)
      NTHREADS="${arg#NTHREADS=}"
      ;;
    ARCH=*)
      ARCH="${arg#ARCH=}"
      ;;
    *)
      echo "Unknown argument: $arg"
      echo "Run 'bash make_compile.sh help' for usage."
      exit 1
      ;;
  esac
done

# ------------------------------------------------------------------------------
# Default threads: 8, or the logical CPU count if smaller
# ------------------------------------------------------------------------------
if [[ -z "$NTHREADS" ]]; then
  NCPU=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
  NTHREADS=$(( NCPU < 8 ? NCPU : 8 ))
fi

# ------------------------------------------------------------------------------
# Checks
# ------------------------------------------------------------------------------
MAKEFILE="$APOST3D_PATH/Makefile"

if [[ ! -f "$MAKEFILE" ]]; then
  echo "ERROR: Makefile not found at $APOST3D_PATH"
  echo "       Make sure APOST3D_PATH is set correctly."
  exit 1
fi

# Read from compile_libxc.sh (the one place it's pinned) rather than
# repeated here, so a version bump only ever needs editing once.
LIBXC_VERSION="$(grep -m1 '^LIBXC_VERSION=' "$APOST3D_PATH/compile_libxc.sh" | sed -E 's/^LIBXC_VERSION="([^"]+)"/\1/')"
BUNDLED_LIBXC_A="$APOST3D_PATH/libxc-${LIBXC_VERSION}/lib/libxcf03.a"

# ------------------------------------------------------------------------------
# libxc preflight/auto-build: probe for an already-usable libxc first, same
# layered order as the Makefile (LIBXC_DIR override, then pkg-config, then
# Homebrew's keg-only prefix), and only fetch+build our own bundled copy if
# none of those are found.
# ------------------------------------------------------------------------------
LIBXC_READY=0
if [[ -n "${LIBXC_DIR:-}" ]]; then
  echo "  libxc: using LIBXC_DIR override ($LIBXC_DIR)"
  LIBXC_READY=1
else
  BREW_LIBXC_PREFIX="$(brew --prefix libxc 2>/dev/null || true)"
  if [[ -n "$BREW_LIBXC_PREFIX" ]]; then
    export PKG_CONFIG_PATH="$BREW_LIBXC_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  fi
  if command -v pkg-config &>/dev/null && pkg-config --exists libxcf03 2>/dev/null; then
    echo "  libxc: found via pkg-config ($(pkg-config --modversion libxcf03))"
    LIBXC_READY=1
  elif [[ -f "$BUNDLED_LIBXC_A" ]]; then
    # gfortran .mod files aren't portable across compiler versions, and a
    # file-existence check alone can't catch a stale bundled build (e.g.
    # after an HPC `module load` swap) -- try actually reading the module,
    # same "verify by using it" approach as the OpenBLAS link-test below,
    # rather than deferring the failure to deep inside the real build with
    # a confusing "different version of GNU Fortran" error.
    LIBXC_MODTEST_DIR="$(mktemp -d)"
    cat > "$LIBXC_MODTEST_DIR/t.f90" <<'EOF'
program t
  use xc_f03_lib_m
end program t
EOF
    if gfortran "$LIBXC_MODTEST_DIR/t.f90" -I"$APOST3D_PATH/libxc-${LIBXC_VERSION}/include" \
        -c -o "$LIBXC_MODTEST_DIR/t.o" &>/dev/null; then
      echo "  libxc: using already-built bundled copy ($APOST3D_PATH/libxc-${LIBXC_VERSION})"
      LIBXC_READY=1
    else
      echo "  libxc: bundled copy exists but its .mod files aren't readable by"
      echo "         the current gfortran (likely built by a different version,"
      echo "         e.g. after an HPC module swap) -- rebuilding."
    fi
    rm -rf "$LIBXC_MODTEST_DIR"
  fi
fi

if [[ "$LIBXC_READY" -eq 0 ]]; then
  echo "  libxc: not found -- fetching and building the pinned ${LIBXC_VERSION} release"
  echo ""
  bash "$APOST3D_PATH/compile_libxc.sh"
  echo ""
  if [[ ! -f "$BUNDLED_LIBXC_A" ]]; then
    echo "ERROR: compile_libxc.sh ran but $BUNDLED_LIBXC_A still doesn't exist."
    exit 1
  fi
fi

# ------------------------------------------------------------------------------
# OpenBLAS preflight: diagonalize() in util.f needs LAPACK's dsyevd. Same
# layered detection as the Makefile's OPENBLAS_LIB -- each a fallback for
# the one above it:
#   1. OPENBLAS_DIR set in the environment -> always wins (nonstandard
#      install, HPC module that doesn't export the right paths, etc).
#   2. pkg-config -- the standards-based mechanism most package managers
#      (apt, dnf, conda, spack) register a .pc file for.
#   3. Homebrew's keg-only prefix on macOS, also fed into pkg-config's own
#      search path so step 2 catches it uniformly.
#   4. Bare -lopenblas, relying on the default linker search path or an
#      HPC `module load` that already exported LIBRARY_PATH/LD_LIBRARY_PATH.
# Links a real test program, not just a file-existence check.
# ------------------------------------------------------------------------------
if [[ -n "${OPENBLAS_DIR:-}" ]]; then
  OPENBLAS_LDFLAGS="-L$OPENBLAS_DIR/lib -lopenblas"
  OPENBLAS_METHOD="OPENBLAS_DIR override ($OPENBLAS_DIR)"
else
  BREW_OPENBLAS_PREFIX="$(brew --prefix openblas 2>/dev/null || true)"
  if [[ -n "$BREW_OPENBLAS_PREFIX" ]]; then
    export PKG_CONFIG_PATH="$BREW_OPENBLAS_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  fi
  if command -v pkg-config &>/dev/null && pkg-config --exists openblas 2>/dev/null; then
    OPENBLAS_LDFLAGS="$(pkg-config --libs openblas)"
    OPENBLAS_METHOD="pkg-config"
  elif [[ -n "$BREW_OPENBLAS_PREFIX" ]]; then
    OPENBLAS_LDFLAGS="-L$BREW_OPENBLAS_PREFIX/lib -lopenblas"
    OPENBLAS_METHOD="Homebrew prefix ($BREW_OPENBLAS_PREFIX)"
  else
    OPENBLAS_LDFLAGS="-lopenblas"
    OPENBLAS_METHOD="default linker search path"
  fi
fi

OPENBLAS_TEST_DIR="$(mktemp -d)"
cat > "$OPENBLAS_TEST_DIR/t.f90" <<'EOF'
program t
  external dsyevd
  print *, "ok"
end program t
EOF
# shellcheck disable=SC2086  # OPENBLAS_LDFLAGS is intentionally unquoted: it
# can hold multiple space-separated flags (e.g. "-L/path -lopenblas") that
# must split into separate arguments.
if ! gfortran "$OPENBLAS_TEST_DIR/t.f90" $OPENBLAS_LDFLAGS -o "$OPENBLAS_TEST_DIR/t" &>/dev/null; then
  echo "ERROR: could not link against OpenBLAS (needed for LAPACK's dsyevd,"
  echo "       used by diagonalize() in sources/util.f)."
  echo "       Tried via: $OPENBLAS_METHOD"
  echo "       Install it with:"
  echo "         macOS:  brew install openblas"
  echo "         Ubuntu: sudo apt install libopenblas-dev"
  echo "         Fedora: sudo dnf install openblas-devel"
  echo "       Nonstandard install location? Point at it directly:"
  echo "         export OPENBLAS_DIR=/path/to/openblas"
  rm -rf "$OPENBLAS_TEST_DIR"
  exit 1
fi
echo "  OpenBLAS found via: $OPENBLAS_METHOD"
rm -rf "$OPENBLAS_TEST_DIR"

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
echo "  Threads      : $NTHREADS (compile jobs, OMP_NUM_THREADS)"
echo "  Target CPU   : ${ARCH:-generic (any CPU of this architecture)}"
echo "  Started      : $(date)"
echo "============================================================"
echo ""

# ------------------------------------------------------------------------------
# Compiler/target drift check. gfortran .mod files aren't portable across
# compiler versions, and objects built for another ARCH would be mixed in:
# the Makefile's timestamp-based deps can't detect either, so stamp the
# compiler and ARCH of the last build and force a clean if they changed.
# ------------------------------------------------------------------------------
COMPILER_STAMP="$APOST3D_PATH/objects/.gfortran_version"
BUILD_ID="$GFORTRAN_VERSION_FULL ARCH=${ARCH:-generic}"

if [[ -f "$COMPILER_STAMP" ]]; then
  PREV_VERSION="$(cat "$COMPILER_STAMP")"
  if [[ "$PREV_VERSION" != "$BUILD_ID" ]] && [[ "$CLEAN" -ne 1 ]]; then
    echo "--- Compiler or target CPU change detected ---"
    echo "  Previous build : $PREV_VERSION"
    echo "  Current        : $BUILD_ID"
    echo "  Forcing 'make clean' so no object of the old build is reused."
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

# Record the compiler/ARCH used for this build (written after clean so a
# failed/interrupted build doesn't falsely mark the stamp as up to date).
mkdir -p "$APOST3D_PATH/objects"
echo "$BUILD_ID" > "$COMPILER_STAMP"

# ------------------------------------------------------------------------------
# Build main binary + utilities
# ------------------------------------------------------------------------------
BINARIES=(apost3d utils/get_energy utils/get_energy_g16 utils/gen_hirsh
          utils/wfn2fchk utils/eos_aom)
echo "--- Building apost3d and the utilities ($NTHREADS jobs) ---"
make -f "$MAKEFILE" -C "$APOST3D_PATH" -j"$NTHREADS" ARCH="$ARCH" all utils
echo ""

# ------------------------------------------------------------------------------
# Ad-hoc code signing (macOS/Apple Silicon: GCC-compiled binaries need a
# valid signature to map the dyld shared cache, or they fail to launch with
# a misleading "Library not loaded" error). codesign doesn't exist on
# Linux, so this block is a no-op there.
# ------------------------------------------------------------------------------
if command -v codesign &>/dev/null; then
  echo "--- Ad-hoc code signing (macOS) ---"
  for bin in "${BINARIES[@]}"; do
    if [[ -x "$APOST3D_PATH/$bin" ]]; then
      # Clear extended attributes first (e.g. a stray com.apple.quarantine
      # flag, which can prevent an ad-hoc signature from being trusted).
      xattr -c "$APOST3D_PATH/$bin" 2>/dev/null || true
      if codesign --force --sign - "$APOST3D_PATH/$bin" 2>&1; then
        echo "  signed: $bin"
      else
        echo "  ERROR: codesign failed for $bin — it will likely fail to run."
        echo "         Try manually:  codesign --force --sign - $APOST3D_PATH/$bin"
      fi
    fi
  done
  echo ""
fi

# ------------------------------------------------------------------------------
# Verify binaries exist
# ------------------------------------------------------------------------------
echo "--- Binaries produced ---"
for bin in "${BINARIES[@]}"; do
  if [[ -x "$APOST3D_PATH/$bin" ]]; then
    echo "  ✓  $APOST3D_PATH/$bin  ($(du -sh "$APOST3D_PATH/$bin" | cut -f1))"
  else
    echo "  ✗  $APOST3D_PATH/$bin  NOT FOUND — build may have failed"
  fi
done
echo ""

# ------------------------------------------------------------------------------
# Smoke test: launch each binary. A binary can exist and be executable yet
# still fail to launch (wrong arch, bad signature, missing shared lib) —
# only shows up at process start. Run without arguments and with stdin
# /dev/null, in a scratch directory: each one stops at once, asking for its
# input file, and any file a utility opens lands there, not here.
# ------------------------------------------------------------------------------
echo "--- Smoke test (launching each binary) ---"
SMOKE_FAILED=0
SMOKE_DIR="$(mktemp -d)"
for bin in "${BINARIES[@]}"; do
  bin_path="$APOST3D_PATH/$bin"
  [[ -x "$bin_path" ]] || continue
  smoke_out="$(cd "$SMOKE_DIR" && "$bin_path" < /dev/null 2>&1 || true)"
  if echo "$smoke_out" | grep -qi "dyld\|Library not loaded\|Killed\|Segmentation fault\|Trace/BPT"; then
    echo "  ✗  $bin failed to launch:"
    echo "$smoke_out" | sed 's/^/       /'
    SMOKE_FAILED=1
  elif [[ -z "$smoke_out" ]]; then
    echo "  ?  $bin produced no output (unexpected, but not a known failure signature)"
  else
    echo "  ✓  $bin launches correctly"
  fi
done
echo ""

rm -rf "$SMOKE_DIR"

if [[ "$SMOKE_FAILED" -eq 1 ]]; then
  echo "============================================================"
  echo "  WARNING: one or more binaries failed to launch"
  echo "============================================================"
  echo "  On macOS this is almost always a code-signing / Gatekeeper"
  echo "  issue, not a compilation problem. Try, then re-run this script:"
  echo "    xattr -cr $APOST3D_PATH"
  echo "    codesign --force --sign - $APOST3D_PATH/apost3d"
  echo "    (and the same for each program in $APOST3D_PATH/utils)"
  echo ""
fi

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
echo "  make -C $APOST3D_PATH test NTHREADS=$NTHREADS"
echo ""
echo "For all available make targets and flags:"
echo "  make -C $APOST3D_PATH help"
