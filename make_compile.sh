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
# Checks, cheapest first: every prerequisite is checked before libxc, whose
# build takes minutes, is started.
# ------------------------------------------------------------------------------
MAKEFILE="$APOST3D_PATH/Makefile"

if [[ ! -f "$MAKEFILE" ]]; then
  echo "ERROR: Makefile not found at $APOST3D_PATH"
  echo "       Make sure APOST3D_PATH is set correctly."
  exit 1
fi

if ! command -v make &>/dev/null; then
  echo "ERROR: make not found on PATH."
  echo "       Install with:"
  echo "         macOS:  xcode-select --install"
  echo "         Ubuntu: sudo apt install make"
  echo "         Fedora: sudo dnf install make"
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
# libxc: use an existing installation if it passes the test below, else build
# the pinned release (compile_libxc.sh). Candidates, in order:
#   1. LIBXC_DIR set in the environment -> must pass, or the build stops.
#   2. pkg-config (libxcf03; Homebrew's keg added to its search path).
#   3. The bundled copy built earlier by compile_libxc.sh.
#   4. Build the bundled copy now.
# The test compiles, links and runs a program that uses the libxc interface
# APOST-3D uses: it catches Fortran modules written by another gfortran
# version, a changed interface, a library that links but cannot be loaded at
# run time, and a version older than LIBXC_MIN_VERSION (newer ones are
# taken). The choice is recorded in objects/libxc.mk for the Makefile.
# ------------------------------------------------------------------------------
# Read from compile_libxc.sh (the one place they're set) rather than
# repeated here, so a version bump only ever needs editing once.
LIBXC_VERSION="$(grep -m1 '^LIBXC_VERSION=' "$APOST3D_PATH/compile_libxc.sh" | sed -E 's/^LIBXC_VERSION="([^"]+)"/\1/')"
LIBXC_MIN_VERSION="$(grep -m1 '^LIBXC_MIN_VERSION=' "$APOST3D_PATH/compile_libxc.sh" | sed -E 's/^LIBXC_MIN_VERSION="([^"]+)"/\1/')"
BUNDLED_LIBXC_DIR="$APOST3D_PATH/libxc-${LIBXC_VERSION}"

LIBXC_TEST_DIR="$(mktemp -d)"
cat > "$LIBXC_TEST_DIR/t.f90" <<'EOF'
program t
  use xc_f03_lib_m
  implicit none
  type(xc_f03_func_t) :: f
  type(xc_f03_func_info_t) :: info
  type(xc_f03_func_reference_t) :: ref
  integer :: vmaj, vmin, vmic, i
  double precision :: rho(1), sigma(1), lapl(1), tau(1), exc(1)
  character(len=256) :: text
  rho = 0.1d0; sigma = 0.01d0; lapl = 0.0d0; tau = 0.05d0
  call xc_f03_version(vmaj, vmin, vmic)
  i = XC_FAMILY_LDA + XC_FAMILY_GGA + XC_FAMILY_MGGA + XC_FAMILY_HYB_MGGA &
    + XC_EXCHANGE + XC_CORRELATION + XC_KINETIC + XC_POLARIZED
! B3LYP: a global hybrid GGA with 20% exact exchange
  call xc_f03_func_init(f, 402, XC_UNPOLARIZED)
  info = xc_f03_func_get_info(f)
  if (xc_f03_func_info_get_family(info) /= XC_FAMILY_HYB_GGA .and. &
      xc_f03_func_info_get_family(info) /= XC_FAMILY_GGA) stop 2
  if (xc_f03_func_info_get_kind(info) /= XC_EXCHANGE_CORRELATION) stop 3
  if (iand(xc_f03_func_info_get_flags(info), XC_FLAGS_HYB_CAM + XC_FLAGS_HYB_CAMY &
      + XC_FLAGS_HYB_LC + XC_FLAGS_HYB_LCY + XC_FLAGS_VV10) /= 0) stop 4
  if (abs(xc_f03_hyb_exx_coef(f) - 0.2d0) > 1.0d-12) stop 5
  text = xc_f03_func_info_get_name(info)
  i = 0
  ref = xc_f03_func_info_get_references(info, i)
  text = xc_f03_func_reference_get_ref(ref)
  call xc_f03_gga_exc(f, int(1,8), rho, sigma, exc)
  call xc_f03_func_end(f)
! Slater exchange (LDA) and TPSS exchange (meta-GGA): the other call shapes
  call xc_f03_func_init(f, 1, XC_UNPOLARIZED)
  call xc_f03_lda_exc(f, int(1,8), rho, exc)
  call xc_f03_func_end(f)
  call xc_f03_func_init(f, 202, XC_UNPOLARIZED)
  call xc_f03_mgga_exc(f, int(1,8), rho, sigma, lapl, tau, exc)
  call xc_f03_func_end(f)
  print '(i0,".",i0,".",i0)', vmaj, vmin, vmic
end program t
EOF

# libxc_try INC LIB: runs the test against one libxc. Sets LIBXC_FOUND_VERSION,
# or LIBXC_WHY (the reason) and returns 1.
libxc_try() {
  local inc="$1" lib="$2" log="$LIBXC_TEST_DIR/log"
  LIBXC_FOUND_VERSION=""
  LIBXC_WHY=""
  rm -f "$LIBXC_TEST_DIR/t.o" "$LIBXC_TEST_DIR/t"
  # shellcheck disable=SC2086  # inc/lib hold several space-separated flags
  if ! gfortran $inc -J"$LIBXC_TEST_DIR" -c "$LIBXC_TEST_DIR/t.f90" \
      -o "$LIBXC_TEST_DIR/t.o" >"$log" 2>&1; then
    if grep -q "Cannot open module file" "$log"; then
      LIBXC_WHY="no Fortran interface (xc_f03_lib_m.mod not found)"
    elif grep -q "different version of GNU Fortran" "$log"; then
      LIBXC_WHY="its Fortran modules were written by another gfortran version"
    else
      LIBXC_WHY="its Fortran interface differs from the one APOST-3D uses ($( (grep -m1 -A1 'Error' "$log" || true) | tail -1 | sed 's/^ *//'))"
    fi
    return 1
  fi
  # shellcheck disable=SC2086
  if ! gfortran "$LIBXC_TEST_DIR/t.o" $lib -o "$LIBXC_TEST_DIR/t" >"$log" 2>&1; then
    LIBXC_WHY="a test program does not link against it"
    return 1
  fi
  if ! LIBXC_FOUND_VERSION="$("$LIBXC_TEST_DIR/t" 2>"$log")"; then
    LIBXC_WHY="a test program linked against it does not run (a shared library not found at run time? check LD_LIBRARY_PATH)"
    return 1
  fi
  if ! printf '%s\n%s\n' "$LIBXC_MIN_VERSION" "$LIBXC_FOUND_VERSION" | sort -C -V; then
    LIBXC_WHY="version $LIBXC_FOUND_VERSION is older than the oldest supported, $LIBXC_MIN_VERSION"
    return 1
  fi
  return 0
}

LIBXC_SOURCE=""
if [[ -n "${LIBXC_DIR:-}" ]]; then
  LIBXC_INC="-I$LIBXC_DIR/include"
  LIBXC_LIB="-L$LIBXC_DIR/lib -lxcf03 -lxc -lm"
  if ! libxc_try "$LIBXC_INC" "$LIBXC_LIB"; then
    echo "ERROR: the libxc in LIBXC_DIR ($LIBXC_DIR) cannot be used:"
    echo "       $LIBXC_WHY."
    echo "       Unset LIBXC_DIR to let this script find or build one."
    rm -rf "$LIBXC_TEST_DIR"
    exit 1
  fi
  LIBXC_SOURCE="LIBXC_DIR ($LIBXC_DIR)"
else
  BREW_LIBXC_PREFIX="$(brew --prefix libxc 2>/dev/null || true)"
  if [[ -n "$BREW_LIBXC_PREFIX" ]]; then
    export PKG_CONFIG_PATH="$BREW_LIBXC_PREFIX/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  fi
  if command -v pkg-config &>/dev/null && pkg-config --exists libxcf03 2>/dev/null; then
    LIBXC_INC="$(pkg-config --cflags libxcf03)"
    LIBXC_LIB="$(pkg-config --libs --static libxcf03) -lm"
    if libxc_try "$LIBXC_INC" "$LIBXC_LIB"; then
      LIBXC_SOURCE="pkg-config"
    else
      echo "  libxc: the one found by pkg-config ($(pkg-config --modversion libxcf03)) cannot be"
      echo "         used: $LIBXC_WHY."
      echo "         Using the bundled copy instead."
    fi
  fi
  if [[ -z "$LIBXC_SOURCE" ]]; then
    LIBXC_INC="-I$BUNDLED_LIBXC_DIR/include"
    LIBXC_LIB="-L$BUNDLED_LIBXC_DIR/lib -lxcf03 -lxc -lm"
    if [[ -f "$BUNDLED_LIBXC_DIR/lib/libxcf03.a" ]]; then
      if libxc_try "$LIBXC_INC" "$LIBXC_LIB"; then
        LIBXC_SOURCE="bundled copy ($BUNDLED_LIBXC_DIR)"
      else
        echo "  libxc: the bundled copy cannot be used ($LIBXC_WHY) -- rebuilding it."
      fi
    else
      echo "  libxc: not found -- fetching and building the pinned ${LIBXC_VERSION} release"
    fi
    if [[ -z "$LIBXC_SOURCE" ]]; then
      echo ""
      NTHREADS="$NTHREADS" bash "$APOST3D_PATH/compile_libxc.sh"
      echo ""
      if ! libxc_try "$LIBXC_INC" "$LIBXC_LIB"; then
        echo "ERROR: the libxc just built cannot be used: $LIBXC_WHY."
        rm -rf "$LIBXC_TEST_DIR"
        exit 1
      fi
      LIBXC_SOURCE="bundled copy, just built ($BUNDLED_LIBXC_DIR)"
    fi
  fi
fi
rm -rf "$LIBXC_TEST_DIR"
echo "  libxc $LIBXC_FOUND_VERSION found via: $LIBXC_SOURCE"
if [[ "$LIBXC_FOUND_VERSION" != "$LIBXC_VERSION" ]]; then
  echo "  note: the test references were made with libxc $LIBXC_VERSION. With"
  echo "        $LIBXC_FOUND_VERSION, 'make test' (numbers) should pass; 'make test-strict'"
  echo "        can differ in libxc's own text, e.g. the functional citations."
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

# Record the libxc checked above for the Makefile, so that a later plain
# `make` (e.g. `make test`) uses the same one instead of detecting it again.
cat > "$APOST3D_PATH/objects/libxc.mk" <<EOF
# Written by make_compile.sh: libxc $LIBXC_FOUND_VERSION, $LIBXC_SOURCE
LIBXC_INC := $LIBXC_INC
LIBXC_LIB := $LIBXC_LIB
EOF

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
