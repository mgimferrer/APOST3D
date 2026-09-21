#!/usr/bin/env bash
# ==============================================================================
# compile_libxc.sh — fetch, verify, and build libxc 7.1.2 with GCC/gfortran
#
# libxc is distributed as source only (no prebuilt binaries, no published
# checksums beyond what we record ourselves here) — see libxc.gitlab.io/download.
# This script fetches the exact pinned release tag archive from upstream,
# verifies it against a recorded SHA256, and builds+installs it via CMake
# into libxc-<version>/ under APOST3D_PATH — the same self-contained-prefix
# convention the old libxc-4.2.3 setup used.
#
# CMake, not Autotools: verified 2026-09-21 that CMake (>= 3.21, required by
# libxc's own CMakeLists.txt) builds libxc 7.1.2's Fortran interface
# correctly — full parity with an Autotools build, confirmed by running
# APOST3D's entire test suite against both. CMake also needs no bootstrap
# step (Autotools' `autoreconf -fi`, and the autoconf/automake/libtool
# prerequisites that implies, since the release tag archive ships no
# pre-generated `configure`) and correctly detects the platform's own
# archiver on its own (no GNU-vs-BSD ar/ranlib workaround needed, unlike
# the old libxc-4.2.3 Autotools setup). Autotools is upstream's own
# recommended default for older libxc releases, and is explicitly flagged
# for removal in libxc 8.0.0 — CMake is also just the forward-compatible
# choice for whenever the pinned version is next bumped.
#
# Usually you don't need to run this directly: make_compile.sh probes for
# an already-usable libxc (an explicit LIBXC_DIR, a system/conda/Homebrew
# install registered with pkg-config) first, and only falls back to this
# script when nothing suitable is found.
#
# Usage:
#   export APOST3D_PATH=/path/to/APOST3D
#   bash compile_libxc.sh
#
# Air-gapped / no internet egress: pre-download the exact tarball below
# (matching the recorded SHA256) and place it at
# $APOST3D_PATH/libxc-7.1.2.tar.gz before running this script — it will
# be used as-is instead of fetching.
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Pinned version. Bumping this is a deliberate, reviewed action (different
# libxc releases can carry different numerics for edge cases) — don't
# auto-track "latest". Update LIBXC_SHA256 together with LIBXC_VERSION if
# the pin ever moves.
# ------------------------------------------------------------------------------
LIBXC_VERSION="7.1.2"
LIBXC_URL="https://gitlab.com/libxc/libxc/-/archive/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
LIBXC_SHA256="c517ce61820ea8114664a4280b6a6bc74a4f22f1fd1ea4ddecd6df0caeeae4f4"
LIBXC_CMAKE_MIN="3.21"

# ------------------------------------------------------------------------------
# Require APOST3D_PATH
# ------------------------------------------------------------------------------
if [[ -z "${APOST3D_PATH:-}" ]]; then
  echo "ERROR: APOST3D_PATH is not set."
  echo "       Run:  export APOST3D_PATH=/path/to/APOST3D"
  exit 1
fi

LIBXCDIR="${APOST3D_PATH}/libxc-${LIBXC_VERSION}"
TARBALL="${APOST3D_PATH}/libxc-${LIBXC_VERSION}.tar.gz"

# ------------------------------------------------------------------------------
# Auto-detect gcc / gfortran
# Homebrew on macOS installs versioned binaries (gcc-14, gfortran-14, …),
# and macOS's own /usr/bin/gcc is actually AppleClang in disguise — prefer
# the versioned name so we don't silently hand CMake the wrong C compiler
# for a Fortran-interoperable build.
# ------------------------------------------------------------------------------
find_compiler() {
  local name="$1"
  for v in 16 15 14 13 12 11 10; do
    if command -v "${name}-${v}" &>/dev/null; then
      echo "${name}-${v}"
      return
    fi
  done
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
# CMake prerequisite (>= 3.21, per libxc's own CMakeLists.txt).
# Not needed by anything else in this codebase, so check explicitly with a
# real version comparison rather than let a cryptic CMake error surface
# mid-configure. If the system cmake is too old (common on conservative
# HPC distros), `pip install --user cmake` or `pipx install cmake` gets a
# modern prebuilt one with no compilation and no root needed.
# ------------------------------------------------------------------------------
if ! command -v cmake &>/dev/null; then
  echo "ERROR: cmake not found (>= ${LIBXC_CMAKE_MIN} required)."
  echo "       Install with:"
  echo "         macOS:  brew install cmake"
  echo "         Ubuntu: sudo apt install cmake"
  echo "         Fedora: sudo dnf install cmake"
  echo "       No root / too old on your system? pip install --user cmake"
  echo "       (ships a modern prebuilt binary, nothing to compile)."
  exit 1
fi

CMAKE_VERSION="$(cmake --version | head -1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+')"
if ! printf '%s\n%s\n' "$LIBXC_CMAKE_MIN" "$CMAKE_VERSION" | sort -C -V; then
  echo "ERROR: cmake ${CMAKE_VERSION} found, but libxc ${LIBXC_VERSION} needs >= ${LIBXC_CMAKE_MIN}."
  echo "       No root / too old on your system? pip install --user cmake"
  echo "       (ships a modern prebuilt binary, nothing to compile)."
  exit 1
fi
echo "Using cmake       : $(command -v cmake)  (${CMAKE_VERSION})"
echo ""

# ------------------------------------------------------------------------------
# Fetch (or reuse a pre-placed tarball) and verify checksum.
# ------------------------------------------------------------------------------
cd "$APOST3D_PATH"

if [[ -f "$TARBALL" ]]; then
  echo "Found existing $TARBALL — reusing it (skipping download)."
else
  echo "Downloading libxc ${LIBXC_VERSION} from upstream..."
  echo "  $LIBXC_URL"
  if ! curl -fL --retry 3 -o "$TARBALL" "$LIBXC_URL"; then
    echo ""
    echo "ERROR: download failed. If this machine has no internet egress"
    echo "       (e.g. an air-gapped HPC node), download the tarball"
    echo "       elsewhere and place it at:"
    echo "         $TARBALL"
    echo "       then re-run this script."
    rm -f "$TARBALL"
    exit 1
  fi
fi

echo "Verifying SHA256 checksum..."
ACTUAL_SHA256=""
if command -v sha256sum &>/dev/null; then
  ACTUAL_SHA256=$(sha256sum "$TARBALL" | awk '{print $1}')
elif command -v shasum &>/dev/null; then
  ACTUAL_SHA256=$(shasum -a 256 "$TARBALL" | awk '{print $1}')
else
  echo "ERROR: neither sha256sum nor shasum found — cannot verify checksum."
  exit 1
fi

if [[ "$ACTUAL_SHA256" != "$LIBXC_SHA256" ]]; then
  echo "ERROR: checksum mismatch for $TARBALL"
  echo "  expected: $LIBXC_SHA256"
  echo "  actual:   $ACTUAL_SHA256"
  echo "The downloaded/pre-placed file does not match the pinned libxc"
  echo "${LIBXC_VERSION} release — refusing to build from it. Delete it and"
  echo "re-run this script to fetch a fresh copy."
  exit 1
fi
echo "  OK: $ACTUAL_SHA256"
echo ""

# ------------------------------------------------------------------------------
# Extract (always start from a clean tree)
# ------------------------------------------------------------------------------
echo "Extracting libxc-${LIBXC_VERSION}.tar.gz ..."
rm -rf "$LIBXCDIR"
tar -xzf "$TARBALL"
echo ""

# ------------------------------------------------------------------------------
# Configure, build, install via CMake — installs into itself (LIBXCDIR) as
# prefix, same self-contained convention as the old libxc-4.2.3 setup.
# Static-only: apost3d links libxc in directly, no need to ship/rpath a
# shared lib. ENABLE_FORTRAN is OFF by default upstream, must opt in.
# ------------------------------------------------------------------------------
cd "$LIBXCDIR"
echo "Running cmake configure ..."
cmake -B build \
  -DCMAKE_INSTALL_PREFIX="$LIBXCDIR" \
  -DCMAKE_C_COMPILER="$CC_CMD" \
  -DCMAKE_Fortran_COMPILER="$FC_CMD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_FORTRAN=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF
echo ""

NCPU=$(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4)
echo "Building with ${NCPU} parallel jobs ..."
cmake --build build -j"$NCPU"
echo ""

echo "Installing into $LIBXCDIR ..."
cmake --install build
echo ""

echo "============================================================"
echo "  libxc-${LIBXC_VERSION} built successfully."
echo "  Headers/modules : $LIBXCDIR/include/"
echo "  Libraries       : $LIBXCDIR/lib/libxc.a, libxcf03.a"
echo "============================================================"
echo ""
echo "Next step:"
echo "  bash \$APOST3D_PATH/make_compile.sh"
