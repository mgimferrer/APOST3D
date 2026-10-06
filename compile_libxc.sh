#!/usr/bin/env bash
# ==============================================================================
# compile_libxc.sh — fetch, verify, and build the pinned libxc release
#
# libxc ships source-only, no prebuilt binaries or published checksums.
# Fetches the pinned release tag archive from upstream, verifies it
# against the recorded SHA256, and builds+installs it via CMake into
# libxc-<version>/ under APOST3D_PATH.
#
# CMake over Autotools: no bootstrap step needed (the tag archive ships
# no pre-generated `configure`), fewer prerequisites, and Autotools is
# deprecated upstream as of libxc 8.0.0.
#
# Normally invoked automatically by make_compile.sh, which probes for an
# already-usable libxc first and only falls back to this script.
#
# Usage:
#   export APOST3D_PATH=/path/to/APOST3D
#   bash compile_libxc.sh
#
# No internet egress (air-gapped HPC): pre-place the tarball at
# $APOST3D_PATH/libxc-<version>.tar.gz (matching LIBXC_SHA256 below) and
# it will be used as-is instead of downloaded.
# ==============================================================================

set -euo pipefail

# Pinned version -- don't auto-track "latest". Update LIBXC_SHA256 together
# with LIBXC_VERSION when the pin moves.
LIBXC_VERSION="7.1.2"
LIBXC_URL="https://gitlab.com/libxc/libxc/-/archive/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
LIBXC_SHA256="c517ce61820ea8114664a4280b6a6bc74a4f22f1fd1ea4ddecd6df0caeeae4f4"
LIBXC_CMAKE_MIN="3.21"
# Oldest libxc an existing installation may have to be used instead of
# building the pinned one (make_compile.sh checks it); newer ones are taken.
# shellcheck disable=SC2034  # read by make_compile.sh
LIBXC_MIN_VERSION="7.0.0"

# APOST3D_PATH defaults to this script's directory, like make_compile.sh
APOST3D_PATH="${APOST3D_PATH:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

LIBXCDIR="${APOST3D_PATH}/libxc-${LIBXC_VERSION}"
TARBALL="${APOST3D_PATH}/libxc-${LIBXC_VERSION}.tar.gz"

# Compilers: the plain gfortran on PATH, the one the Makefile uses (libxc's
# Fortran modules must come from the same gfortran as APOST-3D). The C
# compiler is the gcc of the same version when it exists under a versioned
# name (Homebrew: gcc-16; macOS's own gcc is AppleClang), else plain gcc.
if ! command -v gfortran &>/dev/null; then
  echo "ERROR: gfortran not found. Install with:"
  echo "  macOS:  brew install gcc"
  echo "  Ubuntu: sudo apt install gfortran"
  exit 1
fi
FC_CMD=gfortran
FC_MAJOR="$(gfortran -dumpversion | cut -d. -f1)"
if command -v "gcc-${FC_MAJOR}" &>/dev/null; then
  CC_CMD="gcc-${FC_MAJOR}"
elif command -v gcc &>/dev/null; then
  CC_CMD=gcc
else
  echo "ERROR: gcc not found. Install with:"
  echo "  macOS:  brew install gcc"
  echo "  Ubuntu: sudo apt install gcc"
  exit 1
fi

echo "Using C compiler  : $CC_CMD  ($(${CC_CMD} --version | head -1))"
echo "Using FC compiler : $FC_CMD  ($(${FC_CMD} --version | head -1))"

# CMake >= 3.21 required (libxc's own CMakeLists.txt). Checked explicitly
# so a version mismatch fails clearly instead of mid-configure.
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

# The other tools, also before anything is downloaded.
if ! command -v tar &>/dev/null; then
  echo "ERROR: tar not found."
  exit 1
fi
if command -v sha256sum &>/dev/null; then
  SHA256_CMD="sha256sum"
elif command -v shasum &>/dev/null; then
  SHA256_CMD="shasum -a 256"
else
  echo "ERROR: neither sha256sum nor shasum found -- cannot verify the download."
  exit 1
fi
if [[ ! -f "$TARBALL" ]] && ! command -v curl &>/dev/null; then
  echo "ERROR: curl not found, needed to download libxc ${LIBXC_VERSION}."
  echo "       Install it, or download the tarball elsewhere and place it at:"
  echo "         $TARBALL"
  echo "       (from $LIBXC_URL)"
  exit 1
fi
echo ""

# Fetch (or reuse a pre-placed tarball), then verify checksum.
cd "$APOST3D_PATH"

if [[ -f "$TARBALL" ]]; then
  echo "Found existing $TARBALL — reusing it (skipping download)."
else
  DOWNLOADED=1
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
# shellcheck disable=SC2086  # SHA256_CMD may hold a command and its option
ACTUAL_SHA256=$($SHA256_CMD "$TARBALL" | awk '{print $1}')

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

# Always extract into a clean tree.
echo "Extracting libxc-${LIBXC_VERSION}.tar.gz ..."
rm -rf "$LIBXCDIR"
tar -xzf "$TARBALL"
echo ""

# Configure, build, install via CMake into LIBXCDIR itself. Static-only
# (apost3d links libxc in directly). ENABLE_FORTRAN defaults OFF upstream.
# CMAKE_INSTALL_LIBDIR is pinned to "lib" -- CMake's GNUInstallDirs module
# defaults 64-bit RHEL/CentOS/Rocky-family systems to "lib64" instead, and
# every consumer here (Makefile, make_compile.sh) hardcodes .../lib/.
cd "$LIBXCDIR"
echo "Running cmake configure ..."
cmake -B build \
  -DCMAKE_INSTALL_PREFIX="$LIBXCDIR" \
  -DCMAKE_INSTALL_LIBDIR=lib \
  -DCMAKE_C_COMPILER="$CC_CMD" \
  -DCMAKE_Fortran_COMPILER="$FC_CMD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_FORTRAN=ON \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF
echo ""

# NTHREADS (as for make_compile.sh), else 8 or the CPU count if smaller
if [[ -n "${NTHREADS:-}" ]]; then
  NCPU="$NTHREADS"
else
  NCPU=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
  NCPU=$(( NCPU < 8 ? NCPU : 8 ))
fi
echo "Building with ${NCPU} parallel jobs ..."
cmake --build build -j"$NCPU"
echo ""

echo "Installing into $LIBXCDIR ..."
cmake --install build
echo ""

if [[ ! -f "$LIBXCDIR/lib/libxcf03.a" ]]; then
  echo "ERROR: install finished but $LIBXCDIR/lib/libxcf03.a is missing."
  exit 1
fi

# A tarball this script downloaded is not needed any more; one placed by
# hand (machine without internet access) is kept for a later rebuild.
if [[ "${DOWNLOADED:-0}" -eq 1 ]]; then
  rm -f "$TARBALL"
fi

echo "============================================================"
echo "  libxc-${LIBXC_VERSION} built successfully."
echo "  Headers/modules : $LIBXCDIR/include/"
echo "  Libraries       : $LIBXCDIR/lib/libxc.a, libxcf03.a"
echo "============================================================"
echo ""
echo "Next step:"
echo "  bash \$APOST3D_PATH/make_compile.sh"
