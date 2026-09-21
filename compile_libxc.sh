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

if [[ -z "${APOST3D_PATH:-}" ]]; then
  echo "ERROR: APOST3D_PATH is not set."
  echo "       Run:  export APOST3D_PATH=/path/to/APOST3D"
  exit 1
fi

LIBXCDIR="${APOST3D_PATH}/libxc-${LIBXC_VERSION}"
TARBALL="${APOST3D_PATH}/libxc-${LIBXC_VERSION}.tar.gz"

# macOS's /usr/bin/gcc is AppleClang in disguise -- prefer Homebrew's
# versioned binary (gcc-14 etc.) to avoid handing CMake the wrong compiler.
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
echo ""

# Fetch (or reuse a pre-placed tarball), then verify checksum.
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

# Always extract into a clean tree.
echo "Extracting libxc-${LIBXC_VERSION}.tar.gz ..."
rm -rf "$LIBXCDIR"
tar -xzf "$TARBALL"
echo ""

# Configure, build, install via CMake into LIBXCDIR itself. Static-only
# (apost3d links libxc in directly). ENABLE_FORTRAN defaults OFF upstream.
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
