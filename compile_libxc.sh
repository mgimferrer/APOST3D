#!/usr/bin/env bash
# ==============================================================================
# compile_libxc.sh — fetch, verify, and build libxc 7.1.2 with GCC/gfortran
#
# libxc is distributed as source only (no prebuilt binaries, no published
# checksums beyond what we record ourselves here) — see libxc.gitlab.io/download.
# This script fetches the exact pinned release tag archive from upstream,
# verifies it against a recorded SHA256, bootstraps its autotools build
# (the tag archive ships configure.ac/Makefile.am but no pre-generated
# `configure`), and builds+installs it into libxc-<version>/ under
# APOST3D_PATH — the same self-contained-prefix convention the old
# libxc-4.2.3 setup used.
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
# Homebrew on macOS installs versioned binaries (gcc-14, gfortran-14, …).
# Prefer the versioned name so we avoid picking up Apple's cc wrapper.
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
# Autotools prerequisites.
#
# The GitLab tag archive ships configure.ac/Makefile.am but no pre-generated
# `configure` (unlike a curated release dist tarball) — building it needs
# `autoreconf -fi` first, which needs autoconf/automake/libtool installed.
# Not needed by anything else in this codebase, so check explicitly rather
# than let a cryptic "autoreconf: command not found" surface mid-script.
#
# macOS-specific gotcha: /usr/bin/libtool is Apple's own static-library
# archiver, unrelated to GNU libtool that autoreconf/LT_INIT actually
# needs. Homebrew installs the real one under the keg-only names
# glibtool/glibtoolize precisely to avoid clashing with Apple's — same
# class of name collision as the AR/RANLIB workaround below, just for a
# different tool. Prepending its gnubin dir to PATH makes plain
# `libtoolize` resolve to the GNU one for the rest of this script.
# ------------------------------------------------------------------------------
if [[ "$(uname -s)" == "Darwin" ]]; then
  BREW_LIBTOOL_PREFIX="$(brew --prefix libtool 2>/dev/null || true)"
  if [[ -n "$BREW_LIBTOOL_PREFIX" && -d "$BREW_LIBTOOL_PREFIX/libexec/gnubin" ]]; then
    export PATH="$BREW_LIBTOOL_PREFIX/libexec/gnubin:$PATH"
  fi
fi

MISSING_TOOLS=()
for t in autoconf automake libtoolize; do
  command -v "$t" &>/dev/null || MISSING_TOOLS+=("$t")
done
if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
  echo "ERROR: missing autotools prerequisite(s): ${MISSING_TOOLS[*]}"
  echo "       Needed to bootstrap libxc's build (autoreconf -fi)."
  echo "       Install with:"
  echo "         macOS:  brew install autoconf automake libtool"
  echo "         Ubuntu: sudo apt install autoconf automake libtool"
  echo "         Fedora: sudo dnf install autoconf automake libtool"
  exit 1
fi
echo "Using autoconf    : $(command -v autoconf)  ($(autoconf --version | head -1))"
echo "Using automake     : $(command -v automake)  ($(automake --version | head -1))"
echo "Using libtoolize  : $(command -v libtoolize)"

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
# Bootstrap the autotools build (no pre-generated `configure` in the tag
# archive — see the header comment above).
# ------------------------------------------------------------------------------
cd "$LIBXCDIR"
echo "Running autoreconf -fi ..."
autoreconf -fi
echo ""

# ------------------------------------------------------------------------------
# Configure — installs into itself (LIBXCDIR) as prefix, same
# self-contained convention as the old libxc-4.2.3 setup. Static-only:
# apost3d links libxc in directly, no need to ship/rpath a shared lib.
# ------------------------------------------------------------------------------
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

echo "============================================================"
echo "  libxc-${LIBXC_VERSION} built successfully."
echo "  Headers/modules : $LIBXCDIR/include/"
echo "  Libraries       : $LIBXCDIR/lib/libxc.a, libxcf03.a"
echo "============================================================"
echo ""
echo "Next step:"
echo "  bash \$APOST3D_PATH/make_compile.sh"
