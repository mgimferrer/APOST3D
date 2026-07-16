# Installation

## Prerequisites

GCC/gfortran **10 or newer** (12+ recommended), plus `make`.

**Debian / Ubuntu / Linux Mint**

```bash
sudo apt update
sudo apt install gfortran gcc make
```

**Fedora / RHEL / Rocky Linux**

```bash
sudo dnf install gcc-gfortran gcc make
```

**openSUSE**

```bash
sudo zypper install gcc-fortran gcc make
```

**macOS (via Homebrew)**

```bash
brew install gcc
# gfortran ships bundled with gcc, e.g. as gfortran-14
```

Verify the version before continuing:

```bash
gfortran --version   # must be >= 10.0
```

## Building from source

```bash
# 1. Clone the repository
git clone https://github.com/mgimferrer/APOST3D.git
cd APOST3D

# 2. Set the installation path (add this to your shell profile too)
export APOST3D_PATH=$(pwd)

# 3. Build the bundled libxc-4.2.3 library (once)
bash compile_libxc.sh

# 4. Build apost3d, apost3d-eos, and eos_aom
bash make_compile.sh
```

`make_compile.sh` is the recommended entry point — it checks your gfortran
version up front, calls the `Makefile` for you, and (on macOS) ad-hoc
code-signs the binaries and smoke-tests that they actually launch, catching
the most common install problems immediately with a clear message instead of
a cryptic failure on your first real calculation.

Useful flags:

```bash
bash make_compile.sh --clean        # force a full rebuild
bash make_compile.sh --nthreads 8   # OMP_NUM_THREADS to use for test runs
bash make_compile.sh --help
```

If you'd rather drive `make` directly (custom build setups, CI, etc.):

```bash
make -C $APOST3D_PATH all      # build apost3d, apost3d-eos, eos_aom
make -C $APOST3D_PATH clean    # remove all objects and binaries
```

## Verify the install

```bash
make -C $APOST3D_PATH test
```

See [Running the test suite](testing.md) — a clean pass across all active
tests is the best confirmation your build is sound.
