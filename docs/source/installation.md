# Installation

## Prerequisites

GCC/gfortran **10 or newer** (12+ recommended), `make`, and **OpenBLAS**
(BLAS/LAPACK — `diagonalize()` in `sources/util.f` uses LAPACK's `dsyevd`
for every diagonalization in the program). Free and open-source
throughout — no Intel compiler, no MKL, no license of any kind required,
matching the same reasoning behind the ifort → gfortran move itself.

**Debian / Ubuntu / Linux Mint**

```bash
sudo apt update
sudo apt install gfortran gcc make libopenblas-dev
```

**Fedora / RHEL / Rocky Linux**

```bash
sudo dnf install gcc-gfortran gcc make openblas-devel
```

**openSUSE**

```bash
sudo zypper install gcc-fortran gcc make openblas-devel
```

**macOS (via Homebrew)**

```bash
brew install gcc openblas
# gfortran ships bundled with gcc, e.g. as gfortran-14
```

```{admonition} macOS: openblas is keg-only
:class: note

Homebrew won't symlink `openblas` into the default search path on macOS,
since Apple's Accelerate framework already provides a system BLAS/LAPACK.
`make_compile.sh`/the `Makefile` auto-detect the Homebrew prefix via
`brew --prefix openblas` and link against that directly — no manual
`LDFLAGS`/`CPPFLAGS` exports needed.
```

**HPC clusters**: OpenBLAS is close to universally available as a module
(`module load openblas` or similar) alongside or instead of vendor math
libraries — check `module avail` on your system. As long as `-lopenblas`
resolves (module-provided `LIBRARY_PATH`/`LD_LIBRARY_PATH` is normally
enough), no other setup is required.

Verify the versions before continuing:

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

`make_compile.sh` takes the same arguments a `make` invocation would — bare
words for actions, `KEY=value` for parameters — rather than its own set of
dashed flags:

```bash
bash make_compile.sh clean          # force a full rebuild
bash make_compile.sh NTHREADS=8     # OMP_NUM_THREADS to report/export
bash make_compile.sh clean NTHREADS=8
bash make_compile.sh help           # list all available arguments
```

```{admonition} Same flag for building and testing
:class: tip

`NTHREADS=<n>` means the same thing and is spelled the same way whether
you're building (`bash make_compile.sh NTHREADS=4`) or running the test
suite (`make test NTHREADS=4`) — see [Running the test suite](testing.md).
```

If you'd rather drive `make` directly (custom build setups, CI, etc.):

```bash
make -C $APOST3D_PATH all      # build apost3d, apost3d-eos, eos_aom
make -C $APOST3D_PATH clean    # remove all objects and binaries
make -C $APOST3D_PATH help     # list all available targets and flags
```

## Verify the install

```bash
make -C $APOST3D_PATH test
```

See [Running the test suite](testing.md) — a clean pass across the whole
suite is the best confirmation your build is sound.
