# Installation

## Prerequisites

GCC/gfortran **10 or newer** (12+ recommended), `make`, and **OpenBLAS**
(BLAS/LAPACK — `diagonalize()` in `sources/util.f` uses LAPACK's `dsyevd`
for every diagonalization in the program). Free and open-source
throughout — no Intel compiler, no MKL, no license of any kind required,
matching the same reasoning behind the ifort → gfortran move itself.

libxc (exchange-correlation functionals — see the admonition below) is
handled automatically by `make_compile.sh`, but if it needs to build its
own copy from source, that build needs **CMake ≥ 3.21** too — install it
up front to avoid a mid-build stop.

**Debian / Ubuntu / Linux Mint**

```bash
sudo apt update
sudo apt install gfortran gcc make libopenblas-dev cmake
```

**Fedora / RHEL / Rocky Linux**

```bash
sudo dnf install gcc-gfortran gcc make openblas-devel cmake
```

**openSUSE**

```bash
sudo zypper install gcc-fortran gcc make openblas-devel cmake
```

**macOS (via Homebrew)**

```bash
brew install gcc openblas cmake
# gfortran ships bundled with gcc, e.g. as gfortran-14
```

```{admonition} CMake too old on your system?
:class: tip

Common on conservative HPC distros. `pip install --user cmake` (or
`pipx install cmake`) gets a modern prebuilt binary with nothing to
compile and no root needed — same spirit as the OpenBLAS/libxc
"lives somewhere nonstandard" escape hatches below.
```

```{admonition} macOS: openblas is keg-only
:class: note

Homebrew won't symlink `openblas` into the default search path on macOS,
since Apple's Accelerate framework already provides a system BLAS/LAPACK.
`make_compile.sh`/the `Makefile` auto-detect the Homebrew prefix via
`brew --prefix openblas` — no manual `LDFLAGS`/`CPPFLAGS` exports needed.
```

**HPC clusters**: OpenBLAS is close to universally available as a module
(`module load openblas` or similar) alongside or instead of vendor math
libraries — check `module avail` on your system. As long as `-lopenblas`
resolves (module-provided `LIBRARY_PATH`/`LD_LIBRARY_PATH` is normally
enough), no other setup is required.

```{admonition} OpenBLAS lives somewhere nonstandard?
:class: tip

Both `make_compile.sh` and the `Makefile` try, in order: an explicit
`OPENBLAS_DIR` you set yourself, then `pkg-config` (covers most
conda/spack/package-manager installs automatically, wherever they
actually are), then Homebrew's prefix on macOS, then a bare `-lopenblas`
relying on the default linker search path. If none of those resolve —
a custom-built copy in a one-off location, say — just point at it
directly and every step above is skipped:

    export OPENBLAS_DIR=/path/to/openblas   # expects lib/ and include/ under it
    bash make_compile.sh

`make_compile.sh` actually links a test program against whichever method
it picks before touching the rest of the build, so a misconfigured
`OPENBLAS_DIR` fails immediately with a clear message rather than deep
into compilation.
```

```{admonition} libxc: fetched and built automatically
:class: tip

Unlike OpenBLAS, libxc has no near-universal system package, so
`make_compile.sh` probes for an already-usable install first — an
explicit `LIBXC_DIR` you set yourself, then `pkg-config` (covers a
distro package, a conda environment, or Homebrew's `libxc` formula on
macOS), in that order — and only if none of those are found does it
fetch and build its own pinned copy under `libxc-<version>/` via CMake
(`compile_libxc.sh`, called automatically). That fetch step needs
network access to `gitlab.com` and verifies the download against a
checksum recorded in the script before building; on an air-gapped
machine, pre-download the matching tarball named in `compile_libxc.sh`
and place it at `$APOST3D_PATH/libxc-<version>.tar.gz` first.

Either way, this only happens once — every subsequent
`make_compile.sh`/`make apost3d` call reuses whatever was found or built
the first time, same as OpenBLAS above:

    export LIBXC_DIR=/path/to/libxc   # expects lib/ and include/ under it
    bash make_compile.sh
```

Verify the versions before continuing:

```bash
gfortran --version   # must be >= 10.0
```

## Building from source

```bash
# 1. Clone the repository
git clone https://github.com/mgimferrer/APOST3D.git
cd APOST3D

# 2. Set the installation path (add this to your shell profile too: it is
#    how your job scripts find $APOST3D_PATH/apost3d)
export APOST3D_PATH=$(pwd)

# 3. Build apost3d, apost3d-eos and the utilities in utils/
#    (fetches + builds libxc automatically first, if needed — see above)
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
bash make_compile.sh NTHREADS=4     # 4 compile jobs (default: 8)
bash make_compile.sh ARCH=native    # optimize for this machine's CPU only
bash make_compile.sh help           # list all available arguments
```

```{admonition} Same flag for building and testing
:class: tip

`NTHREADS=<n>` means the same thing and is spelled the same way whether
you're building (`bash make_compile.sh NTHREADS=4`, the number of parallel
compile jobs) or running the test suite (`make test NTHREADS=4`, the threads
of each test) — see [Running the test suite](testing.md). The default is 8,
or fewer if the machine has fewer CPUs. On a shared cluster login node, use a
small value.
```

### Which CPUs the build runs on (`ARCH`)

By default the code is compiled for the generic architecture (x86-64 or
arm64), so the same binary runs on every CPU of that family. This is the
right choice for a cluster whose nodes are of different ages, or when you
compile on a login node and run on compute nodes.

`ARCH=native` compiles for the CPU of the machine doing the build, using all
its instructions (AVX2, AVX-512, ...). That can be somewhat faster, but the
binary may stop with `Illegal instruction` on an older CPU, so use it only if
the program runs on the same kind of machine that compiled it. Any other
`gcc -march` value also works, e.g. `ARCH=x86-64-v3` for CPUs from about 2015
on. Changing `ARCH` triggers a full rebuild automatically.

### Utilities

`make_compile.sh` also builds the programs in `utils/`:

| Program | Purpose |
|---------|---------|
| `get_energy_g16`, `get_energy` | Append the reference energies of a Gaussian 16 / 09 `.log` to its `.fchk` (needed by the zero-error strategy of [ENPART](methods/enpart.md)) |
| `gen_hirsh` | Build the atomic densities file (`densoutput`) for Hirshfeld and Hirshfeld-I |
| `wfn2fchk` | Convert a `.wfn` file to `.fchk` |
| `group_frag` | Group ENPART energy terms by fragment |
| `eos_aom` | EOS from atomic overlap matrices of Multiwfn or AIMAll |
| `eos_alt` | EOS variants from an APOST-3D output |

`utils/apost3d.py` writes `.fchk` files from pySCF (see [ENPART](methods/enpart.md)).

### Using `make` directly

If you'd rather drive `make` directly (custom build setups, CI, etc.):

```bash
make -C $APOST3D_PATH -j8 all    # build apost3d and apost3d-eos
make -C $APOST3D_PATH -j8 utils  # build the utilities in utils/
make -C $APOST3D_PATH clean      # remove all objects and binaries
make -C $APOST3D_PATH help       # list all available targets and flags
```

Object files go to `objects/`; `apost3d` and `apost3d-eos` are written to
the repository root and the utilities to `utils/`.

## Verify the install

```bash
make -C $APOST3D_PATH test
```

See [Running the test suite](testing.md) — a clean pass across the whole
suite is the best confirmation your build is sound.
