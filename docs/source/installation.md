# Installation

Installing APOST-3D takes three steps: install a few standard packages,
run one build script, and run the test suite to check the result. Everything
used is free and open source; no compiler or library license is needed.

## 1. Install the prerequisites

You need a Fortran compiler (**gfortran 10 or newer**), `make`, the
**OpenBLAS** linear algebra library, **CMake** and **Python 3** (for the
test suite). On most systems one command installs them all:

**Debian, Ubuntu, Linux Mint**

```bash
sudo apt update
sudo apt install gfortran gcc make libopenblas-dev cmake python3 curl
```

**Fedora, RHEL, Rocky Linux, AlmaLinux**

```bash
sudo dnf install gcc-gfortran gcc make openblas-devel cmake python3 curl
```

**openSUSE**

```bash
sudo zypper install gcc-fortran gcc make openblas-devel cmake python3 curl
```

**macOS** (with [Homebrew](https://brew.sh))

```bash
xcode-select --install       # command-line tools (make, python3), if not yet installed
brew install gcc openblas cmake
```

**HPC cluster** (no administrator rights): load the compiler and OpenBLAS
modules, for example

```bash
module load gcc openblas      # names differ between clusters: see `module avail`
```

and, if the cluster's CMake is older than 3.21, get a recent one with
`pip install --user cmake`.

Check the compiler version:

```bash
gfortran --version            # must say 10 or higher
```

The exchange-correlation library [libxc](https://libxc.gitlab.io/) is also
needed, but you don't have to install it: the build script downloads and
compiles it the first time (this is what CMake is for).

## 2. Build the program

```bash
git clone https://github.com/mgimferrer/APOST3D.git
cd APOST3D
bash make_compile.sh
```

The first build takes a few minutes, most of it compiling libxc; later
builds take seconds. If something is missing (CMake, OpenBLAS, a recent
enough gfortran), the script stops with a clear message and the command
to fix it (see [Troubleshooting](troubleshooting.md)). At the end it
launches each program once and lists them:

```text
--- Binaries produced ---
  ✓  /home/user/APOST3D/apost3d  (...)
...
--- Smoke test (launching each binary) ---
  ✓  apost3d launches correctly
...
  Build complete: ...
```

The programs are written into the repository folder: `apost3d` (the main
program) and the [utilities](tools/utilities.md) in `utils/`.

Finally, tell your shell where APOST-3D is, so that job scripts can find it.
Add this line to your `~/.bashrc` (or `~/.zshrc` on macOS), replacing the
path with the folder you cloned:

```bash
export APOST3D_PATH=/home/user/APOST3D
```

Open a new terminal (or run `source ~/.bashrc`) for it to take effect.

## 3. Check the build

```bash
cd $APOST3D_PATH
make test
```

This runs 17 complete calculations and compares every number they print
with stored reference outputs. It takes a few minutes and should end with

```text
  17 PASSED   (...s total)
```

If a test fails, the output says which values differ; see
[Running the test suite](developer/testing.md) for how to read it.

You are ready: continue with the [tutorial](tutorial.md).

## Build options

`make_compile.sh` accepts these options (the same names work with `make`):

| Option | Effect |
|---|---|
| `NTHREADS=<n>` | Number of parallel compile jobs. Default 8, or fewer if the machine has fewer CPUs; use a small value on a shared login node. `make test` takes the same option for its threads (`make test NTHREADS=4`). |
| `ARCH=native` | Optimize for the CPU of the machine that compiles. By default the program runs on every CPU of its family (x86-64 or arm64), which is what you want on a cluster with nodes of different ages, or when compiling on a login node. With `ARCH=native` it may be somewhat faster, but can stop with `Illegal instruction` on older CPUs. Any other `gcc -march` value also works, e.g. `ARCH=x86-64-v3`. |
| `clean` | Remove everything compiled and rebuild from scratch. |
| `help` | List all options. |

For example `bash make_compile.sh NTHREADS=4`. Changing the compiler or
`ARCH` triggers a full rebuild automatically.

## Updating

```bash
cd $APOST3D_PATH
git pull
bash make_compile.sh
make test
```

## Special setups

You only need this section if the build script could not find something
by itself.

**OpenBLAS in a nonstandard location.** The script looks for OpenBLAS in
this order: the `OPENBLAS_DIR` variable, `pkg-config` (which covers most
conda, spack and package-manager installations), Homebrew, and finally
the default linker path. To use a specific copy, point at the folder that
contains its `lib/` and `include/`:

```bash
export OPENBLAS_DIR=/path/to/openblas
bash make_compile.sh
```

**An existing libxc.** If libxc is already installed (a distribution
package, conda, Homebrew's `libxc`), the script finds it through
`pkg-config` and does not build its own. It must be libxc 7, built with the
same gfortran; if an installed copy gives build errors, build the bundled
one with `bash compile_libxc.sh` and point at it with
`LIBXC_DIR=$APOST3D_PATH/libxc-7.1.2`. To use a specific copy:

```bash
export LIBXC_DIR=/path/to/libxc       # the folder with lib/ and include/
bash make_compile.sh
```

**A machine without internet access.** The script downloads libxc from
`gitlab.com`. On a machine that cannot, download
<https://gitlab.com/libxc/libxc/-/archive/7.1.2/libxc-7.1.2.tar.gz>
elsewhere, copy it into the APOST-3D folder as `libxc-7.1.2.tar.gz`, and
run `bash make_compile.sh`: the file is used instead of downloading (and
checked against the same checksum).

**CMake too old.** libxc needs CMake 3.21 or newer. `pip install --user
cmake` installs a recent one without administrator rights.

**macOS.** Homebrew's OpenBLAS is found automatically. The script also
signs the programs so that macOS lets them run; if one is still blocked,
see [Troubleshooting](troubleshooting.md).

To compile with `make` directly instead of the script, see
[Building with make](developer/building.md).
