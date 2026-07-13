<p align="center"><img width=25.0% src="https://github.com/mgimferrer/APOST3D/blob/master/media/logo-apost.png"></p>

## Chemical concepts from wave function analysis

A Fortran-based code developed at the Universitat de Girona (UdG) by P. Salvador and collaborators.

Builds with **GCC/gfortran** — free, open-source, and available on every Linux
distribution (no Intel compiler or license required).

## Shortcuts

* [Installation](#installation)
* [How to use](#how-to-use)
* [Running the test suite](#running-the-test-suite)
* [Troubleshooting](#troubleshooting)
* [Documentation](#documentation)
* [Cite the code](#citations)
* [Bug reports and feature requests](#bug-reports-and-feature-requests)

## Installation

### Prerequisites

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
# gfortran ships bundled with gcc, as e.g. gfortran-14
```

Verify the version before continuing:
```bash
gfortran --version   # must be >= 10.0
```

### Building from source

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

### Verify the install

```bash
make -C $APOST3D_PATH test
```

See [Running the test suite](#running-the-test-suite) below — a clean pass
across all active tests is the best confirmation your build is sound.

## How to use

```bash
export OMP_NUM_THREADS=1     # single-core for now; see note below
ulimit -s unlimited          # the code uses large stack-allocated arrays

$APOST3D_PATH/apost3d jobname > jobname.apost 2>&1
```

`jobname.fchk` and `jobname.inp` must be present in the working directory.
Correlated-wavefunction analyses (CASSCF, DMRG) also need `jobname.dm1`/`.dm2`
(1-/2-RDMs in the MO basis).

| File | Contents |
|---|---|
| `jobname.fchk` | Gaussian formatted checkpoint (wavefunction data) |
| `jobname.inp` | APOST-3D keyword input |
| `jobname.dm1` | 1-RDM in MO basis (CASSCF/DMRG only) |
| `jobname.dm2` | 2-RDM in MO basis (CASSCF/DMRG only) |

A detailed description of the input file format and all available keywords is
in the [Documentation](#documentation).

**On threads**: explicit OpenMP parallelization of the numerical integration
routines is planned but not yet implemented — `OMP_NUM_THREADS=1` is the
correct setting for now regardless of how many cores are available.

## Running the test suite

A regression test suite validates numerical output against reference values
for a handful of representative systems (RKS/UKS DFT, fragment analysis,
OSLO, EOS, QCHEM interface). It's the fastest way to confirm a build is
working correctly, and the main safety net when modifying the code.

```bash
make test          # fast tier (~2 min)
make test-full      # everything, including the slower C2H6-B3LYP case
```

```
════════════════════════════════════════════════════════════════
  APOST-3D Test Suite  ·  4 test(s)  ·  1 thread(s)
════════════════════════════════════════════════════════════════

  [ 1/4]  H2O-T-B3LYP                    dft enpart spin tfvc rks
           (10s)
           ✓  Normal Termination
           ✓  Total KS-DFT energy (au)              got -76.22411   ref -76.22411   Δ 0.0e+00
           ...
           PASSED  (10/10 checks)
  ...
════════════════════════════════════════════════════════════════
  ✓  H2O-T-B3LYP                          10s
  ✓  CH3F                                  0s
  ✓  FeCO2-PBEPBE                          0s
  ✓  FeO4-2                                2s

  4 PASSED   (12s total)
════════════════════════════════════════════════════════════════
```

| Command | Description |
|---|---|
| `make test` | Build + run the fast-tier tests (excludes tests tagged `slow`) |
| `make test-only` | Run fast-tier tests without rebuilding |
| `make test-full` | Build + run every test, including the `slow` tier |
| `make test TAGS=slow` | Run just the `slow` tier (e.g. `C2H6-B3LYP`) |
| `make test FILTER=H2O` | Run only tests whose name contains `H2O` |
| `make test TAGS=enpart` | Run only tests tagged `enpart` |
| `make test VERBOSE=1` | Show check details for passing tests too |
| `make test KEEP=1` | Save each test's raw `.apost` output to `tests/report/outputs/` |
| `make update-ref` | Regenerate reference outputs after an intentional code change |

Or invoke the runner directly for more options:
```bash
python3 tests/run_tests.py --help
```

### Active test cases

| System | Description | Tags |
|--------|-------------|------|
| `H2O-T-B3LYP` | Water, RKS B3LYP — TFVC, ENPART (DFT+IQA), local spin | `dft enpart spin tfvc rks` |
| `CH3F` | Fluoromethane, RKS DFT — TFVC, fragment OSLO | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | Iron dicarbonyl⁺, UKS PBE — TFVC, fragment EOS (open-shell) | `dft eos effao tfvc uks fragments openshell` |
| `FeO4-2` | Ferrate(VI)²⁻, UKS — TFVC, QCHEM interface, OSLO+EOS | `dft eos oslo tfvc uks fragments openshell qchem` |
| `C2H6-B3LYP` | Ethane, RKS B3LYP — full ENPART, THREBOD/MOD-GRIDTWOEL (~85s) | `... slow` — `make test-full` or `make test TAGS=slow` |

### Adding a new test case

1. Place `SystemName.fchk` and `SystemName.inp` in `compiler-testset/`.
2. Generate and sanity-check a reference run:
   ```bash
   cd compiler-testset && ulimit -s unlimited
   ../apost3d SystemName > SystemName.apost 2>&1   # confirm "Normal Termination"
   cp SystemName.apost ../tests/reference/SystemName.apost
   ```
3. Add an entry to `tests/manifest.json` (copy an existing similar test and
   adapt the tags, patterns, and reference values).
4. `python3 tests/run_tests.py --filter SystemName` until the checks pass.
5. Commit the input files, the reference output, and the manifest entry together.

## Troubleshooting

**`STOP The required input filename is missing`**
Normal — the program requires a job name argument (`apost3d jobname`).

**Segmentation fault on large jobs**
Run `ulimit -s unlimited` before launching; the code uses large stack-allocated arrays.

**`cannot find -lxcf90` or `-lxc`**
libxc wasn't built, or `APOST3D_PATH` isn't set. Re-run `bash compile_libxc.sh`
with `APOST3D_PATH` exported.

**Compiler version too old**
`gfortran --version` must report 10 or newer (`-fallow-argument-mismatch`
requires GCC 10+).

**macOS: binary fails to launch (`dyld`, "Library not loaded", killed on start)**
`make_compile.sh` ad-hoc code-signs and smoke-tests all three binaries
automatically. If it still fails, or you rebuilt without it:
```bash
xattr -cr $APOST3D_PATH
codesign --force --sign - $APOST3D_PATH/apost3d
codesign --force --sign - $APOST3D_PATH/apost3d-eos
codesign --force --sign - $APOST3D_PATH/eos_aom
```

## Documentation

The `APOST-3D` documentation is [here](DOCUMENTATION.md).

## Citations

### Cite the code

The following paper should be cited in publications utilizing `APOST-3D`:

* P. Salvador, E. Ramos-Cordoba, M. Montilla, L. Pujal and M. Gimferrer, *J. Chem. Phys.*, **2024**, 160, 172502
  DOI: [10.1063/5.0206187](https://doi.org/10.1063/5.0206187)

### Cite implemented methods

For atomic and overlap populations, bond orders and valences:

* I. Mayer and P. Salvador, *Chem. Phys. Lett.*, **2004**, 383, 368-375
  DOI: [10.1016/j.cplett.2003.11.048](https://doi.org/10.1016/j.cplett.2003.11.048)

For Hartree-Fock molecular energy decomposition:

* P. Salvador, M. Duran and I.Mayer, *J. Chem. Phys.*, **2001**, 115, 1153-1157
  DOI: [10.1063/1.1381407](https://doi.org/10.1063/1.1381407)
* P. Salvador and I. Mayer, *J. Chem. Phys.*, **2004**, 120, 5046-5052
  DOI: [10.1063/1.1646354](https://doi.org/10.1063/1.1646354)

For KS-DFT molecular energy decomposition:

* P. Salvador and I. Mayer, *J. Chem. Phys.*, **2007**, 126, 234113
  DOI: [10.1063/1.2741258](https://doi.org/10.1063/1.2741258)
* M. Gimferrer and P. Salvador, *J. Chem. Phys.*, **2023**, 158, 234105
  DOI: [10.1063/5.0142778](https://doi.org/10.1063/5.0142778)

For CAS/DMRG molecular energy decomposition:

For effective atomic/fragment orbitals:

* I. Mayer, *J. Phys. Chem.*, **1996**, 100, 6249
  DOI: [10.1021/jp952779i](https://doi.org/10.1021/jp952779i)
* I. Mayer and P. Salvador, *J. Chem. Phys.*, **2009**, 130, 234106
  DOI: [10.1063/1.3153482](https://doi.org/10.1063/1.3153482)
* E. Ramos-Cordoba, P. Salvador and I. Mayer, *J. Chem. Phys.*, **2013**, 138, 214107
  DOI: [10.1063/1.4807775](https://doi.org/10.1063/1.4807775)

For local spin analysis:

* E. Ramos-Cordoba, E. Matito, I. Mayer and P. Salvador, *J. Chem. Theor. Comput.*, **2012**, 8, 1270-1279
  DOI: [10.1021/ct300050c](https://doi.org/10.1021/ct300050c)
* E. Ramos-Cordoba, E. Matito, P. Salvador and I. Mayer, *Phys. Chem. Chem. Phys.*, **2012**, 14, 15291-15298
  DOI: [10.1039/C2CP42513K](https://doi.org/10.1039/C2CP42513K)

For effective oxidation states analysis:

* E. Ramos-Cordoba, V. Postils and P. Salvador, *J. Chem. Theor. Comput.*, **2015**, 11, 1501-1508
  DOI: [10.1021/ct501088v](https://doi.org/10.1021/ct501088v)
* M. Gimferrer and P. Salvador, _submitted_, **2024**
  DOI: [XX](XX)

For oxidation states from localized orbitals:

* M. Gimferrer, G. Comas-Vila and P. Salvador, *Molecules*, **2020**, 25, 234
  DOI: [10.3390/molecules25010234](https://doi.org/10.3390/molecules25010234)
* M. Gimferrer, J. Van der Mynsbrugge, A. T. Bell, P. Salvador and M. Head-Gordon *Inorg. Chem.*, **2020**, 59, 15410-15420
  DOI: [10.1021/acs.inorgchem.0c02405](https://doi.org/10.1021/acs.inorgchem.0c02405)
* M. Gimferrer, A. Aldossary, P. Salvador and M. Head-Gordon, *J. Chem. Theor. Comput.*, **2022**, 18, 309-322
  DOI: [10.1021/acs.jctc.1c01011](https://doi.org/10.1021/acs.jctc.1c01011)

For decomposition of EDA quantities into one- and two-center IQA terms:

* M. Gimferrer, S. Danes, D. M. Andrada and P. Salvador, *J. Chem. Theory Comput.*, **2023**, 19, 3469-3485
  DOI: [10.1021/acs.jctc.3c00143](https://doi.org/10.1021/acs.jctc.3c00143)

For origin-independent decomposition of static polarizabilities:

* M. Montilla, J. M. Luis and P. Salvador, *J. Chem. Theor. Comput.*, **2021**, 17, 1098-1105
  DOI: [10.1021/acs.jctc.0c00926](https://doi.org/10.1021/acs.jctc.0c00926)

## Bug reports and feature requests

Please submit tickets on the [issues](https://github.com/mgimferrer/APOST3D/issues) page, and/or send an email to mgimferrer18@gmail.com and pedro.salvador@udg.edu
