# CLAUDE.md — APOST-3D Developer Reference

This file is the primary reference for AI-assisted development of APOST-3D. It documents architecture, conventions, known issues, and planned work. Keep it updated as the codebase evolves.

---

## Project Overview

**APOST-3D** (Version 4, last updated March 2025) is an open-source Fortran program for wave function analysis using both real-space and Hilbert-space approaches. It was developed at the Universitat de Girona (UdG) by P. Salvador and collaborators.

- **Language**: Fortran (fixed-form `.f`, free-form `.f90`)
- **Size**: ~100k+ lines of code across ~25 source files
- **Repository**: https://github.com/mgimferrer/APOST3D
- **Primary paper**: P. Salvador, E. Ramos-Cordoba, M. Montilla, L. Pujal and M. Gimferrer, *J. Chem. Phys.*, 2024, 160, 172502. DOI: 10.1063/5.0206187
- **Contact**: mgimferrer18@gmail.com, psalse@gmail.com, eloy.raco@gmail.com

---

## Repository Structure

```
APOST3D/
├── sources/            # All Fortran source files
├── objects/            # Compiled object files (.o) — not tracked in git
├── utils/              # Auxiliary standalone programs and scripts
├── lebedev/            # Lebedev-Laikov angular quadrature routines
├── libxc-4.2.3/        # Bundled libxc library (version 4.2.3, pinned; rebuilt by compile_libxc.sh)
├── libxc-4.2.3.tar.gz  # Tarball for re-extraction
├── compiler-testset/   # Input files (.fchk, .inp, .dm1, .dm2) for regression tests
├── tests/              # Regression test suite (runner, manifest, reference outputs)
├── media/              # Logo and images for GitHub/documentation
├── Makefile            # gfortran build (Phase 0)
├── compile_libxc.sh    # Script to build the bundled libxc library with gfortran
├── make_compile.sh     # Master compilation script (single-step gfortran build)
├── README.md           # GitHub README
├── README-gfortran.md  # Full build and installation guide
├── DOCUMENTATION.md    # User documentation (input keywords, output format)
└── LICENSE
```

---

## Source Files (`sources/`)

| File | Lines | Role |
|------|-------|------|
| `modules.f90` | 566 | F90 modules: `basis_set`, `ao_matrices`, `integration_grid`. Defines all global data structures (basis, MO coefficients, density matrices, integration grid). |
| `main.f` | 1582 | Program entry point. Reads arguments, processes `.inp` and `.fchk` files, dispatches to analysis routines. Contains the `kiir` startup banner subroutine and `readchar`/`readint` input parsing calls. |
| `input2.f` | 1139 | `input` subroutine: reads the `.fchk` file (Gaussian formatted checkpoint), builds the basis set, reads MO coefficients and density matrices for all wavefunction types. Also contains `readchar`, `readint`, `int_locate`, `real_locate` helper functions for keyword-based `.inp` parsing. |
| `enpart.f` | 2770 | Energy partitioning (IQA/ENPART): one- and two-electron energy decomposition into atomic and diatomic contributions for HF, KS-DFT, and CAS wavefunctions. Contains `numint_one`, `numint_two` (RHF), `numint_one_uhf`, `numint_two_uhf` (UHF) subroutines plus polar/multipolar analysis. |
| `enpart_dft.f` | 1225 | DFT exchange-correlation energy decomposition routines. Works with libxc. |
| `enpart_phf.f` | 919 | Energy partitioning for post-Hartree-Fock (PHF) wavefunctions (CAS, CISD, CCSD). |
| `effao.f` | 1521 | Effective Atomic Orbitals (EFAOs) from 3D integration: `ueffao3d_frag` (fragment-based), `ueffao3d` (atom-based), `ueffaolow_frag`, `ueffaomull_frag`, `uefomo`, `eos_analysis`. Core of the EOS (Effective Oxidation States) analysis. |
| `oslo.f` | 2017 | Orbital Symmetry-adapted Localized Orbitals (OSLOs): `rwf_iterative_oslo` (RHF) and `uwf_iterative_oslo` (UHF). Fragment orbital decomposition and population analysis. |
| `numint.f` | 555 | Numerical integration core: `prenumint` builds the grid, evaluates AO values at grid points (`fpoints`), assembles integration weights (`rpoints`), evaluates electron density. |
| `wat.f` | 1476 | Atomic weight functions for real-space partitioning schemes: Becke, Hirshfeld, iterative Hirshfeld (Hirshfeld-I), TFVC. Contains `wat`, `wathirsh`, `wathirshit3`. |
| `qtaim.f` | 1398 | QTAIM (Quantum Theory of Atoms in Molecules) integration via gradient-path topology. Nuclear critical point search and basin integration. |
| `pop.f` | 555 | Atomic/overlap populations and bond orders: `fborder` (Mayer bond order / overlap population matrix), `tomull`, `tolow` (Mulliken, Löwdin population matrices). |
| `mulliken.f` | 292 | Mulliken population analysis routines. |
| `loba.f` | 215 | LOBA (Localized Orbital Bonding Analysis): `eos_loba` assigns oxidation states from localized orbital fragment populations using the Clarity Index (CI). |
| `ueos.f` | 551 | UEOS (Unpaired EOS): `effao3d_u` computes EFAOs from paired and unpaired densities (Takatsuka definition), enabling EOS analysis for open-shell systems. |
| `corr.f` | 285 | `spincorr` subroutine: spin correlation matrices and local spin analysis for correlated wavefunctions using 1-RDM and 2-RDM. |
| `print.f` | 2044 | All formatted output routines: `VPRINT` (vector printing), `MPRINT_NLOP`, and many result printers for populations, energies, EOS, OSLO outputs. |
| `util.f` | 1383 | General utility subroutines: linear algebra (`diagonalize`, `build_Smp` for S^±1/2), integral helpers, DIIS, file I/O helpers, `kiir` banner. |
| `quad.f` | 106 | Interface to Lebedev-Laikov quadrature. |
| `top_iso.f` | 685 | Topology of isosurfaces and orbital topology analysis: `top_3d` subroutine for computing exchange/correlation topology functions on a real-space grid. |
| `subroutines_mmo.f` | 904 | Additional subroutines contributed by M. Gimferrer. |
| `scatt_fact.f` | 146 | X-ray scattering factor computation. |
| `devel.f` | 1292 | Development code (deprecated routines, old implementations being replaced). Contains notes between developers. |
| `parameter.h` | 33 | Global PARAMETER definitions (array size limits, constants). Included in all source files. |

### Files to be resolved

- `devel.f` — contains old/deprecated code and developer notes; should be cleaned up
- `modules.f90.to_do` — future module restructuring work
- `dafh.f.to_do` — Domain-based analysis of Fermi holes (planned)

---

## Global Parameters (`sources/parameter.h`)

```fortran
nmax   = 8000     ! Max number of basis functions
maxat  = 350      ! Max number of atoms
maxp   = 10000    ! Max number of Cartesian primitives
maxg   = nmax
maxc   = 36       ! Max primitives per basis function
maxnna = 20       ! Max number of NNAs (non-nuclear attractors)
maxfrag = maxat   ! Max fragments = max atoms
maxgrid = 10000   ! Max grid size for cubegen
thresh  = 1e-8    ! General numerical threshold
pi      = 3.14159265358979d0
angtoau = 0.52917721067d0   ! Angstrom to Bohr
tokcal  = 627.5096d0        ! a.u. to kcal/mol
```

---

## F90 Modules (`sources/modules.f90`)

Three modules hold all global state. They replace many of the old `COMMON` blocks (migration is ongoing — many `COMMON` blocks still exist in the `.f` files).

### `basis_set`
Stores the Gaussian basis set: primitive exponents, contraction coefficients, angular momentum indices, atom-to-basis maps, overlap matrix, and S^±1/2 matrices. The `build_basis()` subroutine reads all of this from the `.fchk` file and supports angular momenta up to g-type (l=4). Handles both Cartesian and pure (spherical) basis sets.

### `ao_matrices`
Stores MO coefficient matrices (`c`, `cb` for alpha/beta), density matrices (`p`, `pa`, `pb`, `ps`), and natural orbital matrices (`c_no`, `occ_no`). Allocated dynamically in `build_ao_matrices(igr)`.

### `integration_grid`
Stores numerical integration parameters: `Nrad` (radial points, default 40 or 150), `Nang` (angular points from Lebedev grid, default 146 or 590), grid rotation angles (`pha`, `phb`), and the actual grid arrays (`th`, `ph`, `w`, `wr`, `xr`). Lebedev grid orders are tabulated in the `leved` array (32 levels from 6 to 5810 points).

---

## COMMON Blocks (legacy global state in `.f` files)

The main program and most subroutines communicate through a set of named COMMON blocks. These are declared redundantly in every file that uses them — a major source of inconsistency to address in the MAJOR-UPDATE.

| COMMON block | Contents |
|---|---|
| `/nat/` | `nat` (atoms), `igr` (basis functions), `ifg`, `nocc`, `nalf`, `nb`, `kop` |
| `/cas/` | `icas`, `ncasel`, `ncasorb`, `nspinorb`, `norb`, `icisd`, `icass` |
| `/coord/` | `coord2(3,maxat)`, `zn(maxat)`, `iznuc(maxat)` |
| `/iops/` | `iopt(100)` — the master option array (all user-selected options stored here) |
| `/atlist/` | `iatlist(maxat)`, `icuat` — atom list for selected atoms |
| `/frlist/` | `ifrlist`, `nfrlist`, `icufr`, `jfrlist` — fragment lists |
| `/loba/` | `oxi`, `errsav`, `elec`, `effpop` — LOBA/EOS results |
| `/qat/` | `qat(maxat,2)`, `qsat(maxat,2)` — atomic charges |
| `/ovpop/` | `op`, `bo`, `di` (overlap populations, bond orders, delocalization indices), `totq` |
| `/localspin/` | `xlsa(maxat,maxat)`, `ua(maxat)` — local spin matrices |
| `/exchg/` | `exch(maxat,maxat)`, `xmix` — exchange-correlation contributions |
| `/modgrid/` | `nrad22`, `nang22`, `rr0022`, `phb12`, `phb22` — secondary grid for two-electron integrals |
| `/edaiqa/` | `xen`, `xcoul`, `xnn` — EDA-IQA energy terms |
| `/efield/` | `field(4)`, `edipole` — external electric field |
| `/printout/` | `iaccur` — output precision flag |
| `/filename/` | `name0` — base name of the input files |

---

## The `iopt(100)` Option Array

The input parser stores all user-selected keywords as integers in `iopt(100)`. Key entries (from `main.f`):

| Index | Keyword | Meaning |
|-------|---------|---------|
| 5 | MULLI | Mulliken/Löwdin/Davidson populations |
| 6 | HIRSH | Hirshfeld scheme |
| 7 | ALLPOINTS | Use all grid points |
| 9 | DENS | Which density to use from fchk |
| 12 | EFFAO | Effective AO analysis |
| 13 | CUBE | Output cube files |
| 14 | IBCP | Becke with rho correction |
| 16 | QTAIM | QTAIM integration |
| 24 | EFF_THRESH | Occupation threshold for EFAOs |
| 34 | LAPLACIAN | Compute Laplacian |
| 40 | DOFRAGS | Fragment analysis |
| 43 | ORCA | ORCA wavefunction interface |
| 47 | FIELD | External field |
| 85 | PAIRS | Orbital pair analysis (topology) |
| 86 | ETOP | Energy topology |
| 87 | PYSCF | PySCF wavefunction interface |
| 88 | SPINSEP | Spin separation |

---

## Input/Output Format

### Running the program

```bash
export OMP_NUM_THREADS=1          # Phase 0: single core
ulimit -s unlimited               # prevent stack overflow
$APOST3D_PATH/apost3d name-input > name-output.apost
```

Requires two files in the working directory:
- `name-input.fchk` — Gaussian formatted checkpoint file (wavefunction data)
- `name-input.inp` — APOST-3D keyword input file

For correlated calculations (CASSCF, DMRG), optionally:
- `name-input.dm1` — 1-RDM in MO basis
- `name-input.dm2` — 2-RDM in MO basis

### Input file format (`.inp`)

Section-based, keyword-driven. Sections start with `# SECTIONNAME`. Example:

```
# METHOD
TFVC
SPIN
ENPART
#
# ENPART
B3LYP
THREBOD -1
MOD-GRIDTWOEL
#
# GRID
RADIAL 40
ANGULAR 146
phb1 0.162
phb2 0.182
#
```

Key `# METHOD` keywords: `MULLI`, `LOWDIN`, `HIRSH`, `HIRSH-IT`, `BECKE-RHO`, `TFVC`, `QTAIM`, `EFFAO`, `EOS`, `EOS-U`, `OSLO`, `LOBA`, `SPIN`, `ENPART`, `EDAIQA`, `POLAR`, `TOPOLOGY`, `DM`, `DOATOMS`, `DOFRAGS`, `SCATT-FACT`, `QCHEM`, `MOKIT`, `WFN`.

### Supported wavefunction interfaces (via `.fchk`)
- **Gaussian**: native `.fchk` format
- **ORCA** (`ORCA` keyword): ORCA-generated `.fchk`
- **Q-Chem** (`QCHEM` keyword): Q-Chem-generated `.fchk`
- **MOKIT** (`MOKIT` keyword): MOKIT-generated `.fchk`
- **PySCF** (`pySCF` in `# DM` section): PySCF 1-RDM/2-RDM
- **ORCA** (`ORCA` in `# DM` section): ORCA-generated RDMs

---

## Analysis Methods

### Population Analysis (Hilbert space)
- **Mulliken** (`MULLI`) — standard Mulliken atomic populations, overlap populations, bond orders (Mayer), delocalization indices
- **Löwdin** (`LOWDIN`) — Löwdin symmetric orthogonalization
- **Davidson-Löwdin** (`LOWDIN` with value 3) — Davidson's modification

### Atomic Partitioning (real space)
- **Becke** (`BECKE-RHO`) — Becke fuzzy-cell scheme, J. Chem. Phys. 88, 2547 (1988)
- **Hirshfeld** (`HIRSH`) — stockholder partitioning, Theor. Chim. Acta 44, 129 (1977)
- **Hirshfeld-Iterative** (`HIRSH-IT`) — iterative Hirshfeld, J. Chem. Phys. 126, 144111 (2007)
- **TFVC** (`TFVC`) — Topological Fuzzy Voronoi Cells, J. Chem. Phys. 139, 071103 (2013)
- **QTAIM** (`QTAIM`) — quantum theory of atoms in molecules via gradient-path topology

### Orbital Analysis
- **EFFAO** (`EFFAO`) — Effective Atomic Orbitals from 3D integration (Mayer-Salvador-Ramos-Cordoba)
- **EOS** (`EOS`) — Effective Oxidation States from EFAOs (Ramos-Cordoba-Salvador)
- **OSLO** (`OSLO`) — Orbital Symmetry-adapted Localized Orbitals, fragment-based
- **LOBA** (`LOBA`) — Localized Orbital Bonding Analysis (Gimferrer-Salvador)
- **EOS-U** (`EOS-U`) / **UEOS** — EOS for open-shell systems from unpaired density

### Energy Partitioning (IQA-like)
- **ENPART** (`ENPART`) — molecular energy partitioning into one- and two-center terms
  - HF (`HF`), LDA, BP86, B3LYP and other DFT functionals via libxc
  - CAS/DMRG (`CASSCF`, `DM=2` with external RDMs)
  - Supports CISD and CCSD
- **EDAIQA** (`EDAIQA`) — EDA-IQA: decomposition of EDA interaction energies into IQA one- and two-center contributions

### Other
- **SPIN** (`SPIN`) — local spin analysis (Ramos-Cordoba et al.)
- **POLAR** (`POLAR`) — origin-independent static polarizability decomposition (Montilla-Luis-Salvador)
- **TOPOLOGY** (`TOPOLOGY`) — topology analysis of exchange/correlation functions
- **SCATT-FACT** (`SCATT-FACT`) — X-ray scattering factors
- **DAFH** — Domain-based analysis of Fermi holes (in development, `dafh.f.to_do`)

---

## Numerical Integration

The integration engine uses an atom-centered numerical grid:

- **Radial grids**: Gauss-Chebyshev type, `Nrad` points (default 40 for populations, 150 for ENPART)
- **Angular grids**: Lebedev-Laikov spherical grids, `Nang` points (default 146 for populations, 590 for ENPART). 32 levels from 6 to 5810 points, tabulated in `leved` array.
- **Grid rotation**: Small random rotations (`pha`, `phb`) applied to avoid systematic errors at special angles (typically `phb1=0.162`, `phb2=0.182`)
- **Two-electron integrals** (`MOD-GRIDTWOEL`): uses a second, potentially different grid for the two-electron part of ENPART (controlled via `nrad22`, `nang22`)
- **FINEGRID**: activates `nang=974` instead of 590 for higher-accuracy one-electron integrations

The main numerical flow:
1. `build_integration_grid()` — set up grid parameters
2. `prenumint()` — compute grid coordinates, AO values at all points, electron density, and atomic weight functions
3. `numint_one()` / `numint_two()` — contract grid quantities into atomic/diatomic properties

---

## Compilation System

### Phase 0: gfortran build (current, MAJOR-UPDATE branch)

The code now builds with GCC/gfortran — no Intel toolchain required.
See `README-gfortran.md` for full installation instructions.

**Prerequisites**
- GCC/gfortran ≥ 10 (GCC 12+ recommended; GCC 15 tested)
  - macOS: `brew install gcc`
  - Linux: `sudo apt install gfortran gcc make`
- Environment variable `APOST3D_PATH` must be set

**Build sequence**
```bash
export APOST3D_PATH=/path/to/APOST3D

# 1. Build libxc-4.2.3 (once)
bash $APOST3D_PATH/compile_libxc.sh

# 2. Build APOST-3D
bash $APOST3D_PATH/make_compile.sh          # or with --clean
```

**Build targets** (`Makefile`)
- `apost3d` — main executable
- `apost3d-eos` — standalone EOS executable (subset of sources)
- `eos_aom` — utility executable (`utils/eos_aom.f90`)

**Key compiler flags**
| Flag | Reason |
|------|--------|
| `-O3 -ffast-math` | High optimisation (matches old `-Ofast`) |
| `-march=native` | Tune for build CPU; remove for cross-compilation |
| `-fopenmp` | Required by `use OMP_LIB`; code runs single-threaded in Phase 0 |
| `-fbacktrace` | Stack trace on runtime errors |
| `-ffixed-line-length-132` | Allow legacy fixed-form lines up to 132 chars |
| `-fallow-argument-mismatch` | Accept implicit-interface rank mismatches (ifort silently tolerated these; required for UHF `chp2b` scalar/array mismatch in RHF call paths never reached at runtime) |
| `-O1` for `input2.f` | Avoid floating-point parsing issues at high optimisation |

**Parallelization (Phase 0)**
- `OMP_NUM_THREADS=1` (single-core; explicit OpenMP directives planned for Phase 1)
- Always run with `ulimit -s unlimited` to prevent stack overflow

**macOS note**
On macOS 26 (Tahoe) and later, GCC-compiled binaries require an ad-hoc code
signature to access the dyld shared cache. `make_compile.sh`
runs `codesign --force --sign -` automatically after the build.

---

## Utilities (`utils/`)

| File | Description |
|------|-------------|
| `apost3d.py` | Python wrapper script |
| `main_eos.f` | Main program for the standalone `apost3d-eos` executable |
| `eos_aom.f90` | Standalone EOS-AOM (atom/orbital mapping) utility |
| `eos_alt.f90` | Alternative EOS calculation |
| `gen_hirsh.f` | Generate Hirshfeld promolecular densities |
| `get_energy.f` | Extract energy values from output files |
| `get_energy_g16.f` | Same, for Gaussian 16 outputs |
| `group_frag.f` | Fragment grouping utility |
| `wfn2fchk.f90` | Convert `.wfn` files to `.fchk` format |
| `memchk.f` | Memory checking utility |
| `insert.sh` | Shell script helper |

---

## Test Suite

### Overview

The regression test suite lives in `tests/` and is driven by a Python runner
(`tests/run_tests.py`) that reads a structured manifest (`tests/manifest.json`).
It runs the `apost3d` binary on each test case, extracts named numerical
quantities from the output using regex patterns, and compares them against
reference values with configurable tolerances.

```
tests/
├── run_tests.py      # Python test runner (pure stdlib, no pip install needed)
├── manifest.json     # Registry of all test cases and their checks
├── reference/        # Validated reference .apost output files
│   ├── H2O-T-B3LYP.apost
│   ├── CH3F.apost
│   ├── FeCO2-PBEPBE.apost
│   └── FeO4-2.apost
└── report/           # Auto-generated after each run (gitignored)
    ├── last_run.txt
    └── last_run.html
```

Input files (`.fchk`, `.inp`, auxiliary `.fchk`) remain in `compiler-testset/`.

### Running the tests

```bash
# Full build + fast-tier tests (recommended after any code change):
make test

# Build + tests with more threads:
make test TEST_NTHREADS=4

# Run without rebuilding:
make test-only

# Everything, including the 'slow' tier (e.g. C2H6-B3LYP):
make test-full

# Just the slow tier:
make test TAGS=slow

# Run a single test by name:
make test FILTER=H2O

# Run all tests tagged "enpart":
make test TAGS=enpart

# Verbose output (show check details even for passing checks):
make test VERBOSE=1

# Invoke the runner directly with full options:
python3 tests/run_tests.py --filter CH3F --verbose
python3 tests/run_tests.py --tags oslo,eos
python3 tests/run_tests.py --exclude-tags slow
python3 tests/run_tests.py --no-color 2>&1 | tee test.log
```

`make test`/`make test-only` exclude tests tagged `slow` by default (via
`--exclude-tags slow`), *unless* `TAGS=...` is given explicitly — so
`make test TAGS=slow` runs just the slow tier, not fast+slow. `make
test-full` always runs everything, unfiltered.

The runner exits with code 0 if all tests pass, 1 if any fail — suitable for
CI systems (GitHub Actions, etc.).

### Active test cases

| System | Wavefunction | Features tested | Tags |
|--------|-------------|-----------------|------|
| `H2O-T-B3LYP` | RKS B3LYP | TFVC, ENPART (DFT + IQA matrix), SPIN | `dft enpart spin tfvc rks` |
| `CH3F` | RKS DFT | TFVC, fragment OSLO analysis | `dft oslo tfvc rks fragments` |
| `FeCO2-PBEPBE` | UKS PBE | TFVC, fragment EOS (open-shell) | `dft eos effao tfvc uks fragments openshell` |
| `FeO4-2` | UKS, Q-Chem | TFVC, QCHEM interface, OSLO+EOS (open-shell) | `dft eos oslo tfvc uks fragments openshell qchem` |
| `C2H6-B3LYP` | RKS B3LYP | TFVC, full ENPART (DFT+IQA), THREBOD/MOD-GRIDTWOEL, 8 atoms/160 basis fns | `dft enpart tfvc rks threbod slow` — **slow tier**, excluded from default `make test` (~85s single-threaded on Apple Silicon, confirmed working July 2026 — the ">5 min" estimate from the ifort/PGO era no longer applies to the gfortran build). Run via `make test TAGS=slow` or `make test-full`. |
| `O2-CASSCF` | CASSCF | TFVC, SPIN, ENPART CASSCF, DM=2 | excluded: segfault (see Known Issues #12) |

**Validation notes:** All four fast-tier tests (plus `C2H6-B3LYP` in the slow tier) produce Normal Termination and
numerical values within floating-point rounding of the ifort/PGO reference
outputs. Differences are confined to the last 1–4 digits of 7-decimal
quantities — attributable to compiler and architecture differences
(x86-64 ifort PGO vs aarch64 gfortran).

### Tolerance tiers

Checks in `manifest.json` use `tol_abs` (absolute tolerance). Three tiers are
used consistently:

| Tier | `tol_abs` | Used for |
|------|-----------|---------|
| tight | `1e-4` or `1e-5` | Total energies, integration errors |
| normal | `1e-3` | Atomic charges, spin populations, bond orders, IQA energies |
| loose | `5e-3` or `0.05` | OSLO FOLI values, EOS electron counts, reliability indices |

### `manifest.json` — check format

Each check in the `"checks"` list of a test entry is a JSON object:

```json
{
  "label":       "Human-readable name shown in terminal output",
  "type":        "float",        // "float" (default), "present", or "absent"
  "section":     "REGEX",        // optional: search only after this marker in output
  "pattern":     "REGEX",        // regex with exactly one capture group (float type)
                                 // or just a search pattern (present/absent type)
  "match_index": 1,              // optional: which findall() match to use (1-based)
  "ref":         -76.2241116,    // reference value (float type only)
  "tol_abs":     1e-4            // absolute tolerance (float type only)
}
```

**Check types:**
- `"float"` (default): extract a number via `re.findall(pattern, text)[match_index-1]`, compare with `ref ± tol_abs`.
- `"present"`: check that `pattern` appears anywhere in the output (or section). Fail if not found.
- `"absent"`: check that `pattern` does NOT appear. Fail if found.

**`section` field:** When present, the runner finds the first occurrence of `section`
in the output text, then restricts the search to the text that follows it. This
is essential when the same pattern appears multiple times (e.g. atom tables in
multiple analysis sections).

### How to add a new test

1. **Prepare inputs.** Place `SystemName.fchk` and `SystemName.inp` in
   `compiler-testset/`. For OSLO calculations also add the `*-OSLOs.fchk`
   and/or `*-OSLOs-preortho.fchk` files.

2. **Generate a reference output.** Run the code once:
   ```bash
   cd compiler-testset
   ulimit -s unlimited
   ../apost3d SystemName > SystemName.apost 2>&1
   ```
   Verify it terminates normally, then copy the output:
   ```bash
   cp compiler-testset/SystemName.apost tests/reference/SystemName.apost
   ```

3. **Add an entry to `tests/manifest.json`.** Copy an existing entry and adapt:
   ```json
   {
     "name": "SystemName",
     "description": "One-line description of what this tests",
     "tags": ["dft", "enpart"],          // pick from existing tags or add new ones
     "timeout": 300,                      // seconds; be generous
     "extra_input_files": [],             // e.g. ["SystemName-OSLOs.fchk"]
     "checks": [
       { "label": "Normal Termination", "type": "present",
         "pattern": "Normal Termination" },
       { "label": "Total KS-DFT energy (au)",
         "pattern": "Total KS-DFT energy\\s*:\\s*([-+]?\\d+\\.\\d+)",
         "ref": -123.4567890, "tol_abs": 1e-4 },
       ...
     ]
   }
   ```
   Add the entry to the `"tests"` array. The order determines execution order.

4. **Verify the patterns match.** Run:
   ```bash
   python3 tests/run_tests.py --filter SystemName
   ```
   All checks should pass against the reference. If a pattern fails, open the
   `.apost` file and locate the line — then adjust the pattern or add a
   `section` anchor.

5. **Commit everything.** Add and commit:
   - `compiler-testset/SystemName.fchk`
   - `compiler-testset/SystemName.inp`
   - `tests/reference/SystemName.apost`
   - Updated `tests/manifest.json`

### Updating reference outputs

After an intentional code change that alters numerical output (bug fix, new
algorithm), regenerate all references:

```bash
make update-ref
```

This re-runs all tests, writes new `.apost` files to `tests/reference/`, and
updates the `"ref"` values in `manifest.json` automatically. Commit the updated
files as part of the same PR that contains the code change.

### Keyword coverage tracking

The test suite tracks which APOST-3D input keywords are exercised by the active
tests. This is stored in two places:

- **`tests/keywords.json`** — master registry of every keyword APOST-3D accepts,
  organised by category (`partitioning`, `population`, `orbital`, `energy`,
  `analysis`, `selection`, `interface`, `grid`, `output`). Each entry carries:
  - `id` — unique identifier used in manifest `keywords` lists (sub-keywords use
    `SECTION/keyword` notation, e.g. `ENPART/B3LYP`)
  - `key` — exact string as it appears in the `.inp` file
  - `section` — the `.inp` section it belongs to (`# METHOD`, `# ENPART`, …)
  - `description` — one-line explanation
  - `priority` — `high` / `medium` / `low` (whether this capability urgently
    needs a test)

- **`manifest.json` `"keywords"` field** — manually curated list of keyword ids
  that each test case exercises. This is the authoritative source; update it
  whenever you add or modify a test.

Generate the coverage report:

```bash
make coverage                          # full coloured report
make coverage UNCOVERED=1              # only untested keywords
make coverage PRIORITY=high UNCOVERED=1  # high-priority gaps only
make coverage CATEGORY=energy          # single category
make coverage FORMAT=json              # machine-readable JSON
python3 tests/coverage.py --help       # all options
```

**Current status (4 active tests):** 17 / 82 keywords covered (20%).
High-priority gaps include: `BECKE-RHO`, `HIRSH`, `HIRSH-IT`, `QTAIM`,
`MULLI`, `LOWDIN`, `EFFAO`, `EOS-U`, `LOBA`, `POLAR`, `EDAIQA`,
`ENPART/HF`, `ENPART/CASSCF`, `DM/PYSCF`, `DM/ORCA`, `MOKIT`, `ORCA`.

#### Adding a keyword to the registry

When a new keyword is added to the source code:

1. Add an entry to `tests/keywords.json`:
   ```json
   { "id": "NEWKEYWORD", "key": "NEWKEYWORD", "section": "# METHOD",
     "category": "analysis",
     "description": "One-line description",
     "priority": "high" }
   ```
2. If you add a test that exercises it, add its `id` to that test's `keywords`
   list in `manifest.json`.
3. Run `make coverage` to verify the registry is consistent.

---

## Known Technical Debt / Areas for Improvement

These are issues identified during code review that should be addressed in the MAJOR-UPDATE branch:

1. **COMMON blocks everywhere**: Global state is passed via named COMMON blocks redundantly declared in each source file. These should be migrated to F90 modules (the transition is already started with `modules.f90` but is incomplete).

2. **Mixed Fortran standards**: Source files mix old fixed-form Fortran 77 style (`.f`) with some modern F90 constructs. Conventions are inconsistent (implicit typing `IMPLICIT REAL*8(A-H,O-Z)` used throughout, `REAL*8` instead of `REAL(KIND=8)`, etc.).

3. **Compiler dependency**: ~~Hard-wired to Intel `ifort`.~~ **Resolved in Phase 0** — `Makefile` and helper scripts replace the ifort/PGO build with a portable gfortran build.

4. **Auto-parallelization**: ~~Relies on Intel's auto-parallelizer (`-parallel`).~~ **Partially resolved in Phase 0** — gfortran build uses `-fopenmp` with `OMP_NUM_THREADS=1`. Explicit `!$OMP PARALLEL DO` directives on the numerical integration loops are planned for Phase 1.

5. **`parameter.h` as include file**: Array size limits are set at compile time via `parameter.h`. Migrating to dynamic allocation throughout (and removing the fixed-size COMMON arrays) would remove all these limits.

6. **`devel.f`**: Contains deprecated routines and developer discussion comments (e.g., `!! MG: Pedro, aquesta subrutina ja la pots borrar !!`). Should be cleaned up — dead code removed, active code moved to appropriate files.

7. **`print.f`**: The element symbol array (`mend`) encoding H through U (92 elements) using 4-character Hollerith data — stops at uranium, no support for transuranic elements.

8. **Incomplete comments/documentation**: Many routines have minimal or no inline documentation. Variable names follow old conventions (`igr` = number of basis functions, `kop` = open-shell flag, etc.) that are not self-explanatory.

9. **Input parser**: `readchar`/`readint` are custom keyword parsers that operate by linear search through the `.inp` file on every call (re-reading from the beginning each time via `rewind`). Not a bottleneck but architecturally fragile.

10. **libxc pinned at 4.2.3**: Cannot upgrade to newer libxc versions due to interface changes. This limits access to newer functionals.

11. **`TO CHANGE` / `TO DO` markers**: Several `!! TO CHANGE !!` and `! TO DO:` annotations exist in the code marking incomplete implementations (e.g., kinetic energy density for meta-GGA, `xkdens` allocation in `main.f`).

12. **O2-CASSCF segfault**: The CASSCF test case (`O2-CASSCF`, TFVC + SPIN + ENPART + DM=2 with PySCF RDMs) crashes with a segfault in `main.f:898` immediately after printing "POST-HARTREE-FOCK CALCULATION / Number of core + active spin-orbitals: 26". No reference output exists for this case. The crash occurs with both `-O3` and `-O0 -fbounds-check`, and with unlimited stack, suggesting a pointer/allocation issue in the CASSCF ENPART path rather than a stack overflow. Needs investigation.

13. **`istart`/`iend` out-of-bounds bug (fixed in Phase 0)**: In `enpart.f`, the parallel grid-slicing code in `numint_two`, `numint_two_uhf`, and `calc_coul` contained `istart(ithreads)=iend(ithreads-1)+1` which accesses `iend(0)` when `ithreads=1`. ifort coincidentally read zero from heap memory; gfortran does not, causing integration errors for multi-atom systems. Fixed by replacing with `istart(ithreads)=ioffset+((ithreads-1)*ispace)+1` in all 8 occurrences.

14. **libxc archive format mismatch on macOS (fixed in Phase 0, July 2026)**: `compile_libxc.sh` let autotools pick whatever `ar` it found on PATH. On macOS this can resolve to a GNU-format `ar` (e.g. bundled with a Homebrew GCC toolchain, or from Homebrew binutils), which writes archives with a GNU-style symbol table (a member literally named `/`). Apple's system `ld` only understands the BSD archive format (`__.SYMDEF`) and fails at *link time* — not compile or archive time — with `ld: archive member '/' not a mach-o file in '.../libxc-4.2.3/lib/libxc.a'`. This is a toolchain mismatch, not a code bug; `libxcf90.a`/`libxcf03.a` (built via libtool) were unaffected, only the core `libxc.a` (archived directly). Fixed by forcing `AR=/usr/bin/ar RANLIB=/usr/bin/ranlib` (Apple's native, BSD-format, always present via Xcode CLT) into libxc's `./configure` on Darwin, plus a post-`make install` `ranlib` pass on every `.a` as a safety net. Linux is unaffected (GNU ar/ranlib is correct there and unchanged by this fix). Not independently verified end-to-end on macOS — needs confirmation after a real rebuild.

15. **dyld launch failure on macOS — codesigning hardened, but NOT the root cause here (Phase 0, July 2026)**: A gfortran-compiled binary that links successfully can still fail to *launch* on macOS with `dyld[...]: Library not loaded: /usr/lib/libSystem.B.dylib ... (no such file, no dyld cache)`, which looks like a missing system library. `make_compile.sh` was hardened to clear extended attributes (`xattr -c`), surface real `codesign` errors instead of swallowing them, and actually launch all three binaries as a build-time smoke test. **However**, on the actual failing machine (macOS 26.5.1, Apple Silicon), re-signing an already-built `apost3d` did *not* fix the launch failure — ruling out codesigning/quarantine as the cause in this case. See Known Issue #16 for the actual confirmed root cause. The codesign hardening is still worth keeping (it's a real, separate failure mode this project's binaries can hit), just wasn't what was happening here.

16. **dyld launch failure on macOS — actual root cause: ~2.6GB static BSS from `nmax=8000` (Phase 0, July 2026)**: `effao.f`, `mulliken.f`, and `qtaim.f` declare several fixed-size `nmax×nmax` `REAL*8` arrays in old-style `COMMON` blocks (`/effao/`, `/nao/`, `/stv/`) — with `nmax=8000` (`parameter.h`) these alone total ~2.6GB of static BSS compiled directly into the executable. Confirmed by bisection on the reporting machine (macOS 26.5.1, Apple Silicon, SIP enabled): a trivial `clang`/`gfortran` "hello world" and the small `eos_aom` utility (no large COMMON blocks) all launch fine; `apost3d` and `apost3d-eos` (which both pull in `effao.f`/`mulliken.f`) fail with the shared-region dyld error; lowering `nmax` from 8000 to 2000 made `apost3d` launch successfully. Root cause: a Position-Independent Executable with a multi-gigabyte fixed BSS segment can collide with the fixed address window macOS reserves for the shared dyld cache — arm64 macOS mandates PIE, so the classic `-no_pie` workaround isn't available. **Interim fix**: `nmax` lowered from 8000 to 3000 (Known Issue comment in `parameter.h`) — confirmed working, still ~9x more basis functions than any current test case uses. **Not a permanent fix**: a real system needing more than `nmax` basis functions will silently overrun these fixed arrays (undefined behavior, not a bounds-checked error) rather than failing loudly, and the underlying dyld collision could in principle recur at a large enough `nmax` on other machines/macOS versions. **Proper fix (planned)**: convert `/effao/`, `/nao/`, `/stv/` (and any other `nmax`-sized COMMON blocks) to `ALLOCATABLE` module arrays, following the pattern `modules.f90` already uses for `c`, `p`, `ps`, `pa`, `pb` — sized at runtime from the actual basis set, on the heap, not in the binary's static BSS. This removes the compile-time cap entirely (not just raises it) and eliminates the macOS launch failure at its source rather than working around it.

---

## Planned Work (MAJOR-UPDATE Branch)

Tasks to be worked on progressively in this branch:

- [x] **Compiler migration (Phase 0, May 2026)**: `Makefile`, `compile_libxc.sh`, `make_compile.sh` replace the ifort/PGO build. GCC 15 / gfortran tested on Linux (aarch64) and macOS 26 (arm64). Removed `use IFPORT` from `enpart.f`. Fixed latent `istart`/`iend(0)` out-of-bounds bug in two-electron ENPART routines (see Known Issues #13).
- [x] **Test infrastructure (Phase 0, May 2026)**: `tests/run_tests.py` + `tests/manifest.json` provide structured regression testing with per-quantity tolerance-based comparison, colored terminal output, HTML/text reports, `--filter`/`--tags`/`--update-ref` flags, and `make test` / `make update-ref` targets. 4 of 6 systems active (H2O-T-B3LYP, CH3F, FeCO2-PBEPBE, FeO4-2); C2H6-B3LYP excluded (runtime > 5 min), O2-CASSCF excluded (segfault, see Known Issues #12).
- [x] **Solid-base verification (July 2026)**: Full clean rebuild and full test run repeated independently with GCC 11 / gfortran on Linux (aarch64). All 4 active tests pass with every checked value matching the recorded reference exactly (Δ 0.0e+00), confirming numerical portability across compiler versions/vendors and OS. O2-CASSCF segfault and C2H6-B3LYP slow runtime both reproduced as documented.
- [ ] **O2-CASSCF segfault**: Investigate and fix the crash in the CASSCF ENPART path (main.f:898).
- [ ] **Explicit OpenMP parallelization (Phase 1)**: Replace single-threaded Phase 0 build with explicit `!$OMP PARALLEL DO` directives on the key numerical integration loops in `enpart.f`, `numint.f`, `wat.f`.
- [ ] **Module migration**: Convert remaining COMMON blocks to F90 module variables.
- [ ] **HIGH PRIORITY — convert `/effao/`, `/nao/`, `/stv/` to allocatable (July 2026)**: Not just cleanup — these `nmax×nmax` fixed-size COMMON arrays (`effao.f`, `mulliken.f`, `qtaim.f`) are the confirmed root cause of a macOS launch failure (Known Issue #16). `nmax` is currently lowered to 3000 (from 8000) as an interim workaround, which silently caps the max supported basis-set size. Converting these to `ALLOCATABLE`, sized at runtime like `modules.f90` already does for `c`/`p`/`ps`/`pa`/`pb`, removes the cap entirely and fixes the macOS issue at its source.
- [ ] **Code homogenization**: Consistent style, naming conventions, indentation, `IMPLICIT NONE` throughout.
- [ ] **Comment and documentation**: Add subroutine headers documenting purpose, arguments, and references.
- [ ] **Dead code removal**: Clean up `devel.f`, `.to_do` files, `wat.f.old`.
- [ ] **libxc upgrade**: Investigate and implement compatibility with libxc >= 5.x.
- [x] **Makefile module dependency fix (July 2026)**: `$(SRCDIR)/modules.o` had no prerequisites and per-file object rules didn't depend on it, so editing `modules.f90` (or rebuilding under a different gfortran version) silently failed to trigger recompilation of files that `use` those modules. Fixed: `modules.o` now depends on `modules.f90`, and every object rule (main sources, `input2.o`, utils) depends on `modules.o`. Verified: touching `modules.f90` now correctly triggers a full dependent rebuild instead of "Nothing to be done".
- [x] **Build preflight + compiler-drift check (July 2026)**: `make_compile.sh` now checks `gfortran` is present and >= 10 up front with a clear actionable error (was previously a raw `Error 127` mid-build-log). It also stamps the compiler identity used for each build in `objects/.gfortran_version` and automatically forces `make clean` if the compiler changes since the last build (catches the exact "module file created by a different version of GNU Fortran" failure mode — file-mtime-based Make dependencies can't detect a compiler swap on their own). Verified with a simulated version change.
- [x] **CI (GitHub Actions) (July 2026)**: `.github/workflows/test.yml` builds libxc + APOST-3D and runs `make test-only` on a fresh Ubuntu runner on every push/PR to `master`/`MAJOR-UPDATE`, plus manual dispatch. Caches the libxc build (keyed to invalidate on compiler/tarball change). Uploads `tests/report/` as a build artifact and runs `make coverage`. By design the workflow never needs editing to add/remove tests — it just runs whatever `tests/manifest.json` defines; adding a test is purely a data change (input files + reference output + manifest entry). Not yet exercised on real GitHub infrastructure (needs a push to verify).
- [x] **"slow" test tier (July 2026)**: `C2H6-B3LYP` confirmed working (~85s single-threaded on the reporting machine, Normal Termination, 10/10 checks pass) — added to `tests/manifest.json` tagged `slow`, with `tests/reference/C2H6-B3LYP.apost`. `run_tests.py` gained `--exclude-tags`; `make test`/`test-only` now pass `--exclude-tags slow` by default (only when `TAGS` isn't explicitly given), and a new `make test-full` target runs everything unfiltered. `make test TAGS=slow` runs just the slow tier. Keyword coverage unchanged (17/82) — `C2H6-B3LYP` exercises already-covered keywords at larger scale, doesn't add new ones. `O2-CASSCF` remains excluded entirely (segfault, Known Issue #12), not just tagged slow.
- [ ] **Hosted documentation (Read the Docs)**: Once `DOCUMENTATION.md`/`CLAUDE.md` grow further, publish them via Sphinx + MyST-parser (consumes existing Markdown near-verbatim) on readthedocs.org — free for open-source, auto-rebuilds on push via GitHub webhook, gives versioned docs (useful once `master` and `MAJOR-UPDATE` diverge further). Low incremental cost given docs are already Markdown-first.

---

## Git Workflow

- **Main branch**: `master` — stable, matches public GitHub release
- **Development branch**: `MAJOR-UPDATE` — all major refactoring work (current branch)
- **Past branches**: `IQA-CASSCF`, `psalse` (merged or in progress)
- **Convention**: Work incrementally on `MAJOR-UPDATE`, with frequent commits. Run `make test` after every significant change. When output changes intentionally, run `make update-ref` and commit the updated references together with the code change.

---

## Authors

- **Pedro Salvador** (psalse@gmail.com) — lead developer, Universitat de Girona
- **Eloy Ramos-Cordoba** (eloy.raco@gmail.com) — EOS, local spin, EFFAO
- **Marc Montilla** — polarizability decomposition
- **Lluís Pujal** — contributions
- **Martí Gimferrer** (mgimferrer18@gmail.com) — OSLO, LOBA, EOS extensions, EDA-IQA, MOKIT/ORCA interfaces, ongoing development

Original APOST code by I. Mayer and A. Hamza (Budapest, 2000–2003).
