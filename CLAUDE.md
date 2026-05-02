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
├── libxc-4.2.3/        # Bundled libxc library (version 4.2.3, pinned)
├── libxc-4.2.3.tar.gz  # Tarball for re-extraction
├── compiler-testset/   # Input files (.fchk, .inp, .dm1, .dm2) for test runs
├── compiler-runtest    # Script: runs test suite (step 1, generates .dyn files)
├── compiler-runtest2   # Script: runs test suite (step 2, validates PGO build)
├── media/              # Logo and images for GitHub/documentation
├── Makefile_profgen    # Makefile for step 1 of PGO compilation
├── Makefile_profuse    # Makefile for step 2 of PGO compilation
├── make_compile.sh     # Master compilation script (runs both steps + utils)
├── compile_libxc.sh    # Script to build the bundled libxc library
├── README.md           # GitHub README
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
- `wat.f.old` — old version of weight function routines

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
| `/iops/` | `iopt(100)` — 100-element integer option array |
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
export OMP_NUM_THREADS=48
export KMP_STACKSIZE=100m
ulimit -s unlimited
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

### Prerequisites
- Intel oneAPI toolkits (2024 recommended; 2025 drops `ifort` by default)
  - Intel oneAPI Base Toolkit
  - Intel oneAPI HPC Toolkit
- Environment: `source /opt/intel/oneapi/setvars.sh intel64`
- Environment variable `APOST3D_PATH` must be set

### Build targets (both Makefiles)
- `apost3d` — main executable
- `apost3d-eos` — standalone EOS executable (subset of sources)
- `eos_aom` — utility executable (`utils/eos_aom.f90`)

### Two-step PGO compilation (`make_compile.sh`)
1. **Step 1** (`Makefile_profgen`): compile with `-prof-gen` flag; run test suite to generate `.dyn` profiling files
2. **Step 2** (`Makefile_profuse`): recompile with `-prof-use -Ofast` using the profiling data; run test suite again to validate

### Key compiler flags
- `-parallel -qopenmp` — Intel auto-parallelization + OpenMP (both enabled)
- `-unroll-aggressive` — aggressive loop unrolling
- `-xHost` — optimize for the current CPU architecture
- `-extend-source 132` — allow fixed-form lines up to 132 characters
- `-mcmodel=medium` — for large binary/data
- `FOR_INP = ifort` (without `-Ofast`) — `input2.f` compiled without aggressive optimization to avoid parsing issues

### Parallelization
- Runtime threads: `OMP_NUM_THREADS` (recommended: set to max cores on node)
- Stack size: `KMP_STACKSIZE=100m`, `ulimit -s unlimited`
- Intel auto-parallelization + OpenMP are both used. The `-parallel` flag lets ifort auto-parallelize loops.

### Environment variable
```bash
export APOST3D_PATH="/path/to/APOST3D"
```

### Utils compilation
```bash
make -f Makefile_profuse util
```

---

## Utilities (`utils/`)

| File | Description |
|------|-------------|
| `apost3d.py` | Python wrapper script |
| `main_eos.f` | Main program for the standalone `apost3d-eos` executable |
| `eos_aom.f90` | Standalone EOS-AOM (atom/orbital mapping) utility |
| `eos_alt.f90` | Alternative EOS calculation |
| `eos_aom.f90` | EOS from atomic overlap matrix |
| `gen_hirsh.f` | Generate Hirshfeld promolecular densities |
| `get_energy.f` | Extract energy values from output files |
| `get_energy_g16.f` | Same, for Gaussian 16 outputs |
| `group_frag.f` | Fragment grouping utility |
| `wfn2fchk.f90` | Convert `.wfn` files to `.fchk` format |
| `memchk.f` | Memory checking utility |
| `insert.sh` | Shell script helper |

---

## Test Suite (`compiler-testset/`)

Used for PGO profiling and regression testing. Five test systems:

| System | Wavefunction | Features tested |
|--------|-------------|-----------------|
| `H2O-T-B3LYP` | RKS B3LYP | TFVC, ENPART with DFT, SPIN, MOD-GRIDTWOEL |
| `C2H6-B3LYP` | RKS B3LYP | TFVC, ENPART, fragments, THREBOD |
| `CH3F` | (DFT) | TFVC, ENPART |
| `FeCO2-PBEPBE` | UKS PBE | TFVC, DOFRAGS, EOS (open-shell, fragments) |
| `FeO4-2` | UKS | TFVC, DOFRAGS, OSLO, QCHEM, fragment-based OSLO |
| `O2-CASSCF` | CASSCF | TFVC, SPIN, ENPART CASSCF, DM=2 (external RDMs from PySCF) |

---

## Known Technical Debt / Areas for Improvement

These are issues identified during code review that should be addressed in the MAJOR-UPDATE branch:

1. **COMMON blocks everywhere**: Global state is passed via named COMMON blocks redundantly declared in each source file. These should be migrated to F90 modules (the transition is already started with `modules.f90` but is incomplete).

2. **Mixed Fortran standards**: Source files mix old fixed-form Fortran 77 style (`.f`) with some modern F90 constructs. Conventions are inconsistent (implicit typing `IMPLICIT REAL*8(A-H,O-Z)` used throughout, `REAL*8` instead of `REAL(KIND=8)`, etc.).

3. **Compiler dependency**: Currently hard-wired to Intel `ifort` with Intel-specific flags (`-prof-gen`, `-prof-use`, `-parallel`, `-qopenmp`, `-xHost`). Intel oneAPI 2025 drops `ifort`. Migration to `ifx` (Intel's new compiler) or a standard compiler (`gfortran` + OpenMP) is required.

4. **Auto-parallelization**: The current parallelization relies entirely on Intel's auto-parallelizer (`-parallel`). This should be replaced with explicit OpenMP directives for portability and predictability.

5. **`parameter.h` as include file**: Array size limits are set at compile time via `parameter.h`. Migrating to dynamic allocation throughout (and removing the fixed-size COMMON arrays) would remove all these limits.

6. **`devel.f`**: Contains deprecated routines and developer discussion comments (e.g., `!! MG: Pedro, aquesta subrutina ja la pots borrar !!`). Should be cleaned up — dead code removed, active code moved to appropriate files.

7. **`print.f`**: The element symbol array (`mend`) encoding H through U (92 elements) using 4-character Hollerith data — stops at uranium, no support for transuranic elements.

8. **Incomplete comments/documentation**: Many routines have minimal or no inline documentation. Variable names follow old conventions (`igr` = number of basis functions, `kop` = open-shell flag, etc.) that are not self-explanatory.

9. **Input parser**: `readchar`/`readint` are custom keyword parsers that operate by linear search through the `.inp` file on every call (re-reading from the beginning each time via `rewind`). Not a bottleneck but architecturally fragile.

10. **libxc pinned at 4.2.3**: Cannot upgrade to newer libxc versions due to interface changes. This limits access to newer functionals.

11. **`TO CHANGE` / `TO DO` markers**: Several `!! TO CHANGE !!` and `! TO DO:` annotations exist in the code marking incomplete implementations (e.g., kinetic energy density for meta-GGA, `xkdens` allocation in `main.f`).

---

## Planned Work (MAJOR-UPDATE Branch)

Tasks to be worked on progressively in this branch:

- [ ] **Compiler migration**: Replace `ifort` with `ifx` or `gfortran`, remove PGO and Intel-specific flags, create a portable Makefile
- [ ] **Explicit OpenMP parallelization**: Replace Intel auto-parallelizer with explicit `!$OMP PARALLEL DO` directives on the key numerical integration loops
- [ ] **Module migration**: Convert remaining COMMON blocks to F90 module variables
- [ ] **Code homogenization**: Consistent style, naming conventions, indentation, implicit-none throughout
- [ ] **Comment and documentation**: Add subroutine headers documenting purpose, arguments, and references
- [ ] **Test infrastructure**: Convert `compiler-testset` into a proper regression test suite with reference outputs and automated comparison
- [ ] **Dead code removal**: Clean up `devel.f`, `.to_do` files, `wat.f.old`
- [ ] **libxc upgrade**: Investigate and implement compatibility with libxc >= 5.x

---

## Git Workflow

- **Main branch**: `master` — stable, matches public GitHub release
- **Development branch**: `MAJOR-UPDATE` — all major refactoring work (current branch)
- **Past branches**: `IQA-CASSCF`, `psalse` (merged or in progress)
- **Convention**: Work incrementally on `MAJOR-UPDATE`, with frequent commits. Test after each significant change using the `compiler-testset` inputs.

---

## Authors

- **Pedro Salvador** (psalse@gmail.com) — lead developer, Universitat de Girona
- **Eloy Ramos-Cordoba** (eloy.raco@gmail.com) — EOS, local spin, EFFAO
- **Marc Montilla** — polarizability decomposition
- **Lluís Pujal** — contributions
- **Martí Gimferrer** (mgimferrer18@gmail.com) — OSLO, LOBA, EOS extensions, EDA-IQA, MOKIT/ORCA interfaces, ongoing development

Original APOST code by I. Mayer and A. Hamza (Budapest, 2000–2003).
