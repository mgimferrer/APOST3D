# Changelog

What changed for users in each release of APOST-3D. Each entry links to
the documentation page with the details. The full history is in
`git log`.

## [5.0] — Unreleased

Version 5.0 brings two oxidation-state analyses for open-shell systems
(GEOS, and OSLO for unrestricted wavefunctions), wavefunctions from ORCA,
meta-GGA functionals in ENPART and a machine-readable results file. The
program builds with free compilers and libraries only, comes with a test
suite and a [documentation site](https://apost3d.readthedocs.io), and
stops with a clear message where version 4 crashed or gave wrong numbers.
Some version 4 inputs need small changes, listed first.

### Changes you must make to your inputs

- **One atomic definition** in `# METHOD`. Without one the run uses
  `TFVC` and warns (version 4 used Becke atoms with fixed radii, now the
  `BECKE` keyword), and more than one stops the run. See
  [# METHOD](docs/source/input/method.md#atomic-definition).
- **EOS and GEOS need fragments**: `DOFRAGS` and a `# FRAGMENTS` block.
  For one oxidation state per atom, make every atom a fragment. See
  [Fragments](docs/source/guide/fragments.md).
- **Grids are set in `# GRID`**: `RADIAL`/`ANGULAR` for every analysis,
  `RADIAL_2E`/`ANGULAR_2E` for ENPART. The grid on the command line
  (`apost3d job 70 434`), `MOD-GRIDTWOEL`, `FINEGRID` (now `ANGULAR 974`),
  `phb1`/`phb2`, `rr00` (now `R0_2E`) and `THRESH2` stop the run with a
  message. See [# GRID](docs/source/input/grid.md).
- **Keywords are whole words**, and three `# OSLO` keywords are written
  with an underscore: `FOLI_TOLERANCE`, `BRANCH_ITERATION` and
  `PRINT_NONORTHO` (the old spellings stop). See
  [# OSLO](docs/source/input/oslo.md).
- **ENPART takes one method**: more than one of `HF`, a functional
  keyword, `LIBRARY`, `CASSCF` and `CISD` stops. `LDA` is gone (it meant
  Slater exchange only): use `SVWN`, `SVWN5`, or `LIBRARY` with
  `EX_FUNCTIONAL 1`. Range-separated hybrids, VV10 functionals and
  Hilbert-space atoms (`MULLIKEN`, `LOWDIN`) stop before any energy is
  computed. See [# ENPART](docs/source/input/enpart.md).
- **ENPART results change with two new defaults.** In KS-DFT, atom pairs
  with a bond order below `THREBOD` (now 0.005, `THREBOD 50`, was 100)
  take their whole exchange-correlation term from the multipolar
  expansion. In the systems checked the results stay within 0.04 kcal/mol
  of computing every pair, but the interaction and self terms of
  fragments joined by many weak pairs can change by equal and opposite
  amounts (about 11 kcal/mol between the two rings of ferrocene). `THREBOD -1` computes every pair. The
  zero-error correction is applied only when the two-electron error
  exceeds 0.25 kcal/mol (`TWOELTOLER 0.25`, `TWOELTOLER 0.00` applies it
  always, as version 4 did). See [ENPART](docs/source/methods/enpart.md).
- **Analyses that version 4 read but did not compute stop**: `EDAIQA`,
  `TOPOLOGY`, `DAFH`, `SHANNON`, `WFN` and `READINT` in `# METHOD`, and
  `DMRG` in `# DM`.
- **Results of fixed bugs**: the entries marked "Results change" under
  Fixed give different numbers than version 4.

### New

- **Results file** `<jobname>.apost.json`: every run writes its results
  (wavefunction, molecule, fragments, populations, bond orders, local
  spin, EFFAO, EOS, GEOS, OSLO and ENPART) in a machine-readable file,
  for the APOST3Dview visualizer and for scripts. See
  [Results file (JSON)](docs/source/developer/json-output.md).
- **GEOS**: oxidation states of open-shell systems from the paired and
  unpaired densities. Paired orbitals with a negative occupation are
  reported apart. GEOS is not published yet: see the warning on the
  [GEOS](docs/source/methods/geos.md) page.
- **OSLO for unrestricted wavefunctions** (UHF, UKS), and branching:
  `BRANCH_ITERATION` (`BRANCH_ALPHA`, `BRANCH_BETA`) follows an
  alternative assignment, and a warning names close alternatives. See
  [OSLO](docs/source/methods/oslo.md).
- **EFO orbital files**: EOS and GEOS write their effective fragment
  orbitals to `<jobname>-EOS-EFOs.fchk` (`-GEOS-`) for any orbital
  viewer, and with `BECKE`, `BECKE-RHO` or `TFVC` also the orbitals
  before the multiplication by the fragment's weight
  (`-EFOs-unweighted.fchk`). See
  [Visualizing orbitals](docs/source/guide/visualization.md).
- **ORCA wavefunctions**: `utils/orca2fchk.py` converts ORCA 6
  calculations (RHF, UHF, ROHF, RKS, UKS, ROKS and single-root CASSCF,
  with pseudopotentials), checked against pySCF (and Gaussian 16 for
  ENPART) on transition-metal complexes. See
  [Preparing the wavefunction](docs/source/guide/wavefunctions.md#orca).
- **ENPART needs a calculation without density fitting**, in every
  program (ORCA: `NoRI`). `orca2fchk.py` and `apost3d.py` print a note
  when it was used. See
  [Density fitting and ENPART](docs/source/guide/wavefunctions.md#density-fitting-and-enpart).
- **ENPART with meta-GGA functionals** and their global hybrids,
  restricted and unrestricted, and functional keywords that reproduce the
  Gaussian 16 functionals of the same name (`SVWN`, `SVWN5`, `BLYP`,
  `BP86`, `PBE`, `B3LYP`, `B3PW91`, `B3P86`, `PBE0`, `BHANDHLYP`,
  `M062X`). See
  [# ENPART](docs/source/input/enpart.md#predefined-functionals).
- **ENPART atomic and classical interaction energies**: `ATOMIC
  ENERGIES` (self and additive energy of each atom) and `CLASSICAL
  INTERACTION ENERGIES`, also per fragment. See
  [ENPART](docs/source/methods/enpart.md#reading-the-output).
- **Cube options**: `SPACING` and `PADDING` in `# CUBE` set the grid and
  the box, `NEG_EFOS` writes the negative paired orbitals of GEOS, and by
  default every orbital gets a cube file. See
  [# CUBE](docs/source/input/cube.md).
- **`VERBOSE 1`** in `# METHOD` also prints the intermediate matrices
  (for now of ENPART). See [# METHOD](docs/source/input/method.md#options).
- **Exit code 1** when a run stops on an error, so job scripts can
  detect it. See [Running a calculation](docs/source/guide/running.md).
- **Test suite**: `make test` checks a build against reference results,
  and `make test-strict` also compares the whole output. See
  [Installation](docs/source/installation.md#3-check-the-build).
- **Documentation site** at
  [apost3d.readthedocs.io](https://apost3d.readthedocs.io).

### Changed

- **Build with free tools only**: gfortran and OpenBLAS replace the Intel
  compiler and MKL, and libxc 7.1.2 (was 4.2.3) is downloaded and built
  on the first compilation when no suitable one is installed.
  `make_compile.sh` checks every prerequisite first. The build is
  portable by default (`ARCH=native` optimizes for the build machine
  only) and uses 8 threads. See [Installation](docs/source/installation.md).
- **Faster**: diagonalization through LAPACK, and OpenMP in the
  numerical integration, ENPART, EFFAO/EOS, OSLO and the cube files.
- **Size limits**: 8000 basis functions and 500 atoms (were 3000 and
  350), and a larger file stops with a message. See
  [Preparing the wavefunction](docs/source/guide/wavefunctions.md#what-the-file-must-contain).
- **Output**: one aligned layout throughout, sections that start with
  the atomic definition, and value columns named after it (`TFVC`,
  `Lowdin`, ...). The orbital energies of `<jobname>-OSLOs.fchk` are the
  FOLI values. No number changes. See
  [Output files](docs/source/guide/output.md).
- **ENPART output**: every term printed once with its final value, the
  two-electron integration error in its own section, the exact exchange
  already scaled (`EXACT-EXCHANGE SHARE`), and the final matrix named
  `IQA ENERGY TERMS (TFVC)` (was `FUZZY ATOMS ... ENERGY COMPONENTS`). No
  number changes. See
  [ENPART](docs/source/methods/enpart.md#reading-the-output).
- **Input checks**: a misspelled job name, a missing file or block, a
  wrong fragment list, a value that can't be read, a grid that is not a
  Lebedev grid or a `.nao` file of another calculation stop the run with
  a message naming the problem, instead of a Fortran error, a silent
  default or wrong results. See
  [Troubleshooting](docs/source/troubleshooting.md).
- **`utils/apost3d.py`** (pySCF) rewritten: pySCF 2.7 and newer, g
  functions, CASCI, the RDMs of FCI and closed-shell CCSD and CISD, and the
  functional's name in the `.fchk`. Existing scripts work unchanged. See
  [Preparing the wavefunction](docs/source/guide/wavefunctions.md#pyscf).

### Fixed

- **OSLO**: an orbital found from several equivalent fragments is shared
  among them (ethane as two CH₃: 0/0, was −1/−1), incomplete fragments
  stop the run, Mulliken fragment populations no longer give NaN FOLI
  values, and few electrons per channel no longer write past the arrays.
  Results change.
- **OSLO runs that failed**: with Gaussian 09 or pySCF files (writing
  `<jobname>-OSLOs.fchk`), together with ENPART, and without beta
  electrons. A false close-alternative warning is gone.
- **EFFAO and EOS with `LOWDIN-DAVIDSON`** used plain Löwdin orbitals.
  Results change.
- **EOS and GEOS**: the frontier orbitals behind R(%) follow the
  electron assignment (the last digits change in pseudo-degenerate
  cases), the net occupation of a fragment sums its listed orbitals only,
  an `EFF_THRESH` that leaves too few orbitals stops, `EFFAO-U` alone no
  longer crashes, and the negative paired EFOs reach the `.fchk` file
  also with a linearly dependent basis.
- **Orbital `.fchk` and cube files**: values below 1e-99 lost their
  exponent letter (`1.27052287-118`) and are now written as 0. The EFO
  and OSLO files of ROHF/ROKS wavefunctions get Beta blocks and keep the
  density blocks of the input.
- **ENPART with pseudopotentials** no longer stops (`ECP Matrix not
  found`), and the reference electron-nuclear energy of pySCF files
  includes the ECP energy (it distorted the zero-error correction).
  Results change.
- **POLAR field energies** counted the charge-transfer term twice: the
  total is now **μ**·**F**. Results change.
- **Kohn-Sham `.fchk` files** whose functional name contains HF, MP or CI
  (M06HF, HFS, mPW1PW91, the KCIS functionals) were reported as
  Hartree-Fock, MP or CI, and the KCIS ones were analysed as correlated
  wavefunctions. Results change for them.
- **`DENS n`** counts the densities the same way in Gaussian, Q-Chem and
  pySCF files (with Q-Chem, `DENS 2` read the third density). Results
  change.
- **Bond orders of correlated wavefunctions**: the exchange-only bond
  orders of the 1-RDM include the spin-density term, and Mayer's
  improved bond orders are printed with the populations. A correlated
  non-singlet without `# DM` gets no bond orders. Results change for non-singlets.
- **Iterative Hirshfeld** converges to 10⁻⁶ (was 5×10⁻⁴, which left
  charges up to 7×10⁻⁴ off). Results change.
- **Cube files**: no orbital-index line (strict readers shifted the
  data), the printed norm is right (it was 5–10% low), and the file
  names keep the whole job name.
- **Element labels and `.int` files**: Sm (was printed as Sn), every
  element up to Z = 118, and the `NET CHARGE` of atoms with an ECP.
- **`get_energy_g16`** stops on an incomplete Gaussian output instead of
  appending wrong energies, and has no limit of 3000 basis functions.

### Removed

- **Davidson spin decomposition matrix** of `SPIN`: its spin-density
  term had the wrong sign.
- **POLAR, and ENPART in an electric field, for correlated
  wavefunctions**: untested code, they now stop.
- **Programs** `apost3d-eos` (`apost3d` does everything it did),
  `group_frag` (ENPART with `DOFRAGS` prints its fragment terms) and
  `eos_alt`.
- **`efo_occ.dat` and `efo_coeff.dat`**: the EFOs of every scheme are in
  `<jobname>-EOS-EFOs.fchk`.
- **`get_energy`** (Gaussian 09) is deprecated: still built, no longer
  documented.
