# Changelog

This file summarizes what changed in APOST-3D **Version 5.0** relative to
the previous public release (`master`). It is meant for users of the
program, not a full development history — see `git log` for that.

## [5.0] — Unreleased

### New implementations

- **GEOS (Generalized Effective Oxidation States)** — a new open-shell
  oxidation-state analysis (`# METHOD GEOS`), built from paired/unpaired
  densities separately (Takatsuka's definition) rather than the total
  density alone, complementing the existing EOS method.
- **DFT-DM1** (in development, not documented yet) — a new analysis method (`# METHOD DFT-DM1`) for
  UHF/UKS-DFT wavefunctions: builds an approximate one-particle density
  matrix from the local exchange-energy density and prints per-atom-pair
  exchange energies and bond orders, with an optional natural-orbital
  analysis (`NATORB`).
- **OSLO extended to open-shell systems** — OSLO orbital analysis now
  fully supports unrestricted (UHF/UKS) wavefunctions, not just
  restricted ones.
- **Orbitals exportable to `.fchk`** — OSLO orbitals, and now also the
  effective fragment orbitals (EFOs) from EOS (with any population
  scheme) and GEOS, are written out to a standard `.fchk` file, so they
  can be viewed directly in any common orbital-viewer instead of only as
  cube files.
- **New cube-file options** — `SPACING`/`RADIUS_SCALE` keywords under
  `# CUBE` for finer control over the generated grid, backed by one
  unified, faster cube-generation routine.
- **GEOS negative paired EFOs** — paired EFOs with a clearly negative net
  occupation are reported, left out of the electron assignment and
  exported to the `.fchk` file; `# CUBE NEG_EFOS` writes their cube files.
  The exported EFOs approximate the fragment-weighted orbital shown in the
  cube files, with a per-EFO `FIT %` telling how closely.
- **Predefined functional keywords for ENPART** — `SVWN`, `SVWN5`, `BLYP`,
  `BP86`, `PBE` (or `PBEPBE`), `B3LYP`, `B3PW91`, `B3P86`, `PBE0` (or
  `PBE1PBE`) and `BHANDHLYP`, each checked to reproduce the Gaussian 16
  functional of the same name. Any other supported functional can be
  given by its libxc ids (`LIBRARY` with `EXC_FUNCTIONAL`, or
  `EX_FUNCTIONAL`/`EC_FUNCTIONAL`).

### Changes to ENPART that affect inputs or results

- **Weakly bonded atom pairs in KS-DFT** — for pairs with a bond order
  below the `THREBOD` threshold, the whole exchange-correlation term now
  comes from the multipolar expansion, as it already did for HF and
  CASSCF, instead of from BODEN. The default threshold is lowered from
  0.01 to 0.005 (`THREBOD 50`, was 100); in every system checked, the
  results stay within 0.04 kcal/mol of computing all pairs. Use
  `THREBOD -1` to compute every pair.
- **Zero-error strategy only when needed** — the new default
  `TWOELTOLER 0.25` applies it only when the two-electron integration
  error exceeds 0.25 kcal/mol (it used to be applied always);
  `TWOELTOLER 0.00` restores the old behaviour.
- **`LDA` keyword removed** — it meant Slater exchange only. Use `SVWN` or
  `SVWN5`, or `LIBRARY` with `EX_FUNCTIONAL 1` for exchange only.
- **Clear stops instead of wrong numbers** — meta-GGA functionals (not
  yet supported), range-separated hybrids, VV10 functionals, invalid libxc
  ids, and ENPART with a Hilbert-space atom definition (`MULLIKEN`, `LOWDIN`)
  now stop with a message before any energy is computed.
- **Clearer output** — for hybrid functionals, the DFT-part tables now
  say that the HF-type exchange share is added in the two-electron part,
  and point to the complete tables. The zero-error strategy says when it
  is applied, and warns when it extrapolates instead of interpolating.
- `utils/get_energy` and `utils/get_energy_g16`, which append the
  reference energies to a Gaussian `.fchk`, now handle up to 3000 basis
  functions.

### Other changes that affect inputs

- **EOS and GEOS need fragments** — like OSLO, they now stop without
  `DOFRAGS` and a `# FRAGMENTS` block instead of giving every atom its own
  oxidation state by default. For one oxidation state per atom, define
  every atom as a fragment of its own.
- **Atomic definition** — always write one in `# METHOD`. An input
  without any now runs with `TFVC` and a warning (it used plain Becke
  atoms with fixed radii); those are selected with the new `BECKE`
  keyword.

### Code optimization / parallelization

- **Fully open-source build** — the program now builds with gfortran and
  OpenBLAS instead of the Intel `ifort`/MKL/PGO toolchain, removing any
  dependency on licensed compilers.
- **Faster diagonalization** — matrix diagonalization now goes through
  LAPACK/OpenBLAS instead of the previous in-house routine.
- **Broader multi-core parallelization** — OpenMP parallelization was
  extended across the major computational bottlenecks (numerical
  integration, energy partitioning/ENPART, EFFAO/EOS, OSLO, DFT-DM1,
  cube-file generation), giving real speedups on multi-core machines
  beyond what was already parallelized before.

### Bug fixes

- **OSLO with Gaussian 09 or pySCF wavefunctions** — OSLO no longer stops
  with an end-of-file error (before printing the oxidation states) when
  writing `<jobname>-OSLOs.fchk`; the written file now keeps every block
  of the input `.fchk`.
- **OSLO together with ENPART** — an input asking for both no longer
  crashes after the energy decomposition.
- **Checked fragment definitions** — `# FRAGMENTS` with an atom that does
  not exist, an atom in two fragments, or `-1` before the last fragment
  now stops with a message naming the problem; atoms left out of the
  fragments stop EOS, GEOS, OSLO and ENPART (other analyses warn). OSLO
  and GEOS used to run on incomplete fragments and print wrong oxidation
  states.
- **OSLO input checks** — OSLO with a Mulliken, Löwdin or NAO scheme in
  `# METHOD` stops with a message (those schemes go in `# OSLO`), and
  open-shell systems with no beta electrons now run.
- **OSLO linear-dependency check** — with more candidate orbitals than
  electrons in a spin channel (few electrons per channel), it wrote past
  its arrays: the run could crash, or print a false linear-dependency
  warning or garbage oxidation states.
- **EOS with a changed orbital cutoff** — an `EFF_THRESH` that leaves
  fewer orbitals than electrons now stops instead of assigning electrons
  wrongly, very low or negative values no longer crash, and the overall
  reliability index R(%) is defined for systems without beta electrons.
- **EFFAO/EOS/GEOS fragment sums** — `Net occupation for fragment` is now
  the sum over the listed orbitals (above the cutoff) in every scheme, like
  the gross one (in real space it summed all orbitals, negative GEOS ones
  included), and a new `Left out by the cutoff` line gives what the listed
  orbitals miss of the fragment's population. It replaces the `Deviation
  from net/gross population` lines, which were not computed correctly.
  With Mulliken or Löwdin atoms, a fragment with a single basis function
  no longer loses its only orbital.
- **EFFAO/EOS with `LOWDIN-DAVIDSON`** — the orbitals are now built in the
  Davidson-Löwdin basis, as its populations are; they silently used plain
  Löwdin, so oxidation states and R(%) were Löwdin ones. `EFFAO`/`EOS`
  with the experimental `LOWDIN-W` now stop.
- **`efo_occ.dat`/`efo_coeff.dat` removed** — the EFOs of Löwdin,
  Davidson-Löwdin and NAO EOS now go to `<jobname>-EOS-EFOs.fchk` like any
  other scheme (Mulliken included); plain `EFFAO` with `LOWDIN` no longer
  leaves stray `fort.44`/`fort.45` files.
- **Input reading** — a misspelled job name, a missing `.fchk` or `# DM`
  file, a job name given with its extension, or a missing `# METHOD`,
  `# ENPART` or `# OSLO` block now stop with a message naming the file or
  block (a typo used to create empty files and end in a Fortran error).
  `KEY = value` with spaces is read, a value that can't be read stops the
  run showing the line (real keywords used to fall back silently to their
  default), and the last block may end without its closing `#`.
- **Hirshfeld atomic densities** — `HIRSH`/`HIRSH-IT` without a
  `densoutput` file, or with an incomplete one, stop with a message (a
  missing file used to be created empty and end in a Fortran error);
  samarium can now be matched in it. `gen_hirsh` stops on a spherical or
  g-function basis, or one larger than it can hold, instead of writing
  wrong densities.
- **`DM 1`/`DM 2` with a single-determinant `.fchk`** — now stops with a
  message instead of failing in the diagonalization of an empty matrix.
- **`DENS n`** — densities are now counted the same way in Gaussian,
  Q-Chem and pySCF files (with Q-Chem files `DENS 2` read the third
  density), the output lists them, and an unrestricted or ROHF
  wavefunction stops instead of silently using the SCF density.
- **Integration grids are checked** — an angular number of points that is
  not a Lebedev grid (in `# GRID` or on the command line), or a radial one
  outside 1–500, now stops the run. In `# GRID` it used to give wrong
  ENPART numbers, hidden by the zero-error interpolation; on the command
  line it was silently rounded down.

### Others

- **New compilation and testing setup** — the build is now driven by two
  scripts (`compile_libxc.sh`, `make_compile.sh`) that check for a
  compatible compiler/OpenBLAS/libxc automatically; an automated regression
  test suite now runs on every change, checked against independently
  computed reference values, and runs automatically (via GitHub Actions)
  on every update to the code. `make test` checks selected values and
  every printed number of each test; `make test-strict` also requires the
  whole output to match, to validate a build on another machine.
- **Portable build by default** — the code is compiled for the generic
  architecture, so one binary runs on every node of a cluster with CPUs
  of different ages; `ARCH=native` optimizes for the build machine only.
  Compilation and tests use 8 threads by default (`NTHREADS=<n>` to
  change it), and `make_compile.sh` (or `make utils`) also builds the
  utilities in `utils/` (`get_energy`, `get_energy_g16`, `gen_hirsh`,
  `wfn2fchk`, `group_frag`, `eos_aom`, `eos_alt`).
- **libxc upgraded to 7.1.2** (from the previously bundled 4.2.3) — no
  change to any computed energy or property, verified bit-for-bit against
  reference outputs. Fetched and built automatically on first compile if
  no suitable existing install is found, rather than shipped as a bundled
  copy.
- **New documentation website** — a full documentation site is now
  available at [apost3d.readthedocs.io](https://apost3d.readthedocs.io).
- **More consistent program output** — formatting of the printed output
  was made more uniform and readable throughout.
