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
  effective fragment orbitals (EFOs) from EOS and GEOS, can be written out
  to a standard `.fchk` file, so they can be viewed directly in any
  common orbital-viewer instead of only as cube files.
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
