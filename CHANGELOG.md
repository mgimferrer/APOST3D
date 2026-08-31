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
- **DFT-DM1** — a new analysis method (`# METHOD DFT-DM1`) for
  UHF/UKS-DFT wavefunctions: builds an approximate one-particle density
  matrix from the local exchange-energy density and prints per-atom-pair
  exchange energies and bond orders, with an optional natural-orbital
  analysis (`NATORB`).
- **OSLO extended to open-shell systems** — OSLO orbital analysis now
  fully supports unrestricted (UHF/UKS) wavefunctions, not just
  restricted ones.
- **Orbitals exportable to `.fchk`** — OSLO orbitals, and now also the
  Effective Atomic Orbitals (EFOs) from EOS and GEOS, can be written out
  to a standard `.fchk` file, so they can be viewed directly in any
  common orbital-viewer instead of only as cube files.
- **New cube-file options** — `SPACING`/`RADIUS_SCALE` keywords under
  `# CUBE` for finer control over the generated grid, backed by one
  unified, faster cube-generation routine.

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
  compatible compiler/OpenBLAS automatically; an automated regression
  test suite now runs on every change, checked against independently
  computed reference values, and runs automatically (via GitHub Actions)
  on every update to the code.
- **New documentation website** — a full documentation site is now
  available at [apost3d.readthedocs.io](https://apost3d.readthedocs.io).
- **More consistent program output** — formatting of the printed output
  was made more uniform and readable throughout.
