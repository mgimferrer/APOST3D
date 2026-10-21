# Changelog

This file summarizes what changed in APOST-3D **Version 5.0** relative to
the previous public release (`master`). It is meant for users of the
program, not a full development history — see `git log` for that.

## [5.0] — Unreleased

### New implementations

- **Results file `<job>.apost.json`** — every run also writes its
  results in a machine-readable file, for the APOST3Dview visualizer and
  scripts: the wavefunction summary, the input, the molecule and the
  fragments, the populations, bond orders, valences, local spin
  analysis, effective fragment orbitals, EOS, GEOS and OSLO (with the
  orbitals of the EFO and OSLO `.fchk` files mapped to their fragments)
  and the ENPART energy decomposition (HF, KS-DFT, CASSCF/CISD), with
  the other analyses to be added. A run that
  stops on an error still writes a complete file that says so. See
  "Results file (JSON)" in the developer guide.
- **GEOS (Generalized Effective Oxidation States)** — a new open-shell
  oxidation-state analysis (`# METHOD GEOS`), built from paired/unpaired
  densities separately (Takatsuka's definition) rather than the total
  density alone, complementing the existing EOS method.
- **DFT-DM1** (in development, not documented yet) — a new analysis method (`DFTDM1` in `# METHOD`) for
  UHF/UKS-DFT wavefunctions: builds an approximate one-particle density
  matrix from the local exchange-energy density and prints per-atom-pair
  exchange energies and bond orders, with an optional natural-orbital
  analysis (`NATORB`).
- **OSLO extended to open-shell systems** — OSLO orbital analysis now
  fully supports unrestricted (UHF/UKS) wavefunctions, not just
  restricted ones.
- **OSLO branching** — `BRANCH_ITERATION` (and `BRANCH_ALPHA`/`BRANCH_BETA`
  for open shells) in `# OSLO` selects the next FOLI value at the given
  iterations, to follow an alternative assignment; it used to stop the
  run. A warning now names the iteration when the next candidates are
  nearly the selected orbitals with a Δ-FOLI below 0.5, and each run
  prints the sum of its FOLI values to compare the alternatives.
- **Orbitals exportable to `.fchk`** — OSLO orbitals, and now also the
  effective fragment orbitals (EFOs) from EOS (with any population
  scheme) and GEOS, are written out to a standard `.fchk` file, so they
  can be viewed directly in any common orbital-viewer instead of only as
  cube files.
- **New cube-file options** — `SPACING`/`RADIUS_SCALE` keywords under
  `# CUBE` for finer control over the generated grid, backed by one
  unified, faster cube-generation routine. `MAX_OCC 1000` (the default)
  sets no upper limit, so fully occupied orbitals always get a cube file.
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
- **Meta-GGA functionals in ENPART** — meta-GGAs and their global hybrids
  (TPSS, M06-2X, ...), restricted and unrestricted. The pair terms take
  the kinetic energy density of the bond order density, built like its
  gradient. The total exchange-correlation energy reproduces pySCF's to
  0.02 kcal/mol.

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
- **Clear stops instead of wrong numbers** — range-separated hybrids,
  VV10 functionals, invalid libxc ids, and ENPART with a Hilbert-space
  atom definition (`MULLIKEN`, `LOWDIN`) now stop with a message before any energy is computed. So do ENPART of a
  correlated (CASSCF/CISD) wavefunction computed in an electric field, and
  `POLAR` with a correlated wavefunction, whose untested code was removed.
- **Clearer output** — the section starts with the wavefunction, the
  atomic definition and how to read the matrices; every term is printed
  once, with its final values: the two-electron matrices follow a
  `TWO-ELECTRON INTEGRATION ERROR` section (error, tolerance, whether the
  zero-error correction is applied, damping parameter, corrected
  energy), instead of being printed before and again after the
  correction. For hybrid functionals the exact exchange is printed
  already multiplied by the exact-exchange fraction (`EXACT-EXCHANGE
  SHARE`), so `DFT EXCHANGE-CORRELATION ENERGY` + `EXACT-EXCHANGE SHARE` =
  `EXCHANGE-CORRELATION ENERGY`. The kinetic energy is a column per atom
  (it has no pair terms), the KS-DFT one-center terms are one table per
  atom, and with `DOFRAGS` every printed matrix is followed by its
  fragment matrix (the HF exchange and the hybrid exchange-correlation
  had none). The final matrix is `IQA ENERGY TERMS (TFVC)` (or the
  atomic definition used) instead of `FUZZY ATOMS ... ENERGY
  COMPONENTS`. The intermediate matrices are printed with `VERBOSE 1`.
  CASSCF/CISD wavefunctions follow the same layout; with `CORRELATION`,
  exchange and correlation are printed separately only when the
  zero-error correction (which applies to their sum) is not applied, or
  with `VERBOSE 1` before it. The pair count of the correlated part, which
  counted grid points, is now the number of atom pairs. No number changes.
- **Additive atomic energies and classical interaction energies** —
  printed after the final decomposition, from its final matrices:
  `ATOMIC ENERGIES` (self-energy and additive energy of each atom, which
  add up to the total energy) and `CLASSICAL INTERACTION ENERGIES` (the
  electron-nuclear, nuclear-nuclear and Coulomb part of each interaction
  energy), also per fragment.

### Other changes that affect inputs

- **`VERBOSE`** — new option in `# METHOD`: `VERBOSE 1` also prints the
  intermediate matrices of an analysis (for now ENPART: the BODEN of each
  pair, the DFT pair and one-center terms, the unscaled HF-type
  exchange, and the matrices before the zero-error correction). The
  default, `VERBOSE 0`, prints the final results only.
- **EOS and GEOS need fragments** — like OSLO, they now stop without
  `DOFRAGS` and a `# FRAGMENTS` block instead of giving every atom its own
  oxidation state by default. For one oxidation state per atom, define
  every atom as a fragment of its own.
- **Atomic definition** — always write exactly one in `# METHOD`. An
  input without any now runs with `TFVC` and a warning (it used plain
  Becke atoms with fixed radii); those are selected with the new `BECKE`
  keyword. An input with more than one now stops with a message (one of
  them was used without saying which).
- **Integration grids in the input** — the optional `# GRID` block now
  sets every grid and is always read: `RADIAL`/`ANGULAR` for the
  one-electron grid of every analysis (before only on the command line,
  `apost3d job 70 434`, which now stops with a message), `RADIAL_2E`/
  `ANGULAR_2E` for ENPART's two-electron grid. `MOD-GRIDTWOEL` and
  `FINEGRID` are gone (`ANGULAR 974` replaces the latter), the rotation
  angles `phb1`/`phb2` are chosen by the program for the angular grid
  (expert override: `ROTATION_2E`), and `rr00` is now `R0_2E`; the old
  spellings stop with a message. `THRESH2`, which had no effect, is gone
  too. Defaults unchanged.
- **Analyses not available in this version** — `EDAIQA`, `TOPOLOGY`,
  `DAFH`, `SHANNON`, `WFN` and `READINT` in `# METHOD`, and `DMRG` in
  `# DM`, were read but nothing was computed (the input summary still
  listed them); they now stop with a message saying so.
- **Keywords are whole words** — a keyword is no longer recognized inside a
  longer word (e.g. `EOS_THRESH` switched on `EOS`), a block header must be
  alone on its line, and the file names in `# DM` are never read as
  keywords. `FOLI TOLERANCE`, `BRANCH ITERATION` and `PRINT NON-ORTHO`
  (`# OSLO`) are now written `FOLI_TOLERANCE`, `BRANCH_ITERATION` and
  `PRINT_NONORTHO`; the old spellings stop with a message. `MULLIKEN` and
  `MULLI` both select Mulliken atoms.
- **`apost3d-eos`, `group_frag` and `eos_alt` removed** — `apost3d` does
  everything `apost3d-eos` did; ENPART with `DOFRAGS` prints the fragment
  terms `group_frag` computed; `eos_alt` did not work with version 5
  outputs.

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
- **OSLO with an orbital found from several fragments** — when the
  selected orbitals are the same orbital (e.g. the bond between two
  equivalent fragments), it is now assigned once and shared among those
  fragments, as in the OSLO paper. It used to be given whole to each of
  them (ethane as two CH₃ fragments: −1/−1 instead of 0/0), writing past
  the end of its arrays.
- **EOS with a changed orbital cutoff** — an `EFF_THRESH` that leaves
  fewer orbitals than electrons now stops instead of assigning electrons
  wrongly, very low or negative values no longer crash, and the overall
  reliability index R(%) is defined for systems without beta electrons.
- **EOS/GEOS reliability index R(%)** — the frontier orbitals behind it
  ("Last occ."/"First unocc.") are now taken in the same order as the
  electron assignment (gross occupation); one of them could be read from
  an empty slot. In pseudo-degenerate cases R(%) changes in the last
  digits; the R(%) formulas are unchanged.
- **`EFFAO-U` alone** — no longer crashes at the end of the run, and on a
  restricted single determinant its occupations are no longer empty.
- **EFFAO/EOS/GEOS fragment sums** — `Net occupation for fragment` is now
  the sum over the listed orbitals (above the cutoff) in every scheme, like
  the gross one (in real space it summed all orbitals, negative GEOS ones
  included), and a new `Left out by the cutoff` line gives what the listed
  orbitals miss of the fragment's population. It replaces the `Deviation
  from net/gross population` lines, which were not computed correctly.
  With Mulliken or Löwdin atoms, a fragment with a single basis function
  no longer loses its only orbital.
- **GEOS with a linearly dependent basis** — the negative paired EFOs
  were left out of `<jobname>-GEOS-EFOs.fchk` when the basis had near
  linear dependencies (fewer independent functions than basis
  functions, e.g. with diffuse functions); they are now its last alpha
  orbitals, as in any other basis.
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
  default), and the last block may end without its closing `#`. The
  "MO formalism" line of pySCF CASSCF, FCI and CCSD files no longer
  prints unreadable characters.
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
  not a Lebedev grid, or a radial one outside 1–500, now stops the run. In
  `# GRID` it used to give wrong ENPART numbers, hidden by the zero-error
  interpolation; on the command line it was silently rounded down.
- **ENPART with pseudopotentials** — no longer stops with `ECP Matrix not
  found`. `.fchk` files written from pySCF with `utils/apost3d.py` now
  carry the ECP matrix, and their reference electron-nuclear energy
  includes the ECP energy: it used to be counted in the electron-electron
  energy, which made the zero-error strategy distort the ENPART terms.
  `get_energy_g16`/`get_energy` append the ECP matrix for Gaussian
  `.fchk` files.
- **`get_energy_g16`/`get_energy`** — an incomplete Gaussian output file
  (calculation not ended normally) or one without the lines they need now
  gives an error message and adds nothing to the `.fchk` file; they used
  to crash or append wrong energies. They no longer have a limit of 3000
  basis functions. `get_energy` (Gaussian 09) will be removed in a future
  version.
- **Cube file names** — the `.cube` files of `# CUBE` keep the whole job
  name (it was cut to 30 characters in the file name and to 8 in the
  title line inside the file).
- **Samarium** — element 62 is printed as `Sm` (it was printed as `Sn`
  in the output tables, the cube and `.int` files, and in `eos_aom`).
- **Elements beyond uranium** — atom labels in the output, cube and
  `.int` files now cover every element up to oganesson (Z = 118); the
  tables stopped at uranium.
- **`.int` files of atoms with an ECP** — the `NET CHARGE` written by
  `DOINT` uses the nuclear charge left by the pseudopotential (it counted
  the core electrons the ECP replaces as a positive charge).
- **OSLO close-alternative warning** — no longer printed when the
  selection takes every orbital still left: every candidate then lies in
  the selected space, so the overlap was 1 by construction and the
  warning fired without a real alternative (e.g. an H2 triplet whose
  suggested branch gives the same oxidation states). The output says
  `Close alternatives: not checked, no orbital left after this
  selection`.
- **OSLO with Mulliken populations** — `# OSLO MULLIKEN` gave NaN FOLI
  values for orbitals with a slightly negative Mulliken population on a
  fragment, and selected them. An orbital whose population on a fragment
  is 10⁻⁵ or less is now not a candidate of that fragment.
- **libxc references** — functional names and references are printed in
  full (they were cut at 80 and 120 characters).
- **`.nao` file** — `NAO-BASIS` without the `jobname.nao` file stops with
  a message (it ended with a runtime error). The file is now checked
  against the wavefunction: its second line must read `NAOs in the AO
  basis:`, it must hold the whole matrix of the basis set, and the NAOs
  must be orthonormal in that basis (the largest error is printed; above
  10⁻⁵ the run stops). The file of another calculation used to end in a
  runtime error (smaller basis set) or give meaningless populations with
  a normal termination (larger basis set).
- **Iterative Hirshfeld** — the populations are converged to 10⁻⁶ (was
  5×10⁻⁴, which left the charges up to 7×10⁻⁴ from their converged
  values; CH3F takes 32 iterations instead of 13, at most 100 instead of
  40), and a run that does not converge says so. An atom whose charge goes
  beyond ±2, the last charge state of `densoutput`, no longer stops the
  run ("Absolute convergence unlikely"): its free-atom density is held at
  that state, with one warning. The iterations are printed as one table.
- **Cube file layout** — cube files no longer carry a line of orbital
  indices after the atoms. The cube format has that line only with a
  negative atom count, so strict readers (ASE, VMD) shifted the data by
  two values.
- **Cube norm** — the printed norm of each cube (now `Norm in cube`) used
  a voxel volume that was too small, making it 5–10% low. It is the share
  of the orbital inside the box, and a warning now flags values below
  0.9. Fragment and orbital numbers of any size go into the file names.
- **`DOINT` file names and check** — `.int` and `.files` names keep the
  whole job name (they were cut to 55 characters, the lines of the
  `.files` list to 40). The check that the atomic overlap matrices add up
  to the identity covers exactly the orbitals written: for correlated
  wavefunctions it printed false deviations. The output now gives the
  largest deviation, with a warning when it is too large.
- **Davidson spin matrix removed** — the `"FUZZY ATOMS" DAVIDSON SPIN DEC.
  MATRIX` printed by `SPIN` for open-shell single determinants had the
  spin-density term with the wrong sign. It is no longer printed; the
  local spin analysis (effectively unpaired electrons and the a=3/4
  decomposition of ⟨S²⟩) is unchanged.
- **Bond orders of correlated wavefunctions** — with `# DM`, the
  population analysis prints the exchange-only bond orders of the 1-RDM
  (`EXCHANGE-ONLY BOND ORDER MATRIX`), now with its spin-density term
  (missing for non-singlets), and the bond orders from the natural
  orbitals (`IMPROVED BOND ORDER MATRIX`, I. Mayer's improved bond
  orders, before printed only by `SPIN` as `NEW BOND ORDER MATRIX`). The
  local spin analysis no longer repeats the 1-RDM matrix; its LI/DI
  matrix is now `EXCHANGE-CORRELATION BOND ORDER MATRIX (DI)`. Without
  `# DM`, a non-singlet correlated wavefunction gets no bond orders or
  valences (they lacked the spin density), with a note. A non-singlet no
  longer gets the note "closed-shell correlated wavefunction".
- **POLAR field energies** — the "intrinsic" dipole energies used the
  origin-dependent atomic dipoles (so they already gave the whole
  **μ**·**F**) and the charge-transfer term was added on top, counting it
  twice (water RHF in a 0.005 au field: −0.0109182 au instead of
  **μ**·**F** = −0.0041738 au). The intrinsic terms now use the intrinsic
  dipoles, each atom also gets its charge-transfer term, and the total is
  **μ**·**F**. The POLAR output now has the program's layout (tables with
  totals, the bond charges, a molecular dipole summary, a note for
  charged molecules), and POLAR with a Hilbert-space atomic definition
  stops (its dipoles need real-space atoms).

### Others

- **`CUBE` without an EFO analysis** — the input summary warns that no
  cube file is written.
- **Thread count** — printed once, after the banner (`OpenMP threads`),
  instead of in the ENPART section together with the number of cores of
  the machine.
- **Start of the output** — the wavefunction and input summaries and the
  setup blocks (grids, atomic definition, density matrices of correlated
  wavefunctions) use one aligned `label (unit) : value` layout. The
  input summary shows values in the units used later in the run, with
  the keyword value next to them (e.g. `0.2000 (THREBOD 2000)`), and
  names every predefined functional. `SCF/DFT energy` reads `SCF energy`
  (and `not in the .fchk` instead of zero when the file has none); a
  method name the program does not recognise is shown next to the
  assumed KS-DFT. Real-space runs print how accurately the grid
  integrates the atomic overlaps.
- **Tables and matrices** — every atom and fragment matrix has rules as
  wide as the matrix and column headers over their numbers (the S^2
  decomposition and the correlated local-spin matrices now have rules
  too); population, valence and fragment tables line up with their
  headers; the total of a fragment matrix says what it adds up
  (diagonal and one triangle, or all elements). Blank lines are empty
  lines, and a value that rounds to zero is printed as `0.00`, never
  `-0.00`. No number changes.
- **Population analysis and local spin output** — the section starts
  with the atomic definition used, and the value columns are named after
  it (`TFVC`, `Lowdin`, `Hirshfeld`, ...) instead of `apost3d`/`3D-space`;
  titles no longer say `"FUZZY ATOMS"` or `APOST3D` (`BOND ORDER MATRIX`,
  `FRAGMENT ELECTRON POPULATIONS`, ...). A Mulliken run prints one
  column instead of two identical ones. The sum checks of the local spin
  analysis are `Total unpaired electrons N_D` and `Total <S^2>`, with 6
  decimals like the tables, and the fragment tables follow the atom
  tables they sum. No number changes.
- **EFFAO and EOS output** — the section starts with the atomic
  definition, the densities the orbitals come from and the occupation
  cutoff (printed once instead of in every fragment); each fragment lists
  its `Net` and `Gross` occupations and `Fit (%)` under clear labels; the
  assignment of each spin channel is one table with a `Reliability index
  R (%)` line, a fractional assignment is a warning line instead of a
  box, and the oxidation-state table ends with a `Total` row and the
  overall reliability index. A closed-shell run says that the beta
  electrons are assigned as the alpha ones, and a run without beta
  electrons says so. No number changes.
- **GEOS and EFFAO-U output** — the same layout as EFFAO and EOS
  (section header, `Net`/`Gross`/`Fit (%)` rows, one table per density
  and the shared oxidation-state table). The negative paired EFOs are a
  table after a warning line, each pair moved from the paired to the
  unpaired density is one row of a table with its RMSD and whether it
  was kept, the repeated list of electrons per fragment is gone, and the
  final RMSD has 5 decimals like the others. No number changes.
- **OSLO output and orbital files** — the OSLO section starts with the
  atomic definition, the atoms of the fragment populations, the FOLI
  tolerance and the charge centres (now labelled in Å); each iteration
  lists the candidate orbitals of every fragment, the selection, the
  Δ-FOLI and the close-alternative check in aligned tables, with
  warnings and branching as plain lines; the summary of the selected
  OSLOs gives the iteration that selected each, and the oxidation-state
  table a total. In `<jobname>-OSLOs.fchk` (and `-OSLOs-preortho.fchk`)
  the orbital energies are now the FOLI values of the OSLOs (they were
  the SCF orbital energies of the input file, unrelated to the OSLOs).
  No number changes.
- **`utils/apost3d.py` rewritten** — works with current pySCF (2.7 to 2.14
  checked; before, only up to 2.7) and writes basis sets with g functions,
  single atoms, linearly dependent bases, CASCI, the RDMs of FCI and
  closed-shell CCSD, and closed-shell CISD (with its RDMs, read as a CAS
  with every orbital active). An empirical dispersion energy (D3/D4) is removed
  from the reference energies for ENPART. Correlated wavefunctions on UHF
  orbitals, which APOST-3D cannot analyse correctly, are refused with a
  message, as is any other unsupported case, without leaving a partial
  file; the `.fchk` no longer contains bare header lines that strict
  readers reject. Existing scripts work unchanged.
- **Errors end with exit code 1** — a run that stops on an error now
  exits with code 1 (0 after a normal end), so job scripts and workflow
  tools can detect it, and its `STOP` message is the last line of the
  output (it used to appear above the output printed just before it).
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
  `wfn2fchk`, `eos_aom`).
- **libxc upgraded to 7.1.2** (from the previously bundled 4.2.3) — no
  change to any computed energy or property, verified bit-for-bit against
  reference outputs. Fetched and built automatically on first compile if
  no suitable existing install is found, rather than shipped as a bundled
  copy. An installed libxc is used only if it is version 7.0 or newer
  (newer releases need no change) and its Fortran interface works with the
  gfortran in use; `make_compile.sh` checks this, and every other
  prerequisite, before it starts compiling.
- **New documentation website** — a full documentation site is now
  available at [apost3d.readthedocs.io](https://apost3d.readthedocs.io).
- **More consistent program output** — formatting of the printed output
  was made more uniform and readable throughout.
