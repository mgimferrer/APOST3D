# Results file (JSON)

Every run writes, next to its text output, a machine-readable results
file `<job>.apost.json` (the job name is the one given on the command
line, so `water.inp` and `water.fchk` give `water.apost.json`). It is
meant for programs: the APOST3Dview visualizer, scripts, workflows. The
`.apost` text is meant for people and its layout may change between
versions; the results file changes only as described under
[Versions](#versions).

This page is the reference for its content.

## When it is written

- The file is started right after the input and wavefunction files are
  opened, filled while the run proceeds, and completed at the end. It is
  written as `<job>.apost.json.tmp` and renamed when complete, so a file
  named `<job>.apost.json` is always complete and valid JSON. A results
  file left from an earlier run with the same name is removed first.
- A run that stops on an error (usually a `STOP` line at the end of the
  output) still completes its file, with `run.status` set to `"stopped"`
  and the message in `run.stop_message`; it holds the blocks written
  before the stop. A run that stops while reading the input has no
  `input`, `molecule`, `fragments` or `analyses`.
- A run that stops before the files are opened (missing input or
  wavefunction file) writes no results file.

## Conventions

- Atom and fragment numbers are 1-based, in the order of the `.fchk`
  file, as in the text output.
- Units are part of the key: `_au` (hartree), `_bohr`, `_kcal`
  (kcal/mol), `_percent`. The energies of the `enpart` block are in au
  without a suffix. Other keys without a suffix are dimensionless
  (electrons, bond orders, oxidation states).
- Reals are written with 17 significant digits, which read back as the
  same double-precision value the program computed, whatever the
  `FULLPRECISION` setting of the text output. A value that is not
  available is `null` (also a non-finite number).
- Per-atom lists have one entry per atom; per-fragment lists one entry
  per fragment, in `fragments` order. Matrices are lists of rows.
- Atom and fragment values are both written by APOST-3D: a reader never
  builds fragment values by summing atom values.
- Keys may be added in later versions: a reader must ignore keys it does
  not know.

## Top level

| Key | Content |
|---|---|
| `schema` | always `"apost3d-results"` |
| `schema_version` | `"major.minor"`, see [Versions](#versions) |
| `program` | `name`, `version` (`"5.0"`), `version_date` (release date of that version, ISO 8601: `"2026-09-22"`) |
| `wavefunction` | as read from the `.fchk` file (below) |
| `input` | the input as written, and the atomic definition used (below) |
| `molecule` | atoms and geometry (below) |
| `fragments` | the fragments of `# FRAGMENTS` (below) |
| `analyses` | one object per analysis that was performed (below) |
| `atomic_weights` | the weight functions of the BECKE, BECKE-RHO and TFVC atoms (below) |
| `index` | one entry per quantity that can be shown on atoms, atom pairs, fragments or fragment pairs (below) |
| `run` | how the run went (below) |

### `wavefunction`

| Key | Content |
|---|---|
| `formalism` | `"restricted"` or `"unrestricted"` |
| `method` | as in the output: `"Hartree-Fock"`, `"KS-DFT"`, `"CASSCF"`, ... |
| `n_basis_functions` | number of basis functions |
| `n_electrons`, `n_alpha`, `n_beta` | electron counts (valence electrons only with an ECP) |
| `energy_au` | the `SCF Energy` field of the `.fchk` file, or `null` when the file has none. For a correlated Gaussian calculation it is the energy of the reference SCF; the `.fchk` files written by `utils/apost3d.py` (pySCF) store the energy of the method itself there. Any change of this definition will be announced in the release notes. |

### `input`

| Key | Content |
|---|---|
| `text` | the whole input file, as written |
| `atomic_definition` | `keyword` (`"TFVC"`, `"MULLIKEN"`, ...), `name`, `short_name` (the name in the output's column headings and in the index titles: `"TFVC"`, `"Lowdin"`, `"Lowdin-Dav."`, ...), `space` (`"real"` or `"hilbert"`), `default` (`true` when no atomic definition was given and TFVC was used) |

### `molecule`

| Key | Content |
|---|---|
| `n_atoms` | number of atoms |
| `atomic_numbers` | atomic number of each atom |
| `symbols` | element symbol of each atom |
| `nuclear_charges` | nuclear charge of each atom (smaller than the atomic number for an atom with an ECP) |
| `coordinates_bohr` | x, y, z of each atom, in bohr |

### `fragments`

A list with one object per fragment, `{"id": 1, "atoms": [1, 2, 3]}`, in
the order of `# FRAGMENTS`; empty when the input defines no fragments.
Fragments never share an atom (the run stops). An atom can belong to no
fragment only when no fragment analysis needs every atom (EOS, GEOS,
OSLO and ENPART stop; the others warn in the output); `-1` on the last
fragment of the input means "all remaining atoms", so it is written out
as their list.

### `atomic_weights`

Written with the real-space atomic definitions `BECKE`, `BECKE-RHO` and
`TFVC`; absent with Hilbert-space atoms, `HIRSH` and `HIRSH-IT`, and in a
run that stops before its atoms are set up. With the coordinates of
`molecule` it defines the weight w_A(r) of every atom, so a reader can
rebuild the weight of a fragment, w_F = Σ_{A∈F} w_A, anywhere in space.

| Key | Content |
|---|---|
| `profile` | `"becke"` (`BECKE`, `BECKE-RHO`) or `"tfvc"` (`TFVC`): the cell function below |
| `stiffness` | k, the number of iterations of the cell function (`STIFFNESS`, default 4) |
| `radius_ratios` | χ_AB for every atom pair: row A, column B; the diagonal (1) is not used. With `BECKE` the ratio of the atomic radii, with `BECKE-RHO` and `TFVC` the ratio APOST-3D sets from the density along the A–B axis |

$$
w_A(\mathbf r) = \frac{P_A(\mathbf r)}{\sum_B P_B(\mathbf r)}, \qquad
P_A(\mathbf r) = \prod_{B \neq A} \tfrac12 \left(1 - f^{(k)}(\nu_{AB})\right), \qquad
f(x) = \tfrac12\, x\,(3 - x^2)
$$

with f applied k times, the distances in bohr,
μ_AB = (|r − R_A| − |r − R_B|) / |R_A − R_B| and χ = χ_AB:

- `becke`: ν_AB = μ_AB + a (1 − μ_AB²), a = (1 − χ²) / (4χ), limited to
  the range −½ to ½;
- `tfvc`: ν_AB = (1 + μ_AB − χ (1 − μ_AB)) / (1 + μ_AB + χ (1 − μ_AB)).

The orbitals of the EOS and GEOS `.fchk` files are basis-set fits of the
real-space EFOs, each the EFO times its fragment's weight (see
[Orbital `.fchk` files](../guide/visualization.md#orbital-fchk-files)); the fit puts small
tails on the other fragments. Multiplied once more by w_F, an orbital
loses those tails; it stays practically the same where w_F is close to 1
and comes out somewhat weaker near the fragment's boundary than the cube
file of the same EFO.

### `run`

| Key | Content |
|---|---|
| `job_name`, `input_file`, `wavefunction_file` | file names without folders (`"water"`, `"water.inp"`, `"water.fchk"`); the input, wavefunction and results files are always in the same folder |
| `date` | when the run started, ISO 8601 with time zone (`"2026-10-13T14:05:31+02:00"`) |
| `id` | the run id: the start time and 4 hex digits (`"2026-10-13T14:05:31+02:00/7f3a"`), also written in the title line of every orbital `.fchk` file of the run |
| `threads` | number of threads |
| `status` | `"normal"` or `"stopped"` |
| `stop_message` | the `STOP` message, or `null` |

## What goes into `analyses`

Each analysis writes the results it reports, the values a reader of the
output would quote or a visualizer would show: final numbers only. The
intermediate steps (for example, the ENPART matrices before the
zero-error correction), diagnostics and settings stay in the text
output. A block is added in the same version that adds the analysis to
the results file, and is described on this page then. A key is present
only when the run computes it; the conditions are given with each block.

### Matrix objects

An atom x atom or fragment x fragment matrix is an object:

| Key | Content |
|---|---|
| `values` | the full square matrix, a list of rows (also when symmetric) |
| `symmetric` | `true` when the matrix is symmetric |
| `diagonal`, `off_diagonal` | what the elements are, as a short noun for a legend (`"localization index"`, `"bond order"`, `"local spin"`, `"spin coupling"`): the same layout means different things in different analyses |
| `total`, `total_rule` | when the output prints a total: its value, and how it is formed (`"diagonal and one triangle"` or `"all elements"`) |

Fragment matrices always carry `total` and `total_rule`.

### Tables

A table (several quantities per fragment or per orbital) is an object:

| Key | Content |
|---|---|
| `rows` | `"fragment"` (row i is fragment i) or `"orbital"` (one row per orbital, in the order given) |
| `n_rows` | number of rows |
| `columns` | one object per column: `key`, `title`, `unit` (or `null`), `decimals`, `signed` (as in the [index](#index)) |
| `data` | one list per column, under its `key`; `null` where a row has no value. A column may hold a list per row (`fragments`, `shares`) |
| `fchk_orbitals` | in a table of orbitals written to a `.fchk` file: `"alpha"` or `"beta"`, the block of the file its `fchk_column` values refer to |

**Orbital tables** (`efos`, `negative_efos`, `oslos`) have one row per
orbital and the same two columns for its fragments: `fragments`, a list
of fragment numbers, and `shares`, the share of each. An EFO always has
one fragment and share 1 (`[3]`, `[1.0]`); a shared OSLO has several
(`[1, 2]`, `[0.5, 0.5]`). A table whose orbitals are in a `.fchk` file
also has `fchk_orbitals` and a `fchk_column` column (1-based; `null` for
an orbital the file does not hold). The block of such a table gives the
file name (`orbital_file`) and `orbital_file_id`, equal to the first line
of the file after `run `: a file from another run (same molecule,
other fragments or atomic definition) can be refused by comparing them.

### `populations`

Written by every run (the population analysis is always done).

| Key | Content |
|---|---|
| `scheme` | the atomic definition, as `input.atomic_definition.keyword` |
| `atoms` | `electrons`, `charges` and, for open-shell wavefunctions and correlated ones with a spin density, `spin`: one value per atom |
| `fragments` | the same per fragment (only with fragments) |
| `mulliken` | `atoms` and `fragments` as above with Mulliken's definition, printed for comparison (absent when the atomic definition is Mulliken) |

### `bond_orders`

Written by every run, except for a non-singlet correlated wavefunction
without `# DM` (no bond orders printed).

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `atoms` | `matrix`: bond orders (localization index on the diagonal; with `# DM`, the exchange-only bond orders of the 1-RDM); `valence_total`, `valence_bonded` (sum of the bond orders of the atom), `valence_free`; with `# DM` also `improved_bond_orders` (I. Mayer's improved bond orders, from the natural orbitals) with its total |
| `fragments` | `matrix` (localization index of the fragment on the diagonal, bond orders between fragments off it, total over the diagonal and one triangle), the three valences summed per fragment and, with `# DM`, `improved_bond_orders` (only with fragments) |
| `overlap_populations` | `atoms` and `fragments` (only with fragments), each with a `matrix` of overlap populations; written when the output prints them (EOS and GEOS, except with Lowdin, Lowdin-Davidson or NAO atoms) |

### `local_spin`

Written with `SPIN` (or `DM 2`) for an open-shell single determinant or
a correlated wavefunction.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `wavefunction` | `"single determinant"` or `"correlated"` |
| `total_unpaired_electrons` | number of effectively unpaired electrons N_D |
| `total_s2` | total of the <S^2> decomposition |
| `atoms` | `unpaired_electrons` (one value per atom) and `s2`, the a=3/4 decomposition of <S^2> (atomic terms on the diagonal, diatomic ones off it); for a correlated wavefunction also `xc_bond_orders` (exchange-correlation bond orders from the 2-RDM: delocalization indices, localization indices on the diagonal) with its total |
| `fragments` | `unpaired_electrons`, `s2` and, for a correlated wavefunction, `xc_bond_orders`, per fragment (only with fragments) |

### `effao`

Written with `EFFAO` or `UEFFAO` without `EOS`, and with `EFFAO-U`.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `total`, `alpha` and `beta`, or `paired` and `unpaired` | one object per density the EFOs come from (`total` for `EFFAO`, `paired`/`unpaired` for `EFFAO-U`), each with an `efos` table (below) |

The `efos` table has one row per EFO above the occupation cutoff,
fragment by fragment: `fragments`, `shares`, `net` (net occupation) and, in real
space, `gross` (gross occupation) and, except for `EFFAO-U`,
`fit_percent` (how closely the EFO written to a `.fchk` file reproduces
the fragment-weighted orbital).
EFFAO writes no orbital file, so these tables have no `fchk_column`
and their EFOs cannot be drawn from the results file.

### `eos`

Written with `EOS`.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `alpha` | `reliability_percent` (R of the channel), `assignment` (`"integer"`, or `"fractional"` when electrons are shared among EFOs of equal occupation at the frontier), `frontier` (table per fragment: `last_occupied`, `first_unoccupied`, the gross occupations of its last occupied and first unoccupied EFOs, `null` when it has none), `frontier_efos` (the EFOs that set R, below), `efos` (below) |
| `beta` | as `alpha`; or `{"same_as_alpha": true}` when the beta density equals the alpha one (closed shell), or `{"no_electrons": true}` |
| `fragments` | `oxidation_states`, `electrons_alpha`, `electrons_beta` |
| `reliability_percent` | the overall R, the smaller of the two channels |
| `orbital_file`, `orbital_file_id` | the `.fchk` file with the EFOs (in the same folder) and its id (see [Tables](#tables)) |

The `efos` table of a channel lists every EFO of the assignment, sorted
by gross occupation as they are assigned, which is also the order of the
orbitals of that channel in `orbital_file`: `fragments`, `shares`, `net`, and in
real space `gross` and `fit_percent` as in `effao`, `electrons` (given
to that EFO: 1, 0 or a fraction) and `fchk_column` (its orbital in the
file, `null` past the last one when the basis has linear dependencies).
For a closed shell, the file has only alpha orbitals.

`frontier_efos` gives the frontier EFOs of the whole channel, the ones
that set its R, as two lists, `last_occupied` and `first_unoccupied`.
Each item is an object with `efo_row` (1-based row of the channel's
`efos` table, so its `fchk_column` follows), `fragment` and `occupation`
(gross occupation, as in `frontier`). R compares the fragments' frontier
EFOs of the `frontier` table: λ_LO is the smallest last occupied
occupation, λ_FU the largest first unoccupied occupation of the other
fragments (an EFO of the fragment that gives λ_LO does not count:
moving an electron inside a fragment does not change its oxidation
state). `last_occupied` lists the fragments' last occupied EFOs within
`EOS_THRESH` (default 0.0025) of λ_LO, `first_unoccupied` the first
unoccupied EFOs of the other fragments within `EOS_THRESH` of λ_FU, so
that degenerate EFOs (equivalent fragments) all appear. With λ_LO the
smallest occupation of `last_occupied` and λ_FU the largest of
`first_unoccupied`, R = 100 min(1, λ_LO − λ_FU + 0.5) gives
`reliability_percent`. An empty list means R took its default: λ_LO = 1
(no fragment has an occupied EFO), λ_FU = 0 (no other fragment has an
unoccupied EFO, e.g. H2 with one EFO per atom, each holding half an
electron).

### `geos`

Written with `GEOS`. As `eos`, with the paired and unpaired densities in
place of the alpha and beta electrons:

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `paired`, `unpaired` | `reliability_percent`, `frontier`, `frontier_efos` and `efos` as in `eos` (for `paired`, occupations run from 0 to 2: λ_LO = 2 when `last_occupied` is empty, and R = 100 min(1, (λ_LO − λ_FU)/2 + 0.5)). The `efos` rows follow the orbitals of the channel in `orbital_file`: paired EFOs are its alpha orbitals, unpaired ones its beta orbitals; `electrons` is 2 or 0 for a paired EFO and 1 or 0 for an unpaired one |
| `negative_efos` | table of the paired EFOs with a clearly negative net occupation, left out of the assignment: `fragments`, `shares`, `net`, `gross`, `fit_percent` (as in the EFO tables: how closely the orbital in the file reproduces it) and `fchk_column`, their alpha orbital in `orbital_file` (the last ones of the block); no rows when there are none |
| `fragments` | `oxidation_states`, `electrons_paired`, `electrons_unpaired` |
| `reliability_percent` | the overall R |
| `rmsd` | root-mean-square deviation between the pooled occupations and the electrons assigned to them |
| `orbital_file`, `orbital_file_id` | the `.fchk` file with the EFOs and its id |

To label the orbitals of `orbital_file`: alpha orbital k is the row with
`fchk_column` k of a table with `fchk_orbitals` `"alpha"` (`paired.efos`
or `negative_efos`, the last orbitals of the block); beta orbital k the
one of `unpaired.efos`. Orbitals with no row are empty (zero
coefficients). The same rule holds for every orbital file.

### `oslo`

Written with `OSLO`.

| Key | Content |
|---|---|
| `scheme` | the atomic definition of the spreads, as in `populations` |
| `fragment_populations` | the atoms of the fragment populations (`# OSLO` `MULLIKEN`, `LOWDIN`, ...; else as `scheme`) |
| `total`, or `alpha` and `beta` | one object per channel (`total` for a restricted wavefunction, each orbital doubly occupied): `sum_foli` and the `oslos` table (below) |
| `fragments` | `oxidation_states`, and `electrons` (restricted) or `electrons_alpha`, `electrons_beta` |
| `orbital_file`, `orbital_file_id` | the `.fchk` file with the final OSLOs and its id |
| `orbital_file_preortho` | with `PRINT_NONORTHO`: the `.fchk` file with the OSLOs before orthogonalization, the orbitals the FOLI values were computed for. Same order and run id as `orbital_file`, so the `oslos` rows label it too |

The `oslos` table has one row per OSLO, in the order they were selected,
which is the order of the orbitals of the channel in `orbital_file`
(alpha orbitals for `total` and `alpha`, beta ones for `beta`):
`iteration` (that selected it), `foli` (also stored as its orbital
energy in the file), `fragments` (a list: the fragment it was assigned
to, or every fragment of a shared orbital), `shares` (a list, as
`fragments`: 1 for a single fragment; for a shared orbital the shares
printed in the output, divided equally among the orbitals of the
selection) and `fchk_column`. The `-OSLOs-preortho.fchk` file, when
written, has the same orbitals before orthogonalization in the same
order.

### `enpart`

Written with `ENPART`, for Hartree-Fock, KS-DFT and correlated
(CASSCF/CISD) wavefunctions. Energies in au; every matrix
is a [matrix object](#matrix-objects) in the convention of the text
output: one-center terms on the diagonal, the whole A–B term off it, the
total over the diagonal and one triangle. The external-field analysis
(per-atom terms) and `POLAR` are not written; only the field interaction
energy, needed to check the total.

| Key | Content |
|---|---|
| `wavefunction`, `formalism` | `"Hartree-Fock"`, `"KS-DFT"`, `"CASSCF"` or `"CISD"` (pySCF FCI and CISD wavefunctions are read as a CAS with every orbital active, so they say `"CASSCF"`); `"restricted"` or `"unrestricted"` |
| `functional` | KS-DFT: `keyword` (the predefined keyword of `# ENPART`, or `"LIBRARY"`), `components` (one object per libxc functional: `libxc_id`, `name`, `kind`, `family`) and `exact_exchange_fraction` (0 for a pure functional); `null` for Hartree-Fock and correlated wavefunctions |
| `pair_threshold`, `multipolar_pairs` | the `THREBOD` bond order, and the atom pairs below it (a list of `[A, B]`, 1-based), whose exchange-correlation term comes from the multipolar expansion |
| `electron_nuclear`, `kinetic`, `nuclear_repulsion` | one-electron terms (below); `kinetic` is a diagonal matrix |
| `xc_dft` | KS-DFT: the DFT part of the exchange-correlation (for a hybrid, without the exact exchange) |
| `two_electron_integration` | `energy` and `error_kcal` on the first grid, `tolerance_kcal`, `zero_error_correction` (`"applied"`, `"not needed"` or `"no reference energy"`), and when applied `rotated_grid_error_kcal`, `damping`, `corrected_energy`, `corrected_error_kcal` (`null` otherwise); `error_kcal` and `tolerance_kcal` are `null` with `"no reference energy"` |
| `coulomb` | the Coulomb term |
| `exact_exchange_share` | hybrid functionals: the HF-type exchange multiplied by the exact-exchange fraction, as it enters the energy |
| `xc` | the final exchange-correlation term: the exchange for Hartree-Fock, `xc_dft` for a pure functional, `xc_dft` + `exact_exchange_share` for a hybrid, the cumulant term for a correlated wavefunction |
| `exchange`, `correlation` | correlated wavefunctions with `CORRELATION`, only when the zero-error correction was not applied (it corrects their sum `xc`): the exchange and correlation parts, adding up to `xc` |
| `total` | the IQA matrix: self-energies on the diagonal, interaction energies off it |
| `total_energy`, `reference_energy`, `field_energy`, `integration_error`, `integration_error_kcal` | the total of `total`, the energy of the `.fchk` (which includes the interaction with an external field, when there is one), the field interaction energy (`null` without a field; it is not part of the IQA terms), and `integration_error` = `total_energy` + `field_energy` − `reference_energy` |
| `atomic_energies` | `atoms` and, with `DOFRAGS`, `fragments`: `self_energy` and `additive_energy` (self-energy + half of each interaction energy; they add up to `total_energy`) |
| `classical_interaction` | the electron-nuclear + nuclear-nuclear + Coulomb part of each interaction energy (`total` − `xc` off the diagonal; the diagonal is zero) |

Each term (`electron_nuclear` ... `total`, `exchange`, `correlation`,
`classical_interaction`) is an
object with `atoms` (atom matrix), `fragments` (fragment matrix, with
`DOFRAGS`) and, for `electron_nuclear` and `kinetic` when the `.fchk` has
the reference energies, `integration_error_kcal`. All values are the
final ones, after the zero-error correction when it was applied; the
matrices of `VERBOSE 1` are not written.

## Index

`index` lists every quantity a program can show without knowing the
analysis: single numbers and values on atoms, atom pairs, fragments and
fragment pairs. One
entry per quantity, in the order they are written:

| Key | Content |
|---|---|
| `id` | stable name, the path below `analyses` joined with dots (`"populations.atoms.charges"`) |
| `analysis` | the key of its block in `analyses` |
| `spin` | the channel of a quantity given per channel: `"alpha"`, `"beta"`, `"paired"`, `"unpaired"`, or `"total"` for the combined result; `null` otherwise |
| `path` | JSON Pointer to the value (`"/analyses/populations/atoms/charges"`); for a matrix, to its matrix object |
| `role` | `"comparison"` for a value printed only to compare with the result (the Mulliken values of a run with another atomic definition); absent otherwise |
| `kind` | `scalar` (one number, e.g. `local_spin.total_s2`), `atom_vector`, `atom_matrix`, `fragment_vector`, `fragment_matrix` or `table` |
| `title` | a name for a menu, with the atomic definition (`"Atomic charges (TFVC)"`) |
| `unit` | `"au"`, `"%"`, ... or `null` for a dimensionless quantity |
| `decimals` | decimals the text output prints for it (2 to 7; 13 for the 6- and 7-decimal values with `FULLPRECISION`): a default for display; `null` for a table, whose columns give their own |
| `signed` | `true` when the values can be negative by nature (charges, spin populations, the <S^2> decomposition, energies, oxidation states, the occupations of negative EFOs, and Mulliken bond orders and overlap populations): a diverging colour scale centred on zero suits them; `false` when they are not (populations, bond orders, valences), even if a value comes out slightly negative from the numerical integration |

`index` is always present (empty when the run stopped before any
analysis).

## Versions

`schema_version` is `"major.minor"`. A new key or block increases the
minor number; a renamed or removed key, or a change of meaning, increases
the major number. Version 0.x is the draft that becomes 1.0 with
APOST-3D 5.0; while it is a draft a key may still be renamed with only a
minor increase (0.8 replaced the `fragment` column of the EFO tables by
`fragments` and `shares`; 0.13 renamed `recovered_percent` of
`negative_efos` to `fit_percent`; 0.14 moved `new_bond_orders` from
`local_spin` to `bond_orders` as `improved_bond_orders` and renamed `li_di` to `xc_bond_orders`).

## Checking a change

The test suite keeps a reference results file for every test
(`tests/reference/<test>.apost.json`) and compares each run with it
(`tests/compare_json.py`: same structure, text and integers, reals within
a relative 1e-9 or an absolute 1e-8 au, so that values that cancel to
almost zero compare across machines; keys ending in `_kcal` get the same
1e-8 au in kcal/mol, 6.3e-6). `make update-ref` rewrites them together with the
`.apost` references.
