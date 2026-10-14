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
- A run that stops on an error (a `STOP` line at the end of the output)
  still completes its file, with `run.status` set to `"stopped"` and the
  message in `run.stop_message`; it holds the blocks written before the
  stop.
- A run that stops before the files are opened (missing input or
  wavefunction file) writes no results file.

## Conventions

- Atom and fragment numbers are 1-based, in the order of the `.fchk`
  file, as in the text output.
- Units are part of the key: `_au` (hartree), `_bohr`, `_kcal_mol`,
  `_percent`. No suffix means dimensionless (electrons, bond orders,
  oxidation states).
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
| `index` | one entry per quantity that can be shown on atoms, atom pairs, fragments or fragment pairs (below) |
| `run` | how the run went (below) |

### `wavefunction`

| Key | Content |
|---|---|
| `formalism` | `"restricted"` or `"unrestricted"` |
| `method` | as in the output: `"Hartree-Fock"`, `"KS-DFT"`, `"CASSCF"`, ... |
| `n_basis_functions` | number of basis functions |
| `n_electrons`, `n_alpha`, `n_beta` | electron counts (valence electrons only with an ECP) |
| `energy_au` | the `SCF Energy` field of the `.fchk` file (for a correlated wavefunction, the energy of its reference SCF, not the correlated total energy), or `null` when the file has none. Any change of this definition will be announced in the release notes. |

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

### `run`

| Key | Content |
|---|---|
| `job_name`, `input_file`, `wavefunction_file` | file names without folders (`"water"`, `"water.inp"`, `"water.fchk"`); the input, wavefunction and results files are always in the same folder |
| `date` | when the run started, ISO 8601 with time zone (`"2026-10-13T14:05:31+02:00"`) |
| `threads` | number of threads |
| `status` | `"normal"` or `"stopped"` |
| `stop_message` | the `STOP` message, or `null` |

## What goes into `analyses`

Each analysis writes the results it reports, the values a reader of the
output would quote or a visualizer would show: final numbers only. The
intermediate steps (for example, energy terms before the zero-error
interpolation of ENPART), diagnostics and settings stay in the text
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
| `data` | one list per column, under its `key`; `null` where a row has no value |

### `populations`

Written by every run (the population analysis is always done).

| Key | Content |
|---|---|
| `scheme` | the atomic definition, as `input.atomic_definition.keyword` |
| `atoms` | `electrons`, `charges` and, for open-shell wavefunctions and correlated ones with a spin density, `spin`: one value per atom |
| `fragments` | the same per fragment (only with fragments) |
| `mulliken` | `atoms` and `fragments` as above with Mulliken's definition, printed for comparison (absent when the atomic definition is Mulliken) |

### `bond_orders`

Written by every run.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `atoms` | `matrix`: bond orders (localization index on the diagonal); `valence_total`, `valence_bonded` (sum of the bond orders of the atom), `valence_free` |
| `fragments` | `matrix` (localization index of the fragment on the diagonal, bond orders between fragments off it, total over the diagonal and one triangle) and the three valences summed per fragment (only with fragments) |
| `overlap_populations` | `atoms` and `fragments` (only with fragments), each with a `matrix` of overlap populations; written when the output prints them (EOS, GEOS) |

### `local_spin`

Written with `SPIN` (or `DM 2`) for an open-shell single determinant or
a correlated wavefunction.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `wavefunction` | `"single determinant"` or `"correlated"` |
| `total_unpaired_electrons` | number of effectively unpaired electrons N_D |
| `total_s2` | total of the <S^2> decomposition |
| `atoms` | `unpaired_electrons` (one value per atom) and `s2`, the a=3/4 decomposition of <S^2> (atomic terms on the diagonal, diatomic ones off it); for a correlated wavefunction also `new_bond_orders` (bond orders from the natural orbitals, I. Mayer's definition) and `li_di` (localization and delocalization indices from the 2-RDM), both with their totals |
| `fragments` | `unpaired_electrons`, `s2` and, for a correlated wavefunction, `li_di`, per fragment (only with fragments) |

The bond orders from the 1-RDM that a correlated local spin analysis
prints again are those of `bond_orders`.

### `effao`

Written with `EFFAO` or `UEFFAO` without `EOS`.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `total`, or `alpha` and `beta` | one object per density the EFOs come from (`total` for `EFFAO`), each with an `efos` table (below) |

The `efos` table has one row per EFO above the occupation cutoff,
fragment by fragment: `fragment`, `net` (net occupation) and, in real
space, `gross` (gross occupation) and `fit_percent` (how closely the EFO
written to a `.fchk` file reproduces the fragment-weighted orbital).

### `eos`

Written with `EOS`.

| Key | Content |
|---|---|
| `scheme` | as in `populations` |
| `alpha` | `reliability_percent` (R of the channel), `assignment` (`"integer"`, or `"fractional"` when electrons are shared among EFOs of equal occupation at the frontier), `frontier` (table per fragment: `last_occupied`, `first_unoccupied`, the gross occupations of its last occupied and first unoccupied EFOs, `null` when it has none), `efos` (below) |
| `beta` | as `alpha`; or `{"same_as_alpha": true}` when the beta density equals the alpha one (closed shell), or `{"no_electrons": true}` |
| `fragments` | `oxidation_states`, `electrons_alpha`, `electrons_beta` |
| `reliability_percent` | the overall R, the smaller of the two channels |
| `orbital_file` | the `.fchk` file with the EFOs (in the same folder) |

The `efos` table of a channel lists every EFO of the assignment, sorted
by gross occupation as they are assigned, which is also the order of the
orbitals of that channel in `orbital_file`: `fragment`, `net`, and in
real space `gross` and `fit_percent` as in `effao`, `electrons` (given
to that EFO: 1, 0 or a fraction) and `fchk_column` (its orbital in the
file, `null` past the last one when the basis has linear dependencies).
For a closed shell, the file has only alpha orbitals.

## Index

`index` lists every quantity a program can show without knowing the
analysis: single numbers and values on atoms, atom pairs, fragments and
fragment pairs. One
entry per quantity, in the order they are written:

| Key | Content |
|---|---|
| `id` | stable name, the path below `analyses` joined with dots (`"populations.atoms.charges"`) |
| `analysis` | the key of its block in `analyses` |
| `spin` | `"alpha"`, `"beta"` or `"total"` for a quantity given per spin channel; `null` otherwise |
| `path` | JSON Pointer to the value (`"/analyses/populations/atoms/charges"`); for a matrix, to its matrix object |
| `role` | `"comparison"` for a value printed only to compare with the result (the Mulliken values of a run with another atomic definition); absent otherwise |
| `kind` | `scalar` (one number, e.g. `local_spin.total_s2`), `atom_vector`, `atom_matrix`, `fragment_vector`, `fragment_matrix` or `table` |
| `title` | a name for a menu, with the atomic definition (`"Atomic charges (TFVC)"`) |
| `unit` | `"au"`, `"%"`, ... or `null` for a dimensionless quantity |
| `decimals` | decimals the text output prints (6, or 13 with `FULLPRECISION`): a default for display; `null` for a table, whose columns give their own |
| `signed` | `true` when the values can be negative by nature (charges, spin populations, the <S^2> decomposition, and Mulliken bond orders and overlap populations): a diverging colour scale centred on zero suits them; `false` when they are not (populations, bond orders, valences), even if a value comes out slightly negative from the numerical integration |

`index` is always present (empty when the run stopped before any
analysis).

## Versions

`schema_version` is `"major.minor"`. A new key or block increases the
minor number; a renamed or removed key, or a change of meaning, increases
the major number. Version 0.x is the draft that becomes 1.0 with
APOST-3D 5.0.

## Checking a change

The test suite keeps a reference results file for every test
(`tests/reference/<test>.apost.json`) and compares each run with it
(`tests/compare_json.py`: same structure, text and integers, reals within
a relative 1e-9). `make update-ref` rewrites them together with the
`.apost` references.
