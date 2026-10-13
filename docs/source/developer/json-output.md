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
| `analyses` | one object per analysis that was performed; empty in this version |
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
| `atomic_definition` | `keyword` (`"TFVC"`, `"MULLIKEN"`, ...), `name`, `space` (`"real"` or `"hilbert"`), `default` (`true` when no atomic definition was given and TFVC was used) |

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
the results file, and is described on this page then.

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
