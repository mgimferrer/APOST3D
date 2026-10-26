# Output style guide

How the `.apost` output is laid out. New output follows these rules;
sections written before the guide are brought in line by the
[clean-up protocol](cleanup-protocol.md).

## General rules

- Only formatted writes. List-directed output (`write(*,*)`, `print *`)
  adds a leading blank, prints numbers at full precision and cannot be
  aligned.
- Text starts in column 3. Lines have no trailing blanks.
- An empty line is printed as a truly empty line (`write(*,'(a)') ''`),
  never as a line of spaces. Blocks are separated by exactly one empty
  line.
- Lines stay within 80 columns where possible.

## Headings

There are three levels, and nothing else (no `****` lines, no hand-made
boxes):

| Level | Routine | Use |
|---|---|---|
| Section | `print_box('TITLE')` | one analysis or one major step |
| Subsection | `print_box('TITLE')`, or `print_subbox('TITLE')` right after another heading | a part of a section (a spin channel, an iteration, a table) |
| Caption | a plain line, sentence case | the title of a single table or list |

- Section and subsection titles are in capitals. `print_box` prints an
  empty line before and after the box; `print_subbox` prints the same box
  without the empty line before, for a heading that directly follows
  another one. The caller adds none.
- Two headings never follow each other directly with the same level:
  the second one is a subsection.

## Values

A single value is printed as `label (unit) : value`, with the label
padded to 40 characters so that the colons of neighbouring lines line
up:

% source: tests/reference/H2O-Dimer-RHF.apost
```text
  Electron-nuclear attraction (au)        :    -434.8233244
  Integration error (kcal/mol)            :           0.00
  ...
  Kinetic energy (au)                     :     151.9545478
```

- The unit is given for every quantity that has one: `(au)` for
  energies in hartree, `(kcal/mol)`, `(bohr)`, `(%)`. Dimensionless
  quantities (electrons, bond orders, oxidation states) have none.
- A value printed in the input summary and again later in the run is
  printed the same way both times (same units, same decimals).
- Zero is printed as `0.00`, never `-0.00` (the `Makefile` compiles with
  `-fno-sign-zero`).
- In code, such a line is printed with `print_value` (module
  `output_mod`, `modules.f90`): `call print_value('Kinetic energy', x, 7,
  'au')` for a real with 7 decimals, `call print_value('Number of atoms',
  n)` for an integer, `call print_value('Functional', 'B3LYP')` for text.

## Decimals

The same kind of quantity always has the same number of decimals:

| Quantity | Decimals |
|---|---|
| Total and component energies (au) | 7 |
| Energy errors and differences (kcal/mol) | 2 |
| IQA atom and pair terms (au) | 6 |
| Populations, charges, spin populations, bond orders, valences | 6 |
| EFO occupations per fragment (net, gross) | 5 |
| EFO occupations, one per orbital | 4 |
| EOS electrons per fragment | 2 |
| EOS last occupied / first unoccupied occupation | 3 |
| Oxidation states | 2 |
| Reliability index R (%) | 3 |
| OSLO FOLI values | 5 |
| Fit and recovery percentages | 2 |

The SCF energy read from the `.fchk` file is printed with 10 decimals in
the wavefunction summary.

## Tables

% source: tests/reference/CH3F-NAO.apost
```text
  Frag.  Electrons  Last occ.  First unocc.
  -----------------------------------------
      1       5.00      0.720             -
      2       4.00      1.000         0.280
  -----------------------------------------
  Reliability index R (%)                 :         94.026
```

- A header line, a rule, the rows, a closing rule. Rules are exactly as
  wide as the table.
- Column headings are right-aligned over their numbers.
- A missing value is printed as `-`.
- A `Total` row says what it adds up when it is not the plain sum of the
  column (for example, a matrix total over the upper triangle including
  the diagonal).
- Matrices use the shared printers: `MPRINT2` (atom × atom) and
  `group_by_frag_mat` (fragment × fragment, with its total), which split
  wide matrices into blocks of six columns; `VPRINT` prints the rows of a
  per-atom table. All of them print 13 decimals instead of 6 with
  `FULLPRECISION`.
- Element symbols come from `elements_mod` (`element_symbol`,
  `element_label`).

## Messages

- Warnings: `WARNING: ` followed by a sentence saying what happened and
  what the user can do about it.
- Notes that are not problems: a plain sentence, no prefix.
- Errors: through `apost_stop` ([Errors](conventions.md#errors)).

## Words

The same thing has the same name everywhere:

| Use | Not |
|---|---|
| Fragment (`Frag.` in table headings) | Frg., frag, FRAGMENT in running text |
| EFO, EFOs (effective fragment orbitals) | eff-AO-s, eff-AOs, EFFAOs in text |
| electron assignment | assignation |
| neutral titles (`BOND ORDER MATRIX`), the atomic definition named once at the top of its section and as the heading of each value column (`TFVC`, `Mulliken`, `Lowdin-Dav.`) | `"FUZZY ATOMS"`, `APOST3D`, `apost3d`, `3D-space` |
| Lowdin, Lowdin-Davidson (ASCII in the output) | Löwdin, Davidson-Lowdin |
