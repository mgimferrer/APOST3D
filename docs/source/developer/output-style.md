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
| Subsection | `print_subbox('TITLE')` | a part of a section (a spin channel, an iteration, a table) |
| Caption | a plain line, sentence case | the title of a single table or list |

- Section and subsection titles are in capitals; `print_box` and
  `print_subbox` print their own empty lines before and after, so the
  caller adds none.
- Two headings never follow each other directly with the same level:
  the second one is a subsection.

## Values

A single value is printed as `label (unit) : value`, with the label
padded to 40 characters so that the colons of neighbouring lines line
up:

```text
  Electron-nuclear attraction (au)        :   -194.5024157
  Kinetic energy (au)                     :     75.8390458
  Integration error (kcal/mol)            :      0.00
```

- The unit is given for every quantity that has one: `(au)` for
  energies in hartree, `(kcal/mol)`, `(bohr)`, `(%)`. Dimensionless
  quantities (electrons, bond orders, oxidation states) have none.
- A value printed in the input summary and again later in the run is
  printed the same way both times (same units, same decimals).
- Zero is printed as `0.00`, never `-0.00`.

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

## Tables

```text
  Frag.  Electrons  Last occ.  First unocc.
  -----------------------------------------
      1      12.00      0.837         0.321
      2       3.00      0.512             -
  -----------------------------------------
  Reliability index R (%)                 :     96.598
```

- A header line, a rule, the rows, a closing rule. Rules are exactly as
  wide as the table.
- Column headings are right-aligned over their numbers.
- A missing value is printed as `-`.
- A `Total` row says what it adds up when it is not the plain sum of the
  column (for example, a matrix total over the upper triangle including
  the diagonal).
- Matrices (atom × atom, fragment × fragment) use the shared matrix
  printer, which splits wide matrices into blocks of six columns.

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
| the name of the atomic definition in titles (`TFVC BOND ORDER MATRIX`, `MULLIKEN BOND ORDER MATRIX`) | `"FUZZY ATOMS"`, `APOST3D` |
| Lowdin, Davidson-Lowdin (ASCII in the output) | Löwdin |
