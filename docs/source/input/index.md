# Input file

A calculation reads two files with the same name: the wavefunction,
`jobname.fchk` (see [Preparing the wavefunction](../guide/wavefunctions.md)),
and the input file `jobname.inp`, which says what to compute. This section
lists every keyword of the input file, block by block.

## Structure

The input file is plain text made of **blocks**. Each block starts with a
header line (`# METHOD`, `# FRAGMENTS`, ...) and ends with a line holding
only `#`. Inside a block, each keyword goes on its own line, in any order:

```text
# METHOD
TFVC
EOS
DOFRAGS
#
# FRAGMENTS
2
1
1
-1
#
```

`# METHOD` is always required: it holds the [atomic definition](../guide/aim.md),
the analyses to run and the general options. The other blocks are needed
only by the keywords that use them:

| Block | Needed when | Page |
|---|---|---|
| `# METHOD` | always | [# METHOD](method.md) |
| `# FRAGMENTS` | `DOFRAGS` | [# FRAGMENTS](fragments.md) |
| `# ENPART` | `ENPART` | [# ENPART](enpart.md) |
| `# GRID` | `MOD-GRIDTWOEL` in `# ENPART` | [# GRID](grid.md) |
| `# OSLO` | `OSLO` | [# OSLO](oslo.md) |
| `# CUBE` | `CUBE` | [# CUBE](cube.md) |
| `# DM` | `DM 1` or `DM 2` | [# DM](dm.md) |

[Input examples](examples.md) shows complete inputs for the common
analyses.

## Rules

- **Keywords and block names are case-sensitive**: write them exactly as
  in this reference (`TFVC`, `# METHOD`, `rr00`, `pySCF`). The only
  exception are the functional names in `# ENPART` (`B3LYP`, `b3lyp`, ...).
- **Values** follow the keyword after a space or an `=`: `DM 2`, `DM=2`
  and `DM = 2` are the same. A keyword without its value takes the
  default; a value that can't be read (e.g. `THREBOD 2.5`, which must be
  an integer) stops the run, showing the line.
- **No comments inside a block.** Any line that contains a `#` ends the
  block, so a `## comment` line hides every keyword after it. Keywords
  written after the closing `#` of a block are ignored without warning.
- Keywords are recognized anywhere in a line, also inside longer words, so
  don't add free text to keyword lines.
- Check the `INPUT SUMMARY` at the top of the output: it lists the
  analyses and options as the program understood them. An analysis
  missing from it was not read.
