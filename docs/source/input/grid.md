# Block section # GRID

Integration grid for the two-electron terms of [ENPART](../methods/enpart.md).
The block is only read if `MOD-GRIDTWOEL` is given in `# ENPART`; without
it the defaults below are used, and a `# GRID` block present in the input
is ignored (with a warning in the output). The grid of every other
analysis is set on the command line instead (see
[Running a calculation](../guide/running.md)).

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `RADIAL` | integer | 150 | Radial points per atom, 1 to 500. |
| `ANGULAR` | integer | 590 | Angular (Lebedev-Laikov) points per radial shell. Must be one of 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974; any other value stops the run with a message. |
| `rr00` | real | 0.5 | Radial scaling: the distance (in bohr) that contains half of the radial points (*r*ₘ in Eq. 25 of [Becke's scheme](https://doi.org/10.1063/1.454033)). |
| `phb1` | real | 0.169 | Rotation angle (radians) of the grid of the second electron, for the one-center terms. |
| `phb2` | real | 0.170 | Second rotation angle (radians), used by the zero-error strategy. |

```{admonition} The rotation angles depend on the grid
:class: warning

The default `phb1`/`phb2` are calibrated for the 150 × 590 grid. With
another grid, set the angles that belong to it: for 40 × 146, `phb1 0.162`
and `phb2 0.182`.
```

Example, a cheaper two-electron grid for ENPART:

```text
# ENPART
B3LYP
MOD-GRIDTWOEL
#
# GRID
RADIAL 40
ANGULAR 146
phb1 0.162
phb2 0.182
#
```
