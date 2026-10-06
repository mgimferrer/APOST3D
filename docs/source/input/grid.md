# Block section # GRID

Optional block for the numerical integration grids. Without it the
defaults below are used: give it only to change them, and only the
keywords you want to change.

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `RADIAL` | integer | 40 (150 with ENPART) | Radial points per atom of the **one-electron grid**, used by every real-space analysis (populations, EFFAO, EOS, GEOS, OSLO, SPIN) and by the one-electron terms of ENPART. 1 to 500. |
| `ANGULAR` | integer | 146 (590 with ENPART) | Angular (Lebedev-Laikov) points per radial shell of the one-electron grid. Must be one of 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974; any other value stops the run with a message. |
| `RADIAL_2E` | integer | 150 | Radial points per atom of the **two-electron grid** of [ENPART](../methods/enpart.md). |
| `ANGULAR_2E` | integer | 590 | Angular points of the two-electron grid; the same list as `ANGULAR`. |

The defaults of the one-electron grid are enough for populations, bond
orders and effective orbitals with `TFVC` (see
[Running a calculation](../guide/running.md#integration-grid)). The
two-electron grid sets both the accuracy and the cost of ENPART, which
grows with the square of its points per atom: the default 150 × 590 gives
publication-quality two-electron terms; 40 × 146 is much cheaper, for
tests and larger systems (see [ENPART](../methods/enpart.md#grids-and-cost)).

Example, a finer one-electron grid for a population analysis:

```text
# GRID
RADIAL 70
ANGULAR 434
#
```

Example, a cheaper two-electron grid for ENPART:

```text
# GRID
RADIAL_2E 40
ANGULAR_2E 146
#
```

The INPUT SUMMARY at the top of the output shows the grids taken from
`# GRID`, and for ENPART the two-electron grid and its rotation angles
with where they come from.

## Advanced: rotation angles and radial scaling

ENPART integrates the one-center two-electron terms with a copy of the
two-electron grid rotated by a small angle (so that the two electrons
never sit on the same points), and the [zero-error
strategy](../methods/enpart.md#how-it-works) uses a second rotation. The
angles have to suit the angular grid, and the program chooses them:

| `ANGULAR_2E` | Angles (radians) |
|---|---|
| 146 | 0.162, 0.182 |
| 590 | 0.169, 0.170 |

No calibrated angles exist yet for the other angular grids: the 590 pair
is then used, with a warning in the output. Expert users can set them:

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `ROTATION_2E` | two reals | from the table | The two rotation angles in radians: the first for the one-center terms, the second for the zero-error strategy (e.g. `ROTATION_2E 0.162 0.182`). The zero-error interpolation depends on them; change them only with tests on your system. |
| `R0_2E` | real | 0.5 | Radial scaling of the two-electron grid: the distance (in bohr) that contains half of the radial points (*r*ₘ in Eq. 25 of [Becke's scheme](https://doi.org/10.1063/1.454033)). |
