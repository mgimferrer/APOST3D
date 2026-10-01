# Block section # CUBE

Selects the effective orbitals written as cube files when `CUBE` is set
in `# METHOD`. The block is required with `CUBE`, but may be empty (all
defaults). Cube files are written for `EFFAO`, `UEFFAO`, `EFFAO-U`, `EOS`
and `GEOS`; what they contain, how they are named and how to view them is
described in [Visualizing orbitals](../guide/visualization.md).

| Keyword | Value | Default | Description |
| ------- | ----- | ------- | ----------- |
| `MAX_OCC` | integer *n* ≥ 0 | 1000 | Upper bound of the net occupation, *n*/1000. The default leaves out fully occupied orbitals. |
| `MIN_OCC` | integer *n* ≥ 0 | 0 | Lower bound of the net occupation, *n*/1000. |
| `SPACING` | real | 0.25 | Distance between grid points, in bohr. |
| `RADIUS_SCALE` | real | 2.0 | Padding around the fragment's atoms, in units of the covalent radius of the outermost atom on each side. Larger values give bigger boxes. |
| `NEG_EFOS` | integer *n* (optional) | 25 | `GEOS` and `EFFAO-U`: also write the paired orbitals with net occupation ≤ −*n*/1000, whatever `MAX_OCC`/`MIN_OCC` are (see [GEOS](../methods/geos.md)). |

A cube file is written for every orbital of each fragment (or atom) and
each density whose net occupation lies between `MIN_OCC`/1000 and
`MAX_OCC`/1000.

- For densities whose occupations go up to 2 (`EFFAO`, and the paired
  density of `GEOS`/`EFFAO-U`), both bounds are doubled: `MAX_OCC 700`
  then means 1.4. `NEG_EFOS` is never doubled.
- `MAX_OCC` and `MIN_OCC` must be ≥ 0, with `MIN_OCC` ≤ `MAX_OCC`, or the
  run stops. `MAX_OCC 0` with `MIN_OCC 0` selects no orbital, which is
  useful together with `NEG_EFOS`.

Example: orbitals with net occupation between 0.3 and 0.7 (per spin for
`EOS`), on a finer grid than the default:

```text
# CUBE
MAX_OCC 700
MIN_OCC 300
SPACING 0.15
#
```
