# Visualizing orbitals

The orbitals that APOST-3D builds (effective atomic and fragment orbitals,
OSLOs) can be looked at in two ways:

| | Orbital `.fchk` files | Cube files |
|---|---|---|
| Written by | `OSLO`, `GEOS` and `EOS`, always | `EFFAO`, `UEFFAO`, `EFFAO-U`, `EOS`, `GEOS`, with `CUBE` in `# METHOD` and a [`# CUBE` block](../input/cube.md) |
| Content | All orbitals of all fragments, in one file | One file per selected orbital, values on a grid |
| View with | Any program that reads `.fchk` files (GaussView, Avogadro, IQmol, Jmol, Multiwfn, ...) | Any program that reads cube files (VMD, Avogadro, VESTA, Chimera, ...) |

The `.fchk` files are the easiest way to browse all orbitals at once;
cube files show an effective orbital exactly as the program computes it.

```{admonition} APOST3Dview
:class: tip

[APOST3Dview](https://github.com/mgimferrer/APOST3Dview), the companion
viewer of APOST-3D, opens the results file `<jobname>.apost.json`
together with its `.fchk` file and shows the computed values on the
molecule; it also opens cube, `.molden` and `.xyz` files and draws
orbitals directly from the wavefunction. It is under development.
```

## Orbital `.fchk` files

| File | Written by | Orbitals |
|---|---|---|
| `<jobname>-OSLOs.fchk` | `OSLO` | The final (orthogonalized) OSLOs |
| `<jobname>-OSLOs-preortho.fchk` | `OSLO` with `PRINT_NONORTHO` | The OSLOs before orthogonalization |
| `<jobname>-EOS-EFOs.fchk` | `EOS` | Effective fragment orbitals: alpha in the Alpha set, beta in the Beta set (only alpha for closed shells) |
| `<jobname>-GEOS-EFOs.fchk` | `GEOS` | Paired orbitals in the Alpha set, unpaired ones in the Beta set (see [GEOS](../methods/geos.md#orbitals-in-the-fchk-file)) |

Each file is a copy of the input `.fchk` whose molecular orbitals are
replaced by the analysis orbitals; the geometry and basis set are
unchanged. In the OSLO files the orbitals come in the order they were
assigned, and the **"orbital energy" field holds the FOLI value** of
each (0 after the last OSLO). In the EFO files:

- the orbitals of all fragments are sorted together by decreasing gross
  occupation, and the **"orbital energy" field holds the gross
  occupation**, so the viewer's orbital list reads as an occupation list;
- fragment labels are not stored in the file: the results file
  (`<jobname>.apost.json`) maps every orbital to its fragment, and so
  does the main output through the occupations;
- the densities in the file are the original ones, not built from the
  EFOs.

**How faithful the EFOs in the file are.** With a Hilbert-space scheme
(`MULLIKEN`, `LOWDIN`, `LOWDIN-DAVIDSON`, `NAO-BASIS`) the EFOs are
combinations of basis functions and are written exactly. With a
real-space scheme, an effective fragment orbital
is cut to its fragment by the fragment's weight function, and that shape
can't be written exactly in terms of the basis functions, which is all a
`.fchk` file can hold. The file contains its best basis-set
approximation, and the `Fit (%)` row printed under each fragment's
occupations in the main output gives its quality, 100 × (1 − ‖difference‖
/ ‖orbital‖), also for the negative GEOS orbitals. Typical
values:

| EFO net occupation | Typical `Fit (%)` |
|---|---|
| ≳ 0.9 | 86–99.9 (mostly above 94) |
| 0.4–0.8 | 82–99 |
| 0.05–0.1 | 60–80 |
| below ~0.01 | often ≤ 50 |

The measure is strict: 84% still reproduces about 97% of the orbital's
squared norm. The limit comes from the basis set of the original
calculation (a larger one improves it), not from the integration grid.
All numbers in the output are computed from the exact EFOs, never from
the approximation.

**What the surfaces show.** An effective fragment orbital is not zero
outside its fragment: the fragment's weight function goes from 1 to 0
gradually, so the orbital keeps a small part on the neighboring atoms
(a few percent of its squared norm). A surface that reaches onto a
neighbor is therefore part of the orbital itself. The approximation in
the `.fchk` file follows it closely when `Fit (%)` is about 80 or more:
it puts at most about 3% more of its norm outside the fragment, and at an
isovalue of 0.075 it shows the lobes of the exact cube file and no others
(for a diffuse orbital its surface can be somewhat smaller). It shows
lobes that the exact orbital doesn't have:

- for weakly occupied orbitals: small extra lobes with `Fit (%)` between
  about 60 and 80, clear ones below about 60;
- at small isovalues (0.03 and below), where the small differences
  between the two become visible.

To see the exact shape of such an orbital, write its cube file (`CUBE` in
`# METHOD`, see [Cube files](#cube-files)): it holds the exact orbital
times the fragment's weight, on the same scale as the `.fchk` orbital, so
the same isovalue can be used for both.

## Cube files

Selected in the [`# CUBE` block](../input/cube.md) by net occupation: by
default every orbital of each fragment below full occupation is written.

**File names**: `<jobname>_<scheme><density>_<X><n>_<i>[beta].cube`

- `<scheme>`: `tfvc`, `becke`, `beckerho`, `hirsh`, `hirsh-it`,
  `mulliken` or `lowdin` (also for `LOWDIN-DAVIDSON` and `NAO-BASIS`);
- `<density>`: empty for EFFAO/EOS; `_paired`, `_unpaired` or
  `_paired_neg` for GEOS and EFFAO-U;
- `<X><n>`: `FR<n>` for fragment *n* with `DOFRAGS`, otherwise the element
  symbol and atom number (e.g. `O1`);
- `<i>`: orbital number within the fragment, by decreasing net occupation
  (for `_paired_neg`, 1 is the most negative);
- `beta`: added for the beta orbitals of `EOS`/`UEFFAO`.

Example: `LiH-32-FCI_tfvc_paired_neg_FR2_1.cube`. The second line of each
file gives the orbital's gross and net occupation.

**Values**: with a real-space scheme, the orbital multiplied by the
fragment's weight function, i.e. the part of the orbital that belongs to
the fragment; with Mulliken or Löwdin, the plain orbital. Cube files and
`.fchk` orbitals are on the same scale (neither is renormalized), so the
same isovalue can be used for both.

**Box size.** The box encloses the fragment's atoms plus a margin of
`PADDING` bohr (5 by default) on every side. For each fragment, the
output lists the files written with each orbital's net occupation and
its `Norm in cube`, the integral of the squared cube values. Each
effective orbital has norm 1 over all space (in real space, the orbital
times the fragment's weight), so `Norm in cube` is the share of the
orbital inside the box. The output warns when it is below 0.9 or above
1.1:

- below 0.9, the box cuts the orbital: raise `PADDING`. Diffuse
  orbitals need the most room: for the H atom of stretched LiH, the
  norm is 0.93 with the default margin and 0.996 with `PADDING 8`. At
  the default, the orbitals of the test systems stay below 0.03 on the
  box faces, so their surfaces at the usual isovalues are not cut open;
- far from 1, either way, the grid is too coarse for a compact orbital,
  typically a core orbital sampled near its nucleus (3.29 for the F 1s
  of CH3F at the default spacing): lower `SPACING`. The norm, a sum over
  the grid points, suffers much more than the surfaces.

Larger boxes and finer grids make larger files: with the defaults, the
15 cube files of CH3F take 19 MB. The files follow the
standard cube layout (positive atom count, no orbital-index line).
