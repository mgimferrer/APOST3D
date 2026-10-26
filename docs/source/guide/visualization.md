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
| `<jobname>-OSLOs.fchk` | `OSLO` | The final (orthogonalized) OSLOs: alpha in the Alpha set, beta in the Beta set (only alpha for closed shells) |
| `<jobname>-OSLOs-preortho.fchk` | `OSLO` with `PRINT_NONORTHO` | The OSLOs before orthogonalization |
| `<jobname>-EOS-EFOs.fchk` | `EOS` | Effective fragment orbitals: alpha in the Alpha set, beta in the Beta set (only alpha for closed shells) |
| `<jobname>-GEOS-EFOs.fchk` | `GEOS` | Paired orbitals in the Alpha set, unpaired ones in the Beta set (see [GEOS](../methods/geos.md#orbitals-in-the-fchk-file)) |
| `<jobname>-EOS-EFOs-unweighted.fchk`, `<jobname>-GEOS-EFOs-unweighted.fchk` | `EOS`, `GEOS` with `BECKE`, `BECKE-RHO` or `TFVC` atoms | The same EFOs before the multiplication by the fragment's weight (see [below](#the-unweighted-files)) |

Each file is a copy of the input `.fchk` whose molecular orbitals are
replaced by the analysis orbitals; the geometry and basis set are
unchanged (the OSLO files also carry the density built from the
OSLOs). The beta orbitals of an open shell go in the Beta set, also when
the input `.fchk` has a single set of orbitals (ROHF/ROKS: the Beta
blocks are added after the Alpha ones). In the OSLO files the orbitals come in the order they were
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

```{admonition} The shape of an effective fragment orbital
:class: important

With a real-space scheme, judge the shape of an EFO from its cube file
(`CUBE` in `# METHOD`) or, with `BECKE`, `BECKE-RHO` or `TFVC` atoms, in
APOST3Dview: both show the exact orbital. The orbitals of
`<jobname>-EOS-EFOs.fchk` and `<jobname>-GEOS-EFOs.fchk`, which any
program can show, are their least-squares projections onto the basis
set: close to the EFO on its own fragment, with small spurious tails on
the other fragments, and only rough for weakly occupied EFOs, whose fit
is poor. They should not be used to judge an EFO's shape.
```

**Which file holds what.** For an EFO of a real-space scheme, with w_F
the weight function of its fragment:

| File | Holds | Use it for |
|---|---|---|
| `<jobname>-EOS-EFOs.fchk`, `-GEOS-EFOs.fchk` | The projection of the EFO onto the basis set: an approximation, of quality `Fit (%)` | Programs that show the orbitals stored in a `.fchk` file; browsing all EFOs at once |
| `<jobname>-EOS-EFOs-unweighted.fchk`, `-GEOS-EFOs-unweighted.fchk` (`BECKE`, `BECKE-RHO`, `TFVC`) | The basis-set vector φ, exact; the EFO is w_F × φ | APOST3Dview, which multiplies by w_F; programs that apply w_F themselves |
| Cube files (`CUBE`) | w_F × φ on a grid: the EFO, exact | The shape of the EFO, in any program that reads cube files |

With a Hilbert-space scheme the EFOs of the first file are exact, and no
unweighted file is written.

**How faithful the EFOs in the file are.** With a Hilbert-space scheme
(`MULLIKEN`, `LOWDIN`, `LOWDIN-DAVIDSON`, `NAO-BASIS`) the EFOs are
combinations of basis functions and are written exactly. With a
real-space scheme, an effective fragment orbital
is cut to its fragment by the fragment's weight function, and that shape
can't be written exactly in terms of the basis functions, which is all a
`.fchk` file can hold. The file contains its best basis-set
approximation, and the `Fit (%)` row printed under each fragment's
occupations in the main output gives its quality, 100 × (1 − ‖difference‖
/ ‖orbital‖), also for the negative GEOS orbitals. `Fit (%)` concerns
only this file: the cube files and the unweighted file are exact. Typical
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
the `.fchk` file is close to it when `Fit (%)` is about 80 or more: it
puts at most about 3% more of its norm outside the fragment, and at an
isovalue of 0.075 it shows the lobes of the exact cube file and no others
(for a diffuse orbital its surface can be somewhat smaller). Even then,
its values differ from the exact ones by a few percent of the largest
value, so the cube file or APOST3Dview remain the reference for the
shape. The `.fchk` approximation shows lobes that the exact orbital
doesn't have:

- for weakly occupied orbitals: small extra lobes with `Fit (%)` between
  about 60 and 80, clear ones below about 60;
- at small isovalues (0.03 and below), where the small differences
  between the two become visible.

Compared at every grid point with the cube files of the TFVC test
calculations `CH3F-EOS-CUBE` (15 EFOs) and `LiH-32-FCI` (the negative
paired EFO), as a fraction of the cube's largest value:

| Orbital | Largest difference from the cube |
|---|---|
| w_F × φ, from the unweighted file | ≤ 2.9 × 10⁻⁵ (RMS ≤ 1.0 × 10⁻⁶): the cube's 5 printed digits |
| `.fchk` orbital, `Fit (%)` 98–100 | 0.0001–0.003 |
| `.fchk` orbital, `Fit (%)` 86–95 (well-fitted occupied EFOs) | 0.02–0.14 |
| `.fchk` orbital, `Fit (%)` 34–60 (weakly occupied EFOs) | up to 0.62 |
| `.fchk` orbital multiplied by w_F once more | up to 0.67 |

To see the exact shape of an orbital, write its cube file (`CUBE` in
`# METHOD`, see [Cube files](#cube-files)) or, with `BECKE`, `BECKE-RHO`
or `TFVC` atoms, open the results file in APOST3Dview (see
[below](#the-unweighted-files)). The cube file holds the exact orbital
times the fragment's weight, on the same scale as the `.fchk` orbital, so
the same isovalue can be used for both.

### The unweighted files

With `BECKE`, `BECKE-RHO` or `TFVC` atoms, EOS and GEOS also write
`<jobname>-EOS-EFOs-unweighted.fchk` (or `-GEOS-EFOs-unweighted.fchk`),
with the same orbitals, in the same order, as the file above, but taken
**before** the multiplication by the fragment's weight w_F. Each of
these orbitals φ is a combination of basis functions, written exactly,
and the effective fragment orbital is w_F × φ: the cube file of the same
EFO holds exactly that product. On its own, φ is not the EFO (it is not
cut to its fragment and can be large elsewhere), so the file is meant
for programs that multiply each orbital by w_F themselves; the results
file gives everything needed to rebuild w_F
([`atomic_weights`](../developer/json-output.md#atomic_weights)). The
output names both files at the end of EOS and GEOS.

The unweighted file is not written with `HIRSH` and `HIRSH-IT`, whose
weights the results file does not give (the output says so), nor with
Hilbert-space schemes, whose EFOs are already exact in the first file.

**In APOST3Dview.** With `BECKE`, `BECKE-RHO` or `TFVC` atoms,
APOST3Dview reads the unweighted file and multiplies each orbital by its
fragment's weight function, rebuilt from the results file, which
reproduces the cube file of the orbital. With `HIRSH` and `HIRSH-IT` it
shows the orbitals of the first file as stored, as other programs do,
extra lobes included.

## Cube files

Selected in the [`# CUBE` block](../input/cube.md) by net occupation: by
default every orbital of each fragment is written, core ones included
(see `MAX_OCC`).

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
