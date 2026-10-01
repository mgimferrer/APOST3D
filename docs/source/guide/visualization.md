# Visualizing orbitals

The orbitals that APOST-3D builds (effective atomic and fragment orbitals,
OSLOs) can be looked at in two ways:

| | Orbital `.fchk` files | Cube files |
|---|---|---|
| Written by | `OSLO`, `GEOS`, and `EOS` with a real-space scheme, always | `EFFAO`, `UEFFAO`, `EFFAO-U`, `EOS`, `GEOS`, with `CUBE` in `# METHOD` and a [`# CUBE` block](../input/cube.md) |
| Content | All orbitals of all fragments, in one file | One file per selected orbital, values on a grid |
| View with | Any program that reads `.fchk` files (GaussView, Avogadro, IQmol, Jmol, Multiwfn, ...) | Any program that reads cube files (VMD, Avogadro, VESTA, Chimera, ...) |

The `.fchk` files are the easiest way to browse all orbitals at once;
cube files show an effective orbital exactly as the program computes it.

## Orbital `.fchk` files

| File | Written by | Orbitals |
|---|---|---|
| `<jobname>-OSLOs.fchk` | `OSLO` | The final (orthogonalized) OSLOs |
| `<jobname>-OSLOs-preortho.fchk` | `OSLO` with `PRINT NON-ORTHO` | The OSLOs before orthogonalization |
| `<jobname>-EOS-EFOs.fchk` | `EOS` with a real-space scheme | Effective fragment orbitals: alpha in the Alpha set, beta in the Beta set (only alpha for closed shells) |
| `<jobname>-GEOS-EFOs.fchk` | `GEOS` | Paired orbitals in the Alpha set, unpaired ones in the Beta set (see [GEOS](../methods/geos.md#orbitals-in-the-fchk-file)) |

Each file is a copy of the input `.fchk` whose molecular orbitals are
replaced by the analysis orbitals; the geometry and basis set are
unchanged. In the EFO files:

- the orbitals of all fragments are sorted together by decreasing gross
  occupation, and the **"orbital energy" field holds the gross
  occupation**, so the viewer's orbital list reads as an occupation list;
- fragment labels are not stored: match an orbital to its fragment
  through its occupation and the main output;
- the densities in the file are the original ones, not built from the
  EFOs.

**How faithful the EFOs in the file are.** An effective fragment orbital
is cut to its fragment by the fragment's weight function, and that shape
can't be written exactly in terms of the basis functions, which is all a
`.fchk` file can hold. The file contains its best basis-set
approximation, and the `FIT %` row printed under each fragment's
occupations in the main output gives its quality, 100 × (1 − ‖difference‖
/ ‖orbital‖) (`% recovered` for the negative GEOS orbitals). Typical
values:

| EFO net occupation | Typical `FIT %` |
|---|---|
| ≳ 0.9 | 98–99.9 |
| 0.4–0.8 | 84–95 |
| 0.05–0.1 | 60–80 |
| below ~0.01 | often ≤ 50 |

The measure is strict: 84% still reproduces about 97% of the orbital's
squared norm. The limit comes from the basis set of the original
calculation (a larger one improves it), not from the integration grid.
All numbers in the output are computed from the exact EFOs, never from
the approximation. With Mulliken or Löwdin, `EOS` writes no `.fchk` file
(`LOWDIN` and `NAO-BASIS` write the orbitals as text instead, see
[EOS](../methods/eos.md#files-written)).

## Cube files

Selected in the [`# CUBE` block](../input/cube.md) by net occupation: by
default every orbital of each fragment below full occupation is written.

**File names**: `<jobname>_<scheme><density>_<X><n>_<i>[beta].cube`

- `<scheme>`: `tfvc`, `becke`, `beckerho`, `hirsh`, `hirsh-it`,
  `mulliken` or `lowdin`;
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

**Box size.** The box encloses the fragment's atoms plus a margin
(`RADIUS_SCALE` times the covalent radius of the outermost atom on each
side). For each file the output prints `Normalization from cube`, which
approaches 1 (real-space schemes) when the box holds the whole orbital.
Small fragments, especially single H atoms, get small boxes that cut
diffuse orbitals: for the H atom of stretched LiH, the normalization is
0.29 with the default `RADIUS_SCALE 2` and 0.99 with `RADIUS_SCALE 12`.
Larger boxes make larger files; a coarser `SPACING` compensates.
