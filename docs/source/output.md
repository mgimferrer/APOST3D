# Output files

A run of `apost3d jobname` reads `jobname.fchk` and `jobname.inp` and writes
the files below into the working directory. Which ones appear depends on
the keywords used.

| File | Written when | Contents |
|---|---|---|
| main output (standard output) | always | all printed results; redirect it, e.g. `apost3d jobname > jobname.apost` |
| `jobname_<scheme>..._<i>.cube` | `CUBE` | one cube file per selected effective orbital |
| `jobname-GEOS-EFOs.fchk` | `GEOS` | GEOS effective fragment orbitals as viewable "MOs" |
| `jobname-EOS-EFOs.fchk` | `EOS` with a real-space scheme | EOS effective fragment orbitals as viewable "MOs" |
| `jobname-OSLOs.fchk` | `OSLO` | the final OSLOs |
| `jobname-OSLOs-preortho.fchk` | `OSLO` with `PRINT NON-ORTHO` | the OSLOs before orthogonalization |
| `jobname<ext>.files`, `jobname<ext>_<El><n>.int` | `DOINT` | atomic overlap matrices in the MO basis, one `.int` file per atom, for the ESI program |
| `efo_occ.dat`, `efo_coeff.dat` | `EOS`/`UEFFAO` with `LOWDIN` | effective fragment orbital occupations and coefficients as plain text |

In the `DOINT` file names, `<ext>` is `fuz` (TFVC and Becke), `hir`
(Hirshfeld), `ihi` (iterative Hirshfeld), `mul` (Mulliken) or `low`
(Löwdin), and `<El><n>` is the element symbol and atom number, e.g.
`jobnamefuz_O1.int`.

## Main output

The main output is printed in this order:

1. program header and citations;
2. `WAVEFUNCTION SUMMARY`: what was read from the `.fchk` file;
3. `INPUT SUMMARY`: how every keyword was understood; check it first when a
   result looks unexpected;
4. the integration grid and atomic definition (real-space schemes);
5. the population analysis: charges, bond orders and valences, per atom
   and, with `DOFRAGS`, per fragment;
6. the requested analyses, in the order they are run;
7. `...Normal Termination of APOST-3D...`.

The `TIMING CPU` / `TIMING WALL` lines after each step give its cost.

A run that ends without the `Normal Termination` line did not finish;
the last lines of the output usually say why. See
[Troubleshooting](troubleshooting.md).

How to read the results of a specific analysis is described on its page
(e.g. [GEOS](methods/geos.md)).

## Cube files

Written with `CUBE` in `# METHOD`; the orbitals are selected in the
[# CUBE block section](input/other-blocks.md).

**File names**: `jobname_<scheme><density>_<X><n>_<i>[beta].cube`

- `<scheme>`: `tfvc`, `becke`, `beckerho`, `hirsh`, `hirsh-it`,
  `mulliken` or `lowdin`;
- `<density>`: empty for EFFAO/EOS, `_paired`, `_unpaired` or
  `_paired_neg` for GEOS;
- `<X><n>`: `FR<n>` for fragment *n* (with `DOFRAGS`), otherwise element
  symbol and atom number;
- `<i>`: orbital number within that fragment, by decreasing net
  occupation (for `_paired_neg`: 1 is the most negative);
- `beta`: added for the beta orbitals of `EOS`/`UEFFAO`.

The second line of each file gives the orbital's gross and net occupation.

**Values**: for real-space schemes, the orbital multiplied by the
fragment's weight function, i.e. the part of the orbital that belongs to
that fragment; for Mulliken and Löwdin, the plain orbital. Values are not
renormalized. The `Normalization from cube` value printed in the main
output for each file approaches 1 (real-space schemes) when the box contains
the whole orbital; if it is clearly lower, enlarge the box with
`RADIUS_SCALE`.

## Effective orbitals in `.fchk` files

`GEOS` and real-space `EOS` runs always write the effective fragment
orbitals of all fragments into a copy of the input `.fchk`, as its
molecular orbitals, so they can be viewed in any program that reads `.fchk`
files:

- orbitals are sorted by decreasing gross occupation, and the "orbital
  energy" field holds that gross occupation;
- `EOS`: alpha orbitals in the Alpha set, beta orbitals in the Beta set
  (a single set for closed-shell systems);
- `GEOS`: paired orbitals in the Alpha set, unpaired in the Beta set;
  see [GEOS](methods/geos.md) for details, including negative-occupation
  orbitals;
- the exported orbitals are the best basis-set approximation of each
  orbital restricted to its fragment; the `FIT %` line under each
  fragment's occupations in the main output gives how close it is (see
  [GEOS](methods/geos.md));
- everything else in the file (geometry, basis set, densities) is copied
  unchanged.

`EOS` with `MULLIKEN`/`LOWDIN` does not write this file.
