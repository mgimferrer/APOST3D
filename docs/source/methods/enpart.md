# ENPART - Energy partitioning into one- and two-center terms (IQA)

`ENPART` decomposes the molecular energy into **one-center** (atomic
self-energy) and **two-center** (interatomic interaction) terms, in the
spirit of the Interacting Quantum Atoms (IQA) approach, using fuzzy
real-space atoms (e.g. `TFVC`):

*E* = Σ<sub>A</sub> *E*<sub>self</sub>(A) + Σ<sub>A&lt;B</sub> *E*<sub>int</sub>(A,B)

Each interaction energy collects the nuclear repulsion, the two
electron-nuclear attractions, the classical (Coulomb) electron repulsion
and the exchange-correlation (xc) energy between the two atoms. It works for
Hartree-Fock, Kohn-Sham DFT (LDA, GGA and global-hybrid GGA functionals)
and correlated (CASSCF/CISD) wavefunctions.

```{admonition} Cite
:class: note

Hartree-Fock: P. Salvador, M. Duran and I. Mayer, *J. Chem. Phys.*,
**2001**, 115, 1153-1157; P. Salvador and I. Mayer, *J. Chem. Phys.*,
**2004**, 120, 5046-5052. KS-DFT (bond order density approach and
zero-error strategy): P. Salvador and I. Mayer, *J. Chem. Phys.*,
**2007**, 126, 234113; M. Gimferrer and P. Salvador, *J. Chem. Phys.*,
**2023**, 158, 234105. See [Citations](../citations.md).
```

## How it works

1. **One-electron terms.** Kinetic energy (one-center only),
   electron-nuclear attraction and nuclear repulsion are integrated on the
   main (one-electron) grid.
2. **Coulomb term.** The classical electron-electron repulsion is a 6-D
   integral, computed numerically. For the one-center terms the second
   electron's coordinates use a copy of the grid rotated by `phb1`, so
   that the two electrons never sit on the same points.
3. **Exchange-correlation term**, depending on the wavefunction:
   - **Hartree-Fock**: the exact exchange, integrated like the Coulomb
     term.
   - **KS-DFT**: the one-center terms come from the functional's energy
     density integrated within each atom ("exact" one-center terms). The
     interatomic terms come from evaluating the functional on the **bond
     order density** (BODEN) of each atom pair. The one-center terms are
     then rearranged so that all terms add up exactly to the KS-DFT
     exchange-correlation energy. For a global hybrid, the exact-exchange
     fraction is computed as in Hartree-Fock and added in the
     two-electron part.
   - **CASSCF/CISD**: the xc density built from the 1- and 2-RDMs (read
     from the `# DM` block) is diagonalized, and the resulting terms are
     integrated numerically. With `CORRELATION` the exchange and
     correlation parts are given separately.
4. **Pairs with a small bond order (`THREBOD`).** For pairs of atoms
   whose bond order is below the `THREBOD` threshold, the integration
   (6-D for the exchange, BODEN for KS-DFT) is skipped, and the pair's
   **whole** exchange-correlation term is estimated by a multipolar
   expansion of the xc density, from charge-charge up to
   quadrupole-quadrupole terms (E. Francisco *et al.*, *J. Comput.
   Chem.*, **2017**, 38, 816-829). For Hartree-Fock and KS-DFT the xc
   density is that of the (Kohn-Sham) determinant; for CASSCF/CISD, the
   one built from the RDMs. For a hybrid functional, the (1 − *a*) share
   appears in the DFT part and the *a* share, *a* being the
   exact-exchange fraction, with the HF-type exchange. The one-center
   terms are adjusted so that the total exchange-correlation energy does
   not change.

   The default threshold, a bond order of 0.005, was chosen for a good
   accuracy/cost ratio; `THREBOD -1` computes every pair.
5. **Zero-error strategy.** The 6-D integrations carry a numerical error
   of the order of 1 kcal/mol even with good grids. When the `.fchk` file
   contains the reference electron-electron energy (see below), the
   two-electron energy is compared with it. If the error exceeds
   `TWOELTOLER` (0.25 kcal/mol by default), the one-center two-electron
   terms are recomputed with a second rotation (`phb2`), and the two
   results are interpolated so that the error vanishes. Interatomic terms
   are not modified.

## Requirements and input

- A **real-space** atomic definition (`TFVC` recommended). ENPART stops
  with an error for Hilbert-space schemes (`MULLIKEN`, `LOWDIN`, ...).
- A `.fchk` file of the wavefunction. For correlated wavefunctions, the
  1- and 2-RDMs in a [`# DM` block](../input/dm.md).
- **Reference energies**, recommended: the kinetic, electron-nuclear and
  electron-electron energies of the calculation, appended to the `.fchk`
  file, enable the integration-error checks and the zero-error strategy.
  For Gaussian they are added with `utils/get_energy_g16` (or
  `get_energy` for Gaussian 09); `.fchk` files written from pySCF with
  `utils/apost3d.py` already contain them. See
  [Preparing the wavefunction](../guide/wavefunctions.md#gaussian).
- Only the **electronic energy** is decomposed: wavefunctions from
  calculations with an empirical dispersion correction (e.g. Grimme's
  GD3) cannot be decomposed directly.

```text
# METHOD
TFVC
ENPART
#
# ENPART
LIBRARY
EX_FUNCTIONAL 106
EC_FUNCTIONAL 131
#
```

The full keyword list is in [Block section # ENPART](../input/enpart.md).

### Choosing the functional

The functional must be the one used to compute the wavefunction. Give it
with a predefined keyword (`SVWN`, `BLYP`, `BP86`, `PBE`, `B3LYP`,
`B3PW91`, `B3P86`, `PBE0`, `BHANDHLYP`, ...), or with `LIBRARY` and its
[libxc](https://libxc.gitlab.io/functionals/) identifiers: either
`EXC_FUNCTIONAL` alone (a combined exchange-correlation functional), or
`EX_FUNCTIONAL` and/or `EC_FUNCTIONAL` (separate exchange and correlation
parts), never both kinds together. Each predefined keyword was checked to
reproduce the Gaussian 16 functional of the same name; the full list,
with libxc ids, is in [Block section # ENPART](../input/enpart.md).

| Functional family | Supported |
|---|---|
| LDA, GGA | yes |
| Global hybrid GGA (B3LYP, PBE0, ...) | yes |
| Meta-GGA and hybrid meta-GGA (TPSS, M06-2X, ...) | not yet (implementation in progress) |
| Range-separated hybrids (CAM-B3LYP, ωB97X, LC-ωPBE, HSE, ...) | no |
| Functionals with VV10 nonlocal correlation | no |

Unsupported functionals stop the run with a message before any energy is
computed.

```{admonition} Same name, different functional
:class: warning

Programs do not agree on what some names mean: ORCA's and Turbomole's
B3LYP is not Gaussian's (libxc 475 instead of 402, about 23 kcal/mol
apart for water), and libxc's own "B3P86" (403) is not Gaussian's either
(the `B3P86` keyword uses 315). A mismatched functional is not always
obvious in the output, because the zero-error strategy can absorb it into
the one-center terms (see *Checking the result* below). Details in
[Block section # ENPART](../input/enpart.md).
```

## Reading the output

The excerpts below come from water, RKS SVWN/cc-pVDZ from Gaussian 16,
with the `SVWN` keyword, `THREBOD -1` and a 40 × 146 two-electron grid
(test `H2O-SVWN`).

**Functional information.** One box per libxc component, with its name,
references, type and family.

**One-electron part.** One matrix per term (`ELECTRON-NUCLEAR
ATTRACTION`, `KINETIC ENERGY`, `NUCLEAR-NUCLEAR REPULSION`). With
reference energies in the `.fchk`, each total is followed by its
integration error:

```text
  Electron-nuclei energy:   -198.9264417
  Integration error (kcal/mol):     0.01
```

**KS-DFT exchange-correlation.** The BODEN of each pair (`multipolar`
for pairs below `THREBOD`, followed by a note on how they are treated),
then the `DIATOMIC PURE KS-DFT XC TERMS (BODEN)`, the `PURE KS-DFT XC
ONE-CENTER TERMS (EXACT)` and, after rearranging, the final matrix:

```text
    FINAL PURE KS-DFT EXCHANGE-CORRELATION ENERGY COMPONENTS
              1  O        2  H        3  H
    1  O    -8.498405   -0.181686   -0.181686
    2  H    -0.181686   -0.059858   -0.000594
    3  H    -0.181686   -0.000594   -0.059858
  Sum of pure KS-DFT exchange-correlation energy:     -8.9820876
```

For a hybrid functional, this matrix is only the DFT part, and a note
says so. The exact-exchange part is added in the two-electron section,
which prints `HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS` (the full HF-type
exchange, before scaling by the exact-exchange fraction, also noted
there) and then the complete `HYBRID KS-DFT XC TERMS`.

**Two-electron part.** The Coulomb matrix, the electron-electron energy
and its integration error, followed by the zero-error strategy:

```text
  KS-DFT electron-electron energy (au):     37.8661299
  Integration error (kcal/mol):    -3.09

  Max error accepted on the two-electron part (kcal/mol):     0.25

  Zero-error strategy applied: the INTERPOLATED tables below replace
  the ones above (only their one-center terms change).
  Rotating for angles:   0.000000   0.182000
  New error after rotation:   -10.41
  Same-sign error on both grids: extrapolating, not interpolating
  Damping parameter       :      1.4227849
```

When the zero-error strategy is applied, a note says so, and the
`INTERPOLATED ...` matrices that follow replace the first ones (they
differ only on the diagonal). Without reference energies, a note says
that the error is neither estimated nor corrected.

**Final decomposition.**

```text
    FUZZY ATOMS KS-DFT ENERGY COMPONENTS
              1  O        2  H        3  H
    1  O   -74.525615   -0.529236   -0.529236
    2  H    -0.529236   -0.302629    0.138786
    3  H    -0.529236    0.138786   -0.302629
  Total KS-DFT energy      :    -76.0505604
  Total energy in Fchk file:    -76.0505886
  Integration error (au):      0.0000282
  Integration error (kcal/mol):     0.02
```

- The diagonal elements are the atomic self-energies *E*<sub>self</sub>(A)
  and the off-diagonal ones the interaction energies
  *E*<sub>int</sub>(A,B). The total is the sum of the upper triangle
  including the diagonal.
- The same convention holds for every matrix printed before: each
  off-diagonal element is the full A–B term. The exchange-correlation
  part of an interaction is the off-diagonal element of the **last**
  exchange-correlation matrix printed (for hybrids `HYBRID KS-DFT XC
  TERMS`, or its `INTERPOLATED` version); its
  classical part is the rest, *E*<sub>int</sub>(A,B) − *V*<sub>xc</sub>(A,B).

Hartree-Fock and CASSCF/CISD outputs follow the same pattern, with
`HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS` or `POST-HARTREE-FOCK-TYPE
EXCHANGE-CORRELATION ENERGY TERMS`, and a final `FUZZY ATOMS Hartree-Fock
ENERGY COMPONENTS` or `FUZZY ATOMS Post-Hartree-Fock ENERGY COMPONENTS`
matrix. With `DOFRAGS`, every matrix is also condensed to fragments
(`FRAGMENT ANALYSIS: ...`).

## Checking the result

- **The final integration error is not a quality check by itself** when
  the zero-error strategy has been applied: it is close to zero by
  construction, even if the functional does not match the wavefunction.
- The informative numbers are the one-electron integration errors, the
  two-electron error **before** the zero-error strategy (a few kcal/mol
  at most with sensible grids), and the damping parameter (about 1.2 to
  1.4 in the tests, with the 40 × 146 two-electron grid). A much larger error or damping parameter
  points to a functional that does not match the wavefunction.
- Without reference energies, the final integration error is the raw
  numerical error and depends strongly on the two-electron grid,
  especially for heavier atoms.

## Grids and cost

- **One-electron grid** (also used for the DFT exchange-correlation):
  150 radial × 590 angular points per atom by default for ENPART.
- **Two-electron grid**: set in the `# GRID` block, which is only read
  if `MOD-GRIDTWOEL` is given in `# ENPART` (see
  [Block section # GRID](../input/grid.md)). Defaults: 150/590
  with `phb1 0.169`, `phb2 0.170`. The rotation angles are calibrated for
  the grid: for 40/146 use `phb1 0.162`, `phb2 0.182`.
- The 6-D two-electron integrations dominate the cost, which grows with
  the square of the number of points per atom: 40/146 is much cheaper
  than 150/590 and, with the zero-error strategy, often enough.
  `THREBOD` saves time by skipping weakly bonded pairs.
- The DFT part stores the orbital gradients on the one-electron grid:
  about 2 MB per occupied orbital and atom (twice that for unrestricted
  wavefunctions).
