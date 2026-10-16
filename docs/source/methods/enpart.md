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
   electron's coordinates use a copy of the grid rotated by a small
   angle, so that the two electrons never sit on the same points.
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
   terms are recomputed with a second rotation of the grid, and the two
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
- **Pseudopotentials** (ECPs) are supported. The reference
  electron-nuclear energy then includes the ECP energy, and the ECP
  matrix appended to the `.fchk` file (by the same utilities) is split
  among the atoms, Mulliken-like, and added to the one-atom
  electron-nuclear terms. With Gaussian the matrix carries a small
  rounding error, see
  [get_energy_g16](../tools/utilities.md#get_energy_g16-and-get_energy).
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
(test `H2O-SVWN`), unless said otherwise.

**Section header.** The wavefunction, the atomic definition and how to
read the matrices:

% source: tests/reference/H2O-SVWN.apost
```text
  Wavefunction                            : KS-DFT, restricted
  Atomic definition                       : TFVC (real space)

  Each matrix below holds the one-center (self) terms on the diagonal and
  the whole A-B pair term off it; its total adds the diagonal and one
  triangle. Energies in au.
```

Every matrix follows this convention: the diagonal element of atom A is
its one-center term, and each off-diagonal element is the **whole** A–B
term (not half of it), so the total is the sum of the upper triangle
including the diagonal.

**Functional** (KS-DFT). `DENSITY FUNCTIONAL` lists each libxc component
with its id, kind, family, exact-exchange fraction (hybrids) and
references.

**One-electron terms.** `ELECTRON-NUCLEAR ATTRACTION` and
`NUCLEAR-NUCLEAR REPULSION` are matrices; `KINETIC ENERGY` is a column
per atom, since the kinetic energy has no pair terms. With reference
energies in the `.fchk`, each total is followed by its integration error:

% source: tests/reference/H2O-SVWN.apost
```text
    ELECTRON-NUCLEAR ATTRACTION
                 1  O        2  H        3  H
    1  O  -184.537618   -6.402883   -6.402883
    2  H    -6.402883   -0.662774   -0.257508
    3  H    -6.402883   -0.257508   -0.662774
  Electron-nuclear attraction (au)        :    -198.9264417
  Integration error (kcal/mol)            :           0.01
```

**DFT exchange-correlation** (KS-DFT). `DFT EXCHANGE-CORRELATION (BODEN)`
counts the atom pairs computed from the bond order density and those
below `THREBOD`, taken from the multipolar expansion (with a note on how
they are treated). `DFT EXCHANGE-CORRELATION BY ATOM` gives, for each
atom, the functional integrated within the atom and its self term (the
same minus half of each of its pair terms), and `DFT
EXCHANGE-CORRELATION ENERGY` is the final matrix:

% source: tests/reference/H2O-SVWN.apost
```text
    DFT EXCHANGE-CORRELATION ENERGY
  The whole exchange-correlation energy (pure functional).
                 1  O        2  H        3  H
    1  O    -8.498405   -0.181686   -0.181686
    2  H    -0.181686   -0.059858   -0.000594
    3  H    -0.181686   -0.000594   -0.059858
  DFT exchange-correlation (au)           :      -8.9820876
```

For a pure functional this is the whole exchange-correlation energy; for
a hybrid one it is only the DFT part, as the line under the title says.

**Two-electron integration error.** With reference energies, the
electron-electron energy integrated on the first grid is compared with the
reference, and the zero-error strategy is applied when the error exceeds
`TWOELTOLER`:

% source: tests/reference/H2O-SVWN.apost
```text
    TWO-ELECTRON INTEGRATION ERROR
  Electron-electron energy (au)           :      37.8661299
  Integration error (kcal/mol)            :          -3.09
  Tolerance (kcal/mol)                    :           0.25
  Zero-error correction                   : applied to the one-center terms
  Error on the rotated grid (kcal/mol)    :         -10.41
  Same-sign error on both grids: extrapolating, not interpolating
  Damping parameter                       :       1.4227849
  Corrected electron-electron energy (au) :      37.8710578
  Corrected integration error (kcal/mol)  :           0.00
```

Below the tolerance the section says the correction is not needed;
without reference energies it says that the error is neither estimated
nor corrected. Either way, **the matrices printed after this section hold
the final values**: when the correction is applied, their one-center
Coulomb and exchange terms are the corrected ones (pair terms do not
change).

**Two-electron terms.** `COULOMB ENERGY`, then:

- Hartree-Fock: `EXCHANGE ENERGY (HARTREE-FOCK)`, the exchange matrix.
- Hybrid functional: `EXACT-EXCHANGE SHARE` and `EXCHANGE-CORRELATION
  ENERGY`, the complete exchange-correlation terms (DFT part +
  exact-exchange share).
- Pure functional: nothing more, its exchange-correlation is the DFT part
  printed before.

```{admonition} Hybrid functionals: the exact exchange is printed already scaled
:class: important

`EXACT-EXCHANGE SHARE` is the HF-type exchange **multiplied by the
exact-exchange fraction** *a* of the functional (0.20 for B3LYP), i.e. the
part of the exchange that enters the energy. The matrices add up directly:

DFT EXCHANGE-CORRELATION ENERGY + EXACT-EXCHANGE SHARE =
EXCHANGE-CORRELATION ENERGY

The full, unscaled HF-type exchange of the Kohn-Sham orbitals is printed
only with `VERBOSE 1`, as `HF-TYPE EXCHANGE, UNSCALED` (see below).
```

For water, triplet UKS B3LYP (test `H2O-T-B3LYP`):

% source: tests/reference/H2O-T-B3LYP.apost
```text
    EXACT-EXCHANGE SHARE
  0.200 x the HF-type exchange: the exact-exchange part of the
  functional, already scaled, as it enters the energy.
                 1  O        2  H        3  H
    1  O    -1.632116   -0.035618   -0.035618
    2  H    -0.035618   -0.022472   -0.002558
    3  H    -0.035618   -0.002558   -0.022472
  Exact-exchange share (au)               :      -1.7508548
```

**IQA energy terms.** The final decomposition, titled with the atomic
definition:

% source: tests/reference/H2O-SVWN.apost
```text
    IQA ENERGY TERMS (TFVC)
  Self-energies on the diagonal, interaction energies off it.
                 1  O        2  H        3  H
    1  O   -74.525615   -0.529236   -0.529236
    2  H    -0.529236   -0.302629    0.138786
    3  H    -0.529236    0.138786   -0.302629
  Total energy (au)                       :     -76.0505604
  Energy in the .fchk file (au)           :     -76.0505886
  Integration error (au)                  :       0.0000282
  Integration error (kcal/mol)            :           0.02
```

The diagonal elements are the atomic self-energies *E*<sub>self</sub>(A)
and the off-diagonal ones the interaction energies *E*<sub>int</sub>(A,B).

**Atomic energies and classical interactions.** Built from the final
matrices:

- `ATOMIC ENERGIES`: each atom's self-energy and its additive energy,
  *E*<sub>A</sub> = *E*<sub>self</sub>(A) + ½ Σ<sub>B</sub>
  *E*<sub>int</sub>(A,B); the additive energies add up to the total
  energy.
- `CLASSICAL INTERACTION ENERGIES`: the electron-nuclear, nuclear-nuclear
  and Coulomb part of each interaction energy. The rest of it,
  *E*<sub>int</sub>(A,B) − *V*<sub>cl</sub>(A,B), is the
  exchange-correlation part: the off-diagonal element of `EXCHANGE
  ENERGY (HARTREE-FOCK)`, `EXCHANGE-CORRELATION ENERGY` (hybrids) or `DFT
  EXCHANGE-CORRELATION ENERGY` (pure functionals).

% source: tests/reference/H2O-SVWN.apost
```text
    ATOMIC ENERGIES
  Additive energy = self-energy + half of each interaction energy of the atom;
  the additive energies add up to the total energy.
  Atom    Self-energy    Additive
    1  O   -74.525615  -75.054851
    2  H    -0.302629   -0.497855
    3  H    -0.302629   -0.497855
  Total    -75.130874  -76.050560
```

**Fragments.** With `DOFRAGS`, every matrix and per-atom table is followed
by its fragment version (`FRAGMENT ...`): the diagonal of a fragment
matrix is the fragment's one-center term (those of its atoms plus the
pair terms between them), the off-diagonal elements the terms between
fragments.

### External electric field

When the `.fchk` file comes from a calculation in a static electric field
(Gaussian's `Field` keyword; APOST-3D reads it from the `External
E-field` block), the wavefunction summary shows the field, and ENPART
adds an `EXTERNAL ELECTRIC FIELD` section after the one-electron terms.
This works for Hartree-Fock and KS-DFT wavefunctions; for a correlated
(CASSCF/CISD) wavefunction in a field, ENPART stops with a message.
For water, RHF/cc-pVDZ in a field of 0.005 au along *z* (test
`H2O-RHF-FIELD`, which also uses `POLAR` and therefore prints 13
decimals):

% source: tests/reference/H2O-RHF-FIELD.apost
```text
    EXTERNAL ELECTRIC FIELD
  Field x, y, z (au)                      : 0.000000  0.000000  0.005000
  Electronic dipole term (au)             :      -0.0041737862763
  Nuclear dipole term (au)                :      -0.0000000000200
  Field interaction energy (au)           :      -0.0041737862963
```

- The field interaction energy follows Gaussian's convention,
  *E* = +**μ**·**F** (the energy goes down when the dipole points along
  the field), and equals the molecular dipole times the field.
- It is **not** part of the IQA terms: the IQA matrices decompose the
  energy without it, and the energy check of `IQA ENERGY TERMS` adds it
  to their total before comparing with the energy of the `.fchk` (which
  includes it). The integration errors are therefore those of a
  field-free calculation.
- `FIELD INTERACTION PER ATOM` splits it into an intrinsic dipole term
  per atom (origin independent) and a charge-transfer term (origin
  dependent); the two columns add up to the field interaction energy.

### More detail: VERBOSE 1

`VERBOSE 1` in `# METHOD` also prints the intermediate matrices, for
checks and development (the results do not change):

- KS-DFT: the BODEN integral of each atom pair (`multipolar` for the pairs
  below `THREBOD`), `DFT EXCHANGE-CORRELATION, PAIR TERMS` and `DFT
  EXCHANGE-CORRELATION, ONE-CENTER TERMS` (the functional integrated
  within each atom, before half of each pair term is moved out of it);
- hybrid functionals: `HF-TYPE EXCHANGE, UNSCALED`, the full HF-type
  exchange before the exact-exchange fraction is applied;
- when the zero-error correction is applied, the two-electron matrices
  integrated on the first grid, titled `..., BEFORE THE ZERO-ERROR
  CORRECTION` (for a correlated wavefunction with `CORRELATION`, also the
  separate exchange and correlation).

Each comes with its fragment version under `DOFRAGS`.

### Correlated wavefunctions

CASSCF and CISD wavefunctions (also FCI) print the same sections, with
`EXCHANGE-CORRELATION ENERGY` built from the cumulant of the 2-RDM in
place of the exchange or the DFT terms, and the same
`TWO-ELECTRON INTEGRATION ERROR` section. With `CORRELATION` in
`# ENPART`, the exchange-correlation is also split into `EXCHANGE ENERGY`
and `CORRELATION ENERGY`. The zero-error correction applies to the
exchange-correlation as a whole, never to exchange or correlation
separately, so:

- without the correction (error within `TWOELTOLER`, or no reference
  energies), the separate exchange and correlation matrices are printed,
  and they add up to the exchange-correlation;
- with the correction, only the corrected `EXCHANGE-CORRELATION ENERGY`
  is printed, with a note; `VERBOSE 1` also prints the exchange and
  correlation before the correction.

## Checking the result

- **The final integration error is not a quality check by itself** when
  the zero-error strategy has been applied: it is close to zero by
  construction, even if the functional does not match the wavefunction.
- The informative numbers are the one-electron integration errors, the
  two-electron error **before** the zero-error strategy (the first
  `Integration error` of `TWO-ELECTRON INTEGRATION ERROR`, a few kcal/mol
  at most with sensible grids), and the damping parameter (about 1.2 to
  1.4 in the tests, with the 40 × 146 two-electron grid). A much larger error or damping parameter
  points to a functional that does not match the wavefunction.
- Without reference energies, the final integration error is the raw
  numerical error and depends strongly on the two-electron grid,
  especially for heavier atoms.

## Grids and cost

- **One-electron grid** (also used for the DFT exchange-correlation):
  150 radial × 590 angular points per atom by default for ENPART.
- **Two-electron grid**: 150 × 590 by default, set with `RADIAL_2E` and
  `ANGULAR_2E` in the optional `# GRID` block (see
  [Block section # GRID](../input/grid.md)); the rotation angles of the
  grid are chosen by the program for the angular grid.
- The 6-D two-electron integrations dominate the cost, which grows with
  the square of the number of points per atom: 40 × 146 is much cheaper
  than 150 × 590 and, with the zero-error strategy, useful for tests and
  larger systems; for publication-quality ENPART numbers keep the
  default.
  `THREBOD` saves time by skipping weakly bonded pairs.
- The DFT part stores the orbital gradients on the one-electron grid:
  about 2 MB per occupied orbital and atom (twice that for unrestricted
  wavefunctions).
