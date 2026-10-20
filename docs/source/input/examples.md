# Input examples

Complete inputs for the most common analyses. Several analyses can be
combined in the same input, as long as they use the same atomic
definition. More complete inputs, with their outputs, are the test
cases in `tests/inputs/` and `tests/reference/` of the repository (see
[Running the test suite](../developer/testing.md#active-test-cases)).

Each example says which wavefunction of `tests/inputs/` it runs on. To
try it, copy that file (and the other files named) into an empty folder
under one job name, e.g. `mol.fchk`, save the input as `mol.inp`, and run
`$APOST3D_PATH/apost3d mol`.

## Example 1 - EOS with NAO atoms and cube files

```text
# METHOD
NAO-BASIS
EOS
DOFRAGS
CUBE
#
# FRAGMENTS
2
1
1
-1
#
# CUBE
MAX_OCC 700
MIN_OCC 300
#
```

Effective oxidation states with natural atomic orbitals as atoms (which
needs the file `jobname.nao`, see
[Preparing the wavefunction](../guide/wavefunctions.md#gaussian)), and cube
files of the effective fragment orbitals with net occupations between 0.3
and 0.7. Fragment 1 is atom 1; fragment 2 all the other atoms.

Try it on `CH3F-NAO.fchk` with `CH3F-NAO.nao` (copied as `mol.nao`): the
fragments are C and FH₃, with oxidation states −2 and +2 and six cube
files.

## Example 2 - ENPART with a functional given by libxc ids and a custom grid

```text
# METHOD
TFVC
ENPART
DOFRAGS
#
# ENPART
LIBRARY
EX_FUNCTIONAL 106
EC_FUNCTIONAL 131
THREBOD -1
#
# FRAGMENTS
2
1
1
-1
#
# GRID
RADIAL_2E 150
ANGULAR_2E 590
#
```

Energy partitioning of a KS-DFT wavefunction. The libxc ids 106 (B88
exchange) and 131 (LYP correlation) make up BLYP, so the `BLYP` keyword
gives the same result; the functional must be the one of the
wavefunction. `THREBOD -1` computes the exchange-correlation term
of every atom pair. The optional `# GRID` block here
repeats the default two-electron grid. See [ENPART](../methods/enpart.md).

Try it on `H2O-T-BLYP.fchk` (triplet water, UKS BLYP; under a minute on 8
threads): the fragments are O and H₂.

## Example 3 - Local spin of a correlated wavefunction from pySCF

```text
# METHOD
TFVC
SPIN
DM 2
#
# DM
mol.dm1
mol.dm2
pySCF
#
```

Local spin analysis of a CASSCF wavefunction. The 1- and 2-RDM files,
written by `utils/apost3d.py` together with the `.fchk` file (see
[Preparing the wavefunction](../guide/wavefunctions.md#pyscf)), are given
in the `# DM` block. See [SPIN](../methods/spin.md).

Try it on `LiH-35-CAS22.fchk` with `LiH-35-CAS22.dm1` and
`LiH-35-CAS22.dm2` (copied as `mol.dm1` and `mol.dm2`): LiH at 3.5 Å,
CASSCF(2,2).

## Example 4 - OSLO oxidation states

```text
# METHOD
TFVC
OSLO
DOFRAGS
#
# OSLO
LOWDIN
FOLI_TOLERANCE 3
PRINT_NONORTHO
#
# FRAGMENTS
3
1
1
2
2 3
-1
#
```

Oxidation states from localized orbitals, with three fragments: atom 1,
atoms 2 and 3, and the remaining atoms. The orbitals are localized on the
`TFVC` grid, while the fragment populations that enter the FOLI use
Löwdin atoms (`LOWDIN` in `# OSLO`). `PRINT_NONORTHO` also writes the
OSLOs before orthogonalization to a second `.fchk` file. See
[OSLO](../methods/oslo.md).

Try it on `FeCN5NO3--UBLYP-t2.fchk`: the fragments are Fe, one CN and the
rest, with oxidation states +2, −1 and −4.

## Example 5 - EOS and local spin together

```text
# METHOD
TFVC
EOS
SPIN
DOFRAGS
#
# FRAGMENTS
2
1
1
-1
#
```

Effective oxidation states and local spins of an open-shell single
determinant, with two fragments: atom 1, and the remaining atoms. For a
correlated wavefunction, add `DM 2` and a `# DM` block as in Example 3.

Try it on `NaBH3--UHF.fchk` (NaBH₃⁻, UHF): the fragments are Na and BH₃,
with oxidation states 0 and −1.

## Example 6 - GEOS with cube files of the negative-occupation orbitals

```text
# METHOD
TFVC
GEOS
DOFRAGS
CUBE
#
# FRAGMENTS
2
1
1
-1
#
# CUBE
MAX_OCC 0
MIN_OCC 0
NEG_EFOS
RADIUS_SCALE 12.0
#
```

Generalized effective oxidation states of LiH at 3.2 Å, from a pySCF
FCI/cc-pVTZ `.fchk` file, with Li and H as fragments. Try it on
`LiH-32-FCI.fchk`.

`MAX_OCC 0` / `MIN_OCC 0` select none of the regular orbitals, so the only
cube files written are those requested by `NEG_EFOS`: one per paired
orbital with net occupation at or below −0.025 (the default), in a box
large enough for the diffuse orbital on H. Here there is one:

% source: tests/reference/LiH-32-FCI.apost
```text
  WARNING: 1 paired EFO(s) with net occupation below -0.0250 left out of the assignment

  Frag.  Net occ.  Gross occ.        Fit (%)
  ------------------------------------------
      2   -0.0856     -0.0978          84.36
  ------------------------------------------
```

written as `LiH-32-FCI_tfvc_paired_neg_FR2_1.cube`. The same orbital is the
last Alpha orbital of `LiH-32-FCI-GEOS-EFOs.fchk`. The resulting oxidation
states are +1 (Li) and −1 (H). See [GEOS](../methods/geos.md).
