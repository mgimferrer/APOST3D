# Input examples

Sample input files are provided below for the most common types of analyses
currently implemented in APOST-3D. More than one analysis tool using the
same AIM scheme can be requested in the same input file.

## Example 1 — EFFAO/EOS with NAO-BASIS, fragments, and cube files

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
MAX_OCC=700
MIN_OCC=300
#
```

The evaluation of the effective oxidation states is requested, together
with the cube generation of the effective fragment orbitals with net
occupancies within the [0.300, 0.700] range.

The AIM requested is `NAO-BASIS`, which requires an additional
`inputname.nao` file with the transformation matrix from the original AO to
NAO basis. This file corresponds to the `FILE.33` file generated, for
instance, with Gaussian coupled to the NBO software. For this AIM, in a
Gaussian input file the user must add the `pop=(full,nboread)` keyword,
together with the `$NBO AONAO=W $END` additional line at the end of the
file. This generates the `FILE.33` file, to be renamed as `inputname.nao`.

Fragments have been defined: fragment 1 consists of atom 1, and fragment 2
gathers the rest of the atoms of the molecule.

## Example 2 — ENPART (KS-DFT/IQA decomposition) with TFVC and a custom grid

```text
# METHOD
TFVC
ENPART
DOFRAGS
#
# ENPART
LIBRARY
EX_FUNCTIONAL=106
EC_FUNCTIONAL=132
THREBOD=-1
MOD-GRIDTWOEL
#
# FRAGMENTS
2
1
1
-1
#
# GRID
RADIAL 150
ANGULAR 590
rr00 0.5
phb1 0.169
phb2 0.170
#
```

The IQA decomposition of the KS-DFT energy is requested. Functional
ID=106 (B88) and ID=132 (P86) correspond to the well-known BP86 functional,
so the `BP86` keyword gives the same result. See
[ENPART](../methods/enpart.md) for how to read the output.

The AIM requested is TFVC and fragments have been defined: fragment 1
consists of atom 1, and fragment 2 gathers the rest of the atoms of the
molecule. The diatomic XC energy components are computed exactly for all
atomic pairs. User-specific options for the decomposition of the
two-electron energy are given in the `# GRID` block section (in this case
they coincide with the default, *strongly recommended*, values).

## Example 3 — Local spin analysis from a correlated wavefunction

```text
# METHOD
TFVC
SPIN
DM 2
#
# DM
filename.dm1
filename.dm2
pySCF
#
```

The local spin analysis of a correlated WF (e.g. CASSCF) from a pySCF run
is requested.

The files containing the RDM1 and RDM2 information (`.dm1` and `.dm2`) are
mandatory for this type of analysis. They can be obtained from a pySCF run
using the `apost3d.py` utility, as described in
[Extracting .fchk, .dm1 and .dm2 files from pySCF](pyscf.md).

The `# DM` section contains the name of the two files (with extension),
together with the keyword of the code they were created with.

The AIM requested is TFVC and fragments have not been defined.

## Example 4 — OSLO oxidation states

```text
# METHOD
TFVC
OSLO
DOFRAGS
#
# OSLO
LOWDIN
FOLI TOLERANCE 3
BRANCH ITERATION 0
PRINT NON-ORTHO
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

The evaluation of oxidation states using the OSLO method is requested.

`TFVC` in `# METHOD` builds the numerical-integration grid OSLO always
needs to localize orbitals onto fragments — this is mandatory
independently of the AIM scheme chosen below in `# OSLO`; see
[Block section # OSLO](oslo.md) for why. Fragments have been defined:
fragment 1 consists of atom 1, fragment 2 consists of two atoms (2 and
3), and fragment 3 gathers the rest of the atoms of the molecule.

The AIM requested for FOLI values (i.e. the fragment-population step, not
the numerical integration above) is LOWDIN, with a tolerance value of 3
(real(10^-3), the default). Printing of the selected OSLOs
pre-orthogonalization (non-orthogonal if two or more OSLOs are selected in
the same iteration AND belong to different fragments) is invoked. As a
result, a new `.fchk` file is created containing these OSLOs.

## Example 5 — EOS + local spin analysis

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
5
2 3 4 5 6
#
```

The evaluation of the effective oxidation states and the local spin
analysis is requested. For correlated (multireference) WFs, the local spin
analysis requires the RDM1 and RDM2 information (see Example 3).

The AIM requested is TFVC and fragments have been defined: fragment 1
consists of atom 1, and fragment 2 of atoms 2-6 (the system has 6 atoms).

## Example 6 — GEOS with cube files of the negative-occupation orbitals

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
#
```

The generalized effective oxidation states (GEOS) of LiH at 3.2 Å are
requested, from a pySCF FCI/cc-pVTZ `.fchk` file (the files are in
`tests/inputs/LiH-32-FCI.*` in the repository). Fragment 1 is the Li
atom and fragment 2 the H atom.

`MAX_OCC 0` / `MIN_OCC 0` select none of the regular orbitals, so the only
cube files written are those requested by `NEG_EFOS`: one per paired
orbital with net occupation at or below −0.025 (the default). Here there is
one, on the H atom:

```text
  Found 1 EFO(s) with significant negative occupation in the paired channel
  (net occupation < -0.0250), excluded from oxidation-state assignment:
    Frag. 2    Net occ.  -0.0856    Gross occ.  -0.0978    % recovered  84.36
```

written as `LiH-32-FCI_tfvc_paired_neg_FR2_1.cube`. The same orbital is the
last Alpha orbital of `LiH-32-FCI-GEOS-EFOs.fchk`, whose "orbital energy"
(−0.0978) is its gross occupation. The resulting oxidation states are +1
(Li) and −1 (H). See [GEOS](../methods/geos.md) for how to read the full
output.

## Real-case examples

The `tests/inputs` folder of the repository contains the `.fchk` and
`.inp` files of the test suite (see [Testing](../testing.md)), which are
complete, working inputs for a range of analyses; the matching outputs are
in `tests/reference`. A curated set of examples with their outputs will be
added in a future release.
