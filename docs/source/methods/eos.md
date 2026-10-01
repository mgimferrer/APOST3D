# EOS - Effective oxidation states

The **effective oxidation states** (EOS) analysis assigns formal oxidation
states to atoms or fragments (metal centers, ligands) directly from the
wavefunction. It makes the fragments "compete" for the electrons through
their [effective fragment orbitals](effao.md) (EFOs), and gives a
**reliability index** R(%) that says how clear-cut the assignment is. It
works for any wavefunction and any atomic definition, and needs no Lewis
structure or bond assignment.

```{admonition} Cite
:class: note

E. Ramos-Cordoba, V. Postils and P. Salvador, *J. Chem. Theory Comput.*,
**2015**, 11, 1501-1508. Review: M. Gimferrer, *Theor. Chem. Acc.*,
**2026**, 145, 74. See [Citations](../citations.md).
```

## How it works

1. **Spin-resolved EFOs.** For each fragment, the EFOs of the alpha and of
   the beta density are obtained ([EFFAO](effao.md)), with occupations
   between 0 and 1. The gross occupations are used, so that all fragments
   are compared on the same footing.
2. **Electron assignment.** For each spin, the EFOs of all fragments are
   pooled and sorted by decreasing occupation, and the $N_\alpha$ (or
   $N_\beta$) electrons are assigned one by one to the most occupied ones.
   Occupations are not rounded: a fragment can win an electron with an
   EFO occupied less than 0.5 if no other fragment has a more occupied one.
3. **Oxidation states.** A fragment that received $n_F$ electrons (alpha
   plus beta) has the oxidation state $\mathrm{OS}(F) = Z_F - n_F$, where
   $Z_F$ is the sum of the nuclear charges of its atoms (the effective
   ones if pseudopotentials are used).
4. **Reliability index.** The last EFO that received an electron (LO,
   last occupied) and the most occupied EFO of any other fragment that did
   not (FU, first unoccupied) are the frontier EFOs. For each spin σ,

   $$
   R_\sigma(\%) = 100\,\min\!\left(1,\ \lambda^\sigma_{\mathrm{LO}} - \lambda^\sigma_{\mathrm{FU}} + \tfrac12\right),
   $$

   and the overall index is $R = \min(R_\alpha, R_\beta)$. A gap of half an
   electron or more gives 100; two frontier EFOs with the same occupation
   on different fragments give 50, i.e. two equally plausible assignments.

   As a rule of thumb, $R > 80$ is a clear-cut ionic assignment; 60 to 80
   means significant covalency, but the assignment still holds; below 60
   the bond is highly covalent and the result may depend on the
   wavefunction and the atomic definition, so look at the frontier EFOs
   and compare with other descriptors.
5. **Near-degenerate frontier EFOs.** If EFOs of different fragments are
   within 0.0025 of the last occupied one, the electron(s) are shared
   among them, the program prints `EOS: WARNING, PSEUDO-DEGENERACIES
   DETECTED`, and the oxidation states become fractional (e.g. two
   equivalent ligands sharing one electron get −1.5 each).

For a closed-shell wavefunction the beta assignment is the same as the
alpha one, so only the alpha EFOs are computed.

## Requirements and input

- Any atomic definition. Real-space schemes (`TFVC` recommended) give the
  most robust results and also write the EFOs to a `.fchk` file; Mulliken,
  Löwdin and NAO are available too (see [Atoms in molecules](../guide/aim.md)).
- Any wavefunction: restricted or unrestricted single determinants, and
  correlated wavefunctions (their alpha and beta densities, from the
  `.fchk` file or from the 1-RDM in a [`# DM` block](../input/dm.md)).
- Fragments are optional: without `DOFRAGS`, every atom gets its own
  oxidation state. When fragments are defined, every atom must belong to
  exactly one of them. A good choice for a complex is one fragment per
  metal center and one per ligand.

```text
# METHOD
TFVC
EOS
DOFRAGS
#
# FRAGMENTS
3
1
1
2
2 4
-1
#
```

```{admonition} Open-shell singlets and diradicals
:class: tip

For a correlated **restricted** wavefunction (e.g. CASSCF) of a singlet
diradical, the alpha and beta densities are identical, so EOS assigns the
two electrons of the stretched bond together, to one fragment, instead of
one to each. [GEOS](geos.md), which works with the paired and unpaired
densities instead, is designed for these cases.
```

## Reading the output

The excerpts come from the FeCO₂²⁺ complex, closed-shell RKS PBE, TFVC,
with Fe as fragment 1 and the two CO ligands as fragments 2 and 3 (test
`FeCO2-PBEPBE`, input above). The run first prints the EFOs of each
fragment ([EFFAO](effao.md)), then the assignment for each spin:

```text
  Total number of eff-AO-s for analysis:   41

  -------------------------------------------------
    EOS: Unambiguous integer electron assignation
  -------------------------------------------------

  ------------------------------------
    EOS ANALYSIS FOR ALPHA ELECTRONS
  ------------------------------------

   Frag.  Elect.  Last occ.  First unocc.
  ----------------------------------------
     1    12.00     0.837     0.321
     2     7.00     0.787     0.084
     3     7.00     0.787     0.084
  ----------------------------------------
   RELIABILITY INDEX R(%) =  96.598
```

- `Elect.`: alpha electrons assigned to the fragment.
- `Last occ.` / `First unocc.`: gross occupations of the fragment's last
  EFO that received an electron and of its first one that did not.
- R(%) is computed from the smallest `Last occ.` (0.787, a CO EFO) and the
  largest `First unocc.` of the other fragments (0.321, on Fe):
  100 × (0.787 − 0.321 + 0.5) = 96.6.

Since the wavefunction is closed-shell, the beta part is skipped
(`SKIPPING EFFAOs FOR BETA ELECTRONS`) and the final table follows:

```text
  -----------------------------
    FRAGMENT OXIDATION STATES
  -----------------------------

   Frag.  Oxidation State
  ------------------------
     1          2.00
     2          0.00
     3          0.00
  ------------------------
   Total oxidation state:    2.0

  OVERALL RELIABILITY INDEX R(%) =  96.598
```

Iron receives 12 alpha and 12 beta electrons, 24 of its 26: Fe(II) with
two neutral CO ligands, an unambiguous assignment. The total oxidation
state is the charge of the molecule.

A more covalent case: for the doublet [Fe(CN)₅NO]³⁻ (UKS BLYP, Löwdin,
test `FeCN5NO3--UBLYP`), the alpha frontier EFOs are 0.650 (a cyanide)
and 0.505 (on Fe), giving R(%) = 64.5: the assignment Fe(II), five CN⁻
and a neutral NO still holds, but with significant covalency.

## Files written

- `<jobname>-EOS-EFOs.fchk` (real-space schemes): the EFOs of all
  fragments as orbitals, for any viewer that reads `.fchk` files; see
  [Visualizing orbitals](../guide/visualization.md).
- `efo_occ.dat`, `efo_coeff.dat` (`LOWDIN`, `LOWDIN-DAVIDSON`,
  `NAO-BASIS`): occupations and coefficients as text (see
  [Output files](../guide/output.md#other-files)).
- Cube files of selected EFOs with `CUBE`.

## Checking the result

- **Look at R(%)**, and for values below about 80, at the frontier EFOs
  themselves: their shape shows which bond is being split, and how
  polarized it is.
- **Compare atomic definitions** in borderline cases (e.g. `TFVC` against
  `NAO-BASIS`): a robust assignment doesn't change.
- **Fractional oxidation states** mean near-degenerate frontier EFOs on
  different fragments (symmetry-equivalent ligands, mixed valence), not an
  error.
- The total oxidation state must equal the charge of the molecule.
