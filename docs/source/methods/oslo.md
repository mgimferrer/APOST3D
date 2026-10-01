# OSLO - Oxidation states from localized orbitals

`OSLO` assigns oxidation states to fragments by localizing the occupied
orbitals onto the fragments, one at a time, and giving each electron pair
(or electron, for open shells) to the fragment its orbital is localized
on. The **oxidation state localized orbitals** (OSLOs) themselves are a
useful picture of the bonding between fragments, and the **fragment
orbital localization index** (FOLI) measures how clearly each orbital
belongs to one fragment.

```{admonition} Cite
:class: note

M. Gimferrer, A. Aldossary, P. Salvador and M. Head-Gordon, *J. Chem.
Theory Comput.*, **2022**, 18, 309-322. See [Citations](../citations.md).
```

## How it works

1. **Fragment-centered localization.** For each fragment F, the occupied
   orbitals are combined into orbitals of minimal radial spread around a
   reference point $\mathbf R_F$, the fragment's center of nuclear charge
   (printed in `CHARGE CENTER (R_F) FOR EACH FRAGMENT`, in Å). This is a
   single diagonalization of the spread matrix
   $L^F_{ij} = \int \psi_i(\mathbf r)\,|\mathbf r-\mathbf R_F|^2\,\psi_j(\mathbf r)\,d\mathbf r$
   in the occupied space: no iterations, no multiple minima. For an atom,
   the orbitals of smallest spread are its core orbitals.
2. **FOLI.** For each of these orbitals, the populations $N^i_G$ on all
   fragments G give Pipek's delocalization measure
   $D_i = 1/\sum_G (N^i_G)^2$, and the FOLI on its own fragment F,

   $$
   D^F_i = \sqrt{D_i / N^i_F}.
   $$

   A FOLI of 1 means an orbital entirely on fragment F; 2 means an orbital
   shared equally between two fragments; larger values mean more
   delocalization, or an orbital centered on F but mostly elsewhere.
3. **Iterative selection.** Among the candidates of all fragments, the one
   with the lowest FOLI is selected, together with every other candidate
   within the tolerance (0.001 by default) of it. The selected orbitals
   are orthogonalized and removed from the occupied space, and step 1 is
   repeated in the space left, until every occupied orbital has been
   assigned. Core and lone-pair orbitals go first; the orbitals of the
   bonds between fragments, the ones that decide the oxidation states, come
   last, and are better localized for having the others removed first.
4. **Oxidation states.** Each OSLO gives its electrons (two in a
   restricted calculation, one per spin orbital in an unrestricted one) to
   the fragment it was localized on. The oxidation state of a fragment is
   the sum of the nuclear charges of its atoms minus its electrons.

For an unrestricted wavefunction the procedure is done separately for the
alpha and the beta orbitals.

The **Δ-FOLI** of an iteration is the gap between the selected FOLI and the
next candidate. In the **last iteration** it measures how clear the
assignment of the least localized orbital, and so the whole result, is:
above about 0.5 the assignment is usually clear (above 1, clearly ionic);
smaller values mean an increasingly covalent bond between the fragments.
A last FOLI close to 1 means the last orbital is well localized.

## Requirements and input

- A **single-determinant** wavefunction (HF or KS-DFT, restricted or
  unrestricted). The run stops for CASSCF, CISD or FCI wavefunctions.
- **Fragments** (`DOFRAGS` and a [`# FRAGMENTS` block](../input/fragments.md)),
  with every atom in a fragment. OSLO is not meaningful for single atoms,
  and the program does not check that fragments are defined.
- A **real-space atomic definition** in `# METHOD` (use `TFVC`): the
  spread matrix is integrated on its grid.
- A [`# OSLO` block](../input/oslo.md). It can choose another scheme for
  the fragment populations that enter the FOLI (`MULLIKEN`, `LOWDIN`,
  `LOWDIN-DAVIDSON`, `NAO-BASIS`); without one, the real-space scheme of
  `# METHOD` is used.

```text
# METHOD
TFVC
OSLO
DOFRAGS
#
# OSLO
FOLI TOLERANCE 3
#
# FRAGMENTS
2
1
2
-1
#
```

Here, for fluoromethane, fragment 1 is the F atom (atom 2) and fragment 2
the CH₃ group (all the other atoms).

## Reading the output

The excerpts come from the input above (CH₃F, RKS, test `CH3F`).

**Iterations.** Each iteration lists, per fragment, its candidate
orbitals with their spread (bohr²), population on the fragment and FOLI,
then the selected ones:

```text
  ------------------------
    ITERATION NUMBER   1
  ------------------------

  ----------------------------------------
    ORBITAL INFORMATION FOR FRAGMENT   1
  ----------------------------------------

   Orb.    Spread      Frg. Pop.      FOLI
  ------------------------------------------
     1    14.20213      0.02652      6.30549
  ...
     8     1.27930      0.99746      1.00381
     9     0.04119      0.99999      1.00002
  ------------------------------------------
  ...
  ---------------------
    SELECTED ORBITALS
  ---------------------

   Orb.   Frag.   FOLI
  ---------------------
     9      1   1.00002
     9      2   1.00067
  ---------------------

  Number of OSLOs selected:   2
  delta-FOLI value:    0.00380

  Orbitals left to assign:   7
```

In the first iteration the F 1s (FOLI 1.00002) and the C 1s (1.00067,
within the tolerance) are selected together. In the ninth and last
iteration the only orbital left, the C–F bond, is assigned to fluorine
with a FOLI of 1.355 against 2.714 for the CH₃ candidate, a Δ-FOLI of
1.359: a clear ionic assignment.

**Summary.** The selected OSLOs, in the order they were assigned, with
their FOLI and fragment populations before (`pre-ortho`) and after
orthogonalization (`final`), followed by the oxidation states:

```text
  ---------------------------------------------
    Summary of the selected OSLOs (pre-ortho)
  ---------------------------------------------

  OSLO Number :      1         2         3         4         5
  FOLI Value  :   1.00002   1.00067   1.00381   1.03691   1.03688
  Frg. Pop.  1:   0.99999   0.00043   0.99746   0.97593   0.97595
  Frg. Pop.  2:   0.00001   0.99956   0.00253   0.02408   0.02401

  OSLO Number :      6         7         8         9
  FOLI Value  :   1.04773   1.06153   1.06204   1.35476
  Frg. Pop.  1:   0.03102   0.03984   0.03982   0.80051
  Frg. Pop.  2:   0.96906   0.96042   0.96012   0.19950
  ...
  -----------------------------
    FRAGMENT OXIDATION STATES
  -----------------------------

   Frag.  Oxidation State
  ------------------------
     1         -1.00
     2          1.00
  ------------------------
```

The last OSLO (number 9) is the C–F bond, 80% on fluorine: F(−1) and
CH₃(+1). For an unrestricted wavefunction the iterations and summaries
appear twice (`ALPHA PART`, `BETA PART`).

**Linear dependencies.** After each selection, the program checks whether
the selected orbitals and the next candidates are nearly linearly
dependent, and prints `WARNING : LINEAR DEPENDENCY FOUND` if so. This is a
diagnostic for borderline cases (two fragments competing for the same
orbital); the program goes on with the lowest FOLI. In the last
iteration, the warning with a `Lowest eigenvalue` of 0.00000 only means
that more candidates were checked than orbitals were left, which is common.
The alternative assignment that the warning suggests (branching) is not
available in this version.

## Files written

- `<jobname>-OSLOs.fchk`: a copy of the input `.fchk` with the final OSLOs
  as its orbitals (and the matching density), for any orbital viewer. Every
  OSLO run writes it.
- `<jobname>-OSLOs-preortho.fchk`: the OSLOs before orthogonalization,
  with `PRINT NON-ORTHO` in `# OSLO`.

See [Visualizing orbitals](../guide/visualization.md).

## Checking the result

- **Look at the last iteration**: its FOLI and Δ-FOLI tell how clear the
  assignment is. For a Δ-FOLI below about 0.5, look at the last OSLOs in
  the `.fchk` file: a visibly shared orbital means a covalent bond, where
  the ionic assignment is a close call.
- **Several orbitals selected together** in a late iteration, in a molecule
  without symmetry, can change the result: rerun with a tighter tolerance
  (`FOLI TOLERANCE 4`), so that orbitals are selected one by one.
- **Compare population schemes** in borderline cases, e.g. the default
  real-space one against `NAO-BASIS` in `# OSLO`.
