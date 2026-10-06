# GEOS - Generalized effective oxidation states

**Generalized effective oxidation states** (`GEOS`) assigns electrons, and
from them oxidation states, to atoms or fragments using two density
functions instead of one: the **paired** and the **unpaired** electron
density. It extends the [EOS](eos.md) analysis to open-shell,
broken-symmetry and multiconfigurational (correlated) wavefunctions, where
the alpha/beta picture used by EOS is not the natural one.

```{admonition} Cite
:class: note

M. Gimferrer and P. Salvador, *manuscript in preparation*, **2026**,
together with the EOS reference: E. Ramos-Cordoba, V. Postils and
P. Salvador, *J. Chem. Theory Comput.*, **2015**, 11, 1501-1508. See
[Citations](../citations.md).
```

## How it works

1. **Natural orbitals.** The one-particle density of the wavefunction is
   diagonalized into natural orbitals with occupations *n*.
2. **Paired and unpaired densities.** Each natural orbital contributes
   *u* = *n*(2 − *n*) to the unpaired density (Takatsuka's definition).
   The paired density is the total density minus the unpaired one. For a
   closed-shell single determinant every *n* is 0 or 2, so the unpaired
   density is zero and GEOS reduces to a closed-shell EOS analysis.
3. **Effective fragment orbitals (EFOs).** For each fragment and each of
   the two densities, the effective fragment orbitals and their net and
   gross occupations are obtained, as in EOS. EFOs with net occupation
   below 0.0001 are discarded.
4. **Electron assignment.** The EFOs of all fragments are pooled and
   sorted by gross occupation. A paired EFO takes two electrons, an
   unpaired one takes one. The starting assignment places the |N_α − N_β|
   electrons that must be unpaired, and the rest as electron pairs. The
   program then repeatedly tries to break the least-occupied pair into
   the two most-occupied unpaired EFOs, keeping the move only while it
   lowers the RMSD between the assigned and the actual occupations.
5. **Oxidation states.** Each fragment's oxidation state follows from the
   number of electrons assigned to it, and a reliability index R(%) is
   given per density and overall, as in EOS.

## Requirements and input

- A **real-space** atomic definition (`TFVC` recommended). GEOS stops with
  an error for Hilbert-space schemes (`MULLIKEN`, `LOWDIN`, ...).
- Any wavefunction whose `.fchk` provides the total density: restricted,
  unrestricted (including broken-symmetry), or correlated (e.g. a pySCF
  CASSCF/FCI `.fchk`; no `.dm1`/`.dm2` needed).
- Fragments are required (`DOFRAGS` and a
  [`# FRAGMENTS` block](../input/fragments.md)), with every atom in exactly
  one of them; the run stops otherwise. For one oxidation state per atom,
  make every atom a fragment of its own.

```text
# METHOD
TFVC
GEOS
DOFRAGS
#
# FRAGMENTS
2
1
1
-1
#
```

`EFFAO-U` computes the same paired/unpaired EFOs **without** the electron
assignment, the oxidation states and the `.fchk` export described below.

## Reading the output

The excerpts below come from LiH at 3.2 Å, pySCF FCI/cc-pVTZ, with Li and
H as fragments 1 and 2 (the full input is
[Example 6](../input/examples.md)).

**EFO blocks.** One block per density (`EFFAOs FROM THE PAIRED DENSITY`,
then `... UNPAIRED DENSITY`), with one entry per fragment:

```text
  ** FRAGMENT   2 **

  Net occupation for fragment      2    1.05137
  Net occupation using >    0.00010
  OCCUP.   1.0512   0.0002

  Gross occupation for fragment    2    1.08681
  OCCUP.   1.0861   0.0007
  FIT %     93.01    33.07
  Left out by the cutoff (net / gross):   -0.09431   -0.10682
```

- The first `OCCUP.` rows are the net occupations of the fragment's EFOs,
  and the second ones their gross occupations, 8 per line, in the same
  order. Paired occupations go up to 2, unpaired ones up to 1.
- `FIT %` gives, for each EFO, how faithfully the orbital written to the
  `.fchk` file reproduces it (see *Orbitals in the .fchk file* below).
- `Left out by the cutoff` is what the listed EFOs miss of the fragment's
  net and gross population in that density (see [EFFAO](effao.md)). In
  the paired density it can be negative, as here: the EFOs below the
  cutoff include the negative ones described next.

**Negative paired EFOs.** The paired density is the difference of two
densities, so unlike a normal density it can yield EFOs with a *negative*
occupation. Those with net occupation below −0.025 are listed
separately:

```text
  Found 1 EFO(s) with significant negative occupation in the paired channel
  (net occupation < -0.0250), excluded from oxidation-state assignment:
    Frag. 2    Net occ.  -0.0856    Gross occ.  -0.0978    % recovered  84.36
```

They never take part in the electron assignment. They are written to the
`.fchk` file (see below), and cube files of them can be requested with
`NEG_EFOS` (see *Cube files* below).

**Assignment and oxidation states.**

```text
  Ideal assignment: 0 forced-unpaired electron(s), 2 electron pair(s)
  Initial RMSD: 0.53195

  Trying: paired EFO (frag 2, occ  1.0861) -> unpaired (frag 2, occ  0.4272 / frag 1, occ  0.2352)
    RMSD =  0.73340 -- rejected, keeping RMSD =  0.53195

  Electrons assigned per fragment (paired / unpaired):
    Fragment 1:  2.00 /  0.00
    Fragment 2:  2.00 /  0.00
```

Then come the usual EOS tables, one per density (`EOS ANALYSIS FOR PAIRED
ELECTRONS` / `... UNPAIRED ELECTRONS`: electrons per fragment, last
occupied and first unoccupied gross occupation, and R(%); paired
occupations are halved when computing R(%)), followed by:

```text
   Frag.  Oxidation State
  ------------------------
     1          1.00
     2         -1.00
  ------------------------
   Sum:  0.0

  OVERALL RELIABILITY INDEX R(%) =  96.638
  ELECTRONIC ASSIGNMENT RMSD VALUE =   0.532
```

## Orbitals in the .fchk file

Every GEOS run writes `<jobname>-GEOS-EFOs.fchk`, with the EFOs of all
fragments as its orbitals (see [Visualizing orbitals](../guide/visualization.md)
for the general layout and the meaning of `FIT %`):

- **Alpha orbitals:** paired EFOs of all fragments, sorted by decreasing
  gross occupation; the negative paired EFOs, if any, are the **last**
  Alpha orbitals (most negative last).
- **Beta orbitals:** unpaired EFOs, same ordering. For a restricted
  input `.fchk` with a non-empty unpaired density, a Beta set is added.
- **"Orbital energy"** of each orbital: its gross occupation.

In the LiH example, the negative EFO of the H atom is the last Alpha
orbital, with "orbital energy" −0.0978 (its gross occupation), and its
`% recovered` of 84.36 is the `FIT %` of that orbital.

## Cube files

With `CUBE` in `# METHOD` and a [`# CUBE` block](../input/cube.md):

- `<jobname>_<scheme>_paired_FR<n>_<i>.cube` and `..._unpaired_...` for
  the EFOs selected by `MAX_OCC`/`MIN_OCC`. Paired occupations go up to 2,
  so both bounds are doubled for the paired density.
- With `NEG_EFOS [val]`, also `..._paired_neg_FR<n>_<i>.cube` for every
  paired EFO with net occupation ≤ −val/1000 (default val = 25), index 1
  being the most negative. This selection is independent of the −0.025
  used for the output and the `.fchk`: a smaller value plots more EFOs,
  some of which then have no `.fchk` counterpart.

For only the negative EFOs, use an empty occupation window, and a larger
box (negative EFOs are diffuse, and the default box of a single H atom
cuts them):

```text
# CUBE
MAX_OCC 0
MIN_OCC 0
NEG_EFOS 10
RADIUS_SCALE 12.0
#
```

## Checking the result

- R(%) is read as in [EOS](eos.md#how-it-works), per density and overall.
- `ELECTRONIC ASSIGNMENT RMSD VALUE` is the root-mean-square deviation
  between the assigned and the actual EFO occupations: the smaller, the
  closer the electron distribution is to the integer assignment.
- For negative paired EFOs, `% recovered` tells how faithfully the
  orbital in the `.fchk` file represents them.
