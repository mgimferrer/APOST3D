# EFFAO - Effective atomic and fragment orbitals

The **effective atomic orbitals** (eff-AOs) of an atom in a molecule are
the orbitals, with their occupations, that describe the part of the
electron density belonging to that atom. For a group of atoms (a ligand,
a metal center, a molecule in a complex) they are called **effective
fragment orbitals** (EFOs). They need no reference state, Lewis structure
or threshold, they don't depend on the basis set, and they work for any
wavefunction, from Hartree-Fock and KS-DFT to CASSCF.

Two properties make them useful:

- **They recover the minimal basis.** Whatever the basis set of the
  calculation, an atom has only as many eff-AOs with significant occupation
  as orbitals in its classical minimal basis: one for H, five for C, N or O
  (1s, 2s, 2p), and so on. The occupied ones are core orbitals, lone pairs
  and valence hybrids; the rest have negligible occupations.
- **Fragment orbitals look like the orbitals of the free fragment**,
  polarized by the environment, with fractional occupations. For a CO
  ligand bound to a metal, for instance, the σ lone-pair EFO has an
  occupation below 1 (σ donation to the metal) and the formally empty π\*
  EFOs a non-zero one (back-donation from the metal). The occupations
  quantify donation, back-donation and covalency directly.

The EFOs are the basis of the [EOS](eos.md) and [GEOS](geos.md)
oxidation-state analyses.

```{admonition} Cite
:class: note

I. Mayer, *J. Phys. Chem.*, **1996**, 100, 6249-6257; I. Mayer and
P. Salvador, *J. Chem. Phys.*, **2009**, 130, 234106; E. Ramos-Cordoba,
P. Salvador and I. Mayer, *J. Chem. Phys.*, **2013**, 138, 214107.
Review of the method and its applications: M. Gimferrer, *Theor. Chem.
Acc.*, **2026**, 145, 74. See [Citations](../citations.md).
```

## How it works

With a real-space atomic definition, the weight function of fragment F is
the sum of the weight functions of its atoms,
$w_F(\mathbf r)=\sum_{A\in F} w_A(\mathbf r)$. From the natural orbitals
$\phi_i$ and occupations $n_i$ of the wavefunction, the matrix

$$
Q^F_{ij} = \sqrt{n_i n_j}\int w_F(\mathbf r)\,\phi_i^*(\mathbf r)\,\phi_j(\mathbf r)\,w_F(\mathbf r)\,d\mathbf r
$$

is diagonalized. Its eigenvalues $\lambda^F_\mu$ are the **net
occupations** of the EFOs, and the EFOs are the corresponding combinations
of natural orbitals cut to the fragment by $w_F$. They are orthonormal and
together reproduce the net density of the fragment.

When atoms overlap (as fuzzy atoms do), the net occupations don't add up
to the electron population of the fragment: the overlap part is missing.
APOST-3D therefore also gives **gross occupations**, the population on
the fragment of each EFO before it is cut by $w_F$. Gross occupations add
up to the fragment's gross population, and over all fragments to the
number of electrons.

With a Hilbert-space definition (`MULLIKEN`, `LOWDIN`, `LOWDIN-DAVIDSON`,
`NAO-BASIS`), the same construction uses the fragment's block of the
density matrix in the basis of that scheme (Mayer's original formulation):
the original basis functions for `MULLIKEN`, the orthogonalized ones for
the others, the same basis as their populations. There is a single set of
occupations.

Three variants differ in the density that is analyzed:

| Keyword | Density | Occupations |
|---|---|---|
| `EFFAO` | Total density | 0 to 2 |
| `UEFFAO` | Alpha and beta densities, separately | 0 to 1 |
| `EFFAO-U` | Paired and unpaired densities, separately (the orbitals of [GEOS](geos.md), without its oxidation states). Real-space schemes only | 0 to 2 (paired), 0 to 1 (unpaired) |

For a closed-shell wavefunction, `UEFFAO` computes the alpha orbitals
only, since the beta ones are the same.

## Requirements and input

- Any atomic definition (`TFVC` recommended, see
  [Atoms in molecules](../guide/aim.md)); `EFFAO-U` needs a real-space one.
- Any wavefunction.
- Fragments are optional. Without `DOFRAGS` every atom is a fragment, and
  the result is the set of effective atomic orbitals.

```text
# METHOD
TFVC
UEFFAO
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

## Reading the output

The excerpt below comes from the FeCO₂²⁺ complex, RKS PBE, TFVC, with Fe
as fragment 1 and the two CO ligands as fragments 2 and 3 (test
`FeCO2-PBEPBE`, which runs `EOS`; `UEFFAO` prints the same blocks). Each
fragment gets one entry per density:

```text
  ** FRAGMENT   1 **

  Net occupation for fragment      1   12.11689
  Net occupation using >    0.00100
  OCCUP.   1.0000   1.0000   1.0000   1.0000   1.0000   0.9963   0.9963   0.9957
  OCCUP.   0.9854   0.9816   0.8156   0.8144   0.2759   0.1644   0.0609   0.0120
  OCCUP.   0.0086   0.0081   0.0018

  Gross occupation for fragment    1   12.27814
  OCCUP.   1.0000   1.0000   1.0000   1.0000   1.0000   0.9978   0.9978   0.9974
  OCCUP.   0.9908   0.9886   0.8381   0.8374   0.3210   0.1676   0.0981   0.0160
  OCCUP.   0.0124   0.0117   0.0034

  FIT %     99.99    99.99    99.93    99.93    99.90    97.92    97.92    98.89
  FIT %     94.65    95.36    90.95    90.85    86.78    93.03    59.76    84.76
  FIT %     70.62    70.47    58.34
  Left out by the cutoff (net / gross):    0.00007    0.00024
```

- The first `OCCUP.` rows list the net occupations of the fragment's EFOs,
  in decreasing order, 8 per line. Only EFOs with net occupation above
  0.001 are listed, and `Net occupation for fragment` is their sum.
- `Gross occupation for fragment` and the second `OCCUP.` rows are the
  same for the gross occupations, in the same order.
- `Left out by the cutoff` is what the listed EFOs miss of the fragment's
  net and gross (alpha) population: the occupations of the EFOs below
  0.001. Small values mean the listed EFOs describe the fragment's density
  well. With Mulliken, Löwdin or NAO atoms only the net value is printed.
- `FIT %` (real-space schemes) measures how well each EFO can be
  written in terms of the basis functions, which is how `EOS` and `GEOS`
  write them to a `.fchk` file (see
  [Visualizing orbitals](../guide/visualization.md)).

Here the iron alpha EFOs show five core-like orbitals at 1.000, five
more above 0.98 (the rest of the core and the 3d-type orbitals that keep
their electrons), and then a gap: two EFOs at 0.82, one at 0.28 and small
occupations after it (0.84 and 0.32 as gross occupations). Where that gap falls, compared with the EFOs of the
ligands, is what the [EOS](eos.md) analysis turns into an oxidation
state.

With `MULLIKEN`, `LOWDIN` or `NAO-BASIS` only one set of occupations is
printed (there is no separate gross occupation and no `FIT %`). Mulliken
occupations can slightly exceed 1.

## Files written

- Cube files of selected EFOs, with `CUBE` and a `# CUBE` block (see
  [Visualizing orbitals](../guide/visualization.md)).
- `EFFAO`, `UEFFAO` and `EFFAO-U` on their own write no `.fchk` file of the
  orbitals: that is written by `EOS` and `GEOS`.

## Checking the result

- The number of EFOs with a significant occupation should match the
  minimal basis of the fragment (e.g. 5 per C, N or O atom, 1 per H). Many
  more points to a very diffuse basis set with a Hilbert-space scheme; use
  a real-space one.
- With a real-space scheme, the gross occupations of all fragments add up
  to the number of electrons (of each spin, for `UEFFAO`).
