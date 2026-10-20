# Population analysis - Atomic populations, bond orders and valences

Every run starts with a population analysis in the chosen
[atomic definition](../guide/aim.md): electron and spin populations and
charges of each atom, the bond order between every pair of atoms, and
atomic valences. No keyword is needed. With a real-space scheme these are
the "fuzzy atoms" generalizations of the familiar Mulliken quantities,
much less dependent on the basis set.

```{admonition} Cite
:class: note

I. Mayer and P. Salvador, *Chem. Phys. Lett.*, **2004**, 383, 368-375.
See [Citations](../citations.md).
```

## How it works

With $\mathbf P$ the density matrix, $\mathbf P^s$ the spin-density matrix
and $\mathbf S^A$ the overlap matrix of the basis functions integrated
over atom A (with its weight function in real space, or the Hilbert-space
equivalent):

| Quantity | Definition |
|---|---|
| Electron population | $N_A = \mathrm{Tr}(\mathbf P\mathbf S^A)$ |
| Atomic charge | $q_A = Z_A - N_A$ |
| Spin population | $N^s_A = \mathrm{Tr}(\mathbf P^s\mathbf S^A)$ |
| Bond order | $B_{AB} = \mathrm{Tr}(\mathbf P\mathbf S^A\mathbf P\mathbf S^B) + \mathrm{Tr}(\mathbf P^s\mathbf S^A\mathbf P^s\mathbf S^B)$ |
| Total valence | $V_A = 2N_A - \mathrm{Tr}(\mathbf P\mathbf S^A\mathbf P\mathbf S^A)$ |
| Valence used in bonds | $V^B_A = \sum_{B\neq A} B_{AB}$ |
| Free valence | $F_A = V_A - V^B_A$ |

$Z_A$ is the nuclear charge (the effective one if a pseudopotential is
used). The free valence is zero for a closed-shell single determinant;
a non-zero value points to unpaired electrons, or to radical character in
a correlated wavefunction.

## Requirements and input

None: the population analysis is done in every run, with any atomic
definition and any wavefunction (the [atomic definition](../guide/aim.md)
is the only keyword it uses). With `DOFRAGS`, it is also given for the
fragments.

## Reading the output

The excerpts come from fluoromethane, RKS, TFVC (test `CH3F`). The
section starts with the atomic definition used
(`Atomic definition : TFVC (real space)`). The population, charge and
spin tables have a column named after it (`TFVC` here, `Lowdin`,
`Hirshfeld`, ... in other runs) and, for comparison, the Mulliken value;
the other tables, and every table of a Mulliken run, have a single
column:

% source: tests/reference/CH3F.apost
```text
  ------------------------
    ELECTRON POPULATIONS
  ------------------------

  Atom           TFVC    Mulliken
  -------------------------------
    1  C     5.126631    5.845837
    2  F     9.721908    9.259127
    3  H     1.050605    0.965021
    4  H     1.050654    0.965037
    5  H     1.050704    0.964978
  -------------------------------
  Total     18.000503   18.000000
```

The sum of the real-space populations differs from the number of electrons
by the numerical integration error (here 0.0005). `TOTAL ATOMIC CHARGES`
follows in the same layout, and `SPIN POPULATIONS` for open-shell
wavefunctions. With `EOS` or `GEOS` and a real-space or Mulliken atomic
definition, the overlap populations between atoms are also printed
(`OVERLAP POPULATION MATRIX`).

% source: tests/reference/CH3F.apost
```text
  ---------------------
    BOND ORDER MATRIX
  ---------------------

  -------------------------------------------------------------------
                 1  C        2  F        3  H        4  H        5  H
  -------------------------------------------------------------------
    1  C     3.346704    0.812947    0.915475    0.915636    0.916267
    2  F     0.812947    9.142936    0.115002    0.115018    0.114955
    3  H     0.915475    0.115002    0.486215    0.049132    0.049102
  ...
```

- Off-diagonal elements are the bond orders: 0.81 for C–F and 0.92 for
  C–H, close to single bonds; the small F···H and H···H values are
  through-space contributions.
- Diagonal elements are the electrons localized on each atom; for a
  closed-shell molecule, the diagonal element plus half the off-diagonal
  elements of its row add up to the population of the atom.

Then `TOTAL VALENCES` ($V_A$, 3.56 for carbon), `VALENCES USED IN BONDS`
($V^B_A$) and `FREE VALENCES` ($F_A$, zero up to the integration error in
this closed-shell molecule).

With `DOFRAGS`, every table has its fragment version
(`FRAGMENT ELECTRON POPULATIONS`, ...; `FRAGMENT BOND ORDERS` comes last,
after the valences), where the populations and charges are summed over the atoms of each
fragment and the bond order between two fragments is the sum of the bond
orders between their atoms (see [Fragments](../guide/fragments.md)).

**Correlated wavefunctions.** With a [`# DM` block](../input/dm.md)
(`DM 1` or `DM 2`), $\mathbf P$ and $\mathbf P^s$ are those of the 1-RDM
and the matrix is titled `EXCHANGE-ONLY BOND ORDER MATRIX`. After the
valences, `NEW BOND ORDER MATRIX` gives the bond orders from the natural
orbitals (I. Mayer's definition for correlated wavefunctions),

$$
B_{AB} = \sum_{ij} \left(n_i n_j + \sqrt{u_i u_j}\right) S^A_{ij} S^B_{ji},
\qquad u_i = 2n_i - n_i^2,
$$

with $n_i$ the occupations of the natural orbitals and $\mathbf S^A$ in
their basis; its total (diagonal and one triangle) is the number of
electrons. With `DM 2`, the [local spin analysis](spin.md) adds the
exchange-correlation bond orders from the 2-RDM. Without `# DM`, the bond
orders come from the density of the `.fchk` file, which is right for a
singlet only: for any other spin state of a correlated wavefunction the
bond orders and valences are not printed, and a note asks for the 1-RDM.

## Checking the result

- With a real-space scheme, the sum of the populations should equal the
  number of electrons (and the sum of the charges, the charge of the
  molecule) up to a small integration error; a larger difference calls
  for a larger grid (see [Running a calculation](../guide/running.md#integration-grid)).
- For a closed-shell single determinant, the free valences are zero up to
  the integration error.
- Mulliken values (also printed for comparison) depend strongly on the
  basis set; the real-space ones much less.

## Atomic overlap matrices for other programs

`DOINT` in `# METHOD` writes the atomic overlap matrices
$\mathbf S^A$ in the molecular-orbital basis, in the `.int` format of
AIMPAC, so that other programs (e.g. ESI-3D, for aromaticity indices)
can use the atoms of APOST-3D. See [Output files](../guide/output.md) for
the file names. The basis is that of the occupied orbitals (alpha, then
beta, for an unrestricted wavefunction) or, for a correlated
wavefunction, of the natural orbitals.

The output gives the number of files, the file that lists them, the
number of natural orbitals (correlated wavefunctions) and, when every
atom is included, the largest deviation of the summed matrices from the
identity. They add up to the identity up to the integration error; a
warning (above 0.01 on the diagonal, 0.001 off it) points to the
integration grid.
