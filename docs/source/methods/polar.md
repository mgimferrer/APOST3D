# POLAR - Dipole moment decomposition

`POLAR` splits the dipole moment of the molecule into atomic
contributions: an **intrinsic** dipole of each atom (the polarization of
its own electron density) and a **charge-transfer** dipole (the charge
the atom has gained or lost from its neighbors). Both are independent of
the origin of the coordinates. Computed at several field strengths
(finite fields), their changes give the atomic contributions to the
static polarizability.

```{admonition} Cite
:class: note

M. Montilla, J. M. Luis and P. Salvador, *J. Chem. Theory Comput.*,
**2021**, 17, 1098-1105. See [Citations](../citations.md).
```

## How it works

With the weight function $w_A$ of each atom, its electron population
$N_A = \int w_A\,\rho$ and its charge $q_A = Z_A - N_A$, the dipole moment
of the molecule is a sum of atomic terms:

$$
\boldsymbol\mu = \sum_A \boldsymbol\mu_A, \qquad
\boldsymbol\mu_A = Z_A\mathbf R_A - \int w_A(\mathbf r)\,\rho(\mathbf r)\,\mathbf r\,d\mathbf r
= \boldsymbol\mu_A^{\rm intr} + q_A\mathbf R_A,
$$

with the intrinsic dipole
$\boldsymbol\mu_A^{\rm intr} = -\int w_A\,\rho\,(\mathbf r - \mathbf R_A)\,d\mathbf r$,
measured from the nucleus of the atom. The atomic dipoles
$\boldsymbol\mu_A$ depend on the origin through $q_A\mathbf R_A$. Keith's
bond charges remove that dependence: each pair of neighbor atoms gets a
charge $Q(A|B) = -Q(B|A)$, the part of the charge of A that comes from B,
with $\sum_B Q(A|B) = q_A$. The charge-transfer dipole of each atom is then

$$
\boldsymbol\mu_A^{\rm CT} = \sum_B Q(A|B)\,\frac{\mathbf R_A - \mathbf R_B}{2},
\qquad
\sum_A \boldsymbol\mu_A^{\rm CT} = \sum_A q_A\mathbf R_A .
$$

APOST-3D generalizes Keith's scheme to overlapping atoms: two atoms are
neighbors when no other atom is closer to the midpoint between them than
half their distance. The charges are renormalized to the number of
electrons before solving for the bond charges, so that they add up to the
charge of the molecule.

When the `.fchk` file comes from a calculation in a static electric field
(see [ENPART](enpart.md#external-electric-field)), `POLAR` also gives the
energies of the intrinsic and charge-transfer dipoles of each atom in the
field, $\boldsymbol\mu_A^{\rm intr}\cdot\mathbf F$ and
$\boldsymbol\mu_A^{\rm CT}\cdot\mathbf F$, and that of the molecular
dipole, $\boldsymbol\mu\cdot\mathbf F$.

## Requirements and input

- A real-space atomic definition (e.g. `TFVC`): with a Hilbert-space
  scheme the run stops.
- A Hartree-Fock or KS-DFT wavefunction: with a correlated wavefunction
  the run stops.
- `POLAR` switches the whole output to 13 decimals (as `FULLPRECISION`),
  for the finite-field derivatives.

```text
# METHOD
TFVC
POLAR
#
```

**Rings.** In a ring of neighbor atoms there are more bond charges than
independent equations, one more per ring. The run then stops, printing
the neighbor matrix, unless the job directory has a file `ichi.inp`:

```text
title line
-1
2 5
1
3
1 2 3
```

- line 2: minus the number of neighbor pairs to remove (here one, the
  pair 2-5, given on the next line), or 0 to remove none;
- then the number of ring equations, and for each ring the number of its
  atoms and the atoms in order around the ring. Each ring equation sets
  the sum of the bond charges around the ring to zero.

The removed pairs and the ring equations together must make up for the
extra bond charges; otherwise the run stops.

## Reading the output

The excerpts come from water, RHF/cc-pVDZ in a field of 0.005 au along
*z*, TFVC atoms (test `H2O-RHF-FIELD`). All dipoles are in au
(e·bohr). `ORIGIN-DEPENDENT ATOMIC DIPOLES` ($\boldsymbol\mu_A$) and
`INTRINSIC ATOMIC DIPOLES` come first, each with its total (the
intrinsic ones here):

% source: tests/reference/H2O-RHF-FIELD.apost
```text
  Atom                      X                   Y                   Z
  -------------------------------------------------------------------
    1  O      0.0000000000000    -0.0000000013242     0.3282710863653
    2  H     -0.1117396137649    -0.0000000000528     0.0929311062892
    3  H      0.1117396137649    -0.0000000000528     0.0929311062892
  -------------------------------------------------------------------
  Total       0.0000000000000    -0.0000000014299     0.5141332989437
```

`BOND CHARGES` gives the number of neighbor pairs and $Q(A|B)$ for each
of them; here each hydrogen gives 0.61 electrons to oxygen:

% source: tests/reference/H2O-RHF-FIELD.apost
```text
  Atoms A - B                    Q(A|B)
  -------------------------------------
    1  O  -   2  H     -0.6085271294953
    1  O  -   3  H     -0.6085271294953
  -------------------------------------
```

`CHARGE-TRANSFER ATOMIC DIPOLES (BOND CHARGES)` follows
($\boldsymbol\mu_A^{\rm CT}$), and `MOLECULAR DIPOLE MOMENT` sums it all
up, with the charge-transfer part as $\sum_A q_A\mathbf R_A$:

% source: tests/reference/H2O-RHF-FIELD.apost
```text
                                      X                   Y                   Z
  -----------------------------------------------------------------------------
  Intrinsic             0.0000000000000    -0.0000000014299     0.5141332989437
  Charge transfer       0.0000000000000     0.0000000000000    -1.3488905582135
  -----------------------------------------------------------------------------
  Total                 0.0000000000000    -0.0000000014299    -0.8347572592698
```

The total is the dipole moment of the molecule. In a field,
`FIELD ENERGY OF THE ATOMIC DIPOLES` gives
$\boldsymbol\mu_A^{\rm intr}\cdot\mathbf F$ and
$\boldsymbol\mu_A^{\rm CT}\cdot\mathbf F$ for each atom, then the energy
of the molecular dipole:

% source: tests/reference/H2O-RHF-FIELD.apost
```text
  Atom              Intrinsic          Charge tr.
  -----------------------------------------------
    1  O      0.0016413554318    -0.0033722272432
    2  H      0.0004646555314    -0.0016861136216
    3  H      0.0004646555314    -0.0016861136216
  -----------------------------------------------
  Total       0.0025706664947    -0.0067444544864

  Field x, y, z (au)                      : 0.000000  0.000000  0.005000
  Field energy of the dipole, mu.F (au)   :      -0.0041737862963
```

The energies follow Gaussian's convention, *E* = +**μ**·**F**. The
intrinsic column is the same as in ENPART's `FIELD INTERACTION PER ATOM`;
the charge-transfer column differs per atom, because ENPART uses the
origin-dependent $q_A\mathbf R_A\cdot\mathbf F$.

## Checking the result

- The total of `MOLECULAR DIPOLE MOMENT` is the dipole moment of the
  wavefunction (`Dipole Moment` in a Gaussian `.fchk`) up to the
  integration error: here −0.8347573 against −0.8347549 au.
- For a neutral molecule, the total charge-transfer dipole of the bond
  charges equals $\sum_A q_A\mathbf R_A$ up to the integration error
  (here 3e-7 au); a difference above 1e-4 au is flagged with a warning
  (integration grid, or `ichi.inp`).
- For an ion, the bond charges cannot carry the charge of the molecule:
  a note says that their charge-transfer dipoles do not add up to
  $\sum_A q_A\mathbf R_A$. Their differences between two fields, and so
  the polarizabilities, are not affected.
- In a field, `Field energy of the dipole` equals ENPART's
  `Field interaction energy`.
