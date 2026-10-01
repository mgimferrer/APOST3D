# SPIN - Local spin analysis

`SPIN` decomposes the expectation value of the total spin into
**atomic** (local spin) and **diatomic** (spin coupling) contributions:

$$
\langle \hat S^2\rangle = \sum_A \langle \hat S^2\rangle_A + \sum_{A}\sum_{B\neq A} \langle \hat S^2\rangle_{AB}
$$

A local spin $\langle \hat S^2\rangle_A$ of 3/4 means one unpaired electron
localized on atom (or fragment) A, a value of 2 means two parallel unpaired
electrons, and so on. The sign of $\langle \hat S^2\rangle_{AB}$ says how
the local spins of A and B couple: positive for ferromagnetic, negative for
antiferromagnetic coupling. Unlike a spin population, this tells apart the
electron pairing of an ordinary covalent bond (no local spins) from two
antiferromagnetically coupled spins (local spins that cancel through a
negative coupling term). This is what makes the analysis useful for
diradicals and other singlet states described by correlated
wavefunctions.

```{admonition} Cite
:class: note

E. Ramos-Cordoba, E. Matito, I. Mayer and P. Salvador, *J. Chem. Theory
Comput.*, **2012**, 8, 1270-1279; E. Ramos-Cordoba, E. Matito,
P. Salvador and I. Mayer, *Phys. Chem. Chem. Phys.*, **2012**, 14,
15291-15298. See [Citations](../citations.md).
```

## How it works

The decomposition follows Mayer's approach: $\langle \hat S^2\rangle$ is
first written in terms of the first- and second-order density matrices,
and the resulting one- and two-electron integrals are split into one- and
two-center terms with the atomic weight functions of the chosen
[AIM scheme](../guide/aim.md). Several such splittings satisfy the basic
physical requirements (zero local spins for a closed-shell restricted
determinant, the free-atom values at dissociation). APOST-3D uses the one
singled out in the paper above (parameter $a = 3/4$, printed in the box
titles): it is the only one that gives the correct local spin, 3/4 of the
electron population, for a single electron, and it keeps local spins small
in ordinary closed-shell molecules described with correlated wavefunctions.

- **Single determinant** (UHF, UKS): everything follows from the spin
  density matrix $\mathbf P^s = \mathbf P^\alpha - \mathbf P^\beta$ and the
  atomic overlap matrices $\mathbf S^A$:

  $$
  \langle \hat S^2\rangle_A = \tfrac34\,\mathrm{Tr}(\mathbf P^s\mathbf S\,\mathbf P^s\mathbf S^A)
  - \tfrac14\,\mathrm{Tr}(\mathbf P^s\mathbf S^A\mathbf P^s\mathbf S^A)
  + \tfrac14\,[\mathrm{Tr}(\mathbf P^s\mathbf S^A)]^2
  $$

  $$
  \langle \hat S^2\rangle_{AB} = -\tfrac14\,\mathrm{Tr}(\mathbf P^s\mathbf S^A\mathbf P^s\mathbf S^B)
  + \tfrac14\,\mathrm{Tr}(\mathbf P^s\mathbf S^A)\,\mathrm{Tr}(\mathbf P^s\mathbf S^B)
  $$

  For a restricted closed-shell determinant all terms are zero, and the
  analysis is skipped.
- **Correlated wavefunctions** (CASSCF, CISD, ...): the local spin contains
  the density of effectively unpaired electrons $u(\mathbf r)$ (Takatsuka's
  definition, built from the natural orbitals as $\sum_i n_i(2-n_i)|\phi_i|^2$)
  and the cumulant of the second-order density matrix. This needs the
  1- and 2-RDMs, given in a [`# DM` block](../input/dm.md).

## Requirements and input

- Any atomic definition. The tests use `TFVC`; Hilbert-space schemes
  (`MULLIKEN`, `LOWDIN`, ...) are accepted too.
- An open-shell single determinant, or a correlated wavefunction with its
  1- and 2-RDMs (`DM 2` and a `# DM` block). With `DM 1` only, the run
  stops with `Local Spin needs dm1 and dm2 for correlated WFs`.
- Fragments (`DOFRAGS`) are optional: the matrices are then also summed
  over the atoms of each fragment.

Single determinant:

```text
# METHOD
TFVC
SPIN
#
```

Correlated wavefunction (here from pySCF, see
[Preparing the wavefunction](../guide/wavefunctions.md)):

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

`DM 2` switches the local spin analysis on by itself, so a run that
reads the 2-RDM for another analysis (e.g. ENPART) also prints it.

For a restricted closed-shell single determinant, `SPIN` only prints
`No Local Spin Analysis needed for Restricted SD WFs` and is otherwise
ignored.

## Reading the output

**Single determinant.** The excerpts come from water in its triplet state,
UKS B3LYP, TFVC (test `H2O-T-B3LYP`). The spin populations are printed
earlier, with the [population analysis](population.md). The local spin
section starts with the effectively unpaired electrons of each atom:

```text
  ----------------------------------
    EFFECTIVELY UNPAIRED ELECTRONS
  ----------------------------------

    Atom         u_A
  ------------------
   1  O     1.287250
   2  H     0.361616
   3  H     0.361616
  ------------------
  Sum check N_D =    2.01048
```

followed by the decomposition. Diagonal elements are the local spins
$\langle \hat S^2\rangle_A$, off-diagonal ones the coupling terms
$\langle \hat S^2\rangle_{AB}$ (each pair appears twice, as AB and BA):

```text
  -------------------------------------------
    "FUZZY ATOMS" S^2 DECOMPOSITION (a=3/4)
  -------------------------------------------

              1  O        2  H        3  H

   1  O     1.079130    0.091036    0.091036
   2  H     0.091036    0.272245    0.008738
   3  H     0.091036    0.008738    0.272245

  Sum check <S^2> =    2.00524
```

The two unpaired electrons are mostly on oxygen (1.08, compared with 0.75
for one and 2 for two electrons localized on it), partly delocalized onto
the hydrogens, and all couplings are positive (ferromagnetic), as expected for
a triplet. The sum check is the total $\langle \hat S^2\rangle$ of the
wavefunction: 2.00524 here, slightly above the exact 2 because of the spin
contamination of the UKS determinant (it matches the `S**2` value of the
`.fchk` file).

A second matrix, `"FUZZY ATOMS" DAVIDSON SPIN DEC. MATRIX`, follows. It is
under revision: use the $a=3/4$ matrix above.

**Correlated wavefunction.** For LiH at 3.5 Å, CASSCF(2,2) (test
`LiH-35-CAS22`), the section first prints bond order and
localization/delocalization index matrices built from the 2-RDM, then:

```text
    Atom         u_A
  ------------------
   1 Li     0.607882
   2  H     0.793956
  ------------------
  Sum check N_D =    1.40184

  -------------------------------------
    APOST3D S^2 DECOMPOSITION (a=3/4)
  -------------------------------------

              1 Li        2  H

   1 Li     0.452573   -0.452573
   2  H    -0.452573    0.452573

  Sum check <S^2> =   -0.00000
```

The total is zero (a singlet), yet each atom carries a sizeable local spin
that cancels through a negative coupling: the stretched bond has strong
diradical character, two antiferromagnetically coupled electrons on its
way to two free atoms (where each local spin would be 3/4). At the
equilibrium distance the same analysis gives small local spins.

With `DOFRAGS`, `FRAGMENT ANALYSIS : Local Spin Analysis` and
`FRAGMENT ANALYSIS : Num. eff. unpaired elec.` give the same quantities
summed over the atoms of each fragment.

## Checking the result

- The sum check must equal the total $\langle \hat S^2\rangle$ of the
  wavefunction (`S**2` in a Gaussian `.fchk`; 0 for a correlated singlet).
  A clear difference points to an integration problem: use a larger grid
  (see [Running a calculation](../guide/running.md)).
- For a correlated wavefunction, check that the 1- and 2-RDMs correspond to
  the orbitals in the `.fchk` file (written together by `apost3d.py`).
