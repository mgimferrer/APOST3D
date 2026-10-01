# Glossary

```{glossary}
AIM scheme
  Atoms-in-molecules scheme: the rule that decides which part of the
  wavefunction belongs to each atom. See [Atoms in molecules](guide/aim.md).

Real-space scheme
  An AIM scheme that splits the three-dimensional space with atomic
  weight functions (TFVC, Hirshfeld, Becke). Quantities are obtained by
  numerical integration.

Hilbert-space scheme
  An AIM scheme that assigns basis functions to atoms (Mulliken, Löwdin,
  NAO). No integration grid is needed.

Weight function
  In a real-space scheme, the function $w_A(\mathbf r)$ (between 0 and 1,
  summing to 1 over all atoms) that tells how much of each point of space
  belongs to atom A. For a fragment, the sum of the weight functions of its
  atoms.

TFVC
  Topological fuzzy Voronoi cells: the recommended real-space scheme,
  close to Bader's QTAIM atoms at a fraction of the cost.

Fragment
  A group of atoms treated as a unit (`DOFRAGS`). See [Fragments](guide/fragments.md).

Atomic overlap matrix
  The overlap of the orbitals integrated over one atom, $S^A_{ij}$. Every
  analysis of APOST-3D is built from these matrices.

Effective atomic orbitals
EFFAO
  The orbitals, with their occupations, that describe the part of the
  density belonging to one atom. See [EFFAO](methods/effao.md).

EFO
Effective fragment orbitals
  The effective orbitals of a fragment.

Net occupation
  The occupation of an effective orbital within its own fragment
  (eigenvalue of the fragment's density matrix).

Gross occupation
  The occupation of an effective orbital including its share of the
  overlap with other fragments. Gross occupations add up to the number of
  electrons; they are the ones used by EOS.

EOS
  Effective oxidation states: oxidation states from assigning the
  electrons to the most occupied EFOs of all fragments. See [EOS](methods/eos.md).

GEOS
  Generalized effective oxidation states, from the paired and unpaired
  densities. See [GEOS](methods/geos.md).

Reliability index
R(%)
  How clear-cut an EOS or GEOS assignment is, from the occupation gap
  between the last occupied and the first unoccupied EFO: 100 for a gap of
  half an electron or more, 50 for two equally plausible assignments.

Paired density
Unpaired density
  The total density split into the part of electrons that form pairs and
  the part of effectively unpaired electrons (Takatsuka's definition,
  from the natural occupations $n_i$ as $n_i(2-n_i)$). Used by GEOS and
  EFFAO-U.

Natural orbitals
  The eigenvectors of the one-particle density matrix; their eigenvalues
  are the natural occupations (0 or 2 for a closed-shell determinant,
  fractional for correlated wavefunctions).

RDM
  Reduced density matrix. The 1-RDM gives the density; the 2-RDM is needed
  for local spins and energy partitioning of correlated wavefunctions.

OSLO
  Oxidation state localized orbitals: orbitals localized onto fragments,
  one at a time, used to assign oxidation states. See [OSLO](methods/oslo.md).

FOLI
  Fragment orbital localization index: 1 for an orbital entirely on one
  fragment, 2 for one shared equally by two fragments.

Δ-FOLI
  In OSLO, the FOLI gap between the selected orbital and the next
  candidate; in the last iteration, it measures how clear the assignment
  is.

Local spin
  The atomic contribution $\langle \hat S^2\rangle_A$ to the total spin;
  3/4 for one unpaired electron on the atom. See [SPIN](methods/spin.md).

IQA
  Interacting quantum atoms: the decomposition of the energy into atomic
  self-energies and interatomic interaction energies, computed by
  [ENPART](methods/enpart.md).

BODEN
  Bond order density: the density of an atom pair from which ENPART
  computes the interatomic exchange-correlation energy in KS-DFT.

Zero-error strategy
  ENPART's correction of the numerical error of the two-electron
  integrations, using the exact two-electron energy of the calculation.

Multipolar expansion
  The approximation used by ENPART for the exchange-correlation energy of
  weakly bonded atom pairs (bond order below `THREBOD`).
```
