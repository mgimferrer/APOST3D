# APOST-3D

**Chemical concepts from wave function analysis**

APOST-3D reads a converged wavefunction (a formatted checkpoint file from
Gaussian, Q-Chem, pySCF, ...) and translates it into chemical language:
atomic charges and bond orders, effective atomic and fragment orbitals,
oxidation states, local spins, and a decomposition of the molecular
energy into atomic and interatomic terms. It is developed at the
Universitat de Girona by P. Salvador, M. Gimferrer and collaborators, and
is free and open source ([GitHub](https://github.com/mgimferrer/APOST3D)).

## What can I compute?

Every run starts by choosing an **atoms-in-molecules (AIM) scheme**, the
rule that splits the molecule into atoms ([Atoms in molecules](guide/aim.md)).
On top of it, one or more analyses can be requested in the same run:

| Analysis | Keyword | What you get |
|---|---|---|
| [Population analysis](methods/population.md) | *(always)* | Atomic charges, spin populations, bond orders and valences |
| [Effective atomic orbitals](methods/effao.md) | `EFFAO`, `UEFFAO`, `EFFAO-U` | The orbitals and occupations that describe each atom or fragment in the molecule |
| [Effective oxidation states](methods/eos.md) | `EOS` | Oxidation states of atoms or fragments for any wavefunction, with a reliability index |
| [Generalized EOS](methods/geos.md) | `GEOS` | Oxidation states from the paired and unpaired densities, for open-shell and correlated wavefunctions |
| [Oxidation states from localized orbitals](methods/oslo.md) | `OSLO` | Fragment-localized orbitals (OSLOs) and the oxidation states they imply |
| [Local spin](methods/spin.md) | `SPIN` | Atomic and diatomic contributions to ⟨*S*²⟩ |
| [Energy partitioning](methods/enpart.md) | `ENPART` | The molecular energy split into atomic self-energies and interatomic interactions (IQA), for HF, KS-DFT and CASSCF/CISD |

Atoms can be grouped into [fragments](guide/fragments.md) (ligands, metal
centers, molecules) for any of these analyses.

## Where to start

1. [Install](installation.md) the program: one command builds everything.
2. Follow the [tutorial](tutorial.md): a first calculation on a small
   molecule, from the input file to the results.
3. Look up what you need: the [user guide](guide/wavefunctions.md) (how to
   prepare and run a calculation), one page per
   [analysis method](methods/population.md), and the
   [input reference](input/index.md) for every keyword.

```{admonition} Cite this work
:class: note

If you use APOST-3D in published work, please cite the program paper and
the papers of the methods you used; see [Citations](citations.md).
```

```{toctree}
:hidden:
:caption: Getting started

installation
tutorial
```

```{toctree}
:hidden:
:caption: User guide

guide/wavefunctions
guide/running
guide/aim
guide/fragments
guide/output
guide/visualization
```

```{toctree}
:hidden:
:caption: Analysis methods

methods/population
methods/effao
methods/eos
methods/geos
methods/oslo
methods/spin
methods/enpart
```

```{toctree}
:hidden:
:caption: Input reference

input/index
input/method
input/enpart
input/oslo
input/cube
input/grid
input/dm
input/fragments
input/examples
```

```{toctree}
:hidden:
:caption: Tools

tools/utilities
```

```{toctree}
:hidden:
:caption: Help and reference

troubleshooting
glossary
citations
changelog
```

```{toctree}
:hidden:
:caption: Developer guide

developer/building
developer/testing
```

## Acknowledgements

The program uses parts of the program APOST by I. Mayer and A. Hamza
(Budapest, 2000-2003). The numerical integration uses the Lebedev-Laikov
angular quadrature routines
([source](http://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F);
V. I. Lebedev and D. N. Laikov, *Doklady Mathematics*, **1999**, 59,
477-481), and exchange-correlation functionals come from the
[libxc](https://libxc.gitlab.io/) library through its Fortran 2003
interface. We are grateful for the possibility of using these routines,
and thank R. Oswald for technical support with the parallelization and
the build setup.

## Bug reports and questions

Open an issue on [GitHub](https://github.com/mgimferrer/APOST3D/issues),
or write to mgimferrer18@gmail.com or pedro.salvador@udg.edu.
