# APOST-3D

**Chemical concepts from wave function analysis**

APOST-3D is a Fortran-based code developed at the Universitat de Girona (UdG)
by P. Salvador and collaborators. It takes a converged wavefunction (from
Gaussian, Q-Chem, pySCF, or any other source that can produce a formatted
checkpoint file) and extracts chemically meaningful, real-space or
basis-set-based descriptors from it: how the electron density and orbitals
are shared out among atoms and fragments, what that implies about bonding,
oxidation states, and spin, and how the molecular energy itself decomposes
into atomic and interatomic contributions.

It builds with **GCC/gfortran** — free, open-source, and available on every
Linux distribution, with no Intel compiler or license required.

```{admonition} Source code
:class: tip

The source and issue tracker live on GitHub:
[github.com/mgimferrer/APOST3D](https://github.com/mgimferrer/APOST3D)
```

## What it does

Every calculation starts by choosing an **atom-in-a-molecule (AIM) scheme**
— the rule used to partition the molecule into atomic contributions. APOST-3D
supports two families: **real-space** schemes that partition 3D space via
numerical integration (`TFVC`, `HIRSH`, `HIRSH-IT`, the deprecated
`BECKE-RHO`), and **Hilbert-space** schemes that partition the basis-set
overlap instead (`MULLIKEN`, `LOWDIN`, `LOWDIN-DAVIDSON`, `NAO-BASIS`). See
the [AIM definitions table](input/method.md) for the full list.

On top of whichever AIM scheme is selected, one or more **analysis tools**
can be requested in the same run:

- **EFFAO / UEFFAO / EFFAO-U** — effective atomic and fragment orbitals from
  the total, spin-resolved, or paired/unpaired electron density.
- **EOS / EOS-U** — effective oxidation states, from the regular or the
  paired/unpaired density functions.
- **OSLO / OS-CENTROID** — oxidation states from localized orbitals, or from
  the centroids of localized orbitals.
- **SPIN** — local spin analysis (LSA), including for correlated
  wavefunctions via 1-/2-RDM input.
- **ENPART / EDAIQA** — real-space molecular energy decomposition (IQA) for
  HF, DFT, and correlated (CASSCF/CISD) wavefunctions, and decomposition of
  Energy Decomposition Analysis (EDA) terms.
- **POLAR** — Bader-Keith decomposition of the molecular dipole moment.

Fragments of atoms (rather than individual atoms) can be defined for any of
the above via `DOFRAGS`, and cube files of the resulting orbitals can be
generated directly. The full keyword-by-keyword reference is in
[Input reference](input/index.md), with five worked examples in
[Input examples](input/examples.md).

```{admonition} Cite this work
:class: note

If you use APOST-3D in published work, please cite the code and the papers
behind the specific analysis tools you used — see [Citations](citations.md)
for the full reference list.
```

## Getting started

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
quickstart
testing
troubleshooting
```

## Input reference

```{toctree}
:maxdepth: 2
:caption: Input reference

input/index
input/method
input/oslo
input/enpart
input/other-blocks
input/examples
input/pyscf
```

## Citations

```{toctree}
:maxdepth: 1
:caption: Citations

citations
```

## Acknowledgements

The program has been written using parts of the program APOST by I. Mayer
and A. Hamza, Budapest, 2000-2003.

The numerical integration utilizes the Lebedev quadrature subroutines
[available here](http://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/Lebedev-Laikov.F).
The appropriate reference is: V.I. Lebedev and D.N. Laikov, "A quadrature
formula for the sphere of the 131st algebraic order of accuracy," *Doklady
Mathematics*, 59, 477-481 (1999).

The program makes use of the `Libxc` library when necessary, using the F90
interfaces provided by the authors (see https://libxc.gitlab.io/).

We are extremely grateful for the possibility of using these routines. We
also acknowledge R. Oswald for technical support in code parallelization and
compilation setup preparation.

## Bug reports and feature requests

Please submit tickets on the [issues page](https://github.com/mgimferrer/APOST3D/issues),
or send an email to mgimferrer18@gmail.com or pedro.salvador@udg.edu.
