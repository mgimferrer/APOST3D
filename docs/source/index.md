# APOST-3D

**Chemical concepts from wave function analysis**

APOST-3D is a Fortran-based code developed at the Universitat de Girona (UdG)
by P. Salvador and collaborators. It computes atoms-in-molecules partitions,
population and bond-order analyses, effective atomic/fragment orbitals,
oxidation-state assignments, energy decompositions, and related
wave-function-based descriptors from a Gaussian-type formatted checkpoint
file.

It builds with **GCC/gfortran** — free, open-source, and available on every
Linux distribution, with no Intel compiler or license required.

```{admonition} Source code
:class: tip

The source and issue tracker live on GitHub:
[github.com/mgimferrer/APOST3D](https://github.com/mgimferrer/APOST3D)
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

## About

```{toctree}
:maxdepth: 1
:caption: About

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
interfaces provided by the authors (see http://www.tddft.org/programs/libxc).

We are extremely grateful for the possibility of using these routines. We
also acknowledge R. Oswald for technical support in code parallelization and
compilation setup preparation.

## Bug reports and feature requests

Please submit tickets on the [issues page](https://github.com/mgimferrer/APOST3D/issues),
or send an email to mgimferrer18@gmail.com or pedro.salvador@udg.edu.
