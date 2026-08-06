<p align="center"><img width=25.0% src="https://github.com/mgimferrer/APOST3D/blob/master/media/logo-apost.png"></p>

## Chemical concepts from wave function analysis

A Fortran-based code developed at the Universitat de Girona (UdG) by the group of P. Salvador, M. Gimferrer and collaborators.

Builds with **GCC/gfortran** — free, open-source, and available on every Linux
distribution (no Intel compiler or license required).

📖 **Full documentation:** https://apost3d.readthedocs.io

## Shortcuts

* [Installation](#installation)
* [How to use](#how-to-use)
* [Running the test suite](#running-the-test-suite)
* [Documentation](#documentation)
* [Cite the code](#citations)
* [Bug reports and feature requests](#bug-reports-and-feature-requests)

## Installation

Requires GCC/gfortran **10 or newer** (12+ recommended) and `make`.

```bash
# 1. Clone the repository
git clone https://github.com/mgimferrer/APOST3D.git
cd APOST3D

# 2. Set the installation path (add this to your shell profile too)
export APOST3D_PATH=$(pwd)

# 3. Build the bundled libxc-4.2.3 library (once)
bash compile_libxc.sh

# 4. Build apost3d, apost3d-eos, and eos_aom
bash make_compile.sh
```

Both `make_compile.sh` and `make test` accept the same `NTHREADS=<n>` flag for the number of cores to use (e.g. `bash make_compile.sh NTHREADS=4`). Run `bash make_compile.sh help` or `make help` to see all available flags.

Per-distro prerequisite commands, driving `make` directly, and verifying the install are covered in the [Installation](https://apost3d.readthedocs.io/en/latest/installation.html) page of the full documentation.

## How to use

```bash
ulimit -s unlimited          # the code uses large stack-allocated arrays
export OMP_NUM_THREADS=4     # set to the number of cores you want to use

$APOST3D_PATH/apost3d jobname > jobname.apost 2>&1
```

`jobname.fchk` and `jobname.inp` must be present in the working directory. Correlated-wavefunction analyses (CASSCF, DMRG) also need `jobname.dm1`/`.dm2`.

The full input-file format and keyword reference is in the
[Documentation](#documentation).

## Running the test suite

```bash
make test               # build (if needed) + run the entire suite, 1 thread
make test NTHREADS=4    # same, using 4 threads
```

`make test` always runs every case. The active test list, the check format, and how to add a new test case are covered in the [Running the test suite](https://apost3d.readthedocs.io/en/latest/testing.html) page of the full documentation.

## Documentation

Full documentation — installation, usage, the complete input-file keyword reference, worked examples, and troubleshooting — is hosted at **https://apost3d.readthedocs.io**.

## Citations

### Cite the code

The following paper should be cited in publications utilizing `APOST-3D`:

* P. Salvador, E. Ramos-Cordoba, M. Montilla, L. Pujal and M. Gimferrer, *J. Chem. Phys.*, **2024**, 160, 172502
  DOI: [10.1063/5.0206187](https://doi.org/10.1063/5.0206187)

### Cite implemented methods

For atomic and overlap populations, bond orders and valences:

* I. Mayer and P. Salvador, *Chem. Phys. Lett.*, **2004**, 383, 368-375
  DOI: [10.1016/j.cplett.2003.11.048](https://doi.org/10.1016/j.cplett.2003.11.048)

For Hartree-Fock molecular energy decomposition:

* P. Salvador, M. Duran and I.Mayer, *J. Chem. Phys.*, **2001**, 115, 1153-1157
  DOI: [10.1063/1.1381407](https://doi.org/10.1063/1.1381407)
* P. Salvador and I. Mayer, *J. Chem. Phys.*, **2004**, 120, 5046-5052
  DOI: [10.1063/1.1646354](https://doi.org/10.1063/1.1646354)

For KS-DFT molecular energy decomposition:

* P. Salvador and I. Mayer, *J. Chem. Phys.*, **2007**, 126, 234113
  DOI: [10.1063/1.2741258](https://doi.org/10.1063/1.2741258)
* M. Gimferrer and P. Salvador, *J. Chem. Phys.*, **2023**, 158, 234105
  DOI: [10.1063/5.0142778](https://doi.org/10.1063/5.0142778)

For CAS/DMRG molecular energy decomposition:

For effective atomic/fragment orbitals:

* I. Mayer, *J. Phys. Chem.*, **1996**, 100, 6249
  DOI: [10.1021/jp952779i](https://doi.org/10.1021/jp952779i)
* I. Mayer and P. Salvador, *J. Chem. Phys.*, **2009**, 130, 234106
  DOI: [10.1063/1.3153482](https://doi.org/10.1063/1.3153482)
* E. Ramos-Cordoba, P. Salvador and I. Mayer, *J. Chem. Phys.*, **2013**, 138, 214107
  DOI: [10.1063/1.4807775](https://doi.org/10.1063/1.4807775)

For local spin analysis:

* E. Ramos-Cordoba, E. Matito, I. Mayer and P. Salvador, *J. Chem. Theor. Comput.*, **2012**, 8, 1270-1279
  DOI: [10.1021/ct300050c](https://doi.org/10.1021/ct300050c)
* E. Ramos-Cordoba, E. Matito, P. Salvador and I. Mayer, *Phys. Chem. Chem. Phys.*, **2012**, 14, 15291-15298
  DOI: [10.1039/C2CP42513K](https://doi.org/10.1039/C2CP42513K)

For effective oxidation states analysis:

* E. Ramos-Cordoba, V. Postils and P. Salvador, *J. Chem. Theor. Comput.*, **2015**, 11, 1501-1508
  DOI: [10.1021/ct501088v](https://doi.org/10.1021/ct501088v)
* M. Gimferrer and P. Salvador, _submitted_, **2024**
  DOI: [XX](XX)

For oxidation states from localized orbitals:

* M. Gimferrer, G. Comas-Vila and P. Salvador, *Molecules*, **2020**, 25, 234
  DOI: [10.3390/molecules25010234](https://doi.org/10.3390/molecules25010234)
* M. Gimferrer, J. Van der Mynsbrugge, A. T. Bell, P. Salvador and M. Head-Gordon *Inorg. Chem.*, **2020**, 59, 15410-15420
  DOI: [10.1021/acs.inorgchem.0c02405](https://doi.org/10.1021/acs.inorgchem.0c02405)
* M. Gimferrer, A. Aldossary, P. Salvador and M. Head-Gordon, *J. Chem. Theor. Comput.*, **2022**, 18, 309-322
  DOI: [10.1021/acs.jctc.1c01011](https://doi.org/10.1021/acs.jctc.1c01011)

For decomposition of EDA quantities into one- and two-center IQA terms:

* M. Gimferrer, S. Danes, D. M. Andrada and P. Salvador, *J. Chem. Theory Comput.*, **2023**, 19, 3469-3485
  DOI: [10.1021/acs.jctc.3c00143](https://doi.org/10.1021/acs.jctc.3c00143)

For origin-independent decomposition of static polarizabilities:

* M. Montilla, J. M. Luis and P. Salvador, *J. Chem. Theor. Comput.*, **2021**, 17, 1098-1105
  DOI: [10.1021/acs.jctc.0c00926](https://doi.org/10.1021/acs.jctc.0c00926)

## Bug reports and feature requests

Please submit tickets on the [issues](https://github.com/mgimferrer/APOST3D/issues) page, and/or send an email to mgimferrer18@gmail.com and pedro.salvador@udg.edu
