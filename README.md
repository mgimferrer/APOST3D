# APOST-3D
----------

Real-space and Hilbert-space tools for wave function analysis

Computational Chemistry Software developed by P. Salvador's research group from the University of Girona (UdG)

15/09/2023


## Shortcuts
------------

* [Installation](#installation)
* [Documentation](#documentation)
* [Cite the code](#citations)
* [Bug reports and feature requests](#bug-reports-and-feature-requests)


## Installation
---------------

* Building from source (MG: update to MKL soon)

Prerequisites for manual installation:

        Intel oneAPI Base Toolkit
        Intel® oneAPI HPC Toolkit
        Available free of charge from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit

Program compilation:

        Download source code from Github repository
        Set variable PROG in Makefile to the destination folder
        Modify apost's Makefile to set the variable LIBXCDIR to point to where libxc the library has been installed
        Compilation of libxc libraries (see details provided below) 
        Copy the F90 interfaces provided by libxc (libxc_funcs.f90 and libxc.f90 files on /src folder) on $LIBXCDIR
        cd $PROG ; make

* Compilation Libxc Libraries

After loading the modules of the intel compilers above mentioned, one has to follow the next instructions:

        tar -xzvf libxc-4.2.3.tar.gz
        cd $LIBXCDIR 
        export CC=icc
        export FC=ifort
        export FCFLAGS="-u -fpp1 -nbs -pc80 -pad -align -unroll-aggressive -O3 -ip -no-fp-port -mno-ieee-fp -vec-report0 -no-prec-div -parallel -qopenmp"
        export LDFLAGS="-qopenmp"
        export CFLAGS="-O3"
        ./configure --prefix=$LIBXCDIR
        make
        make install

Important: To date it is not possible to couple APOST3D with newer libxc libraries than the provided, due to internal changes on the modules of libxc. We will work on that as soon as possible


## Documentation
----------------

The `APOST-3D` documentation is hosted [here](DOCUMENTATION.md).


## Citations
------------

### Cite `APOST-3D`

The following paper should be cited in publications utilizing `APOST-3D`:

- XX


### Cite implemented methods

For atomic and overlap populations, bond orders and valences:

- I. Mayer and P. Salvador, *Chem. Phys. Lett.*, **2004**, 383, 368-375
  DOI: [10.1016/j.cplett.2003.11.048](https://doi.org/10.1016/j.cplett.2003.11.048)


For Hartree-Fock molecular energy decomposition:

- P. Salvador, M. Duran and I.Mayer, *J. Chem. Phys.*, **2001**, 115, 1153-1157
  DOI: [10.1063/1.1381407](https://doi.org/10.1063/1.1381407)

- P. Salvador and I. Mayer, *J. Chem. Phys.*, **2004**, 120, 5046-5052
  DOI: [10.1063/1.1646354](https://doi.org/10.1063/1.1646354)


For KS-DFT molecular energy decomposition:

- P. Salvador and I. Mayer, *J. Chem. Phys.*, **2007**, 126, 234113
  DOI: [10.1063/1.2741258](https://doi.org/10.1063/1.2741258)

- M. Gimferrer and P. Salvador, *J. Chem. Phys.*, **2023**, 158, 234105
  DOI: [10.1063/5.0142778](https://doi.org/10.1063/5.0142778)


For CAS/DMRG molecular energy decomposition:


For origin-independent decomposition of static polarizabilities:

- M. Montilla, J. M. Luis and P. Salvador, *J. Chem. Theor. Comput.*, **2021**, 17, 1098-1105
  DOI: [10.1021/acs.jctc.0c00926](https://doi.org/10.1021/acs.jctc.0c00926)


For effective atomic/fragment orbitals:

- I. Mayer, *J. Phys. Chem.*, **1996**, 100, 6249
  DOI: [10.1021/jp952779i](https://doi.org/10.1021/jp952779i)

- I. Mayer and P. Salvador, *J. Chem. Phys.*, **2009**, 130, 234106
  DOI: [10.1063/1.3153482](https://doi.org/10.1063/1.3153482)

- E. Ramos-Cordoba, P. Salvador and I. Mayer, *J. Chem. Phys.*, **2013**, 138, 214107
  DOI: [10.1063/1.4807775](https://doi.org/10.1063/1.4807775)


For local spin analysis:

- E. Ramos-Cordoba, E. Matito, I. Mayer and P. Salvador, *J. Chem. Theor. Comput.*, **2012**, 8, 1270-1279
  DOI: [10.1021/ct300050c](https://doi.org/10.1021/ct300050c)

- E. Ramos-Cordoba, E. Matito, P. Salvador and I. Mayer, *Phys. Chem. Chem. Phys.*, **2012**, 14, 15291-15298
  DOI: [10.1039/C2CP42513K](https://doi.org/10.1039/C2CP42513K)


For effective oxidation states analysis:

- E. Ramos-Cordoba, V. Postils and P. Salvador, *J. Chem. Theor. Comput.*, **2015**, 11, 1501-1508
  DOI: [10.1021/ct501088v](https://doi.org/10.1021/ct501088v)

- M. Gimferrer and P. Salvador, *submitted*, **2024**
  DOI: [XX](XX)


For oxidation states from localized orbitals:

- M. Gimferrer, G. Comas-Vila and P. Salvador, *Molecules*, **2020**, 25, 234
  DOI: [10.3390/molecules25010234](https://doi.org/10.3390/molecules25010234)

- M. Gimferrer, A. Aldossary, P. Salvador and M. Head-Gordon, *J. Chem. Theor. Comput.*, **2022**, 18, 309-322
  DOI: [10.1021/acs.jctc.1c01011](https://doi.org/10.1021/acs.jctc.1c01011)


For decomposition of EDA quantities into one- and two-center IQA terms:

- M. Gimferrer, S. Danes, D. M. Andrada and P. Salvador, *J. Chem. Theory Comput.*, **2023**, 19, 3469-3485
  DOI: [10.1021/acs.jctc.3c00143](https://doi.org/10.1021/acs.jctc.3c00143)


## Bug reports and feature requests
-----------------------------------

Please submit tickets on the [issues](https://github.com/mgimferrer/APOST3D/issues) page.

