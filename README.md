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

The following paper should be cited in publications utilizing `APOST-3D`:

- XX

For atomic and overlap populations, bond orders and valences:

- I. Mayer and P. Salvador, *Chem. Phys. Lett.*, **2004**, 383, 368-375
  DOI: [10.1016/j.cplett.2003.11.048](https://doi.org/10.1016/j.cplett.2003.11.048)

For Hartree-Fock molecular energy decomposition:

        P. Salvador, M. Duran, I.Mayer, J. Chem. Phys. 115 1153-1157 2001
        P. Salvador and I. Mayer, J. Chem. Phys. 120 5046-5052 2004

KS-DFT molecular energy decomposition:

        P. Salvador, I. Mayer, J. Chem. Phys. 126 234113 2007
        M. Gimferrer, P. Salvador, J. Chem. Phys. 158 234105 2023

Molecular energy decomposition for CAS/DMRG wavefunctions:

Effective atomic orbitals:

        I. Mayer, J. Phys. Chem. 100 6249 1996
        I. Mayer and P. Salvador, J. Chem. Phys. 130 234106 2009
        E. Ramos-Cordoba et al., J. Chem. Phys. 138 214107 2013

Local spin analysis:

        E. Ramos-Cordoba et al., J. Chem. Theor. Comput. 8, 1270-1279 2012
        E. Ramos-Cordoba et al., Phys. Chem. Chem. Phys. 14 15291-15298 2012

Effective oxidation states analysis:

        E. Ramos-Cordoba et al., J. Chem. Theor. Comput. 11 1501-1508 2015

Oxidation states from localized orbitals:

        M. Gimferrer et al., J. Chem. Theor. Comput. 18 309-322 2022

Decomposition of EDA quantities into one- and two-center IQA terms:

        M. Gimferrer et al., J. Chem. Theory Comput. 19 3469-3485 2023



## Bug reports and feature requests
-----------------------------------

Please submit tickets on the [issues](https://github.com/mgimferrer/APOST3D/issues) page.

