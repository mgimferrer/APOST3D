# APOST-3D
----------

Real-space and Hilbert-space tools for wave function analysis

Computational Chemistry Software developed by P. Salvador's research group from the University of Girona (UdG)

15/09/2023


## Shortcuts
------------

* [Installation](#installation)
* [Documentation](#documentation)
* [Cite the code](#citing-apost3d)
* [Bug reports](#bug-reports-and-feature-requests)


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
-------------

The program has been written by using parts of the program APOST by I. Mayer and A. Hamza, Budapest, 2000-2003

The numerical integration utilizes the subroutines for Lebedev quadrature downloaded from CCL. The appropriate reference is: V.I. Lebedev, and D.N. Laikov "A quadrature formula for the sphere of the 131st algebraic order of accuracy" Doklady Mathematics, 59 477-481 1999

The program makes use of libxc library when necessary, using the F90 interfaces provided by the authors (see http://www.tddft.org/programs/libxc)

We are extremely grateful for the possibility of using these routines!

* Available atomic definitions

Real-space:

        Becke, J. Chem. Phys. 88 2547 1988
        Hirshfeld, Theor. Chim. Acta 44  129 1977
        Hirshfeld-Iterative, J Chem Phys 126 144111 2007
        Topological fuzzy Voronoi cells (TFVC), J Chem Phys, 139 071103 2013
        QTAIM, J. Comput. Chem 30 1082 2009

Hilbert-space:

        Mulliken
        Lowdin
        Davidson-Lowdin
        Natural Atomic Orbitals (NAO)

* Calculating

Atomic and overlap populations, bond orders and valences:

        I. Mayer and P. Salvador, Chem. Phys. Lett. 383 368-375 2004

Hartree-Fock molecular energy decomposition:

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


### Manual
----------

**ONGOING**


* Molecular energy decomposition in the real-space (IQA)

The IQA decomposition requires the keyword ENPART within the # METODE section. This keyword is independent of the WF-type. ENPART keyword needs (mandatory) of the fragments definition (DOFRAGS keyword). Then, the definition of a new section within the .inp file is required, namely # ENPART, including the combination of keywords mentioned below. The keywords to add will depend on the user purpose

Specific keywords for single-determinant wavefunction (Hartree-Fock) decomposition:

        HF --- Decomposition of the Hartree-Fock energy
  
Specific keywords for single-determinant wavefunction (KS-DFT) decomposition:

        LDA --- Decomposition of the LDA energy
        BP86 --- Decomposition of the BP86 energy
        B3LYP --- Decomposition of the B3LYP energy

A general way to define the functional is possible. It requires entering the exchange and correlation (or exchange-correlation) functional ID as provided by the libxc libraries. Specific keywords:

        LIBRARY --- Decomposition of the molecular energy using the combination of exchange and correlation functional desired (if supported by libxc libraries). Mandatory to introduce the ID of the functional using the EXC_FUNCTIONAL, EX_FUNCTIONAL and/or EC_FUNCTIONAL KEYWORDS. For the full list of functionals (and IDs), we guide the user to check libxc libraries
        EXC_FUNCTIONAL --- Exchange-correlation functional ID if functional considered of exchange-correlation type
        EX_FUNCTIONAL --- Exchange functional ID if functional considered of exchange type
        EC_FUNCTIONAL --- Correlation functional ID if functional considered of correlation type

Specific keywords for multi-determinant wavefunction (CASSCF, DMRG) decomposition:

        CASSCF --- Decomposition of the CASSCF energy. Requires of the .dm1 and .dm2 files introduced in the # DM section (see XX)
        CISD --- Decomposition of the CISD energy. Requires of the .dm1 and .dm2 files introduced in the # DM section (see XX)
        CORRELATION --- Decomposition of the exchange and correlation energies, separately. By default the code decompose exchange-correlation altogether

Extra options:

        THREBOD --- Bond order density threshold which decides if evaluating the KS-DFT exchange-correlation between a pair of atoms or not. Expecting an integer number. Value lower than 0 sets the threshold to zero. Default value of 100 (actual threshold = THREBOD/10000.0d0)
        EXACT --- **TODO**
        HOMO --- **TODO**
        DEKIN --- **TODO**
        IONIC --- **TODO**
        TWOELTOLER --- Two-electron integration error threshold to control when the zero-error strategy is invoked. Expects a double-precision number after the keyword. Default value selected as 0.25 (in kcal/mol)
        ANALYTIC --- **TODO**
        MOD-GRIDTWOEL --- Modify the integration grid for the two-electron numerical integrations. Required to define an extra section (# GRID)

Modification of the two-electron integration grid (# GRID section):

        RADIAL --- Number of radial points. Expecting an integer number. Default value of 150
        ANGULAR --- Number of angular points according to the Levedev-Laikov spherical grids. Expecting an integer number. Default value of 590
        rr00 --- Distance value where half of the radial points have been distributed. Expecting a double precision number. Default value of 0.5
        phb1 --- First rotation angle in radians of the second-electron grid points (two-electron zero-error strategy). Expecting a double precision number. Default value of 0.169
        phb2 --- Second rotation angle in radians of the second-electron grid points (two-electron zero-error strategy). Expecting a double precision number. Default value of 0.170

Additional information:

For Gaussian calculations, it is required an extra action to include the analytic energies from the .log file to the .fchk if desired to make use of the two-electron zero-error strategy. In particular:

        INFO: #p must be included in input file from the Gaussian calculation 
        G09 -> $PROG/utils/get_energy mol.log >> mol.fchk
        G16 -> $PROG/utils/get_energy_g16 mol.log >> mol.fchk

For .fchk files obtained from pySCF using the apost3d.py extension, no action is required 

Example input for IQA calculation:

        # METODE
        TFVC
        DOFRAGS
        ENPART
        #
        # ENPART
        LIBRARY
        EX_FUNCTIONAL 106
        EC_FUNCTIONAL 132
        THREBOD -1
        MOD-GRIDTWOEL
        #
        # FRAGMENTS
        2
        1
        1
        -1
        #
        # GRID
        RADIAL 150
        ANGULAR 590
        rr00 0.5
        phb1 0.169
        phb2 0.170
        #

* EDAIQA

**TODO**

* Oxidation states localized orbitals (OSLO)

The calculation of the oxidation states localized orbitals (OSLO) requires to include the OSLO keyword within the # METODE section. It is **mandatory** to include the TFVC keyword in # METODE, independently of the AIM desired to use for the OSLO calculation (technical reasons)

Extracting OSLOs using other AIMs is controlled by definition of a new section (# OSLO) in the .inp file. The AIMs supported to date are: 

        MULLIKEN --- Mulliken AIM
        LOWDIN --- Lowdin AIM
        LOWDIN-DAVIDSON --- Lowdin-Davidson AIM
        NAO-BASIS -- Natural Atomic Orbital (NAO) AIM

The OSLOs can be evaluated from wavefunctions obtained with the Q-Chem software. One simply requires to obtain the .fchk file from Q-Chem and add the QCHEM keyword within # METODE

Example input for OSLO calculation:

        # METODE
        TFVC
        DOFRAGS
        OSLO
        #
        # OSLO
        LOWDIN
        #
        # FRAGMENTS
        2
        1
        1
        -1
        #


## Citing APOST3D
-----------------

The following paper should be cited in publications utilizing the APOST3D program:

P. Salvador, E. Ramos-Cordoba, M. Gimferrer, M. Montilla, Program APOST-3D, Version 4, Girona, 2023


## Bug reports and feature requests
-----------------------------------

Please submit tickets to psalse@gmail.com
