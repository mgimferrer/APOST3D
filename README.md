# APOST-3D

Real-space and Hilbert-space tools for wave function analysis

Computational Chemistry Software developed by P. Salvador's research group from the University of Girona (UdG)

15/09/2023
----------

* [Installation](#installation)
* [Documentation](#documentation)

Installation
------------

* Building from source

Prerequisites for manual installation:

        1) Intel oneAPI Base Toolkit
        2) Intel® oneAPI HPC Toolkit
        Available free of charge from: https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit

Program compilation:

        1) Download source code from Github repository
        2) Set variable PROG in Makefile to the destination folder
        3) Modify apost's Makefile to set the variable LIBXCDIR to point to where libxc the library has been installed
        4) Compilation of libxc libraries (see details provided below) 
        5) Copy the F90 interfaces provided by libxc (libxc_funcs.f90 and libxc.f90 files on /src folder) on $LIBXCDIR
        6) cd $PROG
        7) make

Compilation Libxc Libraries:

        After loading the modules of the intel compilers above mentioned, one has to follow the next instructions:
        1) cd $LIBXCDIR 
        2) export CC=icc
        3) export FC=ifort
        4) export FCFLAGS="-u -fpp1 -nbs -pc80 -pad -align -unroll-aggressive -O3 -ip -no-fp-port -mno-ieee-fp -vec-report0 -no-prec-div -parallel -qopenmp"
        5) export LDFLAGS="-qopenmp"
        6) export CFLAGS="-O3"
        7) ./configure --prefix=$LIBXCDIR
        8) make
        9) make install

Important: To date it is not possible to couple APOST3D with newer libxc libraries than the provided, due to internal changes on the modules of libxc. We will work on that as soon as possible.

Documentation
-------------

TO DO

Citing APOST3D
--------------

The following paper should be cited in publications utilizing the APOST3D program:

TO DO

Bug reports and feature requests
--------------------------------

Please submit tickets to XX.