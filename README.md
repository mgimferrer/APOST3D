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

Documentation
-------------

A manual of the program, together with examples for each type of calculation can be found here(TO DO: readthedocs?)

Citing APOST3D
--------------

The following paper should be cited in publications utilizing the APOST3D program:

P. Salvador, E. Ramos-Cordoba, M. Gimferrer, M. Montilla, Program APOST-3D, Version 4, Girona, 2023

Bug reports and feature requests
--------------------------------

Please submit tickets to psalse@gmail.com