# APOST-3D

Real-space and Hilbert-space tools for wave function analysis

Computational Chemistry Software developed by P. Salvador's research group from the University of Girona (UdG)

15/09/2023
----------

* [Documentation](#documentation)
* [Installation](#installation)


Installation
------------

* Building from source

Prerequisites for manual installation:

        XXX

Program compilation:

        Download source code from Github repository
        Set variable PROG in Makefile to the destination folder
        Modify apost's Makefile to set the variable LIBXCDIR to point to where libxc the library has been installed
        Compilation details of the libxc libraries for APOST3D is provided [here](#compilation-libxc-libraries) 
        Copy the F90 interfaces provided by libxc (libxc_funcs.f90 and libxc.f90 files on /src folder) on $LIBXCDIR
        cd $PROG ; make

Compilation Libxc Libraries
---------------------------

TO DO 

Important: To date it is not possible to couple APOST3D with newer libxc libraries than the provided, due to internal changes on the modules of libxc. We will work on that as soon as possible.

Citing APOST3D
--------------

The following paper should be cited in publications utilizing the APOST3D program:

TO DO

Documentation
-------------

TO DO

Bug reports and feature requests
--------------------------------
Please submit tickets to XX.