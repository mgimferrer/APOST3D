#!/bin/bash -l

## EXECUTING THIS SCRIPT ONE WILL COMPILE THE LIBXC LIBRARIES REQUIRED           ##
## WITH THE FLAGS, AND COPY OF FILES REQUIRED FOR LINKING AFTERWARDS TO APOST-3D ##

LIBXCDIR="${APOST3D_PATH}/libxc-4.2.3"

tar -xzvf libxc-4.2.3.tar.gz
cd $LIBXCDIR

## EXPORTS REQUIRED FOR CREATING PROPERLY THE MAKEFILES ##
export CC=icx
export FC=ifort
export FCFLAGS="-u -fpp1 -nbs -pc80 -pad -align -unroll-aggressive -O3 -ip -no-fp-port -mno-ieee-fp -vec-report0 -no-prec-div -parallel -qopenmp"
export LDFLAGS="-qopenmp"
export CFLAGS="-O3"

## CREATING MAKEFILES AND COMPILING ##
./configure --prefix="${LIBXCDIR}"
make
make install

## COPYING F90 INTERFACES REQUIRED TO LIBXCDIR FOLDER ##
cp ./src/libxc_funcs.f90 . 
cp ./src/libxc.f90 .
