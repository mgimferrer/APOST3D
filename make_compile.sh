#!/bin/bash -l

# THIS SCRIPT RUN SOME TESTING CALCULATIONS AS FIRST STEP OF THE PGO #

PROGDIR="/home/mgimferrer/APOST3D"
APOST="${PROGDIR}/apost3d"

make -f Makefile_profgen clean
make -f Makefile_profgen
./compiler-runtest
make -f Makefile_profuse clean
make -f Makefile_profuse
./compiler-runtest2

