#!/bin/bash -l

## EXECUTING THIS SCRIPT ONE WILL COMPILE THE APOST3D CODE IN HIGH-PERFORMANCE MODE ##
## THIS INVOLVES TWO STEPS: 1) GENERATION OF THE .dyn FILES (PROFILING) AND 2) USE  ##
## OF THEM IN ORDER TO COMPILE WITH PROFILE GUIDED OPTIMIZATION (PGO)               ##

PROGDIR="/home/mgimferrer/APOST3D"
APOST="${PROGDIR}/apost3d"

## STEP 1 ##
make -f Makefile_profgen clean
make -f Makefile_profgen
./compiler-runtest

## STEP 2 ##
make -f Makefile_profuse clean
make -f Makefile_profuse
./compiler-runtest2

