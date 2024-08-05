#!/bin/bash -l

## THIS SCRIPT WILL COMPILE APOST3D CODE IN HIGH-PERFORMANCE MODE IN TWO STEPS      ##
## OF THEM IN ORDER TO COMPILE WITH PROFILE GUIDED OPTIMIZATION (PGO)               ##

## EXPORTS REQUIRED FOR PARALLELIZATION PURPOSES ##
## NUMBER OF CORES, SET TO THE MAXIMUM OF THE NODE (RECOMMENDATION) ##
export NUM_THREADS=$(grep -c ^processor /proc/cpuinfo)
export KMP_STACKSIZE=100m
ulimit -s unlimited

## STEP 1 ##
## COMPILING WITH PROFILING. THE TESTS WILL GENERATE .dyn FILES ##
make -f Makefile_profgen clean
make -f Makefile_profgen
## TESTS WILL RUN ON SINGLE-PROCESSOR FOR ABOUT 10 MIN ##
./compiler-runtest

## STEP 2 ##
## RECOMPILING USING PROFILING INFO GENERATED IN PREVIOUS STEP  ##
make -f Makefile_profuse clean
make -f Makefile_profuse
# TESTS WILL NOW RUN IN PARALLEL FOR COMPARISON  ##
./compiler-runtest2

## COMPILING THE UTILS ##
make -f Makefile_profuse util
