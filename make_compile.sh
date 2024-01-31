#!/bin/bash -l

## EXECUTING THIS SCRIPT ONE WILL COMPILE THE APOST3D CODE IN HIGH-PERFORMANCE MODE ##
## THIS INVOLVES TWO STEPS: 1) GENERATION OF THE .dyn FILES (PROFILING) AND 2) USE  ##
## OF THEM IN ORDER TO COMPILE WITH PROFILE GUIDED OPTIMIZATION (PGO)               ##

## EXPORTS REQUIRED FOR PARALLELIZATION PURPOSES ##
export OMP_NUM_THREADS=48 ## NUMBER OF CORES, SET TO THE MAXIMUM OF THE NODE (RECOMMENDATION) ##
export KMP_STACKSIZE=100m
ulimit -s unlimited

## STEP 1 ##
make -f Makefile_profgen clean
make -f Makefile_profgen
./compiler-runtest

## STEP 2 ##
make -f Makefile_profuse clean
make -f Makefile_profuse
./compiler-runtest2

