#!/bin/bash
# output is generated on screen. must be redirected to *.fchk file for visualization
# will write until Alpha Orbital but excluding the header, then insert efo_occ.dat and efo_coef.dat 
# and then write from Total SCF until the end
#
# assumes correct ordering of the sections in fchk... those coming from Qchem do not follow!!
#
sed -n '1,/Alpha Orbi/p' $1 |sed '$d'
cat efo_occ.dat efo_coeff.dat
sed -n '/Total SCF/,$p' $1
