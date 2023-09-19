#!/bin/bash
# will write until Alpha Orbital but excluding the header, then insert $2
# and then write from Total SCF until the end
#
# assumes correct ordering of the sections in fchk... those coming from Qchem do not follow!!
#
sed -n '1,/Alpha Orbi/p' $1 |sed '$d'
cat efo_occ.dat
cat efo_coeff.dat
#cat $2
sed -n '/Total SCF/,$p' $1
