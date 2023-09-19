%nproc=4
%chk=ruco2-pbepbe.chk
%mem=4000mb
#p pbepbe/gen pop=full scf=xqc 5d 7f gfinput pseudo=read iop(3/33=3)

dsfdsf

+2 1
44       0.000000000      0.000000000      0.000000000
6        1.337845834     -1.351865030     -0.005225752
6       -1.344878776      1.344878776      0.000000000
8        2.456508277     -2.482249861     -0.009595353
8       -2.469421932      2.469421932      0.000000000

Ru 0
sdd
****
O C 0
SVP
****

Ru 0
sdd ECP

