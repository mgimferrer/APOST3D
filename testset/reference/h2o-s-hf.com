%nproc=1
%chk=h2o-s-hf.chk
%mem=2000mb
#p HF/cc-pvtz pop=full field=z+016 gfinput iop(3/33=1)

dsfdsf

0 1
O
H 1 r
H 1 r 2 a

r 0.940
a 105.0

