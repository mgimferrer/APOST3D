#!/bin/bash
#Test set calculations
#reference="/scratch/SOFT/APOST3D.3.1/apost3d"
reference="/home/psalvador/APOST3D.3.1.devel/apost3d"
for mol in n2-cas o2_cas feco2-pbepbe h2o-t-b3lyp ruco2-pbepbe h2o-s-hf ; do
echo "Doing for " $mol
$reference $mol >${mol}.apost_devel
done

