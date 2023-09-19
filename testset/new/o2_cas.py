from pyscf import gto, scf, tools, mcscf, fci
import math
import numpy
import scipy.linalg
import apost3d as apost

##main program##

molname = 'o2_cas'

mol = gto.M(atom = [
            [8,( 0, 0, 0)],
	    [8, (0, 0, 1.20)]
                   ], basis = 'def2svp')
mol.spin = 2 # number of unpaired electrons
mol.charge = 0
mol.cart= True
mol.max_memory = 4000

# start with refernce RHF calc.
myhf = scf.RHF(mol)
myhf.kernel()
#then full-valence CAS with 6 orbitals, 6 electrons
mycas = mcscf.CASSCF(myhf, 8,6)
mycas.fix_spin(shift=0.2,ss=2)
mycas.kernel()
print(mcscf.spin_square(mycas))

# Generate FCHK and RDMs
apost.write_fchk(mol, mycas, molname, myhf.get_ovlp())
apost.write_dm(mol, mycas, molname) 

