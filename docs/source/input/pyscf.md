# Extracting .fchk, .dm1 and .dm2 files from pySCF

The `apost3d.py` file provided in the `$APOST3D_PATH/utils` folder collects
several functions for creating `.fchk`, `.dm1` and `.dm2` files from a
pySCF run. The location of the file must be added to the `PYTHONPATH`
variable, for instance by doing:

```bash
export PYTHONPATH=$PYTHONPATH:$APOST3D_PATH/utils
```

As an example, here is an input for a CASSCF(2,2) calculation on the HF
molecular system in the singlet spin state.

```python
from pyscf import gto, scf, mcscf
import apost3d as apost

molname = 'HF-CASSCF'
mol=gto.M()
mol.atom='''
H    0.0000000    0.0000000    0.0000000
F    0.0000000    0.0000000    0.9500000
'''
mol.basis='def2tzvp'
mol.spin = 0
mol.charge = 0
mol.cart= False
mol.symmetry = False
mol.verbose=4
mol.max_memory=40000
mol.build()

# Calculate RHF reference (guess) #
mf = scf.RHF(mol)
mf.kernel()

# Once our guess (converged HF) is created, time to create and converge the CASSCF  #
mycas = mcscf.CASSCF(mf,2,2) # Active space size: Orbitals, Electrons #
cas_list = [5,6] # Pick orbitals for CAS space, 1-based indices #
mo = mcscf.sort_mo(mycas, mf.mo_coeff, cas_list) # Just orders the MOs #
mycas.natorb = True

# Run CAS #
mycas.kernel(mo)[0]

# Resulting orbitals printing in .fchk file #
s = mf.get_ovlp()
apost.write_fchk(mol, mycas, molname, s)

# .dm1 and .dm2 files creation ONLY for multiconfig WFs #
apost.write_dm12(mol, mycas, molname)
```

For single-determinant calculations (Hartree-Fock or KS-DFT), only the
`.fchk` file is required for running APOST-3D. The `write_fchk` function is
independent of the WF type, and is called the same way as in the CASSCF
example above.

For illustrative purposes, here's how to create the `.fchk` file from a
Hartree-Fock calculation:

```python
from pyscf import gto, scf
import apost3d as apost

molname = 'HF-Hartree-Fock'
mol=gto.M()
mol.atom='''
H    0.0000000    0.0000000    0.0000000
F    0.0000000    0.0000000    0.9500000
'''
mol.basis='def2tzvp'
mol.spin = 0
mol.charge = 0
mol.cart= False
mol.symmetry = False
mol.verbose=4
mol.max_memory=40000
mol.build()

# Calculate RHF  #
mf = scf.RHF(mol)
mf.kernel()

# Resulting orbitals printing in .fchk file #
s = mf.get_ovlp()
apost.write_fchk(mol, mf, molname, s)
```
