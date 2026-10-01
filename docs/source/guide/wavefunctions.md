# Preparing the wavefunction

APOST-3D reads the wavefunction from a **formatted checkpoint file**
(`jobname.fchk`), the text format of Gaussian, also written by Q-Chem, by
pySCF through the script shipped with APOST-3D, and by converters such as
[IOData](https://github.com/theochem/iodata) for other programs. This page
explains how to produce it, and the extra data some analyses need.

## What the file must contain

The geometry, the basis set and the molecular orbitals (and their
density), which every program above writes by default. Spherical and
Cartesian basis functions up to g are supported. The program is compiled
for up to 3000 basis functions and 350 atoms.

| Analysis | Extra requirement |
|---|---|
| All, single determinant (HF, KS-DFT) | none |
| All, correlated wavefunction (CASSCF, CISD, ...) | the density of the correlated wavefunction in the `.fchk` file (natural orbitals are computed from it) |
| `SPIN`, `ENPART` with a correlated wavefunction | the 1- and 2-RDM files, in a [`# DM` block](../input/dm.md) |
| `ENPART` | the reference energies appended to the `.fchk` file (recommended; see below) |
| `NAO-BASIS` | the NAO transformation file `jobname.nao` |
| `HIRSH`, `HIRSH-IT` | the atomic densities file `densoutput` ([utilities](../tools/utilities.md)) |

## Gaussian

Convert the checkpoint file with `formchk`:

```bash
formchk mol.chk mol.fchk
```

**For ENPART**, add `#P`, `iop(3/33=3)` and `Pop=Full` to the route
section, so that the output file contains the energy components, and then
append them to the `.fchk` file with the utility that matches your
Gaussian version:

```bash
$APOST3D_PATH/utils/get_energy_g16 mol.log >> mol.fchk    # Gaussian 16
$APOST3D_PATH/utils/get_energy mol.log >> mol.fchk        # Gaussian 09
```

**For NAO-BASIS**, the transformation to natural atomic orbitals comes
from the NBO program: add `pop=(full,nboread)` to the route section and
the line `$NBO AONAO=W $END` at the end of the input file. NBO writes the
matrix to `FILE.33`; rename it to `mol.nao`.

**Correlated densities** (MP2, CI, ...): Gaussian writes them after the SCF
density in the `.fchk` file when the calculation uses `density=current`;
select them with `DENS 2` in `# METHOD` (restricted wavefunctions).

## Q-Chem

Ask Q-Chem for the formatted checkpoint file with `GUI = 2` in the
`$rem` section, and add
`QCHEM` to `# METHOD` when the analysis writes `.fchk` files of orbitals
(OSLO, EOS, GEOS).

## pySCF

`utils/apost3d.py` writes the `.fchk` file from a pySCF calculation,
including the reference energies needed by ENPART, and for CASSCF also the
1- and 2-RDM files. Make it importable:

```bash
export PYTHONPATH=$PYTHONPATH:$APOST3D_PATH/utils
```

**Single determinant** (HF or KS-DFT):

```python
from pyscf import gto, scf
import apost3d as apost

molname = 'HF-RHF'
mol = gto.M(atom='''
H    0.0000000    0.0000000    0.0000000
F    0.0000000    0.0000000    0.9500000
''', basis='def2tzvp', charge=0, spin=0)

mf = scf.RHF(mol)
mf.kernel()

apost.write_fchk(mol, mf, molname, mf.get_ovlp())     # writes HF-RHF.fchk
```

**CASSCF**, with the RDMs for SPIN and ENPART:

```python
from pyscf import gto, scf, mcscf
import apost3d as apost

molname = 'HF-CASSCF'
mol = gto.M(atom='''
H    0.0000000    0.0000000    0.0000000
F    0.0000000    0.0000000    0.9500000
''', basis='def2tzvp', charge=0, spin=0)

mf = scf.RHF(mol)
mf.kernel()

mycas = mcscf.CASSCF(mf, 2, 2)                  # 2 orbitals, 2 electrons
mo = mcscf.sort_mo(mycas, mf.mo_coeff, [5, 6])  # active orbitals (1-based)
mycas.natorb = True
mycas.kernel(mo)

apost.write_fchk(mol, mycas, molname, mf.get_ovlp())  # HF-CASSCF.fchk
apost.write_dm12(mol, mycas, molname)                 # HF-CASSCF.dm1, .dm2
```

and in the APOST-3D input:

```text
# METHOD
TFVC
SPIN
DM 2
#
# DM
HF-CASSCF.dm1
HF-CASSCF.dm2
pySCF
#
```

`apost3d.py` does not support basis sets with g functions.
