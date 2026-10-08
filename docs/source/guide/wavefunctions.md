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
$APOST3D_PATH/utils/get_energy mol.log >> mol.fchk        # Gaussian 09 (to be removed)
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

**Correlated densities**: in a gradient job (`JOBTYPE force` or `opt`)
with a correlated method such as MP2, Q-Chem writes the relaxed
correlated density in place of the SCF one, so the default `DENS 1`
analyses that density, while the orbitals in the file are still the SCF
ones. The other density blocks Q-Chem writes are not read, so `DENS 2`
does not apply to Q-Chem files.

## pySCF

`utils/apost3d.py` writes the `.fchk` file from a pySCF calculation,
including the reference energies needed by ENPART, and for CASSCF and FCI
also the 1- and 2-RDM files. It works with pySCF 2.7 and newer (checked up
to 2.14). Make it importable:

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

**FCI**: pass the FCI solver, and the mean-field calculation it started
from as `myhf`:

```python
from pyscf import fci

cisolver = fci.FCI(mf)
cisolver.kernel()
apost.write_fchk(mol, cisolver, 'HF-FCI', mf.get_ovlp(), myhf=mf)
apost.write_dm12(mol, cisolver, 'HF-FCI')               # HF-FCI.dm1, .dm2
```

The last argument of `write_fchk`, the overlap matrix, can be left out.

The `.fchk` file of a closed-shell CASSCF or FCI wavefunction is enough for
the population analysis, EFFAO, EOS and GEOS. For SPIN and ENPART, and for
an open-shell wavefunction in any analysis (its alpha and beta densities
come from the RDMs), give the RDM files in the APOST-3D input:

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

**What `apost3d.py` takes:**

- **SCF**: RHF, UHF, ROHF, RKS, UKS and ROKS, also with density fitting.
- **CASSCF/CASCI and FCI** on RHF or ROHF orbitals, with `write_dm12` for
  the RDMs.
- **CCSD** for closed shells (`cc.CCSD` on RHF): its unrelaxed 1-RDM.
  `write_dm12` writes only the `.dm1` file: there is no CCSD 2-RDM, so no
  SPIN or ENPART.
- **Basis sets**: spherical or Cartesian, up to g functions, with
  pseudopotentials, and with linearly dependent functions removed by pySCF.
- **Dispersion corrections** (D3, D4): their energy is removed from the
  reference electron-electron energy (APOST-3D decomposes the electronic
  energy only), with a note printed.

Correlated wavefunctions built on UHF orbitals (UHF-based CASSCF or FCI,
UCCSD) are refused with a message: APOST-3D takes the alpha and beta
densities of a correlated wavefunction from the RDMs, written in one set
of orbitals, so their spin-resolved analyses would be wrong. Use an RHF-
or ROHF-based CASSCF or FCI for open shells. Any other problem (an object
whose calculation has not been run, h functions) also stops with a
message, and no file is written.
